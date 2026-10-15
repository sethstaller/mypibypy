"""Local, pair-aware genome-wide specificity screening.

What makes this different from every general probe-design tool: a piFISH primary
pair only produces signal when **both** halves of the split probe bind adjacent
target. So a genome hit for one half is not an off-target -- it is a hit that
cannot make a spot by itself. The question worth asking is:

    is there anywhere else in the genome where BOTH halves land close together,
    and would those duplexes still be bound after the washes?

This module screens the two 25-nt target footprints (`left25` / `right25`, as
they read on the target's plus strand -- see `design`), collects every hit for
each half within a mismatch budget, and then joins them: a site is reported only
when both halves hit the same chromosome within a window. A hit for one half
alone lowers the score; it is not a rejection.

Two local backends, both shelled out to with explicit arguments:

* ``blastn -task blastn-short``, with ``-dust no -soft_masking false`` -- low
  complexity filtering destroys 25-nt queries, and a masked query silently
  returns nothing.
* ``bowtie`` (v1) ``-v 3 -a``, which is *exhaustive* end-to-end within 3
  mismatches. When a bowtie index is configured this is the better choice: it
  cannot miss a site inside its mismatch budget, whereas any seed-and-extend
  heuristic can.

Measured recovery of planted true sites, 40 x 25-nt queries per cell, from
``scripts/benchmark_specificity.py`` (rerun it; do not trust this table blind):

    backend                                 0mm    1mm    2mm    3mm    4mm   kept
    blastn-short -dust no -word_size 5    40/40  40/40  40/40  40/40  18/40   1316
    blastn-short -dust no -word_size 7    40/40  40/40  40/40  39/40  16/40    765
    blastn-short -dust no -word_size 11   40/40  40/40  27/40  14/40   5/40    151
    bowtie v1 -v 3 -a                     40/40  40/40  40/40  40/40   0/40    160
    bowtie2 --very-sensitive-local        40/40   0/40   0/40   0/40   0/40     40

Three things set the defaults here. ``bowtie -v 3 -a`` is perfect to the edge of
its budget and returns only the planted sites -- 160 kept alignments, no
background -- so it is preferred whenever an index exists. blastn-short's own
default ``-word_size 7`` already loses a site at 3 mismatches, so the default
here is 5; that costs runtime and ~8x more alignments to sift, which is the
price of not coming back falsely clean. And ``bowtie2 --very-sensitive-local``,
the obvious first choice, is unusable at this query length: a 25-nt query with
one mismatch often has no seed left, which is why it is kept in the benchmark as
a baseline rather than offered as a backend.

Filtering is on alignment geometry (aligned length, substitutions, gaps), never
on E-value. E-value scales with database size, so an E-value cutoff tuned on one
database silently changes meaning on another; with short queries it also throws
away true sites. The BLAST call therefore passes a deliberately permissive
``-evalue`` and does the real filtering here.

DEGRADED MODE IS A FIRST-CLASS RESULT. With no database configured the report
comes back ``NOT_SCREENED`` with a reason, and every consumer must render that
honestly. An earlier version of the order workbook asserted a specificity result
unconditionally, which put an unverified "repeat-free, single-copy" claim onto a
real order sheet. `SpecificityReport.screened` is the only thing that licenses
such a claim.

Configuration (first that resolves wins):

1. arguments to `resolve_backend()` / `screen_pair()`
2. ``$MYPIBYPY_BOWTIE_INDEX`` -- basename of a bowtie v1 index
3. ``$MYPIBYPY_BLAST_DB`` -- basename of a ``makeblastdb`` nucleotide database
"""
from __future__ import annotations

import bisect
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from . import thermo
from .sequence import GenomicInterval

#: Report statuses. `SCREENED` is the only one that licenses a specificity claim.
SCREENED = "SCREENED"
NOT_SCREENED = "NOT_SCREENED"
ERROR = "ERROR"

ENV_BOWTIE_INDEX = "MYPIBYPY_BOWTIE_INDEX"
ENV_BLAST_DB = "MYPIBYPY_BLAST_DB"

#: Mismatch budget per half. bowtie v1 is exhaustive up to 3; beyond that it is
#: not supported by the aligner at all, which is why 3 is the ceiling here too.
DEFAULT_MAX_MISMATCHES = 3

#: Two halves closer than this (bases between their footprints) are treated as
#: one candidate site. 200 nt is far wider than the 1-nt design gap on purpose:
#: a pair landing even loosely together is worth a human look.
DEFAULT_WINDOW_NT = 200

#: Gap at or below which the two halves are "juxtaposed" -- close enough to
#: reconstitute the P2 landing site, i.e. close enough to make a real spot.
#: The design gap is 1 nt; a few nt of slack covers indel-containing alignments.
DEFAULT_JUXTAPOSED_GAP_NT = 10

#: Shortest aligned block worth carrying forward. NOT an arbitrary round
#: number: `thermo` puts a perfect 14-bp duplex at roughly -1 C margin against
#: the most stringent wash at ~50% GC, and a 16-bp one at about +2 C. So below
#: ~15 bp a duplex comes apart regardless of how many mismatches it has, and a
#: shorter alignment cannot become an off-target. This is only a tractability
#: pre-filter -- BLAST emits tens of thousands of 7-12 nt scraps per query that
#: would otherwise go into an O(n^2) pairing step. The real decision is still
#: `thermo`'s, on the aligned sub-duplex. Raise it for a GC-rich target, where a
#: short duplex holds better than this default assumes; the count of hits it
#: removed is reported either way.
DEFAULT_MIN_ALIGNED_NT = 15

#: Score deductions. Starting points calibrated against nothing but judgement --
#: there are three designs in hand, far too few to fit these. Treat the hit
#: lists as the result and the score as a sort key.
PENALTY_SURVIVING_JUXTAPOSED = 1.0
PENALTY_SURVIVING_COLOCALIZED = 0.5
PENALTY_NONSURVIVING_PAIR = 0.1
PENALTY_SURVIVING_SINGLE = 0.02
MAX_SINGLE_PENALTY = 0.3


class SpecificityNotConfigured(RuntimeError):
    """No local database/index is available to screen against."""


# --- Hits ------------------------------------------------------------------

@dataclass(frozen=True)
class Hit:
    """One alignment of one half-probe, in plus-strand genomic coordinates."""

    query: str          # "left" or "right"
    chrom: str
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    strand: str         # "+" or "-"
    mismatches: int     # substitutions inside the aligned block
    gaps: int
    aligned: int        # query bases included in the alignment
    query_len: int
    q_start: int = 1    # 1-based inclusive, on the query
    q_end: int = 0      # 0 means "to the end of the query"

    @property
    def unaligned(self) -> int:
        """Query bases the aligner left out of the alignment.

        BLAST aligns locally and will happily trim a mismatched end, so a 25-nt
        query can come back as a clean 20-nt hit. Those five bases did not
        vanish -- they are five unpaired positions, and the duplex that forms is
        20 bp, not 25. `duplex_seq` is what to hand to `thermo`; this is how much
        got left behind.
        """
        return max(0, self.query_len - self.aligned)

    @property
    def effective_mismatches(self) -> int:
        """Severity on one scale across backends: substitutions, gaps, trimming.

        Used for *ranking and reporting*. It is deliberately not the admission
        filter -- dropping a trimmed alignment would turn a real site into a
        clean report, which is the one direction a specificity screen must never
        fail in. See `screen_pair`.
        """
        return self.mismatches + self.gaps + self.unaligned

    def duplex_seq(self, query_seq: str) -> str:
        """The part of the query that actually paired, for thermodynamics."""
        end = self.q_end or self.query_len
        return query_seq[max(0, self.q_start - 1) : end]

    @property
    def interval(self) -> tuple[str, int, int]:
        return (self.chrom, self.start, self.end)

    def overlaps(self, iv: GenomicInterval) -> bool:
        return (self.chrom == iv.chrom and self.start <= iv.end and iv.start <= self.end)

    def __str__(self) -> str:
        trim = f"+{self.unaligned}nt trimmed" if self.unaligned else "full-length"
        return (f"{self.query} {self.chrom}:{self.start:,}-{self.end:,}({self.strand}) "
                f"{self.mismatches}mm/{self.gaps}gap, {trim}")


@dataclass(frozen=True)
class PairSite:
    """A place where both halves hit close together -- the only thing that signals."""

    chrom: str
    left: Hit
    right: Hit
    gap_nt: int             # bases between the two footprints (design value: 1)
    same_strand: bool
    correct_order: bool     # halves in the order the target geometry requires
    survives_wash: bool     # both duplexes predicted still bound after every step
    min_margin_c: float     # worst per-step Tm margin across both halves
    juxtaposed_gap_nt: int = DEFAULT_JUXTAPOSED_GAP_NT

    @property
    def juxtaposed(self) -> bool:
        """Close enough, and oriented right, to reconstitute the P2 site."""
        return (self.same_strand and self.correct_order
                and 0 <= self.gap_nt <= self.juxtaposed_gap_nt)

    @property
    def total_mismatches(self) -> int:
        return self.left.effective_mismatches + self.right.effective_mismatches

    def __str__(self) -> str:
        kind = "JUXTAPOSED" if self.juxtaposed else "co-localized"
        wash = "survives wash" if self.survives_wash else "melts in wash"
        return (f"{kind} {self.chrom}:{min(self.left.start, self.right.start):,} "
                f"gap {self.gap_nt} nt, {self.total_mismatches} mm total, "
                f"{wash} (worst margin {self.min_margin_c:+.1f} C)")


# --- Backends --------------------------------------------------------------

def _write_fasta(path: Path, queries: dict[str, str]) -> None:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in queries.items()))


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


@dataclass(frozen=True)
class Backend:
    """A local aligner plus the database it searches."""

    tool: str            # "blastn" | "bowtie"
    database: str        # db/index basename
    executable: str
    max_mismatches: int = DEFAULT_MAX_MISMATCHES
    word_size: int = 5   # blastn only; 5, not blastn-short's own 7 -- measured
    threads: int = 1

    def describe(self) -> str:
        if self.tool == "bowtie":
            return (f"bowtie v1 -v {self.max_mismatches} -a (exhaustive within "
                    f"{self.max_mismatches} mismatches) against {self.database}")
        return (f"blastn -task blastn-short -dust no -word_size {self.word_size} "
                f"against {self.database}")

    def search(self, queries: dict[str, str]) -> list[Hit]:
        if self.tool == "blastn":
            return self._blastn(queries)
        if self.tool == "bowtie":
            return self._bowtie(queries)
        raise ValueError(f"unknown backend tool {self.tool!r}")

    # -- blastn ------------------------------------------------------------
    _FMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend "
            "sstart send evalue bitscore qlen")

    def _blastn(self, queries: dict[str, str]) -> list[Hit]:
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "query.fa"
            _write_fasta(fa, queries)
            cmd = [
                self.executable,
                "-task", "blastn-short",
                "-db", self.database,
                "-query", str(fa),
                "-outfmt", self._FMT,
                # Short queries are largely low-complexity by BLAST's reckoning;
                # with DUST on, 25-nt pads come back with no hits at all.
                "-dust", "no",
                "-soft_masking", "false",
                # Permissive on purpose -- see module docstring. Filtering is on
                # alignment geometry, which does not move with database size.
                "-evalue", "1000",
                "-word_size", str(self.word_size),
                "-penalty", "-2",        # -3 (the default) prices mismatched
                "-reward", "1",          # 25mers out of the result set
                "-max_target_seqs", "5000",
                "-num_threads", str(self.threads),
            ]
            proc = _run(cmd)
        if proc.returncode != 0:
            raise RuntimeError(f"blastn failed ({proc.returncode}): {proc.stderr.strip()}")
        return self.parse_blast(proc.stdout, queries)

    @staticmethod
    def parse_blast(tabular: str, queries: dict[str, str]) -> list[Hit]:
        """Parse ``-outfmt 6`` with the column list in ``Backend._FMT``."""
        hits: list[Hit] = []
        for line in tabular.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 13:
                continue
            qid, sid = f[0], f[1]
            length, mism, gaps = int(f[3]), int(f[4]), int(f[5])
            qstart, qend = int(f[6]), int(f[7])
            sstart, send = int(f[8]), int(f[9])
            qlen = int(f[12])
            strand = "+" if sstart <= send else "-"
            start, end = (sstart, send) if strand == "+" else (send, sstart)
            if qstart > qend:            # minus-strand query coordinates
                qstart, qend = qend, qstart
            hits.append(Hit(query=qid, chrom=sid, start=start, end=end, strand=strand,
                            mismatches=mism, gaps=gaps, aligned=qend - qstart + 1,
                            query_len=len(queries.get(qid, "")) or qlen,
                            q_start=qstart, q_end=qend))
        return hits

    # -- bowtie v1 ---------------------------------------------------------
    def _bowtie(self, queries: dict[str, str]) -> list[Hit]:
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "query.fa"
            _write_fasta(fa, queries)
            cmd = [
                self.executable,
                "-f",                                # FASTA input
                "-v", str(self.max_mismatches),      # end-to-end, <= N mismatches
                "-a",                                # every alignment, not just one
                "--sam",
                "--quiet",
                "--threads", str(self.threads),
                "-x", self.database,
                str(fa),
            ]
            proc = _run(cmd)
        if proc.returncode != 0:
            raise RuntimeError(f"bowtie failed ({proc.returncode}): {proc.stderr.strip()}")
        return self.parse_sam(proc.stdout, queries)

    @staticmethod
    def parse_sam(sam: str, queries: dict[str, str]) -> list[Hit]:
        """Parse bowtie's SAM output. Mismatches come from the NM tag."""
        hits: list[Hit] = []
        for line in sam.splitlines():
            if not line or line.startswith("@"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            qid, flag, ref, pos = f[0], int(f[1]), f[2], int(f[3])
            if flag & 0x4 or ref == "*":
                continue
            seq_len = len(queries.get(qid, "")) or len(f[9])
            nm = 0
            for tag in f[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
                    break
            # bowtie v1 -v N is end-to-end: the whole read aligns, always.
            hits.append(Hit(query=qid, chrom=ref, start=pos, end=pos + seq_len - 1,
                            strand="-" if flag & 0x10 else "+",
                            mismatches=nm, gaps=0, aligned=seq_len, query_len=seq_len,
                            q_start=1, q_end=seq_len))
        return hits


def resolve_backend(bowtie_index: str | None = None, blast_db: str | None = None,
                    prefer: str | None = None,
                    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
                    threads: int = 1) -> Backend:
    """Find a usable local backend, or raise `SpecificityNotConfigured`.

    bowtie v1 is preferred when both are available because ``-v N -a`` is
    exhaustive: within its mismatch budget it cannot miss a site. Pass
    ``prefer="blastn"`` to override.
    """
    bowtie_index = bowtie_index or os.environ.get(ENV_BOWTIE_INDEX)
    blast_db = blast_db or os.environ.get(ENV_BLAST_DB)

    options = []
    if bowtie_index:
        options.append(("bowtie", bowtie_index))
    if blast_db:
        options.append(("blastn", blast_db))
    if prefer:
        options.sort(key=lambda o: o[0] != prefer)
    if not options:
        raise SpecificityNotConfigured(
            f"no local database configured: set ${ENV_BOWTIE_INDEX} (bowtie v1 "
            f"index basename) or ${ENV_BLAST_DB} (makeblastdb basename), or pass "
            f"one explicitly. No genome-wide screen has been run."
        )
    problems = []
    for tool, db in options:
        exe = shutil.which(tool)
        if not exe:
            problems.append(f"{tool} not on PATH")
            continue
        if not _database_present(tool, db):
            problems.append(f"{tool} database {db!r} not found")
            continue
        return Backend(tool=tool, database=db, executable=exe,
                       max_mismatches=max_mismatches, threads=threads)
    raise SpecificityNotConfigured("; ".join(problems))


def _database_present(tool: str, db: str) -> bool:
    """Cheap existence check on an index/database basename."""
    p = Path(db)
    if tool == "bowtie":
        return bool(list(p.parent.glob(p.name + "*.ebwt"))
                    or list(p.parent.glob(p.name + "*.ebwtl")))
    return bool(list(p.parent.glob(p.name + ".n*"))
                or list(p.parent.glob(p.name + ".*.n*")))


# --- Pairing ---------------------------------------------------------------

def _wash_outcome(hit: Hit, query_seq: str,
                  steps: tuple[thermo.Conditions, ...] = thermo.ALL_STEPS,
                  penalty_c: float = thermo.DEFAULT_MISMATCH_PENALTY_C,
                  ) -> tuple[bool, float]:
    """(survives every step, worst margin) for the duplex this hit would form."""
    duplex = hit.duplex_seq(query_seq)
    if len(duplex) < 2:
        return False, float("-inf")
    pr = thermo.profile(duplex, hit.mismatches + hit.gaps, steps, penalty_c)
    return pr["survives"], pr["min_margin_c"]


def pair_hits(left_hits: list[Hit], right_hits: list[Hit], left_seq: str,
              right_seq: str, window_nt: int = DEFAULT_WINDOW_NT,
              juxtaposed_gap_nt: int = DEFAULT_JUXTAPOSED_GAP_NT,
              steps: tuple[thermo.Conditions, ...] = thermo.ALL_STEPS,
              penalty_c: float = thermo.DEFAULT_MISMATCH_PENALTY_C) -> list[PairSite]:
    """Join per-half hits into candidate signal-producing sites.

    Both halves hybridize to the target's plus strand with a 1-nt gap, so on an
    off-target locus they must still land on the *same* strand and in the right
    order: left before right when both align to the plus strand, right before
    left when the locus is the reverse complement. Order and strand are recorded
    rather than used as a filter -- a co-localized pair in the wrong orientation
    makes no spot, but it is still worth seeing.
    """
    # Indexed by chromosome and sorted by start so each left hit only looks at
    # the right hits that could possibly be in range. A permissive blastn run on
    # a real genome returns thousands of hits per chromosome, and the naive
    # all-against-all is quadratic in that.
    by_chrom: dict[str, tuple[list[int], list[Hit]]] = {}
    for h in sorted(right_hits, key=lambda h: (h.chrom, h.start)):
        starts, hs = by_chrom.setdefault(h.chrom, ([], []))
        starts.append(h.start)
        hs.append(h)
    # Widest a right hit can be, so the search bound is exact rather than a guess.
    longest_right = max((h.end - h.start + 1 for h in right_hits), default=0)

    sites: list[PairSite] = []
    for lh in left_hits:
        starts, chrom_hits = by_chrom.get(lh.chrom, ([], []))
        lo = bisect.bisect_left(starts, lh.start - window_nt - longest_right - 2)
        hi = bisect.bisect_right(starts, lh.end + window_nt + 2)
        for rh in chrom_hits[lo:hi]:
            if lh.strand == rh.strand == "+":
                gap = rh.start - lh.end - 1
                correct_order = rh.start > lh.start
            elif lh.strand == rh.strand == "-":
                gap = lh.start - rh.end - 1
                correct_order = rh.start < lh.start
            else:
                # Opposite strands: no consistent geometry, but measure the
                # separation anyway so the site can be reported.
                gap = max(rh.start - lh.end - 1, lh.start - rh.end - 1)
                correct_order = False
            if abs(gap) > window_nt:
                continue
            # Evaluate the duplex that actually forms -- the aligned block,
            # carrying its own substitutions -- not the full 25-nt pad. A hit
            # BLAST trimmed to 20 nt makes a 20-bp duplex, and asking thermo
            # about the 25-mer would overstate it.
            l_ok, l_margin = _wash_outcome(lh, left_seq, steps, penalty_c)
            r_ok, r_margin = _wash_outcome(rh, right_seq, steps, penalty_c)
            sites.append(PairSite(
                chrom=lh.chrom, left=lh, right=rh, gap_nt=gap,
                same_strand=lh.strand == rh.strand, correct_order=correct_order,
                survives_wash=l_ok and r_ok, min_margin_c=min(l_margin, r_margin),
                juxtaposed_gap_nt=juxtaposed_gap_nt,
            ))
    sites.sort(key=lambda s: (not s.juxtaposed, not s.survives_wash, s.total_mismatches))
    return sites


# --- Report ----------------------------------------------------------------

@dataclass
class SpecificityReport:
    """The outcome of a screen -- including the outcome "no screen happened"."""

    status: str
    reason: str = ""
    backend: str = ""
    database: str = ""
    max_mismatches: int = DEFAULT_MAX_MISMATCHES
    window_nt: int = DEFAULT_WINDOW_NT
    min_aligned_nt: int = DEFAULT_MIN_ALIGNED_NT
    n_raw_hits: int = 0          # before any filtering
    n_dropped_short: int = 0     # removed by the min_aligned_nt pre-filter
    n_dropped_mismatch: int = 0  # over the mismatch budget
    left_seq: str = ""
    right_seq: str = ""
    left_hits: list[Hit] = field(default_factory=list)
    right_hits: list[Hit] = field(default_factory=list)
    pair_sites: list[PairSite] = field(default_factory=list)
    on_target: GenomicInterval | None = None

    @property
    def screened(self) -> bool:
        """True only when a screen actually ran. Gate every claim on this."""
        return self.status == SCREENED

    @property
    def surviving_sites(self) -> list[PairSite]:
        return [s for s in self.pair_sites if s.survives_wash]

    @property
    def juxtaposed_sites(self) -> list[PairSite]:
        return [s for s in self.pair_sites if s.juxtaposed]

    @property
    def single_half_hits(self) -> list[Hit]:
        """Hits that are not part of any reported pair site: no spot on their own."""
        paired = {id(s.left) for s in self.pair_sites} | {id(s.right) for s in self.pair_sites}
        return [h for h in (*self.left_hits, *self.right_hits) if id(h) not in paired]

    @property
    def surviving_single_half_hits(self) -> list[Hit]:
        """Lone-half hits whose duplex is predicted to hold through the washes."""
        seqs = {"left": self.left_seq, "right": self.right_seq}
        return [h for h in self.single_half_hits
                if seqs.get(h.query) and _wash_outcome(h, seqs[h.query])[0]]

    @property
    def score(self) -> float:
        """1.0 = nothing found. Lower is worse. A sort key, not a probability."""
        if not self.screened:
            return float("nan")
        score = 1.0
        for s in self.pair_sites:
            if s.survives_wash and s.juxtaposed:
                score -= PENALTY_SURVIVING_JUXTAPOSED
            elif s.survives_wash:
                score -= PENALTY_SURVIVING_COLOCALIZED
            else:
                score -= PENALTY_NONSURVIVING_PAIR
        # A lone half cannot recruit P2, so it cannot make a spot. It still costs
        # something -- it is probe sequestered elsewhere, and it is where a
        # second nearby hit would turn into a real site -- but only a little, and
        # only when the duplex would actually hold through the washes.
        singles = PENALTY_SURVIVING_SINGLE * len(self.surviving_single_half_hits)
        score -= min(singles, MAX_SINGLE_PENALTY)
        return max(0.0, score)

    @property
    def verdict(self) -> str:
        """``pass`` / ``review`` / ``fail`` / ``not screened``."""
        if not self.screened:
            return "not screened"
        if any(s.survives_wash and s.juxtaposed for s in self.pair_sites):
            return "fail"
        if self.pair_sites or self.score < 0.9:
            return "review"
        return "pass"

    def summary(self) -> str:
        if self.status == NOT_SCREENED:
            return f"NOT SCREENED -- {self.reason}"
        if self.status == ERROR:
            return f"SCREEN FAILED -- {self.reason}"
        n_surv = len(self.surviving_sites)
        n_jux = len([s for s in self.pair_sites if s.juxtaposed and s.survives_wash])
        return (
            f"{self.backend}: {len(self.left_hits)} left / {len(self.right_hits)} right "
            f"half-hits at <={self.max_mismatches} mismatches; "
            f"{len(self.pair_sites)} off-target site(s) with both halves within "
            f"{self.window_nt} nt, of which {n_surv} predicted to survive the wash "
            f"and {n_jux} juxtaposed closely enough to reconstitute the P2 site "
            f"-- verdict {self.verdict} "
            f"[{self.n_raw_hits} raw alignments; {self.n_dropped_mismatch} over "
            f"budget, {self.n_dropped_short} shorter than {self.min_aligned_nt} nt]"
        )

    def to_dict(self) -> dict:
        return {
            "status": self.status,
            "screened": self.screened,
            "reason": self.reason,
            "backend": self.backend,
            "database": self.database,
            "max_mismatches": self.max_mismatches,
            "window_nt": self.window_nt,
            "n_left_hits": len(self.left_hits),
            "n_right_hits": len(self.right_hits),
            "n_pair_sites": len(self.pair_sites),
            "n_surviving_sites": len(self.surviving_sites),
            "n_single_half_hits": len(self.single_half_hits),
            "n_raw_hits": self.n_raw_hits,
            "n_dropped_short": self.n_dropped_short,
            "n_dropped_mismatch": self.n_dropped_mismatch,
            "min_aligned_nt": self.min_aligned_nt,
            "n_juxtaposed_sites": len(self.juxtaposed_sites),
            "score": self.score,
            "verdict": self.verdict,
            "sites": [str(s) for s in self.pair_sites],
            "summary": self.summary(),
        }


def not_screened(reason: str) -> SpecificityReport:
    """The honest no-op result. Use this instead of omitting a report."""
    return SpecificityReport(status=NOT_SCREENED, reason=reason)


def screen_pair(left_seq: str, right_seq: str, backend: Backend | None = None,
                on_target: GenomicInterval | None = None,
                window_nt: int = DEFAULT_WINDOW_NT,
                max_mismatches: int = DEFAULT_MAX_MISMATCHES,
                juxtaposed_gap_nt: int = DEFAULT_JUXTAPOSED_GAP_NT,
                min_aligned_nt: int = DEFAULT_MIN_ALIGNED_NT,
                steps: tuple[thermo.Conditions, ...] = thermo.ALL_STEPS,
                penalty_c: float = thermo.DEFAULT_MISMATCH_PENALTY_C,
                **backend_kwargs) -> SpecificityReport:
    """Screen one split-probe footprint. Never raises for a missing database.

    `left_seq` / `right_seq` are the two 25-nt footprints as they read on the
    target's plus strand (``Candidate.left.seq`` / ``Candidate.right.seq``), not
    the synthesized probes -- the arms are constant and belong to `crosshyb`.

    `on_target` is the intended footprint; hits overlapping it are dropped so the
    design does not flag itself.
    """
    if backend is None:
        try:
            backend = resolve_backend(max_mismatches=max_mismatches, **backend_kwargs)
        except SpecificityNotConfigured as exc:
            return not_screened(str(exc))

    queries = {"left": left_seq.upper(), "right": right_seq.upper()}
    try:
        hits = backend.search(queries)
    except (RuntimeError, OSError) as exc:
        return SpecificityReport(status=ERROR, reason=str(exc),
                                 backend=backend.describe(), database=backend.database)

    # Admission is on substitutions and gaps only. Trimming is NOT grounds for
    # dropping a hit: a locally-aligned 20/25 nt hit is a real 20-bp duplex, and
    # discarding it would turn a genuine off-target into a clean report. The
    # short-alignment pre-filter exists for tractability and is reported.
    n_raw = len(hits)
    kept = [h for h in hits if h.mismatches + h.gaps <= max_mismatches]
    n_dropped_mismatch = n_raw - len(kept)
    long_enough = [h for h in kept if h.aligned >= min_aligned_nt]
    n_dropped_short = len(kept) - len(long_enough)
    if on_target is not None:
        long_enough = [h for h in long_enough if not h.overlaps(on_target)]
    left_hits = [h for h in long_enough if h.query == "left"]
    right_hits = [h for h in long_enough if h.query == "right"]

    sites = pair_hits(left_hits, right_hits, left_seq, right_seq,
                      window_nt=window_nt, juxtaposed_gap_nt=juxtaposed_gap_nt,
                      steps=steps, penalty_c=penalty_c)
    return SpecificityReport(
        status=SCREENED, backend=backend.describe(), database=backend.database,
        max_mismatches=max_mismatches, window_nt=window_nt,
        min_aligned_nt=min_aligned_nt, n_raw_hits=n_raw,
        n_dropped_short=n_dropped_short, n_dropped_mismatch=n_dropped_mismatch,
        left_seq=queries["left"], right_seq=queries["right"],
        left_hits=left_hits, right_hits=right_hits, pair_sites=sites,
        on_target=on_target,
    )
