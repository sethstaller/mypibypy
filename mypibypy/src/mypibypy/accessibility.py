"""Is the target site on the transcript actually reachable?

An RNA probe can be perfectly specific and still never bind, because its site is
folded into a stem. This failure mode does not exist for DNA targets, which the
protocol denatures, and it is the thing primer-shaped tools get wrong when they
are pointed at FISH: they fold the oligo, which is a different question.

WHAT IS SCORED. The accessibility of the **target site on the transcript**, in
the context of its flanking sequence -- not the folding of the probe oligo. Only
the first predicts whether the probe can reach its site. A probe that folds on
itself is a real but separate problem, and `crosshyb.best_duplex` already looks
at self-complementarity.

WHY A CONTIGUOUS-STRETCH PROBABILITY AND NOT A MEAN. The obvious metric is the
mean per-base unpaired probability over the pad. It is misleading, measurably
so. On a synthetic hairpin with an 8-nt loop, RNAplfold gives the loop bases a
per-base unpaired probability of 0.89 -- they really are unpaired -- while the
probability that any 8 consecutive bases there are *simultaneously* unpaired is
0.008. A mean would call that loop accessible. A probe cannot nucleate in it.
So the headline number here is the probability of an unpaired run long enough to
nucleate, and the mean is reported alongside as context rather than as the
verdict.

    measured on the synthetic control in the tests:
                      P(base unpaired)   P(8-nt run unpaired)
      open flank            0.95                 0.86
      stem                  0.28                 0.046
      hairpin loop          0.89                 0.008   <- the trap

THE TWO PADS ARE SCORED SEPARATELY. They bind independently, and the pair needs
both: one buried pad kills the design however open the other is. So the pair's
accessibility is the *minimum* of the two, never the average.

NO CALIBRATED THRESHOLD EXISTS. Nobody has measured, for this assay, the
accessibility below which a piFISH pair fails. Every cutoff in this module is a
starting point for ranking and reporting, and `ACCESSIBILITY_IS_UNCALIBRATED`
says so where a caller will see it. This is deliberately a ranking input and a
reported attribute, never a hard reject.

Requires ViennaRNA's ``RNAplfold`` on PATH. Absent it, results come back
``NOT_EVALUATED`` with a reason, exactly as `specificity` does with no database
configured -- a fresh clone still runs.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

EVALUATED = "EVALUATED"
NOT_EVALUATED = "NOT_EVALUATED"

#: Said out loud wherever a number from this module is reported.
ACCESSIBILITY_IS_UNCALIBRATED = (
    "accessibility thresholds are UNCALIBRATED for this assay: no piFISH "
    "failure has been tied to a measured accessibility, so treat these as a "
    "ranking input, not a pass/fail criterion"
)

#: RNAplfold's sliding window (-W): the span over which the local partition
#: function is computed. 80 nt is ViennaRNA's own default for mRNA-scale input.
DEFAULT_WINDOW_W = 80

#: Maximum base-pair span (-L). No base may pair with a partner further away
#: than this, which is what bounds how much flanking sequence can matter.
DEFAULT_SPAN_L = 40

#: Length of unpaired run treated as "enough to nucleate". Short oligo
#: hybridisation needs a toehold of roughly this order; 8 nt is a reasonable
#: middle and is deliberately exposed as a parameter, because it is a guess.
DEFAULT_NUCLEATION_NT = 8

#: Flanking transcript included either side of the footprint before folding.
#: Local structure depends on surroundings, so the footprint is never folded in
#: isolation. 150 nt is comfortably more than both W (80) and L (40): every
#: footprint base therefore sees a full sliding window, and every base that
#: could pair with it at all -- nothing can pair beyond L -- is present. Larger
#: flanks cannot change a footprint's local structure under these parameters,
#: they only cost runtime.
DEFAULT_FLANK_NT = 150

#: Only used to put a word next to a number in reports. NOT a threshold.
_BANDS = ((0.50, "open"), (0.20, "partly open"), (0.05, "mostly paired"))


def describe_accessibility(p: float) -> str:
    """A word for a probability, for reports. Explicitly not a cutoff."""
    for floor, label in _BANDS:
        if p >= floor:
            return label
    return "occluded"


class AccessibilityNotAvailable(RuntimeError):
    """RNAplfold is not installed or could not be run."""


@dataclass(frozen=True)
class PadAccessibility:
    """Accessibility of one 25-nt pad's site on the transcript."""

    pad: str                   # "left" or "right"
    start: int                 # 0-based offset into the transcript
    length: int
    nucleation_nt: int
    #: P(some run of `nucleation_nt` bases inside this pad is unpaired), taking
    #: the best such run. The headline number.
    best_nucleation_p: float
    #: Mean P(base unpaired) across the pad. Context only -- see module docstring.
    mean_unpaired: float
    #: Lowest per-base unpaired probability in the pad.
    min_unpaired: float

    @property
    def label(self) -> str:
        return describe_accessibility(self.best_nucleation_p)

    def __str__(self) -> str:
        return (f"{self.pad} pad: nucleation P={self.best_nucleation_p:.3f} "
                f"({self.label}), mean per-base unpaired {self.mean_unpaired:.2f}")


@dataclass
class AccessibilityResult:
    """Accessibility of a 51-nt footprint, or an explicit statement that it was
    not evaluated. `evaluated` is the only thing that licenses using a number."""

    status: str
    reason: str = ""
    left: PadAccessibility | None = None
    right: PadAccessibility | None = None
    gap_unpaired: float | None = None
    tool: str = ""
    window_w: int = DEFAULT_WINDOW_W
    span_l: int = DEFAULT_SPAN_L
    flank_nt: int = DEFAULT_FLANK_NT
    nucleation_nt: int = DEFAULT_NUCLEATION_NT
    flanks_available: tuple[int, int] = (0, 0)
    caveat: str = field(default=ACCESSIBILITY_IS_UNCALIBRATED)

    @property
    def evaluated(self) -> bool:
        return self.status == EVALUATED and self.left is not None and self.right is not None

    @property
    def pair_nucleation_p(self) -> float | None:
        """The pair's accessibility: the WORSE pad, not the average.

        Both halves must bind for the pair to signal, so an open left pad cannot
        compensate for a buried right one.
        """
        if not self.evaluated:
            return None
        return min(self.left.best_nucleation_p, self.right.best_nucleation_p)

    @property
    def limiting_pad(self) -> PadAccessibility | None:
        if not self.evaluated:
            return None
        return min((self.left, self.right), key=lambda p: p.best_nucleation_p)

    def summary(self) -> str:
        if not self.evaluated:
            return f"accessibility NOT EVALUATED -- {self.reason}"
        return (
            f"{self.tool}: target-site accessibility, limiting pad is the "
            f"{self.limiting_pad.pad} at P={self.pair_nucleation_p:.3f} "
            f"({self.limiting_pad.label}) for a {self.nucleation_nt}-nt unpaired "
            f"run; left {self.left.best_nucleation_p:.3f} / right "
            f"{self.right.best_nucleation_p:.3f} "
            f"(-W {self.window_w} -L {self.span_l}, {self.flanks_available[0]}/"
            f"{self.flanks_available[1]} nt flanks). {self.caveat}"
        )

    def to_dict(self) -> dict:
        return {
            "status": self.status,
            "evaluated": self.evaluated,
            "reason": self.reason,
            "tool": self.tool,
            "nucleation_nt": self.nucleation_nt,
            "pair_nucleation_p": self.pair_nucleation_p,
            "left_nucleation_p": self.left.best_nucleation_p if self.left else None,
            "right_nucleation_p": self.right.best_nucleation_p if self.right else None,
            "left_mean_unpaired": self.left.mean_unpaired if self.left else None,
            "right_mean_unpaired": self.right.mean_unpaired if self.right else None,
            "gap_unpaired": self.gap_unpaired,
            "limiting_pad": self.limiting_pad.pad if self.limiting_pad else None,
            "window_w": self.window_w,
            "span_l": self.span_l,
            "flank_nt": self.flank_nt,
            "caveat": self.caveat,
        }


def not_evaluated(reason: str) -> AccessibilityResult:
    """The honest no-op. Use instead of omitting a result or defaulting to 1.0."""
    return AccessibilityResult(status=NOT_EVALUATED, reason=reason)


def rnaplfold_path() -> str | None:
    """Where RNAplfold is, honouring ``$MYPIBYPY_RNAPLFOLD`` first."""
    override = os.environ.get("MYPIBYPY_RNAPLFOLD")
    if override:
        return override if Path(override).is_file() else None
    return shutil.which("RNAplfold")


def parse_lunp(text: str) -> list[list[float | None]]:
    """Parse an RNAplfold ``*_lunp`` table.

    Row *i* (1-based) holds the probability that the run of length *l* **ending
    at** position *i* is unpaired, for l = 1..u. ``NA`` appears where l > i.
    Returned 0-indexed by position: ``rows[i-1][l-1]``.
    """
    rows: list[list[float | None]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        try:
            index = int(fields[0])
        except ValueError:
            continue
        values: list[float | None] = []
        for cell in fields[1:]:
            cell = cell.strip()
            if not cell or cell == "NA":
                values.append(None)
            else:
                try:
                    values.append(float(cell))
                except ValueError:
                    values.append(None)
        while len(rows) < index:
            rows.append([])
        rows[index - 1] = values
    return rows


def _run_plfold(seq: str, window_w: int, span_l: int, nucleation_nt: int,
                executable: str) -> list[list[float | None]]:
    """Fold `seq` locally and return the unpaired-probability table.

    RNAplfold writes its output files into the working directory, named after
    the FASTA header, so it is run in a temporary directory with a fixed header
    rather than being trusted to put things anywhere in particular.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        (tmpdir / "win.fa").write_text(f">win\n{seq}\n")
        proc = subprocess.run(
            [executable, "-W", str(min(window_w, len(seq))), "-L", str(span_l),
             "-u", str(nucleation_nt)],
            cwd=tmpdir, input=f">win\n{seq}\n", capture_output=True, text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise AccessibilityNotAvailable(
                f"RNAplfold failed ({proc.returncode}): {proc.stderr.strip()[:200]}")
        lunp = tmpdir / "win_lunp"
        if not lunp.is_file():
            raise AccessibilityNotAvailable(
                f"RNAplfold produced no unpaired-probability table "
                f"(stderr: {proc.stderr.strip()[:200]})")
        return parse_lunp(lunp.read_text())


def _best_run(rows, start: int, length: int, u: int) -> float:
    """Best P(run of `u` unpaired) among runs lying wholly inside [start, start+length).

    Row *i* column *u* describes the run ending at *i*, i.e. covering
    ``i-u+1 .. i``; a run is only counted when both ends are inside the pad, so
    accessibility is never borrowed from the sequence next door.
    """
    best = 0.0
    for end in range(start + u, start + length + 1):      # 1-based inclusive end
        row = rows[end - 1] if end - 1 < len(rows) else []
        if len(row) < u:
            continue
        value = row[u - 1]
        if value is not None:
            best = max(best, value)
    return best


def _per_base(rows, start: int, length: int) -> tuple[float, float]:
    """(mean, min) per-base unpaired probability over a 1-based inclusive span."""
    values = []
    for pos in range(start, start + length):
        row = rows[pos - 1] if pos - 1 < len(rows) else []
        if row and row[0] is not None:
            values.append(row[0])
    if not values:
        return float("nan"), float("nan")
    return sum(values) / len(values), min(values)


def site_accessibility(transcript: str, offset: int, left_len: int = 25,
                       gap_len: int = 1, right_len: int = 25,
                       window_w: int = DEFAULT_WINDOW_W,
                       span_l: int = DEFAULT_SPAN_L,
                       nucleation_nt: int = DEFAULT_NUCLEATION_NT,
                       flank_nt: int = DEFAULT_FLANK_NT,
                       executable: str | None = None) -> AccessibilityResult:
    """Accessibility of the footprint at `offset` within `transcript`.

    `offset` is 0-based into the transcript. Never raises for a missing
    RNAplfold: the result comes back `NOT_EVALUATED` with a reason, because a
    design pipeline must not stop because an optional tool is absent, and must
    equally not silently score an unevaluated site as accessible.
    """
    executable = executable or rnaplfold_path()
    if not executable:
        return not_evaluated(
            "RNAplfold not found. Install ViennaRNA (conda install -c conda-forge "
            "-c bioconda viennarna) or set $MYPIBYPY_RNAPLFOLD. No accessibility "
            "has been evaluated for this site")

    footprint = left_len + gap_len + right_len
    if offset < 0 or offset + footprint > len(transcript):
        return not_evaluated(
            f"footprint at offset {offset} does not fit inside a "
            f"{len(transcript)} nt transcript")

    # Fold the footprint in context. Flanks are clipped at the transcript ends,
    # and how much was actually available is reported -- a site near the 5' or 3'
    # end genuinely has less context, and that should be visible rather than
    # quietly padded.
    lo = max(0, offset - flank_nt)
    hi = min(len(transcript), offset + footprint + flank_nt)
    region = transcript[lo:hi].upper().replace("T", "U")
    if any(b not in "ACGUN" for b in region):
        return not_evaluated("transcript region contains non-RNA characters")

    try:
        rows = _run_plfold(region, window_w, span_l, nucleation_nt, executable)
    except (AccessibilityNotAvailable, OSError) as exc:
        return not_evaluated(str(exc))

    # Positions within the folded region, 1-based for the lunp table.
    left_start = offset - lo + 1
    right_start = left_start + left_len + gap_len
    gap_pos = left_start + left_len

    pads = []
    for name, start, length in (("left", left_start, left_len),
                                ("right", right_start, right_len)):
        mean_p, min_p = _per_base(rows, start, length)
        pads.append(PadAccessibility(
            pad=name, start=(offset if name == "left" else offset + left_len + gap_len),
            length=length, nucleation_nt=nucleation_nt,
            best_nucleation_p=_best_run(rows, start, length, nucleation_nt),
            mean_unpaired=mean_p, min_unpaired=min_p,
        ))

    gap_row = rows[gap_pos - 1] if gap_pos - 1 < len(rows) else []
    return AccessibilityResult(
        status=EVALUATED, tool=f"RNAplfold -W {window_w} -L {span_l} -u {nucleation_nt}",
        left=pads[0], right=pads[1],
        gap_unpaired=(gap_row[0] if gap_row and gap_row[0] is not None else None),
        window_w=window_w, span_l=span_l, flank_nt=flank_nt,
        nucleation_nt=nucleation_nt,
        flanks_available=(offset - lo, hi - (offset + footprint)),
    )


#: Weight applied to accessibility when it is used as a ranking input. Chosen so
#: that a fully occluded site is penalised comparably to a badly unbalanced
#: GC/Tm pair, which is the worst thing `design._score` already knows how to
#: express. Like every number here, uncalibrated.
DEFAULT_ACCESSIBILITY_WEIGHT = 2.0


def accessibility_penalty(result: AccessibilityResult,
                          weight: float = DEFAULT_ACCESSIBILITY_WEIGHT) -> float:
    """Ranking penalty from accessibility. 0.0 when not evaluated.

    Zero, not a guess, when RNAplfold did not run: an unevaluated site must rank
    exactly as it did before this module existed, so that installing ViennaRNA
    is the only thing that changes an ordering. The penalty is driven by the
    limiting pad, for the reason in the module docstring.
    """
    p = result.pair_nucleation_p
    if p is None:
        return 0.0
    return weight * (1.0 - p)
