#!/usr/bin/env python3
"""Measure how many *true* off-target sites each specificity backend recovers.

A specificity screen that misses sites is worse than no screen, because it comes
back clean. So the backends are measured rather than trusted: plant sites with a
known number of mismatches at known coordinates in a synthetic background, run
each backend, and count how many planted sites come back.

This exists because the obvious choice failed badly. ``bowtie2
--very-sensitive-local`` recovered 7/40 true sites for 25-nt queries carrying a
single mismatch, and 0/40 at two or more: bowtie2 is built for reads long enough
to seed, and a 25-nt query with a mismatch often has no seed left. That number is
the baseline to beat.

Usage (needs the aligners on PATH -- see README for the conda environment):

    python scripts/benchmark_specificity.py                 # synthetic queries
    python scripts/benchmark_specificity.py --pads FILE     # one 25-nt pad per line

`--pads` lets the measurement run against real footprints without those
sequences ever entering the repository: pass a file outside version control.

CAVEAT ON THE BACKGROUND. The decoy is random sequence, which has none of a real
genome's repeat structure, and it is megabases rather than gigabases. Recovery
measured here is an upper bound on a real genome for a heuristic backend. It is
*not* an upper bound for ``bowtie -v N -a``, which is exhaustive by construction
and whose recovery is therefore a property of the algorithm, not the database.
The deliberately permissive ``-evalue`` in `specificity.Backend` is what keeps
the BLAST numbers from being flattered by the small database -- filtering is on
alignment geometry, which does not move with database size.
"""
from __future__ import annotations

import argparse
import random
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mypibypy.specificity import (  # noqa: E402
    DEFAULT_MAX_MISMATCHES,
    DEFAULT_MIN_ALIGNED_NT,
    Backend,
)

BASES = "ACGT"
PAD_LEN = 25


def random_seq(n: int, rng: random.Random, gc: float = 0.5) -> str:
    """Random sequence with roughly the requested GC fraction."""
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT")
                   for _ in range(n))


def spike(seq: str, n_mismatches: int, rng: random.Random) -> str:
    """Substitute `n_mismatches` interior bases. Ends are left alone on purpose.

    Terminal mismatches are both easier for an aligner to absorb (local
    alignment just trims them) and far less destabilizing at the bench, so
    putting them interior is the harder, more honest test.
    """
    if n_mismatches == 0:
        return seq
    positions = rng.sample(range(2, len(seq) - 2), n_mismatches)
    out = list(seq)
    for p in positions:
        out[p] = rng.choice([b for b in BASES if b != out[p]])
    return "".join(out)


def build_decoy(pads: list[str], mismatch_levels: list[int], rng: random.Random,
                background_nt: int, n_chroms: int) -> tuple[dict[str, str], list[dict]]:
    """Background sequence with one planted site per (pad, mismatch level).

    Sites are spaced far apart so a recovered hit is unambiguously attributable,
    and so no two plants accidentally form a pair site.
    """
    per_chrom = background_nt // n_chroms
    chroms = {f"decoy{i+1}": list(random_seq(per_chrom, rng)) for i in range(n_chroms)}
    plants: list[dict] = []
    slots = []
    for name in chroms:
        # 2 kb apart: wider than any pairing window used by the screen.
        for pos in range(1000, per_chrom - PAD_LEN - 1000, 2000):
            slots.append((name, pos))
    rng.shuffle(slots)
    need = len(pads) * len(mismatch_levels)
    if need > len(slots):
        raise SystemExit(f"background too small: need {need} slots, have {len(slots)}")
    for i, pad in enumerate(pads):
        for k in mismatch_levels:
            name, pos = slots.pop()
            variant = spike(pad, k, rng)
            chroms[name][pos : pos + PAD_LEN] = list(variant)
            plants.append({"query": f"pad{i:03d}", "pad": pad, "mismatches": k,
                           "chrom": name, "start": pos + 1,
                           "end": pos + PAD_LEN, "variant": variant})
    return {k: "".join(v) for k, v in chroms.items()}, plants


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")


def recovered(hits, plants, slack: int = 5) -> dict[int, tuple[int, int]]:
    """Per mismatch level: (recovered, planted)."""
    index: dict[str, list] = {}
    for h in hits:
        index.setdefault(h.query, []).append(h)
    out: dict[int, list[int]] = {}
    for p in plants:
        got = any(
            h.chrom == p["chrom"]
            and abs(h.start - p["start"]) <= slack + PAD_LEN
            and abs(h.end - p["end"]) <= slack + PAD_LEN
            for h in index.get(p["query"], ())
        )
        tally = out.setdefault(p["mismatches"], [0, 0])
        tally[1] += 1
        tally[0] += int(got)
    return {k: tuple(v) for k, v in sorted(out.items())}


def bowtie2_baseline(index: str, queries: dict[str, str], tmp: Path) -> list:
    """The measured baseline: bowtie2 --very-sensitive-local on 25-nt queries."""
    exe = shutil.which("bowtie2")
    if not exe:
        return []
    fa = tmp / "q.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in queries.items()))
    proc = subprocess.run(
        [exe, "--very-sensitive-local", "-f", "-a", "--quiet", "-x", index, "-U", str(fa)],
        capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        print(f"  bowtie2 failed: {proc.stderr.strip()[:200]}", file=sys.stderr)
        return []
    return Backend.parse_sam(proc.stdout, queries)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pads", help="file of real 25-nt footprints, one per line "
                                   "(kept outside the repo)")
    ap.add_argument("--n-pads", type=int, default=40)
    ap.add_argument("--max-mismatches", type=int, default=4)
    ap.add_argument("--background-nt", type=int, default=4_000_000)
    ap.add_argument("--chroms", type=int, default=4)
    ap.add_argument("--word-sizes", type=int, nargs="+", default=[4, 5, 7, 11])
    ap.add_argument("--seed", type=int, default=20251015)
    ap.add_argument("--workdir", default=None, help="where to build the decoy index")
    ap.add_argument("--report-raw", action="store_true",
                    help="score recovery against the full --max-mismatches budget "
                         "instead of the screen's default budget")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    if args.pads:
        pads = [ln.strip().upper() for ln in Path(args.pads).read_text().splitlines()
                if ln.strip() and not ln.startswith("#")]
        pads = [p for p in pads if len(p) == PAD_LEN][: args.n_pads]
        source = f"{len(pads)} real footprints from {args.pads}"
    else:
        pads = [random_seq(PAD_LEN, rng) for _ in range(args.n_pads)]
        source = f"{len(pads)} synthetic footprints (seed {args.seed})"
    if not pads:
        raise SystemExit("no usable 25-nt pads")

    levels = list(range(args.max_mismatches + 1))
    tmp = Path(args.workdir) if args.workdir else Path("./decoy_bench")
    tmp.mkdir(parents=True, exist_ok=True)
    fasta = tmp / "decoy.fa"

    print(f"Building decoy: {args.background_nt:,} nt in {args.chroms} sequences, "
          f"{source}, mismatch levels {levels}")
    chroms, plants = build_decoy(pads, levels, rng, args.background_nt, args.chroms)
    write_fasta(fasta, chroms)
    queries = {f"pad{i:03d}": p for i, p in enumerate(pads)}

    runs: list[tuple[str, list]] = []

    if shutil.which("makeblastdb") and shutil.which("blastn"):
        subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
                        "-out", str(tmp / "decoy")], capture_output=True, check=True)
        for ws in args.word_sizes:
            b = Backend(tool="blastn", database=str(tmp / "decoy"),
                        executable=shutil.which("blastn"), word_size=ws)
            runs.append((f"blastn-short -dust no -word_size {ws}", b.search(queries)))
    else:
        print("  blastn/makeblastdb not found -- skipped", file=sys.stderr)

    if shutil.which("bowtie-build") and shutil.which("bowtie"):
        subprocess.run(["bowtie-build", "-q", str(fasta), str(tmp / "decoy_bt1")],
                       capture_output=True, check=True)
        for v in (2, 3):
            b = Backend(tool="bowtie", database=str(tmp / "decoy_bt1"),
                        executable=shutil.which("bowtie"), max_mismatches=v)
            runs.append((f"bowtie v1 -v {v} -a", b.search(queries)))
    else:
        print("  bowtie v1 not found -- skipped", file=sys.stderr)

    if shutil.which("bowtie2-build") and shutil.which("bowtie2"):
        subprocess.run(["bowtie2-build", "-q", str(fasta), str(tmp / "decoy_bt2")],
                       capture_output=True, check=True)
        runs.append(("bowtie2 --very-sensitive-local (baseline)",
                     bowtie2_baseline(str(tmp / "decoy_bt2"), queries, tmp)))

    # Recovery is measured on the hits the *screen* keeps, not on raw aligner
    # output: `screen_pair` admits a hit on substitutions and gaps, then drops
    # alignments too short to form a duplex that survives any wash.
    budget = args.max_mismatches if args.report_raw else DEFAULT_MAX_MISMATCHES
    width = max(len(n) for n, _ in runs)
    header = "  ".join(f"{k}mm".rjust(7) for k in levels)
    print(f"\nTrue-site recovery (planted sites found / planted), n={len(pads)} per cell")
    print(f"Hits kept = alignments with <= {budget} substitutions+gaps and at "
          f"least {DEFAULT_MIN_ALIGNED_NT} nt aligned, i.e. exactly what "
          f"specificity.screen_pair carries into pairing.")
    print(f"{'backend'.ljust(width)}  {header}      raw     kept")
    print("-" * (width + len(header) + 20))
    for name, hits in runs:
        kept = [h for h in hits if h.mismatches + h.gaps <= budget
                and h.aligned >= DEFAULT_MIN_ALIGNED_NT]
        tab = recovered(kept, plants)
        cells = []
        for k in levels:
            got, total = tab.get(k, (0, 0))
            cells.append(f"{got}/{total}".rjust(7))
        print(f"{name.ljust(width)}  {'  '.join(cells)}  "
              f"{len(hits):>7}  {len(kept):>7}")
    print("\n'raw' is every alignment the aligner emitted; 'kept' is what survives the "
          "screen's geometry filter. A backend that recovers everything by reporting "
          "everything is no use either, so read recovery and 'kept' together.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
