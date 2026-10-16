#!/usr/bin/env python3
"""Does the agent loop reach the same conclusions as the verified manual path?

The lab's own instruction: "run probe sets you already validated by hand through
the automated loop and check that it reaches your conclusions." A loop that
confidently reaches a *different* answer on a known-good design is a failure, not
a variation -- so this compares, per design, rather than just reporting that the
loop ran.

For each verified design it:
  1. runs the deterministic pipeline, which is the answer of record;
  2. runs the agent loop over the same target with the same starting parameters;
  3. compares the chosen footprint, the pair count, the screen verdict, and
     whether the agent needed the model or relaxed anything.

Reads real targets at run time from the gitignored local files and from an order
workbook supplied out of band. NOTHING identifying is printed: no sequences, no
coordinates, no gene names -- designs are reported by index. Run:

    export SSL_CERT_FILE=/path/to/roots.pem          # if behind a TLS proxy
    MYPIBYPY_FUSION_WORKBOOK=/path/to/fusion.xlsx \
        python scripts/validate_agent.py --bowtie-index ~/genomes/hg38/hg38
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mypibypy.agent import (  # noqa: E402
    Budget,
    DesignAgent,
    DesignParameters,
    ScriptedClient,
    Situation,
)
from mypibypy.agent.loop import classify  # noqa: E402
from mypibypy.agent.tools import filter_candidates, window_stats  # noqa: E402
from mypibypy.design import FOOTPRINT, pairs_from_footprint, scan  # noqa: E402
from mypibypy.loci import load_loci  # noqa: E402
from mypibypy.pipeline import DesignResult, design_locus  # noqa: E402
from mypibypy.sequence import fetch_sequence, mask_repeats, rmsk_hits  # noqa: E402
from mypibypy.specificity import Backend, not_screened, screen_pair  # noqa: E402
from mypibypy.systems import get_system, systems  # noqa: E402


def backend_or_none(index: str | None):
    if not index:
        return None
    exe = shutil.which("bowtie")
    if not exe:
        print("  (bowtie not on PATH -- screens will report NOT_SCREENED)")
        return None
    return Backend(tool="bowtie", database=index, executable=exe, threads=8)


# --- genomic locus ---------------------------------------------------------

def validate_locus(locus, label: str, backend, params: DesignParameters) -> dict:
    """Deterministic design vs agent run, for one genomic locus."""
    print(f"\n=== {label} (genomic) ===")

    # 1. The answer of record.
    deterministic = design_locus(locus, backend=backend, screen=backend is not None)
    det_footprint = deterministic.best.interval
    det_verdict = deterministic.specificity.verdict
    print(f"  deterministic: footprint offset {deterministic.best.offset}, "
          f"{deterministic.n_pairs} pair(s), screen verdict {det_verdict}")

    # 2. The same target through the agent. Sequence is fetched once and reused
    #    so the comparison isolates the loop, not network variation.
    seq = fetch_sequence(locus.interval)
    rmsk = rmsk_hits(locus.interval)
    masked = mask_repeats(seq, locus.interval, rmsk)
    system = get_system(locus.system)

    def design_fn(p: DesignParameters):
        ranked = scan(masked, locus.interval, top=len(masked))
        kept = filter_candidates(ranked, p)
        if not kept:
            empty = DesignResult(locus, system, seq, ranked[:1],
                                 _stub_pair(locus, system, ranked[0]), rmsk,
                                 not_screened("no window met the constraints"))
            empty.pairs = []
            empty.pair_screens = []
            empty.pair_accessibility = []
            return empty, ranked
        result = design_locus(locus, backend=backend, screen=backend is not None,
                              n_pairs=p.n_pairs, min_spacing=p.min_spacing)
        return result, ranked

    # A proceed decision is scripted in so that IF the loop consults the model,
    # it gets a plausible answer rather than falling back. A known-good design
    # should not reach it at all -- model_calls in the table below is the check.
    agent = DesignAgent(design_fn, client=ScriptedClient([
        {"action": "proceed", "reason": "scripted: evidence looks acceptable",
         "primary_signal": "specificity"}]),
                        budget=Budget(max_turns=4, max_relaxations=3),
                        known_systems=tuple(sorted(systems())))
    result = agent.run(params, target_label=label, request="validate against manual")

    agent_offset = result.design.best.offset if result.design else None
    same = agent_offset == deterministic.best.offset
    print(f"  agent        : status {result.status}, "
          f"footprint offset {agent_offset}, "
          f"{result.record.turns_used} turn(s), "
          f"{result.record.relaxations_used} relaxation(s), "
          f"{result.record.model_calls} model call(s)")
    print(f"  SAME FOOTPRINT: {same}")
    print(f"  why: {result.message[:150]}")
    return {
        "design": label, "kind": "genomic",
        "deterministic_offset": deterministic.best.offset,
        "agent_offset": agent_offset,
        "same_footprint": same,
        "deterministic_verdict": det_verdict,
        "agent_status": result.status,
        "model_calls": result.record.model_calls,
        "relaxations": result.record.relaxations_used,
        "agreed": same and result.status == "succeeded",
    }


def _stub_pair(locus, system, cand):
    from mypibypy.design import build_pair
    return build_pair(locus.name, system, cand, "P1L", "P1R")


# --- fusion junctions ------------------------------------------------------

def read_fusion_footprints(path: str) -> dict[str, str]:
    """{sheet: 51-nt footprint} from an order workbook. Structural, not hardcoded."""
    from openpyxl import load_workbook
    wb = load_workbook(path, data_only=True)
    out = {}
    for sheet in wb.sheetnames:
        for row in wb[sheet].iter_rows(values_only=True):
            cells = [c for c in row if c is not None]
            if len(cells) >= 2 and isinstance(cells[0], str) and isinstance(cells[1], str):
                key, val = cells[0].strip().upper(), re.sub(r"\s", "", cells[1]).upper()
                if (len(val) == FOOTPRINT and set(val) <= set("ACGTN")
                        and "SENSE" in key and "ANTISENSE" not in key):
                    out[sheet] = val
    return out


def validate_fusion(index: int, footprint: str, backend,
                    params: DesignParameters) -> dict:
    """A junction is single-pair by construction; the agent must simply agree."""
    label = f"fusion design #{index}"
    print(f"\n=== {label} ===")
    keys = sorted(systems())
    sense_sys, anti_sys = get_system(keys[0]), get_system(keys[-1])
    pairs = pairs_from_footprint(f"junction{index}", footprint, sense_sys, anti_sys)
    sense = pairs[0]
    cand = sense.candidate

    if backend is not None:
        spec = screen_pair(cand.left.seq, cand.right.seq, backend=backend)
    else:
        spec = not_screened("no index supplied to the validator")
    print(f"  deterministic: 1 pair per strand, screen verdict {spec.verdict}")

    locus_like = type("L", (), {"gene": "x", "name": f"junction{index}",
                                "interval": cand.interval, "system": sense_sys.name,
                                "description": ""})()

    def design_fn(p: DesignParameters):
        result = DesignResult(locus_like, sense_sys, footprint, [cand], sense,
                              [], spec)
        return result, [cand]

    agent = DesignAgent(design_fn, client=ScriptedClient([
        {"action": "proceed", "reason": "scripted: evidence looks acceptable",
         "primary_signal": "specificity"}]),
                        budget=Budget(max_turns=3, max_relaxations=3),
                        known_systems=tuple(keys))
    result = agent.run(params, target_label=label, request="validate against manual")
    agreed = result.status == "succeeded" and result.design.pair.p1l == sense.p1l
    print(f"  agent        : status {result.status}, "
          f"{result.record.turns_used} turn(s), "
          f"{result.record.relaxations_used} relaxation(s), "
          f"{result.record.model_calls} model call(s)")
    print(f"  SAME PRIMARIES: {result.design.pair.p1l == sense.p1l}")
    print(f"  why: {result.message[:150]}")
    return {
        "design": label, "kind": "fusion",
        "deterministic_verdict": spec.verdict,
        "agent_status": result.status,
        "same_primaries": result.design.pair.p1l == sense.p1l,
        "model_calls": result.record.model_calls,
        "relaxations": result.record.relaxations_used,
        "agreed": agreed,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bowtie-index", default=os.environ.get("MYPIBYPY_BOWTIE_INDEX"))
    ap.add_argument("--loci", default=None)
    ap.add_argument("--gc", nargs=2, type=float, default=[0.40, 0.60])
    ap.add_argument("--tm", nargs=2, type=float, default=[55.0, 62.0])
    args = ap.parse_args()

    backend = backend_or_none(args.bowtie_index)
    params = DesignParameters(gc_min=args.gc[0], gc_max=args.gc[1],
                              tm_min=args.tm[0], tm_max=args.tm[1])
    rows = []

    # Genomic loci whose build matches the supplied index. A mouse locus
    # screened against a human index would be a meaningless comparison, so it is
    # skipped rather than quietly run.
    want_build = "hg38" if backend and "hg38" in (args.bowtie_index or "") else None
    for i, locus in enumerate(load_loci(args.loci), start=1):
        if want_build and locus.interval.build != want_build:
            print(f"\n=== genomic design #{i} skipped: build "
                  f"{locus.interval.build} does not match the supplied index ===")
            continue
        rows.append(validate_locus(locus, f"genomic design #{i}", backend, params))

    workbook = os.environ.get("MYPIBYPY_FUSION_WORKBOOK")
    if workbook:
        for i, (_, footprint) in enumerate(
                sorted(read_fusion_footprints(workbook).items()), start=1):
            rows.append(validate_fusion(i, footprint, backend, params))
    else:
        print("\n(set $MYPIBYPY_FUSION_WORKBOOK to include the fusion designs)")

    print("\n" + "=" * 72)
    print("AGREEMENT WITH THE VERIFIED MANUAL PATH")
    print("=" * 72)
    print(f"{'design':<22} {'agreed':>7} {'status':>10} {'model':>6} {'relax':>6}  "
          f"same design")
    for r in rows:
        same = r.get("same_footprint", r.get("same_primaries"))
        print(f"{r['design']:<22} {str(r['agreed']):>7} {r['agent_status']:>10} "
              f"{r['model_calls']:>6} {r['relaxations']:>6}  {same}")
    disagreed = [r for r in rows if not r["agreed"]]
    print(f"\n{len(rows) - len(disagreed)}/{len(rows)} designs agree with the "
          f"manual conclusion.")
    if disagreed:
        print("DISAGREEMENTS -- investigate before trusting the loop:")
        for r in disagreed:
            print(f"  {r['design']}: {r}")
    return 1 if disagreed else 0


if __name__ == "__main__":
    raise SystemExit(main())
