"""Command-line interface for mypibypy.

Examples
--------
Design every locus in a loci file (see `loci` for the format and search order):
    python -m mypibypy batch --loci loci.local.json -o out.xlsx

Design one locus given on the command line:
    python -m mypibypy design --build mm10 --region chr5:142906255-142906754 \
        --system cy5 --gene Actb --name promoter -o out.xlsx

With no --loci and no loci.local.json, `batch` runs against the shipped example
(the published mouse Actb positive control) and says so.
"""
from __future__ import annotations

import argparse
import re

from .loci import Locus, is_example, load_loci, resolve_loci_path
from .pipeline import design_locus
from .sequence import GenomicInterval
from .systems import is_example as systems_are_example
from .systems import loaded_from as systems_loaded_from
from .systems import systems
from .workbook import write_workbook

_REGION = re.compile(r"^(\w+):([\d,]+)-([\d,]+)$")


def _parse_region(build: str, s: str) -> GenomicInterval:
    m = _REGION.match(s)
    if not m:
        raise SystemExit(f"bad --region {s!r}; expected chrom:start-end, e.g. chr5:142906255-142906754")
    chrom, a, b = m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if a > b:
        raise SystemExit(f"bad --region {s!r}: start is greater than end")
    return GenomicInterval(build, chrom, a, b)


def _report(results):
    for res in results:
        p = res.pair
        print(f"\n{res.locus.gene} {res.locus.name} -> {res.system.fluorophore}: "
              f"footprint {res.best.interval}  (rmsk hits: {len(res.rmsk)})")
        print(f"  {p.p1l_name}: {p.p1l}")
        print(f"  {p.p1r_name}: {p.p1r}")


def _warn_if_example(loci_path=None) -> list[str]:
    """Placeholder inputs must never be mistaken for an orderable design."""
    warnings = []
    systems()  # force resolution so loaded_from() is populated
    if systems_are_example():
        warnings.append(
            f"amplifier systems are the SYNTHETIC EXAMPLES ({systems_loaded_from()}); "
            "these probes are NOT orderable -- see README on systems.local.json"
        )
    if is_example(loci_path):
        warnings.append(
            f"target loci are the shipped example ({resolve_loci_path(loci_path)}); "
            "pass --loci or create loci.local.json for your own targets"
        )
    return warnings


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mypibypy", description="piFISH split-probe design")
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("batch", help="design every locus in a loci file")
    b.add_argument("--loci", help="path to a loci JSON file (default: search order, see docs)")
    b.add_argument("-o", "--out", required=True)

    d = sub.add_parser("design", help="design one locus given on the command line")
    d.add_argument("--build", required=True)
    d.add_argument("--region", required=True, help="chrom:start-end (1-based inclusive)")
    d.add_argument("--system", required=True, help="amplifier system key, e.g. cy5")
    d.add_argument("--gene", default="gene")
    d.add_argument("--name", default="locus")
    d.add_argument("-o", "--out", required=True)

    args = ap.parse_args(argv)

    if args.cmd == "batch":
        loci = load_loci(args.loci)
        warnings = _warn_if_example(args.loci)
        results = [design_locus(L) for L in loci]
        note = (f"mypibypy batch design of {len(loci)} locus/loci from "
                f"{resolve_loci_path(args.loci).name}.")
    else:
        known = sorted(systems())
        if args.system.lower() not in known:
            raise SystemExit(f"unknown --system {args.system!r}; configured systems: {known}")
        locus = Locus(args.name, args.gene, args.name,
                      _parse_region(args.build, args.region), args.system)
        warnings = _warn_if_example()
        results = [design_locus(locus)]
        note = f"mypibypy design: {locus.interval}."

    for w in warnings:
        note += f" WARNING: {w}."
    write_workbook(args.out, results, note)
    _report(results)
    for w in warnings:
        print(f"\nWARNING: {w}")
    print(f"\nWROTE: {args.out}")


if __name__ == "__main__":
    main()
