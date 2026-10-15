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

from .design import PADS_AS_ORDERED, PADS_DIRECT
from .loci import (
    Locus,
    is_example,
    load_fusions,
    load_loci,
    load_rna_targets,
    resolve_loci_path,
)
from .pipeline import design_fusion, design_locus, design_rna
from .sequence import GenomicInterval
from .specificity import SpecificityNotConfigured, resolve_backend
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
        print(f"  specificity: {res.specificity.summary()}")
        if res.accessibility.evaluated:
            print(f"  accessibility: {res.accessibility.summary()}")


def _specificity_options(ap):
    """Shared screening options. Absent a database the run still succeeds, but
    the workbook says in so many words that no screen was performed."""
    ap.add_argument("--bowtie-index", help="bowtie v1 index basename "
                                          f"(default: $MYPIBYPY_BOWTIE_INDEX)")
    ap.add_argument("--blast-db", help="makeblastdb nucleotide database basename "
                                      f"(default: $MYPIBYPY_BLAST_DB)")
    ap.add_argument("--no-screen", action="store_true",
                    help="skip the genome specificity screen entirely")


def _backend_from(args):
    """Resolve a screening backend, or None to let `screen_pair` degrade."""
    if args.no_screen:
        return None
    try:
        return resolve_backend(bowtie_index=args.bowtie_index, blast_db=args.blast_db)
    except SpecificityNotConfigured as exc:
        print(f"WARNING: no specificity screen will be run -- {exc}")
        return None


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
    _specificity_options(b)

    f = sub.add_parser("fusion", help="design every fusion junction in a loci file")
    f.add_argument("--loci", help="path to a loci JSON file with a 'fusions' key")
    f.add_argument("-o", "--out", required=True)
    f.add_argument("--pad-assignment", default=PADS_DIRECT,
                   choices=[PADS_DIRECT, PADS_AS_ORDERED],
                   help=f"antisense pad-to-arm assignment. {PADS_DIRECT} "
                        f"(default) anchors both arms beside the breakpoint, as "
                        f"the sense pair does; {PADS_AS_ORDERED} reproduces the "
                        f"lab's hand-made convention, which anchors them at the "
                        f"footprint's outer ends so they loop back. Both give "
                        f"signal at the bench; they have not been compared")
    _specificity_options(f)

    r = sub.add_parser("rna", help="design every RNA target in a loci file")
    r.add_argument("--loci", help="path to a loci JSON file with an 'rna_targets' key")
    r.add_argument("-o", "--out", required=True)
    r.add_argument("--transcriptome-index",
                   help="bowtie v1 index of a TRANSCRIPTOME (see pipeline."
                        "RNA_SCREEN_POLICY for why this, not a genome)")
    r.add_argument("--also-genome-index",
                   help="additionally screen against a genome index. Only needed "
                        "when the protocol denatures DNA; a standard RNA-FISH "
                        "protocol does not, so genomic DNA stays double-stranded "
                        "and the probe cannot invade it")
    r.add_argument("--no-screen", action="store_true",
                   help="skip the specificity screen entirely")

    d = sub.add_parser("design", help="design one locus given on the command line")
    d.add_argument("--build", required=True)
    d.add_argument("--region", required=True, help="chrom:start-end (1-based inclusive)")
    d.add_argument("--system", required=True, help="amplifier system key, e.g. cy5")
    d.add_argument("--gene", default="gene")
    d.add_argument("--name", default="locus")
    d.add_argument("-o", "--out", required=True)
    _specificity_options(d)

    args = ap.parse_args(argv)

    if args.cmd == "rna":
        targets = load_rna_targets(args.loci)
        if not targets:
            raise SystemExit(
                f"no RNA targets defined in {resolve_loci_path(args.loci)} "
                "(add an 'rna_targets' key -- see loci.example.json)")
        warnings = _warn_if_example(args.loci)
        backend = genome_backend = None
        if not args.no_screen:
            try:
                backend = resolve_backend(bowtie_index=args.transcriptome_index)
            except SpecificityNotConfigured as exc:
                print(f"WARNING: no transcriptome screen will be run -- {exc}")
            if args.also_genome_index:
                genome_backend = resolve_backend(bowtie_index=args.also_genome_index)
        results = [design_rna(t, backend=backend, screen=not args.no_screen,
                              also_genome_backend=genome_backend) for t in targets]
        note = (f"mypibypy RNA design of {len(targets)} transcript target(s) from "
                f"{resolve_loci_path(args.loci).name}. Screened against a "
                f"transcriptome; see pipeline.RNA_SCREEN_POLICY.")
    elif args.cmd == "fusion":
        fusions = load_fusions(args.loci)
        if not fusions:
            raise SystemExit(
                f"no fusions defined in {resolve_loci_path(args.loci)} "
                "(add a 'fusions' key -- see loci.example.json)")
        warnings = _warn_if_example(args.loci)
        backend = _backend_from(args)
        results = []
        for j in fusions:
            results += design_fusion(j, backend=backend, screen=not args.no_screen,
                                     pad_assignment=args.pad_assignment)
        if args.pad_assignment == PADS_AS_ORDERED:
            warnings.append(
                "built with the as-ordered antisense convention, which anchors "
                "the antisense arms at the footprint's outer ends so they must "
                "loop back to present P2. That is how the lab's hand-made pairs "
                "were built and they do give signal, but the two conventions "
                "have never been compared side by side")
        note = (f"mypibypy fusion design of {len(fusions)} junction(s) from "
                f"{resolve_loci_path(args.loci).name}, "
                f"pad assignment {args.pad_assignment}.")
    elif args.cmd == "batch":
        loci = load_loci(args.loci)
        warnings = _warn_if_example(args.loci)
        backend = _backend_from(args)
        results = [design_locus(L, backend=backend, screen=not args.no_screen)
                   for L in loci]
        note = (f"mypibypy batch design of {len(loci)} locus/loci from "
                f"{resolve_loci_path(args.loci).name}.")
    else:
        known = sorted(systems())
        if args.system.lower() not in known:
            raise SystemExit(f"unknown --system {args.system!r}; configured systems: {known}")
        locus = Locus(args.name, args.gene, args.name,
                      _parse_region(args.build, args.region), args.system)
        warnings = _warn_if_example()
        backend = _backend_from(args)
        results = [design_locus(locus, backend=backend, screen=not args.no_screen)]
        note = f"mypibypy design: {locus.interval}."

    # An unscreened design must carry that fact into the workbook's provenance
    # note, not only into the specificity cell.
    unscreened = [r.locus.name for r in results if not r.specificity.screened]
    if unscreened:
        warnings.append(
            f"no genome-wide specificity screen was run for: {', '.join(unscreened)}")
    for w in warnings:
        note += f" WARNING: {w}."
    write_workbook(args.out, results, note)
    _report(results)
    for w in warnings:
        print(f"\nWARNING: {w}")
    print(f"\nWROTE: {args.out}")


if __name__ == "__main__":
    main()
