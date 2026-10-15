"""Target loci, loaded from a file rather than hardcoded.

A locus is an input to the tool, not part of it: which regions a lab is probing
is experiment-specific and often unpublished. So loci live in a JSON file you
supply, the same way amplifier systems do (see `systems`). Resolution order:

1. An explicit path passed to `load_loci()` (the CLI's ``--loci``).
2. ``$MYPIBYPY_LOCI``.
3. ``loci.local.json`` in the current directory or any parent (gitignored).
4. ``loci.example.json`` next to this package -- a published, citable region so
   the tool runs out of the box.

Create your own from the shipped example:

    cp loci.example.json loci.local.json

RNA targets live under an ``"rna_targets"`` key and are identified by accession,
since a transcript is a molecule rather than an interval:

    {"rna_targets": [
      {"name": "cds", "gene": "GENEA", "accession": "NM_000546.6",
       "system": "cy5", "system_antisense": "cy3",
       "exon_boundaries": [120, 400, 910]}
    ]}

Fusion junctions live in the same file under a ``"fusions"`` key. A junction is
not a genomic interval -- it is a sequence that exists only in the sample -- so
it is supplied as sequence plus a breakpoint offset rather than coordinates:

    {"fusions": [
      {"name": "demo", "five_prime_gene": "GENEA", "three_prime_gene": "GENEB",
       "sequence": "ACGT...", "breakpoint": 120,
       "system": "cy5", "system_antisense": "cy3"}
    ]}

`breakpoint` is the 0-based index of the first base contributed by the 3'
partner, so ``sequence[:breakpoint]`` is 5'-partner sequence. See
`design.fusion_footprint` for how the 1-nt gap is placed on it.

File format -- a list, or an object with a "loci" key holding one:

    {"loci": [
      {"name": "promoter", "gene": "Actb", "description": "positive control",
       "build": "mm10", "chrom": "chr5", "start": 142906255, "end": 142906754,
       "system": "cy5"}
    ]}

`start`/`end` are 1-based inclusive, matching `GenomicInterval` and UCSC.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

from .sequence import GenomicInterval

LOCAL_FILENAME = "loci.local.json"
EXAMPLE_FILENAME = "loci.example.json"
ENV_VAR = "MYPIBYPY_LOCI"


@dataclass(frozen=True)
class Locus:
    name: str
    gene: str
    description: str
    interval: GenomicInterval
    system: str          # amplifier system key, e.g. "cy5" / "cy3"


class LociNotConfigured(RuntimeError):
    """No loci file could be found, or it was unusable."""


def _search_upward(filename: str, start: Path | None = None) -> Path | None:
    here = (start or Path.cwd()).resolve()
    for d in (here, *here.parents):
        if (d / filename).is_file():
            return d / filename
    return None


def resolve_loci_path(path: str | Path | None = None, start: Path | None = None) -> Path:
    """Where loci would be loaded from, without loading them."""
    if path:
        p = Path(path).expanduser()
        if not p.is_file():
            raise LociNotConfigured(f"--loci file not found: {p}")
        return p
    env = os.environ.get(ENV_VAR)
    if env:
        p = Path(env).expanduser()
        if not p.is_file():
            raise LociNotConfigured(f"${ENV_VAR} points at a missing file: {p}")
        return p
    local = _search_upward(LOCAL_FILENAME, start)
    if local:
        return local
    example = Path(__file__).resolve().parent / EXAMPLE_FILENAME
    if example.is_file():
        return example
    raise LociNotConfigured(
        f"no loci file found. Pass --loci, set ${ENV_VAR}, or create "
        f"{LOCAL_FILENAME} (copy {EXAMPLE_FILENAME} to start)."
    )


def _locus_from_dict(d: dict, src: Path) -> Locus:
    missing = {"name", "gene", "build", "chrom", "start", "end", "system"} - set(d)
    if missing:
        raise LociNotConfigured(f"{src}: locus {d.get('name', '?')!r} missing {sorted(missing)}")
    try:
        start, end = int(d["start"]), int(d["end"])
    except (TypeError, ValueError) as exc:
        raise LociNotConfigured(f"{src}: locus {d['name']!r} has non-integer coordinates") from exc
    if start > end:
        raise LociNotConfigured(f"{src}: locus {d['name']!r} has start > end")
    return Locus(
        name=str(d["name"]), gene=str(d["gene"]),
        description=str(d.get("description", "")),
        interval=GenomicInterval(str(d["build"]), str(d["chrom"]), start, end),
        system=str(d["system"]),
    )


def load_loci(path: str | Path | None = None) -> list[Locus]:
    """Load target loci. `path` overrides the search order."""
    src = resolve_loci_path(path)
    try:
        raw = json.loads(src.read_text())
    except json.JSONDecodeError as exc:
        raise LociNotConfigured(f"{src} is not valid JSON: {exc}") from exc
    entries = raw.get("loci") if isinstance(raw, dict) else raw
    if not isinstance(entries, list) or not entries:
        raise LociNotConfigured(f"{src} has no loci")
    return [_locus_from_dict(d, src) for d in entries]


def is_example(path: str | Path | None = None) -> bool:
    """True when the loaded loci are the shipped demonstration region."""
    return resolve_loci_path(path).name == EXAMPLE_FILENAME


# --- Fusion junctions ------------------------------------------------------

@dataclass(frozen=True)
class FusionJunction:
    """A breakpoint, with enough partner sequence either side to probe it.

    `sequence` is the fusion as supplied (transcript or genomic), and
    `breakpoint` is the 0-based index of the first base from the 3' partner.
    Fusions carry two systems because a junction is ordered as both strands --
    the lab puts them on different channels so the strands can be told apart.
    """

    name: str
    five_prime_gene: str
    three_prime_gene: str
    sequence: str
    breakpoint: int
    system: str                          # sense pair
    system_antisense: str | None = None  # antisense pair; defaults to `system`
    gap_side: str = "five_prime"
    description: str = ""

    @property
    def five_prime_seq(self) -> str:
        return self.sequence[: self.breakpoint]

    @property
    def three_prime_seq(self) -> str:
        return self.sequence[self.breakpoint :]

    @property
    def antisense_system(self) -> str:
        return self.system_antisense or self.system

    def __str__(self) -> str:
        return (f"{self.five_prime_gene}::{self.three_prime_gene} ({self.name}), "
                f"breakpoint at {self.breakpoint} of {len(self.sequence)} nt")


def _fusion_from_dict(d: dict, src: Path) -> FusionJunction:
    required = {"name", "five_prime_gene", "three_prime_gene", "sequence",
                "breakpoint", "system"}
    missing = required - set(d)
    if missing:
        raise LociNotConfigured(
            f"{src}: fusion {d.get('name', '?')!r} missing {sorted(missing)}")
    seq = str(d["sequence"]).upper().replace(" ", "")
    try:
        bp = int(d["breakpoint"])
    except (TypeError, ValueError) as exc:
        raise LociNotConfigured(
            f"{src}: fusion {d['name']!r} has a non-integer breakpoint") from exc
    if not 0 < bp < len(seq):
        raise LociNotConfigured(
            f"{src}: fusion {d['name']!r} breakpoint {bp} is outside its "
            f"{len(seq)} nt sequence")
    bad = set(seq) - set("ACGTN")
    if bad:
        raise LociNotConfigured(
            f"{src}: fusion {d['name']!r} sequence has non-DNA characters {sorted(bad)}")
    return FusionJunction(
        name=str(d["name"]), five_prime_gene=str(d["five_prime_gene"]),
        three_prime_gene=str(d["three_prime_gene"]), sequence=seq, breakpoint=bp,
        system=str(d["system"]), system_antisense=d.get("system_antisense"),
        gap_side=str(d.get("gap_side", "five_prime")),
        description=str(d.get("description", "")),
    )


def load_fusions(path: str | Path | None = None) -> list[FusionJunction]:
    """Load fusion junctions. Returns [] when the file defines none."""
    src = resolve_loci_path(path)
    try:
        raw = json.loads(src.read_text())
    except json.JSONDecodeError as exc:
        raise LociNotConfigured(f"{src} is not valid JSON: {exc}") from exc
    entries = raw.get("fusions") if isinstance(raw, dict) else None
    if not entries:
        return []
    if not isinstance(entries, list):
        raise LociNotConfigured(f"{src}: 'fusions' must be a list")
    return [_fusion_from_dict(d, src) for d in entries]


# --- RNA targets -----------------------------------------------------------

@dataclass(frozen=True)
class RnaTarget:
    """A transcript to probe, identified by accession rather than coordinates.

    `exon_boundaries` are 0-based offsets into the transcript at which a new exon
    begins, when they are known. They are optional and purely informational: a
    footprint spanning one is often *desirable* for isoform specificity, so this
    exists to report what a design does, not to forbid anything. Leave it empty
    and junction-spanning simply goes unreported.

    `scan_start` / `scan_end` restrict probe placement to part of the transcript
    (0-based, end-exclusive; `scan_end` of 0 means "to the end").
    """

    name: str
    gene: str
    accession: str
    system: str
    system_antisense: str | None = None
    description: str = ""
    exon_boundaries: tuple[int, ...] = ()
    scan_start: int = 0
    scan_end: int = 0

    @property
    def antisense_system(self) -> str:
        return self.system_antisense or self.system

    def __str__(self) -> str:
        return f"{self.gene} ({self.accession}), target {self.name}"


def _rna_from_dict(d: dict, src: Path) -> RnaTarget:
    missing = {"name", "gene", "accession", "system"} - set(d)
    if missing:
        raise LociNotConfigured(
            f"{src}: RNA target {d.get('name', '?')!r} missing {sorted(missing)}")
    try:
        boundaries = tuple(int(b) for b in d.get("exon_boundaries", ()))
    except (TypeError, ValueError) as exc:
        raise LociNotConfigured(
            f"{src}: RNA target {d['name']!r} has non-integer exon_boundaries") from exc
    if any(b < 0 for b in boundaries):
        raise LociNotConfigured(
            f"{src}: RNA target {d['name']!r} has a negative exon boundary")
    start, end = int(d.get("scan_start", 0)), int(d.get("scan_end", 0))
    if end and end <= start:
        raise LociNotConfigured(
            f"{src}: RNA target {d['name']!r} has scan_end <= scan_start")
    return RnaTarget(
        name=str(d["name"]), gene=str(d["gene"]), accession=str(d["accession"]),
        system=str(d["system"]), system_antisense=d.get("system_antisense"),
        description=str(d.get("description", "")),
        exon_boundaries=tuple(sorted(boundaries)), scan_start=start, scan_end=end,
    )


def load_rna_targets(path: str | Path | None = None) -> list[RnaTarget]:
    """Load RNA targets from the ``rna_targets`` key. [] when none are defined."""
    src = resolve_loci_path(path)
    try:
        raw = json.loads(src.read_text())
    except json.JSONDecodeError as exc:
        raise LociNotConfigured(f"{src} is not valid JSON: {exc}") from exc
    entries = raw.get("rna_targets") if isinstance(raw, dict) else None
    if not entries:
        return []
    if not isinstance(entries, list):
        raise LociNotConfigured(f"{src}: 'rna_targets' must be a list")
    return [_rna_from_dict(d, src) for d in entries]
