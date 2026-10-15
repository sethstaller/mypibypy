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
