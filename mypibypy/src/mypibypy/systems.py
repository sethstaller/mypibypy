"""piFISH amplifier systems, loaded from a local file.

A "system" is a self-contained signal channel: a pair of constant arms (landing
pads appended to the split primary probes), the secondary (P2) and tertiary (P3)
amplifiers, and the fluorophore-labeled signal oligo (P4).

piFISH cascade:  target -> P1L + P1R (split primary) -> P2 -> P3 -> P4 (signal)

Mechanism: P1L = arm_L + revcomp(target_left25); P1R = revcomp(target_right25) +
arm_R. When both bind the target with their 25-nt pads (separated by a 1-nt gap),
the two arms are juxtaposed at the junction and together form the contiguous
landing site for P2's central anneal region -- a coincidence requirement that
gives the system its specificity.

Arm architecture (see `crosshyb` for the checks that depend on it):

    arm_left  = [ top_L: 14 nt ][ mid_L: 8 nt ]
    arm_right = [ mid_R: 8 nt ][ top_R: 14 nt ]
    P2 anneal region (29 nt) = rc(top_R) + 1 nt + rc(top_L)

NO REAL SEQUENCES LIVE IN THIS REPOSITORY. Amplifier sequences are lab-specific
and often unpublished, so they are loaded at runtime from a file you keep
outside version control. Resolution order:

1. ``$MYPIBYPY_SYSTEMS`` -- explicit path to a systems JSON file.
2. ``systems.local.json`` in the current directory or any parent.
3. ``systems.example.json`` next to this package, which contains SYNTHETIC
   placeholder sequences of the correct architecture. They let the tests and a
   dry run work; they will not produce working probes.

Create your own from the shipped example:

    cp systems.example.json systems.local.json   # then paste in real sequences

`systems.local.json` is gitignored. Check `loaded_from()` and `is_example()`
before trusting any output -- the order workbook refuses to mark probes as
orderable while example systems are loaded.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

LOCAL_FILENAME = "systems.local.json"
EXAMPLE_FILENAME = "systems.example.json"
ENV_VAR = "MYPIBYPY_SYSTEMS"


@dataclass(frozen=True)
class AmplifierSystem:
    name: str           # internal/legacy name
    arm_left: str       # appended to the 5' end of P1L
    arm_right: str      # appended to the 3' end of P1R
    p2: str             # secondary amplifier oligo
    p3: str             # tertiary amplifier oligo
    p4_name: str        # signal oligo name
    p4: str             # signal oligo (carries 5' fluorophore mod)
    fluorophore: str    # human label, e.g. "Cy5"
    example: bool = False   # True = synthetic placeholder, not orderable


class SystemsNotConfigured(RuntimeError):
    """No systems file could be found."""


def _search_upward(filename: str, start: Path | None = None) -> Path | None:
    here = (start or Path.cwd()).resolve()
    for d in (here, *here.parents):
        candidate = d / filename
        if candidate.is_file():
            return candidate
    return None


def _example_path() -> Path:
    return Path(__file__).resolve().parent / EXAMPLE_FILENAME


def resolve_systems_path(start: Path | None = None) -> Path:
    """Where systems would be loaded from, without loading them."""
    env = os.environ.get(ENV_VAR)
    if env:
        p = Path(env).expanduser()
        if not p.is_file():
            raise SystemsNotConfigured(f"${ENV_VAR} points at a missing file: {p}")
        return p
    local = _search_upward(LOCAL_FILENAME, start)
    if local:
        return local
    example = _example_path()
    if example.is_file():
        return example
    raise SystemsNotConfigured(
        f"no systems file found. Set ${ENV_VAR}, or create {LOCAL_FILENAME} "
        f"in this directory or a parent (copy {EXAMPLE_FILENAME} to start)."
    )


def load_systems(path: str | Path | None = None) -> dict[str, AmplifierSystem]:
    """Load the system registry. `path` overrides the search order."""
    src = Path(path) if path else resolve_systems_path()
    try:
        raw = json.loads(src.read_text())
    except json.JSONDecodeError as exc:
        raise SystemsNotConfigured(f"{src} is not valid JSON: {exc}") from exc
    entries = raw.get("systems")
    if not isinstance(entries, dict) or not entries:
        raise SystemsNotConfigured(f"{src} has no 'systems' object")
    is_example = src.name == EXAMPLE_FILENAME
    out: dict[str, AmplifierSystem] = {}
    for key, d in entries.items():
        missing = {f for f in ("name", "arm_left", "arm_right", "p2", "p3", "p4_name",
                               "p4", "fluorophore")} - set(d)
        if missing:
            raise SystemsNotConfigured(f"{src}: system {key!r} missing {sorted(missing)}")
        out[key] = AmplifierSystem(
            name=d["name"], arm_left=d["arm_left"].upper(), arm_right=d["arm_right"].upper(),
            p2=d["p2"].upper(), p3=d["p3"].upper(), p4_name=d["p4_name"], p4=d["p4"],
            fluorophore=d["fluorophore"], example=bool(d.get("example", is_example)),
        )
    return out


_CACHE: dict[str, AmplifierSystem] | None = None
_CACHE_PATH: Path | None = None


def systems(reload: bool = False) -> dict[str, AmplifierSystem]:
    """The system registry, loaded once and cached."""
    global _CACHE, _CACHE_PATH
    if _CACHE is None or reload:
        _CACHE_PATH = resolve_systems_path()
        _CACHE = load_systems(_CACHE_PATH)
    return _CACHE


def loaded_from() -> Path | None:
    """Path the cached registry came from, or None if nothing is loaded yet."""
    return _CACHE_PATH


def is_example() -> bool:
    """True when the loaded systems are synthetic placeholders.

    Callers that produce orderable output MUST check this.
    """
    return any(s.example for s in systems().values())


def get_system(key: str) -> AmplifierSystem:
    reg = systems()
    k = key.lower()
    if k not in reg:
        raise KeyError(f"unknown system {key!r}; choose from {sorted(reg)}")
    return reg[k]


def __getattr__(name: str):
    """Lazy module attributes so importing this module never needs a systems file.

    `SYSTEMS` is the registry dict. `SYSTEM_1` / `SYSTEM_2` are kept as aliases
    for the first two registered systems, for older call sites.
    """
    if name == "SYSTEMS":
        return systems()
    if name in ("SYSTEM_1", "SYSTEM_2"):
        reg = list(systems().values())
        idx = 0 if name.endswith("1") else 1
        if idx >= len(reg):
            raise AttributeError(
                f"{name} requires at least {idx + 1} configured systems; "
                f"{loaded_from()} defines {len(reg)}"
            )
        return reg[idx]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
