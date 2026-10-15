"""mypibypy -- piFISH split-probe design.

Design split-probe (P1L/P1R) primary probes for piFISH targets and emit an order
workbook that drops into an existing validated P2/P3/P4 amplifier pool.
"""
from .design import Candidate, ProbePair, build_pair, scan
from .loci import Locus, load_loci
from .sequence import GenomicInterval, fetch_sequence, revcomp
from .systems import AmplifierSystem, get_system, is_example, loaded_from, systems

__version__ = "0.1.0"

# `SYSTEM_1` / `SYSTEM_2` are deliberately not re-exported here: importing them
# would load the systems file at import time, and `import mypibypy` must work
# before any systems are configured. Reach them via `systems()` or
# `get_system("cy5")`, or as lazy attributes of `mypibypy.systems`.
__all__ = [
    "AmplifierSystem", "Candidate", "GenomicInterval", "Locus", "ProbePair",
    "build_pair", "fetch_sequence", "get_system", "is_example", "load_loci",
    "loaded_from",
    "revcomp", "scan", "systems",
]
