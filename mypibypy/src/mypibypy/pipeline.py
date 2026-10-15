"""End-to-end design orchestration: locus -> fetched sequence -> ranked windows
-> chosen split-probe pair."""
from __future__ import annotations

from dataclasses import dataclass

from .design import Candidate, ProbePair, build_pair, scan
from .loci import Locus
from .sequence import fetch_sequence, mask_repeats, rmsk_hits
from .systems import AmplifierSystem, get_system


@dataclass
class DesignResult:
    locus: Locus
    system: AmplifierSystem
    region_seq: str
    candidates: list[Candidate]
    pair: ProbePair
    rmsk: list[dict]

    @property
    def best(self) -> Candidate:
        return self.candidates[0]


def design_locus(locus: Locus, top: int = 5) -> DesignResult:
    system = get_system(locus.system)
    seq = fetch_sequence(locus.interval)
    rmsk = rmsk_hits(locus.interval)
    # Soft-mask repeats so probe pads are never placed in repetitive DNA.
    masked = mask_repeats(seq, locus.interval, rmsk)
    cands = scan(masked, locus.interval, top=top)
    if not cands:
        raise ValueError(f"no valid repeat-free 51-nt window found in {locus.name}")
    pair = build_pair(
        locus.name, system, cands[0],
        p1l_name=f"P1L_{locus.gene}_{locus.name}",
        p1r_name=f"P1R_{locus.gene}_{locus.name}",
    )
    return DesignResult(locus, system, seq, cands, pair, rmsk)
