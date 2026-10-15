"""End-to-end design orchestration: locus -> fetched sequence -> ranked windows
-> chosen split-probe pair -> specificity screen.

Every `DesignResult` carries a `SpecificityReport`, even when no database is
configured -- in that case it is an explicit NOT_SCREENED, never a missing
field. Consumers that make specificity claims (the order workbook) read
`report.screened` rather than assuming one happened.
"""
from __future__ import annotations

from dataclasses import dataclass

from .design import Candidate, ProbePair, build_pair, scan
from .loci import Locus
from .sequence import fetch_sequence, mask_repeats, rmsk_hits
from .specificity import Backend, SpecificityReport, not_screened, screen_pair
from .systems import AmplifierSystem, get_system


@dataclass
class DesignResult:
    locus: Locus
    system: AmplifierSystem
    region_seq: str
    candidates: list[Candidate]
    pair: ProbePair
    rmsk: list[dict]
    specificity: SpecificityReport

    @property
    def best(self) -> Candidate:
        return self.candidates[0]


def design_locus(locus: Locus, top: int = 5, backend: Backend | None = None,
                 screen: bool = True, **screen_kwargs) -> DesignResult:
    """Design one locus end to end.

    `screen=False` skips the genome screen outright; the result still carries a
    NOT_SCREENED report saying so, because "we chose not to screen" and "we
    screened and it was clean" must never render the same way.
    """
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
    if screen:
        spec = screen_pair(cands[0].left.seq, cands[0].right.seq, backend=backend,
                           on_target=cands[0].interval, **screen_kwargs)
    else:
        spec = not_screened("genome screen skipped by request (screen=False)")
    return DesignResult(locus, system, seq, cands, pair, rmsk, spec)
