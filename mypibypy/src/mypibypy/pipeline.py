"""End-to-end design orchestration: locus -> fetched sequence -> ranked windows
-> chosen split-probe pair -> specificity screen.

Every `DesignResult` carries a `SpecificityReport`, even when no database is
configured -- in that case it is an explicit NOT_SCREENED, never a missing
field. Consumers that make specificity claims (the order workbook) read
`report.screened` rather than assuming one happened.
"""
from __future__ import annotations

from dataclasses import dataclass

from .design import (
    PADS_DIRECT,
    Candidate,
    ProbePair,
    arm_separation,
    arms_meet_at_junction,
    build_fusion_pairs,
    build_pair,
    check_partner_specificity,
    fusion_footprint,
    scan,
)
from .loci import FusionJunction, Locus
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


def design_fusion(junction: FusionJunction, backend=None, screen: bool = True,
                  pad_assignment: str = PADS_DIRECT,
                  **screen_kwargs) -> list[DesignResult]:
    """Design a fusion junction: one sense result and one antisense result.

    The junction has no coordinates in the reference, so there is no
    RepeatMasker step -- `rmsk` is empty, and the workbook says 0 repeats
    because none were looked for.

    On what the genome screen means here: each pad is, by construction, a real
    25-mer from its own parent gene, so single-half hits are guaranteed and
    expected. They are not a defect. What matters is that no site anywhere has
    *both* halves together -- which is exactly the pair-aware question
    `specificity` asks, and exactly why a per-probe screen is useless for a
    junction. `check_partner_specificity` adds the complementary local check
    that neither pad is present in the *wrong* partner.
    """
    sense_system = get_system(junction.system)
    anti_system = get_system(junction.antisense_system)
    _, footprint = fusion_footprint(junction.sequence, junction.breakpoint,
                                    junction.gap_side)
    partner_problems = check_partner_specificity(
        footprint, junction.five_prime_seq, junction.three_prime_seq)

    pairs = build_fusion_pairs(
        junction.name, junction.sequence, junction.breakpoint,
        sense_system, anti_system, gap_side=junction.gap_side,
        pad_assignment=pad_assignment,
    )

    results: list[DesignResult] = []
    for pair in pairs:
        cand = pair.candidate
        strand = "antisense" if "antisense" in pair.p1l_name else "sense"
        flags = list(cand.flags)
        flags += [f"partner: {p}" for p in partner_problems]
        if not arms_meet_at_junction(pair, footprint) and not any(
                "loop back" in f for f in flags):
            sep = arm_separation(pair, footprint)
            flags.append(
                f"arms anchored {sep} nt apart rather than either side of the "
                f"breakpoint, so they must loop back to present P2 together. "
                f"The lab's hand-made pairs are built this way and do give "
                f"signal; relative efficiency has not been compared")
        cand.flags = flags
        locus = Locus(
            name=f"{junction.name}_{strand}",
            gene=f"{junction.five_prime_gene}::{junction.three_prime_gene}",
            description=(junction.description
                         or f"fusion junction, {strand} strand").strip(),
            interval=cand.interval, system=pair.system.name,
        )
        if screen:
            spec = screen_pair(cand.left.seq, cand.right.seq, backend=backend,
                               **screen_kwargs)
        else:
            spec = not_screened("genome screen skipped by request (screen=False)")
        results.append(DesignResult(locus, pair.system, junction.sequence,
                                    [cand], pair, rmsk=[], specificity=spec))
    return results
