"""End-to-end design orchestration: locus -> fetched sequence -> ranked windows
-> chosen split-probe pair -> specificity screen.

Every `DesignResult` carries a `SpecificityReport`, even when no database is
configured -- in that case it is an explicit NOT_SCREENED, never a missing
field. Consumers that make specificity claims (the order workbook) read
`report.screened` rather than assuming one happened.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .accessibility import (
    DEFAULT_ACCESSIBILITY_WEIGHT,
    DEFAULT_NUCLEATION_NT,
    AccessibilityResult,
    accessibility_penalty,
    not_evaluated,
    site_accessibility,
)
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
    scan_transcript,
)
from .loci import FusionJunction, Locus, RnaTarget
from .sequence import fetch_sequence, fetch_transcript, mask_repeats, rmsk_hits
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
    #: Target-site accessibility. Only meaningful for RNA: a DNA target is
    #: denatured by the protocol, so the question does not arise there, and a
    #: genomic design carries an explicit "not applicable" rather than a number.
    accessibility: AccessibilityResult = field(
        default_factory=lambda: not_evaluated(
            "DNA target: the protocol denatures it, so target-site secondary "
            "structure does not gate binding"))

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


#: Which reference an RNA design is screened against, and why.
#:
#: TRANSCRIPTOME ONLY, by default. Two reasons, and the second is the one that
#: decides it.
#:
#: 1. The transcriptome is the population of molecules actually present. An RNA
#:    probe competes with every other transcript in the cell, which is precisely
#:    what a transcriptome index represents and what a genome index does not: a
#:    gene's expressed isoforms are separate molecules at separate abundances,
#:    and a genome collapses them into one locus.
#: 2. A footprint may legitimately span an exon-exon junction -- often the point,
#:    since it is what makes a probe isoform-specific. Such a footprint does not
#:    map contiguously to the genome at all. Screened against a genome its two
#:    halves land in different exons, thousands of nt apart, so the real
#:    on-target site is invisible and the output is actively misleading.
#:
#: Genomic DNA is in the sample too, but a standard RNA-FISH protocol has no
#: denaturation step, so that DNA stays double-stranded and a single-stranded
#: probe cannot invade it. If your protocol *does* denature -- a combined
#: DNA/RNA experiment, for instance -- then the genome matters and you should
#: screen against both: pass `also_genome_backend`. Doing so is not the default
#: because for the usual protocol it adds a reference the probe never meets, and
#: junction-spanning designs would be flagged for a geometry that is correct.
RNA_SCREEN_POLICY = "transcriptome-only by default; see RNA_SCREEN_POLICY"


def design_rna(target: RnaTarget, backend=None, screen: bool = True,
               top: int = 5, transcript_seq: str | None = None,
               also_genome_backend=None, check_accessibility: bool = True,
               accessibility_top: int = 25,
               accessibility_weight: float = DEFAULT_ACCESSIBILITY_WEIGHT,
               nucleation_nt: int = DEFAULT_NUCLEATION_NT,
               **screen_kwargs) -> DesignResult:
    """Design a split-probe pair against a transcript.

    `backend` should index a TRANSCRIPTOME -- see `RNA_SCREEN_POLICY` for why,
    and for when to also pass `also_genome_backend`.

    Target-site accessibility is evaluated for each candidate and folded into
    the ranking -- see `accessibility`. It never rejects a window, and when
    RNAplfold is absent every candidate scores 0 for it, so the ordering is
    exactly what it was before accessibility existed.

    The screen is given `target_gene`, so the probe's own gene is excluded
    instead of being reported as an off-target in each of its isoforms. Note the
    co-localisation window is in transcript coordinates here: two half-sites 50
    nt apart in a spliced mRNA may be far apart in the genome, and the spliced
    distance is the one that matters, because the spliced molecule is what the
    probe binds.
    """
    system = get_system(target.system)
    seq = transcript_seq if transcript_seq is not None else fetch_transcript(
        target.accession)
    # Scan wider than requested so accessibility has candidates to choose
    # between: the most accessible window is often not the best-scoring one on
    # GC and Tm alone, and re-ranking a list of one achieves nothing.
    cands = scan_transcript(seq, target.accession, top=max(top, accessibility_top),
                            scan_start=target.scan_start, scan_end=target.scan_end,
                            exon_boundaries=target.exon_boundaries)
    if not cands:
        raise ValueError(f"no valid 51-nt window found in {target.accession}")

    # Accessibility is RNA-only and is applied HERE, not in `design._score`.
    # That keeps `_score` byte-comparable with every design made before this
    # existed -- adding a term there would retroactively reorder genomic
    # designs that have already been ordered. See accessibility.py on why any
    # threshold is uncalibrated; this is a ranking input, never a reject.
    access: dict[int, AccessibilityResult] = {}
    if check_accessibility:
        for cand in cands:
            access[cand.offset] = site_accessibility(
                seq, cand.offset, nucleation_nt=nucleation_nt)
        cands = sorted(
            cands,
            key=lambda c: c.penalty + accessibility_penalty(
                access[c.offset], weight=accessibility_weight),
        )
    cands = cands[:top]
    pair = build_pair(
        target.name, system, cands[0],
        p1l_name=f"P1L_{target.gene}_{target.name}",
        p1r_name=f"P1R_{target.gene}_{target.name}",
    )
    locus = Locus(name=target.name, gene=target.gene,
                  description=(target.description
                               or f"RNA target on {target.accession}"),
                  interval=cands[0].interval, system=system.name)
    if screen:
        spec = screen_pair(cands[0].left.seq, cands[0].right.seq, backend=backend,
                           target_gene=target.gene, **screen_kwargs)
        if also_genome_backend is not None:
            genome = screen_pair(cands[0].left.seq, cands[0].right.seq,
                                 backend=also_genome_backend, **screen_kwargs)
            # Kept as a note rather than merged: the two screens answer different
            # questions against different coordinate spaces, and averaging their
            # verdicts would obscure both.
            spec.reason = (spec.reason + " | additional genome screen: "
                           + genome.summary()).strip(" |")
    else:
        spec = not_screened("genome screen skipped by request (screen=False)")
    return DesignResult(locus, system, seq, cands, pair, rmsk=[], specificity=spec,
                        accessibility=access.get(
                            cands[0].offset,
                            not_evaluated("accessibility check not requested")))
