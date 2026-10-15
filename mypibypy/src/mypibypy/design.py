"""Split-probe (P1L/P1R) design over a target interval, or across a fusion junction.

A target footprint is 51 nt on the plus strand: a 25-nt left pad, a 1-nt gap, and
a 25-nt right pad. The probes hybridize to the plus strand:

    P1L = arm_left  + revcomp(left25)
    P1R = revcomp(right25) + arm_right

FUSION MODE puts that 1-nt gap exactly on the breakpoint, which turns the
coincidence requirement into a junction test: the left pad lies wholly in the 5'
partner, the right pad wholly in the 3' partner, and neither parent gene alone
can bring both halves together. The skipped gap base is the one position that
would otherwise have to be assigned to one partner or the other.

Fusion targets are screened on a rule that normal targets have no use for --
each pad must be found in its own partner and NOT in the other
(`check_partner_specificity`) -- and they are designed as a sense/antisense
quartet, because the lab orders both strands for a junction.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

from .sequence import (
    GenomicInterval,
    gc_fraction,
    max_dinucleotide_repeat,
    max_homopolymer,
    nn_tm,
    revcomp,
)
from .systems import AmplifierSystem

LEFT_LEN = 25
GAP_LEN = 1
RIGHT_LEN = 25
FOOTPRINT = LEFT_LEN + GAP_LEN + RIGHT_LEN  # 51


@dataclass
class PadQC:
    seq: str
    gc: float
    tm: float
    max_homopolymer: int
    max_dinuc: int


@dataclass
class Candidate:
    offset: int                 # 0-based offset within the region sequence
    interval: GenomicInterval   # genomic coords of the 51-nt footprint
    left: PadQC
    right: PadQC
    gap_base: str
    penalty: float
    flags: list[str] = field(default_factory=list)

    @property
    def gc_mean(self) -> float:
        return (self.left.gc + self.right.gc) / 2

    @property
    def tm_delta(self) -> float:
        return abs(self.left.tm - self.right.tm)


@dataclass
class ProbePair:
    locus: str
    system: AmplifierSystem
    candidate: Candidate
    p1l_name: str
    p1r_name: str
    p1l: str
    p1r: str


def _pad_qc(seq: str) -> PadQC:
    return PadQC(
        seq=seq,
        gc=gc_fraction(seq),
        tm=nn_tm(seq),
        max_homopolymer=max_homopolymer(seq),
        max_dinuc=max_dinucleotide_repeat(seq),
    )


def _score(left: PadQC, right: PadQC) -> tuple[float, list[str]]:
    """Lower penalty is better. Returns (penalty, flags)."""
    flags: list[str] = []
    penalty = 0.0
    # distance of each pad's GC from 0.5
    penalty += abs(left.gc - 0.5) + abs(right.gc - 0.5)
    # mean GC away from ideal band, extra weight
    penalty += abs((left.gc + right.gc) / 2 - 0.5)
    # imbalance between the two pads (want matched Tm so they bind together)
    penalty += abs(left.tm - right.tm) / 20.0
    for tag, pad in (("L", left), ("R", right)):
        if pad.max_homopolymer >= 5:
            penalty += 1.0
            flags.append(f"{tag}:homopolymer{pad.max_homopolymer}")
        if pad.max_dinuc >= 8:
            penalty += 1.0
            flags.append(f"{tag}:dinuc{pad.max_dinuc}")
        if not 0.30 <= pad.gc <= 0.70:
            flags.append(f"{tag}:gc{pad.gc*100:.0f}%")
    return penalty, flags


def scan(region_seq: str, region: GenomicInterval, top: int = 5) -> list[Candidate]:
    """Score every 51-nt window in the region; return the best `top`, ranked."""
    if region.length != len(region_seq):
        raise ValueError("region interval length != sequence length")
    cands: list[Candidate] = []
    for i in range(len(region_seq) - FOOTPRINT + 1):
        window = region_seq[i : i + FOOTPRINT]
        if any(b not in "ACGT" for b in window):
            continue
        left = _pad_qc(window[:LEFT_LEN])
        right = _pad_qc(window[LEFT_LEN + GAP_LEN :])
        penalty, flags = _score(left, right)
        fp = GenomicInterval(
            region.build, region.chrom,
            region.start + i, region.start + i + FOOTPRINT - 1,
        )
        cands.append(
            Candidate(
                offset=i, interval=fp, left=left, right=right,
                gap_base=window[LEFT_LEN], penalty=penalty, flags=flags,
            )
        )
    cands.sort(key=lambda c: c.penalty)
    return cands[:top]


def build_pair(
    locus: str, system: AmplifierSystem, cand: Candidate,
    p1l_name: str, p1r_name: str,
) -> ProbePair:
    p1l = system.arm_left + revcomp(cand.left.seq)
    p1r = revcomp(cand.right.seq) + system.arm_right
    return ProbePair(
        locus=locus, system=system, candidate=cand,
        p1l_name=p1l_name, p1r_name=p1r_name, p1l=p1l, p1r=p1r,
    )


# --- Fusion junctions ------------------------------------------------------

#: Which partner gives up the skipped base. The gap sits *on* the breakpoint
#: either way; this says whether the sacrificed base is the last of the 5'
#: partner or the first of the 3' partner. There is no biophysical reason to
#: prefer one, so it is explicit rather than implied -- but it must be recorded,
#: since it shifts the whole footprint by one base.
GAP_FIVE_PRIME = "five_prime"
GAP_THREE_PRIME = "three_prime"


class FusionGeometryError(ValueError):
    """The junction cannot carry a 51-nt footprint with the gap on the breakpoint."""


def fusion_footprint(sequence: str, breakpoint: int,
                     gap_side: str = GAP_FIVE_PRIME) -> tuple[int, str]:
    """The 51-nt footprint centred on a breakpoint.

    `breakpoint` is the 0-based index of the first base contributed by the 3'
    partner, i.e. ``sequence[:breakpoint]`` is 5'-partner sequence and
    ``sequence[breakpoint:]`` is 3'-partner sequence.

    Returns ``(offset, footprint)`` where `offset` indexes into `sequence`.
    """
    if gap_side == GAP_FIVE_PRIME:
        # Skip the last base of the 5' partner.
        offset = breakpoint - 1 - LEFT_LEN
    elif gap_side == GAP_THREE_PRIME:
        # Skip the first base of the 3' partner.
        offset = breakpoint - LEFT_LEN
    else:
        raise ValueError(f"gap_side must be {GAP_FIVE_PRIME!r} or {GAP_THREE_PRIME!r}")
    if offset < 0 or offset + FOOTPRINT > len(sequence):
        need_left = LEFT_LEN + (1 if gap_side == GAP_FIVE_PRIME else 0)
        need_right = RIGHT_LEN + (1 if gap_side == GAP_THREE_PRIME else 0)
        raise FusionGeometryError(
            f"need {need_left} nt of 5' partner and {need_right} nt of 3' partner "
            f"around the breakpoint; breakpoint {breakpoint} in a {len(sequence)} nt "
            f"sequence leaves {breakpoint} / {len(sequence) - breakpoint}"
        )
    return offset, sequence[offset : offset + FOOTPRINT]


def split_footprint(footprint: str) -> tuple[str, str, str]:
    """(left25, gap_base, right25). The one place the 51-nt layout is written down."""
    if len(footprint) != FOOTPRINT:
        raise FusionGeometryError(
            f"footprint must be {FOOTPRINT} nt (25 + 1 + 25), got {len(footprint)}")
    return (footprint[:LEFT_LEN], footprint[LEFT_LEN],
            footprint[LEFT_LEN + GAP_LEN :])


def check_partner_specificity(footprint: str, five_prime_seq: str,
                              three_prime_seq: str) -> list[str]:
    """Each pad must come from its own partner and be absent from the other.

    This is the check that makes a junction probe a junction probe. If the left
    pad can also be found in the 3' partner, then that partner alone presents
    both halves and the pair reports the unfused gene -- which is the entire
    thing a fusion probe exists to avoid.

    Returns a list of problems; empty means the junction test holds. Searches
    both strands, because a pad matching the reverse complement of the wrong
    partner is just as capable of recruiting the probe.
    """
    left, _, right = split_footprint(footprint)
    five, three = five_prime_seq.upper(), three_prime_seq.upper()
    problems: list[str] = []

    def present(pad: str, hay: str) -> bool:
        return pad in hay or revcomp(pad) in hay

    if not present(left, five):
        problems.append("left pad is not found in the 5' partner sequence")
    if not present(right, three):
        problems.append("right pad is not found in the 3' partner sequence")
    if present(left, three):
        problems.append(
            "left pad also occurs in the 3' partner -- that partner alone would "
            "present both halves, so the pair is not junction-specific")
    if present(right, five):
        problems.append(
            "right pad also occurs in the 5' partner -- that partner alone would "
            "present both halves, so the pair is not junction-specific")
    return problems


def fusion_candidate(junction_name: str, sequence: str, breakpoint: int,
                     gap_side: str = GAP_FIVE_PRIME,
                     antisense: bool = False) -> Candidate:
    """A `Candidate` for a junction footprint, on either strand.

    Coordinates are offsets into the supplied junction sequence, not genomic
    ones: a fusion has no single place in the reference. The build is recorded
    as the junction name so anything that prints the interval still says what it
    is looking at.
    """
    offset, footprint = fusion_footprint(sequence, breakpoint, gap_side)
    if antisense:
        footprint = revcomp(footprint)
    left, gap, right = split_footprint(footprint)
    l_qc, r_qc = _pad_qc(left), _pad_qc(right)
    penalty, flags = _score(l_qc, r_qc)
    strand = "antisense" if antisense else "sense"
    iv = GenomicInterval(f"{junction_name}:{strand}", "junction",
                         offset + 1, offset + FOOTPRINT)
    return Candidate(offset=offset, interval=iv, left=l_qc, right=r_qc,
                     gap_base=gap, penalty=penalty, flags=flags)


#: How the antisense pair's two pads are assigned to the two arms.
#:
#: ``juxtaposed`` applies the unchanged P1L/P1R rule to revcomp(footprint). Both
#: arms then meet at the gap, exactly as in the sense pair, and P2's landing site
#: can form.
#:
#: ``as_ordered`` reproduces the convention used in the lab's two hand-made
#: fusion designs: the sense footprint's left25 used *forward* rather than
#: reverse complemented, keeping the left25 -> P1L mapping. Because the antisense
#: probes bind the opposite strand, that points both arms AWAY from the junction
#: -- simulated at 50 nt apart on a 51-nt footprint, versus 2 nt for a working
#: pair. P2's site is never reconstituted, so such a pair should give no signal.
#: It exists here only so the historical designs can be regenerated and audited.
#: See `verify_against`.
PADS_JUXTAPOSED = "juxtaposed"
PADS_AS_ORDERED = "as_ordered"


def build_fusion_pairs(junction_name: str, sequence: str, breakpoint: int,
                       sense_system: AmplifierSystem,
                       antisense_system: AmplifierSystem | None = None,
                       gap_side: str = GAP_FIVE_PRIME,
                       name_prefix: str | None = None,
                       pad_assignment: str = PADS_JUXTAPOSED) -> list[ProbePair]:
    """The sense and antisense pairs for one junction -- the lab's unit of order.

    A junction is always ordered as both strands. `antisense_system` defaults to
    `sense_system`, but the lab's own fusion designs put the two strands on
    *different* channels so the strands can be told apart in one image; pass the
    second system to reproduce that.

    Note the identity worth knowing: the antisense footprint is the reverse
    complement of the sense one, so the antisense probes carry the *sense*
    footprint's bases directly -- P1L(antisense) pad is right25 read forward,
    and P1R(antisense) pad is left25. The two pairs therefore have perfectly
    complementary pads and will dimerize if pooled in one tube; `crosshyb`
    classifies that as expected-by-design rather than flagging it.
    """
    if pad_assignment not in (PADS_JUXTAPOSED, PADS_AS_ORDERED):
        raise ValueError(f"pad_assignment must be {PADS_JUXTAPOSED!r} or "
                         f"{PADS_AS_ORDERED!r}")
    antisense_system = antisense_system or sense_system
    prefix = name_prefix or junction_name
    pairs = []
    for anti, system in ((False, sense_system), (True, antisense_system)):
        cand = fusion_candidate(junction_name, sequence, breakpoint, gap_side,
                                antisense=anti)
        tag = "antisense" if anti else "sense"
        if anti and pad_assignment == PADS_AS_ORDERED:
            # Swap which pad gets which arm, reproducing the hand-made designs.
            cand = replace(cand, left=cand.right, right=cand.left,
                           flags=[*cand.flags, "pads swapped: arms point away "
                                  "from the junction, P2 cannot bind"])
        pairs.append(build_pair(
            f"{junction_name}_{tag}", system, cand,
            p1l_name=f"P1L_{prefix}_{tag}", p1r_name=f"P1R_{prefix}_{tag}",
        ))
    return pairs


def pairs_from_footprint(name: str, footprint: str, sense_system: AmplifierSystem,
                         antisense_system: AmplifierSystem | None = None,
                         name_prefix: str | None = None,
                         pad_assignment: str = PADS_JUXTAPOSED) -> list[ProbePair]:
    """Sense + antisense pairs from a 51-nt footprint that is already chosen.

    The breakpoint-placement step is skipped: this takes the footprint as given.
    It is how a hand-made design is reproduced, and how one is checked.
    """
    if len(footprint) != FOOTPRINT:
        raise FusionGeometryError(
            f"footprint must be {FOOTPRINT} nt, got {len(footprint)}")
    # A footprint supplied directly is treated as its own sequence, with the
    # breakpoint at the gap -- which reduces to the same geometry.
    return build_fusion_pairs(name, footprint, LEFT_LEN + 1, sense_system,
                              antisense_system, gap_side=GAP_FIVE_PRIME,
                              name_prefix=name_prefix, pad_assignment=pad_assignment)


def arm_offsets(pair: ProbePair, footprint: str) -> tuple[int, int] | None:
    """Where each arm attaches, as an offset into the 51-nt footprint.

    This is the geometry check behind `PADS_AS_ORDERED`: a probe's pad anneals
    somewhere on the footprint, and the arm extends outward from whichever end of
    the pad it is attached to. For a working pair both arms land next to the gap.
    Returns None if either pad does not anneal to this footprint at all.
    """
    system = pair.system
    out = []
    for probe, arm, arm_at_5prime in (
            (pair.p1l, system.arm_left, True), (pair.p1r, system.arm_right, False)):
        pad = probe[len(arm):] if arm_at_5prime else probe[: -len(arm)]
        placed = None
        for a in range(len(footprint) - len(pad) + 1):
            window = footprint[a : a + len(pad)]
            if pad == revcomp(window):
                # Binds the plus strand: pad reads 3'->5' left to right.
                placed = (a + len(pad) - 1) if arm_at_5prime else a
                break
            if pad == window:
                # Binds the minus strand: pad reads 5'->3' left to right.
                placed = a if arm_at_5prime else (a + len(pad) - 1)
                break
        if placed is None:
            return None
        out.append(placed)
    return out[0], out[1]


def arms_are_juxtaposed(pair: ProbePair, footprint: str, slack: int = 1) -> bool:
    """Do both arms reach the gap, so P2's landing site can form?

    The whole coincidence mechanism depends on this and nothing in the probe
    sequences alone reveals it -- a pair with its pads on the wrong arms looks
    perfectly well-formed until you ask where the arms point.
    """
    offsets = arm_offsets(pair, footprint)
    if offsets is None:
        return False
    gap = LEFT_LEN  # 0-based index of the skipped base
    return all(abs(o - gap) <= slack + 1 for o in offsets)


def verify_against(built: list[ProbePair], ordered: dict[str, str]) -> list[dict]:
    """Compare freshly built primaries with what was actually ordered.

    `ordered` maps a label to a sequence; the labels are matched positionally
    against `built` in (sense P1L, sense P1R, antisense P1L, antisense P1R)
    order. Returns one row per primary describing agreement and, where they
    differ, by how much -- so a hand-made design can be audited without the
    sequences going anywhere near version control.
    """
    flat = []
    for pair in built:
        flat.append((pair.p1l_name, pair.p1l))
        flat.append((pair.p1r_name, pair.p1r))
    rows = []
    for (name, seq), (label, want) in zip(flat, ordered.items()):
        same = seq == want
        row = {"built_name": name, "ordered_name": label, "match": same,
               "built_len": len(seq), "ordered_len": len(want)}
        if not same:
            if len(seq) == len(want):
                row["note"] = (f"same length, "
                               f"{sum(a != b for a, b in zip(seq, want))} positions differ")
            else:
                row["note"] = f"length differs by {len(seq) - len(want):+d} nt"
        rows.append(row)
    return rows
