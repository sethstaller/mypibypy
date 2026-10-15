"""Split-probe (P1L/P1R) design over a target interval.

A target footprint is 51 nt on the plus strand: a 25-nt left pad, a 1-nt gap, and
a 25-nt right pad. The probes hybridize to the plus strand:

    P1L = arm_left  + revcomp(left25)
    P1R = revcomp(right25) + arm_right
"""
from __future__ import annotations

from dataclasses import dataclass, field

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
