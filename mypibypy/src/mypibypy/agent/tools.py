"""The deterministic tool surface the loop drives.

THE INVARIANT: no tool here accepts a sequence. Every oligo comes out of
`design`/`pipeline`, which read sequence from the genome, the transcriptome or a
loci file -- never from a decision. The model's output schema has no sequence
field, these functions have no sequence parameter, and `TOOL_SIGNATURES` plus
`test_agent.py` assert both. That is what makes it structural rather than a
matter of remembering.

Tools take a `DesignParameters` and a target identifier, and return evidence.
They are thin: the real work is in the modules they call, and keeping them thin
is what stops a second, divergent design path growing inside the agent layer.
"""
from __future__ import annotations

import inspect
from dataclasses import dataclass, field

from .. import thermo
from ..design import FOOTPRINT, PadQC, select_tiled
from ..loci import FusionJunction, Locus, RnaTarget
from ..sequence import GenomicInterval
from .actions import DesignParameters


@dataclass
class WindowStats:
    """How many windows survive the current constraints, and what the near
    misses look like. The evidence behind a relax decision."""

    total: int = 0
    passing: int = 0
    failed_gc: int = 0
    failed_tm: int = 0
    best_penalty: float = float("nan")
    #: The tightest constraint, by how many windows it alone removes.
    limiting_constraint: str = ""

    @property
    def pass_fraction(self) -> float:
        return self.passing / self.total if self.total else 0.0

    def __str__(self) -> str:
        return (f"{self.passing}/{self.total} windows pass "
                f"({self.pass_fraction * 100:.1f}%); "
                f"{self.failed_gc} fail GC, {self.failed_tm} fail Tm; "
                f"limiting constraint: {self.limiting_constraint or 'none'}")


def _pad_ok(pad: PadQC, params: DesignParameters) -> tuple[bool, bool]:
    """(gc_ok, tm_ok) for one pad under the current parameters."""
    gc_ok = params.gc_min <= pad.gc <= params.gc_max
    tm_ok = params.tm_min <= pad.tm <= params.tm_max
    return gc_ok, tm_ok


def window_stats(candidates, params: DesignParameters) -> WindowStats:
    """Count how the current constraints bite across every scanned window.

    Counts each failure cause independently, so "how many does GC alone remove"
    is answerable -- which is the question a relax decision actually turns on.
    """
    stats = WindowStats(total=len(candidates))
    gc_only = tm_only = 0
    for cand in candidates:
        lg, lt = _pad_ok(cand.left, params)
        rg, rt = _pad_ok(cand.right, params)
        gc_ok, tm_ok = lg and rg, lt and rt
        if gc_ok and tm_ok:
            stats.passing += 1
        if not gc_ok:
            stats.failed_gc += 1
            if tm_ok:
                gc_only += 1
        if not tm_ok:
            stats.failed_tm += 1
            if gc_ok:
                tm_only += 1
    if candidates:
        stats.best_penalty = min(c.penalty for c in candidates)
    if gc_only or tm_only:
        stats.limiting_constraint = "GC" if gc_only >= tm_only else "Tm"
    elif stats.passing == 0 and stats.total:
        stats.limiting_constraint = "GC and Tm together"
    return stats


def filter_candidates(candidates, params: DesignParameters) -> list:
    """Windows meeting the current GC and Tm constraints, best first.

    With `enforce_window` off -- the default -- nothing is filtered and the full
    ranked list comes back, which is exactly what the deterministic pipeline
    does. `window_stats` still reports how the window *would* bite, so the
    evidence is available without the constraint being applied.
    """
    if not params.enforce_window:
        return list(candidates)
    kept = []
    for cand in candidates:
        lg, lt = _pad_ok(cand.left, params)
        rg, rt = _pad_ok(cand.right, params)
        if lg and rg and lt and rt:
            kept.append(cand)
    return kept


@dataclass
class ScreenEvidence:
    """What the screens said, reduced to what a decision turns on."""

    verdict: str = "not screened"
    screened: bool = False
    signal_capable_sites: int = 0
    offtarget_genes: list[str] = field(default_factory=list)
    single_arm_hits: int = 0
    paralog_like: bool = False       # one arm hits several loci, partner hits none
    accessibility_evaluated: bool = False
    accessibility_p: float | None = None
    limiting_pad: str = ""

    def __str__(self) -> str:
        bits = [f"specificity verdict {self.verdict}"]
        if self.screened:
            bits.append(f"{self.signal_capable_sites} signal-capable off-target site(s)")
            if self.offtarget_genes:
                bits.append(f"off-target genes {self.offtarget_genes[:5]}")
            bits.append(f"{self.single_arm_hits} lone-arm hit(s)")
            if self.paralog_like:
                bits.append("one arm decorates several loci while its partner has "
                            "nothing there (paralog-like)")
        if self.accessibility_evaluated:
            bits.append(f"accessibility P={self.accessibility_p:.3f} "
                        f"(limiting pad {self.limiting_pad})")
        return "; ".join(bits)


#: A lone arm hitting at least this many places is worth a human decision rather
#: than a silent pass: the pair still cannot make a spot there, but an arm that
#: decorates several paralogous loci may or may not be acceptable depending on
#: the experiment. Starting point, uncalibrated.
PARALOG_ARM_HIT_THRESHOLD = 3


def screen_evidence(specificity, accessibility=None) -> ScreenEvidence:
    """Reduce a `SpecificityReport` (+ accessibility) to decision evidence."""
    ev = ScreenEvidence(verdict=specificity.verdict, screened=specificity.screened)
    if specificity.screened:
        ev.signal_capable_sites = len(
            [s for s in specificity.pair_sites if s.signal_capable and s.survives_wash])
        ev.offtarget_genes = list(specificity.offtarget_genes)
        singles = specificity.single_half_hits
        ev.single_arm_hits = len(singles)
        # Paralog-like: many hits concentrated on ONE arm, with the partner
        # absent. That is the KAT6B shape -- geometry says harmless, the
        # experiment says maybe not.
        left = sum(1 for h in singles if h.query == "left")
        right = len(singles) - left
        ev.paralog_like = (max(left, right) >= PARALOG_ARM_HIT_THRESHOLD
                           and min(left, right) == 0)
    if accessibility is not None and accessibility.evaluated:
        ev.accessibility_evaluated = True
        ev.accessibility_p = accessibility.pair_nucleation_p
        ev.limiting_pad = accessibility.limiting_pad.pad
    return ev


def widened(interval: GenomicInterval, by_nt: int) -> GenomicInterval:
    """A genomic interval extended each side. Clamped at 1 on the left."""
    if by_nt <= 0:
        return interval
    return GenomicInterval(interval.build, interval.chrom,
                           max(1, interval.start - by_nt), interval.end + by_nt)


def capacity(region_length: int, params: DesignParameters) -> int:
    """How many tiled pairs the region can hold at the current spacing.

    Answers "is reduce_pairs even the problem here", which otherwise gets
    guessed at.
    """
    stride = FOOTPRINT + max(0, params.min_spacing)
    if region_length < FOOTPRINT:
        return 0
    return 1 + max(0, (region_length - FOOTPRINT) // stride)


def thermo_margin(pad_seq: str) -> float:
    """Worst wash margin for a pad, in C. Positive means it stays bound."""
    return thermo.profile(pad_seq, 0)["min_margin_c"]


#: Every callable the loop may invoke, by name. The test suite walks this and
#: asserts that not one of them takes a sequence, so a tool added later that
#: accepts model-supplied bases fails the build rather than shipping.
TOOL_SIGNATURES = {
    "window_stats": window_stats,
    "filter_candidates": filter_candidates,
    "screen_evidence": screen_evidence,
    "widened": widened,
    "capacity": capacity,
    "select_tiled": select_tiled,
}

#: Parameter names that would let a sequence in. Checked against every tool.
FORBIDDEN_TOOL_PARAMETERS = frozenset({
    "seq", "sequence", "region_seq", "transcript_seq", "footprint", "pad",
    "left_seq", "right_seq", "p1l", "p1r", "probe", "oligo", "bases", "dna",
})


def sequence_accepting_tools() -> dict[str, list[str]]:
    """Any tool that takes something a sequence could be passed as.

    Used by the test that enforces the invariant. `thermo_margin` is excluded
    from `TOOL_SIGNATURES` deliberately -- it does take a sequence, and it is
    called by the harness with a pad the deterministic code produced, never with
    anything a decision supplied.
    """
    offenders: dict[str, list[str]] = {}
    for name, fn in TOOL_SIGNATURES.items():
        params = inspect.signature(fn).parameters
        bad = [p for p in params if p.lower() in FORBIDDEN_TOOL_PARAMETERS]
        if bad:
            offenders[name] = bad
    return offenders


def target_label(target) -> str:
    """A short name for a locus, fusion junction or RNA target."""
    if isinstance(target, Locus):
        return f"{target.gene} {target.name} ({target.interval})"
    if isinstance(target, FusionJunction):
        return f"{target.five_prime_gene}::{target.three_prime_gene} ({target.name})"
    if isinstance(target, RnaTarget):
        return f"{target.gene} {target.name} ({target.accession})"
    return str(target)
