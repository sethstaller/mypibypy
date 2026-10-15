"""Buffer-aware duplex thermodynamics.

`sequence.nn_tm` is a *fixed-reference* Tm (0.33 M Na+, no formamide) used to rank
candidate windows against each other; it is deliberately left alone so rankings
stay comparable with earlier designs. This module answers a different question:

    at the bench, in this buffer, at this temperature -- does the duplex hold?

Everything here is reported as a **margin**: Tm minus the step temperature.
A positive margin means the duplex is expected to survive that step; negative
means it should come apart. That framing matters for off-target screening, where
"the aligner found it" is not the same as "it survives the wash".

Formamide correction: linear, 0.65 C per % (v/v). The widely used range is
0.6-0.75 C/% for DNA duplexes (McConaughy et al. 1969, Blake & Delcourt 1996);
0.65 is mid-range. It is a first-order correction, not an exact model.

Salt: SantaLucia's log[Na+] entropy correction. 1x SSC is taken as 0.195 M Na+
(0.15 M NaCl contributes 0.15; 0.015 M trisodium citrate contributes 0.045).

MISMATCH HANDLING IS AN APPROXIMATION -- see `tm_mismatched`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

from .sequence import _NN_DH, _NN_DS, _R, gc_fraction

#: Na+ molarity contributed by 1x SSC.
SSC_NA = 0.195

#: Formamide depression, degrees C per % (v/v).
FORMAMIDE_C_PER_PCT = 0.65

#: Default conservative penalty per mismatch, in degrees C. See `tm_mismatched`.
DEFAULT_MISMATCH_PENALTY_C = 8.0


@dataclass(frozen=True)
class Conditions:
    """One hybridization or wash step."""

    name: str
    temp_c: float
    formamide_pct: float
    ssc_x: float
    probe_nM: float = 10.0

    @property
    def na_molar(self) -> float:
        return self.ssc_x * SSC_NA

    def __str__(self) -> str:
        return (f"{self.name} ({self.temp_c:g} C, {self.formamide_pct:g}% formamide, "
                f"{self.ssc_x:g}x SSC)")


# --- The lab's piFISH conditions -------------------------------------------
# Primary (P1) hybridization. Matches the published pi-FISH rainbow protocol
# (Tao et al., Nat Commun 2023): solution A, 25% formamide / 6x SSC, 40 C.
HYB_PRIMARY = Conditions("primary hyb", temp_c=40.0, formamide_pct=25.0, ssc_x=6.0, probe_nM=10.0)

# Wash series, also from the published protocol (elution buffers A/B/C, 40 C).
# NOTE: transcribed from the paper, NOT confirmed against this lab's bench
# protocol. Override with your own if they differ.
WASH_A = Conditions("wash A", temp_c=40.0, formamide_pct=20.0, ssc_x=2.0, probe_nM=10.0)
WASH_B = Conditions("wash B", temp_c=40.0, formamide_pct=10.0, ssc_x=1.0, probe_nM=10.0)
WASH_C = Conditions("wash C", temp_c=40.0, formamide_pct=0.0, ssc_x=0.1, probe_nM=10.0)

WASH_SERIES = (WASH_A, WASH_B, WASH_C)

#: The step a probe must survive to stay bound: the most stringent wash.
#: Stringency here is dominated by the formamide/salt combination, so the
#: limiting step is computed rather than assumed -- see `most_stringent`.
ALL_STEPS = (HYB_PRIMARY, *WASH_SERIES)


def nn_thermo(seq: str) -> tuple[float, float]:
    """(dH kcal/mol, dS cal/mol/K) at 1 M Na+ for a perfect DNA/DNA duplex."""
    s = seq.upper()
    if len(s) < 2 or any(b not in "ACGT" for b in s):
        return float("nan"), float("nan")
    dh = ds = 0.0
    for end in (s[0], s[-1]):
        if end in "GC":
            dh += 0.1
            ds += -2.8
        else:
            dh += 2.3
            ds += 4.1
    for i in range(len(s) - 1):
        dh += _NN_DH[s[i : i + 2]]
        ds += _NN_DS[s[i : i + 2]]
    return dh, ds


def tm_buffer(seq: str, cond: Conditions) -> float:
    """Melting temperature (C) of a perfect duplex under `cond`.

    Probe in excess over target (the FISH case), so Ct/4 is used as in
    `sequence.nn_tm`.
    """
    dh, ds = nn_thermo(seq)
    if math.isnan(dh):
        return float("nan")
    ds += 0.368 * (len(seq) - 1) * math.log(cond.na_molar)
    tm_k = (dh * 1000.0) / (ds + _R * math.log(cond.probe_nM * 1e-9 / 4.0))
    return tm_k - 273.15 - FORMAMIDE_C_PER_PCT * cond.formamide_pct


def margin(seq: str, cond: Conditions) -> float:
    """Tm - step temperature. Positive = expected to stay hybridized."""
    return tm_buffer(seq, cond) - cond.temp_c


def tm_mismatched(seq: str, n_mismatches: int, cond: Conditions,
                  penalty_c: float = DEFAULT_MISMATCH_PENALTY_C) -> float:
    """Approximate Tm of a duplex carrying `n_mismatches`.

    APPROXIMATION. A proper treatment needs the mismatch identities, their
    positions and internal-loop NN parameters, which a seed-and-extend aligner
    does not give us. Published single-internal-mismatch penalties for 20-25mers
    span roughly 4-12 C depending on identity and position; `penalty_c` defaults
    to 8 C, mid-range and deliberately on the destabilizing side so that
    off-targets are not dismissed too readily.

    Use this to *rank* and to *triage* off-target sites, never as a measured Tm.
    Terminal mismatches are far less destabilizing than internal ones, so a site
    whose mismatches sit at the very ends will be under-estimated here.
    """
    if n_mismatches < 0:
        raise ValueError("n_mismatches must be >= 0")
    return tm_buffer(seq, cond) - penalty_c * n_mismatches


def most_stringent(seq: str, steps: tuple[Conditions, ...] = ALL_STEPS) -> Conditions:
    """The step in `steps` giving this sequence its smallest margin."""
    return min(steps, key=lambda c: margin(seq, c))


def survives(seq: str, n_mismatches: int = 0, steps: tuple[Conditions, ...] = ALL_STEPS,
             penalty_c: float = DEFAULT_MISMATCH_PENALTY_C) -> bool:
    """Would a duplex with this many mismatches still be bound after every step?"""
    return all(
        tm_mismatched(seq, n_mismatches, c, penalty_c) > c.temp_c for c in steps
    )


def profile(seq: str, n_mismatches: int = 0,
            steps: tuple[Conditions, ...] = ALL_STEPS,
            penalty_c: float = DEFAULT_MISMATCH_PENALTY_C) -> dict:
    """Per-step Tm and margin, plus the limiting step. For reports and QC sheets."""
    rows = []
    for c in steps:
        tm = tm_mismatched(seq, n_mismatches, c, penalty_c)
        rows.append({"step": c.name, "conditions": str(c), "tm_c": tm,
                     "temp_c": c.temp_c, "margin_c": tm - c.temp_c})
    worst = min(rows, key=lambda r: r["margin_c"]) if rows else None
    return {
        "seq": seq,
        "length": len(seq),
        "gc": gc_fraction(seq),
        "n_mismatches": n_mismatches,
        "steps": rows,
        "limiting_step": worst["step"] if worst else None,
        "min_margin_c": worst["margin_c"] if worst else float("nan"),
        "survives": all(r["margin_c"] > 0 for r in rows),
        "mismatch_penalty_c": penalty_c,
        "approximate": n_mismatches > 0,
    }


def at_temperature(cond: Conditions, temp_c: float) -> Conditions:
    """Same buffer, different temperature."""
    return replace(cond, temp_c=temp_c)
