"""Buffer-aware thermodynamics.

Checks internal consistency and physical direction, not absolute accuracy: the
formamide term is a first-order correction and the mismatch penalty is an
explicit approximation, so these tests pin down behaviour and sign, not values
to the degree.
"""
from __future__ import annotations

import math

import pytest

from mypibypy import thermo
from mypibypy.sequence import nn_tm
from mypibypy.thermo import (
    ALL_STEPS,
    HYB_PRIMARY,
    SSC_NA,
    WASH_A,
    WASH_B,
    WASH_C,
    Conditions,
    margin,
    most_stringent,
    profile,
    survives,
    tm_buffer,
    tm_mismatched,
)

# Neutral synthetic 25-nt pads, matched in length and GC to a real target pad.
# No real design sequences are committed to this repository.
PAD_A = "ACTGACAGTCATGCAGTACGTACAG"
PAD_B = "GTCAGTCATGACGTCAGTACTGACT"


class TestConditions:
    def test_ssc_to_sodium(self):
        assert HYB_PRIMARY.na_molar == pytest.approx(6 * SSC_NA)
        assert WASH_C.na_molar == pytest.approx(0.1 * SSC_NA)

    def test_published_protocol_values(self):
        assert (HYB_PRIMARY.temp_c, HYB_PRIMARY.formamide_pct, HYB_PRIMARY.ssc_x) == (40.0, 25.0, 6.0)
        assert [c.formamide_pct for c in (WASH_A, WASH_B, WASH_C)] == [20.0, 10.0, 0.0]


class TestTmBuffer:
    def test_formamide_lowers_tm_linearly(self):
        base = Conditions("x", 40, 0, 6)
        with_fa = Conditions("x", 40, 25, 6)
        drop = tm_buffer(PAD_A, base) - tm_buffer(PAD_A, with_fa)
        assert drop == pytest.approx(25 * thermo.FORMAMIDE_C_PER_PCT)

    def test_more_salt_raises_tm(self):
        low = Conditions("low", 40, 0, 0.1)
        high = Conditions("high", 40, 0, 6)
        assert tm_buffer(PAD_A, high) > tm_buffer(PAD_A, low)

    def test_gc_rich_melts_higher(self):
        gc = "GCGCGCGCGCGCGCGCGCGCGCGCG"
        at = "ATATATATATATATATATATATATA"
        assert tm_buffer(gc, HYB_PRIMARY) > tm_buffer(at, HYB_PRIMARY)

    def test_invalid_sequence_is_nan(self):
        assert math.isnan(tm_buffer("ACGTN", HYB_PRIMARY))
        assert math.isnan(tm_buffer("A", HYB_PRIMARY))

    def test_agrees_with_reference_tm_when_conditions_match(self):
        """Same physics as sequence.nn_tm, modulo salt/formamide/concentration."""
        same = Conditions("ref", 0, 0, 0.33 / SSC_NA, probe_nM=250.0)
        assert tm_buffer(PAD_A, same) == pytest.approx(nn_tm(PAD_A), abs=0.2)

    def test_reference_tm_is_not_a_buffer_tm(self):
        """The whole reason this module exists: the ranking Tm overstates things.

        For a typical 25-nt pad the gap is >10 C: the reference Tm reads near
        68 C where the primary hybridization sits nearer 55 C.
        """
        assert nn_tm(PAD_A) - tm_buffer(PAD_A, HYB_PRIMARY) > 10


class TestMargins:
    def test_typical_pads_survive_every_step(self):
        for pad in (PAD_A, PAD_B):
            assert survives(pad), profile(pad)

    def test_margin_is_tm_minus_step_temperature(self):
        assert margin(PAD_A, WASH_A) == pytest.approx(
            tm_buffer(PAD_A, WASH_A) - WASH_A.temp_c)

    def test_limiting_step_is_identified(self):
        worst = most_stringent(PAD_A)
        assert worst in ALL_STEPS
        assert margin(PAD_A, worst) == min(margin(PAD_A, c) for c in ALL_STEPS)

    def test_short_probe_fails_the_wash(self):
        assert not survives("ACGTACGTAC")


class TestMismatches:
    def test_mismatches_lower_tm(self):
        t0 = tm_mismatched(PAD_A, 0, HYB_PRIMARY)
        t1 = tm_mismatched(PAD_A, 1, HYB_PRIMARY)
        t2 = tm_mismatched(PAD_A, 2, HYB_PRIMARY)
        assert t0 > t1 > t2
        assert t0 - t1 == pytest.approx(thermo.DEFAULT_MISMATCH_PENALTY_C)

    def test_zero_mismatches_equals_perfect_duplex(self):
        assert tm_mismatched(PAD_A, 0, WASH_B) == tm_buffer(PAD_A, WASH_B)

    def test_negative_mismatch_count_rejected(self):
        with pytest.raises(ValueError):
            tm_mismatched(PAD_A, -1, HYB_PRIMARY)

    def test_enough_mismatches_stop_survival(self):
        """The point of the module: an aligner hit is not a surviving duplex."""
        assert survives(PAD_A, n_mismatches=0)
        assert not survives(PAD_A, n_mismatches=4)

    def test_penalty_is_configurable(self):
        lenient = tm_mismatched(PAD_A, 2, HYB_PRIMARY, penalty_c=4.0)
        strict = tm_mismatched(PAD_A, 2, HYB_PRIMARY, penalty_c=12.0)
        assert lenient > strict


class TestProfile:
    def test_profile_reports_every_step(self):
        p = profile(PAD_A)
        assert len(p["steps"]) == len(ALL_STEPS)
        assert p["limiting_step"] in {c.name for c in ALL_STEPS}
        assert p["survives"] is True
        assert p["approximate"] is False

    def test_profile_marks_mismatched_results_approximate(self):
        assert profile(PAD_A, n_mismatches=2)["approximate"] is True

    def test_profile_min_margin_matches_steps(self):
        p = profile(PAD_B)
        assert p["min_margin_c"] == pytest.approx(min(s["margin_c"] for s in p["steps"]))
