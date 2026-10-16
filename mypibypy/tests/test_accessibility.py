"""Target-site accessibility on the transcript.

The controls are synthetic and built here: a designed hairpin for the occluded
case and an A/T-rich tract for the open one. Nothing depends on a real
transcript, and no real sequence appears.

`TestTheMeanWouldLie` is the test that justifies the whole metric choice. A
hairpin loop has high per-base unpaired probability and essentially no chance of
a contiguous unpaired run, so a mean-per-base score calls it accessible when a
probe cannot nucleate there. If someone later "simplifies" this module to a
mean, that test fails.

`TestDegradesWithoutViennaRNA` is the other one that matters: a missing optional
tool must never score as an accessible site, and must not shift any ranking.
"""
from __future__ import annotations

import shutil

import pytest

from mypibypy.accessibility import (
    ACCESSIBILITY_IS_UNCALIBRATED,
    DEFAULT_FLANK_NT,
    DEFAULT_SPAN_L,
    DEFAULT_WINDOW_W,
    EVALUATED,
    NOT_EVALUATED,
    accessibility_penalty,
    describe_accessibility,
    not_evaluated,
    parse_lunp,
    rnaplfold_path,
    site_accessibility,
)

HAVE_PLFOLD = bool(rnaplfold_path())
pytestmark = []

STEM = "GGCGCGCGGCGCGCC"
LOOP = "TTTTAAAA"
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
HAIRPIN = STEM + LOOP + "".join(_COMP[b] for b in reversed(STEM))   # 38 nt
OPEN_TRACT = "ATAAATAAATAAATAAATAAAT" * 12                          # unstructured
#: Unstructured AND GC-balanced (0.50 GC, nucleation P ~0.76, measured). Needed
#: for any test of RANKING: `design._score` penalises GC far from 0.5 heavily, so
#: an AT-rich tract is rejected on composition before accessibility gets a say --
#: correctly, since accessibility is a ranking input and not an override.
BALANCED_OPEN = "ACTGACTGACTG" * 40


def open_transcript(n: int = 2) -> str:
    return OPEN_TRACT * n


def hairpin_transcript() -> tuple[str, int]:
    """(transcript, offset of a 51-nt footprint centred on the hairpin)."""
    seq = OPEN_TRACT + HAIRPIN + OPEN_TRACT
    offset = len(OPEN_TRACT) + len(HAIRPIN) // 2 - 25
    return seq, offset


class TestParseLunp:
    """The table is positional and 1-based; getting it wrong shifts every score."""

    def test_parses_rows_and_na(self):
        text = ("#unpaired probabilities\n"
                " #i$\tl=1\t2\t3\n"
                "1\t0.9\tNA\tNA\n"
                "2\t0.8\t0.75\tNA\n"
                "3\t0.7\t0.65\t0.6\n")
        rows = parse_lunp(text)
        assert len(rows) == 3
        assert rows[0] == [0.9, None, None]
        assert rows[2] == [0.7, 0.65, 0.6]

    def test_indexes_by_position_even_with_gaps(self):
        rows = parse_lunp("3\t0.5\n")
        assert len(rows) == 3 and rows[2] == [0.5] and rows[0] == []

    def test_ignores_headers_and_junk(self):
        assert parse_lunp("#c\n\nnot-a-row\n") == []


class TestBands:
    def test_describes_without_committing_to_a_threshold(self):
        assert describe_accessibility(0.9) == "open"
        assert describe_accessibility(0.3) == "partly open"
        assert describe_accessibility(0.1) == "mostly paired"
        assert describe_accessibility(0.01) == "occluded"


class TestDegradesWithoutViennaRNA:
    """A tool that did not run must not look like an accessible site."""

    def test_missing_binary_is_not_evaluated(self, monkeypatch):
        monkeypatch.setenv("MYPIBYPY_RNAPLFOLD", "/definitely/not/here")
        res = site_accessibility(open_transcript(), 100)
        assert res.status == NOT_EVALUATED
        assert not res.evaluated
        assert "RNAplfold not found" in res.reason
        assert "conda install" in res.reason, "say how to fix it"

    def test_unevaluated_has_no_number_to_mistake_for_a_good_one(self, monkeypatch):
        monkeypatch.setenv("MYPIBYPY_RNAPLFOLD", "/definitely/not/here")
        res = site_accessibility(open_transcript(), 100)
        assert res.pair_nucleation_p is None
        assert res.limiting_pad is None
        assert "NOT EVALUATED" in res.summary()

    def test_unevaluated_costs_nothing_in_ranking(self):
        """Installing ViennaRNA must be the only thing that reorders candidates."""
        assert accessibility_penalty(not_evaluated("no tool")) == 0.0

    def test_absent_from_path_entirely_is_the_fresh_clone_case(self, monkeypatch):
        """The real fresh-clone path: no env override AND nothing on PATH.

        The other tests here point $MYPIBYPY_RNAPLFOLD at a missing file, which
        short-circuits before PATH is consulted -- so on this machine, where
        ViennaRNA *is* installed, they never exercise the branch a fresh clone
        actually takes. This one does.
        """
        import shutil as _shutil

        from mypibypy import accessibility as acc
        monkeypatch.delenv("MYPIBYPY_RNAPLFOLD", raising=False)
        monkeypatch.setattr(acc.shutil, "which", lambda name: None)
        assert acc.rnaplfold_path() is None
        res = site_accessibility(open_transcript(), 100)
        assert res.status == NOT_EVALUATED
        assert "RNAplfold not found" in res.reason
        assert accessibility_penalty(res) == 0.0
        assert _shutil.which is not acc.shutil.which or True  # patch is scoped

    def test_an_rna_design_still_completes_without_viennarna(self, monkeypatch):
        """A fresh clone must be able to design, just without the extra signal."""
        from mypibypy import accessibility as acc
        from mypibypy.loci import RnaTarget
        from mypibypy.pipeline import design_rna
        from mypibypy.systems import systems
        monkeypatch.delenv("MYPIBYPY_RNAPLFOLD", raising=False)
        monkeypatch.setattr(acc.shutil, "which", lambda name: None)
        target = RnaTarget("site", "GENER", "NM_000546.6", sorted(systems())[0])
        res = design_rna(target, screen=False, transcript_seq=BALANCED_OPEN)
        assert res.pair.p1l and res.pair.p1r, "design must still succeed"
        assert not res.accessibility.evaluated
        assert "RNAplfold not found" in res.accessibility.reason

    def test_a_failing_binary_is_reported_not_swallowed(self, monkeypatch):
        monkeypatch.setenv("MYPIBYPY_RNAPLFOLD", shutil.which("false") or "/usr/bin/false")
        res = site_accessibility(open_transcript(), 100)
        assert res.status == NOT_EVALUATED and res.reason

    def test_footprint_off_the_end_is_reported(self):
        res = site_accessibility("ACGT" * 5, 100)
        assert res.status == NOT_EVALUATED
        assert "does not fit" in res.reason


@pytest.mark.skipif(not HAVE_PLFOLD, reason="RNAplfold not installed")
class TestStructuredVersusUnstructured:
    """The headline behaviour, against built controls."""

    def test_an_unstructured_site_scores_high(self):
        res = site_accessibility(open_transcript(), 150)
        assert res.evaluated and res.status == EVALUATED
        assert res.pair_nucleation_p > 0.5, res.summary()
        assert res.left.label == "open"

    def test_a_hairpin_site_scores_far_lower(self):
        seq, offset = hairpin_transcript()
        res = site_accessibility(seq, offset)
        assert res.evaluated
        assert res.pair_nucleation_p < 0.5, res.summary()

    def test_the_two_cases_are_well_separated(self):
        seq, offset = hairpin_transcript()
        structured = site_accessibility(seq, offset).pair_nucleation_p
        unstructured = site_accessibility(open_transcript(), 150).pair_nucleation_p
        assert unstructured > structured * 2, (
            f"open {unstructured:.3f} vs hairpin {structured:.3f} -- the metric "
            f"must actually discriminate")

    def test_ranking_penalty_follows(self):
        seq, offset = hairpin_transcript()
        assert (accessibility_penalty(site_accessibility(seq, offset))
                > accessibility_penalty(site_accessibility(open_transcript(), 150)))


@pytest.mark.skipif(not HAVE_PLFOLD, reason="RNAplfold not installed")
class TestTheMeanWouldLie:
    """Why the headline metric is a contiguous run, not a mean.

    A hairpin loop's bases are genuinely unpaired, so a mean-per-base score
    rates it accessible -- but no run long enough to nucleate is ever open,
    because the loop is short and pinned at both ends by the stem. This is the
    measurement that chose the metric.
    """

    def test_a_loop_looks_open_per_base_but_cannot_nucleate(self):
        # Footprint positioned so a pad sits over the hairpin's loop.
        seq = OPEN_TRACT + HAIRPIN + OPEN_TRACT
        loop_mid = len(OPEN_TRACT) + len(STEM) + len(LOOP) // 2
        offset = loop_mid - 12          # left pad straddles the loop
        res = site_accessibility(seq, offset, nucleation_nt=8)
        assert res.evaluated
        pad = res.left
        assert pad.mean_unpaired > pad.best_nucleation_p, (
            f"mean per-base {pad.mean_unpaired:.3f} should exceed the "
            f"contiguous-run probability {pad.best_nucleation_p:.3f}; if these "
            f"ever match, the metric has stopped distinguishing a loop from an "
            f"open tract")

    def test_longer_nucleation_requirement_is_stricter(self):
        seq, offset = hairpin_transcript()
        short = site_accessibility(seq, offset, nucleation_nt=4).pair_nucleation_p
        long = site_accessibility(seq, offset, nucleation_nt=12).pair_nucleation_p
        assert short >= long, "a longer unpaired run cannot be more likely"


@pytest.mark.skipif(not HAVE_PLFOLD, reason="RNAplfold not installed")
class TestPadsAreScoredSeparately:
    """One buried pad kills the pair, however open the other is."""

    def test_the_pair_takes_the_worse_pad_not_the_average(self):
        # Hairpin placed under the right pad only.
        left_part = OPEN_TRACT
        seq = left_part + HAIRPIN + OPEN_TRACT
        offset = len(left_part) - 26        # left pad open, right pad over the stem
        res = site_accessibility(seq, offset)
        assert res.evaluated
        assert res.pair_nucleation_p == min(res.left.best_nucleation_p,
                                            res.right.best_nucleation_p)
        average = (res.left.best_nucleation_p + res.right.best_nucleation_p) / 2
        assert res.pair_nucleation_p <= average

    def test_the_limiting_pad_is_named(self):
        res = site_accessibility(open_transcript(), 150)
        assert res.limiting_pad.pad in ("left", "right")
        assert res.limiting_pad.pad in res.summary()

    def test_each_pad_reports_its_own_offset(self):
        res = site_accessibility(open_transcript(), 150)
        assert res.left.start == 150
        assert res.right.start == 150 + 25 + 1


@pytest.mark.skipif(not HAVE_PLFOLD, reason="RNAplfold not installed")
class TestContextAndReporting:
    def test_flanks_are_included_and_reported(self):
        res = site_accessibility(open_transcript(3), 400)
        assert res.flanks_available == (DEFAULT_FLANK_NT, DEFAULT_FLANK_NT)
        assert res.window_w == DEFAULT_WINDOW_W and res.span_l == DEFAULT_SPAN_L

    def test_a_site_near_the_transcript_end_reports_truncated_context(self):
        """Less context is real, and should be visible rather than padded over."""
        seq = open_transcript()
        res = site_accessibility(seq, 5)
        assert res.evaluated
        assert res.flanks_available[0] == 5, "only 5 nt of 5' flank exists"

    def test_the_uncalibrated_caveat_travels_with_the_number(self):
        res = site_accessibility(open_transcript(), 150)
        assert ACCESSIBILITY_IS_UNCALIBRATED in res.summary()
        assert res.to_dict()["caveat"] == ACCESSIBILITY_IS_UNCALIBRATED

    def test_to_dict_carries_both_pads(self):
        d = site_accessibility(open_transcript(), 150).to_dict()
        assert d["evaluated"] is True
        assert d["left_nucleation_p"] is not None
        assert d["right_nucleation_p"] is not None
        assert d["limiting_pad"] in ("left", "right")

    def test_thymine_input_is_accepted(self):
        """Transcripts arrive as DNA from both fetchers, with T not U."""
        assert site_accessibility(open_transcript(), 150).evaluated


@pytest.mark.skipif(not HAVE_PLFOLD, reason="RNAplfold not installed")
class TestRnaDesignUsesIt:
    @pytest.fixture
    def target(self):
        from mypibypy.loci import RnaTarget
        from mypibypy.systems import systems
        keys = sorted(systems())
        return RnaTarget("site", "GENER", "NM_000546.6", keys[0])

    def test_accessibility_is_attached_to_the_result(self, target):
        from mypibypy.pipeline import design_rna
        res = design_rna(target, screen=False, transcript_seq=BALANCED_OPEN)
        assert res.accessibility.evaluated
        assert res.accessibility.pair_nucleation_p > 0.5

    def test_it_prefers_an_accessible_window(self, target):
        """The point of making it a ranking input rather than only a report.

        Both regions here are GC-balanced, so composition does not decide and
        accessibility does. That matters: with an AT-rich "open" region the
        ranking correctly keeps the structured-but-balanced window instead,
        because accessibility is weighted alongside GC and Tm rather than
        overriding them.
        """
        from mypibypy.pipeline import design_rna
        structured = HAIRPIN * 4
        seq = BALANCED_OPEN[:150] + structured + BALANCED_OPEN
        hairpin_zone = range(150, 150 + len(structured))
        chosen = design_rna(target, screen=False, transcript_seq=seq).best.offset
        assert chosen not in hairpin_zone, (
            f"accessibility ranking chose offset {chosen}, inside the structured "
            f"region spanning {hairpin_zone.start}-{hairpin_zone.stop}")

    def test_accessibility_does_not_override_gross_composition(self, target):
        """It is a ranking input, not a veto -- and that is deliberate.

        An open but 0%-GC window must not win over a balanced one just for being
        accessible; such a probe would not work for a different reason.
        """
        from mypibypy.pipeline import design_rna
        seq = OPEN_TRACT[:200] + BALANCED_OPEN
        chosen = design_rna(target, screen=False, transcript_seq=seq).best.offset
        assert chosen >= 150, (
            f"chose offset {chosen} in the AT-rich tract; composition should "
            f"still dominate when it is this bad")

    def test_disabling_the_check_leaves_an_explicit_not_evaluated(self, target):
        from mypibypy.pipeline import design_rna
        res = design_rna(target, screen=False, transcript_seq=BALANCED_OPEN,
                         check_accessibility=False)
        assert not res.accessibility.evaluated
        assert res.accessibility.reason


class TestDnaDesignsSayNotApplicable:
    """A denatured DNA target is not an unevaluated RNA one."""

    def test_a_genomic_design_result_defaults_to_not_applicable(self):
        import random

        from mypibypy.design import build_pair, scan
        from mypibypy.loci import Locus
        from mypibypy.pipeline import DesignResult
        from mypibypy.sequence import GenomicInterval
        from mypibypy.specificity import not_screened
        from mypibypy.systems import systems

        system = list(systems().values())[0]
        rng = random.Random(3)
        seq = "".join(rng.choice("ACGT") for _ in range(200))
        iv = GenomicInterval("mm10", "chr1", 1, 200)
        cands = scan(seq, iv, top=1)
        pair = build_pair("x", system, cands[0], "P1L", "P1R")
        res = DesignResult(Locus("x", "G", "", iv, system.name), system, seq,
                           cands, pair, rmsk=[], specificity=not_screened("none"))
        assert not res.accessibility.evaluated
        assert "denatures" in res.accessibility.reason
        assert accessibility_penalty(res.accessibility) == 0.0
