"""Presentation helpers for the Streamlit app.

Streamlit itself is not imported here -- `mypibypy.ui` holds the logic precisely
so it can be tested without a browser or a server. What the app does with these
rows (layout, widgets) was verified by running it.

The point worth guarding: a missing screen must render as words, not as a blank
cell. A table is exactly where "not screened" would otherwise be mistaken for
"passed".
"""
from __future__ import annotations

import random

import pytest

from mypibypy.design import build_tiled_pairs, scan, select_tiled
from mypibypy.loci import Locus
from mypibypy.pipeline import DesignResult
from mypibypy.sequence import GenomicInterval
from mypibypy.specificity import SCREENED, SpecificityReport, not_screened
from mypibypy.systems import systems
from mypibypy.ui import candidate_rows, unscreened_pairs, workbook_bytes


@pytest.fixture(scope="module")
def system():
    reg = systems()
    assert reg, "no amplifier systems configured"
    return list(reg.values())[0]


def make_result(system, n_pairs=1, screened=False, accessibility=None):
    rng = random.Random(77)
    seq = "".join(rng.choice("ACGT") for _ in range(600))
    iv = GenomicInterval("hg38", "chr1", 1, 600)
    locus = Locus("target", "GENE", "ui test", iv, system.name)
    cands = scan(seq, iv, top=len(seq))
    pairs = build_tiled_pairs("target", "GENE", system, select_tiled(cands, n_pairs))
    if screened:
        screens = [SpecificityReport(status=SCREENED, backend="bowtie v1 -v 3 -a",
                                     database="hg38") for _ in pairs]
    else:
        screens = [not_screened("no index supplied") for _ in pairs]
    access = [accessibility or _not_evaluated() for _ in pairs]
    return DesignResult(locus, system, seq, cands[:5], pairs[0], [], screens[0],
                        accessibility=access[0], pairs=pairs,
                        pair_screens=screens, pair_accessibility=access)


def _not_evaluated():
    from mypibypy.accessibility import not_evaluated
    return not_evaluated("RNAplfold not installed")


class TestCandidateRows:
    def test_one_row_per_ordered_pair(self, system):
        assert len(candidate_rows(make_result(system, n_pairs=4))) == 4

    def test_columns_cover_what_the_brief_asked_for(self, system):
        row = candidate_rows(make_result(system))[0]
        for column in ("GC L", "GC R", "Tm L", "Tm R", "accessibility",
                       "specificity", "P1L", "P1R"):
            assert column in row, f"{column} missing from the candidate table"

    def test_an_unscreened_pair_says_so_in_words(self, system):
        """A blank cell would read as a pass."""
        row = candidate_rows(make_result(system, screened=False))[0]
        assert row["specificity"] == "NOT SCREENED"

    def test_a_screened_pair_shows_its_verdict(self, system):
        row = candidate_rows(make_result(system, screened=True))[0]
        assert row["specificity"] == "pass"

    def test_unevaluated_accessibility_says_so_in_words(self, system):
        row = candidate_rows(make_result(system))[0]
        assert row["accessibility"] == "not evaluated"

    def test_evaluated_accessibility_is_a_number(self, system):
        from mypibypy.accessibility import (
            EVALUATED,
            AccessibilityResult,
            PadAccessibility,
        )
        pad = PadAccessibility("left", 0, 25, 8, 0.734, 0.9, 0.8)
        right = PadAccessibility("right", 26, 25, 8, 0.812, 0.9, 0.8)
        access = AccessibilityResult(status=EVALUATED, left=pad, right=right,
                                     tool="RNAplfold")
        row = candidate_rows(make_result(system, accessibility=access))[0]
        assert row["accessibility"] == pytest.approx(0.734, abs=1e-3)

    def test_rows_are_numbered_from_one(self, system):
        rows = candidate_rows(make_result(system, n_pairs=3))
        assert [r["#"] for r in rows] == [1, 2, 3]

    def test_names_are_distinct_across_tiled_pairs(self, system):
        rows = candidate_rows(make_result(system, n_pairs=3))
        names = [r["P1L"] for r in rows] + [r["P1R"] for r in rows]
        assert len(set(names)) == len(names)


class TestUnscreenedCount:
    def test_counts_unscreened_pairs(self, system):
        assert unscreened_pairs(make_result(system, n_pairs=3)) == 3
        assert unscreened_pairs(make_result(system, n_pairs=3, screened=True)) == 0


class TestWorkbookBytes:
    def test_produces_a_real_xlsx_in_memory(self, system):
        data = workbook_bytes([make_result(system)], "ui test")
        assert data[:2] == b"PK", "an xlsx is a zip"
        assert len(data) > 4000

    def test_nothing_is_written_into_the_working_tree(self, system, tmp_path,
                                                      monkeypatch):
        """The workbook carries probe sequences; the download is the only exit."""
        monkeypatch.chdir(tmp_path)
        workbook_bytes([make_result(system)], "ui test")
        assert list(tmp_path.iterdir()) == [], (
            "workbook_bytes must leave no file behind")

    def test_tiled_results_render(self, system):
        data = workbook_bytes([make_result(system, n_pairs=5)], "tiled")
        assert data[:2] == b"PK"
