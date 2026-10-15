"""The design-to-result record.

Mostly these tests pin down that the store is append-only and that nothing in it
pretends to be a score. With three designs in hand, the thing most worth
guarding against is this module growing a threshold.

zFISHer is treated as read-only: the fixtures write their own CSV in zFISHer's
column format rather than reaching into a real session folder.
"""
from __future__ import annotations

import csv
import json
import random

import pytest

from mypibypy.design import build_pair, scan
from mypibypy.loci import Locus
from mypibypy.outcomes import (
    CALLS,
    DesignFacts,
    ImagingResult,
    Outcome,
    OutcomeStore,
    facts_from_result,
    read_zfisher_puncta,
    resolve_store_path,
)
from mypibypy.pipeline import DesignResult
from mypibypy.sequence import GenomicInterval
from mypibypy.specificity import SCREENED, SpecificityReport, not_screened
from mypibypy.systems import systems

# zFISHer's own column order -- zfisher/constants.py PUNCTA_CSV_COLUMNS.
ZFISHER_COLUMNS = ["puncta_id", "Z", "Y", "X", "Nucleus_ID", "Intensity", "SNR", "Source"]


@pytest.fixture
def store(tmp_path):
    return OutcomeStore(tmp_path / "outcomes.jsonl")


@pytest.fixture(scope="module")
def design_result():
    reg = systems()
    assert reg, "no amplifier systems configured"
    system = list(reg.values())[0]
    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(200))
    iv = GenomicInterval("mm10", "chr3", 500_000, 500_199)
    locus = Locus("prom", "Gene9", "synthetic", iv, system.name)
    cands = scan(seq, iv, top=2)
    pair = build_pair(locus.name, system, cands[0], "P1L_x", "P1R_x")
    return DesignResult(locus, system, seq, cands, pair, rmsk=[],
                        specificity=not_screened("no database in tests"))


def write_puncta(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ZFISHER_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    return path


def punctum(i, nucleus, intensity=1000.0, snr=5.0):
    return {"puncta_id": i, "Z": 1.0, "Y": 2.0, "X": 3.0, "Nucleus_ID": nucleus,
            "Intensity": intensity, "SNR": snr, "Source": "detected"}


class TestDesignFacts:
    def test_flattens_a_design_result(self, design_result):
        f = facts_from_result(design_result)
        assert f.design_id == "Gene9_prom"
        assert f.p1l and f.p1r and f.left_pad and f.right_pad
        assert f.left_gc is not None and f.left_tm is not None
        assert f.specificity_status == "NOT_SCREENED"
        assert f.specificity_verdict == "not screened"

    def test_an_unscreened_design_is_recorded_as_unscreened(self, design_result):
        """The store must not launder a missing screen into a blank field."""
        f = facts_from_result(design_result)
        assert f.specificity_status and "NOT SCREENED" in f.specificity_summary

    def test_a_screened_design_keeps_its_verdict(self, design_result):
        import dataclasses
        screened = dataclasses.replace(
            design_result,
            specificity=SpecificityReport(status=SCREENED, backend="bowtie v1 -v 3 -a",
                                          database="mm10"))
        f = facts_from_result(screened)
        assert f.specificity_verdict == "pass"
        assert "bowtie" in f.specificity_summary

    def test_crosshyb_is_optional_and_recorded_when_given(self, design_result):
        from mypibypy.crosshyb import check_pair
        rep = check_pair(design_result.pair.p1l, design_result.pair.p1r,
                         design_result.system)
        f = facts_from_result(design_result, crosshyb_report=rep)
        assert f.crosshyb_summary
        assert f.crosshyb_worst_bp is not None and f.crosshyb_worst_bp > 0
        assert facts_from_result(design_result).crosshyb_worst_bp is None

    def test_custom_design_id_is_used(self, design_result):
        assert facts_from_result(design_result, "ORDER_2026_01").design_id == "ORDER_2026_01"


class TestReadingZFisher:
    def test_summarizes_a_puncta_table(self, tmp_path):
        rows = [punctum(1, "1", 100, 4.0), punctum(2, "1", 200, 6.0),
                punctum(3, "2", 300, 8.0)]
        res = read_zfisher_puncta(write_puncta(tmp_path / "a_puncta.csv", rows))
        assert res.n_puncta == 3
        assert res.n_nuclei == 2
        assert res.puncta_per_nucleus == pytest.approx(1.5)
        assert res.median_snr == pytest.approx(6.0)
        assert res.median_intensity == pytest.approx(200.0)

    def test_puncta_outside_a_nucleus_do_not_count_as_nuclei(self, tmp_path):
        """zFISHer marks unassigned puncta with 0 / -1 / blank."""
        rows = [punctum(1, "0"), punctum(2, "-1"), punctum(3, ""), punctum(4, "7")]
        res = read_zfisher_puncta(write_puncta(tmp_path / "b_puncta.csv", rows))
        assert res.n_puncta == 4
        assert res.n_nuclei == 1

    def test_an_empty_table_is_a_real_result_not_an_error(self, tmp_path):
        """No puncta is exactly the outcome worth recording for a failed probe."""
        res = read_zfisher_puncta(write_puncta(tmp_path / "c_puncta.csv", []),
                                  call="no_signal")
        assert res.n_puncta == 0
        assert res.n_nuclei is None and res.median_snr is None
        assert res.call == "no_signal"

    def test_missing_optional_columns_are_tolerated(self, tmp_path):
        """That CSV comes from another project on its own release schedule."""
        p = tmp_path / "d_puncta.csv"
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["puncta_id", "Z", "Y", "X"])
            w.writeheader()
            w.writerows([{"puncta_id": 1, "Z": 0, "Y": 0, "X": 0}])
        res = read_zfisher_puncta(p)
        assert res.n_puncta == 1 and res.median_snr is None

    def test_unparseable_numbers_do_not_abort_the_read(self, tmp_path):
        rows = [punctum(1, "1", 100, 4.0), punctum(2, "1", 200, 6.0)]
        rows[0]["SNR"] = "NA"
        res = read_zfisher_puncta(write_puncta(tmp_path / "e_puncta.csv", rows))
        assert res.n_puncta == 2 and res.median_snr == pytest.approx(6.0)

    def test_a_missing_file_says_so(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_zfisher_puncta(tmp_path / "nope.csv")

    def test_an_invalid_call_is_rejected(self, tmp_path):
        p = write_puncta(tmp_path / "f_puncta.csv", [punctum(1, "1")])
        with pytest.raises(ValueError, match="call must be one of"):
            read_zfisher_puncta(p, call="great")


class TestStore:
    def test_round_trips_a_record(self, store, design_result):
        store.record_design(design_result)
        got = store.get("Gene9_prom")
        assert got is not None
        assert got.design.p1l == design_result.pair.p1l
        assert got.imaging.call == "not_imaged"

    def test_empty_store_reads_as_empty(self, store):
        assert store.all() == []
        assert "no records yet" in store.summary()

    def test_attach_imaging_to_a_recorded_design(self, store, design_result, tmp_path):
        store.record_design(design_result)
        p = write_puncta(tmp_path / "g_puncta.csv",
                         [punctum(i, "1", 100, 5.0) for i in range(12)])
        store.attach_zfisher_report("Gene9_prom", p, call="signal",
                                    notes="clear nuclear spots")
        rec = store.get("Gene9_prom")
        assert rec.imaging.call == "signal" and rec.imaging.n_puncta == 12
        assert rec.imaging.notes == "clear nuclear spots"
        assert rec.design.p1l, "attaching imaging must not drop the design"

    def test_attaching_to_an_unknown_design_lists_what_is_known(self, store, design_result):
        store.record_design(design_result)
        with pytest.raises(KeyError, match="Gene9_prom"):
            store.attach_imaging("does_not_exist", ImagingResult())

    def test_the_store_is_append_only(self, store, design_result, tmp_path):
        """Nothing is ever overwritten -- the earlier version stays on disk."""
        store.record_design(design_result)
        p = write_puncta(tmp_path / "h_puncta.csv", [punctum(1, "1")])
        store.attach_zfisher_report("Gene9_prom", p, call="weak")
        store.attach_zfisher_report("Gene9_prom", p, call="signal")
        assert len(store.path.read_text().splitlines()) == 3
        assert store.get("Gene9_prom").imaging.call == "signal", "latest wins"
        calls = [r.imaging.call for r in store.history("Gene9_prom")]
        assert calls == ["not_imaged", "weak", "signal"]

    def test_created_at_survives_an_update(self, store, design_result, tmp_path):
        store.record_design(design_result)
        created = store.get("Gene9_prom").created_at
        p = write_puncta(tmp_path / "i_puncta.csv", [punctum(1, "1")])
        store.attach_zfisher_report("Gene9_prom", p, call="signal")
        rec = store.get("Gene9_prom")
        assert rec.created_at == created

    def test_several_designs_coexist(self, store, design_result):
        store.record_design(design_result, "d1")
        store.record_design(design_result, "d2")
        assert {r.design.design_id for r in store.all()} == {"d1", "d2"}

    def test_a_corrupt_line_is_reported_with_the_path(self, store, design_result):
        store.record_design(design_result)
        with store.path.open("a") as fh:
            fh.write("{not json\n")
        with pytest.raises(ValueError, match="unreadable record"):
            store.all()

    def test_comments_and_blank_lines_are_skipped(self, store, design_result):
        store.record_design(design_result)
        with store.path.open("a") as fh:
            fh.write("\n# written by hand\n")
        assert len(store.all()) == 1

    def test_stored_json_is_plain_and_inspectable(self, store, design_result):
        """This file must still be readable by hand in five years."""
        store.record_design(design_result)
        line = store.path.read_text().splitlines()[0]
        d = json.loads(line)
        assert set(d) == {"design", "imaging", "created_at", "updated_at"}
        assert isinstance(d["design"]["p1l"], str)
        assert all(not isinstance(v, dict) for v in d["design"].values())


class TestItDoesNotPretendToScore:
    """The constraint that matters: three designs cannot support a threshold."""

    def test_summary_counts_and_says_the_n_is_too_small(self, store, design_result, tmp_path):
        p = write_puncta(tmp_path / "j_puncta.csv", [punctum(1, "1")])
        for i in range(3):
            store.record_design(design_result, f"d{i}")
            store.attach_zfisher_report(f"d{i}", p, call="signal")
        text = store.summary()
        assert "3 design(s) recorded, 3 imaged" in text
        assert "far too few to calibrate" in text

    def test_there_is_no_score_anywhere_in_the_record(self, store, design_result):
        store.record_design(design_result)
        d = json.loads(store.path.read_text().splitlines()[0])
        flat = {**d["design"], **d["imaging"]}
        assert not [k for k in flat if "score" in k.lower() or "rank" in k.lower()], (
            "outcomes records facts; scoring belongs nowhere near three data points")

    def test_table_puts_design_next_to_result(self, store, design_result, tmp_path):
        store.record_design(design_result)
        p = write_puncta(tmp_path / "k_puncta.csv",
                         [punctum(i, "1", 100, 7.0) for i in range(5)])
        store.attach_zfisher_report("Gene9_prom", p, call="signal")
        (row,) = store.table()
        assert row["specificity"] == "not screened"
        assert row["call"] == "signal" and row["n_puncta"] == 5
        assert row["median_snr"] == pytest.approx(7.0)

    def test_csv_export_round_trips(self, store, design_result, tmp_path):
        store.record_design(design_result)
        out = store.write_csv(tmp_path / "out.csv")
        with out.open() as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 1 and rows[0]["design_id"] == "Gene9_prom"

    def test_every_call_label_is_accepted(self, store, design_result, tmp_path):
        p = write_puncta(tmp_path / "l_puncta.csv", [punctum(1, "1")])
        for call in CALLS:
            assert read_zfisher_puncta(p, call=call).call == call


class TestStorePath:
    def test_explicit_path_wins(self, tmp_path):
        assert resolve_store_path(tmp_path / "x.jsonl") == tmp_path / "x.jsonl"

    def test_env_var_is_used(self, tmp_path, monkeypatch):
        monkeypatch.setenv("MYPIBYPY_OUTCOMES", str(tmp_path / "env.jsonl"))
        assert resolve_store_path() == tmp_path / "env.jsonl"

    def test_default_is_the_gitignored_local_name(self, monkeypatch):
        monkeypatch.delenv("MYPIBYPY_OUTCOMES", raising=False)
        assert resolve_store_path().name == "outcomes.local.jsonl"
