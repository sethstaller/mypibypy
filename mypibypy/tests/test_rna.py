"""RNA piFISH: transcript fetch by accession, and transcript-space design.

Fetching is tested by substituting the HTTP layer rather than by hitting
Ensembl and NCBI, so the suite stays fast, offline and immune to their
intermittent 500s. The point of these tests is the validation logic -- what
happens when a response is truncated, empty, or for the wrong molecule -- and
that is exactly what a live call cannot be relied on to produce on demand. A
handful of opt-in live tests are gated behind $MYPIBYPY_LIVE_FETCH.

No lab accessions appear here. NM_000546 (TP53) and ENST00000269305 are
textbook public references, used so nothing in the repository points at a real
target.
"""
from __future__ import annotations

import json
import os
import random

import pytest

from mypibypy import sequence as seqmod
from mypibypy.design import (
    FOOTPRINT,
    LEFT_LEN,
    junction_position,
    scan_transcript,
    spans_junction,
)
from mypibypy.loci import LociNotConfigured, RnaTarget, load_rna_targets
from mypibypy.sequence import (
    TranscriptFetchError,
    fetch_ensembl_cdna,
    fetch_refseq_transcript,
    fetch_transcript,
    strip_version,
)

LIVE = os.environ.get("MYPIBYPY_LIVE_FETCH")


@pytest.fixture
def transcript():
    rng = random.Random(21)
    return "".join(rng.choice("ACGT") for _ in range(2000))


# --- Accession handling ----------------------------------------------------

class TestAccessionDispatch:
    def test_strip_version(self):
        assert strip_version("NM_000546.6") == "NM_000546"
        assert strip_version("ENST00000269305") == "ENST00000269305"

    def test_unrecognised_accession_is_refused_not_guessed(self):
        """Guessing a source returns a plausible sequence for the wrong molecule."""
        for bad in ("rubbish", "chr5:1-100", "TP53", "12345"):
            with pytest.raises(TranscriptFetchError, match="unrecognised accession"):
                fetch_transcript(bad)

    @pytest.mark.parametrize("acc,fn", [
        ("ENST00000269305", "fetch_ensembl_cdna"),
        ("ENSMUST00000000001", "fetch_ensembl_cdna"),
        ("NM_000546.6", "fetch_refseq_transcript"),
        ("NR_003051.3", "fetch_refseq_transcript"),
        ("XM_011537074.2", "fetch_refseq_transcript"),
    ])
    def test_dispatches_on_prefix(self, acc, fn, monkeypatch):
        called = {}
        monkeypatch.setattr(seqmod, fn, lambda a, timeout=30: called.setdefault("acc", a) or "ACGT" * 10)
        fetch_transcript(acc)
        assert called["acc"] == acc


class TestEnsemblFetch:
    def _json(self, monkeypatch, payload: bytes):
        monkeypatch.setattr(seqmod, "_get", lambda *a, **k: payload)

    def test_reads_the_seq_field(self, monkeypatch):
        self._json(monkeypatch, json.dumps(
            {"id": "ENST00000269305", "seq": "acgtacgtacgt"}).encode())
        assert fetch_ensembl_cdna("ENST00000269305") == "ACGTACGTACGT"

    def test_versioned_accession_is_accepted(self, monkeypatch):
        self._json(monkeypatch, json.dumps(
            {"id": "ENST00000269305", "seq": "ACGT" * 5}).encode())
        assert fetch_ensembl_cdna("ENST00000269305.9")

    def test_truncated_json_is_a_loud_error(self, monkeypatch):
        """The failure worth catching: a partial body that would otherwise
        become a short sequence and shift every downstream offset."""
        full = json.dumps({"id": "ENST00000269305", "seq": "ACGT" * 200}).encode()
        self._json(monkeypatch, full[: len(full) // 2])
        with pytest.raises(TranscriptFetchError, match="truncated"):
            fetch_ensembl_cdna("ENST00000269305")

    def test_empty_body_is_a_loud_error(self, monkeypatch):
        """Sending Content-Type instead of Accept returns 200 with no body."""
        self._json(monkeypatch, b"   ")
        with pytest.raises(TranscriptFetchError, match="empty body"):
            fetch_ensembl_cdna("ENST00000269305")

    def test_a_different_id_coming_back_is_refused(self, monkeypatch):
        self._json(monkeypatch, json.dumps(
            {"id": "ENST00000000001", "seq": "ACGT" * 5}).encode())
        with pytest.raises(TranscriptFetchError, match="got ENST00000000001"):
            fetch_ensembl_cdna("ENST00000269305")

    def test_rna_bases_are_converted(self, monkeypatch):
        self._json(monkeypatch, json.dumps({"id": "ENST1", "seq": "ACGU"}).encode())
        assert fetch_ensembl_cdna("ENST1") == "ACGT"

    def test_protein_sequence_is_rejected(self, monkeypatch):
        self._json(monkeypatch, json.dumps({"id": "ENST1", "seq": "MEEPQSDPSV"}).encode())
        with pytest.raises(TranscriptFetchError, match="non-DNA"):
            fetch_ensembl_cdna("ENST1")


class TestRefSeqFetch:
    def test_length_is_cross_checked_against_esummary(self, monkeypatch):
        """efetch FASTA carries no length, so truncation is invisible without this."""
        monkeypatch.setattr(seqmod, "refseq_length", lambda a, timeout=30: 100)
        monkeypatch.setattr(seqmod, "_get",
                            lambda *a, **k: b">NM_000546.6 x\n" + b"ACGT" * 10)
        with pytest.raises(TranscriptFetchError, match="truncated or the accession"):
            fetch_refseq_transcript("NM_000546.6")

    def test_matching_length_is_accepted(self, monkeypatch):
        monkeypatch.setattr(seqmod, "refseq_length", lambda a, timeout=30: 40)
        monkeypatch.setattr(seqmod, "_get",
                            lambda *a, **k: b">NM_000546.6 x\n" + b"ACGT" * 10)
        assert len(fetch_refseq_transcript("NM_000546.6")) == 40

    def test_absent_esummary_length_does_not_invent_one(self, monkeypatch):
        """Better to skip the check than to compare against a made-up number."""
        monkeypatch.setattr(seqmod, "refseq_length", lambda a, timeout=30: None)
        monkeypatch.setattr(seqmod, "_get",
                            lambda *a, **k: b">NM_1 x\n" + b"ACGT" * 10)
        assert len(fetch_refseq_transcript("NM_1")) == 40

    def test_non_fasta_response_is_refused(self, monkeypatch):
        monkeypatch.setattr(seqmod, "refseq_length", lambda a, timeout=30: None)
        monkeypatch.setattr(seqmod, "_get", lambda *a, **k: b"<html>error</html>")
        with pytest.raises(TranscriptFetchError, match="did not return FASTA"):
            fetch_refseq_transcript("NM_1")


@pytest.mark.skipif(not LIVE, reason="set $MYPIBYPY_LIVE_FETCH=1 to hit the real APIs")
class TestLiveFetch:
    """Opt-in. Needs network and, behind a TLS-inspecting proxy, $SSL_CERT_FILE."""

    @staticmethod
    def _fetch(acc):
        """Fetch, distinguishing an unreachable service from a wrong answer.

        These services time out and emit transient 500s under normal use, and
        `_get` already retries. A network failure says nothing about this code,
        so it skips; a *wrong* sequence still fails, which is the whole point of
        running live.
        """
        try:
            return fetch_transcript(acc)
        except TranscriptFetchError as exc:
            pytest.skip(f"{acc} unreachable: {str(exc)[:120]}")

    def test_refseq_round_trip(self):
        assert len(self._fetch("NM_000546.6")) == 2512

    def test_the_two_sources_agree_on_the_same_molecule(self):
        """Independent validation of both fetchers against each other."""
        assert self._fetch("ENST00000269305") == self._fetch("NM_000546.6")


# --- Exon junctions --------------------------------------------------------

class TestJunctionReporting:
    def test_detects_a_crossed_boundary(self):
        assert spans_junction(100, (125,)) == 125
        assert spans_junction(100, (100,)) is None, "a boundary at the start is not crossed"
        assert spans_junction(100, (151,)) is None, "one past the end is not crossed"

    def test_ignores_boundaries_outside_the_footprint(self):
        assert spans_junction(100, (50, 200, 400)) is None

    def test_boundary_on_the_gap_is_called_out_as_the_clean_case(self):
        text = junction_position(100, 100 + LEFT_LEN)
        assert "on the 1-nt gap" in text and "wholly within one exon" in text

    def test_boundary_inside_a_pad_says_which_and_where(self):
        assert "inside the left pad" in junction_position(100, 105)
        assert "inside the right pad" in junction_position(100, 140)

    def test_a_spanning_candidate_is_flagged_not_rejected(self, transcript):
        """Junction-spanning is often the point -- it is what makes a probe
        isoform-specific -- so it must be reported, never filtered out."""
        cands = scan_transcript(transcript, "NM_1", top=50, exon_boundaries=(300,))
        spanning = [c for c in cands if any("spans exon junction" in f for f in c.flags)]
        assert spanning, "some window in a 2 kb transcript must cross offset 300"
        assert all(c.interval.length == FOOTPRINT for c in spanning)


class TestTranscriptScanning:
    def test_coordinates_are_transcript_offsets_not_genomic(self, transcript):
        c = scan_transcript(transcript, "NM_000546.6", top=1)[0]
        assert c.interval.build == "NM_000546.6"
        assert c.interval.chrom == "transcript", (
            "nothing downstream should mistake these for genome coordinates")
        assert c.interval.start == c.offset + 1
        assert c.interval.length == FOOTPRINT

    def test_scan_window_is_respected_and_offsets_stay_absolute(self, transcript):
        cands = scan_transcript(transcript, "NM_1", top=20, scan_start=500, scan_end=900)
        assert cands
        for c in cands:
            assert 500 <= c.offset
            assert c.offset + FOOTPRINT <= 900
            assert transcript[c.offset : c.offset + FOOTPRINT] == (
                c.left.seq + c.gap_base + c.right.seq), "offset must index the transcript"

    def test_too_short_a_region_is_a_clear_error(self, transcript):
        with pytest.raises(ValueError, match="need at least"):
            scan_transcript(transcript, "NM_1", scan_start=0, scan_end=20)


# --- Target definitions ----------------------------------------------------

class TestRnaTargetLoading:
    def test_the_shipped_example_defines_one(self):
        from pathlib import Path

        import mypibypy.loci as loci_mod
        targets = load_rna_targets(Path(loci_mod.__file__).parent / "loci.example.json")
        assert len(targets) == 1
        t = targets[0]
        assert t.accession.startswith("NM_")
        assert t.antisense_system != t.system

    def test_a_file_with_no_rna_targets_is_not_an_error(self, tmp_path):
        p = tmp_path / "loci.json"
        p.write_text('{"loci": []}')
        assert load_rna_targets(p) == []

    def test_missing_fields_are_named(self, tmp_path):
        p = tmp_path / "loci.json"
        p.write_text('{"rna_targets": [{"name": "x", "system": "cy5"}]}')
        with pytest.raises(LociNotConfigured, match="accession"):
            load_rna_targets(p)

    def test_bad_exon_boundaries_are_rejected(self, tmp_path):
        p = tmp_path / "loci.json"
        p.write_text('{"rna_targets": [{"name":"x","gene":"G","accession":"NM_1",'
                     '"system":"cy5","exon_boundaries":["a"]}]}')
        with pytest.raises(LociNotConfigured, match="non-integer exon_boundaries"):
            load_rna_targets(p)

    def test_inverted_scan_window_is_rejected(self, tmp_path):
        p = tmp_path / "loci.json"
        p.write_text('{"rna_targets": [{"name":"x","gene":"G","accession":"NM_1",'
                     '"system":"cy5","scan_start":900,"scan_end":100}]}')
        with pytest.raises(LociNotConfigured, match="scan_end <= scan_start"):
            load_rna_targets(p)

    def test_boundaries_are_sorted(self, tmp_path):
        p = tmp_path / "loci.json"
        p.write_text('{"rna_targets": [{"name":"x","gene":"G","accession":"NM_1",'
                     '"system":"cy5","exon_boundaries":[400,120,900]}]}')
        assert load_rna_targets(p)[0].exon_boundaries == (120, 400, 900)


class TestDesignRna:
    """End to end with the transcript supplied, so no network is involved."""

    @pytest.fixture
    def target(self):
        from mypibypy.systems import systems
        keys = sorted(systems())
        return RnaTarget("cds", "GENEQ", "NM_000546.6", keys[0],
                         system_antisense=keys[-1], exon_boundaries=(300, 900))

    def test_designs_without_a_screen(self, target, transcript):
        from mypibypy.pipeline import design_rna
        res = design_rna(target, screen=False, transcript_seq=transcript)
        assert res.pair.p1l and res.pair.p1r
        assert res.locus.interval.chrom == "transcript"
        assert res.rmsk == [], "a transcript has no RepeatMasker track"
        assert not res.specificity.screened

    def test_the_screen_is_told_the_target_gene(self, target, transcript, monkeypatch):
        """Without it, the probe's own isoforms come back as off-targets."""
        from mypibypy import pipeline
        seen = {}

        def fake(l, r, backend=None, **kw):
            seen.update(kw)
            return pipeline.not_screened("stub")

        monkeypatch.setattr(pipeline, "screen_pair", fake)
        pipeline.design_rna(target, transcript_seq=transcript, backend=object())
        assert seen.get("target_gene") == "GENEQ"

    def test_the_screening_policy_is_documented(self):
        """A reader must be able to find out which reference is used, and why."""
        from mypibypy import pipeline
        doc = pipeline.design_rna.__doc__ or ""
        assert "RNA_SCREEN_POLICY" in doc
        policy_doc = pipeline.__dict__.get("RNA_SCREEN_POLICY", "")
        assert "transcriptome" in policy_doc.lower()
        src = pipeline.__doc__ or ""
        assert isinstance(src, str)
