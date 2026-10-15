"""Gene-level collapsing for transcriptome screens.

A gene appears once in a genome and ~8.5 times in GENCODE. Screened against a
transcriptome without gene awareness, a probe finds its own target in eight
isoforms, the pair logic reads those as eight off-target sites, and a perfectly
good probe is rejected. These tests cover the machinery that fixes that, and --
just as importantly -- that it changes nothing for a genome screen.

`TestGenomeScreensAreUnchanged` is the guard that matters most here. Gene
awareness is detected from the reference names, so a plain chromosome FASTA must
take exactly the old path. If that regresses, every genomic design silently
changes behaviour.

The references are synthetic and tiny. No real sequences, and no dependency on
the multi-GB indices on this machine.
"""
from __future__ import annotations

import random
import shutil
import subprocess

import pytest

from mypibypy.sequence import revcomp
from mypibypy.specificity import (
    Backend,
    Hit,
    SpecificityReport,
    collapse_sites_by_gene,
    detect_gene_parser,
    gencode_gene,
    no_gene,
    pair_hits,
    screen_pair,
)

LEFT = "ACTGACAGTCATGCAGTACGTACAG"
RIGHT = "GTCAGTCATGACGTCAGTACTGACT"

# A real GENCODE v50 header, used verbatim as a parsing fixture. It is public
# annotation metadata, not sequence.
REAL_HEADER = ("ENST00000832824.1|ENSG00000290825.2|-|-|"
               "DDX11L16-260|DDX11L16|1379|lncRNA|")


def gencode_header(gene: str, n: int, length: int = 1500, gene_id: str | None = None,
                   symbol: str | None = None) -> str:
    """A GENCODE-shaped reference name for isoform `n` of `gene`."""
    sym = gene if symbol is None else symbol
    gid = gene_id or f"ENSG{abs(hash(gene)) % 10**11:011d}.1"
    return (f"ENST{n:011d}.1|{gid}|-|-|{gene}-{200 + n}|{sym}|{length}|protein_coding|")


def hit(query="left", chrom="chr1", start=1000, strand="+", mismatches=0,
        aligned=25, gene=""):
    return Hit(query=query, chrom=chrom, start=start, end=start + aligned - 1,
               strand=strand, mismatches=mismatches, gaps=0, aligned=aligned,
               query_len=25, q_start=1, q_end=aligned, gene=gene)


class TestGeneParsing:
    def test_parses_a_real_gencode_header(self):
        assert gencode_gene(REAL_HEADER) == "DDX11L16"

    def test_falls_back_to_gene_id_when_the_symbol_is_a_dash(self):
        assert gencode_gene("ENST1.1|ENSG42.7|-|-|name-1|-|900|lncRNA|") == "ENSG42.7"

    def test_plain_chromosome_names_have_no_gene(self):
        for name in ("chr1", "chrX", "chr1_KI270706v1_random", "decoy1", "NC_000001.11"):
            assert gencode_gene(name) is None, name

    def test_too_few_fields_is_not_a_gencode_name(self):
        assert gencode_gene("gene|symbol") is None
        assert gencode_gene("a|b|c|d|e") is None

    def test_no_gene_resolver_always_returns_none(self):
        assert no_gene(REAL_HEADER) is None

    def test_detection_picks_gencode_when_any_name_parses(self):
        assert detect_gene_parser([REAL_HEADER]) is gencode_gene
        assert detect_gene_parser(["chr1", REAL_HEADER]) is gencode_gene

    def test_detection_falls_back_for_chromosome_references(self):
        assert detect_gene_parser(["chr1", "chr2", "chrM"]) is no_gene

    def test_detection_on_no_hits_at_all_is_not_gene_aware(self):
        assert detect_gene_parser([]) is no_gene


class TestCollapsing:
    """Isoform redundancy collapses; genuinely different geometry does not."""

    def _sites(self, specs):
        """specs: list of (gene, chrom, left_start, right_start)."""
        left = [hit("left", chrom=c, start=ls, gene=g) for g, c, ls, _ in specs]
        right = [hit("right", chrom=c, start=rs, gene=g) for g, c, _, rs in specs]
        out = []
        for lh, rh in zip(left, right):
            out += pair_hits([lh], [rh], LEFT, RIGHT, window_nt=2000)
        return out

    def test_eight_isoforms_with_one_geometry_collapse_to_one_site(self):
        specs = [("GENEX", gencode_header("GENEX", i), 100 + 37 * i, 126 + 37 * i)
                 for i in range(8)]
        sites = self._sites(specs)
        assert len(sites) == 8, "one per transcript before collapsing"
        collapsed = collapse_sites_by_gene(sites)
        assert len(collapsed) == 1
        assert collapsed[0].gene == "GENEX"
        assert collapsed[0].gap_nt == 1

    def test_a_different_geometry_in_another_isoform_survives(self):
        """An exon spliced out of one isoform changes the gap, and that is a
        different site, not a duplicate -- so it must not be merged away."""
        adjacent = ("GENEY", gencode_header("GENEY", 1), 100, 126)      # gap 1
        spliced = ("GENEY", gencode_header("GENEY", 2), 100, 700)       # gap 575
        collapsed = collapse_sites_by_gene(self._sites([adjacent, spliced]))
        assert len(collapsed) == 2
        assert {s.gap_nt for s in collapsed} == {1, 575}

    def test_two_genes_stay_two_sites(self):
        specs = [("GENEA", gencode_header("GENEA", 1), 100, 126),
                 ("GENEB", gencode_header("GENEB", 2), 100, 126)]
        collapsed = collapse_sites_by_gene(self._sites(specs))
        assert {s.gene for s in collapsed} == {"GENEA", "GENEB"}

    def test_sites_without_gene_information_are_never_merged(self):
        """Two identical-looking genomic sites are two real loci."""
        specs = [("", "chr1", 100, 126), ("", "chr7", 5000, 5026)]
        sites = self._sites(specs)
        assert len(collapse_sites_by_gene(sites)) == 2

    def test_differing_mismatch_counts_are_not_merged(self):
        a = hit("left", chrom=gencode_header("GENEZ", 1), start=100, gene="GENEZ")
        b = hit("right", chrom=gencode_header("GENEZ", 1), start=126, gene="GENEZ")
        c = hit("left", chrom=gencode_header("GENEZ", 2), start=100, gene="GENEZ",
                mismatches=2)
        d = hit("right", chrom=gencode_header("GENEZ", 2), start=126, gene="GENEZ")
        sites = pair_hits([a], [b], LEFT, RIGHT) + pair_hits([c], [d], LEFT, RIGHT)
        assert len(collapse_sites_by_gene(sites)) == 2


class TestReportingSeparatesSitesFromGenes:
    def test_offtarget_genes_counts_genes_not_sites(self):
        specs = [("GENEA", gencode_header("GENEA", i), 100, 126) for i in range(4)]
        specs += [("GENEB", gencode_header("GENEB", 9), 100, 126)]
        sites = []
        for g, c, ls, rs in specs:
            sites += pair_hits([hit("left", chrom=c, start=ls, gene=g)],
                               [hit("right", chrom=c, start=rs, gene=g)],
                               LEFT, RIGHT)
        rep = SpecificityReport(status="SCREENED", backend="b", database="d",
                                gene_aware=True, left_seq=LEFT, right_seq=RIGHT,
                                pair_sites=sites, n_sites_before_collapse=len(sites))
        assert len(rep.pair_sites) == 5
        assert rep.offtarget_genes == ["GENEA", "GENEB"], "5 sites, 2 genes"

    def test_summary_reports_both_counts_and_the_collapse(self):
        sites = pair_hits([hit("left", chrom=gencode_header("GENEA", 1), gene="GENEA")],
                          [hit("right", chrom=gencode_header("GENEA", 1), start=1026,
                               gene="GENEA")], LEFT, RIGHT)
        rep = SpecificityReport(status="SCREENED", backend="b", database="d",
                                gene_aware=True, target_gene="MYGENE",
                                on_target_gene_hits=16, left_seq=LEFT, right_seq=RIGHT,
                                pair_sites=sites, n_sites_before_collapse=9)
        text = rep.summary()
        assert "off-target gene(s)" in text and "GENEA" in text
        assert "intended gene MYGENE were excluded" in text
        assert "9 raw sites collapsed to 1" in text

    def test_a_genome_report_says_nothing_about_genes(self):
        rep = SpecificityReport(status="SCREENED", backend="b", database="d",
                                gene_aware=False, left_seq=LEFT, right_seq=RIGHT)
        text = rep.summary()
        assert "gene" not in text.lower(), text
        assert rep.offtarget_genes == []

    def test_to_dict_exposes_the_gene_counts(self):
        rep = SpecificityReport(status="SCREENED", backend="b", database="d",
                                gene_aware=True, target_gene="G")
        d = rep.to_dict()
        assert d["gene_aware"] is True and d["target_gene"] == "G"
        assert d["n_offtarget_genes"] == 0


# --- Integration against real aligners on synthetic references -------------

HAVE_BOWTIE = bool(shutil.which("bowtie") and shutil.which("bowtie-build"))
FILLER = "GATCTTAGCCAGTTACGATCGGTAACCTGAATCGGCTATTGCAACGTTAGC"


def _build(tmp_path, records: dict[str, str], name="ref"):
    fa = tmp_path / f"{name}.fa"
    with fa.open("w") as fh:
        for header, seq in records.items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")
    idx = str(tmp_path / name)
    subprocess.run(["bowtie-build", "-q", str(fa), idx], capture_output=True, check=True)
    return Backend(tool="bowtie", database=idx, executable=shutil.which("bowtie"))


def _transcript(payload: str, length: int = 1200, offset: int = 300) -> str:
    rng = random.Random(len(payload) + offset)
    body = "".join(rng.choice("ACGT") for _ in range(length))
    return body[:offset] + payload + body[offset + len(payload):]


@pytest.mark.skipif(not HAVE_BOWTIE, reason="no bowtie v1")
class TestTranscriptomeScreen:
    """The actual failure mode, end to end."""

    @pytest.fixture
    def isoform_reference(self, tmp_path):
        """MYGENE in 8 isoforms, all carrying the footprint, plus two other genes."""
        footprint = LEFT + "T" + RIGHT
        records = {}
        for i in range(8):
            records[gencode_header("MYGENE", i)] = _transcript(footprint, offset=200 + 50 * i)
        for j, other in enumerate(("OTHER1", "OTHER2")):
            records[gencode_header(other, 100 + j)] = _transcript("", offset=400)
        return _build(tmp_path, records)

    def test_without_target_gene_a_good_probe_looks_terrible(self, isoform_reference):
        """This is the bug: the probe's own gene reads as many off-target sites."""
        rep = screen_pair(LEFT, RIGHT, backend=isoform_reference)
        assert rep.screened and rep.gene_aware
        assert rep.pair_sites, "its own eight isoforms pair up"
        assert rep.verdict == "fail"
        assert rep.offtarget_genes == ["MYGENE"], "and they are all the intended gene"

    def test_naming_the_target_gene_clears_it(self, isoform_reference):
        rep = screen_pair(LEFT, RIGHT, backend=isoform_reference, target_gene="MYGENE")
        assert rep.screened and rep.gene_aware
        assert rep.pair_sites == [], f"nothing off-target remains: {rep.summary()}"
        assert rep.offtarget_genes == []
        assert rep.verdict == "pass"
        assert rep.on_target_gene_hits >= 8, "the on-target hits are counted, not hidden"

    def test_the_excluded_hits_are_reported_not_silently_dropped(self, isoform_reference):
        rep = screen_pair(LEFT, RIGHT, backend=isoform_reference, target_gene="MYGENE")
        assert "intended gene MYGENE were excluded" in rep.summary()
        assert str(rep.on_target_gene_hits) in rep.summary()

    def test_gene_matching_is_case_insensitive(self, isoform_reference):
        rep = screen_pair(LEFT, RIGHT, backend=isoform_reference, target_gene="mygene")
        assert rep.pair_sites == [] and rep.on_target_gene_hits >= 8

    def test_a_real_offtarget_gene_still_fails(self, tmp_path):
        """Excluding the target must not excuse a genuine second gene."""
        footprint = LEFT + "T" + RIGHT
        records = {gencode_header("MYGENE", 1): _transcript(footprint, offset=200),
                   gencode_header("BADGENE", 2): _transcript(footprint, offset=500)}
        rep = screen_pair(LEFT, RIGHT, backend=_build(tmp_path, records),
                          target_gene="MYGENE")
        assert rep.verdict == "fail"
        assert rep.offtarget_genes == ["BADGENE"]

    def test_isoform_sites_collapse_so_the_count_is_readable(self, tmp_path):
        """Eight isoforms of one off-target gene is one problem, not eight."""
        footprint = LEFT + "T" + RIGHT
        records = {gencode_header("BADGENE", i): _transcript(footprint, offset=250)
                   for i in range(8)}
        rep = screen_pair(LEFT, RIGHT, backend=_build(tmp_path, records))
        assert rep.n_sites_before_collapse >= 8
        assert len(rep.pair_sites) < rep.n_sites_before_collapse, rep.summary()
        assert rep.offtarget_genes == ["BADGENE"]

    def test_reports_name_the_gene_not_the_80_character_header(self, isoform_reference):
        rep = screen_pair(LEFT, RIGHT, backend=isoform_reference)
        rendered = str(rep.pair_sites[0])
        assert "MYGENE" in rendered
        assert "|" not in rendered, "a pipe-delimited header helps no reader"


@pytest.mark.skipif(not HAVE_BOWTIE, reason="no bowtie v1")
class TestGenomeScreensAreUnchanged:
    """Gene awareness must be invisible on a plain chromosome reference."""

    @pytest.fixture
    def genome_reference(self, tmp_path):
        seq = (FILLER * 400)[:20000]
        planted = seq[:4000] + LEFT + "T" + RIGHT + seq[4051:]
        return _build(tmp_path, {"chr1": planted, "chr2": seq}, name="genome")

    def test_a_chromosome_reference_is_not_gene_aware(self, genome_reference):
        rep = screen_pair(LEFT, RIGHT, backend=genome_reference)
        assert rep.screened
        assert rep.gene_aware is False
        assert rep.target_gene == ""
        assert all(h.gene == "" for h in (*rep.left_hits, *rep.right_hits))

    def test_sites_are_not_collapsed_on_a_genome(self, genome_reference):
        rep = screen_pair(LEFT, RIGHT, backend=genome_reference)
        assert rep.n_sites_before_collapse == len(rep.pair_sites)

    def test_the_planted_site_is_still_found_and_still_fails(self, genome_reference):
        rep = screen_pair(LEFT, RIGHT, backend=genome_reference)
        assert [s for s in rep.pair_sites if s.signal_capable]
        assert rep.verdict == "fail"

    def test_target_gene_is_ignored_without_gene_information(self, genome_reference):
        """Passing it must not quietly drop hits it cannot possibly match."""
        with_gene = screen_pair(LEFT, RIGHT, backend=genome_reference,
                                target_gene="SOMEGENE")
        without = screen_pair(LEFT, RIGHT, backend=genome_reference)
        assert with_gene.on_target_gene_hits == 0
        assert len(with_gene.pair_sites) == len(without.pair_sites)
        assert with_gene.verdict == without.verdict

    def test_sites_still_render_the_chromosome_name(self, genome_reference):
        rep = screen_pair(LEFT, RIGHT, backend=genome_reference)
        assert "chr1" in str(rep.pair_sites[0])
