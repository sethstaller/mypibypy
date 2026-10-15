"""Cross-hybridization checks, including synthetic positive controls.

No real probe or amplifier sequences live here. Primaries are built from
whichever systems are configured (`systems.local.json`, or the synthetic
`systems.example.json` fallback) plus neutral synthetic target pads, so these
tests exercise the real architecture without committing anyone's sequences.

`TestNucleationDiscriminates` matters most: a check that always returns "ok"
would pass every other test in this file.
"""
from __future__ import annotations

import pytest

from mypibypy.crosshyb import (
    MAX_NUCLEATION_PAD_CONTRIBUTION,
    TOP_LEN,
    arm_regions,
    best_duplex,
    check_order,
    check_pair,
    duplex_dg,
    p2_anneal,
    pool_parts,
    repeat_unit,
)
from mypibypy.sequence import revcomp
from mypibypy.systems import systems

# Neutral synthetic target pads: 25 nt, balanced GC, no long runs. They stand in
# for genomic sequence and are not derived from any real locus.
LEFT_PAD = "ACTGACAGTCATGCAGTACGTACAG"
RIGHT_PAD = "GTCAGTCATGACGTCAGTACTGACT"


@pytest.fixture(scope="module")
def registry():
    reg = systems()
    assert reg, "no amplifier systems configured"
    return reg


@pytest.fixture(scope="module")
def sys_a(registry):
    return list(registry.values())[0]


@pytest.fixture
def pair(sys_a):
    """A well-formed primary pair on the first configured system."""
    return (sys_a.arm_left + revcomp(LEFT_PAD), revcomp(RIGHT_PAD) + sys_a.arm_right)


class TestBestDuplex:
    def test_perfect_revcomp_spans_full_length(self):
        a = "ACGTACGTACGT"
        d = best_duplex(a, revcomp(a))
        assert d.length == len(a)
        assert d.matched == a

    def test_no_complementarity(self):
        assert best_duplex("AAAAAAAA", "AAAAAAAA").length < 2

    def test_reports_coordinates_into_b(self):
        core = "GGGGTTTT"
        b = "AAAA" + revcomp(core) + "AAAA"
        d = best_duplex(core, b)
        assert d.length == len(core)
        assert b[d.b_start : d.b_start + d.length] == revcomp(core)

    def test_prefers_more_stable_run_on_ties(self):
        a = "GCGCGC"
        b = revcomp("GCGCGC") + "TTTTTT" + revcomp("ATATAT")
        assert best_duplex(a, b).matched == "GCGCGC"

    def test_strips_modification_tags(self):
        probe = "ACGTTGCAATCGGATCCAGT"
        assert best_duplex(probe, "/5Cy5/" + revcomp(probe)).length == 20

    def test_dg_is_negative_for_stable_duplex(self):
        assert duplex_dg("GCGCGCGCGC") < duplex_dg("ATATATATAT") < 0


class TestArchitecture:
    """The geometry every other check depends on must actually hold."""

    def test_arms_split_into_14nt_top_and_8nt_middle(self, registry):
        for s in registry.values():
            r = arm_regions(s)
            assert len(r["top_L"]) == TOP_LEN and len(r["top_R"]) == TOP_LEN, s.name
            assert len(r["mid_L"]) == 8 and len(r["mid_R"]) == 8, s.name

    def test_p2_anneal_is_29nt_and_built_from_both_tops(self, registry):
        for s in registry.values():
            anneal, offset = p2_anneal(s)
            r = arm_regions(s)
            assert len(anneal) == 29, f"{s.name}: pi-FISH specifies a 29-nt P2 middle"
            assert anneal.startswith(revcomp(r["top_R"])), s.name
            assert anneal.endswith(revcomp(r["top_L"])), s.name
            assert s.p2[offset : offset + 29] == anneal

    def test_middles_form_a_short_pi_bond(self, registry):
        for s in registry.values():
            d = best_duplex(arm_regions(s)["mid_L"], arm_regions(s)["mid_R"])
            assert 2 <= d.length <= 4, f"{s.name}: paper specifies 2-4 bp between halves"

    def test_p2_and_p3_are_tandem_repeats(self, registry):
        for s in registry.values():
            assert repeat_unit(s.p2) is not None, s.name
            assert repeat_unit(s.p3) is not None, s.name

    def test_p4_binds_the_p3_repeat(self, registry):
        """The last rung of the cascade: P4 must be complementary to P3's repeat."""
        for s in registry.values():
            unit = repeat_unit(s.p3)
            assert best_duplex(s.p4, unit).length >= 20, s.name

    def test_pool_covers_every_configured_channel(self, registry):
        names = {p.name for p in pool_parts(list(registry.values()))}
        for s in registry.values():
            assert any(s.name in n for n in names), s.name
            assert s.p4_name in names, s.p4_name


class TestWellFormedPairPasses:
    def test_clean_pair_has_no_flags(self, pair, sys_a):
        rep = check_pair(*pair, sys_a)
        assert rep.ok, rep.summary()

    def test_lone_primary_pad_adds_little_to_the_p2_site(self, pair, sys_a):
        rep = check_pair(*pair, sys_a)
        nuc = [f for f in rep.findings if f.kind == "lone_primary_nucleation"]
        assert len(nuc) == 2
        for f in nuc:
            assert not f.flagged, f
            assert f.duplex.length <= MAX_NUCLEATION_PAD_CONTRIBUTION

    def test_top_region_contributes_exactly_the_designed_half_site(self, pair, sys_a):
        anneal, _ = p2_anneal(sys_a)
        for primary in pair:
            assert best_duplex(primary, anneal).length == TOP_LEN

    def test_designed_pi_bond_is_not_flagged(self, pair, sys_a):
        partner = check_pair(*pair, sys_a).worst("primary_vs_partner")
        assert partner is not None and not partner.flagged


class TestNucleationDiscriminates:
    """Prove the checks fire. Without this, an always-ok check looks perfect."""

    def test_nucleation_catches_pad_that_also_binds_p2(self, pair, sys_a):
        # A pad complementary to P2's rc(top_R) half: this primary offers its own
        # 14-nt top PLUS a second long duplex, so P2 no longer needs the partner.
        anneal, _ = p2_anneal(sys_a)
        bad_p1l = sys_a.arm_left + revcomp(anneal[:14]) + "ACGTACGTACG"
        f = [x for x in check_pair(bad_p1l, pair[1], sys_a).findings
             if x.kind == "lone_primary_nucleation" and "left" in x.subject][0]
        assert f.duplex.length > MAX_NUCLEATION_PAD_CONTRIBUTION
        assert f.flagged, "a pad that independently binds P2 must be flagged"

    def test_pad_complementary_to_amplifier_is_flagged(self, pair, sys_a):
        unit = repeat_unit(sys_a.p3)
        bad_p1l = sys_a.arm_left + revcomp(unit[:20])
        rep = check_pair(bad_p1l, pair[1], sys_a)
        assert [f for f in rep.flagged if f.kind == "pad_vs_amplifier"], rep.summary()

    def test_the_two_tops_together_cover_the_whole_p2_site(self, pair, sys_a):
        """The coincidence mechanism itself: neither primary alone suffices."""
        anneal, _ = p2_anneal(sys_a)
        left = best_duplex(pair[0], anneal).length
        right = best_duplex(pair[1], anneal).length
        assert left == right == TOP_LEN
        assert left + right == len(anneal) - 1, "14 + 14, 1-nt spacer unpaired"
        assert max(left, right) < len(anneal)


class TestOrderLevel:
    def test_sense_antisense_classified_not_flagged(self, pair, sys_a, registry):
        other = list(registry.values())[1] if len(registry) > 1 else sys_a
        anti = other.arm_left + revcomp(revcomp(LEFT_PAD))
        rep = check_order({"P1L_sense": (pair[0], sys_a), "P1L_anti": (anti, other)})
        assert "sense_antisense_pair" in {f.kind for f in rep.findings}
        assert rep.ok, "expected sense/antisense complementarity is informational"

    def test_unrelated_primaries_do_not_flag(self, pair, sys_a):
        rep = check_order({"a": (pair[0], sys_a), "b": (pair[1], sys_a)})
        assert rep.ok, rep.summary()
