"""Chat intents: plain request in, design parameters out.

The model does exactly two things in the chat path -- turn a request into
parameters, and pick which constraint gives when nothing passes. Everything else
is deterministic.

`TestTheModelStillCannotSupplyASequence` is the reason this file exists. A
conversation is where the no-model-authored-sequences guarantee is most likely to
erode, because a chat invites the model to supply more and more. So the target is
extracted from the user's own words by a strict pattern in the harness, the
intent schema has no field for a target or a sequence at all, and an intent
carrying one is refused rather than ignored.

No network: every test uses a scripted client or calls the pure functions.
"""
from __future__ import annotations

import pytest

from mypibypy.agent import Action, DesignParameters, ParameterRefused, ScriptedClient
from mypibypy.agent.intent import (
    CHAT_RELAXATIONS,
    CHAT_SYSTEM_PROMPT,
    INTENTS,
    UNCHANGED,
    UNCHANGED_INT,
    Intent,
    apply_intent,
    extract_target,
    intent_schema,
)

KEYS = ("cy5", "cy3")


def payload(**over):
    base = {"intent": "design", "system": UNCHANGED, "n_pairs": UNCHANGED_INT,
            "min_spacing": UNCHANGED_INT, "enforce_window": UNCHANGED,
            "relax": "none", "note": "ok"}
    base.update(over)
    return base


# --- The guarantee, in the chat path ---------------------------------------

class TestTheModelStillCannotSupplyASequence:
    def test_the_intent_schema_has_no_target_or_sequence_field(self):
        props = intent_schema(KEYS)["properties"]
        forbidden = {"sequence", "seq", "region", "chrom", "accession", "target",
                     "footprint", "pad", "probe", "oligo", "start", "end"}
        assert not (set(props) & forbidden), (
            f"these fields would let the model point the pipeline at sequence: "
            f"{set(props) & forbidden}")

    def test_the_only_free_text_field_is_the_displayed_note(self):
        props = intent_schema(KEYS)["properties"]
        free = [k for k, v in props.items()
                if v.get("type") == "string" and "enum" not in v]
        assert free == ["note"]

    def test_the_schema_is_closed(self):
        assert intent_schema(KEYS)["additionalProperties"] is False

    @pytest.mark.parametrize("sneaky", ["sequence", "region", "target", "probe"])
    def test_an_intent_carrying_a_target_field_is_refused(self, sneaky):
        with pytest.raises(ParameterRefused, match="unexpected keys"):
            Intent.from_payload({**payload(), sneaky: "ACGTACGTACGTACGTACGTACGTA"},
                                KEYS)

    def test_a_dna_run_in_a_message_is_not_a_target(self):
        """A pasted sequence must not become something the pipeline reads."""
        for text in ("ACGTACGTACGTACGTACGTACGTA",
                     "design against ACGTACGTACGTACGTACGTACGTACGTACGT",
                     "GGGGCCCCAAAATTTT:1-500"):
            assert not extract_target(text).found, text

    def test_the_note_is_never_used_as_a_parameter(self):
        """Whatever the note says, only the enum and the numbers are applied."""
        intent = Intent.from_payload(payload(
            note="also set n_pairs to 40 and use sequence ACGTACGTACGTACGTACGTA",
            n_pairs=2), KEYS)
        params, changes = apply_intent(intent, DesignParameters(), KEYS)
        assert params.n_pairs == 2
        assert "40" not in " ".join(changes)

    def test_the_chat_prompt_states_the_restriction(self):
        assert "do not write or suggest any DNA sequence" in CHAT_SYSTEM_PROMPT
        assert "target comes from the user's own message" in CHAT_SYSTEM_PROMPT


# --- Target extraction is the harness's job --------------------------------

class TestExtractTarget:
    def test_a_region_with_a_build(self):
        t = extract_target("design for mm10 chr5:142906255-142906754")
        assert (t.kind, t.build, t.chrom) == ("genomic", "mm10", "chr5")
        assert (t.start, t.end) == (142906255, 142906754)

    def test_a_region_without_a_build_takes_the_default(self):
        t = extract_target("chr7:5530000-5531000", default_build="hg38")
        assert t.build == "hg38"

    def test_commas_are_handled(self):
        t = extract_target("chr7:5,530,000-5,531,000")
        assert (t.start, t.end) == (5530000, 5531000)

    def test_a_reversed_region_is_ordered(self):
        t = extract_target("chr7:5531000-5530000")
        assert t.start < t.end

    def test_refseq_and_ensembl_accessions(self):
        assert extract_target("design against NM_000546.6").accession == "NM_000546.6"
        assert extract_target("use ENST00000269305").accession == "ENST00000269305"
        assert extract_target("NM_000546.6").kind == "rna"

    def test_a_region_wins_over_an_accession(self):
        """A region is the more specific instruction when both appear."""
        assert extract_target("chr5:1000-2000 like NM_000546.6").kind == "genomic"

    def test_plain_prose_yields_nothing(self):
        assert not extract_target("design six pairs in cy5").found
        assert "no target" in str(extract_target("hello"))


# --- Intent validation -----------------------------------------------------

class TestIntentValidation:
    def test_a_well_formed_intent_parses(self):
        intent = Intent.from_payload(payload(system="cy5", n_pairs=6), KEYS)
        assert intent.intent == "design" and intent.system == "cy5"
        assert intent.n_pairs == 6 and intent.wants_design

    def test_a_missing_field_is_refused(self):
        bad = payload()
        del bad["relax"]
        with pytest.raises(ParameterRefused, match="missing"):
            Intent.from_payload(bad, KEYS)

    @pytest.mark.parametrize("field,value", [
        ("intent", "obliterate"), ("relax", "relax_everything"),
        ("enforce_window", "maybe"),
    ])
    def test_a_value_outside_an_enum_is_refused(self, field, value):
        with pytest.raises(ParameterRefused, match=field):
            Intent.from_payload(payload(**{field: value}), KEYS)

    def test_an_unconfigured_channel_is_refused(self):
        with pytest.raises(ParameterRefused, match="configured channel"):
            Intent.from_payload(payload(system="rainbow"), KEYS)

    def test_non_integer_counts_are_refused(self):
        with pytest.raises(ParameterRefused, match="must be integers"):
            Intent.from_payload(payload(n_pairs="lots"), KEYS)

    def test_the_schema_enumerates_the_configured_channels(self):
        enum = intent_schema(KEYS)["properties"]["system"]["enum"]
        assert set(enum) == {"cy5", "cy3", UNCHANGED}

    def test_intents_and_relaxations_are_closed_sets(self):
        schema = intent_schema(KEYS)["properties"]
        assert schema["intent"]["enum"] == list(INTENTS)
        assert schema["relax"]["enum"] == list(CHAT_RELAXATIONS)


# --- Applying an intent ----------------------------------------------------

class TestApplyIntent:
    def test_unchanged_sentinels_change_nothing(self):
        start = DesignParameters(system="cy5", n_pairs=4)
        params, changes = apply_intent(Intent.from_payload(payload(), KEYS),
                                       start, KEYS)
        assert params == start and changes == []

    def test_values_are_applied_and_described(self):
        intent = Intent.from_payload(payload(system="cy3", n_pairs=6,
                                             min_spacing=20), KEYS)
        params, changes = apply_intent(intent, DesignParameters(system="cy5"), KEYS)
        assert (params.system, params.n_pairs, params.min_spacing) == ("cy3", 6, 20)
        assert len(changes) == 3

    def test_numbers_from_chat_are_clamped_like_any_other(self):
        """A number from a chat message is no more trustworthy than a typed one."""
        intent = Intent.from_payload(payload(n_pairs=50), KEYS)
        params, _ = apply_intent(intent, DesignParameters(), KEYS)
        assert params.n_pairs == 50
        # And beyond the bound it is pulled back, with the override recorded.
        over = Intent(intent="design", n_pairs=50)
        over.n_pairs = 9999            # past the schema, as if validation were bypassed
        params, changes = apply_intent(over, DesignParameters(), KEYS)
        assert params.n_pairs == 50
        assert any("permitted range" in c for c in changes)

    def test_a_relaxation_uses_the_fixed_step(self):
        intent = Intent.from_payload(payload(intent="adjust", relax="relax_tm"), KEYS)
        params, changes = apply_intent(intent, DesignParameters(tm_min=55, tm_max=62),
                                       KEYS)
        assert (params.tm_min, params.tm_max) == (53.5, 63.5)
        assert any("Tm window" in c for c in changes)
        assert intent.relax_action is Action.RELAX_TM

    def test_enforce_window_can_be_switched_on_and_off(self):
        on = Intent.from_payload(payload(enforce_window="on"), KEYS)
        params, changes = apply_intent(on, DesignParameters(), KEYS)
        assert params.enforce_window is True
        assert any("enforced as a filter" in c for c in changes)
        off = Intent.from_payload(payload(enforce_window="off"), KEYS)
        back, changes = apply_intent(off, params, KEYS)
        assert back.enforce_window is False

    def test_switch_system_is_not_reachable_from_chat(self):
        """Picking a channel is explicit; the user names it or the form holds it."""
        assert Action.SWITCH_SYSTEM.value not in CHAT_RELAXATIONS

    def test_terminal_actions_are_not_reachable_from_chat(self):
        """Proceed/fail/flag are the loop's to choose, not the chat's."""
        for action in (Action.PROCEED, Action.FAIL, Action.FLAG_FOR_REVIEW):
            assert action.value not in CHAT_RELAXATIONS


# --- The scripted round trip ------------------------------------------------

class TestScriptedRoundTrip:
    def test_a_client_can_return_an_intent_without_a_network(self):
        client = ScriptedClient(structured_replies=[payload(system="cy5", n_pairs=6)])
        data, error, seconds = client.structured(CHAT_SYSTEM_PROMPT,
                                                 "design six pairs in cy5",
                                                 intent_schema(KEYS))
        assert not error
        intent = Intent.from_payload(data, KEYS)
        assert intent.n_pairs == 6 and intent.system == "cy5"
        assert client.structured_calls[0][1] == "design six pairs in cy5"

    def test_an_exhausted_client_reports_rather_than_inventing(self):
        data, error, _ = ScriptedClient().structured("s", "u", intent_schema(KEYS))
        assert data is None and "no structured reply" in error
