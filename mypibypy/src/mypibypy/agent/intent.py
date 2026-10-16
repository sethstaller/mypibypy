"""Turning a plain-English request into design parameters.

This is the model's *second* job -- the first being "which constraint gives when
nothing passes". Together those are the only two things it does; everything else
in the chat path is deterministic.

THE GUARANTEE HOLDS HERE TOO, and this is where it is most likely to erode,
because a conversation invites the model to supply more and more. Two rules keep
it structural rather than hopeful:

* **The model never supplies the target.** Not the region, not the accession, not
  a junction. It returns parameters only. The target is pulled out of the user's
  own text by `extract_target`, a strict regex in the harness, so there is no
  path by which a model-authored string becomes something the pipeline reads
  sequence from. `INTENT_SCHEMA` has no target field at all.
* **The only free-text field is `note`**, which is displayed to the user and
  never parsed. Every other field is an enum or a bounded number, and the
  numbers go through `actions.clamp` exactly like parameters typed into the form.

`UNCHANGED` sentinels rather than nullable fields: a strict JSON schema with
`additionalProperties: false` is easier to constrain when every property is
always present, and a sentinel is unambiguous in the run record.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, replace

from .actions import Action, DesignParameters, ParameterRefused, clamp

#: Sentinel values meaning "the user did not ask to change this".
UNCHANGED_INT = -1
UNCHANGED = "unchanged"

#: What the user is asking for. Closed.
INTENTS = ("design", "adjust", "explain", "unclear")

#: Relaxations reachable from chat. A subset of `Action`: `switch_system` is
#: excluded because picking a channel is explicit (the user names it or the form
#: holds it), and the terminal actions are the loop's to choose, not the chat's.
CHAT_RELAXATIONS = ("none", Action.RELAX_TM.value, Action.RELAX_GC.value,
                    Action.WIDEN_REGION.value, Action.REDUCE_PAIRS.value,
                    Action.ACCEPT_SINGLE_ARM_HITS.value)


def intent_schema(system_keys: tuple[str, ...] = ()) -> dict:
    """The structured-output schema, with the configured channels enumerated.

    Built per call so `system` is a closed enum of channels that actually exist
    -- a model cannot name a dye the lab does not have, and `clamp` would refuse
    it anyway.
    """
    return {
        "type": "object",
        "properties": {
            "intent": {
                "type": "string", "enum": list(INTENTS),
                "description": ("design = run a design now; adjust = change a "
                                "constraint on the existing one; explain = the "
                                "user is asking a question; unclear = ask them."),
            },
            "system": {
                "type": "string", "enum": [*system_keys, UNCHANGED],
                "description": "Amplifier channel, or unchanged.",
            },
            "n_pairs": {
                "type": "integer", "minimum": UNCHANGED_INT, "maximum": 50,
                "description": f"Pairs to tile, or {UNCHANGED_INT} for unchanged.",
            },
            "min_spacing": {
                "type": "integer", "minimum": UNCHANGED_INT, "maximum": 10000,
                "description": f"Minimum nt between footprints, or "
                               f"{UNCHANGED_INT} for unchanged.",
            },
            "enforce_window": {
                "type": "string", "enum": ["on", "off", UNCHANGED],
                "description": ("Whether the GC/Tm window should REJECT windows. "
                                "Off by default; the pipeline normally ranks on "
                                "GC/Tm distance rather than rejecting."),
            },
            "relax": {
                "type": "string", "enum": list(CHAT_RELAXATIONS),
                "description": ("Which constraint to loosen one step, if the user "
                                "asked for that. 'none' otherwise."),
            },
            "note": {
                "type": "string",
                "description": ("One sentence back to the user. Displayed as your "
                                "comment; never interpreted by the code."),
            },
        },
        "required": ["intent", "system", "n_pairs", "min_spacing",
                     "enforce_window", "relax", "note"],
        "additionalProperties": False,
    }


@dataclass
class Intent:
    """A validated request. Carries no sequence and no target, by construction."""

    intent: str
    system: str = UNCHANGED
    n_pairs: int = UNCHANGED_INT
    min_spacing: int = UNCHANGED_INT
    enforce_window: str = UNCHANGED
    relax: str = "none"
    note: str = ""

    @property
    def wants_design(self) -> bool:
        return self.intent in ("design", "adjust")

    @property
    def relax_action(self) -> Action | None:
        return None if self.relax == "none" else Action(self.relax)

    @staticmethod
    def from_payload(payload: dict, system_keys: tuple[str, ...] = ()) -> Intent:
        """Validate against the closed schema. Refuses anything outside it."""
        schema = intent_schema(system_keys)
        allowed = set(schema["properties"])
        missing = allowed - set(payload)
        if missing:
            raise ParameterRefused(f"intent is missing {sorted(missing)}")
        extra = set(payload) - allowed
        if extra:
            raise ParameterRefused(
                f"intent carries unexpected keys {sorted(extra)}; the schema is "
                f"closed, and a conversation is exactly where that matters")
        for field, values in (("intent", INTENTS), ("relax", CHAT_RELAXATIONS),
                              ("enforce_window", ("on", "off", UNCHANGED))):
            if payload[field] not in values:
                raise ParameterRefused(
                    f"{field}={payload[field]!r} is not one of {list(values)}")
        if system_keys and payload["system"] not in (*system_keys, UNCHANGED):
            raise ParameterRefused(
                f"system={payload['system']!r} is not a configured channel "
                f"{sorted(system_keys)}")
        try:
            n_pairs = int(payload["n_pairs"])
            min_spacing = int(payload["min_spacing"])
        except (TypeError, ValueError) as exc:
            raise ParameterRefused("n_pairs and min_spacing must be integers") from exc
        return Intent(intent=payload["intent"], system=payload["system"],
                      n_pairs=n_pairs, min_spacing=min_spacing,
                      enforce_window=payload["enforce_window"],
                      relax=payload["relax"], note=str(payload["note"]))


def apply_intent(intent: Intent, params: DesignParameters,
                 known_systems: tuple[str, ...] = ()) -> tuple[DesignParameters, list[str]]:
    """Fold an intent into the current parameters. Returns (parameters, changes).

    Every number still goes through `clamp`, because a number from a chat message
    is no more trustworthy than one typed into a form -- less, if anything.
    """
    changes: list[str] = []
    new = params

    if intent.system != UNCHANGED and intent.system != params.system:
        new = replace(new, system=intent.system)
        changes.append(f"channel → {intent.system}")
    if intent.n_pairs != UNCHANGED_INT and intent.n_pairs != params.n_pairs:
        new = replace(new, n_pairs=intent.n_pairs)
        changes.append(f"pairs → {intent.n_pairs}")
    if intent.min_spacing != UNCHANGED_INT and intent.min_spacing != params.min_spacing:
        new = replace(new, min_spacing=intent.min_spacing)
        changes.append(f"minimum spacing → {intent.min_spacing} nt")
    if intent.enforce_window != UNCHANGED:
        wanted = intent.enforce_window == "on"
        if wanted != params.enforce_window:
            new = replace(new, enforce_window=wanted)
            changes.append("GC/Tm window " + ("enforced as a filter" if wanted
                                              else "no longer enforced"))

    action = intent.relax_action
    if action is not None:
        from .actions import apply_action
        new, change = apply_action(action, new, known_systems)
        changes.append(change)

    new, clamps = clamp(new, known_systems)
    changes.extend(str(c) for c in clamps)
    return new, changes


# --- Target extraction: harness-side, never the model's ---------------------

#: A genomic region, strictly. Chromosome names are letters, digits and
#: underscores only, and the coordinates are digits -- so nothing shaped like a
#: sequence can come through this, whoever typed it.
REGION_RE = re.compile(r"\b(?:(?P<build>hg38|hg19|mm10|mm39|t2t)\s+)?"
                       r"(?P<chrom>chr[0-9A-Za-z_]{1,20})\s*:\s*"
                       r"(?P<start>[\d,]{1,15})\s*[-–]\s*(?P<end>[\d,]{1,15})\b")

#: An accession, by the prefixes `sequence.fetch_transcript` actually dispatches
#: on. Anchored and bounded; a run of bases cannot match it.
ACCESSION_RE = re.compile(r"\b(?P<acc>(?:ENS[A-Z]*T|NM_|NR_|XM_|XR_)[0-9]{4,12}"
                          r"(?:\.[0-9]{1,3})?)\b")


@dataclass(frozen=True)
class ExtractedTarget:
    """A target found in the user's own words, by the harness."""

    kind: str                     # "genomic" | "rna" | "none"
    build: str = ""
    chrom: str = ""
    start: int = 0
    end: int = 0
    accession: str = ""

    @property
    def found(self) -> bool:
        return self.kind != "none"

    def __str__(self) -> str:
        if self.kind == "genomic":
            return f"{self.build or '?'} {self.chrom}:{self.start:,}-{self.end:,}"
        if self.kind == "rna":
            return self.accession
        return "no target in the message"


def extract_target(message: str, default_build: str = "hg38") -> ExtractedTarget:
    """Pull a region or accession out of the user's text. Deterministic, no model.

    This exists so the model never supplies a target. It is the one piece of the
    request that points the pipeline at sequence, and it comes from the user's
    own words through a strict pattern -- not from anything generated.
    """
    region = REGION_RE.search(message)
    if region:
        start = int(region.group("start").replace(",", ""))
        end = int(region.group("end").replace(",", ""))
        if start > end:
            start, end = end, start
        return ExtractedTarget(kind="genomic",
                               build=region.group("build") or default_build,
                               chrom=region.group("chrom"), start=start, end=end)
    accession = ACCESSION_RE.search(message)
    if accession:
        return ExtractedTarget(kind="rna", accession=accession.group("acc"))
    return ExtractedTarget(kind="none")


CHAT_SYSTEM_PROMPT = """\
You help a molecular biology laboratory set up piFISH probe designs. This is \
routine laboratory work: short DNA oligonucleotides that label a genomic locus \
or a transcript for imaging under a microscope.

Your ONLY job is to turn the user's request into design PARAMETERS. You do not \
design anything, you do not write or suggest any DNA sequence, and you do not \
choose the target -- the target comes from the user's own message and is \
extracted by the software. You have no field in which to put any of those.

Fill in only what the user actually asked for and leave everything else at its \
unchanged sentinel. If they ask to loosen a constraint ("too few passed, relax \
the Tm"), set `relax` to the matching value. If you cannot tell what they want, \
use intent "unclear" and ask them in `note`.

Note that the GC/Tm window is NOT applied as a filter by default: the pipeline \
ranks windows by distance from it. Only set enforce_window to "on" if the user \
clearly wants candidates rejected for falling outside it.

`note` is one sentence shown to the user as your comment. It is never \
interpreted by the software.
"""
