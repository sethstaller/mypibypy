"""Model clients for the decision step.

The loop depends on `DecisionClient`, not on any vendor, so the endpoint is
swappable: `AnthropicDecisionClient` talks to the Claude API (or any
Anthropic-compatible gateway via `base_url`), and `ScriptedClient` returns
canned decisions for tests. Nothing in `loop` imports `anthropic`.

Three things in the Anthropic path are not optional and are easy to get wrong:

1. **`stop_reason == "refusal"` is checked BEFORE `response.content` is read.**
   A refusal returns HTTP 200 with content that may not hold what you expect, so
   reading content first and checking afterwards is the wrong order. Probe design
   is routine laboratory work, but the prompt talks about oligonucleotides
   targeting human genes, which is close enough to the biology classifier's
   territory that refusals must be handled as a normal branch rather than an
   exception.
2. **Server-side fallback is opted into** (`fallbacks="default"` with the
   matching beta), so a policy decline is re-run on a fallback model inside the
   same call rather than failing the run.
3. **The action enum arrives as structured output**, not as prose to be parsed.

A refusal that survives the fallback chain is not an error here: it returns
`None`, and the loop falls back to its deterministic default. The model is an
advisor in this design, so losing it degrades the run rather than ending it.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Protocol

from .actions import ACTION_SCHEMA, Decision, ParameterRefused

#: Default model. Pinned rather than "latest" so a run record names something
#: reproducible; override per-client when the lab wants a cheaper or local one.
DEFAULT_MODEL = "claude-opus-5-5"

#: Opus 5.5 defaults to `medium` effort, unlike its predecessor. Set explicitly
#: so a model change does not silently alter how hard it thinks about a routing
#: decision. `low` is plenty for picking one of nine actions from a prepared
#: evidence summary.
DEFAULT_EFFORT = "low"

#: Opt into server-side fallback on a policy decline. The scalar "default" form
#: routes by refusal category, so no model list has to be maintained here.
FALLBACK_BETA = "server-side-fallback-2026-07-01"


@dataclass
class ModelReply:
    """What a decision call returned, including the shape of a non-answer."""

    decision: Decision | None
    refused: bool = False
    refusal_category: str = ""
    refusal_explanation: str = ""
    input_tokens: int = 0
    output_tokens: int = 0
    model: str = ""
    error: str = ""

    @property
    def usable(self) -> bool:
        return self.decision is not None


class DecisionClient(Protocol):
    """Anything that can turn an evidence summary into one `Decision`."""

    name: str

    def decide(self, system: str, evidence: str) -> ModelReply:
        ...


class ScriptedClient:
    """A canned client for tests and dry runs. Makes no network calls.

    Every test in this package uses this. The real client is exercised only by
    shape assertions on the request it would build, because no credentials are
    configured in this environment and a test suite that needs them is a test
    suite that stops working.
    """

    name = "scripted"

    def __init__(self, decisions: list[dict] | None = None,
                 refuse: bool = False, error: str = ""):
        self._decisions = list(decisions or [])
        self.refuse = refuse
        self.error = error
        self.calls: list[tuple[str, str]] = []

    def decide(self, system: str, evidence: str) -> ModelReply:
        self.calls.append((system, evidence))
        if self.error:
            return ModelReply(decision=None, error=self.error, model=self.name)
        if self.refuse:
            return ModelReply(decision=None, refused=True,
                              refusal_category="bio",
                              refusal_explanation="scripted refusal",
                              model=self.name)
        if not self._decisions:
            return ModelReply(decision=None, error="scripted client exhausted",
                              model=self.name)
        payload = self._decisions.pop(0)
        return ModelReply(decision=Decision.from_payload(payload, source="model"),
                          input_tokens=len(evidence) // 4, output_tokens=24,
                          model=self.name)


class AnthropicDecisionClient:
    """Claude API client for the decision step.

    `base_url` is passed through, so an Anthropic-compatible gateway or proxy can
    be used instead of the public API. A local OpenAI-compatible server speaks a
    different protocol and needs its own `DecisionClient` implementation rather
    than a base_url override -- the seam is here precisely so that can be added
    without touching the loop.
    """

    def __init__(self, model: str = DEFAULT_MODEL, effort: str = DEFAULT_EFFORT,
                 base_url: str | None = None, api_key: str | None = None,
                 max_tokens: int = 1024, timeout: float = 60.0,
                 use_fallbacks: bool = True, client=None):
        self.model = model
        self.effort = effort
        self.max_tokens = max_tokens
        self.use_fallbacks = use_fallbacks
        self.name = f"anthropic:{model}"
        self._client = client
        if self._client is None:
            try:
                import anthropic
            except ImportError as exc:  # pragma: no cover - depends on install
                raise RuntimeError(
                    "the anthropic package is not installed; install it or use a "
                    "different DecisionClient") from exc
            kwargs = {"timeout": timeout}
            if base_url or os.environ.get("ANTHROPIC_BASE_URL"):
                kwargs["base_url"] = base_url or os.environ["ANTHROPIC_BASE_URL"]
            if api_key:
                kwargs["api_key"] = api_key
            # No api_key argument means the SDK resolves credentials itself:
            # ANTHROPIC_API_KEY, then ANTHROPIC_AUTH_TOKEN, then an `ant auth
            # login` profile. Do not pre-empt that with a key check here.
            self._client = anthropic.Anthropic(**kwargs)

    def request_kwargs(self, system: str, evidence: str) -> dict:
        """The exact request body. Separated out so tests can assert its shape
        without a network call or a credential."""
        kwargs: dict = {
            "model": self.model,
            "max_tokens": self.max_tokens,
            "system": system,
            "messages": [{"role": "user", "content": evidence}],
            # Structured output: the action is an enum the API enforces, not
            # prose this code has to parse out of a paragraph.
            "output_config": {
                "format": {"type": "json_schema", "schema": ACTION_SCHEMA},
                "effort": self.effort,
            },
        }
        if self.use_fallbacks:
            kwargs["betas"] = [FALLBACK_BETA]
            kwargs["fallbacks"] = "default"
        return kwargs

    def decide(self, system: str, evidence: str) -> ModelReply:
        import anthropic

        kwargs = self.request_kwargs(system, evidence)
        try:
            if self.use_fallbacks:
                response = self._client.beta.messages.create(**kwargs)
            else:
                kwargs.pop("betas", None)
                kwargs.pop("fallbacks", None)
                response = self._client.messages.create(**kwargs)
        except anthropic.APIStatusError as exc:
            return ModelReply(decision=None, model=self.model,
                              error=f"API error {exc.status_code}: {exc.message}")
        except anthropic.APIConnectionError as exc:
            return ModelReply(decision=None, model=self.model,
                              error=f"could not reach the API: {exc}")

        usage = getattr(response, "usage", None)
        reply = ModelReply(
            decision=None,
            input_tokens=getattr(usage, "input_tokens", 0) or 0,
            output_tokens=getattr(usage, "output_tokens", 0) or 0,
            model=getattr(response, "model", self.model),
        )

        # CHECK REFUSAL BEFORE READING CONTENT. On a refusal the content is not
        # the decision, and `stop_details` is populated only in this case.
        if getattr(response, "stop_reason", None) == "refusal":
            details = getattr(response, "stop_details", None)
            reply.refused = True
            reply.refusal_category = getattr(details, "category", "") or ""
            reply.refusal_explanation = getattr(details, "explanation", "") or ""
            return reply

        text = next((b.text for b in response.content if b.type == "text"), "")
        if not text.strip():
            reply.error = "model returned no text content"
            return reply
        try:
            payload = json.loads(text)
        except json.JSONDecodeError as exc:
            reply.error = f"structured output was not valid JSON: {exc}"
            return reply
        try:
            reply.decision = Decision.from_payload(payload, source="model")
        except ParameterRefused as exc:
            # The schema is enforced server-side, so this should not happen --
            # but a decision outside the closed set is refused here too rather
            # than trusted because the API said it was valid.
            reply.error = f"decision rejected by the harness: {exc}"
        return reply
