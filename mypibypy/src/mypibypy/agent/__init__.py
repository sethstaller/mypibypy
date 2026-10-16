"""Agent layer: a harness where the model's output decides the next step.

The terms, in the lab's own framing:

    MODEL    stateless text in, text out.
    HARNESS  the code around it -- schema, state, routing, retries, and the
             calls into deterministic tools.
    AGENT    a harness where the model's output decides what happens next.

This package is the agent. The pipeline it drives was built first and works
without it; the model is reserved for the cases where several signals conflict
and no rule settles the question. Thresholds come from the lab as hardcoded
rules, because a language model has no innate sense of what "good" means for
this assay.

What is load-bearing here, in one place:

* The model never writes or edits a sequence. No tool accepts one, and its
  output schema has no field that could carry one. See `tools` and `actions`.
* The action set is closed and enforced by structured output, not by parsing.
* Parameter bounds are code, not prompt text. See `actions.clamp`.
* Every tool call, clamp and decision lands in a durable record. See `record`.
* Turn, relaxation and token budgets are hard. Three relaxations, then the run
  reports failure instead of loosening until something passes.
"""
from .actions import (
    ACTION_SCHEMA,
    BOUNDS,
    RELAXING_ACTIONS,
    Action,
    Decision,
    DesignParameters,
    ParameterRefused,
    apply_action,
    clamp,
)
from .client import AnthropicDecisionClient, DecisionClient, ModelReply, ScriptedClient
from .loop import AgentResult, Budget, DesignAgent, Evaluation, Situation, classify
from .record import RunRecord, RunStore
from .tools import ScreenEvidence, WindowStats, screen_evidence, window_stats

__all__ = [
    "ACTION_SCHEMA", "BOUNDS", "RELAXING_ACTIONS", "Action", "AgentResult",
    "AnthropicDecisionClient", "Budget", "Decision", "DecisionClient",
    "DesignAgent", "DesignParameters", "Evaluation", "ModelReply",
    "ParameterRefused", "RunRecord", "RunStore", "ScreenEvidence",
    "ScriptedClient", "Situation", "WindowStats", "apply_action", "clamp",
    "classify", "screen_evidence", "window_stats",
]
