"""Streamlit front end for mypibypy.  Run with:  streamlit run app.py

Deliberately small. This UI is a form, a table and a download, and the measured
runtimes make Streamlit's rerun model fine -- a single-pair hg38 design is about
2 s end to end and six tiled pairs about 4 s, so there is no job queue, no
progress streaming and no async machinery here. The alternative considered was
the React path (React + MUI + Vite, a FastAPI backend, auth, two repos,
thousands of lines); for this surface it would be more framework than feature.

CACHING AND SEQUENCE DATA. `st.cache_data` is memory-only here on purpose --
`persist="disk"` would write real target sequence into Streamlit's cache
directory, and real sequences stay out of anything that could be committed or
synced. The same goes for the workbook: it is offered as a download from memory
rather than written into the repo.
"""
from __future__ import annotations

import os
from pathlib import Path

import streamlit as st

# `systems`/`loci` search the working directory and its PARENTS, which works for
# the CLI (run from the package directory) but not for a UI that may be launched
# from the repository root. Point them at the files next to this script when
# nothing else is configured, so the app behaves the same either way. The sidebar
# always shows which file was actually loaded, so this can never be ambiguous.
_HERE = Path(__file__).resolve().parent
for _var, _name in (("MYPIBYPY_SYSTEMS", "systems.local.json"),
                    ("MYPIBYPY_LOCI", "loci.local.json")):
    if not os.environ.get(_var) and (_HERE / _name).is_file():
        os.environ[_var] = str(_HERE / _name)

from mypibypy.agent import (
    Action,
    Budget,
    DesignAgent,
    DesignParameters,
    build_client,
)
from mypibypy.agent.intent import (
    CHAT_SYSTEM_PROMPT,
    Intent,
    apply_intent,
    extract_target,
    intent_schema,
)
from mypibypy.agent.tools import filter_candidates, window_stats
from mypibypy.design import DEFAULT_MIN_SPACING_NT, scan
from mypibypy.loci import Locus, load_fusions, load_loci, load_rna_targets
from mypibypy.pipeline import design_fusion, design_locus, design_rna
from mypibypy.sequence import GenomicInterval, fetch_sequence, fetch_transcript
from mypibypy.specificity import Backend, resolve_backend, SpecificityNotConfigured
from mypibypy.systems import is_example, resolve_systems_path, systems
from mypibypy.ui import candidate_rows, unscreened_pairs, workbook_bytes

st.set_page_config(page_title="mypibypy", page_icon="🧬", layout="wide")


# --- cached resources ------------------------------------------------------

@st.cache_resource
def registry() -> dict:
    """The amplifier registry. A resource, not data: it is loaded once and held."""
    return systems()


@st.cache_resource
def screening_backend(index: str) -> Backend | None:
    """A bowtie/BLAST handle. Resource-cached so the index is resolved once."""
    if not index:
        return None
    try:
        return resolve_backend(bowtie_index=index)
    except SpecificityNotConfigured as exc:
        st.warning(f"No screen will run: {exc}")
        return None


@st.cache_data(show_spinner="Fetching genomic sequence…")
def cached_region(build: str, chrom: str, start: int, end: int) -> str:
    """Memory-only. Never persist=\"disk\": that writes target sequence to disk."""
    return fetch_sequence(GenomicInterval(build, chrom, start, end))


@st.cache_data(show_spinner="Fetching transcript…")
def cached_transcript(accession: str) -> str:
    return fetch_transcript(accession)


# --- helpers ---------------------------------------------------------------

def run_design(fn, *args, **kwargs):
    """Call a design function and turn its failures into readable messages.

    A design can legitimately fail -- a region that is entirely repeat-masked has
    no valid window, an accession can be wrong, a fetch can be blocked by a TLS
    proxy. None of those is a bug, and none should reach the user as a
    traceback.
    """
    from mypibypy.design import FusionGeometryError
    from mypibypy.loci import LociNotConfigured
    from mypibypy.sequence import TranscriptFetchError
    try:
        return fn(*args, **kwargs)
    except (ValueError, FusionGeometryError) as exc:
        st.error(f"**No design could be produced.** {exc}")
    except TranscriptFetchError as exc:
        st.error(f"**Could not fetch the transcript.** {exc}")
    except LociNotConfigured as exc:
        st.error(f"**Configuration problem.** {exc}")
    except OSError as exc:
        st.error(f"**Could not reach a reference service.** {exc}")
    return None


def chat_panel(backend, system_key: str, client_kind: str,
               system_keys: list[str]) -> None:
    """Converse with the pipeline.

    The model does exactly two things here: turn a request into parameters, and
    pick which constraint gives when nothing passes. Everything else -- the
    target, the design, the screens, the workbook -- is deterministic.

    THE TARGET IS NEVER THE MODEL'S. It is pulled from the user's own words by
    `extract_target`, a strict pattern in the harness. The model has no field in
    which to supply a region, an accession or a sequence, which is what keeps
    the no-model-authored-sequences guarantee intact once a conversation is
    involved.
    """
    state = st.session_state
    state.setdefault("chat", [])
    state.setdefault("params", DesignParameters(system=system_key))
    state.setdefault("target", None)
    state.setdefault("chat_results", [])

    st.caption(
        "Try: *design six pairs for this enhancer in cy5* · *too few passed, "
        "loosen the Tm* · *give me 12 pairs at 20 nt spacing on "
        "chr7:5530000-5531000*"
    )
    params: DesignParameters = state["params"]
    cols = st.columns(4)
    cols[0].metric("Channel", params.system or "—")
    cols[1].metric("Pairs", params.n_pairs)
    cols[2].metric("Spacing", f"{params.min_spacing} nt")
    cols[3].metric("Target", "set" if state["target"] else "none")

    for entry in state["chat"]:
        with st.chat_message(entry["role"]):
            st.markdown(entry["text"])
            if entry.get("rows"):
                st.dataframe(entry["rows"], hide_index=True,
                             use_container_width=True)

    if state.pop("needs_refresh", False):
        st.rerun()

    message = st.chat_input("Ask for a design, or adjust one…")
    if not message:
        return

    state["chat"].append({"role": "user", "text": message})
    with st.chat_message("user"):
        st.markdown(message)

    client = build_client(client_kind)
    if client is None:
        reply = ("No decision model is selected, so I cannot interpret that. "
                 "Pick **ollama** in the sidebar, or use the Form interface.")
        state["chat"].append({"role": "assistant", "text": reply})
        with st.chat_message("assistant"):
            st.markdown(reply)
        return

    # 1. The target comes from the user's words, by pattern, not from the model.
    found = extract_target(message)
    if found.found:
        state["target"] = found

    # 2. The model's first job: request -> parameters.
    with st.chat_message("assistant"):
        with st.spinner("Reading the request…"):
            payload, error, seconds = client.structured(
                CHAT_SYSTEM_PROMPT, message, intent_schema(tuple(system_keys)))
        if error or payload is None:
            reply = f"I could not read that request. `{error}`"
            st.markdown(reply)
            state["chat"].append({"role": "assistant", "text": reply})
            return
        try:
            intent = Intent.from_payload(payload, tuple(system_keys))
        except Exception as exc:
            reply = f"That reply did not fit the allowed shape and was refused: `{exc}`"
            st.markdown(reply)
            state["chat"].append({"role": "assistant", "text": reply})
            return

        new_params, changes = apply_intent(intent, params, tuple(system_keys))
        state["params"] = new_params

        # The note is the MODEL'S COMMENT. Labelled as such, never treated as a
        # justification -- in trials the chosen action was right while the
        # sentence explaining it was backwards.
        lines = [f"_Model comment ({seconds:.1f}s): {intent.note}_"]
        if changes:
            lines.append("**Applied:** " + "; ".join(changes))

        # Loosening a window that is not being applied changes nothing, and
        # saying "Applied: Tm window 55-62 -> 53.5-63.5" without that caveat
        # reads as though it did something. The pipeline ranks on distance from
        # the window by default rather than rejecting on it.
        if intent.relax_action in (Action.RELAX_TM, Action.RELAX_GC) \
                and not new_params.enforce_window:
            lines.append(
                "ℹ️ Note that the GC/Tm window is **not currently enforced**, so "
                "widening it does not change which windows are available — the "
                "pipeline ranks by distance from it instead of rejecting. Say "
                "*enforce the GC/Tm window* first if you want candidates "
                "actually filtered.")

        rows = None
        if intent.intent == "explain":
            lines.append("Nothing was changed — that read as a question rather "
                         "than a request.")
        elif intent.intent == "unclear":
            lines.append("Nothing was changed. Tell me the target and what you "
                         "want designed.")
        elif not state["target"]:
            lines.append("I have the parameters but no target yet. Give me a "
                         "region like `chr7:5530000-5531000` or an accession "
                         "like `NM_000546.6`.")
        else:
            target = state["target"]
            with st.spinner(f"Designing against {target}…"):
                result = design_for_target(target, new_params, backend)
            if result is not None:
                state["chat_results"] = [result]
                rows = candidate_rows(result)
                n_unscreened = unscreened_pairs(result)
                lines.append(
                    f"**{result.n_pairs} pair(s)** designed against `{target}`.")
                if n_unscreened:
                    lines.append(
                        f"⚠️ **No specificity screen ran** for {n_unscreened} of "
                        f"{result.n_pairs} pair(s) — nothing here establishes "
                        f"single-copy status.")
                stats = window_stats(list(result.candidates), new_params)
                lines.append(f"_Windows: {stats}_")

        reply = "\n\n".join(lines)
        st.markdown(reply)
        if rows:
            st.dataframe(rows, hide_index=True, use_container_width=True)
        state["chat"].append({"role": "assistant", "text": reply, "rows": rows})
        # The metric row at the top rendered before this message was processed,
        # so it shows the previous parameters until the next interaction. Rerun
        # once so what it displays matches what was just applied.
        state["needs_refresh"] = True

    if state["chat_results"]:
        st.download_button(
            "⬇️ Download order workbook (.xlsx)",
            data=workbook_bytes(state["chat_results"], "mypibypy chat session"),
            file_name="mypibypy_order.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )


def design_for_target(target, params: DesignParameters, backend):
    """Run the deterministic design for an extracted target. No model involved."""
    from mypibypy.loci import RnaTarget

    if target.kind == "genomic":
        locus = Locus("target", "GENE", "from chat",
                      GenomicInterval(target.build, target.chrom,
                                      target.start, target.end),
                      params.system)
        return run_design(design_locus, locus, backend=backend,
                          screen=backend is not None, n_pairs=params.n_pairs,
                          min_spacing=params.min_spacing)
    if target.kind == "rna":
        rna = RnaTarget("target", "GENE", target.accession, params.system)
        seq = run_design(cached_transcript, target.accession)
        if seq is None:
            return None
        return run_design(design_rna, rna, backend=backend,
                          screen=backend is not None, transcript_seq=seq,
                          n_pairs=params.n_pairs, min_spacing=params.min_spacing)
    return None


def show_run_record(record) -> None:
    """Why, not just what. The decisions and the evidence behind them."""
    st.subheader("Agent run record")
    cols = st.columns(4)
    cols[0].metric("Turns", record.turns_used)
    cols[1].metric("Relaxations", record.relaxations_used)
    cols[2].metric("Model calls", record.model_calls)
    cols[3].metric("Status", record.status)
    if not record.consulted_model:
        st.caption("Settled deterministically — the model was not consulted. "
                   "That is the expected path for a clear case.")
    for event in record.events:
        icon = {"model_decision": "🤖", "deterministic_decision": "⚙️",
                "parameter_clamp": "📏", "budget": "⏳", "refusal": "🚫",
                "error": "⚠️", "outcome": "🏁"}.get(event.kind, "•")
        st.write(f"{icon} **turn {event.turn}** · `{event.kind}` — {event.summary}")
    with st.expander("Full record (JSON)"):
        st.json(record.to_dict())


# --- sidebar ---------------------------------------------------------------

st.sidebar.title("🧬 mypibypy")
st.sidebar.caption("piFISH split-probe design")

reg = registry()
# `loaded_from()` reports only what is currently cached, and Streamlit can
# re-import the module on a hot reload while `registry()` stays cached -- which
# showed "systems: None" even though the real file was loaded. Resolve the path
# directly so the caption always names the file that would be used.
st.sidebar.caption(f"systems: `{resolve_systems_path()}`")
if is_example():
    st.sidebar.error(
        f"Amplifier systems are the SYNTHETIC EXAMPLES "
        f"({resolve_systems_path().name}). "
        "Probes designed now are **not orderable** — see the README on "
        "systems.local.json."
    )

interface = st.sidebar.radio("Interface", ["Form", "Chat"], horizontal=True)
mode = st.sidebar.radio("Target type", ["Genomic region", "RNA transcript",
                                        "Fusion junction"],
                        disabled=interface == "Chat",
                        help="In Chat the target comes from what you type.")
system_key = st.sidebar.selectbox("Amplifier system", sorted(reg),
                                  format_func=lambda k: f"{k} ({reg[k].fluorophore})")

if mode == "Fusion junction":
    st.sidebar.caption("A junction is one position — the gap sits on the "
                       "breakpoint — so tiling does not apply.")
    n_pairs, min_spacing = 1, DEFAULT_MIN_SPACING_NT
else:
    n_pairs = st.sidebar.number_input("Pairs to tile", 1, 50, 1)
    min_spacing = st.sidebar.number_input("Minimum spacing (nt)", 0, 5000,
                                         DEFAULT_MIN_SPACING_NT)

index_path = st.sidebar.text_input(
    "bowtie index (optional)",
    help="Leave blank to skip the genome/transcriptome screen. The workbook "
         "will then say in writing that no screen was run.")
use_agent = st.sidebar.checkbox(
    "Run through the agent loop", value=False, disabled=interface == "Chat",
    help="Adds the deterministic classifier and, for genuinely conflicting "
         "signals, a model decision. Clear cases are settled without a model.")
client_kind = st.sidebar.selectbox(
    "Decision model", ["ollama", "none", "anthropic"], index=0,
    help="Local model by default. 'none' runs the pipeline with no model at "
         "all, flagging anything that would need a judgement call.")

# --- main ------------------------------------------------------------------

st.title("Design piFISH probes")
backend = screening_backend(index_path.strip())
results, note = [], ""

if interface == "Chat":
    chat_panel(backend, system_key, client_kind, sorted(reg))
    st.stop()

if mode == "Genomic region":
    col1, col2 = st.columns([2, 3])
    build = col1.text_input("Build", "hg38")
    region = col2.text_input("Region (chrom:start-end)", "",
                             placeholder="chr5:142906255-142906754")
    gene = col1.text_input("Gene label", "GENE")
    name = col2.text_input("Target name", "target")
    if st.button("Design", type="primary", disabled=not region):
        try:
            chrom, span = region.split(":")
            start, end = (int(x.replace(",", "")) for x in span.split("-"))
        except ValueError:
            st.error("Region must look like chr5:142906255-142906754")
            st.stop()
        locus = Locus(name, gene, "designed in the UI",
                      GenomicInterval(build, chrom, start, end), system_key)
        with st.spinner("Designing…"):
            one = run_design(design_locus, locus, backend=backend,
                             screen=backend is not None, n_pairs=int(n_pairs),
                             min_spacing=int(min_spacing))
        if one is not None:
            results = [one]
            note = f"mypibypy UI: {locus.interval}, system {system_key}."

elif mode == "RNA transcript":
    accession = st.text_input("Accession", "", placeholder="NM_000546.6 or ENST…")
    gene = st.text_input("Gene symbol", "", placeholder="used to exclude on-target hits")
    st.caption("Screened against a **transcriptome** index — see "
               "`pipeline.RNA_SCREEN_POLICY` for why, and why a junction-spanning "
               "footprint does not map contiguously to a genome.")
    if st.button("Design", type="primary", disabled=not (accession and gene)):
        from mypibypy.loci import RnaTarget
        target = RnaTarget("target", gene, accession.strip(), system_key)
        with st.spinner("Fetching and designing…"):
            seq = run_design(cached_transcript, accession.strip())
            one = run_design(design_rna, target, backend=backend,
                             screen=backend is not None, transcript_seq=seq,
                             n_pairs=int(n_pairs),
                             min_spacing=int(min_spacing)) if seq else None
        if one is not None:
            results = [one]
            note = f"mypibypy UI: RNA target {accession}, system {system_key}."

else:
    try:
        fusions = load_fusions()
    except Exception as exc:  # a missing or malformed loci file
        st.error(f"Could not load fusions: {exc}")
        fusions = []
    if not fusions:
        st.info("No fusion junctions are defined. Add a `fusions` entry to your "
                "loci file — see `loci.example.json`.")
    else:
        chosen = st.selectbox("Junction", fusions, format_func=str)
        if st.button("Design", type="primary"):
            with st.spinner("Designing both strands…"):
                both = run_design(design_fusion, chosen, backend=backend,
                                  screen=backend is not None)
            if both:
                results = both
                note = f"mypibypy UI: fusion junction {chosen.name}."

# --- results ---------------------------------------------------------------

if results:
    if use_agent:
        primary = results[0]
        candidates = scan(primary.region_seq, primary.best.interval,
                          top=len(primary.region_seq)) if mode == "Genomic region" \
            else list(primary.candidates)

        def design_fn(params: DesignParameters):
            return primary, candidates

        agent = DesignAgent(design_fn, client=build_client(client_kind),
                            budget=Budget(max_turns=2),
                            known_systems=tuple(sorted(reg)))
        run = agent.run(DesignParameters(system=system_key, n_pairs=int(n_pairs),
                                         min_spacing=int(min_spacing)),
                        target_label=note or "UI run", request="designed from the UI")
        st.info(f"**Agent: {run.status}** — {run.message}")
        show_run_record(run.record)
        st.divider()

    for result in results:
        st.subheader(f"{result.locus.gene} {result.locus.name} → "
                     f"{result.system.fluorophore}")
        n_unscreened = unscreened_pairs(result)
        if n_unscreened:
            st.warning(
                f"**No genome-wide specificity screen was run** for "
                f"{n_unscreened} of {result.n_pairs} pair(s). Nothing here "
                f"establishes single-copy status. Supply a bowtie index and "
                f"re-run before ordering."
            )
        st.dataframe(candidate_rows(result), hide_index=True,
                     use_container_width=True)
        stats = window_stats(list(result.candidates),
                             DesignParameters(system=system_key))
        st.caption(f"Windows scored: {stats}. The GC/Tm window is reported for "
                   f"context — the pipeline ranks on distance from it rather "
                   f"than rejecting on it.")

    st.download_button(
        "⬇️ Download order workbook (.xlsx)",
        data=workbook_bytes(results, note),
        file_name="mypibypy_order.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        type="primary",
    )
    st.caption("Probe sequences are in the workbook you download. They are not "
               "written anywhere inside the repository.")
else:
    st.caption("Pick a target and press **Design**.")
