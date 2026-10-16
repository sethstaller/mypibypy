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

from mypibypy.agent import Budget, DesignAgent, DesignParameters, ScriptedClient
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

mode = st.sidebar.radio("Target type", ["Genomic region", "RNA transcript",
                                        "Fusion junction"])
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
    "Run through the agent loop", value=False,
    help="Adds the deterministic classifier and, for genuinely conflicting "
         "signals, a model decision. Clear cases are settled without a model.")

# --- main ------------------------------------------------------------------

st.title("Design piFISH probes")
backend = screening_backend(index_path.strip())
results, note = [], ""

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

        agent = DesignAgent(design_fn, client=ScriptedClient(),
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
