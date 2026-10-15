# mypibypy

Design **piFISH** split-probe primary probes (P1L / P1R) for genomic targets and
emit an order workbook that drops straight into an existing, validated
P2 / P3 / P4 amplifier pool.

## The piFISH cascade

```
target  ->  P1L + P1R (split primary)  ->  P2 (secondary)  ->  P3 (tertiary)  ->  P4 (Cy signal)
```

- `P1L = arm_left + revcomp(target_left25)`
- `P1R = revcomp(target_right25) + arm_right`
- Target footprint = **25 nt + 1 nt gap + 25 nt = 51 nt** (plus strand).
- The two constant arms meet at the 1-nt junction and together form the single
  landing site for P2's central anneal region — so P2 only lands when **both**
  primaries are bound (coincidence detection = specificity).
- P2 carries ~10 P3 sites; P3 carries ~10 P4 sites → ~100× signal per pair.

A **system** bundles one channel: the two arms + P2 + P3 + the fluorophore-labeled
P4. Designing a new target reuses these systems, so the **only new oligos are the
P1L/P1R primaries**.

## Configuring amplifier systems

**No real amplifier sequences are committed to this repository.** They are
lab-specific and frequently unpublished, so they load at runtime from a file you
keep outside version control. Resolution order:

1. `$MYPIBYPY_SYSTEMS` — explicit path to a systems JSON file
2. `systems.local.json` — in the working directory or any parent (gitignored)
3. `src/mypibypy/systems.example.json` — **synthetic placeholders**, so a fresh
   clone and the test suite work out of the box. They will not produce working
   probes.

To set up your own:

```bash
cp src/mypibypy/systems.example.json systems.local.json
# then replace the placeholder sequences with your real arms / P2 / P3 / P4
```

The example file documents the required architecture, which `crosshyb` and the
test suite both enforce:

```
arm_left  = [ top_L: 14 nt ][ mid_L: 8 nt ]
arm_right = [ mid_R: 8 nt ][ top_R: 14 nt ]
P2 anneal region (29 nt) = rc(top_R) + 1 nt + rc(top_L)
P2 repeat (30 nt) = spacer(10) + P3-binding(20)
P3 anneal (20 nt) = rc(P2's P3-binding region)
P3 repeat (30 nt) = spacer(10) + rc(P4)
mid_L / mid_R pair over 2–4 bp — the "pi bond" holding the pair together
```

Check `mypibypy.is_example()` before treating any output as orderable.

Order workbooks (`*.xlsx`) are gitignored for the same reason: they carry probe
sequences and target/fusion designs. Keep them in a local `designs/` directory.

## Setup

The design core is pure Python, but specificity screening shells out to BLAST
and bowtie, which are not pip-installable. Create the environment with conda:

```bash
conda create -n pifish-tools -c conda-forge -c bioconda \
    python=3.11 blast bowtie bowtie2 samtools viennarna openpyxl pytest
conda activate pifish-tools
pip install -e .          # from this directory
```

Then configure your amplifier systems (see below) and run the tests:

```bash
pytest                    # 273 tests, no PYTHONPATH needed
```

## Usage

```bash
# every locus in a loci file (see "Configuring target loci" below)
python -m mypibypy batch --loci loci.local.json -o order.xlsx

# one locus straight from the command line
python -m mypibypy design --build mm10 --region chr5:142906255-142906754 \
    --system cy5 --gene Actb --name promoter -o order.xlsx

# fusion junctions: sense + antisense pairs, gap placed on the breakpoint
python -m mypibypy fusion --loci loci.local.json -o fusion_order.xlsx

# with a genome screen (otherwise the workbook says no screen was run)
python -m mypibypy batch --loci loci.local.json -o order.xlsx \
    --bowtie-index /data/genomes/mm10
```

```python
from mypibypy.loci import load_loci
from mypibypy.pipeline import design_locus

for locus in load_loci("loci.local.json"):
    res = design_locus(locus)
    print(locus.gene, locus.name, res.pair.p1l, res.pair.p1r)
```

## Configuring target loci

Target loci are **inputs, not source**: which regions you probe is
experiment-specific, so they load from a JSON file rather than being hardcoded.
Resolution order is `--loci` → `$MYPIBYPY_LOCI` → `loci.local.json` (searched
upward, gitignored) → `loci.example.json` (the published mouse *Actb* positive
control, so a fresh clone runs).

```bash
cp src/mypibypy/loci.example.json loci.local.json   # then edit
```

Both `batch` and `design` print a warning, and write one into the workbook, when
example systems or example loci are in use -- those outputs are not orderable.

## Layout

| module           | role                                                            |
|------------------|-----------------------------------------------------------------|
| `sequence.py`    | UCSC sequence fetch, revcomp, GC, nearest-neighbor Tm, motifs    |
| `systems.py`     | amplifier-system registry (arms + P2/P3/P4 + fluorophore)        |
| `design.py`      | 51-nt window scan/scoring + P1L/P1R construction                 |
| `loci.py`        | target loci loaded from JSON (`Locus`, `load_loci`)              |
| `pipeline.py`    | locus → fetch → scan → chosen pair                               |
| `specificity.py` | pair-aware genome screening via local bowtie v1 / blastn         |
| `crosshyb.py`    | cross-hybridization against the amplifier pool and the order     |
| `thermo.py`      | buffer-aware Tm and wash-survival margins                        |
| `accessibility.py` | RNA target-site secondary structure (RNAplfold)                |
| `outcomes.py`    | design → imaging-result record, for calibrating thresholds later |
| `workbook.py`    | render the order workbook (ORDER_ssDNA / Parts list / per-locus) |
| `cli.py`         | `python -m mypibypy`                                             |

## QC and specificity

Window scoring balances GC (per-pad and matched between pads), matched Tm, and
rejects homopolymer / low-complexity runs. RepeatMasker overlap is reported.

### Pair-aware genome screening

A piFISH pair only makes a spot when **both** halves bind adjacent target, so a
genome hit for one half is not an off-target. `specificity.py` screens each half
separately and then joins the hits: a site is reported only where both halves
land on the same chromosome within the signal window (50 nt, see below), on the
same strand and in an order that could rebuild P2's landing site. A lone-half
hit lowers the score
rather than rejecting the pair. Each site is then passed through `thermo` to turn
"the aligner found it" into "it would still be bound after the washes".

Point it at a local database:

```bash
export MYPIBYPY_BOWTIE_INDEX=/data/genomes/mm10     # bowtie v1 index basename
export MYPIBYPY_BLAST_DB=/data/genomes/mm10         # or a makeblastdb basename
```

With neither set, the screen returns a `NOT_SCREENED` status with a reason, and
the workbook prints that verbatim. **A screen that did not run never renders as
a screen that passed** — check `SpecificityReport.screened` before making any
specificity claim.

Backends are measured, not assumed — rerun
`scripts/benchmark_specificity.py`, which plants sites with known mismatch
counts at known coordinates and counts recovery (40 × 25-nt queries per cell):

| backend | 0mm | 1mm | 2mm | 3mm | 4mm | hits kept |
|---|---|---|---|---|---|---|
| `bowtie -v 3 -a` | 40/40 | 40/40 | 40/40 | 40/40 | 0/40 | 160 |
| `blastn-short -word_size 5` | 40/40 | 40/40 | 40/40 | 40/40 | 18/40 | 1316 |
| `blastn-short -word_size 7` | 40/40 | 40/40 | 40/40 | 39/40 | 16/40 | 765 |
| `blastn-short -word_size 11` | 40/40 | 40/40 | 27/40 | 14/40 | 5/40 | 151 |
| `bowtie2 --very-sensitive-local` | 40/40 | 0/40 | 0/40 | 0/40 | 0/40 | 40 |

bowtie v1 is preferred: `-v 3 -a` is exhaustive within its budget and returns
only the planted sites. `bowtie2 --very-sensitive-local`, the obvious first
choice, is unusable at 25 nt — a single mismatch often leaves no seed — and is
kept only as a baseline. blastn needs `-dust no`, since DUST masks short queries
into silence, which looks exactly like a clean result.

## Fusion junctions

`mypibypy fusion` places the 1-nt gap **on the breakpoint**, so the left pad lies
wholly in the 5' partner and the right pad wholly in the 3' partner, and neither
parent gene alone can bring both halves together. `check_partner_specificity`
enforces that each pad is present in its own partner and absent from the other
(both strands). Sense and antisense pairs are emitted together, on two channels.

`design.arm_offsets` / `arm_separation` answer a question no per-probe check can:
where are the two constant arms actually anchored? A pair with its pads on the
other arms has the right lengths, the right arms and the right pads, and differs
only in that its arms sit at the footprint's outer ends rather than either side
of the gap. `--pad-assignment` chooses between the two, defaulting to `direct`
(arms beside the gap, as in the sense pair).

Both conventions give signal at the bench. An earlier version of this tool
claimed the outer-end geometry could not work, reasoning that arms 50 nt apart
could never present one 29-nt P2 site — that treated the arms as rigid rods when
they are flexible single-stranded tethers (22 nt is ~14 nm of contour against the
~17 nm spanned by 50 bp of duplex), and the lab's hand-made antisense pairs are
built that way and do produce signal. `arms_meet_at_junction` reports the
geometry only; it is not a verdict on whether a pair works, and the relative
efficiency of the two has never been measured.

That correction propagates into screening: if a *designed* pair bridges 50 nt,
an accidental one can too, so `specificity.DEFAULT_SIGNAL_GAP_NT` is 50 nt rather
than the 10 nt originally used. It is an observed lower bound, not a measured
limit.

## Screening a transcriptome: gene-level collapsing

A gene appears once in a genome and ~8.5 times in GENCODE. Screened against a
transcriptome without gene awareness, a probe finds its own target in every one
of its isoforms, the pair logic reads those as off-target sites, and a good probe
is rejected. Measured on real GENCODE v50 data, one footprint produced **70 raw
pair sites which collapsed to 1, and that single "off-target gene" was the
intended gene itself — verdict `fail`.** Naming the target gene excluded 140
on-target hits and the same probe passed.

```python
screen_pair(left25, right25, backend=transcriptome, target_gene="MYGENE")
```

- The gene is parsed out of the reference name (`gencode_gene`, field 6 of a
  pipe-delimited GENCODE header), and the parser is **detected from the names
  that actually came back**. A plain chromosome FASTA has no gene field, so
  `gene_aware` is `False` and a genome screen behaves exactly as before —
  guarded by `TestGenomeScreensAreUnchanged`.
- Site count and gene count are reported separately: 14 sites across 8 isoforms
  of 2 genes is 2 problems, not 14. Read `offtarget_genes`.
- Collapsing runs *after* the pair geometry, not before. The gap between the two
  halves decides whether a site can signal, and that gap belongs to one
  transcript's coordinates — two isoforms can place the same pads 1 nt apart in
  one and 800 nt apart in the other, because an exon is spliced out of only one.
  Collapsing first would have to pick a representative coordinate and would
  silently invent or destroy sites.

## RNA piFISH

RNA targets are identified by accession, not coordinates, because a transcript is
a molecule rather than an interval:

```bash
python -m mypibypy rna --loci loci.local.json -o rna_order.xlsx \
    --transcriptome-index /path/to/gencode/transcripts
```

`sequence.fetch_transcript` dispatches on the accession prefix — Ensembl REST for
`ENST`/`ENSMUST`, NCBI E-utilities for `NM_`/`NR_`/`XM_`/`XR_` — and validates
what comes back, because a truncated response that still parses as sequence would
silently shorten the transcript and shift every downstream offset:

- Ensembl is asked for JSON rather than FASTA, so a partial body fails to parse
  instead of becoming a short sequence. It needs an `Accept` header; sending
  `Content-Type` returns 200 with an empty body.
- RefSeq length is fetched independently via `esummary` and compared against the
  downloaded sequence, since efetch FASTA carries no length.

**Which reference an RNA design is screened against — transcriptome only, by
default.** The transcriptome is the population of molecules actually present, and
a footprint may legitimately span an exon-exon junction, which is often the point
since it is what makes a probe isoform-specific. Such a footprint does not map
contiguously to the genome at all: screened against one, its two halves land in
different exons thousands of nt apart, the real on-target site is invisible, and
the output is actively misleading. Genomic DNA is in the sample too, but a
standard RNA-FISH protocol has no denaturation step, so that DNA stays
double-stranded and a single-stranded probe cannot invade it. If your protocol
*does* denature, pass `--also-genome-index`. See `pipeline.RNA_SCREEN_POLICY`.

Note the co-localisation window is in **transcript** coordinates when screening a
transcriptome: two half-sites 50 nt apart in a spliced mRNA may be far apart in
the genome, and the spliced distance is the one that matters.

`exon_boundaries` on an RNA target are optional 0-based transcript offsets. They
are used only to *report* that a footprint spans a junction, and where —
`junction_position` distinguishes a boundary landing on the 1-nt gap, where each
pad lies wholly within one exon, from one landing mid-pad.

### Target-site accessibility

A probe can be perfectly specific and still never bind, because its site is
folded into a stem. This failure mode does not exist for DNA targets, which the
protocol denatures. `accessibility.py` scores the **target site on the
transcript** in the context of its flanks — not the folding of the probe oligo,
which is a different question.

The metric is the probability of an unpaired run long enough to nucleate
(default 8 nt), not the mean per-base unpaired probability. That choice is
measured, not stylistic — on a synthetic hairpin:

| region | P(base unpaired) | P(8-nt run unpaired) |
|---|---|---|
| open flank | 0.95 | 0.86 |
| stem | 0.28 | 0.046 |
| hairpin loop | 0.89 | **0.008** |

A mean would call that loop accessible; a probe cannot nucleate in it.

- The two pads are scored **separately** and the pair takes the **worse** of
  them, never the average — both halves must bind, so one buried pad kills the
  design however open the other is.
- Folded with `RNAplfold -W 80 -L 40` and 150 nt of flanking transcript either
  side. 150 is comfortably more than both `W` and `L`: every footprint base sees
  a full window, and nothing beyond `L` can pair with it at all. Where a site
  sits near a transcript end, the flank actually available is reported rather
  than padded over.
- It is a **ranking input and a reported attribute, never a hard reject**, and
  it does not override gross GC or Tm problems. Applied in the RNA path only, so
  `design._score` stays comparable with designs already ordered — installing
  ViennaRNA is the only thing that changes an ordering.
- **No calibrated threshold exists.** Nothing ties a piFISH failure to a measured
  accessibility, and the caveat travels with every number into the workbook.

Needs ViennaRNA (`conda install -c conda-forge -c bioconda viennarna`). Without
it, results are `NOT_EVALUATED` with a reason and contribute zero to ranking.

### Behind a TLS-inspecting proxy

If `curl` reaches these APIs but Python raises
`CERTIFICATE_VERIFY_FAILED`, the proxy's root is in the macOS keychain but not in
this interpreter's CA bundle. Export the system roots rather than disabling
verification:

```bash
security find-certificate -a -p /System/Library/Keychains/SystemRootCertificates.keychain > roots.pem
security find-certificate -a -p /Library/Keychains/System.keychain >> roots.pem
export SSL_CERT_FILE=$PWD/roots.pem
```

The error message says this too. The same applies to the UCSC fetch used by
genomic designs.

## Tiling: several pairs per target

`design.scan` always ranked many windows, but only the best became a pair. One
pair is defensible — at ~100x amplification per pair, and every validated design
in the lab is single-pair — so **tiling is headroom for marginal signal, not a
fix for a defect.** For scale: the pi-FISH paper used 35–50 probes per DNA locus
per channel, and for RNA roughly 10–15 pairs per transcript at medium-to-low
abundance, 1–5 for a highly abundant one.

```bash
python -m mypibypy batch --loci loci.local.json -o order.xlsx --pairs 6
python -m mypibypy rna   --loci loci.local.json -o rna.xlsx   --pairs 10 \
    --transcriptome-index /path/to/gencode/transcripts
```

- Selection is greedy over the existing ranking: walk candidates best-first and
  keep each that clears every chosen footprint by `--min-spacing` (default 3 nt)
  on top of not overlapping. Greedy rather than optimal on purpose — maximising
  the *count* of packed footprints would start preferring mediocre windows, and
  the reason to tile is more signal from good pairs.
- **Every selected pair is screened**, and for RNA has its accessibility
  evaluated, separately. Screening only the best and ordering ten would be worse
  than not tiling, since the unchecked nine look just as validated on the sheet.
  In a real run one of five tiled pairs came back with 17 off-target half-hits
  where the others had none.
- Asking for more pairs than fit returns what fits and warns, rather than
  failing.
- **Fusion junctions refuse `--pairs > 1`** rather than ignoring it. A junction
  is a single position — the design *is* the 1-nt gap on the breakpoint — and a
  second pair could only lie wholly inside one partner, where it reports that
  gene unfused, which is the one thing a fusion probe exists to rule out.
  (Sense and antisense are a separate axis and still both emitted.)
- Names gain a `_t01` suffix only when tiling. **A single-pair design is
  byte-identical to before tiling existed** — same names, same form-style locus
  sheet, same notes — verified by diffing a regenerated workbook against a
  pre-change snapshot, not assumed. A tiled locus sheet becomes a table with one
  row per pair, and the REUSED amplifier count does not grow with pair count,
  which is the whole point of tiling into an existing pool.

## Recording what happened

`outcomes.py` keeps an append-only JSONL record linking each ordered pair
(sequences, locus, system, design-time QC, screen verdicts) to its imaging result
from [zFISHer](https://github.com/sstaller-sj/zFISHer), read from that project's
`reports/<layer>_puncta.csv`. It computes no score and fits no threshold: with
three designs there is nothing to fit, and this exists so the comparison becomes
possible later. The store holds real sequences, so it is gitignored
(`outcomes.local.jsonl`).

## Credit and license

Developed by Seth Staller, St. Jude Children's Research Hospital.

Released under the MIT License. Copyright is held by St. Jude Children's
Research Hospital; see LICENSE. To cite this software, see CITATION.cff.
