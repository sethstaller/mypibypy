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
    python=3.11 blast bowtie bowtie2 samtools openpyxl pytest
conda activate pifish-tools
pip install -e .          # from this directory
```

Then configure your amplifier systems (see below) and run the tests:

```bash
pytest                    # 158 tests, no PYTHONPATH needed
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
