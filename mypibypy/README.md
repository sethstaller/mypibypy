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
pytest                    # 41 tests, no PYTHONPATH needed
```

## Usage

```bash
# every locus in a loci file (see "Configuring target loci" below)
python -m mypibypy batch --loci loci.local.json -o order.xlsx

# one locus straight from the command line
python -m mypibypy design --build mm10 --region chr5:142906255-142906754 \
    --system cy5 --gene Actb --name promoter -o order.xlsx
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
| `specificity.py` | genome-wide screening (NCBI web BLAST; being replaced by a local aligner) |
| `crosshyb.py`    | cross-hybridization against the amplifier pool and the order     |
| `thermo.py`      | buffer-aware Tm and wash-survival margins                        |
| `workbook.py`    | render the order workbook (ORDER_ssDNA / Parts list / per-locus) |
| `cli.py`         | `python -m mypibypy`                                             |

## QC / specificity

Window scoring balances GC (per-pad and matched between pads), matched Tm, and
rejects homopolymer / low-complexity runs. RepeatMasker overlap is reported.
Genome-wide BLAT/BLAST is the recommended final check (`specificity.py` provides
an NCBI fallback since UCSC web BLAT now requires a captcha).

## Credit and license

Developed by Seth Staller, St. Jude Children's Research Hospital.

Released under the MIT License. Copyright is held by St. Jude Children's
Research Hospital; see LICENSE. To cite this software, see CITATION.cff.
