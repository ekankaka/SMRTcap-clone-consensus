# SMRTcap clone consensus pipeline

Version **1.0.0**

`SMRTcap-clone-consensus` generates one HIV consensus sequence for each pre-assigned clone from randomly sheared HIV fragments.

The pipeline is intentionally independent of any upstream integration-site or non-flanked-read workflow. It does **not** decide which reads belong to a clone. Clone membership is supplied in the input CSV through `participant_id` and `clone_id`.

The same pipeline can therefore be used for:

- flanked reads only;
- flanked reads plus successfully assigned non-flanked reads; or
- any other HIV fragments that already have a clone assignment.

## 1. Required software

- Bash
- Python 3
- MAFFT
- minimap2
- Biopython

Install the Python dependency with:

```bash
python3 -m pip install -r requirements.txt
```

MAFFT and minimap2 must be available on `PATH`.

## 2. Input CSV

The input CSV must contain these six columns:

| Column | Meaning |
| --- | --- |
| `participant_id` | Participant to whom the fragment belongs |
| `sample_id` | Sample/visit containing the fragment |
| `read` | Read identifier |
| `clone_id` | Pre-assigned clone identifier |
| `STRAND` | `plus` or `minus` |
| `HIV_SEQ` | HIV fragment sequence |

Example:

```csv
participant_id,sample_id,read,clone_id,STRAND,HIV_SEQ
101,101-1,read001,chr3_123456,plus,ACTG...
101,101-1,read002,chr3_123456,minus,TTGA...
101,101-2,read003,chr3_123456,plus,GCTA...
```

A biological clone is defined by:

```text
participant_id + clone_id
```

`sample_id` does not define the clone. Fragments from several samples belonging to the same participant can therefore contribute to the same clone consensus.

The combination `participant_id + sample_id + read` must be unique. Identifier fields cannot be blank or contain whitespace, commas, semicolons, or `|`.

## 3. Reference FASTA

Provide one full-length HIV reference sequence, for example HXB2.

The reference is used **only as an alignment scaffold** to place randomly sheared fragments in the correct approximate genome position. It never contributes nucleotides to the clone consensus.

The reference FASTA must contain exactly one sequence.

## 4. Run the pipeline

```bash
bash run_clone_consensus.sh input.csv HXB2.fasta my_output
```

The third argument is the analysis output folder. The pipeline creates two subfolders inside it:

```text
my_output/
├── work/
└── results/
```

`work/` contains intermediate and QC files. `results/` contains the final combined outputs.

## 5. Consensus method

### A. Forward normalization

The pipeline normalizes orientation directly from the input CSV:

```text
STRAND = plus   -> keep HIV_SEQ as supplied
STRAND = minus  -> reverse-complement HIV_SEQ
```

The normalized sequences are saved in:

```text
OUTPUT_DIR/work/normalized_sequences.csv
```

### B. Group fragments into clones

Fragments are grouped by `participant_id + clone_id`.

Each clone receives an internal work identifier such as `clone_000001`. The mapping between this internal identifier and the original participant/clone IDs is retained in `OUTPUT_DIR/work/clone_manifest.csv`.

### C. Reference-scaffolded fragment alignment

For each clone, MAFFT aligns the forward-normalized fragments to the supplied reference using:

```bash
mafft --addfragments clone_fragments.fasta reference.fasta
```

`--keeplength` is deliberately **not** used, so insertions present in clone fragments can be retained.

### D. Coverage-aware majority consensus

The reference row is completely excluded when consensus bases are called.

For each alignment position:

- terminal gaps in a fragment mean **no coverage** and do not vote;
- internal gaps represent deletion evidence and do vote;
- `A`, `C`, `G`, `T`, and internal gaps are counted;
- ambiguous bases such as `N` do not vote;
- a base/deletion is called only when it has **more than 50%** of informative observations;
- a tie or lack of majority is called `N`;
- positions with no fragment coverage are called `N` when they lie between covered regions.

Leading and trailing regions with no fragment coverage are removed. Internal uncovered regions remain as `N`, so the pipeline does not invent sequence between non-overlapping fragments.

If a clone has only one fragment, that fragment still produces a sequence, but the summary labels it `single_fragment` rather than `multi_fragment`.

### E. Fragment-to-consensus QC

Every original forward-normalized fragment is remapped to its clone consensus with minimap2 `map-hifi`.

This remapping is **QC only**. The pipeline reports identity, fragment coverage, alignment length, mapping strand, and MAPQ. It does not automatically remove low-quality or discordant fragments.

## 6. Generated folders

```text
OUTPUT_DIR/
├── work/
│   ├── normalized_sequences.csv
│   ├── clone_manifest.csv
│   ├── reference.fasta
│   └── clones/
│       └── clone_000001/
│           ├── clone_info.json
│           ├── fragments.fasta
│           ├── alignment.fasta
│           ├── consensus.fasta
│           ├── consensus_stats.json
│           ├── depth.csv
│           └── fragments_vs_consensus.paf
└── results/
    ├── clone_consensus.fasta
    ├── clone_consensus_summary.csv
    └── fragment_qc.csv
```

`OUTPUT_DIR/work/` contains intermediate and audit/QC files. `OUTPUT_DIR/results/` contains only the three combined outputs intended for downstream use.

## 7. Main outputs

### `clone_consensus.fasta`

One consensus sequence per `participant_id + clone_id`.

FASTA headers use:

```text
>participant_id|clone_id|n=<number_of_fragments>
```

### `clone_consensus_summary.csv`

One row per clone, including:

- participant and clone identifiers;
- number of fragments and samples;
- sample IDs;
- shortest and longest fragment;
- consensus length and non-`N` length;
- percentage `N`;
- mean, median, minimum and maximum depth;
- number of alignment positions showing disagreement;
- number of fragments remapped to the consensus; and
- remapping identity and coverage summaries.

### `fragment_qc.csv`

One row per input fragment, including:

- participant, sample, read and clone identifiers;
- fragment length;
- whether it remapped to the clone consensus;
- mapping strand;
- alignment length;
- percent identity;
- percent of the fragment covered; and
- MAPQ.

## 8. Important interpretation points

- The pipeline accepts clone membership as given; it does not infer or validate the biological integration site.
- The reference guides alignment only and never supplies consensus bases.
- Internal regions without fragment coverage are represented by `N` rather than filled from the reference.
- Insertions relative to the reference are retained when supported by the fragment alignment.
- A single-fragment sequence is reported but should not be interpreted as having the same support as a consensus reconstructed from several overlapping fragments.
- Remapping metrics are reported for QC rather than used as automatic exclusion thresholds.

## 9. Example data

`example_data/input.csv` and `example_data/reference.fasta` contain small synthetic data for testing the workflow after MAFFT and minimap2 are installed.

Run:

```bash
bash run_clone_consensus.sh example_data/input.csv example_data/reference.fasta example_output
```

The example data are synthetic and are not intended for biological interpretation.

## 10. Files in this repository

```text
run_clone_consensus.sh   Main pipeline runner
clone_consensus.py       Validation, normalization, consensus, and QC logic
requirements.txt         Python dependency
example_data/            Synthetic test data
LICENSE                  MIT license
```

## 11. Reproducibility

For each analysis, record:

- the pipeline version;
- the input CSV used;
- the reference FASTA used; and
- the versions of MAFFT and minimap2.

Current release: **v1.0.0**.
