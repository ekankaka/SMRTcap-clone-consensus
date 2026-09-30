# SMRTcap clone consensus pipeline

Version **1.0.0**

`SMRTcap-clone-consensus` generates one HIV consensus sequence for each pre-assigned clone from randomly sheared HIV fragments.

It is designed as an **independent downstream analysis of `nf-viral-integration`**. It does not depend on `SMRTcap-nonflanked-HIV-analysis` and does not call or import that pipeline.

The repository keeps the workflow deliberately small: one short adapter prepares standard `nf-viral-integration` output, one core Python helper handles normalization/consensus/QC, and the Bash runners keep the MAFFT and minimap2 steps visible.

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

## 2. Recommended: run directly from `nf-viral-integration`

Use the `final_results/` folder produced by `nf-viral-integration`:

```bash
bash run_from_nf.sh FINAL_RESULTS_DIR REFERENCE_FASTA [SAMPLE_MAP_CSV] [OUTPUT_DIR]
```

Example without a sample mapping file:

```bash
bash run_from_nf.sh final_results HXB2.fasta
```

Example with several samples belonging to the same participant:

```bash
bash run_from_nf.sh final_results HXB2.fasta sample_to_participant.csv
```

Example with a custom output folder but no mapping file:

```bash
bash run_from_nf.sh final_results HXB2.fasta - my_clone_consensus
```

`OUTPUT_DIR` is optional. The default is:

```text
clone_consensus_output/
```

The pipeline creates:

```text
OUTPUT_DIR/
├── work/
└── results/
```

### Sample-to-participant mapping

The optional mapping CSV contains exactly the identifiers needed to combine samples from the same participant:

```csv
sample_id,participant_id
101-1,101
101-2,101
102-1,102
```

Without a mapping file, each `sample_id` is treated as a separate `participant_id`.

### What the direct nf mode uses

The runner finds:

```text
FINAL_RESULTS_DIR/<sample_id>/<sample_id>.annotated.csv
```

For standard `nf-viral-integration` output it uses **host-flanked reads only**:

```text
chromosome != HIV
```

and defines:

```text
clone_id = chromosome_INTEGRATION_SITE
```

It then extracts `READ`, `STRAND`, and `HIV_SEQ` and creates the standard clone-consensus input CSV at:

```text
OUTPUT_DIR/work/nf_clone_input.csv
```

Non-flanked reads are not included automatically in this mode because standard `nf-viral-integration` output does not provide a host integration-site clone assignment for them.

## 3. Generic pre-assigned CSV mode

The consensus pipeline remains independent of how clone membership was assigned. This mode can therefore be used for:

- flanked reads only;
- flanked reads plus successfully assigned non-flanked reads; or
- any other HIV fragments with a trusted clone assignment.

Run:

```bash
bash run_clone_consensus.sh INPUT_CSV REFERENCE_FASTA [OUTPUT_DIR]
```

The default output folder is again `clone_consensus_output/`.

The input CSV must contain:

| Column | Meaning |
| --- | --- |
| `participant_id` | Participant containing the clone |
| `sample_id` | Sample/visit containing the fragment |
| `read` | Read identifier |
| `clone_id` | Pre-assigned clone identifier |
| `STRAND` | `plus` or `minus` |
| `HIV_SEQ` | HIV fragment sequence |

A biological clone is defined by:

```text
participant_id + clone_id
```

`sample_id` does not define the clone, so fragments from several samples of the same participant can contribute to one consensus.

The combination `participant_id + sample_id + read` must be unique.

## 4. Reference FASTA

Provide one full-length HIV reference sequence, for example HXB2.

The reference is used **only as an alignment scaffold**. It helps place randomly sheared fragments at approximately the correct HIV genome coordinates but never contributes bases to the final consensus.

The FASTA must contain exactly one reference sequence.

## 5. Consensus method

### A. Forward normalization

The pipeline normalizes every fragment itself:

```text
STRAND = plus   -> keep HIV_SEQ
STRAND = minus  -> reverse-complement HIV_SEQ
```

The normalized sequences are saved in:

```text
OUTPUT_DIR/work/normalized_sequences.csv
```

### B. Clone grouping

Fragments are grouped by:

```text
participant_id + clone_id
```

Each clone receives a simple internal work identifier such as `clone_000001`. The mapping is stored in `clone_manifest.csv`.

### C. Reference-scaffolded alignment

For each clone, MAFFT aligns the forward-normalized fragments using:

```bash
mafft --addfragments clone_fragments.fasta reference.fasta
```

`--keeplength` is deliberately not used, so supported insertions relative to the reference can be retained.

### D. Coverage-aware majority consensus

The reference row is excluded completely when consensus bases are called.

For each alignment position:

- terminal fragment gaps mean no coverage and do not vote;
- internal gaps can support a deletion;
- `A`, `C`, `G`, `T`, and internal gaps vote;
- ambiguous bases such as `N` do not vote;
- a base or deletion is called only when it has **more than 50%** of informative observations;
- a tie or absence of a majority is called `N`;
- uncovered regions between observed fragments remain `N` rather than being filled from the reference.

Leading and trailing positions with no fragment coverage are removed.

A clone with one fragment is still reported, but it is labelled `single_fragment` rather than `multi_fragment`.

### E. Fragment-to-consensus QC

Each forward-normalized fragment is remapped to its clone consensus with minimap2 `map-hifi`.

This is QC only. The pipeline reports mapping identity and fragment coverage but does not automatically remove discordant fragments.

## 6. Generated folders

```text
OUTPUT_DIR/
├── work/
│   ├── nf_clone_input.csv          # direct nf mode only
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

`work/` contains intermediate and audit/QC files. `results/` contains only the combined outputs intended for downstream use.

## 7. Main outputs

### `clone_consensus.fasta`

One consensus sequence per `participant_id + clone_id`.

Headers use:

```text
>participant_id|clone_id|n=<number_of_fragments>
```

### `clone_consensus_summary.csv`

One row per clone with fragment/sample counts, consensus length, percentage `N`, depth metrics, disagreement positions, and remapping summaries.

### `fragment_qc.csv`

One row per input fragment with fragment length, mapping status, mapping strand, alignment length, percent identity, percent fragment coverage, and MAPQ.

## 8. Important interpretation points

- Clone membership is accepted as given; the consensus pipeline does not infer integration sites.
- Direct `nf-viral-integration` mode derives clone membership only for host-flanked reads.
- The generic CSV mode is the route to include non-flanked reads that have already been assigned to a clone.
- The HIV reference guides alignment only and never supplies consensus bases.
- Internal regions without fragment coverage remain `N`.
- Single-fragment sequences are reported but have less supporting evidence than multi-fragment consensuses.
- Remapping metrics are QC measures, not automatic exclusion thresholds.

## 9. Example generic input

Synthetic example files are in `example_data/`.

Run:

```bash
bash run_clone_consensus.sh example_data/input.csv example_data/reference.fasta
```

The example data are synthetic and are not intended for biological interpretation.

## 10. Files in this repository

```text
run_from_nf.sh           Recommended runner for nf-viral-integration final_results/
run_clone_consensus.sh   Core runner for a pre-assigned clone CSV
prepare_nf_input.py      Adapter from nf-viral-integration final_results/ to the standard input CSV
clone_consensus.py       Validation, normalization, consensus, and QC logic
requirements.txt         Python dependency
example_data/            Synthetic generic-input test data
LICENSE                  MIT license
```

## 11. Reproducibility

For each analysis, record:

- the pipeline version;
- the `nf-viral-integration final_results/` folder or generic input CSV used;
- the optional sample mapping CSV;
- the reference FASTA; and
- the versions of MAFFT and minimap2.

Current release: **v1.0.0**.
