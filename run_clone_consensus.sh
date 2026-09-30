#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'TXT'
Usage:
  bash run_clone_consensus.sh INPUT_CSV REFERENCE_FASTA OUTPUT_DIR

Arguments:
  INPUT_CSV        CSV containing participant_id, sample_id, read, clone_id, STRAND, HIV_SEQ
  REFERENCE_FASTA  One full-length HIV reference sequence used only as an alignment scaffold
  OUTPUT_DIR       Analysis output folder; work/ and results/ are created inside it
TXT
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if (( $# != 3 )); then
    usage >&2
    exit 1
fi

INPUT_CSV="$1"
REFERENCE_FASTA="$2"
OUTPUT_DIR="$3"
WORK_DIR="$OUTPUT_DIR/work"
RESULTS_DIR="$OUTPUT_DIR/results"

for command_name in python3 mafft minimap2; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "Missing required program: $command_name" >&2
        exit 1
    }
done

python3 -c 'import Bio' >/dev/null 2>&1 || {
    echo "Missing Python package. Install with: python3 -m pip install -r $SCRIPT_DIR/requirements.txt" >&2
    exit 1
}

[[ -s "$INPUT_CSV" ]] || { echo "Missing input CSV: $INPUT_CSV" >&2; exit 1; }
[[ -s "$REFERENCE_FASTA" ]] || { echo "Missing reference FASTA: $REFERENCE_FASTA" >&2; exit 1; }

mkdir -p "$WORK_DIR" "$RESULTS_DIR"
rm -rf "$WORK_DIR/clones"
rm -f "$WORK_DIR/normalized_sequences.csv" "$WORK_DIR/clone_manifest.csv" "$WORK_DIR/reference.fasta"
rm -f "$RESULTS_DIR/clone_consensus.fasta" \
      "$RESULTS_DIR/clone_consensus_summary.csv" \
      "$RESULTS_DIR/fragment_qc.csv"

echo "SMRTcap-clone-consensus version $(python3 "$SCRIPT_DIR/clone_consensus.py" version)"
echo "Input CSV:       $INPUT_CSV"
echo "Reference FASTA: $REFERENCE_FASTA"
echo "Output folder:   $OUTPUT_DIR"
echo "Work folder:     $WORK_DIR"
echo "Results folder:  $RESULTS_DIR"

echo ""
echo "1. Validate input and forward-normalize HIV fragments"
python3 "$SCRIPT_DIR/clone_consensus.py" prepare \
    --input "$INPUT_CSV" \
    --reference "$REFERENCE_FASTA" \
    --work "$WORK_DIR"

echo ""
echo "2. Align fragments, call clone consensuses, and remap fragments for QC"
for clone_dir in "$WORK_DIR"/clones/clone_*; do
    echo "  $(basename "$clone_dir")"

    # The reference is only an alignment scaffold. --keeplength is deliberately
    # not used, so insertions present in clone fragments are retained.
    mafft --quiet --addfragments "$clone_dir/fragments.fasta" "$WORK_DIR/reference.fasta" \
        > "$clone_dir/alignment.fasta"

    # Consensus bases are called only from clone fragments; the REFERENCE row
    # never contributes bases to the consensus.
    python3 "$SCRIPT_DIR/clone_consensus.py" consensus \
        --alignment "$clone_dir/alignment.fasta"

    # Remapping is QC only. No fragment is automatically removed here.
    minimap2 -x map-hifi -c --cs=long --secondary=no \
        "$clone_dir/consensus.fasta" "$clone_dir/fragments.fasta" \
        > "$clone_dir/fragments_vs_consensus.paf"
done

echo ""
echo "3. Combine final consensus sequences and QC summaries"
python3 "$SCRIPT_DIR/clone_consensus.py" finalize \
    --work "$WORK_DIR" \
    --results "$RESULTS_DIR"

echo ""
echo "Finished."
