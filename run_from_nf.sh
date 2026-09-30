#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'TXT'
Usage:
  bash run_from_nf.sh FINAL_RESULTS_DIR REFERENCE_FASTA [SAMPLE_MAP_CSV] [OUTPUT_DIR]

Arguments:
  FINAL_RESULTS_DIR  nf-viral-integration final_results/ folder
  REFERENCE_FASTA    One full-length HIV reference used only as an alignment scaffold
  SAMPLE_MAP_CSV     Optional CSV with sample_id,participant_id columns. Use - for none.
  OUTPUT_DIR         Optional analysis output folder (default: clone_consensus_output)

Standard nf-viral-integration input mode uses host-flanked reads only
(chromosome != HIV) and defines clone_id as chromosome_INTEGRATION_SITE.
TXT
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if (( $# < 2 || $# > 4 )); then
    usage >&2
    exit 1
fi

FINAL_RESULTS_DIR="$1"
REFERENCE_FASTA="$2"
SAMPLE_MAP="${3:-}"
[[ "$SAMPLE_MAP" == "-" ]] && SAMPLE_MAP=""
OUTPUT_DIR="${4:-clone_consensus_output}"
WORK_DIR="$OUTPUT_DIR/work"
INPUT_CSV="$WORK_DIR/nf_clone_input.csv"

[[ -d "$FINAL_RESULTS_DIR" ]] || { echo "Missing final-results folder: $FINAL_RESULTS_DIR" >&2; exit 1; }
[[ -s "$REFERENCE_FASTA" ]] || { echo "Missing reference FASTA: $REFERENCE_FASTA" >&2; exit 1; }
if [[ -n "$SAMPLE_MAP" && ! -s "$SAMPLE_MAP" ]]; then
    echo "Missing sample mapping CSV: $SAMPLE_MAP" >&2
    exit 1
fi

mkdir -p "$WORK_DIR"

echo "Preparing clone input from nf-viral-integration output"
if [[ -n "$SAMPLE_MAP" ]]; then
    python3 "$SCRIPT_DIR/prepare_nf_input.py" \
        "$FINAL_RESULTS_DIR" "$INPUT_CSV" \
        --sample-map "$SAMPLE_MAP"
else
    python3 "$SCRIPT_DIR/prepare_nf_input.py" \
        "$FINAL_RESULTS_DIR" "$INPUT_CSV"
fi

echo ""
bash "$SCRIPT_DIR/run_clone_consensus.sh" "$INPUT_CSV" "$REFERENCE_FASTA" "$OUTPUT_DIR"
