#!/usr/bin/env bash
# Stage the reads of a Grinder run as nf-metabarcoding inputs.
# For each results/<run>/<sample>-reads.fastq, write
# results/nf_input/<run>/<sample>.fastq with RNA bases converted to DNA
# (U -> T: the SILVA reference, hence the Grinder reads, are RNA, whereas
# vsearch, swarm and cutadapt expect DNA), and list every sample in a
# samplesheet (nf-metabarcoding --input, single-end: fastq_2 left empty).
# Read headers and quality lines are copied unchanged.
#
# Run from the project root (or any directory), once per run. The chimera sweep
# has one run per chimera level:
#   bash scripts/03_prepare_nf_inputs.sh chimera_sweep/chimera_05
#   for run in results/chimera_sweep/chimera_??; do
#       bash scripts/03_prepare_nf_inputs.sh "chimera_sweep/${run##*/}"
#   done

set -euo pipefail

# Work from the project root so the relative paths below resolve (the
# samplesheet stores paths relative to it, nextflow is launched from it).
cd "$(dirname "$0")/.."

RUN="${1:?usage: bash scripts/03_prepare_nf_inputs.sh RUN (e.g. chimera_sweep/chimera_05)}"
IN_DIR="results/$RUN"
OUT_DIR="results/nf_input/$RUN"
SAMPLESHEET="$OUT_DIR/samplesheet.csv"

[[ -d "$IN_DIR" ]] || { echo "Grinder run not found: $IN_DIR" >&2; exit 1; }
shopt -s nullglob
fastqs=("$IN_DIR"/*-reads.fastq)
[[ ${#fastqs[@]} -gt 0 ]] || { echo "no *-reads.fastq in $IN_DIR" >&2; exit 1; }

mkdir -p "$OUT_DIR"
echo "sample,fastq_1,fastq_2" > "$SAMPLESHEET"

printf '%-12s %8s %9s %s\n' sample reads_in reads_out status
failures=0
for fastq in "${fastqs[@]}"; do
    sample=$(basename "$fastq" -reads.fastq)
    # nf-metabarcoding accepts [A-Za-z0-9._-] sample IDs, not starting with
    # '-' or '.', and reserves the 'notmerged' suffix.
    if [[ ! $sample =~ ^[A-Za-z0-9_][A-Za-z0-9._-]*$ || $sample == *notmerged ]]; then
        echo "invalid sample ID for nf-metabarcoding: $sample" >&2
        exit 1
    fi
    out="$OUT_DIR/$sample.fastq"

    # Grinder writes 4-line records: line 2 of each record is the sequence.
    awk 'NR % 4 == 2 { gsub(/U/, "T"); gsub(/u/, "t") } { print }' "$fastq" > "$out"

    # Checks: same number of records, well-formed records, no U left.
    reads_in=$(awk 'END { print NR / 4 }' "$fastq")
    reads_out=$(awk 'END { print NR / 4 }' "$out")
    malformed=$(awk '(NR % 4 == 1 && !/^@/) || (NR % 4 == 3 && !/^\+/) { n++ }
                     END { print n + (NR % 4 != 0) }' "$out")
    rna_left=$(awk 'NR % 4 == 2 && /[Uu]/ { n++ } END { print n + 0 }' "$out")
    status=OK
    if [[ $reads_in != "$reads_out" || $malformed -ne 0 || $rna_left -ne 0 ]]; then
        status=FAIL
        failures=$((failures + 1))
    fi
    printf '%-12s %8s %9s %s\n' "$sample" "$reads_in" "$reads_out" "$status"

    echo "$sample,$out," >> "$SAMPLESHEET"
done

echo
if [[ $failures -gt 0 ]]; then
    echo "$failures sample(s) failed the checks" >&2
    exit 1
fi
echo "${#fastqs[@]} samples staged in $OUT_DIR/, samplesheet: $SAMPLESHEET"
