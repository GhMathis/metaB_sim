#!/usr/bin/env bash
# Run nf-metabarcoding (https://github.com/frederic-mahe/nf-metabarcoding) on
# a staged Grinder run: Part A (per-sample dereplication and clustering),
# Part B (occurrence table: swarm, chimera detection, mumu) and Part C
# (stampa taxonomic assignment against SILVA).
#
# Runs on a Unix machine with nextflow, not on the Windows workstation. A run is
# a folder of Grinder samples under results/; the chimera sweep has one run per
# chimera level (chimera_sweep/chimera_00 ... chimera_sweep/chimera_09, 5 samples
# each). Stage the inputs first:
#   bash scripts/03_prepare_nf_inputs.sh chimera_sweep/chimera_05
#   bash scripts/04_build_stampa_reference.sh                # once
# then, from the project root and inside tmux or screen (the nextflow driver
# must stay alive for the whole run):
#   bash scripts/05_run_nf_metabarcoding.sh chimera_sweep/chimera_05
#   NF_PROFILE=local,singularity bash scripts/05_run_nf_metabarcoding.sh chimera_sweep/chimera_05
# Output files are prefixed with the last part of the run name (chimera_05).
# Running it again resumes from the cached tasks (-resume).
#
# Parameters: config/nf_metabarcoding_params.yml (analysis),
# config/nf_metabarcoding_server.config (resources on the server), and the
# primers of the grinder simulation (-forward_reverse in config/grinder_params.txt).

set -euo pipefail

# Work from the project root: nextflow resolves the relative paths of the
# params file and of the samplesheet against the launch directory.
cd "$(dirname "$0")/.."

RUN="${1:?usage: bash scripts/05_run_nf_metabarcoding.sh RUN (e.g. chimera_sweep/chimera_05)}"
PIPELINE="frederic-mahe/nf-metabarcoding"
PIPELINE_VERSION="v0.6.0"   # pinned release: change it deliberately, then rerun everything
# 'local' feeds the whole representatives fasta to a single vsearch for
# taxonomy (no scatter); the second profile provides vsearch, swarm, cutadapt
# and mumu: conda (bioconda), singularity, apptainer, docker, or none (tools on PATH).
NF_PROFILE="${NF_PROFILE:-local,conda}"
PARAMS="config/nf_metabarcoding_params.yml"
SITE_CONFIG="config/nf_metabarcoding_server.config"
GRINDER_PROFILE="config/grinder_params.txt"
SILVA="data/reference/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta"
SAMPLESHEET="results/nf_input/$RUN/samplesheet.csv"
OUT_DIR="results/nf_metabarcoding/$RUN"
WORK_DIR="results/nf_work/$RUN"
# Build the conda environment once and reuse it across runs and work dirs.
export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-$HOME/.nextflow/conda_cache}"
# Cap the heap of the nextflow driver (it runs on a shared machine).
export NXF_OPTS="${NXF_OPTS:--Xms512m -Xmx4g}"

command -v nextflow > /dev/null || { echo "nextflow not found on PATH" >&2; exit 1; }
[[ -f "$PARAMS" ]] || { echo "params file not found: $PARAMS" >&2; exit 1; }
[[ -f "$SITE_CONFIG" ]] || { echo "site config not found: $SITE_CONFIG" >&2; exit 1; }
[[ -f "$SAMPLESHEET" ]] \
    || { echo "samplesheet not found: $SAMPLESHEET (run scripts/03_prepare_nf_inputs.sh $RUN)" >&2; exit 1; }

# Primers: the pair grinder used to simulate the amplicons (first record
# forward, second reverse), so that simulation and analysis cannot diverge.
PRIMERS=$(awk '$1 == "-forward_reverse" { print $2 }' "$GRINDER_PROFILE")
[[ -f "$PRIMERS" ]] \
    || { echo "primer file not found: '$PRIMERS' (-forward_reverse in $GRINDER_PROFILE)" >&2; exit 1; }
primer_seqs=$(awk '{ sub(/\r$/, "") } /^>/ { n++; next } { seq[n] = seq[n] toupper($0) }
                   END { print seq[1]; print seq[2] }' "$PRIMERS")
FORWARD=$(sed -n 1p <<< "$primer_seqs")
REVERSE=$(sed -n 2p <<< "$primer_seqs")
[[ -n $FORWARD && -n $REVERSE ]] \
    || { echo "$PRIMERS must hold two primers (forward, then reverse)" >&2; exit 1; }

# Taxonomic reference for that primer pair, named as scripts/04 names it.
marker=$(basename "$PRIMERS"); marker=${marker%.fasta}; marker=${marker#primers_}
name=$(basename "$SILVA")
REFERENCE="results/reference/${name%.fasta}_${marker}_stampa.fasta"
[[ -s "$REFERENCE" ]] \
    || { echo "reference not found: $REFERENCE (run scripts/04_build_stampa_reference.sh)" >&2; exit 1; }

echo "[$(date +%H:%M:%S)] $PIPELINE $PIPELINE_VERSION, profile $NF_PROFILE, run $RUN"
echo "primers $FORWARD / $REVERSE, reference $REFERENCE"
nextflow run "$PIPELINE" -r "$PIPELINE_VERSION" \
    -profile "$NF_PROFILE" \
    -c "$SITE_CONFIG" \
    -params-file "$PARAMS" \
    -work-dir "$WORK_DIR" \
    -resume \
    --input "$SAMPLESHEET" \
    --project_name "${RUN##*/}" \
    --outdir "$OUT_DIR" \
    --forward_primer "$FORWARD" \
    --reverse_primer "$REVERSE" \
    --reference_dataset "$REFERENCE"

# Keep the driver log next to the execution reports of this run.
cp .nextflow.log "$OUT_DIR/pipeline_info/nextflow.log"

echo
echo "Occurrence tables in $OUT_DIR/occurrence_table/:"
ls -l "$OUT_DIR/occurrence_table/"
echo "Per-sample read counts: $OUT_DIR/logs/"
echo "work dir (needed by -resume, delete when done): $WORK_DIR"
