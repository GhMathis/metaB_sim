#!/usr/bin/env bash
# Chimera sweep: 10 chimera levels (0 to 9 % chimeric reads) x 5 replicates = 50 samples.
#
# 1. Community: grinder draws it once (config/grinder_community.txt) among the SILVA
#    references that carry the primer pair. Its references and relative abundances
#    are saved in results/chimera_sweep/community/.
# 2. Samples: every sample receives that community (-abundance_file) and, as reference,
#    only those sequences. This is much faster than re-reading SILVA and gives the same
#    reads (checked: grinder skips references missing from the abundance file anyway).
#    Each sample also gets its
#    -chimera_perc and its own -random_seed (100 x level + replicate, so the 50 seeds
#    differ). Every other option comes from config/grinder_params.txt.
# 3. Checks on every sample, and a summary table: results/chimera_sweep/summary.tsv.
#
# Output: one folder per chimera level, each to be processed as one nf-metabarcoding run:
#   results/chimera_sweep/chimera_05/chimera_05_r1-reads.fastq ... chimera_05_r5-reads.fastq
#
# Run from the project root (or any directory). Grinder is installed in WSL:
#   wsl -d Ubuntu -- bash scripts/02_run_grinder_chimera_sweep.sh
# If grinder is not on PATH, set GRINDER, for example:
#   GRINDER="perl /path/to/grinder.pl" bash scripts/02_run_grinder_chimera_sweep.sh

set -euo pipefail

# Work from the project root so the relative paths below resolve.
cd "$(dirname "$0")/.."

# Grinder builds chimeras by walking Perl hashes, whose order Perl shuffles in every
# process. Without a fixed hash seed, the same -random_seed reproduces the
# non-chimeric reads but not the chimeras (checked on chimera_09_r3). Fixing it makes
# every sample reproducible from its seed.
export PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0

GRINDER="${GRINDER:-grinder}"
SILVA="data/reference/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta"
COMMUNITY_PROFILE="config/grinder_community.txt"
PROFILE="config/grinder_params.txt"
OUT_DIR="results/chimera_sweep"
COMMUNITY_DIR="$OUT_DIR/community"
ABUNDANCES="$COMMUNITY_DIR/community-abundances.txt"
REFERENCES="$COMMUNITY_DIR/community-references.fasta"
SUMMARY="$OUT_DIR/summary.tsv"
CHIMERA_LEVELS=$(seq 0 9)   # percentage of chimeric reads, one level per run
REPLICATES=$(seq 1 5)       # samples per level, each with its own seed

command -v "${GRINDER%% *}" > /dev/null || { echo "grinder not found: set GRINDER" >&2; exit 1; }
for file in "$SILVA" "$COMMUNITY_PROFILE" "$PROFILE"; do
    [[ -f "$file" ]] || { echo "not found: $file" >&2; exit 1; }
done
if [[ -e "$OUT_DIR" ]]; then
    echo "$OUT_DIR already exists: move or delete it before simulating a new sweep" >&2
    exit 1
fi

expected_reads=$(awk '$1 == "-total_reads" { print $2 }' "$PROFILE")
expected_richness=$(awk '$1 == "-diversity" { print $2 }' "$COMMUNITY_PROFILE")
primers=$(awk '$1 == "-forward_reverse" { print $2 }' "$PROFILE")
[[ -n "$expected_reads" && -n "$expected_richness" && -f "$primers" ]] \
    || { echo "need -total_reads and -forward_reverse in $PROFILE, -diversity in $COMMUNITY_PROFILE" >&2; exit 1; }

# GRINDER is unquoted on purpose below, so that "perl grinder.pl" splits into words.
# shellcheck disable=SC2086
run_grinder() {
    local log="$1"
    shift
    if ! $GRINDER "$@" > "$log" 2>&1; then
        echo "grinder failed; see $log" >&2
        exit 1
    fi
}

# 1. Community, drawn once among the references carrying both primers
mkdir -p "$COMMUNITY_DIR"
echo "[$(date +%H:%M:%S)] community: $expected_richness references drawn from $SILVA"
run_grinder "$COMMUNITY_DIR/community.log" \
    -reference_file "$SILVA" -profile_file "$COMMUNITY_PROFILE" \
    -forward_reverse "$primers" -total_reads 1 \
    -base_name community -output_dir "$COMMUNITY_DIR"
# ranks file: rank, reference ID, relative abundance (%) -> abundance file: ID, abundance
awk '!/^#/ { print $2, $3 }' "$COMMUNITY_DIR/community-ranks.txt" > "$ABUNDANCES"
# the community's reference sequences, without SILVA's CRLF line endings
awk 'NR == FNR { keep[$1]; next }
     { sub(/\r$/, "") }
     /^>/ { keep_record = (substr($1, 2) in keep) }
     keep_record' "$ABUNDANCES" "$SILVA" > "$REFERENCES"
n_abundances=$(wc -l < "$ABUNDANCES")
n_references=$(grep -c '^>' "$REFERENCES")
if [[ $n_abundances -ne $expected_richness || $n_references -ne $expected_richness ]]; then
    echo "community: $n_abundances abundances, $n_references references (expected $expected_richness)" >&2
    exit 1
fi

# 2. Samples
samples=()
for level in $CHIMERA_LEVELS; do
    level_dir=$(printf '%s/chimera_%02d' "$OUT_DIR" "$level")
    mkdir -p "$level_dir"
    for replicate in $REPLICATES; do
        base=$(printf 'chimera_%02d_r%d' "$level" "$replicate")
        seed=$((100 * level + replicate))
        samples+=("$level_dir/$base")
        echo "[$(date +%H:%M:%S)] $base: -chimera_perc $level -random_seed $seed"
        run_grinder "$level_dir/$base.log" \
            -reference_file "$REFERENCES" -profile_file "$PROFILE" \
            -abundance_file "$ABUNDANCES" \
            -chimera_perc "$level" -random_seed "$seed" \
            -base_name "$base" -output_dir "$level_dir"
    done
done

# 3. Checks on every sample: configured read depth, the community's richness and
# references, and at most one distinct sequence per reference among non-chimeric
# reads (no sequencing errors are simulated). Chimeric reads list several
# references in their header: reference=A,B.
failures=0
echo -e "sample\tchimera_level\treplicate\tseed\treads\tchimeric_reads\tchimeric_pct\tdistinct_amplicons\tdistinct_chimeras" > "$SUMMARY"
echo
echo "expected: reads=$expected_reads richness=$expected_richness"
printf '%-14s %6s %9s %9s %10s %9s %s\n' sample reads chimeric "chim_%" amplicons chimeras status
for sample in "${samples[@]}"; do
    base=$(basename "$sample")
    level=$((10#$(sed 's/^chimera_\([0-9]*\)_r.*/\1/' <<< "$base")))
    replicate=${base##*_r}
    seed=$((100 * level + replicate))
    richness=$(awk '!/^#/ { n++ } END { print n + 0 }' "$sample-ranks.txt")
    same_references=yes
    cmp -s <(awk '!/^#/ { print $2 }' "$sample-ranks.txt" | sort) <(cut -d' ' -f1 "$ABUNDANCES" | sort) \
        || same_references=no
    # -fastq_output 1: 4-line records; header line 1, sequence line 2
    read -r reads chimeric amplicons chimeras < <(paste - - - - < "$sample-reads.fastq" | awk -F'\t' '
        { n++; split($1, f, " "); chimeric = (f[2] ~ /,/) }
        chimeric { c++; if (!($2 in chim)) { chim[$2]; nc++ } }
        !chimeric { if (!($2 in amp)) { amp[$2]; na++ } }
        END { print n + 0, c + 0, na + 0, nc + 0 }')
    pct=$(awk -v c="$chimeric" -v n="$reads" 'BEGIN { printf "%.1f", n ? 100 * c / n : 0 }')

    status=OK
    if [[ $reads -ne $expected_reads || $richness -ne $expected_richness \
          || $same_references != yes || $amplicons -gt $expected_richness ]]; then
        status=FAIL
        failures=$((failures + 1))
    fi
    printf '%-14s %6s %9s %9s %10s %9s %s\n' "$base" "$reads" "$chimeric" "$pct" "$amplicons" "$chimeras" "$status"
    echo -e "$base\t$level\t$replicate\t$seed\t$reads\t$chimeric\t$pct\t$amplicons\t$chimeras" >> "$SUMMARY"
done

echo
if [[ $failures -gt 0 ]]; then
    echo "$failures sample(s) failed the checks" >&2
    exit 1
fi
echo "All ${#samples[@]} samples passed the checks. Outputs: $OUT_DIR/, summary: $SUMMARY"
