#!/usr/bin/env bash
# Build the taxonomic reference of nf-metabarcoding (--reference_dataset) from
# SILVA: the region amplified by the primers that grinder uses, in stampa format.
#
#   1. Primers: the file given to grinder's -forward_reverse in
#      config/grinder_params.txt (first record forward, second reverse).
#   2. For each SILVA sequence (RNA converted to DNA, CRLF line endings removed),
#      find the forward primer, then downstream the reverse complement of the
#      reverse primer, with the same exact IUPAC matching as grinder. Keep the
#      region between them, primers excluded, as nf-metabarcoding trims the
#      primers off the reads. Sequences without both primer sites are dropped.
#   3. Stampa headers, lineage levels separated by '|' and spaces replaced by '_':
#        SILVA   >HG374866.9.1343 Bacteria;Pseudomonadati;...;Rhizoctonia solani
#        stampa  >HG374866.9.1343 Bacteria|Pseudomonadati|...|Rhizoctonia_solani
#
# Why the amplified region and not the full-length sequences: stampa scores hits
# with vsearch --iddef 1, which counts terminal gaps, so a ~250-nt amplicon
# aligned to a ~1,450-nt full-length sequence never reaches --stampa_id 0.5.
#
# The output is named after the primer file and skipped when it already exists:
# delete it to rebuild (for example after editing the primer sequences).
#
# Run from the project root (or any directory):
#   bash scripts/04_build_stampa_reference.sh
#   bash scripts/04_build_stampa_reference.sh path/to/other_silva.fasta[.gz]

set -euo pipefail

# Work from the project root so the relative paths below resolve.
cd "$(dirname "$0")/.."

SILVA="${1:-data/reference/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta}"
GRINDER_PROFILE="config/grinder_params.txt"
OUT_DIR="results/reference"

[[ -f "$SILVA" ]] || { echo "reference not found: $SILVA" >&2; exit 1; }
PRIMERS=$(awk '$1 == "-forward_reverse" { print $2 }' "$GRINDER_PROFILE")
[[ -f "$PRIMERS" ]] \
    || { echo "primer file not found: '$PRIMERS' (-forward_reverse in $GRINDER_PROFILE)" >&2; exit 1; }
primer_seqs=$(awk '{ sub(/\r$/, "") } /^>/ { n++; next } { seq[n] = seq[n] toupper($0) }
                   END { print seq[1]; print seq[2] }' "$PRIMERS")
FORWARD=$(sed -n 1p <<< "$primer_seqs")
REVERSE=$(sed -n 2p <<< "$primer_seqs")
[[ $FORWARD =~ ^[ACGTURYSWKMBDHVN]+$ && $REVERSE =~ ^[ACGTURYSWKMBDHVN]+$ ]] \
    || { echo "$PRIMERS must hold two IUPAC primers (forward, then reverse)" >&2; exit 1; }

# primers_16S_V4.fasta -> ..._16S_V4_stampa.fasta (scripts/05 derives the same name)
marker=$(basename "$PRIMERS"); marker=${marker%.fasta}; marker=${marker#primers_}
name=$(basename "${SILVA%.gz}")
OUT="$OUT_DIR/${name%.fasta}_${marker}_stampa.fasta"
if [[ -s "$OUT" ]]; then
    echo "already built: $OUT (delete it to rebuild)"
    exit 0
fi

mkdir -p "$OUT_DIR"
reader=(cat)
[[ $SILVA == *.gz ]] && reader=(gzip -dc)

echo "[$(date +%H:%M:%S)] extracting the $FORWARD ... $REVERSE region from $SILVA"
"${reader[@]}" "$SILVA" | awk -v fwd="$FORWARD" -v rev="$REVERSE" '
    # IUPAC primer -> regular expression (U is matched as T: sequences are DNA here)
    function to_regex(primer,    i, c, re) {
        re = ""
        for (i = 1; i <= length(primer); i++) {
            c = substr(primer, i, 1)
            re = re (c in iupac ? iupac[c] : c)
        }
        return re
    }
    function revcomp(primer,    i, rc) {
        rc = ""
        for (i = length(primer); i >= 1; i--) rc = rc comp[substr(primer, i, 1)]
        return rc
    }
    function emit(    region) {
        if (seq == "") return
        n_in++
        if (!match(seq, fwd_re)) return
        region = substr(seq, RSTART + RLENGTH)
        if (!match(region, rev_re)) return
        region = substr(region, 1, RSTART - 1)
        # Shortest amplicon, as grinder: start after the last forward primer site
        # found before the reverse primer site.
        while (match(region, fwd_re)) region = substr(region, RSTART + RLENGTH)
        if (region == "") return
        print header
        print region
        n_out++
        len = length(region)
        lengths[len]++
        if (min == "" || len < min) min = len
        if (len > max) max = len
    }
    BEGIN {
        split("A:A C:C G:G T:T U:T R:[AG] Y:[CT] S:[CG] W:[AT] K:[GT] M:[AC] " \
              "B:[CGT] D:[AGT] H:[ACT] V:[ACG] N:[ACGT]", pairs, " ")
        for (i in pairs) { split(pairs[i], kv, ":"); iupac[kv[1]] = kv[2] }
        split("A:T C:G G:C T:A U:A R:Y Y:R S:S W:W K:M M:K B:V V:B D:H H:D N:N", pairs, " ")
        for (i in pairs) { split(pairs[i], kv, ":"); comp[kv[1]] = kv[2] }
        fwd_re = to_regex(fwd)
        rev_re = to_regex(revcomp(rev))
    }
    { sub(/\r$/, "") }
    /^>/ {
        emit()
        seq = ""
        id = substr($1, 2)
        lineage = substr($0, length($1) + 2)
        gsub(/ /, "_", lineage)
        gsub(/;/, "|", lineage)
        header = ">" id " " lineage
        next
    }
    {
        line = toupper($0)
        gsub(/U/, "T", line)
        seq = seq line
    }
    END {
        emit()
        # Length distribution of the extracted regions, to stderr
        for (len = min; len <= max; len++) {
            cum += lengths[len]
            if (p01 == "" && cum >= 0.01 * n_out) p01 = len
            if (p50 == "" && cum >= 0.50 * n_out) p50 = len
            if (p99 == "" && cum >= 0.99 * n_out) p99 = len
        }
        printf "sequences: %d, with both primer sites: %d (%.1f %%)\n",
            n_in, n_out, (n_in > 0 ? 100 * n_out / n_in : 0) > "/dev/stderr"
        printf "region length without primers: min %d, 1 %% %d, median %d, 99 %% %d, max %d\n",
            min, p01, p50, p99, max > "/dev/stderr"
    }
' > "$OUT.tmp"

# Checks: something extracted, one sequence line per header, no U left.
awk '
    /^>/ { headers++; if (prev_header) bad++; prev_header = 1; next }
    { if (!prev_header) bad++; prev_header = 0; if (/U/) rna++ }
    END {
        printf "check: %d sequences, layout errors: %d, U left: %d\n", headers, bad + 0, rna + 0
        exit !(headers > 0 && bad == 0 && rna == 0)
    }' "$OUT.tmp" \
    || { echo "check failed; partial output kept in $OUT.tmp" >&2; exit 1; }

mv "$OUT.tmp" "$OUT"
echo "[$(date +%H:%M:%S)] stampa reference: $OUT"
