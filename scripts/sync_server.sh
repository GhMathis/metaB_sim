#!/usr/bin/env bash
# Copy the project to the analysis server, and fetch the nf-metabarcoding
# results back. The pipeline runs on the server (nextflow does not run on the
# Windows workstation). Files travel as a tar stream over ssh, since rsync is
# not part of Git for Windows.
#
# Run from Git Bash on the workstation:
#   bash scripts/sync_server.sh push [RUN]   # scripts, config, notes, results/<RUN>,
#                                            # and the SILVA reference if not there yet
#   bash scripts/sync_server.sh pull [RUN]   # results/nf_metabarcoding/<RUN>/
# RUN defaults to chimera_sweep. 'push' overwrites the server copies of the
# files it sends: edit on the workstation, not on the server.
#
# Settings (environment variables):
#   REMOTE      user@host                         default mgheno@194.254.138.192
#   REMOTE_DIR  project directory on the server   default metaB_sim (relative to home)
#   SSH_OPTS    extra ssh options                 e.g. "-i ~/.ssh/id_ed25519_metabsim"

set -euo pipefail

# Work from the project root so the relative paths below resolve.
cd "$(dirname "$0")/.."

ACTION="${1:-}"
RUN="${2:-chimera_sweep}"
REMOTE="${REMOTE:-mgheno@194.254.138.192}"
REMOTE_DIR="${REMOTE_DIR:-metaB_sim}"
SILVA="data/reference/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta"
# SSH_OPTS is unquoted on purpose so that it splits into separate options.
# shellcheck disable=SC2086
remote() { ssh $SSH_OPTS "$REMOTE" "$@"; }

case "$ACTION" in
    push)
        [[ -d "results/$RUN" ]] || { echo "Grinder run not found: results/$RUN" >&2; exit 1; }
        echo "[$(date +%H:%M:%S)] code, config and results/$RUN -> $REMOTE:$REMOTE_DIR"
        tar -czf - README.md LICENSE .gitignore .gitattributes config scripts notes "results/$RUN" \
            | remote "mkdir -p '$REMOTE_DIR' && tar -xzf - -C '$REMOTE_DIR'"

        if remote "test -s '$REMOTE_DIR/$SILVA'"; then
            echo "SILVA reference already on the server"
        else
            [[ -f "$SILVA" ]] || { echo "reference not found: $SILVA" >&2; exit 1; }
            echo "[$(date +%H:%M:%S)] $SILVA -> server (1.5 GB, compressed on the fly)"
            tar -cf - "$SILVA" | gzip -1 \
                | remote "tar -xzf - -C '$REMOTE_DIR'"
        fi
        echo "[$(date +%H:%M:%S)] done"
        ;;
    pull)
        echo "[$(date +%H:%M:%S)] $REMOTE:$REMOTE_DIR/results/nf_metabarcoding/$RUN -> here"
        # --hard-dereference: published results are hard links into the work dir.
        remote "tar -czf - --hard-dereference -C '$REMOTE_DIR' 'results/nf_metabarcoding/$RUN'" \
            | tar -xzf -
        echo "[$(date +%H:%M:%S)] results in results/nf_metabarcoding/$RUN/"
        ;;
    *)
        echo "usage: bash scripts/sync_server.sh push|pull [RUN]" >&2
        exit 1
        ;;
esac
