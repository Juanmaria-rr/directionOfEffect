#!/bin/bash
# Runs colocAndCaviarDoE_biosampleNameCombinations_dedupPairs.py twice, sequentially:
#   1) PROPAGATE_DISEASES=true   -> outputs tagged _dedupPairs_propag
#   2) PROPAGATE_DISEASES=false  -> outputs tagged _dedupPairs_noPropag
# Sequential on purpose: the single-node Spark cluster cannot host two sessions
# (OOM, already documented). Each run took ~19h on 2026-08-09.
# Launched with nohup+disown so it survives the terminal and the Claude session.

set -u

REPRO_DIR="/home/jmroldan/jr-doe-temp1810_machineData/reproducibility"
PYTHONPATH_DIR="/home/jmroldan/jr-doe-temp1810_machineData/src/data/analysis"
SCRIPT="colocAndCaviarDoE_biosampleNameCombinations_dedupPairs.py"
STATUS_LOG="$REPRO_DIR/dedup_chain_status.log"

cd "$REPRO_DIR" || exit 1

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" >> "$STATUS_LOG"; }

run_variant() {
    local propagate="$1"      # true | false
    local tag="$2"            # propag | noPropag
    local out="$REPRO_DIR/run_dedupPairs_${tag}_$(date +%Y%m%d_%H%M%S).out"

    log "launching variant '$tag' (PROPAGATE_DISEASES=$propagate) -> $out"
    PYTHONPATH="$PYTHONPATH_DIR" PROPAGATE_DISEASES="$propagate" \
        python3 "$SCRIPT" > "$out" 2>&1
    local rc=$?

    if [ $rc -eq 0 ] && grep -q "Analysis finished" "$out" 2>/dev/null \
       && ! grep -qE "Traceback \(most recent call last\)|ERROR ApplicationMaster" "$out" 2>/dev/null; then
        log "variant '$tag' finished OK"
        return 0
    fi

    log "variant '$tag' FAILED (exit=$rc, no 'Analysis finished' marker or an error was found). Check $out"
    return 1
}

# Optional argument: run a single variant instead of the whole chain.
#   ./run_dedup_chain.sh noPropag   -> only PROPAGATE_DISEASES=false
#   ./run_dedup_chain.sh propag     -> only PROPAGATE_DISEASES=true
#   ./run_dedup_chain.sh            -> both, sequentially (default)
ONLY="${1:-}"

case "$ONLY" in
    propag)
        log "=== dedup chain started (single variant: propag) ==="
        run_variant "true" "propag"
        ;;
    noPropag)
        log "=== dedup chain started (single variant: noPropag) ==="
        run_variant "false" "noPropag"
        ;;
    "")
        log "=== dedup chain started ==="
        if run_variant "true" "propag"; then
            run_variant "false" "noPropag"
        else
            log "propagated variant failed -> non-propagated variant NOT launched"
        fi
        ;;
    *)
        log "unknown variant '$ONLY' (expected 'propag', 'noPropag' or no argument)"
        echo "unknown variant '$ONLY' (expected 'propag', 'noPropag' or no argument)" >&2
        exit 2
        ;;
esac

log "=== dedup chain finished ==="
