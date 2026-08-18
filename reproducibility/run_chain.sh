#!/bin/bash
# Waits for the currently running genEvidDoE_colocAndCaviar_new.py (PID passed as $1)
# to finish, checks its log for success, and if successful launches
# colocAndCaviarDoE_biosampleNameCombinations.py next. Runs detached via nohup so it
# survives the launching terminal/session disconnecting.

set -u

JOB1_PID="$1"
REPRO_DIR="/home/jmroldan/jr-doe-temp1810_machineData/reproducibility"
JOB1_OUT="$2"
STATUS_LOG="$REPRO_DIR/chain_status.log"
PYTHONPATH_DIR="/home/jmroldan/jr-doe-temp1810_machineData/src/data/analysis"

echo "$(date '+%Y-%m-%d %H:%M:%S') chain watcher started, waiting on PID $JOB1_PID" >> "$STATUS_LOG"

while kill -0 "$JOB1_PID" 2>/dev/null; do
    sleep 30
done

echo "$(date '+%Y-%m-%d %H:%M:%S') job1 (PID $JOB1_PID) is no longer running" >> "$STATUS_LOG"

if grep -q "Analysis finished" "$JOB1_OUT" 2>/dev/null && ! grep -qE "Traceback \(most recent call last\)|ERROR ApplicationMaster" "$JOB1_OUT" 2>/dev/null; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') job1 finished successfully -> launching job2" >> "$STATUS_LOG"
    cd "$REPRO_DIR" || exit 1
    JOB2_OUT="$REPRO_DIR/run_colocAndCaviarDoE_$(date +%Y%m%d_%H%M%S).out"
    PYTHONPATH="$PYTHONPATH_DIR" nohup python3 colocAndCaviarDoE_biosampleNameCombinations.py > "$JOB2_OUT" 2>&1 &
    JOB2_PID=$!
    echo "$(date '+%Y-%m-%d %H:%M:%S') job2 launched, PID $JOB2_PID, output -> $JOB2_OUT" >> "$STATUS_LOG"
else
    echo "$(date '+%Y-%m-%d %H:%M:%S') job1 did NOT finish successfully (no 'Analysis finished' marker, or an error was found) -> job2 NOT launched. Check $JOB1_OUT" >> "$STATUS_LOG"
fi
