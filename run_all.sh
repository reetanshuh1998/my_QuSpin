#!/bin/bash
# =============================================================
#  MBL Simulation Suite — Overnight Runner
#  Runs all 6 scripts sequentially and logs output + timing.
#  Usage:  nohup ./run_all.sh &
# =============================================================

PYTHON="./venv/bin/python"
LOGFILE="run_all.log"

echo "=============================================" | tee $LOGFILE
echo "  MBL Simulation Suite — Started: $(date)"     | tee -a $LOGFILE
echo "=============================================" | tee -a $LOGFILE

SCRIPTS=(
    "01_reproduce_fig3.py"
    "02_temperature_effect.py"
    "03_heat_reservoir_effect.py"
    "04_strong_interaction_effect.py"
    "05_entanglement_entropy.py"
    "06_level_spacing.py"
)

TOTAL_START=$SECONDS

for script in "${SCRIPTS[@]}"; do
    echo "" | tee -a $LOGFILE
    echo ">>> Starting $script at $(date)" | tee -a $LOGFILE
    echo "---------------------------------------------" | tee -a $LOGFILE
    
    START=$SECONDS
    $PYTHON "$script" 2>&1 | tee -a $LOGFILE
    EXIT_CODE=$?
    ELAPSED=$(( SECONDS - START ))
    MINS=$(( ELAPSED / 60 ))
    SECS=$(( ELAPSED % 60 ))
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "<<< DONE: $script — ${MINS}m ${SECS}s ✓" | tee -a $LOGFILE
    else
        echo "<<< FAILED: $script (exit code $EXIT_CODE) — ${MINS}m ${SECS}s ✗" | tee -a $LOGFILE
    fi
done

TOTAL_ELAPSED=$(( SECONDS - TOTAL_START ))
TOTAL_MINS=$(( TOTAL_ELAPSED / 60 ))
TOTAL_SECS=$(( TOTAL_ELAPSED % 60 ))

echo "" | tee -a $LOGFILE
echo "=============================================" | tee -a $LOGFILE
echo "  ALL DONE — Total: ${TOTAL_MINS}m ${TOTAL_SECS}s"  | tee -a $LOGFILE
echo "  Finished: $(date)"                              | tee -a $LOGFILE
echo "=============================================" | tee -a $LOGFILE
