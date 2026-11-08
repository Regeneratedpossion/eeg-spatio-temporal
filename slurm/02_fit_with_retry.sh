#!/bin/bash
# Wrapper: retries code/02_fit/02_fit_bgam.R up to 5 times on failure (mitigates
# sporadic GCC internal-compiler-error segfaults, ~12% observed rate, not
# tied to any specific seed or node).
set -uo pipefail

MAX_ATTEMPTS=5
SEED_PADDED=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
SEED_DIR="${PROJECT_ROOT}/results/single_node_cond${COND}_ds_cr_seeds001_100/seed_${SEED_PADDED}"

for attempt in $(seq 1 $MAX_ATTEMPTS); do
    echo "attempt ${attempt}/${MAX_ATTEMPTS} for seed ${SLURM_ARRAY_TASK_ID} at $(date --iso-8601=seconds)"
    Rscript --vanilla code/02_fit/02_fit_bgam.R
    status=$?
    if [ $status -eq 0 ]; then
        echo "succeeded on attempt ${attempt}"
        exit 0
    fi
    echo "attempt ${attempt} failed with status ${status}; clearing partial output and retrying"
    rm -rf "$SEED_DIR"
    sleep 10
done

echo "all ${MAX_ATTEMPTS} attempts failed for seed ${SLURM_ARRAY_TASK_ID}"
exit 1
