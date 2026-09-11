#!/usr/bin/env bash
set -euo pipefail

PREPROC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${PREPROC_DIR}/.." && pwd)"
CONFIG="${PREPROC_DIR}/configs/pbmc.yaml"
LOG_DIR="${PREPROC_DIR}/logs"

mkdir -p "${LOG_DIR}"
cd "${PROJECT_ROOT}"

echo "[ATLAS] Running PBMC RNA-ATAC preprocessing..."

python "${PREPROC_DIR}/scripts/preprocess_pbmc_rna_atac.py" \
  --config "${CONFIG}" \
  2>&1 | tee "${LOG_DIR}/pbmc_rna_atac.log"

echo "[ATLAS] PBMC preprocessing finished."
