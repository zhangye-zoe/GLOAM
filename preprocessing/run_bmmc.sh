#!/usr/bin/env bash
set -euo pipefail

PREPROC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${PREPROC_DIR}/.." && pwd)"
CONFIG="${PREPROC_DIR}/configs/bmmc.yaml"
LOG_DIR="${PREPROC_DIR}/logs"

mkdir -p "${LOG_DIR}"
cd "${PROJECT_ROOT}"

echo "[ATLAS] Running BMMC RNA-ATAC preprocessing..."
python "${PREPROC_DIR}/scripts/preprocess_bmmc_rna_atac.py" --config "${CONFIG}" \
  2>&1 | tee "${LOG_DIR}/bmmc_rna_atac.log"

echo "[ATLAS] Running BMMC RNA-Protein preprocessing..."
python "${PREPROC_DIR}/scripts/preprocess_bmmc_rna_protein.py" --config "${CONFIG}" \
  2>&1 | tee "${LOG_DIR}/bmmc_rna_protein.log"

echo "[ATLAS] All BMMC preprocessing finished."
