#!/usr/bin/env bash
# Rerun the whole workflow. Run from the repository root:  bash scripts/run_all.sh
# Set PYTHON to choose the interpreter, and SKIP_DOWNLOAD=1 to reuse data already in data/raw/.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python}"

if [ "${SKIP_DOWNLOAD:-0}" != "1" ]; then
  "$PY" 01_download_data.py
fi
"$PY" 02_point_cloud_qa.py
"$PY" 03_dtm_dsm_chm.py
"$PY" 04_accuracy_vs_usgs_dem.py
"$PY" 05_terrain_derivatives.py
"$PY" 06_hydrology.py
"$PY" 07_compare_nhd_interpret.py
"$PY" 08_export_products.py
echo "All steps finished. See outputs/figures, outputs/tables, outputs/rasters and outputs/vectors."
