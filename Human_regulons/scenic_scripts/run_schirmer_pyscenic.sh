#!/usr/bin/env bash

set -euo pipefail

PROJECT="/iss-scratch/CoreBioinformatics/rk720/human_sucnr1"

IMAGE="${PROJECT}/scenic/containers/pyscenic_0.12.1.sif"

INPUT="${PROJECT}/scenic/input/Schirmer/Schirmer_microglia_scenic_input.loom"

RESOURCE_DIR="${PROJECT}/scenic/resources"
RESULT_DIR="${PROJECT}/scenic/results/Schirmer"
LOG_DIR="${PROJECT}/scenic/logs/Schirmer"

TF_LIST="${RESOURCE_DIR}/allTFs_hg38.txt"

DB_10KB="${RESOURCE_DIR}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

DB_PROMOTER="${RESOURCE_DIR}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

MOTIF_ANNOTATIONS="${RESOURCE_DIR}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

ADJACENCIES="${RESULT_DIR}/Schirmer_adjacencies.tsv"
REGULONS="${RESULT_DIR}/Schirmer_regulons.csv"
AUC_MATRIX="${RESULT_DIR}/Schirmer_regulon_auc.csv"

N_WORKERS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${RESULT_DIR}" "${LOG_DIR}"

for required_file in \
  "${IMAGE}" \
  "${INPUT}" \
  "${TF_LIST}" \
  "${DB_10KB}" \
  "${DB_PROMOTER}" \
  "${MOTIF_ANNOTATIONS}"
do
  if [[ ! -s "${required_file}" ]]; then
    echo "Missing or empty file: ${required_file}"
    exit 1
  fi
done

RUNNER=(
  apptainer
  exec
  --env
  "HDF5_USE_FILE_LOCKING=FALSE"
  --bind
  "${PROJECT}:/data"
  "${IMAGE}"
)

echo "Dataset: Schirmer microglia"
echo "Workers: ${N_WORKERS}"

echo
echo "Step 1/3: GRNBoost2"

"${RUNNER[@]}" \
  pyscenic grn \
  --num_workers "${N_WORKERS}" \
  --method grnboost2 \
  --seed 42 \
  --output "/data/scenic/results/Schirmer/Schirmer_adjacencies.tsv" \
  "/data/scenic/input/Schirmer/Schirmer_microglia_scenic_input.loom" \
  "/data/scenic/resources/allTFs_hg38.txt" \
  2>&1 | tee "${LOG_DIR}/01_grn.log"

test -s "${ADJACENCIES}"

echo
echo "Step 2/3: motif enrichment and pruning"

"${RUNNER[@]}" \
  pyscenic ctx \
  "/data/scenic/results/Schirmer/Schirmer_adjacencies.tsv" \
  "/data/scenic/resources/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  "/data/scenic/resources/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  --annotations_fname \
  "/data/scenic/resources/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
  --expression_mtx_fname \
  "/data/scenic/input/Schirmer/Schirmer_microglia_scenic_input.loom" \
  --mode custom_multiprocessing \
  --num_workers "${N_WORKERS}" \
  --output "/data/scenic/results/Schirmer/Schirmer_regulons.csv" \
  2>&1 | tee "${LOG_DIR}/02_ctx.log"

test -s "${REGULONS}"

echo
echo "Step 3/3: AUCell"

"${RUNNER[@]}" \
  pyscenic aucell \
  "/data/scenic/input/Schirmer/Schirmer_microglia_scenic_input.loom" \
  "/data/scenic/results/Schirmer/Schirmer_regulons.csv" \
  --num_workers "${N_WORKERS}" \
  --seed 42 \
  --output "/data/scenic/results/Schirmer/Schirmer_regulon_auc.csv" \
  2>&1 | tee "${LOG_DIR}/03_aucell.log"

test -s "${AUC_MATRIX}"

echo
echo "Schirmer pySCENIC completed."

ls -lh \
  "${ADJACENCIES}" \
  "${REGULONS}" \
  "${AUC_MATRIX}"