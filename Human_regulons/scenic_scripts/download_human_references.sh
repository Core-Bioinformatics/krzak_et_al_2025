#!/usr/bin/env bash

set -euo pipefail

PROJECT="/iss-scratch/CoreBioinformatics/rk720/human_sucnr1"
RESOURCE_DIR="${PROJECT}/scenic/resources"

mkdir -p "${RESOURCE_DIR}"
cd "${RESOURCE_DIR}"

DB_BASE="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based"

DB_10KB="hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

DB_PROMOTER="hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

MOTIF_ANNOTATIONS="motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

TF_LIST="allTFs_hg38.txt"

echo "Downloading hg38 ranking databases..."

wget -c \
  "${DB_BASE}/${DB_10KB}"

wget -c \
  "${DB_BASE}/${DB_PROMOTER}"

echo "Downloading motif annotations..."

wget -c \
  "https://resources.aertslab.org/cistarget/motif2tf/${MOTIF_ANNOTATIONS}"

echo "Downloading human transcription-factor list..."

wget -c \
  "https://resources.aertslab.org/cistarget/tf_lists/${TF_LIST}"

echo "Downloading checksums..."

wget -c \
  "https://resources.aertslab.org/cistarget/databases/sha256sum.txt"

echo "Checking ranking-database checksums..."

for database in "${DB_10KB}" "${DB_PROMOTER}"; do
    checksum_line="$(
      awk \
        -v feather_database="${database}" \
        '$2 == feather_database' \
        sha256sum.txt
    )"

    if [[ -n "${checksum_line}" ]]; then
        echo "${checksum_line}" | sha256sum -c -
    else
        echo "Warning: checksum not found for ${database}"
    fi
done

echo
echo "Downloaded resources:"
ls -lh \
  "${DB_10KB}" \
  "${DB_PROMOTER}" \
  "${MOTIF_ANNOTATIONS}" \
  "${TF_LIST}"