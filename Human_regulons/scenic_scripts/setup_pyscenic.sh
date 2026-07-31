#!/usr/bin/env bash

set -euo pipefail

PROJECT="/iss-scratch/CoreBioinformatics/rk720/human_sucnr1"
CONTAINER_DIR="${PROJECT}/scenic/containers"
IMAGE="${CONTAINER_DIR}/pyscenic_0.12.1.sif"

mkdir -p "${CONTAINER_DIR}"

if command -v apptainer >/dev/null 2>&1; then
    ENGINE="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    ENGINE="singularity"
else
    echo "Neither Apptainer nor Singularity was found."
    exit 1
fi

echo "Using container engine: ${ENGINE}"

if [[ ! -f "${IMAGE}" ]]; then
    "${ENGINE}" build \
        "${IMAGE}" \
        docker://aertslab/pyscenic:0.12.1
else
    echo "Container already exists: ${IMAGE}"
fi

"${ENGINE}" exec \
    --bind "${PROJECT}:/data" \
    "${IMAGE}" \
    pyscenic --version

echo "pySCENIC container is ready."