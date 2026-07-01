#!/usr/bin/env bash
set -euo pipefail

cd /iss-corescratch/se524/NatImm

THREADS=12
STAR_INDEX="/iss-corescratch/se524/NatImm/reference/star_GRCm39_115_overhang149"
INDIR="subsampled_4M/large_sample_iterations"
OUTDIR="alignment/star_large_sample_iterations"
MANIFEST="${OUTDIR}/star_large_iterations_manifest.tsv"

mkdir -p "$OUTDIR"

if [[ ! -f "$STAR_INDEX/SAindex" || ! -f "$STAR_INDEX/Genome" || ! -f "$STAR_INDEX/genomeParameters.txt" ]]; then
  echo "ERROR: STAR index is missing or incomplete: $STAR_INDEX"
  exit 1
fi

grep -E "genomeFastaFiles|sjdbGTFfile|sjdbOverhang|genomeSAindexNbases" "$STAR_INDEX/genomeParameters.txt"

echo
df -h /iss-corescratch

echo -e "sample\titeration\tR1\tR2\toutdir\tbam\tlog_final" > "$MANIFEST"

for sample in S7 S19 S43; do
  for iter in 1 2 3 4 5; do

    label="${sample}_iter${iter}"

    R1="${INDIR}/${label}_R1.subsampled_4M.fq.gz"
    R2="${INDIR}/${label}_R2.subsampled_4M.fq.gz"

    SAMPLE_OUT="${OUTDIR}/${label}"
    PREFIX="${SAMPLE_OUT}/${label}."

    BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
    LOG="${PREFIX}Log.final.out"

    mkdir -p "$SAMPLE_OUT"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
      echo "ERROR: missing FASTQ pair for ${label}"
      echo "R1: $R1"
      echo "R2: $R2"
      exit 1
    fi

    if [[ -f "$BAM" && -f "$LOG" ]]; then
      echo "[$(date)] SKIP ${label}: BAM and Log.final.out already exist"
    else
      echo "[$(date)] STAR mapping ${label}"

      STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$PREFIX" \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 8000000000 \
        > "${SAMPLE_OUT}/${label}.STAR.stdout.log" \
        2> "${SAMPLE_OUT}/${label}.STAR.stderr.log"
    fi

    if [[ ! -f "$BAM" || ! -f "$LOG" ]]; then
      echo "ERROR: STAR output missing for ${label}"
      exit 1
    fi

    samtools quickcheck "$BAM"

    echo -e "${sample}\t${iter}\t${R1}\t${R2}\t${SAMPLE_OUT}\t${BAM}\t${LOG}" >> "$MANIFEST"

  done
done

echo
df -h /iss-corescratch

