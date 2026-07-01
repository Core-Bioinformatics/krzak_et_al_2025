#!/usr/bin/env bash
set -euo pipefail

cd /iss-corescratch/se524/NatImm

THREADS=12
STAR_INDEX="/iss-corescratch/se524/NatImm/reference/star_GRCm39_115_overhang149"
OUTDIR="alignment/star_main"
MANIFEST="${OUTDIR}/star_main_manifest.tsv"

mkdir -p "$OUTDIR"

if [[ ! -f "$STAR_INDEX/SAindex" || ! -f "$STAR_INDEX/Genome" || ! -f "$STAR_INDEX/genomeParameters.txt" ]]; then
  echo "ERROR: STAR index is missing or incomplete: $STAR_INDEX"
  exit 1
fi

grep -E "genomeFastaFiles|sjdbGTFfile|sjdbOverhang|genomeSAindexNbases" "$STAR_INDEX/genomeParameters.txt"

echo
df -h /iss-corescratch

echo -e "sample\tR1\tR2\toutdir\tbam\tlog_final" > "$MANIFEST"

while read -r sample; do
  [[ -z "$sample" ]] && continue

  R1="subsampled_4M/main/${sample}_R1.subsampled_4M.fq.gz"
  R2="subsampled_4M/main/${sample}_R2.subsampled_4M.fq.gz"

  SAMPLE_OUT="${OUTDIR}/${sample}"
  PREFIX="${SAMPLE_OUT}/${sample}."

  BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
  LOG="${PREFIX}Log.final.out"

  mkdir -p "$SAMPLE_OUT"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: missing FASTQ pair for $sample"
    exit 1
  fi

  if [[ -f "$BAM" && -f "$LOG" ]]; then
    echo "[$(date)] SKIP ${sample}: BAM and Log.final.out already exist"
  else
    echo "[$(date)] STAR mapping ${sample}"

    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$R1" "$R2" \
      --readFilesCommand zcat \
      --outFileNamePrefix "$PREFIX" \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM 8000000000 \
      > "${SAMPLE_OUT}/${sample}.STAR.stdout.log" \
      2> "${SAMPLE_OUT}/${sample}.STAR.stderr.log"
  fi

  if [[ ! -f "$BAM" || ! -f "$LOG" ]]; then
    echo "ERROR: STAR output missing for $sample"
    exit 1
  fi

  samtools quickcheck "$BAM"

  echo -e "${sample}\t${R1}\t${R2}\t${SAMPLE_OUT}\t${BAM}\t${LOG}" >> "$MANIFEST"

done < cleanup/samples_keep_after_cleanup.txt

echo
df -h /iss-corescratch

