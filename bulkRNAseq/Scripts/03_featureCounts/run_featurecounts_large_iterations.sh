#!/usr/bin/env bash
set -euo pipefail

cd /iss-corescratch/se524/NatImm

GTF="/iss-corescratch/se524/GSE261455/mapping/references/Mus_musculus.GRCm39.115.gtf"
OUTDIR="quantification/featurecounts_large_iterations"
BAM_LIST="${OUTDIR}/large_iterations_bam_files.txt"
OUT_COUNTS="${OUTDIR}/featureCounts_large_iterations_gene_counts.txt"

mkdir -p "$OUTDIR"

if [[ ! -f "$GTF" ]]; then
  echo "ERROR: GTF not found: $GTF"
  exit 1
fi

rm -f "$BAM_LIST"

for sample in S7 S19 S43; do
  for iter in 1 2 3 4 5; do
    label="${sample}_iter${iter}"
    bam="alignment/star_large_sample_iterations/${label}/${label}.Aligned.sortedByCoord.out.bam"

    if [[ ! -f "$bam" ]]; then
      echo "ERROR: missing BAM for $label: $bam"
      exit 1
    fi

    samtools quickcheck "$bam"
    echo "$bam" >> "$BAM_LIST"
  done
done

echo "BAM files:"
wc -l "$BAM_LIST"
cat "$BAM_LIST"

featureCounts \
  -T 12 \
  -p \
  -B \
  -C \
  -t exon \
  -g gene_id \
  -a "$GTF" \
  -o "$OUT_COUNTS" \
  $(cat "$BAM_LIST")

