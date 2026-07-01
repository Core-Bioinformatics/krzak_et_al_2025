#!/usr/bin/env bash
set -euo pipefail

cd /iss-corescratch/se524/NatImm

GTF="/iss-corescratch/se524/GSE261455/mapping/references/Mus_musculus.GRCm39.115.gtf"
OUTDIR="quantification/featurecounts_main"
BAM_LIST="${OUTDIR}/main_bam_files.txt"
OUT_COUNTS="${OUTDIR}/featureCounts_main_gene_counts.txt"

mkdir -p "$OUTDIR"

if [[ ! -f "$GTF" ]]; then
  echo "ERROR: GTF not found: $GTF"
  exit 1
fi

rm -f "$BAM_LIST"

while read -r sample; do
  [[ -z "$sample" ]] && continue

  bam="alignment/star_main/${sample}/${sample}.Aligned.sortedByCoord.out.bam"

  if [[ ! -f "$bam" ]]; then
    echo "ERROR: missing BAM for $sample: $bam"
    exit 1
  fi

  samtools quickcheck "$bam"
  echo "$bam" >> "$BAM_LIST"

done < cleanup/samples_keep_after_cleanup.txt

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

echo "Done"
