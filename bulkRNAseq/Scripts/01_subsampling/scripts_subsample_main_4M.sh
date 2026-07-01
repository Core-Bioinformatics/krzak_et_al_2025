#!/usr/bin/env bash
set -euo pipefail

cd /iss-corescratch/se524/NatImm

N=4000000
OUTDIR="subsampled_4M/main"
MANIFEST="subsampled_4M/main_subsampling_manifest.tsv"

mkdir -p "$OUTDIR"

echo -e "sample\tseed\tinput_R1\tinput_R2\toutput_R1\toutput_R2\ttarget_read_pairs" > "$MANIFEST"

i=0
while read -r sample; do
    [[ -z "$sample" ]] && continue

    i=$((i + 1))
    seed=$((1000 + i))

    r1="fastq/${sample}_R1.fq.gz"
    r2="fastq/${sample}_R2.fq.gz"

    out1="${OUTDIR}/${sample}_R1.subsampled_4M.fq.gz"
    out2="${OUTDIR}/${sample}_R2.subsampled_4M.fq.gz"

    echo "[$(date)] Subsampling ${sample} with seed ${seed}"

    seqtk sample -s"$seed" "$r1" "$N" | pigz -p 4 > "$out1"
    seqtk sample -s"$seed" "$r2" "$N" | pigz -p 4 > "$out2"

    echo -e "${sample}\t${seed}\t${r1}\t${r2}\t${out1}\t${out2}\t${N}" >> "$MANIFEST"

done < cleanup/samples_keep_after_cleanup.txt

echo "[$(date)] Done main subsampling"
