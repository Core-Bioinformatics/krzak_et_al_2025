#!/usr/bin/env bash
set -euo pipefail

cd /iss-corescratch/se524/NatImm

N=4000000
OUTDIR="subsampled_4M/large_sample_iterations"
MANIFEST="subsampled_4M/large_sample_iterations_manifest.tsv"

mkdir -p "$OUTDIR"

echo "[$(date)] Starting large sample iteration setup"

existing=$(find "$OUTDIR" -name "*.fq.gz" | wc -l)
if [[ "$existing" -ne 0 ]]; then
    echo "ERROR: $OUTDIR already contains $existing FASTQ files."
    exit 1
fi

echo -e "sample\titeration\tseed\tinput_R1\tinput_R2\toutput_R1\toutput_R2\ttarget_read_pairs\tnote" > "$MANIFEST"

for sample in S7 S19 S43; do
    main_r1="subsampled_4M/main/${sample}_R1.subsampled_4M.fq.gz"
    main_r2="subsampled_4M/main/${sample}_R2.subsampled_4M.fq.gz"

    out1="${OUTDIR}/${sample}_iter1_R1.subsampled_4M.fq.gz"
    out2="${OUTDIR}/${sample}_iter1_R2.subsampled_4M.fq.gz"

    if [[ ! -f "$main_r1" || ! -f "$main_r2" ]]; then
        echo "missing validated main subsample for $sample"
        exit 1
    fi

    ln -s "../../main/${sample}_R1.subsampled_4M.fq.gz" "$out1"
    ln -s "../../main/${sample}_R2.subsampled_4M.fq.gz" "$out2"

    echo -e "${sample}\t1\tmain_existing\t${main_r1}\t${main_r2}\t${out1}\t${out2}\t${N}\tvalidated_main_subsample" >> "$MANIFEST"
done

for sample in S7 S19 S43; do
  for iter in 2 3 4 5; do

    if [[ "$sample" == "S7" ]]; then
      base_seed=7000
    elif [[ "$sample" == "S19" ]]; then
      base_seed=19000
    elif [[ "$sample" == "S43" ]]; then
      base_seed=43000
    fi

    seed=$((base_seed + iter))

    r1="fastq/${sample}_R1.fq.gz"
    r2="fastq/${sample}_R2.fq.gz"

    out1="${OUTDIR}/${sample}_iter${iter}_R1.subsampled_4M.fq.gz"
    out2="${OUTDIR}/${sample}_iter${iter}_R2.subsampled_4M.fq.gz"

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "missing raw FASTQ pair for $sample"
        exit 1
    fi

    echo "[$(date)] Subsampling ${sample} iteration ${iter} with seed ${seed}"

    seqtk sample -s"$seed" "$r1" "$N" | pigz -p 4 > "$out1"
    seqtk sample -s"$seed" "$r2" "$N" | pigz -p 4 > "$out2"

    echo -e "${sample}\t${iter}\t${seed}\t${r1}\t${r2}\t${out1}\t${out2}\t${N}\tnew_independent_subsample" >> "$MANIFEST"

  done
done

echo "[$(date)] Done large sample iterations 1-5"
