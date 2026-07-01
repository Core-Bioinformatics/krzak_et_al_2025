#!/usr/bin/env bash
set -euo pipefail

cd /iss-corescratch/se524/NatImm

OUT="subsampled_4M/audit/large_iterations_read_counts.tsv"

echo -e "file\tlines\treads\tlines_mod_4\tread_count_ok\tfastq_structure_ok" > "$OUT"

for fq in subsampled_4M/large_sample_iterations/*.fq.gz; do
    echo "Checking $fq"

    lines=$(gzip -cd "$fq" | wc -l)
    reads=$((lines / 4))
    mod=$((lines % 4))

    if [[ "$reads" -eq 4000000 ]]; then
        read_ok="TRUE"
    else
        read_ok="FALSE"
    fi

    if [[ "$mod" -eq 0 ]]; then
        struct_ok="TRUE"
    else
        struct_ok="FALSE"
    fi

    echo -e "${fq}\t${lines}\t${reads}\t${mod}\t${read_ok}\t${struct_ok}" >> "$OUT"
done

echo "DONE large iterations read count check"
echo
echo "Summary:"
cut -f5,6 "$OUT" | sort | uniq -c
