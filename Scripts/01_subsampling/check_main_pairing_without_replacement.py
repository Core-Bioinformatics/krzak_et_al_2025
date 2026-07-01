from pathlib import Path
import gzip
import pandas as pd

samples = [x.strip() for x in Path("cleanup/samples_keep_after_cleanup.txt").read_text().splitlines() if x.strip()]
outdir = Path("subsampled_4M/main")
audit = Path("subsampled_4M/audit")
audit.mkdir(parents=True, exist_ok=True)

def norm_id(header):
    h = header.strip()
    if h.startswith("@"):
        h = h[1:]
    h = h.split()[0]
    if h.endswith("/1") or h.endswith("/2"):
        h = h[:-2]
    return h

rows = []

for s in samples:
    print(f"Checking {s}", flush=True)

    r1 = outdir / f"{s}_R1.subsampled_4M.fq.gz"
    r2 = outdir / f"{s}_R2.subsampled_4M.fq.gz"

    n = 0
    mismatch = 0
    duplicate_ids = 0
    seen = set()

    with gzip.open(r1, "rt", errors="replace") as f1, gzip.open(r2, "rt", errors="replace") as f2:
        while True:
            h1 = f1.readline()
            h2 = f2.readline()

            if not h1 and not h2:
                break

            if not h1 or not h2:
                mismatch += 1
                break

            f1.readline()
            f1.readline()
            f1.readline()

            f2.readline()
            f2.readline()
            f2.readline()

            id1 = norm_id(h1)
            id2 = norm_id(h2)

            if id1 != id2:
                mismatch += 1

            if id1 in seen:
                duplicate_ids += 1
            else:
                seen.add(id1)

            n += 1

    rows.append({
        "sample": s,
        "read_pairs_checked": n,
        "expected_read_pairs": 4_000_000,
        "pair_count_ok": n == 4_000_000,
        "R1_R2_id_mismatches": mismatch,
        "duplicate_read_ids_in_subsample": duplicate_ids,
        "paired_ids_ok": mismatch == 0,
        "no_duplicate_ids": duplicate_ids == 0,
    })

    pd.DataFrame(rows).to_csv(audit / "main_pairing_and_duplicates.tsv", sep="\t", index=False)

df = pd.DataFrame(rows)
df.to_csv(audit / "main_pairing_and_duplicates.tsv", sep="\t", index=False)

bad = df[
    (~df["pair_count_ok"]) |
    (~df["paired_ids_ok"]) |
    (~df["no_duplicate_ids"])
]

print()
print(df.to_string(index=False))

if len(bad):
    print()
    print("BAD SAMPLES:")
    print(bad.to_string(index=False))
    raise SystemExit(1)

print()
