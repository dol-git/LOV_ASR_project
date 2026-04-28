#!/usr/bin/env python3
from Bio import AlignIO
import argparse
import textwrap

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignment", required=True)
parser.add_argument("--ref-id", required=True)
parser.add_argument("--target-id", required=True)
parser.add_argument("--ref-start", type=int, required=True)
parser.add_argument("--ref-end", type=int, required=True)
parser.add_argument("-o", "--output", required=True)
parser.add_argument("--name", required=True)
args = parser.parse_args()

aln = AlignIO.read(args.alignment, "fasta")

records = {r.id: str(r.seq) for r in aln}

if args.ref_id not in records:
    raise SystemExit(f"[ERROR] ref-id not found: {args.ref_id}\nAvailable: {list(records)}")
if args.target_id not in records:
    raise SystemExit(f"[ERROR] target-id not found: {args.target_id}\nAvailable: {list(records)}")

ref = records[args.ref_id]
target = records[args.target_id]

ref_pos = 0
selected_cols = []

for i, aa in enumerate(ref):
    if aa != "-":
        ref_pos += 1
    if args.ref_start <= ref_pos <= args.ref_end and aa != "-":
        selected_cols.append(i)

target_region = "".join(target[i] for i in selected_cols)
target_region = target_region.replace("-", "")

with open(args.output, "w") as f:
    f.write(f">{args.name} mapped_from_{args.ref_id}_{args.ref_start}_{args.ref_end} length={len(target_region)}\n")
    for line in textwrap.wrap(target_region, 60):
        f.write(line + "\n")

print(f"[DONE] {args.output}")
print(f"[INFO] extracted length: {len(target_region)}")
print(f"[INFO] selected alignment columns: {len(selected_cols)}")
