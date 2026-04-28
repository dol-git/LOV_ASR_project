#!/usr/bin/env python3
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
from itertools import combinations
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

records = list(SeqIO.parse(args.input, "fasta"))

def pid(seq1, seq2):
    aln = pairwise2.align.globalxx(str(seq1), str(seq2), one_alignment_only=True)[0]
    a, b = aln.seqA, aln.seqB
    matches = sum(x == y for x, y in zip(a, b) if x != "-" and y != "-")
    aligned = sum(1 for x, y in zip(a, b) if x != "-" and y != "-")
    return 100 * matches / aligned if aligned else 0

rows = []
scores = {r.id: [] for r in records}

for r1, r2 in combinations(records, 2):
    score = pid(r1.seq, r2.seq)
    scores[r1.id].append(score)
    scores[r2.id].append(score)
    rows.append([r1.id, r2.id, round(score, 2)])

with open(args.output, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["seq1", "seq2", "identity_percent"])
    w.writerows(rows)

summary = args.output.replace(".tsv", ".summary.tsv")
with open(summary, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["seq", "mean_identity", "min_identity", "max_identity", "n_compared"])
    for sid, vals in scores.items():
        if vals:
            w.writerow([sid, round(sum(vals)/len(vals), 2), round(min(vals), 2), round(max(vals), 2), len(vals)])

print(f"[DONE] pairwise: {args.output}")
print(f"[DONE] summary : {summary}")
