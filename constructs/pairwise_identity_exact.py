#!/usr/bin/env python3
from Bio import SeqIO, pairwise2
from itertools import combinations
import argparse
import csv
import os

parser = argparse.ArgumentParser(
    description="Calculate pairwise amino-acid identity and mean identity for FASTA sequences."
)
parser.add_argument("-i", "--input", required=True, help="Input FASTA")
parser.add_argument("-o", "--output", required=True, help="Output pairwise TSV")
parser.add_argument("--gap-open", type=float, default=-10)
parser.add_argument("--gap-extend", type=float, default=-0.5)
args = parser.parse_args()

records = list(SeqIO.parse(args.input, "fasta"))

if len(records) < 2:
    raise SystemExit("[ERROR] Need at least two sequences.")

def identity(seq1, seq2):
    aln = pairwise2.align.globalms(
        str(seq1), str(seq2),
        1, 0,                 # match, mismatch
        args.gap_open, args.gap_extend,
        one_alignment_only=True
    )[0]

    a, b = aln.seqA, aln.seqB
    matches = 0
    compared = 0
    gaps = 0

    for x, y in zip(a, b):
        if x == "-" or y == "-":
            gaps += 1
            continue
        compared += 1
        if x == y:
            matches += 1

    pid = 100 * matches / compared if compared else 0
    coverage1 = 100 * compared / len(seq1) if len(seq1) else 0
    coverage2 = 100 * compared / len(seq2) if len(seq2) else 0

    return pid, compared, matches, gaps, coverage1, coverage2

scores = {r.id: [] for r in records}
rows = []

for r1, r2 in combinations(records, 2):
    pid, compared, matches, gaps, cov1, cov2 = identity(r1.seq, r2.seq)
    rows.append([
        r1.id, r2.id,
        round(pid, 2),
        compared,
        matches,
        gaps,
        round(cov1, 2),
        round(cov2, 2),
        len(r1.seq),
        len(r2.seq)
    ])
    scores[r1.id].append(pid)
    scores[r2.id].append(pid)

with open(args.output, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow([
        "seq1", "seq2",
        "identity_percent",
        "compared_positions",
        "matches",
        "gap_positions",
        "coverage_seq1_percent",
        "coverage_seq2_percent",
        "len_seq1",
        "len_seq2"
    ])
    w.writerows(rows)

summary = os.path.splitext(args.output)[0] + ".summary.tsv"

with open(summary, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["seq", "mean_identity", "min_identity", "max_identity", "n_compared"])
    for sid, vals in scores.items():
        w.writerow([
            sid,
            round(sum(vals) / len(vals), 2),
            round(min(vals), 2),
            round(max(vals), 2),
            len(vals)
        ])

print(f"[DONE] pairwise table: {args.output}")
print(f"[DONE] summary table : {summary}")
