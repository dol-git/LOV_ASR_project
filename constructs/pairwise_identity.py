#!/usr/bin/env python3
"""
pairwise_identity.py

Compute global pairwise sequence identity for all pairs in a FASTA file using
global alignment without gap penalty (Bio.pairwise2.align.globalxx).

Outputs two files:
  <output>             — pairwise TSV:  seq1, seq2, identity_percent
  <output stem>.summary.tsv — per-sequence mean/min/max identity

NOTE: Bio.pairwise2 is deprecated in Biopython >= 1.80.
      Future migration: replace with Bio.Align.PairwiseAligner.

Usage:
    python pairwise_identity.py -i sequences.fasta -o pairwise.tsv
"""

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path

from Bio import SeqIO, pairwise2


def pairwise_id(seq1, seq2) -> float:
    """Return percent identity from an ungapped global alignment."""
    aln = pairwise2.align.globalxx(str(seq1), str(seq2), one_alignment_only=True)[0]
    a, b = aln.seqA, aln.seqB
    matches = sum(x == y for x, y in zip(a, b) if x != "-" and y != "-")
    aligned = sum(1 for x, y in zip(a, b) if x != "-" and y != "-")
    return 100 * matches / aligned if aligned else 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute global pairwise sequence identity for all pairs in a FASTA."
    )
    p.add_argument("-i", "--input",  required=True, type=Path,
                   help="Input FASTA file")
    p.add_argument("-o", "--output", required=True, type=Path,
                   help="Output pairwise TSV")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        sys.exit(f"[ERROR] Input file not found: {args.input}")

    records = list(SeqIO.parse(str(args.input), "fasta"))

    summary_path = Path(str(args.output).replace(".tsv", ".summary.tsv"))

    rows   = []
    scores = {r.id: [] for r in records}

    for r1, r2 in combinations(records, 2):
        score = pairwise_id(r1.seq, r2.seq)
        scores[r1.id].append(score)
        scores[r2.id].append(score)
        rows.append([r1.id, r2.id, round(score, 2)])

    with open(args.output, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["seq1", "seq2", "identity_percent"])
        w.writerows(rows)

    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["seq", "mean_identity", "min_identity", "max_identity", "n_compared"])
        for sid, vals in scores.items():
            if vals:
                w.writerow([
                    sid,
                    round(sum(vals) / len(vals), 2),
                    round(min(vals), 2),
                    round(max(vals), 2),
                    len(vals),
                ])

    print(f"[DONE] pairwise: {args.output}")
    print(f"[DONE] summary : {summary_path}")


if __name__ == "__main__":
    main()
