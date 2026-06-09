#!/usr/bin/env python3
"""
pairwise_identity_exact.py

Compute global pairwise sequence identity for all pairs in a FASTA file using
a scored global alignment (match=1, mismatch=0, configurable gap penalties).

Outputs two files:
  <output>                  — pairwise TSV with per-pair detail columns
  <output stem>.summary.tsv — per-sequence mean/min/max identity

Column definitions (pairwise TSV):
  identity_percent     — 100 * matches / compared_positions
  compared_positions   — aligned positions excluding gaps on either side
  matches              — identical compared positions
  gap_positions        — positions where at least one side is a gap
  coverage_seqN_percent— 100 * compared_positions / length of seqN

NOTE: Bio.pairwise2 is deprecated in Biopython >= 1.80.
      Future migration: replace with Bio.Align.PairwiseAligner.

Usage:
    python pairwise_identity_exact.py -i sequences.fasta -o pairwise.tsv
"""

import argparse
import csv
import os
from itertools import combinations

from Bio import SeqIO, pairwise2


def sequence_identity(seq1, seq2, gap_open: float, gap_extend: float):
    """
    Return (pid, compared, matches, gaps, coverage1, coverage2) from a
    scored global alignment.
    """
    aln = pairwise2.align.globalms(
        str(seq1), str(seq2),
        1, 0,
        gap_open, gap_extend,
        one_alignment_only=True,
    )[0]

    a, b = aln.seqA, aln.seqB
    matches  = 0
    compared = 0
    gaps     = 0

    for x, y in zip(a, b):
        if x == "-" or y == "-":
            gaps += 1
            continue
        compared += 1
        if x == y:
            matches += 1

    pid       = 100 * matches / compared  if compared      else 0
    coverage1 = 100 * compared / len(seq1) if len(seq1) > 0 else 0
    coverage2 = 100 * compared / len(seq2) if len(seq2) > 0 else 0

    return pid, compared, matches, gaps, coverage1, coverage2


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Calculate pairwise amino-acid identity and mean identity "
                    "for FASTA sequences."
    )
    p.add_argument("-i", "--input",      required=True,
                   help="Input FASTA")
    p.add_argument("-o", "--output",     required=True,
                   help="Output pairwise TSV")
    p.add_argument("--gap-open",   type=float, default=-10,
                   help="Gap-open penalty (default: -10)")
    p.add_argument("--gap-extend", type=float, default=-0.5,
                   help="Gap-extend penalty (default: -0.5)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    records = list(SeqIO.parse(args.input, "fasta"))

    if len(records) < 2:
        raise SystemExit("[ERROR] Need at least two sequences.")

    summary_path = os.path.splitext(args.output)[0] + ".summary.tsv"

    scores = {r.id: [] for r in records}
    rows   = []

    for r1, r2 in combinations(records, 2):
        pid, compared, matches, gaps, cov1, cov2 = sequence_identity(
            r1.seq, r2.seq, args.gap_open, args.gap_extend
        )
        rows.append([
            r1.id, r2.id,
            round(pid, 2),
            compared,
            matches,
            gaps,
            round(cov1, 2),
            round(cov2, 2),
            len(r1.seq),
            len(r2.seq),
        ])
        scores[r1.id].append(pid)
        scores[r2.id].append(pid)

    with open(args.output, "w", encoding="utf-8", newline="") as f:
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
            "len_seq2",
        ])
        w.writerows(rows)

    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["seq", "mean_identity", "min_identity", "max_identity", "n_compared"])
        for sid, vals in scores.items():
            w.writerow([
                sid,
                round(sum(vals) / len(vals), 2),
                round(min(vals), 2),
                round(max(vals), 2),
                len(vals),
            ])

    print(f"[DONE] pairwise table: {args.output}")
    print(f"[DONE] summary table : {summary_path}")


if __name__ == "__main__":
    main()
