#!/usr/bin/env python3
"""
make_consensus.py

Build a threshold-based consensus sequence from a multiple-sequence alignment
in FASTA format.  For each alignment column:

  - columns where gap fraction >= gap_threshold are skipped (removed_gap_sites)
  - the most-common non-gap residue is written
  - if its frequency among non-gap residues is below --threshold the site is
    counted as ambiguous_sites; with --lowercase-ambiguous it is written in
    lower case

Usage:
    python make_consensus.py -i aligned.fasta -o consensus.fasta -n MyConsensus
"""

import argparse
import sys
import textwrap
from collections import Counter
from pathlib import Path

from Bio import AlignIO


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Make threshold-based consensus sequence from aligned FASTA."
    )
    p.add_argument("-i", "--input", required=True, type=Path,
                   help="Input aligned FASTA")
    p.add_argument("-o", "--output", required=True, type=Path,
                   help="Output consensus FASTA")
    p.add_argument("-n", "--name", required=True,
                   help="Consensus sequence name")
    p.add_argument("-t", "--threshold", type=float, default=0.7,
                   help="Residue frequency threshold among non-gap residues")
    p.add_argument("--gap-threshold", type=float, default=0.5,
                   help="Remove column if gap fraction is >= this value")
    p.add_argument("--lowercase-ambiguous", action="store_true",
                   help="Write low-confidence majority residues as lowercase "
                        "instead of uppercase")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        sys.exit(f"[ERROR] Input file not found: {args.input}")

    aln  = AlignIO.read(str(args.input), "fasta")
    nseq = len(aln)

    if nseq == 0:
        sys.exit("[ERROR] Alignment contains no sequences.")

    consensus         = []
    ambiguous_sites   = 0
    removed_gap_sites = 0

    for i in range(aln.get_alignment_length()):
        col      = aln[:, i]
        gap_frac = col.count("-") / nseq

        if gap_frac >= args.gap_threshold:
            removed_gap_sites += 1
            continue

        residues = [c.upper() for c in col if c != "-"]
        if not residues:
            continue

        counts    = Counter(residues)
        aa, count = counts.most_common(1)[0]
        freq      = count / len(residues)

        if freq >= args.threshold:
            consensus.append(aa)
        else:
            ambiguous_sites += 1
            consensus.append(aa.lower() if args.lowercase_ambiguous else aa)

    seq = "".join(consensus)

    with open(args.output, "w", encoding="utf-8") as f:
        f.write(
            f">{args.name} length={len(seq)} threshold={args.threshold} "
            f"gap_threshold={args.gap_threshold} "
            f"ambiguous_sites={ambiguous_sites} "
            f"removed_gap_sites={removed_gap_sites}\n"
        )
        for line in textwrap.wrap(seq, 60):
            f.write(line + "\n")

    print(f"[DONE] {args.output}")
    print(f"[INFO] consensus length: {len(seq)}")
    print(f"[INFO] ambiguous sites below threshold: {ambiguous_sites}")
    print(f"[INFO] removed high-gap columns: {removed_gap_sites}")


if __name__ == "__main__":
    main()
