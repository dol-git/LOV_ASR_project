#!/usr/bin/env python3
from Bio import AlignIO
from collections import Counter
import argparse
import textwrap

parser = argparse.ArgumentParser(
    description="Make threshold-based consensus sequence from aligned FASTA."
)
parser.add_argument("-i", "--input", required=True, help="Input aligned FASTA")
parser.add_argument("-o", "--output", required=True, help="Output consensus FASTA")
parser.add_argument("-n", "--name", required=True, help="Consensus sequence name")
parser.add_argument("-t", "--threshold", type=float, default=0.7,
                    help="Residue frequency threshold among non-gap residues")
parser.add_argument("--gap-threshold", type=float, default=0.5,
                    help="Remove column if gap fraction is >= this value")
parser.add_argument("--lowercase-ambiguous", action="store_true",
                    help="Write low-confidence majority residues as lowercase instead of uppercase")
args = parser.parse_args()

aln = AlignIO.read(args.input, "fasta")
nseq = len(aln)
consensus = []
ambiguous_sites = 0
removed_gap_sites = 0

for i in range(aln.get_alignment_length()):
    col = aln[:, i]
    gap_frac = col.count("-") / nseq

    if gap_frac >= args.gap_threshold:
        removed_gap_sites += 1
        continue

    residues = [c.upper() for c in col if c != "-"]
    if not residues:
        continue

    counts = Counter(residues)
    aa, count = counts.most_common(1)[0]
    freq = count / len(residues)

    if freq >= args.threshold:
        consensus.append(aa)
    else:
        ambiguous_sites += 1
        consensus.append(aa.lower() if args.lowercase_ambiguous else aa)

seq = "".join(consensus)

with open(args.output, "w") as f:
    f.write(f">{args.name} length={len(seq)} threshold={args.threshold} gap_threshold={args.gap_threshold} ambiguous_sites={ambiguous_sites} removed_gap_sites={removed_gap_sites}\n")
    for line in textwrap.wrap(seq, 60):
        f.write(line + "\n")

print(f"[DONE] {args.output}")
print(f"[INFO] consensus length: {len(seq)}")
print(f"[INFO] ambiguous sites below threshold: {ambiguous_sites}")
print(f"[INFO] removed high-gap columns: {removed_gap_sites}")
