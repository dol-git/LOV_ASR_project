#!/usr/bin/env python3
from collections import defaultdict
import argparse
import csv

# -----------------------------
# Default thresholds
# -----------------------------
DEFAULT_IEVALUE = 1e-3

DECISION_PRIORITY = {"keep": 2, "review": 1, "exclude": 0}

PAS_DOMAINS = {
    "PAS", "PAS_3", "PAS_4", "PAS_8", "PAS_9", "PAS_11", "PAS_12", "PAS_13", "PAC"
}

KINASE_EUK_DOMAINS = {
    "Pkinase", "PK_Tyr_Ser-Thr", "Kinase-like"
}

HK_CORE_DOMAINS = {
    "HisKA", "HisKA_2", "HWE_HK", "HATPase_c", "HATPase_c_2", "Response_reg"
}

STAS_DOMAINS = {"STAS"}

HTH_DOMAINS = {
    "HTH", "LuxR", "GerE", "HTH_3", "TetR_N", "AraC_binding"
}

BAD_PATTERNS = {
    "dioxygenase", "dehydrogenase", "oxidoreductase",
    "monooxygenase", "amine_oxidase",
    "fad_binding", "fad_oxidase", "nad_binding",
    "fer2", "pyr_redox", "aldedh"
}


# -----------------------------
# Utility
# -----------------------------
def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None


def safe_int(x):
    try:
        return int(x)
    except Exception:
        return None


# -----------------------------
# Parsing
# -----------------------------
def parse_domtblout(path, min_ievalue):
    proteins = defaultdict(list)

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split()
            if len(parts) < 23:
                continue

            target_name = parts[0]
            query_name = parts[3]
            i_evalue = safe_float(parts[12])
            ali_from = safe_int(parts[17])
            ali_to = safe_int(parts[18])
            description = " ".join(parts[22:])

            if i_evalue is None or i_evalue > min_ievalue:
                continue

            proteins[query_name].append({
                "domain": target_name,
                "i_evalue": i_evalue,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "description": description
            })

    return proteins


def read_fasta(path):
    seqs = {}
    current_id = None
    current_lines = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_lines)
                current_id = line[1:].split()[0]
                current_lines = []
            else:
                current_lines.append(line.strip())

        if current_id:
            seqs[current_id] = "".join(current_lines)

    return seqs


def write_fasta(path, ids, seqs):
    with open(path, "w", encoding="utf-8") as f:
        for seq_id in ids:
            if seq_id not in seqs:
                continue

            f.write(f">{seq_id}\n")
            seq = seqs[seq_id]
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


# -----------------------------
# Classification
# -----------------------------
def classify_protein(hits):
    domain_set = {h["domain"] for h in hits}

    has_pas = any(d in PAS_DOMAINS for d in domain_set)
    has_euk_kinase = any(d in KINASE_EUK_DOMAINS for d in domain_set)
    has_hk = any(d in HK_CORE_DOMAINS for d in domain_set)
    has_stas = any(d in STAS_DOMAINS for d in domain_set)
    has_hth = any(d in HTH_DOMAINS for d in domain_set)

    has_bad = False
    bad_hits = []

    for h in hits:
        text = f"{h['domain']} {h['description']}".lower()
        for bad in BAD_PATTERNS:
            if bad in text:
                has_bad = True
                bad_hits.append(bad)

    # subtype
    if has_pas and has_euk_kinase:
        subtype = "lov_kinase_euk"
    elif has_pas and has_hk:
        subtype = "lov_hk_bact"
    elif has_pas and has_stas:
        subtype = "lov_stas"
    elif has_pas and has_hth:
        subtype = "lov_hth"
    elif has_pas:
        subtype = "lov_pas_only"
    else:
        subtype = "non_lov_like"

    # score
    score = (
        has_pas * 5 +
        has_euk_kinase * 6 +
        has_hk * 6 +
        has_stas * 5 +
        has_hth * 4 -
        has_bad * 10
    )

    # decision
    if not has_pas:
        decision = "exclude"
    elif has_bad and score < 5:
        decision = "exclude"
    elif subtype in {"lov_kinase_euk", "lov_hk_bact", "lov_stas", "lov_hth"}:
        decision = "keep"
    else:
        decision = "review"

    return {
        "decision": decision,
        "subtype": subtype,
        "score": score,
        "domains": ";".join(sorted(domain_set)),
        "has_bad": has_bad,
        "bad_hits": ",".join(set(bad_hits))
    }


# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--domtblout", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--outprefix", default="lov_arch")
    parser.add_argument("--min_ievalue", type=float, default=DEFAULT_IEVALUE)

    args = parser.parse_args()

    proteins = parse_domtblout(args.domtblout, args.min_ievalue)
    seqs = read_fasta(args.fasta)

    rows = []
    decision_groups = defaultdict(list)
    subtype_groups = defaultdict(list)

    for pid in seqs:
        hits = proteins.get(pid, [])

        if hits:
            result = classify_protein(hits)
        else:
            result = {
                "decision": "exclude",
                "subtype": "no_hit",
                "score": 0,
                "domains": "",
                "has_bad": False,
                "bad_hits": ""
            }

        row = {"protein_id": pid, **result}
        rows.append(row)

        decision_groups[result["decision"]].append(pid)
        subtype_groups[result["subtype"]].append(pid)

    # 🔥 정렬 개선
    rows.sort(
        key=lambda x: (
            DECISION_PRIORITY[x["decision"]],
            x["score"]
        ),
        reverse=True
    )

    # summary
    summary_file = f"{args.outprefix}.summary.tsv"
    with open(summary_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # decision FASTA
    for decision, ids in decision_groups.items():
        if decision != "exclude":
            write_fasta(f"{args.outprefix}.{decision}.fasta", ids, seqs)

    # subtype FASTA (exclude 제거)
    for subtype, ids in subtype_groups.items():
        valid_ids = [
            pid for pid in ids
            if next(r for r in rows if r["protein_id"] == pid)["decision"] != "exclude"
        ]
        if valid_ids:
            write_fasta(f"{args.outprefix}.{subtype}.fasta", valid_ids, seqs)

    print(f"[DONE] {summary_file}")
    for d in ["keep", "review", "exclude"]:
        print(f"[{d.upper()}] {len(decision_groups[d])}")


if __name__ == "__main__":
    main()