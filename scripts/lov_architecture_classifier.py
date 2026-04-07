#!/usr/bin/env python3
from collections import defaultdict
import argparse
import csv
import re

MIN_IEVALUE = 1e-3

PAS_DOMAINS = {
    "PAS", "PAS_3", "PAS_4", "PAS_8", "PAS_9", "PAS_11", "PAS_12", "PAS_13", "PAC"
}

KINASE_EUK_DOMAINS = {
    "Pkinase", "PK_Tyr_Ser-Thr", "Kinase-like"
}

HK_CORE_DOMAINS = {
    "HisKA", "HisKA_2", "HWE_HK", "HATPase_c", "HATPase_c_2", "Response_reg"
}

STAS_DOMAINS = {
    "STAS"
}

HTH_DOMAINS = {
    "HTH", "LuxR", "GerE", "HTH_3", "TetR_N", "AraC_binding"
}

BAD_PATTERNS = {
    "Dioxygenase", "dioxygenase",
    "dehydrogenase", "oxidoreductase",
    "Monooxygenase", "monooxygenase",
    "amine_oxidase", "FAD_binding",
    "FAD_oxidase", "NAD_binding",
    "Fer2", "Pyr_redox", "Aldedh"
}

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

def parse_domtblout(path):
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

            if i_evalue is None or i_evalue > MIN_IEVALUE:
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
                if current_id is not None:
                    seqs[current_id] = "".join(current_lines)
                current_id = line[1:].split()[0]
                current_lines = []
            else:
                current_lines.append(line.strip())

        if current_id is not None:
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

def classify_protein(hits):
    domain_names = [h["domain"] for h in hits]
    domain_set = set(domain_names)

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
            if bad.lower() in text:
                has_bad = True
                bad_hits.append(bad)

    # subtype 분류
    subtype = "ambiguous"

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
    score = 0
    if has_pas:
        score += 5
    if has_euk_kinase:
        score += 6
    if has_hk:
        score += 6
    if has_stas:
        score += 5
    if has_hth:
        score += 4
    if has_bad:
        score -= 10

    # decision
    if not has_pas:
        decision = "exclude"
    elif has_bad and score < 5:
        decision = "exclude"
    elif subtype in {"lov_kinase_euk", "lov_hk_bact", "lov_stas", "lov_hth"}:
        decision = "keep"
    elif subtype == "lov_pas_only":
        decision = "review"
    else:
        decision = "review"

    reasons = []
    for h in hits:
        d = h["domain"]
        if d in PAS_DOMAINS:
            reasons.append(f"{d}[{h['ali_from']}-{h['ali_to']}]: PAS-like")
        elif d in KINASE_EUK_DOMAINS:
            reasons.append(f"{d}[{h['ali_from']}-{h['ali_to']}]: euk_kinase")
        elif d in HK_CORE_DOMAINS:
            reasons.append(f"{d}[{h['ali_from']}-{h['ali_to']}]: hk_core")
        elif d in STAS_DOMAINS:
            reasons.append(f"{d}[{h['ali_from']}-{h['ali_to']}]: STAS")
        elif d in HTH_DOMAINS:
            reasons.append(f"{d}[{h['ali_from']}-{h['ali_to']}]: HTH")

    if has_bad:
        reasons.append("bad_hits=" + ",".join(sorted(set(bad_hits))))

    return {
        "decision": decision,
        "subtype": subtype,
        "score": score,
        "has_pas": has_pas,
        "has_euk_kinase": has_euk_kinase,
        "has_hk": has_hk,
        "has_stas": has_stas,
        "has_hth": has_hth,
        "has_bad": has_bad,
        "domains": "; ".join(sorted(domain_set)),
        "reasons": " | ".join(reasons)
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--domtblout", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--outprefix", default="lov_arch")
    args = parser.parse_args()

    proteins = parse_domtblout(args.domtblout)
    seqs = read_fasta(args.fasta)

    rows = []
    keep_ids, review_ids, exclude_ids = [], [], []
    subtype_to_ids = defaultdict(list)

    for pid in seqs:
        hits = proteins.get(pid, [])
        if not hits:
            result = {
                "decision": "exclude",
                "subtype": "no_hit",
                "score": 0,
                "has_pas": False,
                "has_euk_kinase": False,
                "has_hk": False,
                "has_stas": False,
                "has_hth": False,
                "has_bad": False,
                "domains": "",
                "reasons": "no significant domain hits"
            }
        else:
            result = classify_protein(hits)

        row = {
            "protein_id": pid,
            "decision": result["decision"],
            "subtype": result["subtype"],
            "score": result["score"],
            "has_pas": result["has_pas"],
            "has_euk_kinase": result["has_euk_kinase"],
            "has_hk": result["has_hk"],
            "has_stas": result["has_stas"],
            "has_hth": result["has_hth"],
            "has_bad": result["has_bad"],
            "domains": result["domains"],
            "reasons": result["reasons"],
        }
        rows.append(row)
        subtype_to_ids[result["subtype"]].append(pid)

        if result["decision"] == "keep":
            keep_ids.append(pid)
        elif result["decision"] == "review":
            review_ids.append(pid)
        else:
            exclude_ids.append(pid)

    rows.sort(key=lambda x: (x["decision"], x["subtype"], x["score"]), reverse=True)

    summary_file = f"{args.outprefix}.summary.tsv"
    with open(summary_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "protein_id", "decision", "subtype", "score",
                "has_pas", "has_euk_kinase", "has_hk", "has_stas", "has_hth", "has_bad",
                "domains", "reasons"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    # decision별 fasta
    for fname, ids in [
        (f"{args.outprefix}.keep.fasta", keep_ids),
        (f"{args.outprefix}.review.fasta", review_ids),
        (f"{args.outprefix}.exclude_ids.txt", exclude_ids),
    ]:
        if fname.endswith(".fasta"):
            write_fasta(fname, ids, seqs)
        else:
            with open(fname, "w", encoding="utf-8") as f:
                for x in ids:
                    f.write(x + "\n")

    # subtype별 fasta
    for subtype, ids in subtype_to_ids.items():
        if subtype in {"lov_kinase_euk", "lov_hk_bact", "lov_stas", "lov_hth", "lov_pas_only"}:
            write_fasta(f"{args.outprefix}.{subtype}.fasta", ids, seqs)

    print(f"[DONE] {summary_file}")
    print(f"[KEEP]   {len(keep_ids)}")
    print(f"[REVIEW] {len(review_ids)}")
    print(f"[EXCLUDE] {len(exclude_ids)}")
    print("[SUBTYPE COUNTS]")
    for subtype, ids in sorted(subtype_to_ids.items(), key=lambda x: len(x[1]), reverse=True):
        print(f"  {subtype}: {len(ids)}")

if __name__ == "__main__":
    main()
