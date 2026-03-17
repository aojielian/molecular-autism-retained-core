#!/usr/bin/env python3
import argparse
import json
import math
import os
import re
import time
from typing import List, Tuple

import pandas as pd
import requests
import matplotlib.pyplot as plt

ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"


def log(msg: str):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def safe_mkdir(path: str):
    os.makedirs(path, exist_ok=True)


def read_any_table(path: str) -> pd.DataFrame:
    if path.endswith(".csv") or path.endswith(".csv.gz"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def write_tsv(df: pd.DataFrame, path: str):
    df.to_csv(path, sep="\t", index=False)


def pick_col(cols: List[str], candidates: List[str]) -> str:
    lower_map = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def extract_from_frozen_signature(signature_tsv: str) -> Tuple[List[str], List[str], pd.DataFrame]:
    df = read_any_table(signature_tsv)
    if not {"direction", "gene"}.issubset(df.columns):
        raise ValueError(f"{signature_tsv} must contain columns: direction, gene")

    cand = (
        df.loc[df["direction"].astype(str) == "disc_cluster2_up", "gene"]
        .astype(str).str.upper().str.strip().drop_duplicates().tolist()
    )
    ref = (
        df.loc[df["direction"].astype(str) == "disc_cluster0_up", "gene"]
        .astype(str).str.upper().str.strip().drop_duplicates().tolist()
    )
    return cand, ref, df


def extract_top_de_genes(de_path: str, n_up: int = 150, n_down: int = 150) -> Tuple[List[str], List[str], pd.DataFrame]:
    df = read_any_table(de_path)
    gene_col = pick_col(df.columns.tolist(), ["gene", "Gene", "symbol", "SYMBOL"])
    logfc_col = pick_col(df.columns.tolist(), ["logFC", "avg_log2FC", "avg_logFC", "log2FoldChange"])
    p_col = pick_col(df.columns.tolist(), ["FDR", "adj.P.Val", "padj", "p_val_adj", "PValue", "P.Value", "pvalue"])

    if gene_col is None or logfc_col is None or p_col is None:
        raise ValueError(
            f"Could not identify gene/logFC/P columns in {de_path}. Columns: {list(df.columns)}"
        )

    tmp = df[[gene_col, logfc_col, p_col]].copy()
    tmp.columns = ["gene", "logFC", "priority_p"]
    tmp["gene"] = tmp["gene"].astype(str).str.upper().str.strip()
    tmp = tmp.replace([float("inf"), -float("inf")], pd.NA).dropna()
    tmp = tmp.sort_values(["priority_p", "logFC"], ascending=[True, False])

    up = (
        tmp[tmp["logFC"] > 0]
        .drop_duplicates(subset="gene")
        .head(n_up)["gene"]
        .tolist()
    )
    down = (
        tmp[tmp["logFC"] < 0]
        .sort_values(["priority_p", "logFC"], ascending=[True, True])
        .drop_duplicates(subset="gene")
        .head(n_down)["gene"]
        .tolist()
    )

    return up, down, tmp


def enrichr_add_list(genes: List[str], description: str) -> str:
    payload = {
        "list": (None, "\n".join(genes)),
        "description": (None, description)
    }
    r = requests.post(ENRICHR_ADD_URL, files=payload, timeout=120)
    r.raise_for_status()
    js = r.json()
    if "userListId" not in js:
        raise RuntimeError(f"Unexpected addList response: {js}")
    return str(js["userListId"])


def parse_enrichr_rows(rows: list, library: str) -> pd.DataFrame:
    parsed = []
    for row in rows:
        parsed.append({
            "Rank": row[0] if len(row) > 0 else None,
            "Term": row[1] if len(row) > 1 else None,
            "P_value": row[2] if len(row) > 2 else None,
            "Z_score": row[3] if len(row) > 3 else None,
            "Combined_score": row[4] if len(row) > 4 else None,
            "Overlapping_genes": ";".join(row[5]) if len(row) > 5 and isinstance(row[5], list) else (row[5] if len(row) > 5 else None),
            "Adjusted_p_value": row[6] if len(row) > 6 else None,
            "Old_p_value": row[7] if len(row) > 7 else None,
            "Old_adjusted_p_value": row[8] if len(row) > 8 else None,
            "Library": library
        })
    df = pd.DataFrame(parsed)
    if not df.empty:
        for c in ["P_value", "Z_score", "Combined_score", "Adjusted_p_value", "Old_p_value", "Old_adjusted_p_value"]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def enrichr_query(genes: List[str], library: str, description: str) -> pd.DataFrame:
    user_list_id = enrichr_add_list(genes, description)
    params = {"userListId": user_list_id, "backgroundType": library}
    r = requests.get(ENRICHR_ENRICH_URL, params=params, timeout=180)
    r.raise_for_status()
    js = r.json()

    if library in js:
        rows = js[library]
    elif len(js) == 1:
        rows = list(js.values())[0]
    else:
        raise RuntimeError(f"Unexpected enrich response for {library}: {js}")

    return parse_enrichr_rows(rows, library)


STOPWORDS = {
    "up", "down", "vs", "and", "or", "of", "in", "on", "at", "after", "before",
    "treated", "treatment", "vehicle", "control", "compound", "drug", "cells", "cell",
    "human", "mouse", "mice", "rat", "rats", "line", "lines", "hour", "hours", "day", "days",
    "week", "weeks", "dose", "doses", "mg", "ug", "um", "nm", "mm", "gse", "gsm", "geo"
}


def normalize_perturbagen(term: str) -> str:
    if pd.isna(term):
        return ""
    x = str(term).lower()
    x = re.sub(r"gse\d+", " ", x)
    x = re.sub(r"gsm\d+", " ", x)
    x = re.sub(r"pmid\d+", " ", x)
    x = re.sub(r"[^a-z0-9\s\-_\/]", " ", x)
    x = re.sub(r"[_/]", " ", x)
    x = re.sub(r"\s+", " ", x).strip()

    tokens = x.split()
    kept = []
    for tok in tokens:
        if tok in STOPWORDS:
            break
        if re.fullmatch(r"\d+(\.\d+)?", tok):
            break
        if tok in {"mcf7", "hela", "hek293", "293t", "u2os", "npc", "ipsc", "neuron", "neurons", "microglia", "astrocyte", "astrocytes"}:
            break
        kept.append(tok)
        if len(kept) >= 4:
            break

    return " ".join(kept).strip()


def add_common_metrics(df: pd.DataFrame, side: str) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out[f"neglog10_adjP_{side}"] = -out["Adjusted_p_value"].clip(lower=1e-300).apply(math.log10)
    out[f"neglog10_P_{side}"] = -out["P_value"].clip(lower=1e-300).apply(math.log10)
    out["perturbagen_norm"] = out["Term"].apply(normalize_perturbagen)
    out = out[out["perturbagen_norm"] != ""].copy()
    return out


def combine_reversal(down_df: pd.DataFrame, up_df: pd.DataFrame, top_n: int = 200) -> pd.DataFrame:
    if down_df.empty or up_df.empty:
        return pd.DataFrame()

    d = add_common_metrics(down_df, "down")
    u = add_common_metrics(up_df, "up")

    d2 = d.sort_values(["Adjusted_p_value", "P_value", "Combined_score"], ascending=[True, True, False]).drop_duplicates("perturbagen_norm")
    u2 = u.sort_values(["Adjusted_p_value", "P_value", "Combined_score"], ascending=[True, True, False]).drop_duplicates("perturbagen_norm")

    merged = d2.merge(
        u2,
        on="perturbagen_norm",
        how="inner",
        suffixes=("_candDown", "_refUp")
    )

    if merged.empty:
        return merged

    merged["reversal_score"] = (
        merged["neglog10_adjP_down"] +
        merged["neglog10_adjP_up"] +
        0.01 * merged["Combined_score_candDown"].fillna(0) +
        0.01 * merged["Combined_score_refUp"].fillna(0)
    )
    merged["support_strength"] = merged["neglog10_adjP_down"] + merged["neglog10_adjP_up"]

    merged = merged.sort_values(["reversal_score", "support_strength"], ascending=[False, False]).reset_index(drop=True)
    return merged.head(top_n)


def run_dsigdb(genes: List[str], label: str) -> pd.DataFrame:
    try:
        df = enrichr_query(genes, "DSigDB", label)
        return df.sort_values(["Adjusted_p_value", "P_value", "Combined_score"], ascending=[True, True, False])
    except Exception as e:
        log(f"WARNING: DSigDB query failed for {label}: {e}")
        return pd.DataFrame()


def make_barplot(df: pd.DataFrame, title: str, out_png: str, out_pdf: str, label_col: str = "perturbagen_norm"):
    if df.empty:
        return
    plot_df = df.head(15).copy().iloc[::-1]

    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(plot_df) + 1.5)))
    ax.barh(plot_df[label_col], plot_df["reversal_score"])
    ax.set_xlabel("Integrated reversal score")
    ax.set_ylabel("")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def make_scatter(df: pd.DataFrame, title: str, out_png: str, out_pdf: str):
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    ax.scatter(df["neglog10_adjP_down"], df["neglog10_adjP_up"], s=25, alpha=0.8)

    top = df.head(min(10, len(df)))
    for _, r in top.iterrows():
        ax.text(r["neglog10_adjP_down"], r["neglog10_adjP_up"], str(r["perturbagen_norm"]), fontsize=8)

    ax.set_xlabel("Candidate genes enriched in GEO drug-down signatures\n(-log10 adjusted P)")
    ax.set_ylabel("Reference genes enriched in GEO drug-up signatures\n(-log10 adjusted P)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def summarise_run(primary_combined: pd.DataFrame, sens_combined: pd.DataFrame, out_txt: str):
    lines = []
    lines.append("Perturbation-based prioritization summary")
    lines.append("======================================")
    lines.append("")

    if not primary_combined.empty:
        lines.append("Top primary reversal candidates:")
        for _, r in primary_combined.head(10).iterrows():
            lines.append(
                f"- {r['perturbagen_norm']}: reversal_score={r['reversal_score']:.3f}; "
                f"candDown_adjP={r['Adjusted_p_value_candDown']:.3g}; "
                f"refUp_adjP={r['Adjusted_p_value_refUp']:.3g}"
            )
        lines.append("")

    if not sens_combined.empty:
        lines.append("Top pooled-DE sensitivity reversal candidates:")
        for _, r in sens_combined.head(10).iterrows():
            lines.append(
                f"- {r['perturbagen_norm']}: reversal_score={r['reversal_score']:.3f}; "
                f"pooledUp->drugDown_adjP={r['Adjusted_p_value_candDown']:.3g}; "
                f"pooledDown->drugUp_adjP={r['Adjusted_p_value_refUp']:.3g}"
            )
        lines.append("")

    with open(out_txt, "w") as f:
        f.write("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description="Perturbation-based prioritization for the candidate/reference remodeling axis")
    ap.add_argument(
        "--signature_tsv",
        default="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package4_CrossCohortValidation_Velmeshev_v2/08_discovery_frozen_state_signature.tsv"
    )
    ap.add_argument(
        "--pooled_de_tsv",
        default="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation/13_pseudobulk_all_retained_ASD_vs_Control.tsv.gz"
    )
    ap.add_argument(
        "--outdir",
        default="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package9_PerturbationPrioritization"
    )
    ap.add_argument("--top_primary", type=int, default=100)
    ap.add_argument("--top_sensitivity", type=int, default=150)
    args = ap.parse_args()

    safe_mkdir(args.outdir)
    safe_mkdir(os.path.join(args.outdir, "tables"))
    safe_mkdir(os.path.join(args.outdir, "plots"))

    log("Reading frozen discovery signature ...")
    cand, ref, _ = extract_from_frozen_signature(args.signature_tsv)
    cand = cand[:args.top_primary]
    ref = ref[:args.top_primary]

    with open(os.path.join(args.outdir, "tables", "01_primary_candidate_genes.txt"), "w") as f:
        f.write("\n".join(cand) + "\n")
    with open(os.path.join(args.outdir, "tables", "02_primary_reference_genes.txt"), "w") as f:
        f.write("\n".join(ref) + "\n")

    log(f"Primary candidate genes: {len(cand)}; reference genes: {len(ref)}")

    log("Querying Enrichr: candidate genes vs Drug_Perturbations_from_GEO_down ...")
    primary_geo_down = enrichr_query(cand, "Drug_Perturbations_from_GEO_down", "candidate_genes_vs_GEO_drug_down")
    write_tsv(primary_geo_down, os.path.join(args.outdir, "tables", "03_primary_candidate_vs_GEO_drug_down.tsv"))

    log("Querying Enrichr: reference genes vs Drug_Perturbations_from_GEO_up ...")
    primary_geo_up = enrichr_query(ref, "Drug_Perturbations_from_GEO_up", "reference_genes_vs_GEO_drug_up")
    write_tsv(primary_geo_up, os.path.join(args.outdir, "tables", "04_primary_reference_vs_GEO_drug_up.tsv"))

    log("Combining directional reversal support ...")
    primary_combined = combine_reversal(primary_geo_down, primary_geo_up, top_n=200)
    write_tsv(primary_combined, os.path.join(args.outdir, "tables", "05_primary_reversal_combined.tsv"))

    log("Querying Enrichr: DSigDB for candidate genes ...")
    dsig_cand = run_dsigdb(cand, "candidate_genes_DSigDB")
    write_tsv(dsig_cand, os.path.join(args.outdir, "tables", "06_primary_candidate_DSigDB.tsv"))

    log("Querying Enrichr: DSigDB for reference genes ...")
    dsig_ref = run_dsigdb(ref, "reference_genes_DSigDB")
    write_tsv(dsig_ref, os.path.join(args.outdir, "tables", "07_primary_reference_DSigDB.tsv"))

    log("Reading pooled retained microglia DE for sensitivity analysis ...")
    up2, down2, _ = extract_top_de_genes(args.pooled_de_tsv, n_up=args.top_sensitivity, n_down=args.top_sensitivity)
    with open(os.path.join(args.outdir, "tables", "08_sensitivity_pooled_up_genes.txt"), "w") as f:
        f.write("\n".join(up2) + "\n")
    with open(os.path.join(args.outdir, "tables", "09_sensitivity_pooled_down_genes.txt"), "w") as f:
        f.write("\n".join(down2) + "\n")

    log("Querying Enrichr: pooled ASD-up genes vs Drug_Perturbations_from_GEO_down ...")
    sens_geo_down = enrichr_query(up2, "Drug_Perturbations_from_GEO_down", "pooled_ASD_up_vs_GEO_drug_down")
    write_tsv(sens_geo_down, os.path.join(args.outdir, "tables", "10_sensitivity_pooled_up_vs_GEO_drug_down.tsv"))

    log("Querying Enrichr: pooled ASD-down genes vs Drug_Perturbations_from_GEO_up ...")
    sens_geo_up = enrichr_query(down2, "Drug_Perturbations_from_GEO_up", "pooled_ASD_down_vs_GEO_drug_up")
    write_tsv(sens_geo_up, os.path.join(args.outdir, "tables", "11_sensitivity_pooled_down_vs_GEO_drug_up.tsv"))

    sens_combined = combine_reversal(sens_geo_down, sens_geo_up, top_n=200)
    write_tsv(sens_combined, os.path.join(args.outdir, "tables", "12_sensitivity_reversal_combined.tsv"))

    primary_set = set(primary_combined.get("perturbagen_norm", pd.Series(dtype=str)).dropna().astype(str))
    sens_set = set(sens_combined.get("perturbagen_norm", pd.Series(dtype=str)).dropna().astype(str))
    overlap = sorted(primary_set & sens_set)
    overlap_df = pd.DataFrame({"shared_perturbagen": overlap})
    write_tsv(overlap_df, os.path.join(args.outdir, "tables", "13_primary_vs_sensitivity_overlap.tsv"))

    log("Making plots ...")
    make_barplot(
        primary_combined,
        "Primary perturbation reversal prioritization\n(candidate genes -> GEO drug-down; reference genes -> GEO drug-up)",
        os.path.join(args.outdir, "plots", "01_primary_reversal_barplot.png"),
        os.path.join(args.outdir, "plots", "01_primary_reversal_barplot.pdf")
    )

    make_scatter(
        primary_combined,
        "Primary directional reversal support",
        os.path.join(args.outdir, "plots", "02_primary_reversal_scatter.png"),
        os.path.join(args.outdir, "plots", "02_primary_reversal_scatter.pdf")
    )

    make_barplot(
        sens_combined,
        "Sensitivity perturbation reversal prioritization\n(pooled retained ASD DE)",
        os.path.join(args.outdir, "plots", "03_sensitivity_reversal_barplot.png"),
        os.path.join(args.outdir, "plots", "03_sensitivity_reversal_barplot.pdf")
    )

    summarise_run(
        primary_combined,
        sens_combined,
        os.path.join(args.outdir, "14_summary.txt")
    )

    params = {
        "signature_tsv": args.signature_tsv,
        "pooled_de_tsv": args.pooled_de_tsv,
        "outdir": args.outdir,
        "top_primary": args.top_primary,
        "top_sensitivity": args.top_sensitivity,
        "libraries": [
            "Drug_Perturbations_from_GEO_down",
            "Drug_Perturbations_from_GEO_up",
            "DSigDB"
        ]
    }
    with open(os.path.join(args.outdir, "15_run_parameters.json"), "w") as f:
        json.dump(params, f, indent=2)

    log("Done.")


if __name__ == "__main__":
    main()
