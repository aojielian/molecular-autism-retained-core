#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import math
import os
import time
from typing import List, Tuple

import pandas as pd
import requests
import matplotlib.pyplot as plt

ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"

TF_LIBRARIES = [
    "ChEA_2022",
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
    "TRRUST_Transcription_Factors_2019",
]

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

def parse_enrichr_rows(rows: list, library: str, query_label: str) -> pd.DataFrame:
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
            "Library": library,
            "Query_label": query_label
        })
    df = pd.DataFrame(parsed)
    if not df.empty:
        for c in ["P_value", "Z_score", "Combined_score", "Adjusted_p_value", "Old_p_value", "Old_adjusted_p_value"]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def extract_tf_name(term: str) -> str:
    if pd.isna(term):
        return ""
    x = str(term).strip()
    # 常见 Enrichr TF term 通常首个 token 就是 TF
    x = x.replace("_", " ")
    tf = x.split()[0].upper().strip()
    tf = tf.replace(",", "").replace(";", "")
    return tf

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

    out = parse_enrichr_rows(rows, library, description)
    if not out.empty:
        out["TF"] = out["Term"].apply(extract_tf_name)
        out["neglog10_adjP"] = -out["Adjusted_p_value"].clip(lower=1e-300).apply(math.log10)
    return out

def run_tf_libraries(genes: List[str], label: str) -> pd.DataFrame:
    all_res = []
    for lib in TF_LIBRARIES:
        try:
            log(f"Querying {label} -> {lib}")
            df = enrichr_query(genes, lib, label)
            if not df.empty:
                all_res.append(df)
        except Exception as e:
            log(f"WARNING: failed for {label} / {lib}: {e}")
    if len(all_res) == 0:
        return pd.DataFrame()
    return pd.concat(all_res, axis=0, ignore_index=True)

def summarize_tf_hits(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[
            "TF", "n_libraries_supported", "best_adjusted_p", "best_combined_score",
            "supporting_libraries", "best_term", "aggregate_score"
        ])

    tmp = df.copy()
    tmp = tmp[tmp["TF"].astype(str) != ""].copy()

    rows = []
    for tf, sub in tmp.groupby("TF"):
        sub2 = sub.sort_values(["Adjusted_p_value", "Combined_score"], ascending=[True, False])
        rows.append({
            "TF": tf,
            "n_libraries_supported": sub["Library"].nunique(),
            "best_adjusted_p": sub["Adjusted_p_value"].min(),
            "best_combined_score": sub["Combined_score"].max(),
            "supporting_libraries": ";".join(sorted(sub["Library"].dropna().astype(str).unique().tolist())),
            "best_term": sub2.iloc[0]["Term"],
            "aggregate_score": sub["neglog10_adjP"].sum()
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["n_libraries_supported", "best_adjusted_p", "aggregate_score", "best_combined_score"],
        ascending=[False, True, False, False]
    ).reset_index(drop=True)
    return out

def make_convergence_table(df1: pd.DataFrame, df2: pd.DataFrame, label1: str, label2: str) -> pd.DataFrame:
    if df1.empty or df2.empty:
        return pd.DataFrame()

    a = df1[["TF", "n_libraries_supported", "best_adjusted_p", "aggregate_score"]].copy()
    b = df2[["TF", "n_libraries_supported", "best_adjusted_p", "aggregate_score"]].copy()

    a.columns = ["TF", f"{label1}_n_libs", f"{label1}_best_adjP", f"{label1}_aggregate_score"]
    b.columns = ["TF", f"{label2}_n_libs", f"{label2}_best_adjP", f"{label2}_aggregate_score"]

    m = a.merge(b, on="TF", how="inner")
    if m.empty:
        return m

    m["joint_score"] = (
        m[f"{label1}_aggregate_score"].fillna(0) +
        m[f"{label2}_aggregate_score"].fillna(0)
    )
    m = m.sort_values(
        ["joint_score", f"{label1}_best_adjP", f"{label2}_best_adjP"],
        ascending=[False, True, True]
    ).reset_index(drop=True)
    return m

def make_barplot(df: pd.DataFrame, title: str, out_png: str, out_pdf: str, n: int = 15):
    if df.empty:
        return
    plot_df = df.head(n).copy().iloc[::-1]

    fig, ax = plt.subplots(figsize=(7.8, max(4.5, 0.35 * len(plot_df) + 1.2)))
    ax.barh(plot_df["TF"], plot_df["aggregate_score"])
    ax.set_xlabel("Aggregate TF prioritization score")
    ax.set_ylabel("")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

def summarise_run(
    cand_sum: pd.DataFrame,
    ref_sum: pd.DataFrame,
    up_sum: pd.DataFrame,
    down_sum: pd.DataFrame,
    conv_up: pd.DataFrame,
    conv_down: pd.DataFrame,
    outfile: str
):
    lines = []
    lines.append("TF prioritization summary")
    lines.append("========================")
    lines.append("")

    def add_top_block(title: str, df: pd.DataFrame):
        lines.append(title)
        if df.empty:
            lines.append("- No supported TF hits.")
        else:
            for _, r in df.head(10).iterrows():
                lines.append(
                    f"- {r['TF']}: n_libs={r['n_libraries_supported']}; "
                    f"best_adjP={r['best_adjusted_p']:.3g}; "
                    f"aggregate_score={r['aggregate_score']:.3f}; "
                    f"libs={r['supporting_libraries']}"
                )
        lines.append("")

    add_top_block("Top TFs for primary candidate program:", cand_sum)
    add_top_block("Top TFs for primary reference program:", ref_sum)
    add_top_block("Top TFs for pooled ASD-up genes:", up_sum)
    add_top_block("Top TFs for pooled ASD-down genes:", down_sum)

    lines.append("Convergent TFs: candidate program + pooled ASD-up genes")
    if conv_up.empty:
        lines.append("- No convergent TFs.")
    else:
        for _, r in conv_up.head(15).iterrows():
            lines.append(
                f"- {r['TF']}: joint_score={r['joint_score']:.3f}; "
                f"candidate_adjP={r['candidate_best_adjP']:.3g}; "
                f"pooledUp_adjP={r['pooledUp_best_adjP']:.3g}"
            )
    lines.append("")

    lines.append("Convergent TFs: reference program + pooled ASD-down genes")
    if conv_down.empty:
        lines.append("- No convergent TFs.")
    else:
        for _, r in conv_down.head(15).iterrows():
            lines.append(
                f"- {r['TF']}: joint_score={r['joint_score']:.3f}; "
                f"reference_adjP={r['reference_best_adjP']:.3g}; "
                f"pooledDown_adjP={r['pooledDown_best_adjP']:.3g}"
            )

    with open(outfile, "w") as f:
        f.write("\n".join(lines))


def env_or_none(name: str):
    value = os.environ.get(name, "").strip()
    return value or None


def project_default(*parts: str):
    project_dir = env_or_none("MA_PROJECT_DIR")
    if project_dir is None:
        return None
    return os.path.join(project_dir, *parts)

def main():
    ap = argparse.ArgumentParser(description="TF-target prioritization for the candidate/reference remodeling axis")
    ap.add_argument(
        "--signature_tsv",
        default=env_or_none("MA_FROZEN_SIGNATURE_TSV") or project_default("data", "08_discovery_frozen_state_signature.tsv")
    )
    ap.add_argument(
        "--pooled_de_tsv",
        default=env_or_none("MA_POOLED_DE_TSV") or project_default("data", "13_pseudobulk_all_retained_ASD_vs_Control.tsv.gz")
    )
    ap.add_argument(
        "--outdir",
        default=env_or_none("MA_TF_OUTDIR") or project_default("results", "Package10_TF_prioritization")
    )
    ap.add_argument("--top_primary", type=int, default=100)
    ap.add_argument("--top_sensitivity", type=int, default=150)
    args = ap.parse_args()

    if args.signature_tsv is None or args.pooled_de_tsv is None or args.outdir is None:
        ap.error("Provide --signature_tsv, --pooled_de_tsv, and --outdir, or set MA_PROJECT_DIR/MA_FROZEN_SIGNATURE_TSV/MA_POOLED_DE_TSV/MA_TF_OUTDIR")

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

    log("Reading pooled retained microglia DE ...")
    up2, down2, _ = extract_top_de_genes(args.pooled_de_tsv, n_up=args.top_sensitivity, n_down=args.top_sensitivity)

    with open(os.path.join(args.outdir, "tables", "03_sensitivity_pooled_up_genes.txt"), "w") as f:
        f.write("\n".join(up2) + "\n")
    with open(os.path.join(args.outdir, "tables", "04_sensitivity_pooled_down_genes.txt"), "w") as f:
        f.write("\n".join(down2) + "\n")

    log("Running TF enrichment for primary candidate genes ...")
    cand_raw = run_tf_libraries(cand, "primary_candidate")
    write_tsv(cand_raw, os.path.join(args.outdir, "tables", "05_primary_candidate_tf_raw.tsv"))
    cand_sum = summarize_tf_hits(cand_raw)
    write_tsv(cand_sum, os.path.join(args.outdir, "tables", "06_primary_candidate_tf_summary.tsv"))

    log("Running TF enrichment for primary reference genes ...")
    ref_raw = run_tf_libraries(ref, "primary_reference")
    write_tsv(ref_raw, os.path.join(args.outdir, "tables", "07_primary_reference_tf_raw.tsv"))
    ref_sum = summarize_tf_hits(ref_raw)
    write_tsv(ref_sum, os.path.join(args.outdir, "tables", "08_primary_reference_tf_summary.tsv"))

    log("Running TF enrichment for pooled ASD-up genes ...")
    up_raw = run_tf_libraries(up2, "pooled_ASD_up")
    write_tsv(up_raw, os.path.join(args.outdir, "tables", "09_pooled_up_tf_raw.tsv"))
    up_sum = summarize_tf_hits(up_raw)
    write_tsv(up_sum, os.path.join(args.outdir, "tables", "10_pooled_up_tf_summary.tsv"))

    log("Running TF enrichment for pooled ASD-down genes ...")
    down_raw = run_tf_libraries(down2, "pooled_ASD_down")
    write_tsv(down_raw, os.path.join(args.outdir, "tables", "11_pooled_down_tf_raw.tsv"))
    down_sum = summarize_tf_hits(down_raw)
    write_tsv(down_sum, os.path.join(args.outdir, "tables", "12_pooled_down_tf_summary.tsv"))

    log("Computing convergence tables ...")
    conv_up = make_convergence_table(cand_sum, up_sum, "candidate", "pooledUp")
    write_tsv(conv_up, os.path.join(args.outdir, "tables", "13_candidate_plus_pooledUp_convergent_TFs.tsv"))

    conv_down = make_convergence_table(ref_sum, down_sum, "reference", "pooledDown")
    write_tsv(conv_down, os.path.join(args.outdir, "tables", "14_reference_plus_pooledDown_convergent_TFs.tsv"))

    log("Making plots ...")
    make_barplot(
        cand_sum,
        "Top TFs for primary candidate program",
        os.path.join(args.outdir, "plots", "01_primary_candidate_topTFs.png"),
        os.path.join(args.outdir, "plots", "01_primary_candidate_topTFs.pdf"),
        n=15
    )
    make_barplot(
        up_sum,
        "Top TFs for pooled retained microglia ASD-up genes",
        os.path.join(args.outdir, "plots", "02_pooled_up_topTFs.png"),
        os.path.join(args.outdir, "plots", "02_pooled_up_topTFs.pdf"),
        n=15
    )
    make_barplot(
        conv_up.rename(columns={"joint_score": "aggregate_score"}),
        "Convergent TFs: candidate program + pooled ASD-up genes",
        os.path.join(args.outdir, "plots", "03_candidate_pooledUp_convergentTFs.png"),
        os.path.join(args.outdir, "plots", "03_candidate_pooledUp_convergentTFs.pdf"),
        n=15
    )

    summarise_run(
        cand_sum, ref_sum, up_sum, down_sum, conv_up, conv_down,
        os.path.join(args.outdir, "15_summary.txt")
    )

    params = {
        "signature_tsv": args.signature_tsv,
        "pooled_de_tsv": args.pooled_de_tsv,
        "outdir": args.outdir,
        "top_primary": args.top_primary,
        "top_sensitivity": args.top_sensitivity,
        "tf_libraries": TF_LIBRARIES
    }
    with open(os.path.join(args.outdir, "16_run_parameters.json"), "w") as f:
        json.dump(params, f, indent=2)

    log("Done.")

if __name__ == "__main__":
    main()
