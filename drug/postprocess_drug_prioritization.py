#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
from collections import Counter

import pandas as pd


def log(msg: str):
    print(msg, flush=True)


def safe_mkdir(path: str):
    os.makedirs(path, exist_ok=True)


def write_tsv(df: pd.DataFrame, path: str):
    df.to_csv(path, sep="\t", index=False)


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def classify_hits(
    df: pd.DataFrame,
    cand_strict: float = 0.05,
    ref_relaxed: float = 0.25,
    bidir_relaxed: float = 0.25
):
    """
    划分三类：
    1) bidirectional: candDown 和 refUp 都 < bidir_relaxed
    2) one-sided candidate suppressor: candDown < cand_strict, refUp >= ref_relaxed
    3) weak/mixed: 其余
    """
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    out = df.copy()

    required = ["Adjusted_p_value_candDown", "Adjusted_p_value_refUp"]
    for c in required:
        if c not in out.columns:
            raise ValueError(f"Missing required column: {c}")

    bidirectional = out[
        (out["Adjusted_p_value_candDown"] < bidir_relaxed) &
        (out["Adjusted_p_value_refUp"] < bidir_relaxed)
    ].copy()

    one_sided = out[
        (out["Adjusted_p_value_candDown"] < cand_strict) &
        (out["Adjusted_p_value_refUp"] >= ref_relaxed)
    ].copy()

    weak_mixed = out.drop(bidirectional.index.union(one_sided.index)).copy()

    return bidirectional, one_sided, weak_mixed


def normalize_text(x: str) -> str:
    x = str(x).lower().strip()
    x = re.sub(r"db\d+", " ", x)
    x = re.sub(r"cid\s*\d+", " ", x)
    x = re.sub(r"\b\d+\b", " ", x)
    x = re.sub(r"[_/]", " ", x)
    x = re.sub(r"[^a-z0-9\s\-]", " ", x)
    x = re.sub(r"\s+", " ", x).strip()
    return x


def mechanism_bucket(term: str) -> str:
    t = normalize_text(term)

    # 先做一些最有生物学意义的规则归类
    if any(k in t for k in ["vemurafenib", "plx4032", "plx4720"]):
        return "MAPK/BRAF-related"
    if any(k in t for k in ["metformin", "pioglitazone", "vitamin d3", "calcitriol", "ascorbic acid", "vitamin c", "creatine", "ubiquinol", "nicotinamide riboside", "resveratrol"]):
        return "Metabolic / nutrient-related"
    if any(k in t for k in ["etanercept", "curcumin", "imatinib", "cetuximab", "trastuzumab", "bexarotene", "formoterol"]):
        return "Immune / signaling-targeted"
    if any(k in t for k in ["lipopolysaccharide", "hypochlorous acid", "hydrogen peroxide", "cigarette smoke", "ethanol", "formaldehyde", "bpde"]):
        return "Inflammatory / stress perturbation"
    if any(k in t for k in ["cisplatin", "decitabine", "5-fluorouracil", "trovafloxacin", "menadione", "4-hydroxytamoxifen"]):
        return "Cytotoxic / DNA-stress-like perturbation"
    if any(k in t for k in ["morphine", "halothane", "isoflurane", "phenytoin", "d-serine"]):
        return "Neuroactive / anesthetic-related"

    return "Other / unclassified"


def make_bucket_table(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["mechanism_bucket", "n_hits", "top_examples"])

    tmp = df.copy()
    source_col = "perturbagen_norm" if "perturbagen_norm" in tmp.columns else "Term_candDown"
    tmp["mechanism_bucket"] = tmp[source_col].apply(mechanism_bucket)

    rows = []
    for bucket, sub in tmp.groupby("mechanism_bucket"):
        examples = sub[source_col].dropna().astype(str).head(5).tolist()
        rows.append({
            "mechanism_bucket": bucket,
            "n_hits": sub.shape[0],
            "top_examples": "; ".join(examples)
        })

    out = pd.DataFrame(rows).sort_values(["n_hits", "mechanism_bucket"], ascending=[False, True]).reset_index(drop=True)
    return out.head(top_n)


def summarize_to_txt(
    primary_bidirectional: pd.DataFrame,
    primary_one_sided: pd.DataFrame,
    sens_bidirectional: pd.DataFrame,
    sens_one_sided: pd.DataFrame,
    primary_bucket: pd.DataFrame,
    sens_bucket: pd.DataFrame,
    outfile: str
):
    lines = []
    lines.append("Postprocessed perturbation prioritization summary")
    lines.append("===============================================")
    lines.append("")

    lines.append(f"Primary bidirectional candidates: {len(primary_bidirectional)}")
    if not primary_bidirectional.empty:
        lines.append("Top primary bidirectional candidates:")
        for _, r in primary_bidirectional.head(10).iterrows():
            lines.append(
                f"- {r.get('perturbagen_norm', 'NA')}: "
                f"candDown_adjP={r.get('Adjusted_p_value_candDown', 'NA'):.3g}; "
                f"refUp_adjP={r.get('Adjusted_p_value_refUp', 'NA'):.3g}; "
                f"reversal_score={r.get('reversal_score', 'NA'):.3f}"
            )
    lines.append("")

    lines.append(f"Primary one-sided candidate suppressors: {len(primary_one_sided)}")
    if not primary_one_sided.empty:
        lines.append("Top primary one-sided suppressors:")
        for _, r in primary_one_sided.head(10).iterrows():
            lines.append(
                f"- {r.get('perturbagen_norm', 'NA')}: "
                f"candDown_adjP={r.get('Adjusted_p_value_candDown', 'NA'):.3g}; "
                f"refUp_adjP={r.get('Adjusted_p_value_refUp', 'NA'):.3g}; "
                f"reversal_score={r.get('reversal_score', 'NA'):.3f}"
            )
    lines.append("")

    lines.append(f"Sensitivity bidirectional candidates: {len(sens_bidirectional)}")
    if not sens_bidirectional.empty:
        lines.append("Top sensitivity bidirectional candidates:")
        for _, r in sens_bidirectional.head(10).iterrows():
            lines.append(
                f"- {r.get('perturbagen_norm', 'NA')}: "
                f"candDown_adjP={r.get('Adjusted_p_value_candDown', 'NA'):.3g}; "
                f"refUp_adjP={r.get('Adjusted_p_value_refUp', 'NA'):.3g}; "
                f"reversal_score={r.get('reversal_score', 'NA'):.3f}"
            )
    lines.append("")

    lines.append(f"Sensitivity one-sided suppressors: {len(sens_one_sided)}")
    if not sens_one_sided.empty:
        lines.append("Top sensitivity one-sided suppressors:")
        for _, r in sens_one_sided.head(10).iterrows():
            lines.append(
                f"- {r.get('perturbagen_norm', 'NA')}: "
                f"candDown_adjP={r.get('Adjusted_p_value_candDown', 'NA'):.3g}; "
                f"refUp_adjP={r.get('Adjusted_p_value_refUp', 'NA'):.3g}; "
                f"reversal_score={r.get('reversal_score', 'NA'):.3f}"
            )
    lines.append("")

    lines.append("Primary mechanism buckets:")
    if not primary_bucket.empty:
        for _, r in primary_bucket.iterrows():
            lines.append(f"- {r['mechanism_bucket']} (n={r['n_hits']}): {r['top_examples']}")
    lines.append("")

    lines.append("Sensitivity mechanism buckets:")
    if not sens_bucket.empty:
        for _, r in sens_bucket.iterrows():
            lines.append(f"- {r['mechanism_bucket']} (n={r['n_hits']}): {r['top_examples']}")

    with open(outfile, "w") as f:
        f.write("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--indir",
        default="/Users/aojie/PROJECTs/molecular_autism/drug/Package9_PerturbationPrioritization_local/tables"
    )
    ap.add_argument(
        "--outdir",
        default="/Users/aojie/PROJECTs/molecular_autism/drug/Package9_PerturbationPrioritization_local/postprocessed"
    )
    ap.add_argument("--cand_strict", type=float, default=0.05)
    ap.add_argument("--ref_relaxed", type=float, default=0.25)
    ap.add_argument("--bidir_relaxed", type=float, default=0.25)
    args = ap.parse_args()

    safe_mkdir(args.outdir)

    primary_path = os.path.join(args.indir, "05_primary_reversal_combined.tsv")
    sens_path = os.path.join(args.indir, "10_sensitivity_reversal_combined.tsv")

    log("Reading primary reversal table ...")
    primary = read_tsv(primary_path)

    log("Reading sensitivity reversal table ...")
    sens = read_tsv(sens_path)

    log("Classifying primary hits ...")
    p_bidir, p_one, p_weak = classify_hits(
        primary,
        cand_strict=args.cand_strict,
        ref_relaxed=args.ref_relaxed,
        bidir_relaxed=args.bidir_relaxed
    )

    log("Classifying sensitivity hits ...")
    s_bidir, s_one, s_weak = classify_hits(
        sens,
        cand_strict=args.cand_strict,
        ref_relaxed=args.ref_relaxed,
        bidir_relaxed=args.bidir_relaxed
    )

    # 排序
    sort_cols = [c for c in ["reversal_score", "support_strength"] if c in primary.columns]
    if sort_cols:
        p_bidir = p_bidir.sort_values(sort_cols, ascending=[False] * len(sort_cols))
        p_one = p_one.sort_values(sort_cols, ascending=[False] * len(sort_cols))
        p_weak = p_weak.sort_values(sort_cols, ascending=[False] * len(sort_cols))
        s_bidir = s_bidir.sort_values(sort_cols, ascending=[False] * len(sort_cols))
        s_one = s_one.sort_values(sort_cols, ascending=[False] * len(sort_cols))
        s_weak = s_weak.sort_values(sort_cols, ascending=[False] * len(sort_cols))

    write_tsv(p_bidir, os.path.join(args.outdir, "01_primary_bidirectional_candidates.tsv"))
    write_tsv(p_one, os.path.join(args.outdir, "02_primary_one_sided_candidate_suppressors.tsv"))
    write_tsv(p_weak, os.path.join(args.outdir, "03_primary_weak_or_mixed_hits.tsv"))

    write_tsv(s_bidir, os.path.join(args.outdir, "04_sensitivity_bidirectional_candidates.tsv"))
    write_tsv(s_one, os.path.join(args.outdir, "05_sensitivity_one_sided_candidate_suppressors.tsv"))
    write_tsv(s_weak, os.path.join(args.outdir, "06_sensitivity_weak_or_mixed_hits.tsv"))

    log("Building mechanism bucket summaries ...")
    primary_bucket = make_bucket_table(pd.concat([p_bidir, p_one], axis=0, ignore_index=True))
    sens_bucket = make_bucket_table(pd.concat([s_bidir, s_one], axis=0, ignore_index=True))

    write_tsv(primary_bucket, os.path.join(args.outdir, "07_primary_mechanism_buckets.tsv"))
    write_tsv(sens_bucket, os.path.join(args.outdir, "08_sensitivity_mechanism_buckets.tsv"))

    summarize_to_txt(
        p_bidir, p_one,
        s_bidir, s_one,
        primary_bucket, sens_bucket,
        os.path.join(args.outdir, "09_postprocessed_summary.txt")
    )

    log("Done.")


if __name__ == "__main__":
    main()
