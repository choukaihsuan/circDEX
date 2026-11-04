#!/usr/bin/env python3
"""Item 26 (committee question 12): Type-wise score distribution analysis
across the 15-dataset cross-dataset cohort (not the platform's full 20).

Zero new pipeline computation. Reads existing de_results.tsv,
circbase_annotated.tsv, and biomarker_candidates.tsv per dataset. Reuses
scripts/rank_biomarkers.py's own _minmax() directly rather than
reimplementing the normalization by hand.

Run on the server (needs pandas/scipy in the ciriquant env):
  conda run -n ciriquant python3 analysis/typewise_score_distributions.py

Two datasets (GSE133998, GSE143797) have a biomarker_candidates.tsv
predating the `interaction_queried` column -- handled explicitly (treated
as not-queried, consistent with how the main pipeline treats that state),
not silently defaulted without comment.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/home/choukaihsuan/circRNA_agent/scripts")
import rank_biomarkers as rb  # noqa: E402  (reuse _minmax, not reimplemented)

BASE = "/home3/choukaihsuan"
DATASETS = ["GSE113230", "GSE133998", "SRP156355", "GSE77509", "GSE130078",
            "GSE248612", "GSE221107", "PRJNA553289", "GSE148036", "GSE121842",
            "GSE136569", "GSE143797", "GSE108735", "GSE171011", "GSE97239"]

SIG_CAP = 10.0
FC_CAP = 5.0
SIG_COL = "pvalue"
MIN_N_FOR_TEST = 5


def load_merged(gse):
    d = f"{BASE}/{gse}_results"
    de = pd.read_csv(f"{d}/de/de_results.tsv", sep="\t")
    annot = pd.read_csv(f"{d}/circRNA/circbase_annotated.tsv", sep="\t")
    annot_cols = [c for c in ["circ_id", "confidence_score", "in_circbase",
                               "circbase_id", "circbase_gene"] if c in annot.columns]
    if "circ_id" in annot.columns and "circ_id" in de.columns:
        merged = de.merge(annot[annot_cols], on="circ_id", how="left")
    else:
        merged = de.copy()
    merged["confidence_score"] = pd.to_numeric(merged.get("confidence_score", 0), errors="coerce").fillna(0)
    return merged


def load_bm(gse):
    d = f"{BASE}/{gse}_results"
    bm = pd.read_csv(f"{d}/de/biomarker_candidates.tsv", sep="\t")
    if "interaction_queried" not in bm.columns:
        bm = bm.copy()
        bm["interaction_queried"] = False
    return bm


def compute_components(sig):
    """Mirrors rank_biomarkers._compute_weighted_score's raw/norm construction,
    exposing the intermediate normalized components it only ever returns
    combined (as biomarker_score)."""
    raw = pd.DataFrame(index=sig.index)
    raw["sig"] = sig[SIG_COL].apply(lambda p: min(-np.log10(max(p, 1e-10)), SIG_CAP))
    raw["fc"] = sig["log2FC"].abs().clip(upper=FC_CAP)
    raw["conf"] = sig["confidence_score"].clip(lower=0)
    raw["known"] = sig["in_circbase"].astype(float)
    raw["mirna"] = sig["n_mirna"].where(sig["interaction_queried"])
    raw["rbp"] = sig["n_rbp"].where(sig["interaction_queried"])
    norm = pd.DataFrame(index=sig.index, columns=raw.columns, dtype=float)
    for c in raw.columns:
        if c == "known":
            norm[c] = raw[c]
            continue
        valid = raw[c].notna()
        if valid.any():
            norm.loc[valid, c] = rb._minmax(raw.loc[valid, c])
    return raw, norm


def describe(x):
    x = pd.Series(x).dropna()
    n = len(x)
    if n == 0 or n < MIN_N_FOR_TEST:
        return {"n": n}
    return {"n": n, "median": round(float(x.median()), 4),
            "q1": round(float(x.quantile(0.25)), 4), "q3": round(float(x.quantile(0.75)), 4)}


def per_dataset_two_group(values_by_dataset, group_col, value_col, g1_label, g2_label):
    """For each dataset, (n1, n2, median_I, median_II, diff, MWU p) if both
    groups have n >= MIN_N_FOR_TEST, else just n1/n2. Extracts value_col as
    a Series explicitly -- passing a multi-column frame through to
    .median() silently returns a per-column Series, not a scalar, which is
    the wrong thing and was caught only by actually running this, not by
    reading it."""
    out = {}
    for gse, df in values_by_dataset.items():
        g1 = df.loc[df[group_col] == g1_label, value_col].dropna()
        g2 = df.loc[df[group_col] == g2_label, value_col].dropna()
        n1, n2 = len(g1), len(g2)
        entry = {"n1": n1, "n2": n2}
        if n1 >= MIN_N_FOR_TEST and n2 >= MIN_N_FOR_TEST:
            entry["median_I"] = round(float(g1.median()), 4)
            entry["median_II"] = round(float(g2.median()), 4)
            entry["diff_II_minus_I"] = round(entry["median_II"] - entry["median_I"], 4)
            try:
                _, p = stats.mannwhitneyu(g1, g2, alternative="two-sided")
                entry["mwu_p"] = round(float(p), 6)
            except ValueError:
                entry["mwu_p"] = None
        out[gse] = entry
    return out


def signed_rank_summary(per_dataset):
    diffs = [e["diff_II_minus_I"] for e in per_dataset.values() if "diff_II_minus_I" in e]
    result = {"n_usable_datasets": len(diffs), "diffs": diffs}
    if len(diffs) >= 2 and any(d != 0 for d in diffs):
        try:
            w, p = stats.wilcoxon(diffs)
            result["wilcoxon_W"] = round(float(w), 2)
            result["wilcoxon_p"] = round(float(p), 6)
        except ValueError as e:
            result["wilcoxon_error"] = str(e)
    result["n_II_gt_I"] = sum(1 for d in diffs if d > 0)
    result["n_II_lt_I"] = sum(1 for d in diffs if d < 0)
    result["n_tied"] = sum(1 for d in diffs if d == 0)
    return result


def main():
    merged_cache = {gse: load_merged(gse) for gse in DATASETS}
    bm_cache = {gse: load_bm(gse) for gse in DATASETS}

    print("=" * 70)
    print("(a) Per-dataset Type_I / Type_II / Type_III / NS counts")
    print("=" * 70)
    for gse in DATASETS:
        vc = merged_cache[gse]["Type"].value_counts(dropna=False).to_dict()
        print(f"{gse:15s} Type_I={vc.get('Type_I',0):6d} Type_II={vc.get('Type_II',0):6d} "
              f"Type_III={vc.get('Type_III',0):6d} NS={vc.get('NS',0):6d}")

    print("\n" + "=" * 70)
    print("(b) confidence_score -- 4 groups per dataset; pooled test on Type_I vs Type_II only")
    print("=" * 70)
    pool = {"Type_I": [], "Type_II": []}
    for gse in DATASETS:
        m = merged_cache[gse]
        row = {}
        for t in ["Type_I", "Type_II", "Type_III", "NS"]:
            vals = m.loc[m["Type"] == t, "confidence_score"]
            row[t] = describe(vals)
            if t in pool:
                pool[t].extend(vals.dropna().tolist())
        print(f"{gse:15s} I={row['Type_I']} II={row['Type_II']} III={row['Type_III']} NS={row['NS']}")
    print("\nPooled (15 datasets), confidence_score:")
    for t in ["Type_I", "Type_II"]:
        print(f"  {t}: {describe(pool[t])}")
    u_stat, u_p = stats.mannwhitneyu(pool["Type_I"], pool["Type_II"], alternative="two-sided")
    print(f"  Mann-Whitney U (pooled Type_I vs Type_II): U={u_stat:.1f}, p={u_p:.6g}")

    print("\n" + "=" * 70)
    print("(c) biomarker_score -- per-dataset Type_I vs Type_II, then signed-rank across datasets")
    print("=" * 70)
    bm_result = per_dataset_two_group(bm_cache, "Type", "biomarker_score", "Type_I", "Type_II")
    for gse in DATASETS:
        print(f"  {gse:15s} {bm_result[gse]}")
    print(f"\n  {signed_rank_summary(bm_result)}")

    print("\n" + "=" * 70)
    print("(d) Six components -- same per-dataset + signed-rank treatment (known handled in (e))")
    print("=" * 70)
    for comp in ["sig", "fc", "conf", "mirna", "rbp"]:
        comp_frames = {}
        for gse in DATASETS:
            bm = bm_cache[gse]
            _, norm = compute_components(bm)
            df = bm[["Type"]].copy()
            df["_val"] = norm[comp].values
            comp_frames[gse] = df
        per_ds = per_dataset_two_group(comp_frames, "Type", "_val", "Type_I", "Type_II")
        print(f"\n  --- component: {comp} ---")
        for gse in DATASETS:
            print(f"    {gse:15s} {per_ds[gse]}")
        print(f"    {signed_rank_summary(per_ds)}")

    print("\n" + "=" * 70)
    print("(e) known (in_circbase, 0/1) -- proportion + Fisher exact per dataset")
    print("=" * 70)
    for gse in DATASETS:
        bm = bm_cache[gse]
        g1 = bm.loc[bm["Type"] == "Type_I", "in_circbase"]
        g2 = bm.loc[bm["Type"] == "Type_II", "in_circbase"]
        n1, n2 = len(g1), len(g2)
        if n1 < MIN_N_FOR_TEST or n2 < MIN_N_FOR_TEST:
            print(f"  {gse:15s} n1={n1} n2={n2} (too small, not tested)")
            continue
        k1, k2 = int(g1.sum()), int(g2.sum())
        odds, p = stats.fisher_exact([[k1, n1 - k1], [k2, n2 - k2]])
        print(f"  {gse:15s} Type_I: {k1}/{n1} ({k1/n1:.3f})  Type_II: {k2}/{n2} ({k2/n2:.3f})  "
              f"Fisher p={p:.6g} odds={odds:.3f}")


if __name__ == "__main__":
    main()
