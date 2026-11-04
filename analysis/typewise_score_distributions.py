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
sys.path.insert(0, "/home/choukaihsuan/circRNA_agent/benchmark/scripts")
import rank_biomarkers as rb  # noqa: E402  (reuse _minmax, not reimplemented)
import accuracy_benchmark as ab  # noqa: E402  (reuse _parse_id, not reimplemented)

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


# GSE133998 and GSE143797 predate the interaction_queried column. Verified
# (item 27, not assumed) by reading their actual biomarker_candidates.tsv:
# n_mirna/n_rbp are present for every row (no NaN) and non-zero for the
# large majority (72/79, 74/79 and 50/56, 51/56 respectively) -- a pattern
# consistent with "this file predates the top-N-only query restriction and
# every candidate was queried," not "these were never queried." The column
# is absent, not the data; defaulting to interaction_queried=False here
# would silently zero out the one component (mirna) that turned out to
# carry the analysis's only consistent cross-dataset signal.
_FULLY_QUERIED_LEGACY_DATASETS = {"GSE133998", "GSE143797"}


def load_bm(gse):
    d = f"{BASE}/{gse}_results"
    bm = pd.read_csv(f"{d}/de/biomarker_candidates.tsv", sep="\t")
    if "interaction_queried" not in bm.columns:
        bm = bm.copy()
        if gse in _FULLY_QUERIED_LEGACY_DATASETS:
            bm["interaction_queried"] = True
        else:
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

    # Pooling can be driven by cross-dataset composition (size, depth) rather
    # than the within-dataset Type_I vs Type_II difference -- same risk the
    # per-dataset + signed-rank treatment exists for biomarker_score/
    # components. confidence_score IS pool-comparable (not per-dataset
    # normalized), but that doesn't make the pooled test immune to a
    # composition artifact; run the same paired check as (c) on it too and
    # compare the two answers directly, since they can legitimately disagree.
    conf_result = per_dataset_two_group(merged_cache, "Type", "confidence_score", "Type_I", "Type_II")
    conf_signed_rank = signed_rank_summary(conf_result)
    print("\nPer-dataset confidence_score Type_I vs Type_II (paired check, same method as (c)):")
    for gse in DATASETS:
        print(f"  {gse:15s} {conf_result[gse]}")
    print(f"  {conf_signed_rank}")
    print(f"  Pooled MWU p={u_p:.6g} vs paired Wilcoxon signed-rank p="
          f"{conf_signed_rank.get('wilcoxon_p', 'n/a')} -- compare directly, do not average.")

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

    # ── (f) selection-bias check for the mirna/rbp comparison ─────────────
    # interaction_queried is set by predict_interactions.py's top-N-by-p-value
    # union (see that script's --de-edger/--de-deseq2/--de-limma union mode).
    # Type_II's p-values are structurally smaller (it requires significance
    # on BOTH the BSJ ratio test and the independent FSJ test), so a larger
    # share of Type_II candidates can end up inside that top-N net than of
    # Type_I -- making the mirna/rbp comparison implicitly "typical Type_II"
    # vs "only the most significant slice of Type_I," not "Type_II vs
    # Type_I" in general. Checked directly: per-dataset queried-rate by
    # Type, and the p-value distribution *within* the queried subset by
    # Type (if queried Type_I's p-values are already as small as queried
    # Type_II's, the two groups being compared are not differentially
    # filtered by significance even though the raw pools are).
    print("\n" + "=" * 70)
    print("(f) Selection-bias check: queried-rate and within-queried p-value by Type")
    print("=" * 70)
    for gse in DATASETS:
        bm = bm_cache[gse]
        row = {}
        for t in ["Type_I", "Type_II"]:
            sub = bm.loc[bm["Type"] == t]
            n_total = len(sub)
            n_queried = int(sub["interaction_queried"].sum())
            rate = round(n_queried / n_total, 4) if n_total else None
            p_queried = sub.loc[sub["interaction_queried"], SIG_COL]
            row[t] = {"n_total": n_total, "n_queried": n_queried, "queried_rate": rate,
                      "p_within_queried": describe(p_queried)}
        print(f"  {gse:15s} I={row['Type_I']}")
        print(f"  {'':15s} II={row['Type_II']}")

    # ── (g) length confound check for the mirna signal ─────────────────────
    # n_mirna is a count of predicted binding sites (CircInteractome/ENCORI);
    # predicted site count scales with sequence length -- a longer sequence
    # has more positions for a seed match. mirna_norm only min-maxes within
    # each dataset's own significant set; it is never length-normalized. If
    # Type_II circRNAs are simply longer, and n_mirna correlates with length
    # within a dataset, the mirna direction found in (d) could be a length
    # artifact rather than a regulatory signal -- checked directly, not
    # assumed either way.
    print("\n" + "=" * 70)
    print("(g) Length confound check: Type_I vs Type_II circRNA length, and length-vs-n_mirna correlation")
    print("=" * 70)
    MIRNA_USABLE = ["GSE113230", "SRP156355", "GSE77509", "GSE130078", "GSE248612",
                    "GSE221107", "PRJNA553289", "GSE121842", "GSE136569", "GSE143797"]
    length_diffs = []
    for gse in MIRNA_USABLE:
        bm = bm_cache[gse]
        df = bm.copy()
        parsed = df["circ_id"].apply(ab._parse_id)
        df["length"] = parsed.apply(lambda t: (t[2] - t[1]) if t is not None else np.nan)
        g1 = df.loc[df["Type"] == "Type_I", "length"].dropna()
        g2 = df.loc[df["Type"] == "Type_II", "length"].dropna()
        n1, n2 = len(g1), len(g2)
        entry = {"n1": n1, "n2": n2}
        if n1 >= MIN_N_FOR_TEST and n2 >= MIN_N_FOR_TEST:
            entry["median_I"] = round(float(g1.median()), 1)
            entry["median_II"] = round(float(g2.median()), 1)
            entry["diff_II_minus_I"] = round(entry["median_II"] - entry["median_I"], 1)
            try:
                _, p = stats.mannwhitneyu(g1, g2, alternative="two-sided")
                entry["mwu_p"] = round(float(p), 6)
            except ValueError:
                entry["mwu_p"] = None
            length_diffs.append(entry["diff_II_minus_I"])

        # Spearman(n_mirna, length) within this dataset's queried subset
        # (both types combined -- asking whether length predicts n_mirna at
        # all, not a Type-specific question).
        qd = df.loc[df["interaction_queried"], ["n_mirna", "length"]].dropna()
        rho_entry = {"n_queried": len(qd)}
        if len(qd) >= MIN_N_FOR_TEST:
            rho, rho_p = stats.spearmanr(qd["n_mirna"], qd["length"])
            rho_entry["spearman_rho"] = round(float(rho), 4)
            rho_entry["spearman_p"] = round(float(rho_p), 6)

        print(f"  {gse:15s} length: {entry}")
        print(f"  {'':15s} n_mirna~length Spearman: {rho_entry}")

    print(f"\n  length usable datasets: {len(length_diffs)}/{len(MIRNA_USABLE)}, diffs={length_diffs}")
    if len(length_diffs) >= 2 and any(d != 0 for d in length_diffs):
        try:
            w, p = stats.wilcoxon(length_diffs)
            print(f"  Wilcoxon signed-rank on length diffs: W={w:.2f}, p={p:.6g}")
        except ValueError as e:
            print(f"  Wilcoxon failed: {e}")
    n_pos = sum(1 for d in length_diffs if d > 0)
    n_neg = sum(1 for d in length_diffs if d < 0)
    print(f"  direction: {n_pos} datasets II>I (longer), {n_neg} datasets II<I (shorter)")


if __name__ == "__main__":
    main()
