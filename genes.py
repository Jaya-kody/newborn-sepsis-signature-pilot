# genes.py
# Do genes in the An 2024 sepsis signature change during healthy first-week development?
# Data: GSE111404 (healthy Gambian newborns, blood RNA-seq, day of life 0/1/3/7).
# Needs GSE111404_counts.csv and GSE111404_samples.csv (made by pilot.py).
# Writes tables to results/ and figures to figures/

# ------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------
import os
from datetime import date
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mygene
from scipy.stats import wilcoxon, t

# all outputs go into these two folders
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)


# ------------------------------------------------------------------
# 1. Load the data
# ------------------------------------------------------------------
# counts: one row per gene (Ensembl ID), one column per blood sample (GSM ID)
counts = pd.read_csv("GSE111404_counts.csv", index_col=0)

# samples: one row per blood sample, with columns baby, day, sex
samples = pd.read_csv("GSE111404_samples.csv", index_col=0)

# basic checks before anything else
assert counts.columns.is_unique, "duplicate sample IDs in counts"
assert not counts.isna().any().any(), "missing values in counts"
assert set(counts.columns) == set(samples.index), "counts and samples list different samples"

# put the sample rows in the same order as the count columns
samples = samples.loc[counts.columns]

# Ensembl IDs sometimes end in a version number (ENSG00000123.12); drop it
counts.index = counts.index.str.split(".").str[0]


# ------------------------------------------------------------------
# 2. Turn Ensembl IDs into gene names (e.g. ENSG00000075624 -> ACTB)
# ------------------------------------------------------------------
# The lookup is slow, so do it once and save the result to gene_symbols.csv
if not os.path.exists("gene_symbols.csv"):
    mg = mygene.MyGeneInfo()
    lookup = mg.querymany(
        list(counts.index.unique()),
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=True,
        verbose=False,
    )
    # some IDs return more than one answer; keep the first, but record which ones
    ambiguous = lookup.index[lookup.index.duplicated()].unique()
    pd.Series(ambiguous, name="ensembl_id").to_csv("results/ambiguous_gene_ids.csv", index=False)
    lookup = lookup[~lookup.index.duplicated()]   # keep one answer per ID
    lookup = lookup["symbol"].dropna()            # drop IDs with no gene name
    lookup.to_csv("gene_symbols.csv")
    with open("results/gene_mapping_info.txt", "w") as f:
        f.write(f"mygene {mygene.__version__}, queried {date.today()}\n")
        f.write(f"{len(ambiguous)} Ensembl IDs had more than one symbol; first kept\n")

# sym: a lookup table, Ensembl ID -> gene name
sym = pd.read_csv("gene_symbols.csv", index_col=0)["symbol"]

# keep only genes that have a name, and rename the rows to gene names
named = counts.loc[counts.index.isin(sym.index)].copy()
named.index = sym.loc[named.index].values

# a few names appear twice; add their counts together
n_before = len(named)
named = named.groupby(level=0).sum()
print(f"{len(counts)} Ensembl IDs -> {n_before} with a symbol -> {len(named)} unique symbols")


# ------------------------------------------------------------------
# 3. How much of each sample is haemoglobin (globin) RNA?
# ------------------------------------------------------------------
globin_genes = ["HBB", "HBA1", "HBA2", "HBG1", "HBG2", "HBD", "HBM", "HBQ1", "HBZ", "HBE1"]
globin_present = named.index.intersection(globin_genes)

globin_reads = named.loc[globin_present].sum()   # globin reads per sample
total_reads = counts.sum()                       # ALL gene counts per sample, mapped to a symbol or not
globin_pct = globin_reads / total_reads * 100

print("% of all gene counts from globin genes, by day of life:")
print(globin_pct.groupby(samples["day"]).describe()[["count", "mean", "min", "max"]].round(1))
print()


# ------------------------------------------------------------------
# 4. Normalise (DESeq2-style "median of ratios")
# ------------------------------------------------------------------
def median_of_ratios(table):
    # use only genes detected in every sample to compute the scaling factors
    detected_everywhere = table[(table > 0).all(axis=1)]
    log_counts = np.log(detected_everywhere)
    gene_average = log_counts.mean(axis=1)
    # for each sample: the median of (its value / that gene's average)
    return np.exp(log_counts.sub(gene_average, axis=0).median())

# main analysis: drop globin genes, then normalise
no_globin = named.drop(globin_present)
size_factor = median_of_ratios(no_globin)
expr = np.log2(no_globin.div(size_factor) + 1)
expr.to_csv("results/expr_norm_no_globin.csv")   # used by pca.py

# for comparison only: the two ingredients separately
cpm_log = np.log2(named / named.sum() * 1e6 + 1)                       # CPM, globin kept
mor_with_globin = np.log2(named.div(median_of_ratios(named)) + 1)      # median of ratios, globin kept


# ------------------------------------------------------------------
# 5. Pair each baby's day-0 sample with its later sample
# ------------------------------------------------------------------
pairs = []
for baby, baby_rows in samples.groupby("baby"):
    day0 = baby_rows[baby_rows["day"] == 0]
    later = baby_rows[baby_rows["day"] != 0]
    if len(day0) != 1 or len(later) != 1:
        print(f"EXCLUDED {baby}: {len(day0)} day-0 and {len(later)} later samples")
        continue
    pairs.append({"baby": baby, "day0": day0.index[0], "later": later.index[0],
                  "later_day": later["day"].iloc[0]})
pairs = pd.DataFrame(pairs).set_index("baby")
assert len(pairs) == 25, f"expected 25 pairs, found {len(pairs)}"
print(len(pairs), "babies with a valid pair;", pairs["later_day"].value_counts().sort_index().to_dict())
print()

def paired_change(table, genes):
    # change = later value - day-0 value, one row per baby
    later = table.loc[genes, pairs["later"]].T.set_axis(pairs.index)
    day0 = table.loc[genes, pairs["day0"]].T.set_axis(pairs.index)
    return later - day0


# ------------------------------------------------------------------
# 6. Does removing globin or changing the normalisation matter?
# ------------------------------------------------------------------
# Separate the two effects using the control genes and the signature genes.
check_genes = ["ACTB", "GAPDH", "HSPH1", "BORA", "NCAPG2", "PRIM1"]
compare = pd.DataFrame({
    "CPM, globin kept": paired_change(cpm_log, check_genes).mean(),
    "median-of-ratios, globin kept": paired_change(mor_with_globin, check_genes).mean(),
    "median-of-ratios, globin removed": paired_change(expr, check_genes).mean(),
}).round(2)
print("Mean paired change (log2) under three normalisations:")
print(compare)
compare.to_csv("results/normalisation_comparison.csv")
print()


# ------------------------------------------------------------------
# 7. Statistics
# ------------------------------------------------------------------
# Note: the change is in log2(normalised count + 1). +1 means (count + 1) doubled,
# which is close to "count doubled" only when counts are well above 1 (see baseline table).
signature = ["HSPH1", "BORA", "NCAPG2", "PRIM1"]
others = ["ISG15", "IFI44L", "MMP8", "ACTB", "GAPDH"]
genes = signature + others
changes = paired_change(expr, genes)
changes["later_day"] = pairs["later_day"]
changes.to_csv("results/per_baby_changes.csv")

def bh_adjust(p):
    # Benjamini-Hochberg adjusted p-values
    p = pd.Series(p)
    ranked = p.sort_values()
    m = len(p)
    adj = (ranked * m / np.arange(1, m + 1))[::-1].cummin()[::-1].clip(upper=1)
    return adj.reindex(p.index)

def mean_ci(x):
    # mean and 95% confidence interval (t distribution)
    half = t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return x.mean(), x.mean() - half, x.mean() + half

# 7a. pooled test: is the later-vs-day-0 change centred away from zero?
rows = []
for g in genes:
    m, lo, hi = mean_ci(changes[g])
    rows.append({"gene": g, "family": "signature" if g in signature else "check",
                 "mean_change": m, "ci_low": lo, "ci_high": hi,
                 "n_up": int((changes[g] > 0).sum()), "n": len(changes),
                 "p": wilcoxon(changes[g]).pvalue})
pooled = pd.DataFrame(rows).set_index("gene")
# multiple-testing correction within the primary family (the 4 signature genes)
pooled["p_BH_signature"] = bh_adjust(pooled.loc[signature, "p"])
pooled.to_csv("results/pooled_tests.csv")
print("Pooled change from day 0 to the later sample (all 25 babies):")
print(pooled.round(4).to_string())
print()

# 7b. the same change split by later day, with 95% CIs
rows = []
for g in genes:
    for d, grp in changes.groupby("later_day"):
        m, lo, hi = mean_ci(grp[g])
        rows.append({"gene": g, "later_day": d, "n": len(grp), "mean_change": m,
                     "ci_low": lo, "ci_high": hi, "n_up": int((grp[g] > 0).sum())})
by_day = pd.DataFrame(rows)
by_day.to_csv("results/change_by_later_day.csv", index=False)
print("Change split by later day (signature genes):")
print(by_day[by_day["gene"].isin(signature)].round(2).to_string(index=False))
print()

# 7c. how highly expressed are the signature genes at day 0? (low counts = noisier log changes)
day0_cols = pairs["day0"]
baseline = pd.DataFrame({
    "raw_count_median": named.loc[signature, day0_cols].median(axis=1),
    "normalised_count_median": no_globin.div(size_factor).loc[signature, day0_cols].median(axis=1),
}).round(1)
baseline.to_csv("results/baseline_levels.csv")
print("Day-0 expression of the signature genes:")
print(baseline)
print()


# ------------------------------------------------------------------
# 8. Sensitivity analysis for PRIM1 (exploratory)
# ------------------------------------------------------------------
# Leave each baby out in turn and re-test, so no single baby drives the result.
loo = []
for baby in changes.index:
    rest = changes.drop(index=baby)["PRIM1"]
    loo.append({"left_out": baby, "mean_change": rest.mean(), "p": wilcoxon(rest).pvalue})
loo = pd.DataFrame(loo)
loo.to_csv("results/prim1_leave_one_out.csv", index=False)
print(f"PRIM1 leave-one-baby-out: mean change {loo['mean_change'].min():+.2f} to "
      f"{loo['mean_change'].max():+.2f}, largest p = {loo['p'].max():.2g}")
print()


# ------------------------------------------------------------------
# 9. Figures
# ------------------------------------------------------------------
# 9a. globin share by day
days = [0, 1, 3, 7]
globin_by_day = [globin_pct[samples["day"] == d].values for d in days]
plt.figure(figsize=(4, 3.5))
plt.boxplot(globin_by_day, tick_labels=[f"day {d}" for d in days])
for i, values in enumerate(globin_by_day, start=1):
    plt.scatter(np.full(len(values), i), values, s=12, alpha=0.6)
plt.ylabel("% of gene counts from globin genes")
plt.title("Globin share by day of life")
plt.tight_layout()
plt.savefig("figures/globin_by_day.png", dpi=150)
plt.close()

# 9b. every baby's paired change, split by later day, with mean and 95% CI
def plot_changes(gene_list, filename):
    fig, axes = plt.subplots(1, len(gene_list), figsize=(2.6 * len(gene_list), 3), sharey=True)
    for ax, g in zip(np.atleast_1d(axes), gene_list):
        for i, d in enumerate([1, 3, 7]):
            x = changes.loc[changes["later_day"] == d, g]
            ax.scatter(np.full(len(x), i) + np.random.uniform(-0.1, 0.1, len(x)), x, s=12, alpha=0.6)
            m, lo, hi = mean_ci(x)
            ax.errorbar(i + 0.25, m, yerr=[[m - lo], [hi - m]], fmt="o", color="black", capsize=3)
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xticks([0, 1, 2], ["d1", "d3", "d7"])
        ax.set_title(g)
    np.atleast_1d(axes)[0].set_ylabel("change from day 0 (log2)")
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()

np.random.seed(0)
plot_changes(signature, "figures/signature_change_by_day.png")
plot_changes(["ISG15", "MMP8", "ACTB", "GAPDH"], "figures/check_genes_change_by_day.png")

# 9c. each baby's two samples joined by a line (expression level, not change)
def plot_genes(gene_list, filename):
    fig, axes = plt.subplots(1, len(gene_list), figsize=(3 * len(gene_list), 3), squeeze=False)
    for ax, g in zip(axes[0], gene_list):
        for baby, row in pairs.iterrows():
            ax.plot([0, row["later_day"]], expr.loc[g, [row["day0"], row["later"]]], color="grey", alpha=0.4)
        for d in days:
            day_rows = samples[samples["day"] == d]
            ax.scatter(day_rows["day"], expr.loc[g, day_rows.index], s=15, label=f"day {d}")
        ax.set_title(g)
        ax.set_xlabel("day of life")
        ax.set_xticks(days)
    axes[0][0].set_ylabel("log2(normalised count + 1)")
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()

plot_genes(signature, "figures/an2024_genes_norm.png")
plot_genes(["ACTB", "GAPDH"], "figures/control_genes_norm.png")
plot_genes(["ISG15", "IFI44L", "MMP8"], "figures/known_genes_norm.png")
print("figures saved in figures/, tables saved in results/")
