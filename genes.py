# genes.py
# Do genes in the An 2024 sepsis signature change during healthy first-week development?
# Data: GSE111404 (healthy Gambian newborns, blood RNA-seq, day of life 0/1/3/7).
# Needs GSE111404_counts.csv and GSE111404_samples.csv in the same folder (made by pilot.py).

# ------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mygene
from scipy.stats import wilcoxon


# ------------------------------------------------------------------
# 1. Load the data
# ------------------------------------------------------------------
# counts: one row per gene (Ensembl ID), one column per blood sample (GSM ID)
counts = pd.read_csv("GSE111404_counts.csv", index_col=0)

# samples: one row per blood sample, with columns baby, day, sex
samples = pd.read_csv("GSE111404_samples.csv", index_col=0)

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
    lookup = lookup[~lookup.index.duplicated()]   # keep one answer per ID
    lookup = lookup["symbol"].dropna()            # drop IDs with no gene name
    lookup.to_csv("gene_symbols.csv")

# sym: a lookup table, Ensembl ID -> gene name
sym = pd.read_csv("gene_symbols.csv", index_col=0)["symbol"]

# keep only genes that have a name, and rename the rows to gene names
named = counts.loc[counts.index.isin(sym.index)].copy()
named.index = sym.loc[named.index].values

# a few names appear twice; add their counts together
named = named.groupby(level=0).sum()


# ------------------------------------------------------------------
# 3. How much of each sample is haemoglobin (globin) RNA?
# ------------------------------------------------------------------
globin_genes = ["HBB", "HBA1", "HBA2", "HBG1", "HBG2", "HBD", "HBM", "HBQ1", "HBZ", "HBE1"]
globin_present = named.index.intersection(globin_genes)

globin_reads = named.loc[globin_present].sum()   # globin reads per sample
total_reads = named.sum()                        # all reads per sample
globin_pct = globin_reads / total_reads * 100

print("% of reads from globin genes, by day of life:")
print(globin_pct.groupby(samples["day"]).describe()[["mean", "min", "max"]].round(1))
print()


# ------------------------------------------------------------------
# 4. Normalise without globin (DESeq2-style "median of ratios")
# ------------------------------------------------------------------
# Why: globin falls from ~68% to ~26% of reads over the week. Simple CPM then makes
# every other gene look like it rises. Removing globin and using median-of-ratios
# gives each sample a fair scaling factor.

# drop the globin genes
no_globin = named.drop(globin_present)

# use only genes detected in every sample to compute the scaling factors
detected_everywhere = no_globin[(no_globin > 0).all(axis=1)]
log_counts = np.log(detected_everywhere)

# each gene's average (in log scale) across all samples
gene_average = log_counts.mean(axis=1)

# for each sample: the median of (its value / that gene's average)
size_factor = np.exp(log_counts.sub(gene_average, axis=0).median())

# divide by the size factor, then log2 (the +1 avoids log of zero)
expr = np.log2(no_globin.div(size_factor) + 1)


# ------------------------------------------------------------------
# 5. A function that plots genes against day of life
# ------------------------------------------------------------------
def plot_genes(genes, filename):
    # check which genes exist in the data
    found = [g for g in genes if g in expr.index]
    missing = [g for g in genes if g not in expr.index]
    print(filename, "found:", found, "missing:", missing)
    if not found:
        return

    # one small panel per gene
    fig, axes = plt.subplots(1, len(found), figsize=(3 * len(found), 3), squeeze=False)
    axes = axes[0]

    for ax, gene in zip(axes, found):
        # grey line joining each baby's two samples
        for baby, baby_rows in samples.groupby("baby"):
            ax.plot(baby_rows["day"], expr.loc[gene, baby_rows.index], color="grey", alpha=0.4)

        # one dot per sample, coloured by day
        ax.scatter(samples["day"], expr.loc[gene], c=samples["day"], cmap="viridis", s=15)

        ax.set_title(gene)
        ax.set_xlabel("day of life")
        ax.set_xticks([0, 1, 3, 7])

    axes[0].set_ylabel("log2 normalised (no globin)")
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print("saved", filename)


# ------------------------------------------------------------------
# 6. Make the plots
# ------------------------------------------------------------------
plot_genes(["ACTB", "GAPDH"], "control_genes_norm.png")                              # should be roughly flat
plot_genes(["IFI27", "ISG15", "IFI44L", "LTF", "CAMP", "MMP8"], "known_genes_norm.png")  # Lee 2019 check
plot_genes(["HSPH1", "BORA", "NCAPG2", "PRIM1"], "an2024_genes_norm.png")              # the pilot question
print()


# ------------------------------------------------------------------
# 7. Numbers: how much does each gene change within each baby?
# ------------------------------------------------------------------
# Each baby has a day-0 sample and one later sample (day 1, 3 or 7).
# For every baby: change = later value - day-0 value (in log2 units; +1 = doubled).
genes = ["HSPH1", "BORA", "NCAPG2", "PRIM1", "ISG15", "IFI44L", "MMP8", "ACTB", "GAPDH"]

rows = []
for baby, baby_rows in samples.groupby("baby"):
    day0 = baby_rows[baby_rows["day"] == 0]
    later = baby_rows[baby_rows["day"] != 0]
    if len(day0) != 1 or len(later) != 1:
        continue   # skip babies without exactly one of each

    change = expr.loc[genes, later.index[0]] - expr.loc[genes, day0.index[0]]
    change = change.rename(baby)
    change["later_day"] = later["day"].iloc[0]
    rows.append(change)

changes = pd.DataFrame(rows)

# Wilcoxon signed-rank test: is the change consistently above or below zero?
print("Change from day 0 to the later sample, within each baby:")
for g in genes:
    mean_change = changes[g].mean()
    n_up = (changes[g] > 0).sum()
    p = wilcoxon(changes[g]).pvalue
    print(f"{g:7s} mean change {mean_change:+.2f} log2   up in {n_up}/{len(changes)} babies   p={p:.2g}")
print()

print("Mean change for the 4 signature genes, split by the later day:")
print(changes.groupby("later_day")[genes[:4]].mean().round(2))


# ------------------------------------------------------------------
# 8. Figure: globin share of reads by day of life
# ------------------------------------------------------------------
days = [0, 1, 3, 7]
globin_by_day = [globin_pct[samples["day"] == d].values for d in days]

plt.figure(figsize=(4, 3.5))
plt.boxplot(globin_by_day, tick_labels=[f"day {d}" for d in days])
for i, values in enumerate(globin_by_day, start=1):
    plt.scatter(np.full(len(values), i), values, s=12, alpha=0.6)
plt.ylabel("% of reads from globin genes")
plt.title("Globin share falls over the first week")
plt.tight_layout()
plt.savefig("globin_by_day.png", dpi=150)
plt.close()
print("saved globin_by_day.png")
print()


# ------------------------------------------------------------------
# 9. PRIM1 check: do the few babies with very low day-0 PRIM1 drive the rise?
# ------------------------------------------------------------------
day0_samples = samples[samples["day"] == 0]
prim1_day0 = expr.loc["PRIM1", day0_samples.index]
prim1_day0.index = day0_samples["baby"].values

print("Lowest day-0 PRIM1 values:")
print(prim1_day0.sort_values().head(5).round(2))

low_babies = prim1_day0.sort_values().index[:3]
kept = changes.drop(index=low_babies, errors="ignore")

p_all = wilcoxon(changes["PRIM1"]).pvalue
p_kept = wilcoxon(kept["PRIM1"]).pvalue
print(f"PRIM1 all babies:          mean change {changes['PRIM1'].mean():+.2f}  "
      f"up in {(changes['PRIM1'] > 0).sum()}/{len(changes)}  p={p_all:.2g}")
print(f"PRIM1 without lowest 3 {list(low_babies)}: mean change {kept['PRIM1'].mean():+.2f}  "
      f"up in {(kept['PRIM1'] > 0).sum()}/{len(kept)}  p={p_kept:.2g}")

changes.to_csv("per_baby_changes.csv")
print("saved per_baby_changes.csv")
