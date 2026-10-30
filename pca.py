# pca.py
# A first overview of the data: do the 50 blood samples group by day of life? By sex?
# PCA squeezes thousands of genes into 2 numbers per sample (PC1, PC2) so we can plot them.
# Run genes.py first: it writes results/expr_norm_no_globin.csv (globin removed,
# median-of-ratios normalised, log2), so the PCA uses the same normalisation as the gene tests.

# ------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.decomposition import PCA

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)


# ------------------------------------------------------------------
# 1. Load the data
# ------------------------------------------------------------------
# expr: one row per gene, one column per sample, log2(normalised count + 1)
expr = pd.read_csv("results/expr_norm_no_globin.csv", index_col=0)

# samples: one row per blood sample, with columns baby, day, sex
samples = pd.read_csv("GSE111404_samples.csv", index_col=0)
samples = samples.loc[expr.columns]


# ------------------------------------------------------------------
# 2. Remove genes that are barely detected
# ------------------------------------------------------------------
# Keep a gene only if its normalised count is above 10 in at least 10 of the 50 samples.
# Very rare genes are mostly noise.
keep = (expr > np.log2(10 + 1)).sum(axis=1) >= 10
expr = expr[keep]
print(expr.shape, "(genes kept, samples)")


# ------------------------------------------------------------------
# 3. PCA on the 2,000 most variable genes
# ------------------------------------------------------------------
gene_variance = expr.var(axis=1)
top_genes = gene_variance.sort_values().index[-2000:]
data_for_pca = expr.loc[top_genes].T.values      # one row per sample

pca = PCA(n_components=2)
xy = pca.fit_transform(data_for_pca)              # xy[:, 0] = PC1, xy[:, 1] = PC2

pc1_label = f"PC1 ({pca.explained_variance_ratio_[0]:.0%})"
pc2_label = f"PC2 ({pca.explained_variance_ratio_[1]:.0%})"
print(pc1_label, pc2_label)


# ------------------------------------------------------------------
# 4. How strongly is each PC associated with day and with sex?
# ------------------------------------------------------------------
# These are associations only; they do not say what causes a component.
rows = []
for i in range(2):
    rho, p_day = spearmanr(xy[:, i], samples["day"])
    p_sex = mannwhitneyu(xy[samples["sex"].values == "M", i], xy[samples["sex"].values == "F", i]).pvalue
    rows.append({"PC": f"PC{i + 1}", "variance": pca.explained_variance_ratio_[i],
                 "spearman_rho_day": rho, "p_day": p_day, "p_sex": p_sex})
assoc = pd.DataFrame(rows).set_index("PC")
assoc.to_csv("results/pca_associations.csv")
print(assoc.round(4))


# ------------------------------------------------------------------
# 5. Plots: coloured by day of life, then by sex
# ------------------------------------------------------------------
def pca_plot(groups, column, title, filename):
    plt.figure(figsize=(6, 5))
    for value, name, colour in groups:
        rows = samples[column].values == value
        plt.scatter(xy[rows, 0], xy[rows, 1], c=colour, label=name, s=40)
    plt.xlabel(pc1_label)
    plt.ylabel(pc2_label)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()

pca_plot([(0, "day 0", "tab:blue"), (1, "day 1", "tab:green"), (3, "day 3", "tab:orange"), (7, "day 7", "tab:red")],
         "day", "Newborn blood, first week (globin removed)", "figures/pca_by_day.png")
pca_plot([("M", "Male", "tab:blue"), ("F", "Female", "tab:red")],
         "sex", "Newborn blood by sex (globin removed)", "figures/pca_by_sex.png")
print("saved figures/pca_by_day.png and figures/pca_by_sex.png")
