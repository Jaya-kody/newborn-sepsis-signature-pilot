# pca.py
# A first overview of the data: do the 50 blood samples group by day of life? By sex?
# PCA squeezes thousands of genes into 2 numbers per sample (PC1, PC2) so we can plot them.
# Needs GSE111404_counts.csv and GSE111404_samples.csv in the same folder.

# ------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


# ------------------------------------------------------------------
# 1. Load the data
# ------------------------------------------------------------------
# counts: one row per gene, one column per blood sample
counts = pd.read_csv("GSE111404_counts.csv", index_col=0)

# samples: one row per blood sample, with columns baby, day, sex
samples = pd.read_csv("GSE111404_samples.csv", index_col=0)

# put the sample rows in the same order as the count columns
samples = samples.loc[counts.columns]


# ------------------------------------------------------------------
# 2. Normalise: counts per million (CPM)
# ------------------------------------------------------------------
# Samples were sequenced to different depths (2.6 to 14.6 million reads).
# Dividing by each sample's total and multiplying by 1 million makes them comparable.
# Note: CPM still includes globin genes. We found later (genes.py) that this matters
# for single genes; for this overview plot it is fine as a first look.
reads_per_sample = counts.sum()
cpm = counts / reads_per_sample * 1e6


# ------------------------------------------------------------------
# 3. Remove genes that are barely detected
# ------------------------------------------------------------------
# Keep a gene only if it has CPM above 1 in at least 10 of the 50 samples.
# Very rare genes are mostly noise.
samples_with_gene = (cpm > 1).sum(axis=1)
keep = samples_with_gene >= 10
cpm = cpm[keep]

# log2 scale: +1 doubles the value. The +1 avoids log of zero.
logcpm = np.log2(cpm + 1)
print(logcpm.shape, "(genes kept, samples)")


# ------------------------------------------------------------------
# 4. PCA on the 2,000 most variable genes
# ------------------------------------------------------------------
# Genes that barely change between samples carry no information about grouping,
# so use only the 2,000 genes that vary most.
gene_variance = logcpm.var(axis=1)
top_genes = gene_variance.sort_values().index[-2000:]
top = logcpm.loc[top_genes]

# PCA wants one row per sample, so flip the table (.T = transpose)
data_for_pca = top.T.values

pca = PCA(n_components=2)
pca.fit(data_for_pca)
xy = pca.transform(data_for_pca)   # xy[:, 0] = PC1, xy[:, 1] = PC2 for each sample

# share of the total variation each PC captures, used in the axis labels
pc1_label = f"PC1 ({pca.explained_variance_ratio_[0]:.0%})"
pc2_label = f"PC2 ({pca.explained_variance_ratio_[1]:.0%})"


# ------------------------------------------------------------------
# 5. Plot 1: samples coloured by day of life
# ------------------------------------------------------------------
plt.figure(figsize=(6, 5))

day_colours = {0: "tab:blue", 1: "tab:green", 3: "tab:orange", 7: "tab:red"}
for day, colour in day_colours.items():
    is_this_day = samples["day"].values == day
    plt.scatter(xy[is_this_day, 0], xy[is_this_day, 1], c=colour, label=f"day {day}", s=40)

plt.xlabel(pc1_label)
plt.ylabel(pc2_label)
plt.title("Newborn blood, first week (GSE111404)")
plt.legend()
plt.tight_layout()
plt.savefig("pca_by_day.png", dpi=150)
print("saved pca_by_day.png")


# ------------------------------------------------------------------
# 6. Plot 2: the same points, coloured by sex
# ------------------------------------------------------------------
plt.figure(figsize=(6, 5))

sex_groups = [("M", "Male", "tab:blue"), ("F", "Female", "tab:red")]
for code, name, colour in sex_groups:
    is_this_sex = samples["sex"].values == code
    plt.scatter(xy[is_this_sex, 0], xy[is_this_sex, 1], c=colour, label=name, s=40)

plt.xlabel(pc1_label)
plt.ylabel(pc2_label)
plt.title("Newborn blood by sex (GSE111404)")
plt.legend()
plt.tight_layout()
plt.savefig("pca_by_sex.png", dpi=150)
print("saved pca_by_sex.png")


# ------------------------------------------------------------------
# 7. Show both plots on screen (close the windows to end the script)
# ------------------------------------------------------------------
plt.show()
