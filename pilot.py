# pilot.py
# Step 1 of the pilot: download GSE111404 from GEO and turn it into two tidy tables.
#   GSE111404_samples.csv : one row per blood sample (baby, day of life, sex)
#   GSE111404_counts.csv  : one row per gene, one column per sample (read counts)
# GSE111404 = whole-blood RNA-seq from healthy newborns in The Gambia (Lee et al. 2019).

# ------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------
import os
import gzip
import tarfile
import urllib.request

import pandas as pd
import GEOparse


# ------------------------------------------------------------------
# 1. Download the sample information from GEO
# ------------------------------------------------------------------
# GEOparse downloads a description of the study and every sample.
# If the file is already in this folder, it reuses it instead of downloading again.
gse = GEOparse.get_GEO("GSE111404", destdir=".", silent=True)

# pheno: one row per sample, one column per piece of information
pheno = gse.phenotype_data.copy()


# ------------------------------------------------------------------
# 2. Look at what information we have
# ------------------------------------------------------------------
print(pheno.shape, "(samples, columns)")
print(pheno.columns.tolist())

# Columns starting with "characteristics" usually hold day of life, sex and baby ID.
# Print the first few different values in each one.
for col in pheno.columns:
    if col.startswith("characteristics"):
        print(col, pheno[col].unique()[:10])
print()


# ------------------------------------------------------------------
# 3. Pull out baby ID and day of life from the sample title
# ------------------------------------------------------------------
# Titles look like "RNA_H01_D0": baby H01, day of life 0.
# str.extract keeps the part of the text inside the brackets ( ).
pheno["baby"] = pheno["title"].str.extract(r"_(H\d+)_")        # "H" followed by digits
pheno["day"] = pheno["title"].str.extract(r"_D(\d+)$")         # digits after "_D" at the end
pheno["day"] = pheno["day"].astype(int)                        # text "7" -> number 7

# sex is already in its own column; give it a short name
pheno["sex"] = pheno["characteristics_ch1.2.Sex"]

print(pheno[["title", "baby", "day", "sex"]].head(10))
print(pheno["day"].value_counts(), "\n")                       # samples per day
print(pheno["baby"].nunique(), "babies\n")


# ------------------------------------------------------------------
# 4. Download the gene counts
# ------------------------------------------------------------------
# GEO stores the counts as one archive (.tar) holding one small file per sample.
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111404/suppl/GSE111404_RAW.tar"

if not os.path.exists("GSE111404_RAW.tar"):
    urllib.request.urlretrieve(url, "GSE111404_RAW.tar")

# unpack the archive into a folder called "counts" (only the first time)
if not os.path.exists("counts"):
    with tarfile.open("GSE111404_RAW.tar") as tar:
        tar.extractall("counts")

files = sorted(f for f in os.listdir("counts") if f.endswith(".count.txt.gz"))
print(len(files), "count files, e.g.", files[:3])


# ------------------------------------------------------------------
# 5. Peek inside one file to see its format
# ------------------------------------------------------------------
# Each line is: Ensembl gene ID <tab> number of reads for that gene
first_file = os.path.join("counts", files[0])
with gzip.open(first_file, "rt") as f:
    for _ in range(5):
        print(f.readline().rstrip())
print()


# ------------------------------------------------------------------
# 6. Combine the 50 files into one table
# ------------------------------------------------------------------
columns = {}
for f in files:
    sample_id = f.split("_")[0]                 # file name starts with the GSM sample ID
    path = os.path.join("counts", f)
    one_sample = pd.read_csv(path, sep="\t", header=None, index_col=0, compression="gzip")
    columns[sample_id] = one_sample.iloc[:, 0]  # keep the counts column as a Series

counts = pd.DataFrame(columns)                  # genes x samples

# htseq adds summary rows at the bottom (e.g. "__no_feature"); they are not genes
is_summary_row = counts.index.str.startswith("__")
counts = counts[~is_summary_row]


# ------------------------------------------------------------------
# 7. Quick checks
# ------------------------------------------------------------------
print(counts.shape, "(genes, samples)")
print(counts.iloc[:5, :5])

reads_per_sample = counts.sum()
print(reads_per_sample.describe())              # sequencing depth per sample

genes_detected = (counts.sum(axis=1) > 0).sum()
print(genes_detected, "genes detected in at least one sample")


# ------------------------------------------------------------------
# 8. Save the two tidy tables
# ------------------------------------------------------------------
counts.to_csv("GSE111404_counts.csv")
pheno[["baby", "day", "sex"]].to_csv("GSE111404_samples.csv")
print("saved GSE111404_counts.csv and GSE111404_samples.csv")
