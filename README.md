# Do sepsis-signature genes change with age in healthy newborns?

A small pilot on public data. The An et al. (2024) early-onset sepsis signature uses four genes (HSPH1, BORA, NCAPG2, PRIM1) measured in blood at birth. If any of them also change during normal first-week development, a single cut-off ignores age. This pilot asks whether they move in healthy babies.

## Data
GEO [GSE111404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111404) (Lee et al. 2019, *Nat Commun*): whole-blood RNA-seq from 25 healthy term newborns in The Gambia (EPIC cohort). Every baby was sampled on day 0 and once more on day 1 (n=8), day 3 (n=8) or day 7 (n=9). That makes 50 samples.

## Method
1. `pilot.py` downloads the counts and sample information and builds one counts table (58,051 genes x 50 samples).
2. `pca.py` gives an overview (PCA on the 2,000 most variable genes): PC1 (80%) follows day of life, and PC2 (3%) follows sex.
3. `genes.py` does the following:
   - Measures the share of reads from globin genes. It falls from about 68% at day 0 to 26% at day 7. Plain CPM normalisation therefore makes every other gene look like it rises.
   - Removes the globin genes and normalises with DESeq2-style median-of-ratios size factors. Lee et al. 2019 also removed globin.
   - Computes each baby's change from day 0 to the later sample and tests it with a paired Wilcoxon signed-rank test.

## Results (globin removed, n = 25 babies)

| Gene | Role | Mean change (log2) | Up in | p |
|---|---|---|---|---|
| HSPH1 | signature | -0.04 | 14/25 | 0.87 |
| BORA | signature | +0.11 | 16/25 | 0.21 |
| NCAPG2 | signature | +0.39 | 17/25 | 0.03 |
| PRIM1 | signature | +0.69 | 21/25 | 1e-05 |
| ISG15 | known (Lee 2019) | +1.34 | 23/25 | 4e-06 |
| IFI44L | known (Lee 2019) | +1.58 | 19/25 | 3e-04 |
| MMP8 | known (Lee 2019) | -2.16 | 1/25 | 2e-07 |
| ACTB | control | +0.21 | 20/25 | 1e-04 |
| GAPDH | control | +0.01 | 13/25 | 0.98 |

- **Pipeline check passes:** the interferon genes rise and MMP8 falls, as in Lee 2019. GAPDH stays flat.
- **Most of the signature does not change with age.** HSPH1 and BORA show no change.
- **PRIM1 rises with age, and NCAPG2 possibly does.** For PRIM1, removing the 3 babies with the lowest day-0 values leaves the rise in place (+0.56, p = 8e-05).

![Globin share by day](figures/globin_by_day.png)
![Signature genes after removing globin](figures/an2024_genes_norm.png)

## Limitations
- The sample is small (8-9 babies per later day), and no baby was sampled on all four days.
- All babies are term babies from one population (The Gambia).
- The samples are bulk blood, so a change can reflect a shift in cell mix rather than a change in gene activity. That matters here because these are cell-cycle genes.
- No correction for multiple testing was applied (this is a pilot). The analysis uses Python with simplified methods; the standard approach is DESeq2 in R.

## How to run
```
pip install -r requirements.txt
python pilot.py
python pca.py
python genes.py
```

## References
- Lee AH et al. 2019. Dynamic molecular changes during the first week of human life follow a robust developmental trajectory. *Nat Commun* 10:1092.
- An et al. 2024. Early-onset sepsis host-response signature (HSPH1, BORA, NCAPG2, PRIM1).
