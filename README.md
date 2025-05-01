# Interval basics

## Task 1: Implement functions to count overlaps and calculate Jaccard index

Implement the 3 functions in [intervals.R](src/intervals.R). Descriptions of the functions can be found there. Additional instructions for each function follow:

- `countOverlapsSimple`: your function may use a sequential search to identify overlaps (no need to implement binary search or lower complexity algorithm). You should write this function using base R only, without relying on GenomicRanges or other Bioconductor packages. The inputs should be base R objects (a list for the query, and a data.frame for the database, see examples in the function documentation). (20 pts)

- `calculateJaccardScore`: for this function, you may make use of the Bioconductor GenomicRanges package, which includes fast implementations of interval overlap and union. In particular, you may be interested in the `intersect` and `union` methods that operate on `GenomicRanges::GRanges` objects. The inputs to this function should both be `GenomicRanges::GRanges` objects. (20 pts)

- `pairwiseJaccard`: Input should be a base R list of GRanges objects, and you should return a square matrix with dimension equal to the number of elements in the input list. (10 pts)

## Task 2: Calculate and interpret similarity results

There are 4 narrowPeak files stored under [/data](/data). These are ChIP-seq datasets, for two cell types (HepG2 and HeLaS3) and two transcription factors (Jun and CTCF). I filtered the data to keep things smaller and faster for this assignment.

Use your functions to calculate the pairwise Jaccard similarity for this set of 4 experiments, then answer the following questions:

1. Which two interval sets are the most similar? (10 pts)

The highest Jaccard index (J = 0.607) was observed between the CTCF peak sets from HeLa-S3 and HepG2. In other words, approximately 60.7 % of the combined CTCF-binding footprints are shared between these two cell types, indicating a high degree of conservation in CTCF occupancy across contexts.

2. Which two interval sets are the most different? (10 pts)

The lowest Jaccard index (J = 0.013) occurs between HeLa-S3 Jun peaks and HepG2 CTCF peaks. This near-zero overlap reflects almost completely distinct genomic binding landscapes when comparing Jun binding in one cell type to CTCF binding in another.

3. Based on these results, which factor, CTCF or Jun, would you predict varies more across cell types? (10 pts)

Jun shows substantially more variability across cell types (HeLa-S3 vs. HepG2 Jun: J = 0.166) than does CTCF (HeLa-S3 vs. HepG2 CTCF: J = 0.607). The much lower Jaccard coefficient for Jun implies that its binding repertoire is more cell-type–specific, whereas CTCF sites remain largely consistent.


4. Based on these results, do the genomic locations found by ChIP-seq experiments depend more on the cell-type, or on the transcription factor being assayed? (20 pts)

Comparisons of “same-factor” (CTCF vs. CTCF, Jun vs. Jun) versus “same-cell-type” (CTCF vs. Jun in HeLa-S3 or HepG2) overlaps reveal that transcription factor identity exerts the primary influence on binding-site selection. Same-factor comparisons yield markedly higher Jaccard indices than same-cell-type comparisons, demonstrating that who you profile (CTCF vs. Jun) matters more than where you profile it (HeLa-S3 vs. HepG2).


## Testing:

You can test your work by running within the interval_basics directory:

```
Rscript tests/testDriver.R
```
