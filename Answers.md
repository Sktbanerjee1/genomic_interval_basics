### 1. Which two interval sets are the most similar? (10 pts)

The highest Jaccard index (J = 0.607) was observed between the CTCF peak sets from HeLa-S3 and HepG2. In other words, approximately 60.7 % of the combined CTCF-binding footprints are shared between these two cell types, indicating a high degree of conservation in CTCF occupancy across contexts.

### 2. Which two interval sets are the most different? (10 pts)

The lowest Jaccard index (J = 0.013) occurs between HeLa-S3 Jun peaks and HepG2 CTCF peaks. This near-zero overlap reflects almost completely distinct genomic binding landscapes when comparing Jun binding in one cell type to CTCF binding in another.

### 3. Based on these results, which factor, CTCF or Jun, would you predict varies more across cell types? (10 pts)

Jun shows substantially more variability across cell types (HeLa-S3 vs. HepG2 Jun: J = 0.166) than does CTCF (HeLa-S3 vs. HepG2 CTCF: J = 0.607). The much lower Jaccard coefficient for Jun implies that its binding repertoire is more cell-type–specific, whereas CTCF sites remain largely consistent.


### 4. Based on these results, do the genomic locations found by ChIP-seq experiments depend more on the cell-type, or on the transcription factor being assayed? (20 pts)

Comparisons of “same-factor” (CTCF vs. CTCF, Jun vs. Jun) versus “same-cell-type” (CTCF vs. Jun in HeLa-S3 or HepG2) overlaps reveal that transcription factor identity exerts the primary influence on binding-site selection. Same-factor comparisons yield markedly higher Jaccard indices than same-cell-type comparisons, demonstrating that who you profile (CTCF vs. Jun) matters more than where you profile it (HeLa-S3 vs. HepG2).