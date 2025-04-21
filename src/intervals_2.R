# src/intervals.R

#––– Dependencies
if (!requireNamespace("GenomicRanges", quietly=TRUE)) {
  stop("Please install Bioconductor package 'GenomicRanges'")
}
if (!requireNamespace("IRanges", quietly=TRUE)) {
  stop("Please install Bioconductor package 'IRanges'")
}

#––– Random tester (leave as is) –––
randomGRanges <- function(){
  starts = sample(1:1e7, 500)
  sizes  = sample(1:1e4, 500)
  ends   = starts + sizes
  ret = GenomicRanges::GRanges(
    seqnames = "chrX",
    ranges   = IRanges::IRanges(starts, ends)
  )
  GenomicRanges::reduce(ret)
}

#––– Simple overlap counter (unchanged) –––
countOverlapsSimple <- function(query, database) {
  q_start <- query$start
  q_end   <- query$end
  overlaps <- !(database$ends < q_start | database$starts > q_end)
  sum(overlaps)
}

#––– Robust Jaccard for half‑open BED (narrowPeak) inputs –––
calculateJaccardScore <- function(gr1, gr2) {
  # sanity check
  if (!inherits(gr1, "GRanges") || !inherits(gr2, "GRanges")) {
    stop("Inputs must be GRanges objects")
  }

  # 1) convert loaded [start,end) → IRanges closed [start+1, end]
  gr1_ho <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr1),
    ranges   = IRanges::IRanges(
      start = IRanges::start(gr1) + 1,
      end   = IRanges::end(gr1)
    )
  )
  gr2_ho <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr2),
    ranges   = IRanges::IRanges(
      start = IRanges::start(gr2) + 1,
      end   = IRanges::end(gr2)
    )
  )

  # 2) merge only overlapping (no abutting!)
  gr1r <- GenomicRanges::reduce(gr1_ho,
                                ignore.strand = TRUE,
                                min.gapwidth  = 0)
  gr2r <- GenomicRanges::reduce(gr2_ho,
                                ignore.strand = TRUE,
                                min.gapwidth  = 0)

  # 3) intersect
  intr <- GenomicRanges::intersect(gr1r, gr2r,
                                   ignore.strand = TRUE)

  # 4) measure true half‑open widths: (end – start)
  w1 <- sum(IRanges::end(gr1r) - IRanges::start(gr1r))
  w2 <- sum(IRanges::end(gr2r) - IRanges::start(gr2r))
  wi <- sum(IRanges::end(intr)  - IRanges::start(intr))

  # 5) union by formula
  uw <- w1 + w2 - wi
  if (uw == 0) return(0)

  wi / uw
}

#––– Pairwise matrix using calculateJaccardScore –––
pairwiseJaccard <- function(lst) {
  n <- length(lst)
  mat <- matrix(NA, nrow=n, ncol=n)
  rownames(mat) <- colnames(mat) <- names(lst)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i,j] <- signif(calculateJaccardScore(lst[[i]], lst[[j]]), 3)
    }
  }
  mat
}
