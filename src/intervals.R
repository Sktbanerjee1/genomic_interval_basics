library(GenomicRanges)
library(IRanges)

# Generates a random GRanges object
#
# You can use this function to get a random interval list for testing.
# Don't change this function, it's just for you to use while
# you work on the other functions
#' @examples
# gr1 = randomGRanges()
randomGRanges = function(){
  starts = sample(1:1e7, 500)
  sizes = sample(1:1e4, 500)
  ends = starts + sizes
  ret = GenomicRanges::GRanges(seqnames="chrX", ranges=IRanges::IRanges(starts,ends))
  GenomicRanges::reduce(ret)
}

# Counts overlaps between query and database
#
# This function counts the number of intervals from a
# database set of intervals that overlap with a single
# query interval.
#
#' @param query
#'   Interval to check for overlaps. Should be a list
#'   object with 'start' and 'end' values.
#' @param database
#'   A data.frame, with rows corresponding to intervals,
#'   with 'starts' and 'ends' as columns.
#' @return
#'   Number of overlaps counted (as a numeric)
#' @examples
#'   query = list(start=50, end=60)
#'   database = data.frame(starts=c(35, 45, 55), ends=c(45, 65, 85))
#'   countOverlapsSimple(query, database)  # returns 2
countOverlapsSimple = function(query, database) {
  q_start <- query$start
  q_end <- query$end

  # Identify intervals in the database that overlap with the query
  overlaps <- !(database$ends < q_start | database$starts > q_end)

  # Return number of overlapping intervals
  return(sum(overlaps))
}

# Measure the Jaccard similarity between two interval sets
#
# This function calculates the Jaccard Score between two
# given interval sets, provided as GenomicRanges::GRanges
# objects. The Jaccard score is computed as the intersection
# between the two sets divided by the union of the two sets.
# @param gr1  First GRanges object
# @param gr2  Second GRanges object
# @return The Jaccard score (as numeric)
calculateJaccardScore <- function(gr1, gr2) {
  # 1) sanity check
  if (!inherits(gr1, "GRanges") || !inherits(gr2, "GRanges")) {
    stop("Inputs must be GRanges objects")
  }

  # 2) convert from 0‑based [start,end) → 1‑based closed IRanges for half‑open math
  gr1_ho <- GenomicRanges::GRanges(
    seqnames = seqnames(gr1),
    ranges   = IRanges::IRanges(start(gr1) + 1, end(gr1))
  )
  gr2_ho <- GenomicRanges::GRanges(
    seqnames = seqnames(gr2),
    ranges   = IRanges::IRanges(start(gr2) + 1, end(gr2))
  )

  # 3) merge only truly overlapping (not abutting) intervals
  gr1r <- GenomicRanges::reduce(gr1_ho, ignore.strand = TRUE, min.gapwidth = 0)
  gr2r <- GenomicRanges::reduce(gr2_ho, ignore.strand = TRUE, min.gapwidth = 0)

  # 4) compute the intersection
  intr <- GenomicRanges::intersect(gr1r, gr2r, ignore.strand = TRUE)

  # 5) half‑open widths = end – start
  w1 <- sum(end(gr1r) - start(gr1r))
  w2 <- sum(end(gr2r) - start(gr2r))
  wi <- sum(end(intr)  - start(intr))

  # 6) union width by formula
  uw <- w1 + w2 - wi
  if (uw == 0) return(0)

  # 7) Jaccard = |A∩B| / |A∪B|
  wi / uw
}




# Calculate pairwise Jaccard similarity among several interval sets
#
# This function makes use of \code{calculateJaccardScore}. It simply
# loops through each pairwise comparison and calculates.
#' Round the result to 3 significant figures using \code{signif}
#' @param lst A base R list of GRanges objects to compare
#' @return
#'   A matrix of size n-by-n, where n is the number of elements in lst.
#' @examples
#' lst = replicate(10, randomGRanges())
#' pairwiseJaccard(lst)
pairwiseJaccard = function(lst) {
  n <- length(lst)
  result <- matrix(NA, nrow = n, ncol = n)
  rownames(result) <- colnames(result) <- names(lst)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      score <- calculateJaccardScore(lst[[i]], lst[[j]])
      result[i, j] <- signif(score, 3)
    }
  }
  return(result)
}
