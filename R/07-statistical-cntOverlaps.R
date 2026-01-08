#' Count overlaps between two GRanges objects
#' 
#' @description 
#' Counts the number of ranges in the first GRanges object that overlap with
#' ranges in the second GRanges object. Optionally allows a maximum gap between
#' ranges to be considered as overlapping. This function is a wrapper around
#' \code{\link[regioneR]{numOverlaps}} that provides gap-based expansion
#' functionality by expanding ranges in the first object before counting overlaps.
#' 
#' This function is useful for statistical analyses where you need to count how
#' many peaks or regions from one set overlap with regions from another set,
#' potentially allowing for small gaps between ranges.
#' 
#' @param A,B \code{\link[GenomicRanges:GRanges-class]{GRanges}} objects. The
#'        function counts how many ranges in \code{A} overlap with ranges in
#'        \code{B}. Overlap is determined after optionally expanding ranges in
#'        \code{A} by \code{maxgap}.
#' @param maxgap A single non-negative integer. Maximum gap (in base pairs)
#'        allowed between ranges for them to be considered overlapping. Default
#'        is \code{0L} (ranges must directly overlap with no gap).
#'        \itemize{
#'          \item \code{maxgap = 0L}: Only ranges that directly overlap (share
#'                at least one base pair) are counted
#'          \item \code{maxgap > 0L}: Ranges in \code{A} are expanded by
#'                \code{maxgap} base pairs on both sides before checking for
#'                overlaps. This allows counting overlaps even when there's a
#'                gap up to \code{maxgap} bp between the original ranges
#'        }
#' @param ... Additional parameters passed to
#'        \code{\link[regioneR]{numOverlaps}} for controlling overlap detection
#'        (e.g., \code{count.once}).
#' 
#' @return Returns an integer representing the number of ranges in \code{A} that
#'        overlap with at least one range in \code{B} (within the specified
#'        \code{maxgap}). Each range in \code{A} is counted at most once,
#'        regardless of how many ranges in \code{B} it overlaps with.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item If \code{maxgap > 0L}, ranges in \code{A} are expanded symmetrically
#'         by \code{maxgap} base pairs on both sides using \code{\link{expandGR}}.
#'         Expanded ranges are trimmed to valid genomic coordinates.
#'   \item Overlaps between (expanded) ranges in \code{A} and ranges in \code{B}
#'         are counted using \code{\link[regioneR]{numOverlaps}}.
#'   \item Returns the count of ranges in \code{A} that overlap with at least
#'         one range in \code{B}.
#' }
#' 
#' \strong{Gap-based expansion:}
#' When \code{maxgap > 0}, the function expands each range in \code{A} by
#' \code{maxgap} bp upstream and downstream before checking for overlaps. This
#' means:
#' \itemize{
#'   \item A range in \code{A} at position [100, 200] with \code{maxgap = 20}
#'         becomes [80, 220]
#'   \item If the expanded range [80, 220] overlaps with any range in \code{B},
#'         the original range [100, 200] is counted as overlapping
#'   \item This effectively allows overlaps with gaps up to \code{maxgap} bp
#'         between the original ranges
#' }
#' 
#' \strong{Difference from findOverlaps:}
#' Unlike \code{\link[GenomicRanges]{findOverlaps}}, which returns detailed
#' overlap information, this function only returns a count. It's optimized for
#' cases where you only need to know how many ranges overlap, not which specific
#' ranges overlap.
#' 
#' \strong{Use cases:}
#' \itemize{
#'   \item Statistical testing: Count how many peaks overlap with annotation
#'         features
#'   \item Quality control: Check overlap rates between replicates
#'   \item Permutation testing: Count overlaps in random vs. observed data
#' }
#' 
#' @seealso \code{\link[regioneR]{numOverlaps}} for the underlying overlap
#'          counting function, \code{\link{expandGR}} for range expansion,
#'          \code{\link[GenomicRanges]{findOverlaps}} for detailed overlap
#'          information
#' 
#' @author Internal utility function
#' @keywords misc
#' @export
#' @importFrom regioneR numOverlaps
#' @examples
#' ## Example 1: Count direct overlaps
#' A <- GRanges("chr1", IRanges(c(100, 200, 300), width = 50))
#' B <- GRanges("chr1", IRanges(c(150, 250), width = 50))
#' cntOverlaps(A, B, maxgap = 0L)
#' # Returns: 2 (ranges at 200 and 300 overlap with B)
#' 
#' ## Example 2: Count overlaps allowing 20 bp gap
#' A <- GRanges("chr1", IRanges(c(100, 200, 300), width = 50))
#' B <- GRanges("chr1", IRanges(c(175, 275), width = 50))
#' cntOverlaps(A, B, maxgap = 0L)   # No direct overlaps
#' cntOverlaps(A, B, maxgap = 20L)  # Overlaps with gap allowed
#' 
#' ## Example 3: Count peak overlaps with annotation features
#' peaks <- GRanges("chr1", IRanges(c(1000, 5000, 10000), width = 200))
#' genes <- GRanges("chr1", IRanges(c(800, 4800, 15000), width = 2000))
#' cntOverlaps(peaks, genes, maxgap = 1000L)
#' # Counts peaks within 1kb of gene boundaries
cntOverlaps <- function(A, B, maxgap = 0L, ...) {
    if (maxgap > 0L) {
        A <- expandGR(A, maxgap)
    }
    
    numOverlaps(A, B, ...)
}
