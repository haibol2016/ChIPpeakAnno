#' Get downstream genomic coordinates
#' 
#' @description 
#' Returns a \code{\link[GenomicRanges]{GRanges}} object containing downstream
#' ranges relative to the input ranges. The function extracts regions downstream
#' of the input ranges, taking strand information into account. For
#' minus-strand ranges, "downstream" refers to the direction opposite to the
#' gene's transcription direction (which is upstream in genomic coordinates).
#' 
#' The output range is defined differently based on strand:
#' 
#' \itemize{
#'   \item \strong{Plus strand (+)} and \strong{unstranded (*)}: 
#'         The downstream region extends from the end of the input range.
#'         Coordinates: \code{(end(gr) - upstream)} to \code{(end(gr) + downstream - 1)}
#'   \item \strong{Minus strand (-)}: 
#'         The downstream region extends from the start of the input range
#'         (which is the 3' end for minus-strand features).
#'         Coordinates: \code{(start(gr) - downstream + 1)} to \code{(start(gr) + upstream)}
#' }
#' 
#' @param gr A \code{\link[GenomicRanges]{GRanges}} object. The strand levels
#'        must be exactly \code{c("+", "-", "*")} in that order. The function
#'        will stop with an error if the strand levels differ.
#' @param upstream A non-negative integer (length 1) specifying the number of
#'        base pairs upstream to include in the output range. For plus-strand
#'        ranges, this extends backward from the end; for minus-strand ranges,
#'        this extends forward from the start.
#' @param downstream A non-negative integer (length 1) specifying the number
#'        of base pairs downstream to include in the output range. For
#'        plus-strand ranges, this extends forward from the end; for
#'        minus-strand ranges, this extends backward from the start.
#' 
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object of the same
#'        length as \code{gr}, containing downstream coordinates. The returned
#'        object:
#'        \itemize{
#'          \item Has the same number of ranges as the input
#'          \item Preserves all metadata columns from the input
#'          \item Has strand levels in the order \code{c("+", "-", "*")}
#'          \item May extend beyond chromosome boundaries (ranges are not
#'                automatically trimmed)
#'        }
#'        \strong{Note:} The returned ranges might extend beyond chromosome
#'        boundaries and should be trimmed using \code{\link[GenomicRanges]{trim}}
#'        if needed.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' The function uses a clever trick with the \code{\link[GenomicRanges]{promoters}}
#' function to extract downstream regions:
#' \enumerate{
#'   \item Reverses the strand levels of the input (swaps "+" and "-", keeps "*")
#'   \item Calls \code{promoters()} with reversed parameters:
#'         \code{promoters(gr_rev, upstream = downstream, downstream = upstream)}
#'   \item Reverses the strand levels back to the original order
#'   \item Returns the result
#' }
#' This approach leverages the fact that promoters extract upstream regions, so
#' by reversing strands and swapping parameters, we get downstream regions.
#' 
#' \strong{Strand-aware behavior:}
#' The function correctly handles strand information:
#' \itemize{
#'   \item \strong{Plus strand (+)}: Downstream is in the positive genomic
#'         direction (toward higher coordinates), starting from the end of the
#'         range
#'   \item \strong{Minus strand (-)}: Downstream is in the negative genomic
#'         direction (toward lower coordinates), starting from the start of the
#'         range (which is the 3' end for minus-strand features)
#'   \item \strong{Unstranded (*)}: Treated the same as plus strand
#' }
#' 
#' \strong{Coordinate calculation examples:}
#' For a plus-strand range at positions 100-105:
#' \itemize{
#'   \item With \code{upstream = 2}, \code{downstream = 3}:
#'         Output: 103-107 (from 105-2 to 105+3-1)
#' }
#' 
#' For a minus-strand range at positions 100-105:
#' \itemize{
#'   \item With \code{upstream = 2}, \code{downstream = 3}:
#'         Output: 98-102 (from 100-3+1 to 100+2)
#' }
#' 
#' \strong{Validation:}
#' The function performs strict validation:
#' \itemize{
#'   \item \code{gr} must be a GRanges object
#'   \item \code{upstream} and \code{downstream} must be numeric, length 1,
#'         and non-negative
#'   \item Strand levels must be exactly \code{c("+", "-", "*")} in that order
#'   \item The function stops with an error if any validation fails
#' }
#' 
#' \strong{Use cases:}
#' This function is useful for:
#' \itemize{
#'   \item Extracting downstream regulatory regions (e.g., 3' UTRs, downstream
#'         enhancers)
#'   \item Analyzing sequence features downstream of genomic features
#'   \item Creating flanking regions for downstream analysis
#' }
#' 
#' @note
#' \itemize{
#'   \item The function does not trim ranges that extend beyond chromosome
#'         boundaries; use \code{\link[GenomicRanges]{trim}} if needed
#'   \item Strand levels must be exactly \code{c("+", "-", "*")}; the function
#'         will error if they differ
#'   \item Both \code{upstream} and \code{downstream} can be 0, which would
#'         return a range of width 1 at the boundary
#'   \item The function preserves all metadata columns from the input
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link[GenomicRanges]{promoters}} for upstream region extraction
#'   \item \code{\link[GenomicRanges]{trim}} for trimming ranges to chromosome
#'         boundaries
#'   \item \code{\link[GenomicRanges]{flank}} for extracting flanking regions
#' }
#' 
#' @author Internal utility function
#' @keywords misc
#' @export
#' 
#' @examples
#' 
#' # Example 1: Basic usage with different strands
#' gr <- GRanges("chr1", IRanges(rep(10, 3), width = 6), 
#'               strand = c("+", "-", "*"))
#' 
#' # Get 2 bp upstream and 2 bp downstream of range ends
#' downstreams(gr, upstream = 2L, downstream = 2L)
#' 
#' # Example 2: Extract downstream region only (no upstream)
#' # Plus strand: extends from end(gr) to end(gr) + 5
#' # Minus strand: extends from start(gr) - 5 to start(gr)
#' downstreams(gr, upstream = 0L, downstream = 5L)
#' 
#' # Example 3: Extract region with more upstream than downstream
#' downstreams(gr, upstream = 10L, downstream = 5L)
#' 
#' # Example 4: With metadata columns (preserved in output)
#' gr_with_meta <- GRanges("chr1", IRanges(10, width = 6), 
#'                         strand = "+",
#'                         gene_id = "G1", score = 100)
#' downstreams(gr_with_meta, upstream = 2L, downstream = 2L)
#' 
#' # Example 5: Multiple ranges
#' gr_multi <- GRanges("chr1", 
#'                     IRanges(start = c(100, 200, 300), width = 10),
#'                     strand = c("+", "-", "+"))
#' downstreams(gr_multi, upstream = 5L, downstream = 10L)
downstreams <- function(gr, upstream, downstream) {
    stopifnot(
        inherits(gr, "GRanges"),
        is.numeric(upstream),
        length(upstream) == 1L,
        upstream >= 0L,
        is.numeric(downstream),
        length(downstream) == 1L,
        downstream >= 0L
    )
    
    stopifnot(identical(levels(strand(gr)), c("+", "-", "*")))
    
    gr_rev <- gr
    levels(strand(gr_rev)) <- c("-", "+", "*")
    out <- promoters(gr_rev, upstream = downstream, downstream = upstream)
    levels(strand(out)) <- c("-", "+", "*")
    
    # Auto-correct strand levels
    stopifnot(identical(levels(strand(out)), c("+", "-", "*")))
    
    out
}
