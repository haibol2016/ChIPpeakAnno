#' Slide windows over genomic regions
#' 
#' @description 
#' Creates a set of overlapping or non-overlapping genomic windows (tiles) by
#' sliding a fixed-size window across input regions with a specified step size.
#' This function is useful for binning genomic regions for coverage analysis,
#' signal aggregation, metagene plots, or other window-based operations.
#' 
#' The function creates windows by:
#' \enumerate{
#'   \item Dividing each region into tiles of size \code{step}
#'   \item Expanding each tile to size \code{windowSize} (creating overlaps)
#'   \item Filtering out windows that extend beyond region boundaries
#'   \item Tracking which original region each window came from
#' }
#' 
#' @param targetRegions A \code{\link[GenomicRanges]{GRanges}} object
#'        containing genomic regions to be tiled (e.g., genes, peaks, exons).
#'        Regions smaller than \code{windowSize} or \code{step} are removed
#'        with a warning. All metadata columns are preserved in the output.
#' @param windowSize A positive integer (length 1) specifying the size of each
#'        window (tile) in base pairs. Each window will have this width (except
#'        partial windows at region ends if \code{keepPartialWindow = TRUE}).
#'        Must be > 0.
#' @param step A positive integer (length 1) specifying the step size for
#'        sliding windows in base pairs. This determines the spacing between
#'        window starts. Must be > 0.
#'        \itemize{
#'          \item If \code{step < windowSize}: Windows overlap by
#'                \code{windowSize - step} base pairs
#'          \item If \code{step == windowSize}: Windows are adjacent
#'                (non-overlapping)
#'          \item If \code{step > windowSize}: Windows have gaps between them
#'        }
#' @param keepPartialWindow A logical value. If \code{TRUE}, includes the last
#'        partial window at the end of each region (even if it's smaller than
#'        \code{windowSize}). The partial window is trimmed to the region end.
#'        If \code{FALSE} (default), partial windows that extend beyond region
#'        boundaries are excluded.
#' @param ... Additional arguments (currently not used).
#' 
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object containing
#'        tiled windows. The returned object:
#'        \itemize{
#'          \item Has one window per step position within each input region
#'          \item Preserves chromosome, strand, and other metadata from input
#'                regions
#'          \item Contains a metadata column \code{oid} (original ID) indicating
#'                which input region each window came from (1-based index)
#'          \item Windows have width equal to \code{windowSize} (except partial
#'                windows if \code{keepPartialWindow = TRUE})
#'          \item Partial windows extending beyond region boundaries are excluded
#'                unless \code{keepPartialWindow = TRUE}
#'        }
#'        The number of windows per region depends on the region width,
#'        \code{windowSize}, and \code{step}.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates input: ensures \code{targetRegions} is a GRanges object,
#'         and \code{windowSize} and \code{step} are positive integers
#'   \item Filters regions: removes regions smaller than \code{windowSize} or
#'         \code{step} (with a warning)
#'   \item Creates initial tiles: uses \code{\link[IRanges]{tile}} to divide
#'         each region into tiles of size \code{step}
#'   \item Expands tiles: expands each tile to size \code{windowSize} (this
#'         creates overlaps when \code{step < windowSize})
#'   \item Handles boundaries: identifies windows that extend beyond region
#'         ends:
#'         \itemize{
#'           \item If \code{keepPartialWindow = TRUE}: Trims partial windows to
#'                 region end
#'           \item If \code{keepPartialWindow = FALSE}: Removes partial windows
#'         }
#'   \item Adds metadata: assigns \code{oid} to track original region indices
#'   \item Returns tiled windows
#' }
#' 
#' \strong{Window creation example:}
#' For a 200 bp region with \code{windowSize = 50} and \code{step = 10}:
#' \itemize{
#'   \item Initial tiles (step=10): positions 1-10, 11-20, 21-30, ..., 191-200
#'   \item Expanded to windowSize=50: positions 1-50, 11-60, 21-70, ..., 151-200
#'   \item Windows overlap by 40 bp (50 - 10 = 40)
#'   \item Total: 20 windows (200 / 10 = 20)
#' }
#' 
#' For a 200 bp region with \code{windowSize = 50} and \code{step = 50}:
#' \itemize{
#'   \item Windows: positions 1-50, 51-100, 101-150, 151-200
#'   \item Windows are adjacent (non-overlapping)
#'   \item Total: 4 windows (200 / 50 = 4)
#' }
#' 
#' \strong{Overlap calculation:}
#' When \code{step < windowSize}, windows overlap by:
#' \preformatted{
#' overlap = windowSize - step
#' }
#' For example, with \code{windowSize = 100} and \code{step = 25}, windows
#' overlap by 75 bp.
#' 
#' \strong{Partial window handling:}
#' Partial windows occur when a window would extend beyond a region's end:
#' \itemize{
#'   \item \code{keepPartialWindow = FALSE} (default): Partial windows are
#'         excluded. This ensures all windows have the same size.
#'   \item \code{keepPartialWindow = TRUE}: Partial windows are included but
#'         trimmed to the region end. This ensures complete coverage of the
#'         region but creates windows of different sizes.
#' }
#' 
#' \strong{Region filtering:}
#' Regions smaller than \code{windowSize} or \code{step} are automatically
#' removed with a warning. This ensures that:
#' \itemize{
#'   \item At least one full window can be created per region
#'   \item The tiling process works correctly
#' }
#' 
#' \strong{Use cases:}
#' \itemize{
#'   \item \strong{Metagene plots}: Create windows across genes to visualize
#'         average signal profiles
#'   \item \strong{Coverage analysis}: Bin large regions for efficient coverage
#'         calculation
#'   \item \strong{Signal aggregation}: Summarize signal within fixed-size
#'         windows
#'   \item \strong{Feature alignment}: Align features of different sizes to
#'         common window structure
#' }
#' 
#' @note
#' \itemize{
#'   \item Regions smaller than \code{windowSize} or \code{step} are removed
#'         with a warning
#'   \item The \code{oid} metadata column uses 1-based indexing (first region
#'         has oid=1)
#'   \item Windows preserve all metadata from their parent regions
#'   \item When \code{step > windowSize}, windows have gaps between them
#'   \item The function uses \code{\link[IRanges]{tile}} internally for initial
#'         tiling
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link[IRanges]{tile}} for the underlying tiling function
#'   \item \code{\link{summarizeOverlapsByBins}} for counting reads in binned
#'         features
#'   \item \code{\link{tileCount}} for genome-wide tiling windows
#' }
#' 
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @examples
#' 
#' # Example 1: Basic usage with overlapping windows
#' genes <- GRanges(
#'     seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
#'     ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 
#'                        4000, 7500, 5000, 5400), 
#'                      width=c(rep(500, 3), 600, 900, 500, 300, 900, 
#'                              300, 500, 500),
#'                      names=letters[1:11])) 
#' windows <- tileGRanges(genes, windowSize=50, step=10)
#' # View windows and their original region IDs
#' windows
#' mcols(windows)$oid  # Original region indices
#' 
#' # Example 2: Non-overlapping windows (step == windowSize)
#' windows <- tileGRanges(genes, windowSize=50, step=50)
#' # Windows are adjacent, no overlap
#' 
#' # Example 3: Include partial windows at region ends
#' windows <- tileGRanges(genes, windowSize=50, step=10,
#'                        keepPartialWindow = TRUE)
#' # Last window in each region may be smaller than windowSize
#' 
#' # Example 4: Larger windows with larger step
#' windows <- tileGRanges(genes, windowSize=200, step=100)
#' # Windows overlap by 100 bp
#' 
#' # Example 5: Small windows for fine resolution
#' windows <- tileGRanges(genes, windowSize=20, step=5)
#' # Many overlapping windows for detailed analysis
#' 
#' # Example 6: Regions with metadata (preserved in output)
#' genes_with_meta <- GRanges(
#'     "chr1", IRanges(c(1000, 5000), width=500),
#'     gene_id = c("G1", "G2"),
#'     score = c(10, 20))
#' windows <- tileGRanges(genes_with_meta, windowSize=50, step=25)
#' # Metadata is preserved in windows
#' mcols(windows)
#' 
#' # Example 7: Handling small regions (will be removed with warning)
#' small_regions <- GRanges("chr1", IRanges(1000, width=30))
#' # Warning: region smaller than windowSize/step will be removed
#' windows <- tileGRanges(small_regions, windowSize=50, step=10)
#' 
tileGRanges <- function(targetRegions,
                       windowSize,
                       step,
                       keepPartialWindow = FALSE,
                       ...) {
    if (!inherits(targetRegions, "GRanges")) {
        stop("'targetRegions' must be a GRanges object", call. = FALSE)
    }
    
    if (!is.numeric(windowSize) || length(windowSize) != 1L || windowSize <= 0L) {
        stop("'windowSize' must be a positive integer", call. = FALSE)
    }
    
    if (!is.numeric(step) || length(step) != 1L || step <= 0L) {
        stop("'step' must be a positive integer", call. = FALSE)
    }
    
    if (any(width(targetRegions) < windowSize, na.rm = TRUE) ||
        any(width(targetRegions) < step, na.rm = TRUE)) {
        warning("Some regions are smaller than windowSize or step. ",
                "They will be removed.",
                call. = FALSE,
                immediate. = TRUE)
    }
    targetRegions <- targetRegions[width(targetRegions) >= windowSize]
    targetRegions <- targetRegions[width(targetRegions) >= step]
    
    tileTargetRanges <- tile(x = ranges(targetRegions), width = step)
    nt <- elementNROWS(tileTargetRanges)
    tileTargetRanges_end <- rep(end(targetRegions), nt)
    tileTargetRanges <- unlist(tileTargetRanges)
    width(tileTargetRanges) <- windowSize
    
    tileTargetRegions <- GRanges(
        rep(seqnames(targetRegions), nt),
        tileTargetRanges,
        rep(strand(targetRegions), nt),
        oid = rep(seq_along(targetRegions), nt)
    )
    
    id <- end(tileTargetRanges) > tileTargetRanges_end
    if (keepPartialWindow) {
        end(tileTargetRegions[id]) <- tileTargetRanges_end[id]
    } else {
        tileTargetRegions <- tileTargetRegions[!id]
    }
    
    tileTargetRegions
}
