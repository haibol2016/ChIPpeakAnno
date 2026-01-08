#' Re-center peaks to fixed width
#' 
#' Create a new set of peaks centered on the original peak centers with a
#' specified fixed width. This is useful for standardizing peak sizes for
#' downstream analysis or visualization.
#' 
#' @param peaks A \code{\link[GenomicRanges]{GRanges}} or \code{\link{annoGR}}
#'        object containing peaks to be re-centered.
#' @param width Integer. The width (in base pairs) of the re-centered peaks.
#'        Default is \code{2000L}.
#' @param ... Additional arguments (not currently used).
#' 
#' @return A \code{\link[GenomicRanges]{GRanges}} object with re-centered peaks
#'         of the specified width. The new peaks are centered on the midpoint
#'         of the original peaks.
#' 
#' @details
#' This function:
#' \enumerate{
#'   \item Calculates the center of each peak: \code{start + floor(width/2)}
#'   \item Creates new peaks centered on these positions with the specified width
#'   \item Issues warnings if peaks extend beyond chromosome boundaries
#' }
#' 
#' If sequence length information is available, the function checks for
#' out-of-bound peaks and issues warnings.
#' 
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom GenomeInfoDb seqlengths
#' 
#' @examples
#' peaks <- GRanges("chr1", IRanges(100, 200))
#' 
#' # Re-center to 2 bp width
#' reCenterPeaks(peaks, width = 2L)
#' 
#' # Re-center to 2000 bp width (default)
#' reCenterPeaks(peaks, width = 2000L)
reCenterPeaks <- function(peaks, width = 2000L, ...) {
    stopifnot(
        inherits(peaks, c("annoGR", "GRanges")),
        is.numeric(width),
        length(width) == 1L,
        width > 0L
    )
    
    peaks_center <- start(peaks) + floor(width(peaks) / 2L)
    peaks_recentered <- peaks
    start(peaks_recentered) <- peaks_center - floor(width / 2L)
    width(peaks_recentered) <- as.integer(width)
    
    if (any(start(peaks_recentered) < 1L, na.rm = TRUE)) {
        warning("Some start positions of the peaks are less than 1!", 
                call. = FALSE)
    }
    
    seq_len <- seqlengths(peaks)
    seq_len <- seq_len[!is.na(seq_len)]
    
    if (length(seq_len) > 0L) {
        peaks_subset <- 
            peaks_recentered[seqnames(peaks_recentered) %in% names(seq_len)]
        seqnames_char <- as.character(seqnames(peaks_subset))
        if (any(end(peaks_subset) > seq_len[seqnames_char], na.rm = TRUE)) {
            warning("Some end positions of the peaks are out of bound!", 
                    call. = FALSE)
        }
    }
    
    peaks_recentered
}
