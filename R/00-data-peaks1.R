#' Example ChIP-seq peak dataset 1
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing example ChIP-seq peaks
#' for demonstration purposes. This is one of three example peak datasets
#' (\code{peaks1}, \code{peaks2}, \code{peaks3}) that can be used to demonstrate
#' overlap analysis, Venn diagram generation, and multi-set peak comparisons.
#' 
#' @name peaks1
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Peak identifiers as character vector}
#' }
#' 
#' @details
#' This dataset contains 12 example peaks that can be used with \code{peaks2} and
#' \code{peaks3} to demonstrate:
#' \itemize{
#'   \item Multi-set peak overlap analysis with \code{\link{findOverlapsOfPeaks}}
#'   \item Venn diagram generation with \code{\link{makeVennDiagram}}
#'   \item Statistical testing of overlaps with \code{\link{peakPermTest}}
#' }
#' 
#' @seealso
#' \code{\link{peaks2}}, \code{\link{peaks3}}, \code{\link{findOverlapsOfPeaks}},
#' \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(peaks1)
#' 
#' # Inspect the structure
#' peaks1
#' length(peaks1)
#' 
#' # View first few peaks
#' head(peaks1, n = 2)
#' 
#' # Compare with other peak sets
#' data(peaks2)
#' data(peaks3)
#' 
#' # Find overlaps
#' overlaps <- findOverlapsOfPeaks(peaks1, peaks2, peaks3)
"peaks1"
