#' Example ChIP-seq peak dataset 3
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing example ChIP-seq peaks
#' for demonstration purposes. This is one of three example peak datasets
#' (\code{peaks1}, \code{peaks2}, \code{peaks3}) that can be used to demonstrate
#' overlap analysis, Venn diagram generation, and multi-set peak comparisons.
#' 
#' @name peaks3
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
#' This dataset contains example peaks that can be used with \code{peaks1} and
#' \code{peaks2} to demonstrate:
#' \itemize{
#'   \item Multi-set peak overlap analysis with \code{\link{findOverlapsOfPeaks}}
#'   \item Venn diagram generation with \code{\link{makeVennDiagram}}
#'   \item Statistical testing of overlaps with \code{\link{peakPermTest}}
#' }
#' 
#' @seealso
#' \code{\link{peaks1}}, \code{\link{peaks2}}, \code{\link{findOverlapsOfPeaks}},
#' \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(peaks3)
#' 
#' # Inspect the structure
#' peaks3
#' length(peaks3)
#' 
#' # View first few peaks
#' head(peaks3, n = 2)
"peaks3"
