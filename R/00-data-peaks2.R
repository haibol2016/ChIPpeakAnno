#' Example ChIP-seq peak dataset 2
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing example ChIP-seq peaks
#' for demonstration purposes. This is one of three example peak datasets
#' (\code{peaks1}, \code{peaks2}, \code{peaks3}) that can be used to demonstrate
#' overlap analysis, Venn diagram generation, and multi-set peak comparisons.
#' 
#' @name peaks2
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
#' \code{peaks3} to demonstrate:
#' \itemize{
#'   \item Multi-set peak overlap analysis with \code{\link{findOverlapsOfPeaks}}
#'   \item Venn diagram generation with \code{\link{makeVennDiagram}}
#'   \item Statistical testing of overlaps with \code{\link{peakPermTest}}
#' }
#' 
#' @seealso
#' \code{\link{peaks1}}, \code{\link{peaks3}}, \code{\link{findOverlapsOfPeaks}},
#' \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(peaks2)
#' 
#' # Inspect the structure
#' peaks2
#' length(peaks2)
#' 
#' # View first few peaks
#' head(peaks2, n = 2)
"peaks2"
