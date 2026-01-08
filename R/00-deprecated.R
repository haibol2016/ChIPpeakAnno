#' Deprecated Functions and Features in Package ChIPpeakAnno
#' 
#' This page documents functions and features that are deprecated in
#' \pkg{ChIPpeakAnno}. Deprecated items are maintained for backward
#' compatibility but will be removed in a future release. Users are strongly
#' encouraged to migrate to the recommended alternatives.
#' 
#' @section Deprecated Functions:
#' 
#' \strong{findOverlappingPeaks}
#' 
#' The function \code{findOverlappingPeaks} is deprecated. Use
#' \code{\link{findOverlapsOfPeaks}} instead, which provides the same
#' functionality with improved performance and supports multiple peak sets.
#' 
#' \strong{Migration:}
#' \preformatted{
#' # Old (deprecated):
#' result <- findOverlappingPeaks(peaks1, peaks2, ...)
#' 
#' # New (recommended):
#' result <- findOverlapsOfPeaks(peaks1, peaks2, ...)
#' }
#' 
#' @section Deprecated Data Objects:
#' 
#' \strong{Pre-computed TSS Datasets (Legacy Assemblies)}
#' 
#' Several pre-computed TSS annotation datasets for legacy genome assemblies
#' are deprecated and will be removed in a future release:
#' \itemize{
#'   \item \code{TSS.human.NCBI36} (hg18)
#'   \item \code{TSS.human.GRCh38} (hg38) - kept for examples only
#'   \item \code{TSS.mouse.NCBIM37} (mm9)
#'   \item \code{TSS.rat.RGSC3.4} (rn4)
#'   \item \code{TSS.rat.Rnor_5.0} (rn5)
#'   \item \code{TSS.zebrafish.Zv8} (danRer7)
#'   \item \code{TSS.zebrafish.Zv9} (danRer10)
#' }
#' 
#' \strong{Note:} The datasets \code{TSS.human.GRCh37} and
#' \code{TSS.mouse.GRCm38} are retained for package examples and unit testing
#' only. Users should not use pre-computed datasets for their own analyses.
#' 
#' \strong{Migration:}
#' 
#' Generate annotations matching your genome assembly using one of these
#' approaches:
#' 
#' \preformatted{
#' # Option 1: Using biomaRt (recommended)
#' library(biomaRt)
#' mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#' TSS <- getAnnotation(mart, featureType = "TSS")
#' 
#' # Option 2: Using EnsDb packages (alternative)
#' library(EnsDb.Hsapiens.v86)
#' annoData <- annoGR(EnsDb.Hsapiens.v86)
#' 
#' # Option 3: Using TxDb packages
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' }
#' 
#' @rdname ChIPpeakAnno-deprecated
#' 
#' @param Peaks1 \link[GenomicRanges:GRanges-class]{GRanges} object containing
#' the first set of peaks. See example below.
#' @param Peaks2 \link[GenomicRanges:GRanges-class]{GRanges} object containing
#' the second set of peaks. See example below.
#' @param maxgap,minoverlap Integer values used in the internal call to
#' \code{\link[IRanges:findOverlaps-methods]{findOverlaps}} to detect overlaps.
#' \code{maxgap} is the maximum gap allowed between ranges for them to be
#' considered overlapping. \code{minoverlap} is the minimum overlap required.
#' See \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} in the
#' \pkg{IRanges} package for details.
#' @param multiple Logical. \code{TRUE} may return multiple overlapping peaks
#' in Peaks2 for one peak in Peaks1; \code{FALSE} will return at most one
#' overlapping peak in Peaks2 for one peak in Peaks1. This parameter is kept
#' for backward compatibility only. Please use \code{select} instead.
#' @param NameOfPeaks1 Character string. Name of Peaks1, used for generating
#' column names in the output.
#' @param NameOfPeaks2 Character string. Name of Peaks2, used for generating
#' column names in the output.
#' @param select Character. Controls which overlapping peaks are returned:
#' \describe{
#'   \item{"all"}{Return all overlapping peaks}
#'   \item{"first"}{Return the first overlapping peak}
#'   \item{"last"}{Return the last overlapping peak}
#'   \item{"arbitrary"}{Return one arbitrary overlapping peak}
#' }
#' @param annotate Integer. Whether to include overlapFeature and
#' shortestDistance in the OverlappingPeaks output. \code{1} means yes,
#' \code{0} means no. Default is \code{0}.
#' @param ignore.strand Logical. When \code{TRUE}, strand information is ignored
#' in the overlap calculations. Default is \code{TRUE}.
#' @param connectedPeaks Character. Controls how multiple peaks involved in
#' overlapping groups are handled:
#' \describe{
#'   \item{"merge"}{Count connected peaks as only 1}
#'   \item{"min"}{Count as the minimal number of involved peaks in any
#'         concerned group}
#' }
#' @param \dots Additional \link[GenomicRanges:GRanges-class]{GRanges} objects
#' (for \code{findOverlapsOfPeaks}).
#' 
#' @return For \code{findOverlappingPeaks}:
#' \describe{
#'   \item{OverlappingPeaks}{A data frame containing input peaks information
#'         with added columns: \code{overlapFeature} (spatial relationship
#'         between peaks) and \code{shortestDistance} (shortest distance between
#'         overlapping peaks)}
#'   \item{MergedPeaks}{A \link[GenomicRanges:GRanges-class]{GRanges} object
#'         containing merged overlapping peaks}
#' }
#' 
#' @seealso 
#' \itemize{
#'   \item \code{\link{findOverlapsOfPeaks}} - Recommended replacement for
#'         \code{findOverlappingPeaks}
#'   \item \code{\link{getAnnotation}} - Generate TSS annotations matching your
#'         genome assembly
#'   \item \code{\link{annoGR}} - Convert EnsDb or TxDb objects to GRanges
#'   \item \code{\link{toGRanges}} - Convert peak files to GRanges format
#' }
#' 
#' @references
#' For information about deprecated functions in R, see
#' \code{\link[base:Deprecated]{Deprecated}}.
#' 
#' @name ChIPpeakAnno-deprecated
#' 
#' @examples
#' 
#' ## Example of deprecated function usage (for reference only)
#' \dontrun{
#' if (interactive()) {
#'     peaks1 <- GRanges(
#'         seqnames = c(6, 6, 6, 6, 5),
#'         IRanges(
#'             start = c(1543200, 1557200, 1563000, 1569800, 167889600),
#'             end = c(1555199, 1560599, 1565199, 1573799, 167893599),
#'             names = c("p1", "p2", "p3", "p4", "p5")
#'         ),
#'         strand = "+"
#'     )
#'     peaks2 <- GRanges(
#'         seqnames = c(6, 6, 6, 6, 5),
#'         IRanges(
#'             start = c(1549800, 1554400, 1565000, 1569400, 167888600),
#'             end = c(1550599, 1560799, 1565399, 1571199, 167888999),
#'             names = c("f1", "f2", "f3", "f4", "f5")
#'         ),
#'         strand = "+"
#'     )
#'     
#'     # Deprecated function (shows warning)
#'     result <- findOverlappingPeaks(
#'         peaks1, peaks2,
#'         maxgap = 1000,
#'         NameOfPeaks1 = "TF1",
#'         NameOfPeaks2 = "TF2",
#'         select = "all",
#'         annotate = 1
#'     )
#'     
#'     # Recommended replacement
#'     result <- findOverlapsOfPeaks(
#'         peaks1, peaks2,
#'         maxgap = 1000,
#'         NameOfPeaks = c("TF1", "TF2"),
#'         select = "all"
#'     )
#' }
#' }
#' 
NULL
