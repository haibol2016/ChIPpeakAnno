#' Build binding distribution from peaks
#' 
#' @description 
#' Analyzes the distribution of peak distances from annotation features (e.g.,
#' TSS or gene ends) and creates a \code{\link{bindist}} object representing
#' this distribution. The function calculates the shortest distance from each
#' peak to the nearest annotation feature, bins these distances, and counts
#' peaks in each bin. This distribution is then used to generate matched random
#' background regions for permutation testing, ensuring that random regions
#' follow the same distance distribution as the observed peaks.
#' 
#' The function uses robust statistical methods (boxplot statistics) to
#' automatically determine appropriate bin sizes and ranges, handling outliers
#' and ensuring the distribution captures the main pattern of peak binding
#' relative to features.
#' 
#' @param x A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        representing peaks to analyze. Each range should represent a single
#'        peak region.
#' @param AnnotationData A \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#'        or \code{\link{annoGR}} object containing annotation features (e.g.,
#'        transcripts, exons, or genes). The function will find the nearest
#'        feature to each peak.
#' @param bindingType A character string specifying the feature position to use
#'        for distance calculation. Options:
#'        \itemize{
#'          \item \code{"TSS"} (default): Transcription start site. For
#'                positive strand features, this is the start position. For
#'                negative strand features, this is the end position.
#'          \item \code{"geneEnd"}: Transcription end site. For positive strand
#'                features, this is the end position. For negative strand
#'                features, this is the start position.
#'        }
#' @param featureType A character string specifying the type of annotation
#'        feature to use. Options:
#'        \itemize{
#'          \item \code{"transcript"} (default): Use transcript-level features
#'          \item \code{"exon"}: Use exon-level features
#'        }
#'        Note: This parameter is passed to \code{annotatePeakInBatch} but may
#'        not be directly used if \code{AnnotationData} already contains the
#'        desired feature type.
#' 
#' @return Returns an object of class \code{\link{bindist}} containing:
#'   \itemize{
#'     \item \code{counts}: An integer vector with the number of peaks in each
#'           distance bin. Only bins with non-zero counts are included.
#'     \item \code{mids}: An integer vector of bin midpoints (in base pairs)
#'           corresponding to the distance bins. Same length as \code{counts}.
#'     \item \code{halfBinSize}: An integer specifying half the width of each
#'           bin (in base pairs). The full bin width is \code{2 * halfBinSize}.
#'     \item \code{bindingType}: The \code{bindingType} value used (character)
#'     \item \code{featureType}: The \code{featureType} value used (character)
#'   }
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Annotates each peak in \code{x} to find the nearest feature in
#'         \code{AnnotationData} using \code{annotatePeakInBatch} with
#'         \code{output = "shortestDistance"} and
#'         \code{PeakLocForDistance = "middle"}
#'   \item Calculates the shortest distance from each peak center to the
#'         specified feature position (TSS or geneEnd)
#'   \item Uses \code{boxplot.stats} with coefficient 1.5 to identify the
#'         interquartile range and whiskers, which helps handle outliers
#'   \item Determines bin breaks based on the whisker range:
#'         \itemize{
#'           \item If whisker range > 10,000 bp: Creates 100 bins
#'           \item If whisker range > 100 bp: Creates bins with width 100 bp
#'           \item Otherwise: Creates a single bin centered at the mean
#'         }
#'   \item Extends breaks to include all distances (including outliers)
#'   \item Creates a histogram and extracts counts, midpoints, and bin size
#'   \item Removes bins with zero counts
#'   \item Creates and returns a \code{bindist} object
#' }
#' 
#' \strong{Distance calculation:}
#' The function uses \code{annotatePeakInBatch} with \code{output = "shortestDistance"},
#' which calculates the shortest gap between any boundary of the peak and any
#' boundary of the feature. The peak center is used as the reference point
#' (\code{PeakLocForDistance = "middle"}), and the feature reference point is
#' determined by \code{bindingType} (TSS or geneEnd).
#' 
#' \strong{Binning algorithm:}
#' The binning strategy adapts to the data distribution:
#' \itemize{
#'   \item For wide distributions (>10kb): Uses 100 evenly-spaced bins
#'   \item For medium distributions (100bp-10kb): Uses 100bp-wide bins
#'   \item For narrow distributions (<100bp): Uses a single bin
#' }
#' The minimum bin width is 100 bp. Outliers beyond the whiskers are included
#' by extending the break range to cover all distances.
#' 
#' \strong{Use in permutation testing:}
#' The \code{bindist} object returned by this function is used by
#' \code{\link{preparePool}} to generate random background regions that match
#' the distance distribution of the observed peaks. This ensures that
#' permutation tests compare peaks to appropriately matched random regions,
#' making the statistical test more powerful and appropriate.
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{bindist}} for the class definition,
#'          \code{\link{preparePool}} for generating matched random regions,
#'          \code{\link{peakPermTest}} for performing permutation tests,
#'          \code{\link{annotatePeakInBatch}} for peak annotation
#' @keywords misc
#' @export
#' @importFrom graphics hist
#' @import methods
#' @importFrom grDevices boxplot.stats
#' @examples
#' \dontrun{
#' ## Example 1: Build distribution relative to TSS
#' library(EnsDb.Hsapiens.v75)
#' data("myPeakList")
#' annoGR <- annoGR(EnsDb.Hsapiens.v75, feature = "transcript")
#' dist <- buildBindingDistribution(myPeakList, annoGR, 
#'                                  bindingType = "TSS",
#'                                  featureType = "transcript")
#' dist
#' 
#' ## Example 2: Build distribution relative to gene ends
#' dist <- buildBindingDistribution(myPeakList, annoGR,
#'                                  bindingType = "geneEnd",
#'                                  featureType = "transcript")
#' 
#' ## Example 3: Use with preparePool for permutation testing
#' pool <- preparePool(annoGR, dist)
#' result <- peakPermTest(myPeakList, pool, TxDb = EnsDb.Hsapiens.v75)
#' }
buildBindingDistribution <- function(x, AnnotationData,
                                     bindingType = c("TSS", "geneEnd"), 
                                     featureType = c("transcript", "exon")) {
    if (missing(AnnotationData) || missing(x)) {
        stop("'x' and 'AnnotationData' are required", call. = FALSE)
    }
    if (!inherits(x, "GRanges")) {
        stop("'x' must be a GRanges object", call. = FALSE)
    }
    if (inherits(AnnotationData, "annoGR")) {
        AnnotationData <- as(AnnotationData, "GRanges")
    }
    if (!inherits(AnnotationData, "GRanges")) {
        stop("'AnnotationData' must be a GRanges object", call. = FALSE)
    }
    
    bindingType <- match.arg(bindingType)
    featureType <- match.arg(featureType)
    AnnotationData <- unique(AnnotationData)
    
    suppressWarnings({
        anno <- annotatePeakInBatch(x, AnnotationData = AnnotationData, 
                                    output = "shortestDistance",
                                    FeatureLocForDistance = bindingType, 
                                    PeakLocForDistance = "middle")
    })
    
    anno <- unique(anno)
    if (length(anno) < 1L) {
        stop("Cannot annotate input 'x' with the AnnotationData", 
             call. = FALSE)
    }
    distancetoFeature <- anno$distancetoFeature
    # Calculate the breaks
    # Use boxplot.stats to remove outliers
    box_stats <- boxplot.stats(distancetoFeature, coef = 1.5, 
                               do.conf = FALSE, do.out = FALSE)
    whiskers <- box_stats$stats[5L] - box_stats$stats[1L]
    minWid <- 100L
    if (whiskers > 100L * minWid) {
        breaks <- seq(box_stats$stats[1L], box_stats$stats[5L],
                      length.out = 100L)
        minWid <- diff(breaks)[1L]
    } else {
        if (whiskers > minWid) {
            breaks <- seq(floor(box_stats$stats[1L] / minWid), 
                          ceiling(box_stats$stats[5L] / minWid),
                          by = 1L) * minWid
        } else {
            breaks <- mean(box_stats$stats[c(1L, 5L)])
            breaks <- c(breaks - minWid / 2L, breaks + minWid / 2L)
        }
    }
    if (min(distancetoFeature) < min(breaks)) {
        breaks <- c(min(distancetoFeature), breaks)
    }
    if (max(distancetoFeature) > max(breaks)) {
        breaks <- c(breaks, max(distancetoFeature))
    }
    h <- hist(distancetoFeature, breaks = breaks, plot = FALSE)
    N <- as.integer(h$counts)
    diff <- as.integer(floor(minWid / 2L))
    offset <- as.integer(floor(h$mids))
    offset <- offset[N > 0L]
    N <- N[N > 0L]
    new("bindist", counts = N,
        mids = offset,
        halfBinSize = diff,
        bindingType = bindingType,
        featureType = featureType)
}