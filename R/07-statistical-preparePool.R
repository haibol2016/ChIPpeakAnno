#' Prepare data for permutation test
#' 
#' @description 
#' Prepares a pool of genomic regions for permutation testing. The function
#' generates distance-binned regions around annotation features (e.g., TSS or
#' gene ends) based on the observed binding distribution of input peaks. This
#' pool is used by \code{\link{peakPermTest}} to generate random background
#' regions that match the observed binding pattern, ensuring that permutation
#' tests account for genomic biases in peak distribution.
#' 
#' The function creates regions at specific distances from annotation features,
#' grouped into bins matching the binding distribution. Regions that overlap
#' with original annotation features are filtered out to avoid sampling from
#' annotated regions. The resulting pool can be used to generate random peaks
#' that have the same distance distribution as observed peaks.
#' 
#' @param TxDb An object of class
#'        \code{\link[GenomicFeatures:TxDb-class]{TxDb}} containing annotation
#'        data. Required. Used to extract transcripts or exons for building the
#'        sampling pool.
#' @param template A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        representing the peaks to be tested. Required if
#'        \code{bindingDistribution} is not provided. Used to build the binding
#'        distribution by calculating distances from peaks to annotation features.
#' @param bindingDistribution Optional. An object of class \code{\link{bindist}}
#'        containing the observed binding distribution. If not provided, will be
#'        computed from \code{template} using \code{\link{buildBindingDistribution}}.
#'        Providing a pre-computed distribution can save time when running
#'        multiple tests with the same binding pattern.
#' @param bindingType A character string specifying the feature position to use
#'        for distance calculation when building the binding distribution.
#'        Options:
#'        \itemize{
#'          \item \code{"TSS"} (default): Transcription start site
#'          \item \code{"geneEnd"}: Transcription end site
#'        }
#'        Only used if \code{bindingDistribution} is not provided.
#' @param featureType A character string specifying the type of annotation
#'        feature to use. Options:
#'        \itemize{
#'          \item \code{"transcript"} (default): Use transcript-level features
#'          \item \code{"exon"}: Use exon-level features
#'        }
#'        This determines which features are extracted from \code{TxDb} and used
#'        to build the sampling pool.
#' @param seqn Optional character vector. If provided, restricts the pool to
#'        the specified chromosome names (seqnames). Default is \code{NA} (no
#'        restriction, uses all chromosomes in the annotation). Useful for
#'        testing associations on specific chromosomes only.
#' 
#' @return Returns a list with two elements that can be used to create a
#'        \code{\link{permPool}} object:
#'        \itemize{
#'          \item \code{grs}: A \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}
#'                object containing distance-binned genomic regions for sampling.
#'                Each element of the list corresponds to one distance bin from
#'                the binding distribution. Regions in each bin are positioned at
#'                the appropriate distance from annotation features and filtered
#'                to exclude overlaps with original annotation features.
#'          \item \code{N}: An integer vector specifying the number of regions
#'                to sample from each corresponding element in \code{grs}. The
#'                length of \code{N} equals the length of \code{grs}. Values in
#'                \code{N} correspond to the number of peaks observed in each
#'                distance bin (from the binding distribution), ensuring that
#'                random peaks match the observed distance distribution.
#'        }
#'        This list can be converted to a \code{permPool} object using
#'        \code{new("permPool", grs = result$grs, N = result$N)} or passed
#'        directly to \code{peakPermTest} via the \code{pool} parameter.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Extracts annotation features (transcripts or exons) from \code{TxDb}
#'         and removes duplicates
#'   \item If \code{bindingDistribution} is missing, builds it from
#'         \code{template} using \code{buildBindingDistribution}
#'   \item Filters annotation features by chromosome if \code{seqn} is provided
#'   \item For each distance bin in the binding distribution:
#'         \enumerate{
#'           \item Shifts annotation features by the bin midpoint distance
#'                 (\code{offset}) from their reference point (TSS or geneEnd)
#'           \item Expands regions by \code{halfBinSize} on both sides to create
#'                 the full bin width (\code{2 * halfBinSize})
#'           \item Trims regions to valid genomic coordinates
#'           \item Removes regions with width <= 1 bp
#'           \item Filters out regions that overlap with original annotation
#'                 features:
#'                 \itemize{
#'                   \item For upstream regions (negative offset): Removes all
#'                         overlapping regions
#'                   \item For downstream regions (positive offset): Removes
#'                         regions that overlap with different annotation features
#'                         (keeps overlaps with the same feature, as these are
#'                         expected)
#'                 }
#'         }
#'   \item Returns the list of distance-binned regions and sampling counts
#' }
#' 
#' \strong{Region creation process:}
#' For a distance bin centered at, for example, -5000 bp from TSS:
#' \itemize{
#'   \item Annotation features are shifted by -5000 bp (upstream)
#'   \item Regions are expanded by \code{halfBinSize} on both sides
#'   \item If \code{halfBinSize = 50}, each region spans [-5050, -4950] relative
#'         to TSS
#'   \item Regions are trimmed to valid coordinates
#'   \item Regions overlapping with original annotation features are removed
#' }
#' 
#' \strong{Overlap filtering:}
#' The function filters out regions that overlap with original annotation
#' features to avoid sampling from annotated regions. The filtering logic differs
#' for upstream vs. downstream regions:
#' \itemize{
#'   \item \strong{Upstream regions} (negative offset): All overlapping regions
#'         are removed, as these would be within or very close to annotated
#'         regions
#'   \item \strong{Downstream regions} (positive offset): Only regions
#'         overlapping with \emph{different} annotation features are removed.
#'         Overlaps with the same feature are kept, as these represent valid
#'         downstream positions
#' }
#' 
#' \strong{Use in permutation testing:}
#' The pool created by this function is used by \code{\link{randPeaks}} to
#' generate random peaks:
#' \itemize{
#'   \item For each distance bin i, \code{randPeaks} samples \code{N[i]}
#'         positions uniformly from \code{grs[[i]]}
#'   \item Random widths are assigned matching the width distribution of
#'         observed peaks
#'   \item The resulting random peaks have the same distance distribution as
#'         observed peaks
#' }
#' 
#' \strong{Performance considerations:}
#' \itemize{
#'   \item Building the pool can be time-consuming for large genomes or many
#'         distance bins
#'   \item Consider pre-computing and reusing the pool when running multiple
#'         tests with the same binding distribution
#'   \item Filtering by \code{seqn} can reduce computation time if only
#'         specific chromosomes are needed
#' }
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{peakPermTest}} for using the pool in permutation tests,
#'          \code{\link{buildBindingDistribution}} for building binding
#'          distributions, \code{\link{bindist}} for the distribution class,
#'          \code{\link{permPool}} for the pool class,
#'          \code{\link{randPeaks}} for random peak generation
#' @keywords misc
#' @export
#' @importFrom S4Vectors queryHits subjectHits
#' @examples
#' \dontrun{
#' ## Example 1: Basic usage
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' data("myPeakList")
#' pool <- preparePool(TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                     template = myPeakList,
#'                     bindingType = "TSS",
#'                     featureType = "transcript")
#' 
#' ## Access the pool components
#' pool$grs  # GRangesList of distance-binned regions
#' pool$N    # Number of regions to sample from each bin
#' 
#' ## Example 2: Using pre-computed binding distribution
#' dist <- buildBindingDistribution(myPeakList,
#'                                   TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                   bindingType = "TSS")
#' pool <- preparePool(TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                     bindingDistribution = dist,
#'                     featureType = "transcript")
#' 
#' ## Example 3: Restrict to specific chromosomes
#' pool <- preparePool(TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                     template = myPeakList,
#'                     seqn = c("chr1", "chr2"))
#' 
#' ## Example 4: Use with exon-level features
#' pool <- preparePool(TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                     template = myPeakList,
#'                     featureType = "exon",
#'                     bindingType = "TSS")
#' 
#' ## Example 5: Use in permutation test
#' peaks2 <- GRanges("chr1", IRanges(c(1000000, 2000000), width = 500))
#' result <- peakPermTest(myPeakList, peaks2, pool = pool)
#' }
#' 
preparePool <- function(TxDb, template, bindingDistribution,
                        bindingType = c("TSS", "geneEnd"), 
                        featureType = c("transcript", "exon"), 
                        seqn = NA) {
    if (missing(TxDb)) {
        stop("'TxDb' is required", call. = FALSE)
    }
    if (!inherits(TxDb, "TxDb")) {
        stop("'TxDb' must be a TxDb object", call. = FALSE)
    }
    if (!missing(bindingDistribution)) {
        if (!inherits(bindingDistribution, "bindist")) {
            stop("'bindingDistribution' must be a bindist object", 
                 call. = FALSE)
        }
    }
    
    bindingType <- match.arg(bindingType)
    featureType <- match.arg(featureType)
    
    if (featureType == "transcript") {
        tx <- transcripts(x = TxDb)
    } else {
        tx <- exons(x = TxDb)
    }
    tx <- unique(tx)
    names(tx) <- seq_along(tx)
    
    if (missing(bindingDistribution)) {
        bindingDistribution <- 
            buildBindingDistribution(template, AnnotationData = tx,
                                     bindingType = bindingType,
                                     featureType = featureType)
    }
    
    offset <- bindingDistribution$mids
    diff <- bindingDistribution$halfBinSize
    N <- bindingDistribution$counts
    
    if (!is.na(seqn[1L])) {
        tx <- tx[as.character(seqnames(tx)) %in% seqn]
    }
    
    grs <- lapply(offset, function(.off) {
        x <- NULL
        suppressWarnings({
            x <- shift(tx, .off)
            start(x) <- start(x) - diff
            width(x) <- 2L * diff
            x <- trim(x)
            x <- x[width(x) > 1L]
        })
        ol <- findOverlaps(x, tx)
        if (length(ol) > 0L) {
            if (.off + diff < 0L) {
                x <- x[-(unique(queryHits(ol)))]
            } else {
                x <- x[-(unique(queryHits(ol)[
                    names(x)[queryHits(ol)] != names(tx)[subjectHits(ol)]]))]
            }
        }
        x
    })
    
    list(grs = GRangesList(grs), N = as.integer(N))
}
