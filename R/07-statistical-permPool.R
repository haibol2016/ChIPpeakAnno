#' Class \code{"permPool"}
#' 
#' @description 
#' An S4 class representing a pool of genomic regions for permutation testing.
#' Objects of this class store genomic regions grouped by distance bins (relative
#' to annotation features like TSS) and specify how many regions to sample from
#' each bin. This enables generation of random background regions that match the
#' observed binding distribution of peaks, making permutation tests more
#' appropriate by accounting for genomic biases.
#' 
#' The class is typically created by \code{\link{preparePool}} and used by
#' \code{\link{peakPermTest}} to generate random peaks for permutation testing.
#' Each element in \code{grs} corresponds to a distance bin, and the
#' corresponding element in \code{N} specifies how many regions to sample from
#' that bin to match the observed binding distribution.
#' 
#' @name permPool-class
#' @rdname permPool
#' @aliases permPool permPool-class permPool-method $,permPool-method
#' $<-,permPool-method
#' @docType class
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form:
#' \code{new("permPool", grs = GRangesList(...), N = integer(...))}
#' 
#' However, it is recommended to use \code{\link{preparePool}} to create
#' \code{permPool} objects, as it handles the complex logic of building
#' distance-binned regions from annotation data and binding distributions.
#' 
#' @slot grs An object of class \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}
#'        containing genomic regions grouped by distance bins. Each element of
#'        the list corresponds to one distance bin (e.g., -5000 to -4900 bp from
#'        TSS). Regions in each bin are shifted versions of annotation features
#'        (transcripts or exons) positioned at the appropriate distance from
#'        their reference point (TSS or geneEnd). Regions that overlap with the
#'        original annotation features are filtered out to avoid sampling from
#'        annotated regions.
#' @slot N An integer vector specifying the number of regions to sample from
#'        each corresponding element in \code{grs}. The length of \code{N} must
#'        equal the length of \code{grs}. Values in \code{N} correspond to the
#'        number of peaks observed in each distance bin (from the binding
#'        distribution), ensuring that random peaks match the observed distance
#'        distribution.
#' 
#' @section Validity:
#' The \code{permPool} class has the following validity check:
#' \itemize{
#'   \item The length of \code{grs} must equal the length of \code{N}
#' }
#' 
#' @section Methods:
#' \describe{
#'   \item{\code{x$grs}, \code{x$N}}{Access slots using the \code{$} operator}
#'   \item{\code{x$grs <- value}, \code{x$N <- value}}{Replace slot values using
#'         the \code{$<-} operator}
#' }
#' 
#' @details
#' 
#' \strong{Purpose:}
#' The \code{permPool} class enables matched permutation testing by providing
#' a structured way to sample random genomic regions that match the distance
#' distribution of observed peaks. This is crucial because ChIP-seq peaks are
#' not randomly distributed - they cluster near regulatory elements like
#' promoters. Random regions must match this distribution to make statistical
#' tests appropriate.
#' 
#' \strong{Structure:}
#' \itemize{
#'   \item \code{grs} is a \code{GRangesList} where each element contains
#'         regions at a specific distance from annotation features
#'   \item \code{N} is an integer vector where each element specifies how many
#'         regions to sample from the corresponding element in \code{grs}
#'   \item The i-th element of \code{grs} and \code{N} work together: sample
#'         \code{N[i]} regions from \code{grs[[i]]}
#' }
#' 
#' \strong{How regions are created:}
#' For each distance bin (e.g., centered at -5000 bp from TSS):
#' \enumerate{
#'   \item Annotation features (transcripts/exons) are shifted by the bin
#'         midpoint distance
#'   \item Regions are expanded by \code{halfBinSize} on both sides to create
#'         the bin width
#'   \item Regions are trimmed to valid genomic coordinates
#'   \item Regions that overlap with original annotation features are removed
#'         (to avoid sampling from annotated regions)
#'   \item The resulting regions are stored in \code{grs[[i]]}
#' }
#' 
#' \strong{Use in permutation testing:}
#' When generating random peaks:
#' \enumerate{
#'   \item For each distance bin i, sample \code{N[i]} positions uniformly from
#'         \code{grs[[i]]}
#'   \item Assign random widths matching the width distribution of observed peaks
#'   \item The resulting random peaks have the same distance distribution as
#'         observed peaks
#' }
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{preparePool}} for creating \code{permPool} objects,
#'          \code{\link{peakPermTest}} for using \code{permPool} in permutation
#'          tests, \code{\link{bindist}} for the binding distribution class,
#'          \code{\link{randPeaks}} for random peak generation
#' @keywords classes
#' @exportClass permPool
#' @exportMethod "$" "$<-"
#' @examples
#' \dontrun{
#' ## Example 1: Create permPool using preparePool (recommended)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' data("myPeakList")
#' pool <- preparePool(TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                     template = myPeakList,
#'                     bindingType = "TSS",
#'                     featureType = "transcript")
#' 
#' ## Access slots
#' pool$grs  # GRangesList of distance-binned regions
#' pool$N    # Number of regions to sample from each bin
#' 
#' ## Example 2: Manual creation (not recommended)
#' grs <- GRangesList(
#'     GRanges("chr1", IRanges(1000, 2000)),
#'     GRanges("chr1", IRanges(3000, 4000))
#' )
#' N <- c(10L, 20L)
#' pool <- new("permPool", grs = grs, N = N)
#' 
#' ## Example 3: Use in permutation test
#' peaks2 <- GRanges("chr1", IRanges(c(1500, 3500), width = 500))
#' result <- peakPermTest(myPeakList, peaks2, pool = pool)
#' }

setClass("permPool", representation(grs="GRangesList",
                                    N="integer"),
         validity=function(object){
             re <- TRUE
             if(length(object@grs) != length(object@N)) {
                 re <- "the length of grs and N are not identical"
             }
             re
         })

setMethod("$", "permPool", function(x, name) slot(x, name))
setReplaceMethod("$", "permPool",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })