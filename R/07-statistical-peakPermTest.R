#' Permutation test for two peak lists
#' 
#' @description 
#' Performs a permutation test to assess whether there is a significant
#' association between two peak lists. The test generates random background
#' peaks that match the observed binding distribution of \code{peaks1} (e.g.,
#' distance from TSS) and compares the overlap between \code{peaks1} and
#' \code{peaks2} to the distribution of overlaps with random peaks. This approach
#' accounts for genomic biases in peak distribution, making the test more
#' appropriate than simple random sampling.
#' 
#' The test is one-sided (alternative = "greater"), testing whether the observed
#' overlap between \code{peaks1} and \code{peaks2} is significantly greater than
#' expected by chance, given the binding distribution of \code{peaks1}.
#' 
#' @param peaks1 A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        representing the first peak list. This is the set of peaks that will
#'        be randomized (replaced with random peaks matching their binding
#'        distribution) in each permutation.
#' @param peaks2 A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        representing the second peak list. This is the set of peaks that
#'        overlaps are counted against in each permutation.
#' @param ntimes An integer specifying the number of permutations to perform.
#'        Default is 100. More permutations provide more accurate p-values but
#'        take longer to compute.
#' @param seed An integer specifying the random seed for reproducibility.
#'        Default is \code{as.integer(Sys.time())}, which uses the current
#'        system time. Set a fixed seed for reproducible results.
#' @param mc.cores An integer specifying the number of cores to use for parallel
#'        computation. Default is \code{getOption("mc.cores", 2L)}. This parameter
#'        is passed directly to \code{\link[regioneR]{permTest}}, which uses
#'        \code{\link[parallel]{mclapply}} for parallelization. Set to 1 to disable
#'        parallelization. Note: This function uses \code{mc.cores} instead of
#'        \code{BPPARAM} for direct compatibility with the underlying
#'        \code{regioneR::permTest} function.
#' @param maxgap An integer specifying the maximum gap (in base pairs) allowed
#'        between ranges for them to be considered overlapping. Default is
#'        \code{-1L} (ranges must directly overlap with no gap). See
#'        \code{\link[IRanges]{findOverlaps}} for details. Note: In this function,
#'        \code{maxgap} is passed to \code{\link{cntOverlaps}}, which expands
#'        ranges in \code{peaks1} by \code{maxgap} before counting overlaps.
#' @param pool Optional. An object of class \code{\link{permPool}} containing
#'        pre-computed sampling regions. If not provided, will be generated from
#'        \code{TxDb} and \code{peaks1} using \code{\link{preparePool}}. Providing
#'        a pre-computed pool can save time when running multiple tests with the
#'        same binding distribution.
#' @param TxDb Optional. An object of class
#'        \code{\link[GenomicFeatures:TxDb-class]{TxDb}} containing annotation
#'        data. Required if \code{pool} is not provided. Used to determine the
#'        binding distribution and create the sampling pool.
#' @param bindingDistribution Optional. An object of class \code{\link{bindist}}
#'        containing the observed binding distribution. If not provided, will be
#'        computed from \code{peaks1} using \code{\link{buildBindingDistribution}}.
#'        Providing a pre-computed distribution can save time when running
#'        multiple tests.
#' @param bindingType A character string specifying the feature position to use
#'        for distance calculation when building the binding distribution.
#'        Options:
#'        \itemize{
#'          \item \code{"TSS"} (default): Transcription start site
#'          \item \code{"geneEnd"}: Transcription end site
#'        }
#'        Only used if \code{bindingDistribution} is not provided.
#' @param featureType A character string specifying the type of annotation
#'        feature to use when building the binding distribution. Options:
#'        \itemize{
#'          \item \code{"transcript"} (default): Use transcript-level features
#'          \item \code{"exon"}: Use exon-level features
#'        }
#'        Only used if \code{bindingDistribution} is not provided.
#' @param seqn Optional character vector. If provided, restricts the sampling
#'        pool to the specified chromosome names (seqnames). Default is
#'        \code{NA} (no restriction, uses all chromosomes in the annotation).
#'        Useful for testing associations on specific chromosomes only.
#' @param \dots Further arguments passed to \code{\link[regioneR]{numOverlaps}}
#'        for controlling overlap detection (e.g., \code{count.once}).
#' 
#' @return Returns an object of class \code{permTestResults} from the
#'        \code{regioneR} package. This object contains:
#'        \itemize{
#'          \item \code{observed}: The observed number of overlaps between
#'                \code{peaks1} and \code{peaks2}
#'          \item \code{randomized}: A numeric vector of overlap counts from
#'                \code{ntimes} permutations
#'          \item \code{pval}: The p-value (proportion of random overlaps >=
#'                observed overlap)
#'          \item \code{zscore}: The z-score (standardized difference from mean
#'                of random distribution)
#'          \item \code{alternative}: The alternative hypothesis tested
#'          \item \code{ntimes}: The number of permutations performed
#'        }
#'        See \code{\link[regioneR]{permTest}} for complete details on the
#'        return object structure.
#' 
#' @details
#' 
#' \strong{How the permutation test works:}
#' \enumerate{
#'   \item Builds or uses the binding distribution of \code{peaks1} (distance
#'         from TSS/geneEnd)
#'   \item Creates a sampling pool using \code{\link{preparePool}} that matches
#'         this binding distribution
#'   \item For each permutation:
#'         \itemize{
#'           \item Generates random peaks using \code{\link{randPeaks}} that:
#'                 \enumerate{
#'                   \item Sample positions from genomic regions grouped by
#'                         distance bins (matching the binding distribution)
#'                   \item Assign random widths from an exponential distribution
#'                         matching the width distribution of \code{peaks1}
#'                 }
#'           \item Counts overlaps between random peaks and \code{peaks2} using
#'                 \code{\link{cntOverlaps}}
#'         }
#'   \item Compares the observed overlap count to the distribution of random
#'         overlap counts
#'   \item Calculates p-value as the proportion of random overlaps >= observed
#'         overlap
#' }
#' 
#' \strong{Why match the binding distribution?}
#' ChIP-seq peaks are not randomly distributed in the genome. They tend to
#' cluster near regulatory elements like promoters (TSS). If random peaks were
#' sampled uniformly from the genome, they would have a different distance
#' distribution than observed peaks, making the test less appropriate. By
#' matching the binding distribution, random peaks have the same genomic biases
#' as observed peaks, making the test more powerful and appropriate.
#' 
#' \strong{Random peak generation:}
#' Random peaks are generated by \code{\link{randPeaks}}, which:
#' \itemize{
#'   \item Samples positions uniformly from genomic regions in each distance
#'         bin (from \code{pool$grs})
#'   \item Samples the number of peaks from each bin according to
#'         \code{pool$N} (matching the binding distribution)
#'   \item Assigns widths from an exponential distribution with mean adjusted
#'         to match the width distribution of \code{peaks1}
#'   \item Trims peaks to valid genomic coordinates
#' }
#' 
#' \strong{Overlap counting:}
#' Overlaps are counted using \code{\link{cntOverlaps}}, which:
#' \itemize{
#'   \item If \code{maxgap > 0}: Expands ranges in \code{peaks1} (or random
#'         peaks) by \code{maxgap} bp on both sides before counting
#'   \item Counts how many ranges in \code{peaks1} (or random peaks) overlap
#'         with at least one range in \code{peaks2}
#'   \item Each range is counted at most once, regardless of how many ranges
#'         in \code{peaks2} it overlaps with
#' }
#' 
#' \strong{Interpretation:}
#' \itemize{
#'   \item \code{pval < 0.05}: Significant association - \code{peaks1} and
#'         \code{peaks2} overlap more than expected by chance
#'   \item \code{pval >= 0.05}: No significant association - overlap could be
#'         due to chance
#'   \item \code{zscore > 0}: Observed overlap is above the mean of random
#'         distribution
#'   \item \code{zscore < 0}: Observed overlap is below the mean of random
#'         distribution
#' }
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{preparePool}} for creating sampling pools,
#'          \code{\link{buildBindingDistribution}} for building binding
#'          distributions, \code{\link{bindist}} for the distribution class,
#'          \code{\link{randPeaks}} for random peak generation,
#'          \code{\link{cntOverlaps}} for overlap counting,
#'          \code{\link[regioneR]{permTest}} for the underlying permutation
#'          test framework
#' @references Davison, A. C. and Hinkley, D. V. (1997) Bootstrap methods and
#' their application, Cambridge University Press, United Kingdom, 156-160
#' @keywords misc
#' @export
#' @importFrom regioneR permTest numOverlaps
#' @examples
#' \dontrun{
#' ## Example 1: Basic permutation test
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' data("myPeakList")
#' peaks2 <- GRanges("chr1", IRanges(c(1000000, 2000000), width = 500))
#' result <- peakPermTest(myPeakList, peaks2, 
#'                        TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
#' result$pval  # P-value
#' result$zscore  # Z-score
#' 
#' ## Example 2: Using pre-computed binding distribution
#' dist <- buildBindingDistribution(myPeakList, 
#'                                  TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                  bindingType = "TSS")
#' result <- peakPermTest(myPeakList, peaks2,
#'                        TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                        bindingDistribution = dist)
#' 
#' ## Example 3: Allow gaps in overlap detection
#' result <- peakPermTest(myPeakList, peaks2,
#'                        TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                        maxgap = 1000L)  # Allow 1kb gaps
#' 
#' ## Example 4: Test on specific chromosomes only
#' result <- peakPermTest(myPeakList, peaks2,
#'                        TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                        seqn = c("chr1", "chr2"))
#' 
#' ## Example 5: More permutations for higher accuracy
#' result <- peakPermTest(myPeakList, peaks2,
#'                        TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                        ntimes = 1000L, seed = 12345)
#' }
#' 
peakPermTest <- function(peaks1, peaks2, ntimes = 100L, 
                         seed = as.integer(Sys.time()),
                         mc.cores = getOption("mc.cores", 2L),
                         maxgap = -1L, pool,
                         TxDb, bindingDistribution,
                         bindingType = c("TSS", "geneEnd"), 
                         featureType = c("transcript", "exon"),
                         seqn = NA, ...) {
    if (!inherits(peaks1, "GRanges") || !inherits(peaks2, "GRanges")) {
        stop("'peaks1' and 'peaks2' must be GRanges objects", call. = FALSE)
    }
    if (!missing(bindingDistribution)) {
        if (!inherits(bindingDistribution, "bindist")) {
            stop("'bindingDistribution' must be a bindist object", 
                 call. = FALSE)
        }
    }
    
    bindingType <- match.arg(bindingType)
    featureType <- match.arg(featureType)
    
    if (missing(pool)) {
        if (missing(TxDb)) {
            stop("Either 'pool' or 'TxDb' must be provided", call. = FALSE)
        }
        pool <- preparePool(TxDb, template = peaks1, 
                            bindingDistribution = bindingDistribution, 
                            bindingType = bindingType, 
                            featureType = featureType, seqn = seqn)
    } else {
        if (!inherits(pool, "permPool")) {
            stop("'pool' must be a permPool object", call. = FALSE)
        }
    }
    
    set.seed(seed)
    pt <- permTest(A = peaks1, ntimes = ntimes,
                   randomize.function = randPeaks, 
                   grs = pool$grs, N = pool$N, 
                   maxgap = maxgap,
                   evaluate.function = cntOverlaps, 
                   B = peaks2, verbose = FALSE, alternative = "greater", 
                   mc.set.seed = FALSE, mc.cores = mc.cores,
                   ...)
    pt
}
