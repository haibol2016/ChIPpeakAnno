#' Filter peaks by IDR (irreproducible discovery rate)
#' 
#' @description 
#' Uses the IDR (irreproducible discovery rate) method to assess consistency
#' between replicate ChIP-seq experiments and filter peaks to obtain a
#' high-confidence set. IDR measures the reproducibility of peak rankings
#' between replicates based on signal strength (read coverage density).
#' 
#' The function identifies overlapping peaks between two replicates, calculates
#' normalized read coverage for each overlapping peak, and uses the IDR
#' statistical framework to estimate the probability that a peak is
#' irreproducible. Peaks with IDR values below the threshold are retained as
#' high-confidence peaks.
#' 
#' IDR is particularly useful for ChIP-seq quality control, as it provides a
#' principled way to combine information from replicates and identify peaks
#' that are consistently detected across replicates.
#' 
#' @param peaksA,peaksB \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#'        objects containing peaks from two replicate ChIP-seq experiments.
#'        Peaks should be called independently for each replicate (e.g., using
#'        MACS, SPP, or other peak callers). The function will identify
#'        overlapping peaks between the two sets.
#' @param bamfileA,bamfileB Character strings specifying file paths to BAM
#'        files containing aligned reads for the corresponding replicates. The
#'        BAM files must be indexed (have corresponding .bai files). Read
#'        coverage from these files is used to calculate signal strength for
#'        IDR estimation.
#' @param maxgap An integer specifying the maximum gap (in base pairs) allowed
#'        between ranges for them to be considered overlapping. Default is
#'        \code{-1L} (ranges must actually overlap with no gap allowed). See
#'        \code{\link[IRanges]{findOverlaps}} for details.
#' @param minoverlap An integer specifying the minimum overlap required (in
#'        base pairs). Default is \code{0L} (any overlap is considered). See
#'        \code{\link[IRanges]{findOverlaps}} for details.
#' @param singleEnd A logical value. If \code{TRUE} (default), reads are
#'        treated as single-end. If \code{FALSE}, reads are treated as
#'        paired-end. This affects how coverage is calculated by
#'        \code{summarizeOverlaps}.
#' @param IDRcutoff A numeric value between 0 and 1 specifying the IDR
#'        threshold. Peaks with IDR >= \code{IDRcutoff} are removed (filtered
#'        out). Default is \code{0.01} (1\% IDR threshold), meaning peaks with
#'        IDR < 0.01 are retained. Lower values are more stringent (retain
#'        fewer peaks), higher values are more lenient (retain more peaks).
#'        Common thresholds: 0.01 (1\%), 0.05 (5\%), 0.1 (10\%).
#' @param ... Additional arguments (currently not used).
#' 
#' @return Returns a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        containing filtered peaks that pass the IDR threshold. The returned
#'        object:
#'        \itemize{
#'          \item Contains only overlapping peaks between \code{peaksA} and
#'                \code{peaksB} that have IDR < \code{IDRcutoff}
#'          \item Uses merged peak regions (from \code{findOverlapsOfPeaks})
#'          \item Has names in the format "olp0001", "olp0002", etc.
#'          \item Returns an empty GRanges object if no overlapping peaks are
#'                found or if all peaks fail the IDR threshold
#'        }
#'        Note: The function only processes overlapping peaks. Peaks that don't
#'        overlap between replicates are excluded from the analysis.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Finds overlapping peaks between \code{peaksA} and \code{peaksB} using
#'         \code{\link{findOverlapsOfPeaks}} with the specified \code{maxgap}
#'         and \code{minoverlap} parameters
#'   \item Extracts the merged overlapping peaks (peaks present in both
#'         replicates)
#'   \item If no overlapping peaks are found, returns an empty GRanges object
#'   \item Calculates read coverage for each overlapping peak from both BAM
#'         files using \code{\link[GenomicAlignments]{summarizeOverlaps}} with
#'         \code{mode = Union}
#'   \item Normalizes coverage by peak width to obtain signal density (reads
#'         per base pair) for each replicate
#'   \item Estimates IDR values using \code{idr::est.IDR} with default
#'         parameters (mu = 2.07, sigma = 1.34, rho = 0.89, p = 0.84)
#'   \item Filters peaks: retains only peaks with IDR < \code{IDRcutoff}
#' }
#' 
#' \strong{IDR estimation:}
#' The function uses the IDR package's \code{est.IDR} function with default
#' parameters optimized for ChIP-seq data:
#' \itemize{
#'   \item \code{mu = 2.07}: Mean parameter for the bivariate normal
#'         distribution
#'   \item \code{sigma = 1.34}: Standard deviation parameter
#'   \item \code{rho = 0.89}: Correlation parameter between replicates
#'   \item \code{p = 0.84}: Mixing proportion parameter
#' }
#' These parameters are based on typical ChIP-seq replicate behavior and may
#' need adjustment for specific experimental conditions.
#' 
#' \strong{Signal density calculation:}
#' Signal strength is calculated as normalized read coverage:
#' \code{coverage / peak_width}
#' This gives reads per base pair, which accounts for differences in peak
#' sizes. The IDR framework then compares these normalized signal values
#' between replicates to assess reproducibility.
#' 
#' \strong{Overlap detection:}
#' Only peaks that overlap between the two replicates are analyzed. The
#' overlap detection uses \code{findOverlapsOfPeaks}, which:
#' \itemize{
#'   \item Identifies peaks that overlap based on \code{maxgap} and
#'         \code{minoverlap}
#'   \item Merges overlapping peaks into a single region
#'   \item Returns merged regions in the \code{peaklist} with names containing
#'         "///" (indicating overlap between multiple peak sets)
#' }
#' 
#' \strong{Requirements:}
#' \itemize{
#'   \item The \code{idr} package must be installed (available from
#'         Bioconductor)
#'   \item The \code{DelayedArray} package must be installed (available from
#'         Bioconductor)
#'   \item BAM files must be indexed (have .bai files)
#'   \item BAM files must be readable and contain aligned reads
#' }
#' 
#' \strong{Limitations:}
#' \itemize{
#'   \item Only processes overlapping peaks; non-overlapping peaks are excluded
#'   \item Uses fixed IDR parameters; may need adjustment for specific
#'         experimental conditions
#'   \item Requires both BAM files to be available (not just peak files)
#'   \item Can be computationally intensive for large peak sets
#' }
#' 
#' @author Jianhong Ou
#' @references Li, Qunhua, et al. "Measuring reproducibility of high-throughput
#' experiments." The annals of applied statistics (2011): 1752-1779.
#' @seealso \code{\link{findOverlapsOfPeaks}} for overlap detection,
#'          \code{\link[GenomicAlignments]{summarizeOverlaps}} for coverage
#'          calculation, \code{\link[idr]{est.IDR}} for IDR estimation
#' @keywords misc
#' @export
#' @importFrom GenomicAlignments summarizeOverlaps Union
#' @importFrom SummarizedExperiment assay
#' @importFrom idr est.IDR
#' @importFrom DelayedArray rowRanges
#' @examples
#' \dontrun{
#' ## Example 1: Basic IDR filtering
#' library(idr)
#' library(DelayedArray)
#' peaksA <- toGRanges("replicate1_peaks.bed", format = "BED")
#' peaksB <- toGRanges("replicate2_peaks.bed", format = "BED")
#' filtered_peaks <- IDRfilter(peaksA, peaksB,
#'                              bamfileA = "replicate1.bam",
#'                              bamfileB = "replicate2.bam")
#' 
#' ## Example 2: More stringent IDR threshold (0.5%)
#' filtered_peaks <- IDRfilter(peaksA, peaksB,
#'                              bamfileA = "replicate1.bam",
#'                              bamfileB = "replicate2.bam",
#'                              IDRcutoff = 0.005)
#' 
#' ## Example 3: Allow gaps in overlap detection
#' filtered_peaks <- IDRfilter(peaksA, peaksB,
#'                              bamfileA = "replicate1.bam",
#'                              bamfileB = "replicate2.bam",
#'                              maxgap = 100L)  # Allow 100bp gaps
#' 
#' ## Example 4: Paired-end reads
#' filtered_peaks <- IDRfilter(peaksA, peaksB,
#'                              bamfileA = "replicate1.bam",
#'                              bamfileB = "replicate2.bam",
#'                              singleEnd = FALSE)
#' 
#' ## Example 5: Using MACS output files
#' peaksA <- toGRanges("MACS_rep1_peaks.xls", format = "MACS2")
#' peaksB <- toGRanges("MACS_rep2_peaks.xls", format = "MACS2")
#' filtered_peaks <- IDRfilter(peaksA, peaksB,
#'                              bamfileA = "rep1.bam",
#'                              bamfileB = "rep2.bam")
#' }
#' 
IDRfilter <- function(peaksA, peaksB, bamfileA, bamfileB, 
                      maxgap = -1L, minoverlap = 0L, singleEnd = TRUE,
                      IDRcutoff = 0.01, ...) {
    if (!requireNamespace("idr", quietly = TRUE)) {
        stop("The 'idr' package is required. Install it with: ",
             "BiocManager::install('idr')", call. = FALSE)
    }
    if (!requireNamespace("DelayedArray", quietly = TRUE)) {
        stop("The 'DelayedArray' package is required. Install it with: ",
             "BiocManager::install('DelayedArray')", call. = FALSE)
    }
    
    stopifnot(
        inherits(peaksA, "GRanges"),
        inherits(peaksB, "GRanges"),
        file.exists(bamfileA),
        file.exists(bamfileB)
    )
    
    ol <- findOverlapsOfPeaks(peaksA,
                              peaksB,
                              maxgap = maxgap,
                              minoverlap = minoverlap)
    ol <- ol$peaklist[grepl("\\/\\/\\/", names(ol$peaklist))][[1L]]
    
    if (length(ol) < 1L) {
        return(GRanges())
    }
    
    names(ol) <- paste0("olp",
                        formatC(seq_along(ol),
                                width = nchar(as.character(length(ol))),
                                flag = "0"))
    
    coverage <- summarizeOverlaps(features = ol,
                                  reads = c(bamfileA, bamfileB),
                                  mode = Union,
                                  ignore.strand = FALSE,
                                  singleEnd = singleEnd)
    
    idr <- idr::est.IDR(
        assay(coverage) / width(DelayedArray::rowRanges(coverage)),
        mu = 2.07,
        sigma = 1.34,
        rho = 0.89,
        p = 0.84
    )
    
    ol[idr$IDR < IDRcutoff]
}
