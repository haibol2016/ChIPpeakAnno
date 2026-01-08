#' Estimate fragment length from BAM files
#' 
#' @description 
#' Estimates the average fragment length (insert size) for ChIP-seq or other
#' DNA-seq experiments from BAM files. The function automatically detects
#' whether the data is paired-end or single-end and uses the appropriate
#' method:
#' \itemize{
#'   \item \strong{Paired-end}: Analyzes insert size distribution from proper
#'         pairs and finds the peak
#'   \item \strong{Single-end}: Performs cross-correlation analysis between
#'         plus-strand and minus-strand read coverage to find the optimal shift
#' }
#' 
#' The estimated fragment length is essential for:
#' \itemize{
#'   \item Peak calling parameter optimization (e.g., MACS2 \code{--extsize})
#'   \item Quality control assessment (checking if fragment size matches
#'         expected values)
#'   \item Downstream analysis (e.g., shifting reads by half the fragment
#'         length for visualization)
#' }
#' 
#' @param bamfiles A character vector of file paths to BAM files to be
#'        processed. The BAM files must be indexed (have corresponding .bai
#'        files). Multiple files can be processed in a single call.
#' @param index A character vector of index file names (without '.bai'
#'        extension) corresponding to \code{bamfiles}. Defaults to \code{bamfiles},
#'        meaning the index files are assumed to have the same base name as the
#'        BAM files with a '.bai' extension. If your index files have different
#'        names, specify them explicitly here.
#' @param plot A logical value. If \code{TRUE} (default), plots diagnostic
#'        visualizations:
#'        \itemize{
#'          \item For paired-end: Insert size distribution with smoothed curve
#'                and peak position
#'          \item For single-end: Cross-correlation function (CCF) with peak
#'                position
#'        }
#'        The plots help verify that the estimation is reasonable.
#' @param lag.max An integer specifying the maximum lag (in base pairs) for
#'        cross-correlation calculation in single-end mode. Default is
#'        \code{1000}. Larger values allow detection of longer fragments but
#'        increase computation time. See \code{\link[stats]{ccf}} for details.
#' @param minFragmentSize An integer specifying the minimal fragment size to
#'        consider. Fragments smaller than this are excluded to avoid phantom
#'        peaks in cross-correlation analysis. Default is \code{100} base pairs.
#'        This parameter helps filter out artifacts from very short fragments.
#' @param ... Additional arguments (currently not used).
#' 
#' @return Returns a numeric vector of estimated fragment lengths (in base
#'        pairs), one value per BAM file in \code{bamfiles}. The values
#'        represent:
#'        \itemize{
#'          \item For paired-end: The peak of the insert size distribution
#'                (most common fragment length)
#'          \item For single-end: The lag with maximum cross-correlation plus
#'                the mean read length (fragment length = shift + read length)
#'        }
#' 
#' @details
#' 
#' \strong{Paired-end mode:}
#' For paired-end data, the function:
#' \enumerate{
#'   \item Tests if the BAM file contains paired-end reads using
#'         \code{testPairedEndBam}
#'   \item Extracts insert sizes from proper pairs (both mates mapped, correct
#'         orientation, same chromosome)
#'   \item Filters to first mate reads only, excludes secondary alignments and
#'         low-quality reads
#'   \item Creates a histogram of insert sizes (absolute values)
#'   \item Smooths the distribution using LOESS regression (span = 1/3)
#'   \item Sets values below \code{minFragmentSize} to zero to avoid phantom
#'         peaks
#'   \item Finds the position with maximum smoothed value (peak of distribution)
#'   \item Returns the peak position as the estimated fragment length
#' }
#' 
#' \strong{Single-end mode:}
#' For single-end data, the function uses cross-correlation analysis:
#' \enumerate{
#'   \item Counts reads per chromosome and selects a representative chromosome
#'         (one with median read depth)
#'   \item Reads all alignments from the selected chromosome (excluding
#'         unmapped, secondary, and low-quality reads)
#'   \item Splits reads by strand and computes coverage for plus and minus
#'         strands separately
#'   \item Trims leading and trailing zeros from both strand coverages
#'   \item If the chromosome is very long (>100,000 bp), samples a representative
#'         region (50,000 bp block with highest average signal)
#'   \item Performs cross-correlation analysis between minus-strand and
#'         plus-strand coverage using \code{\link[stats]{ccf}}
#'   \item Finds the lag (shift) with maximum correlation (excluding lags <
#'         \code{minFragmentSize})
#'   \item Calculates fragment length as: lag + mean read width
#' }
#' 
#' \strong{Why cross-correlation for single-end?}
#' In ChIP-seq, reads from the same fragment are sequenced from opposite ends.
#' For a transcription factor binding site, reads cluster on both strands with
#' a characteristic spacing equal to the fragment length. Cross-correlation
#' finds this optimal shift by measuring how well the minus-strand coverage
#' aligns with the plus-strand coverage when shifted.
#' 
#' \strong{Chromosome selection:}
#' For single-end data, the function selects a chromosome with median read
#' depth to ensure:
#' \itemize{
#'   \item Sufficient read coverage for reliable cross-correlation
#'   \item Representative of the overall data (not too sparse or too dense)
#'   \item Reasonable computation time (avoids very long chromosomes)
#' }
#' 
#' \strong{Read filtering:}
#' The function applies strict filtering to ensure quality:
#' \itemize{
#'   \item \strong{Paired-end}: Proper pairs only, first mate, no secondary
#'         alignments, passing quality controls
#'   \item \strong{Single-end}: Unmapped reads excluded, no secondary
#'         alignments, passing quality controls
#' }
#' 
#' \strong{Plotting:}
#' When \code{plot = TRUE}, diagnostic plots are generated:
#' \itemize{
#'   \item \strong{Paired-end}: Shows insert size histogram (raw counts) with
#'         smoothed curve overlay and vertical line at the estimated peak
#'   \item \strong{Single-end}: Shows cross-correlation function (CCF) with
#'         vertical line at the estimated lag
#' }
#' These plots help verify that the estimation is reasonable and identify
#' potential issues (e.g., multiple peaks, no clear peak).
#' 
#' \strong{Performance considerations:}
#' \itemize{
#'   \item Paired-end mode is generally faster (only needs to read insert sizes)
#'   \item Single-end mode requires reading full alignments and computing
#'         coverage, which can be slow for large files
#'   \item For very long chromosomes, the function automatically samples a
#'         representative region to reduce computation time
#' }
#' 
#' @note
#' \itemize{
#'   \item BAM files must be indexed (have .bai files)
#'   \item The function automatically detects paired-end vs single-end data
#'   \item For single-end data, only one chromosome is analyzed (selected by
#'         median depth)
#'   \item The estimation may be inaccurate if:
#'         \itemize{
#'           \item Read coverage is too low
#'           \item Fragment size distribution is multimodal
#'           \item Data quality is poor
#'         }
#'   \item The function uses absolute insert sizes for paired-end (handles both
#'         forward and reverse pairs)
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link[Rsamtools]{testPairedEndBam}} for detecting paired-end
#'         data
#'   \item \code{\link[stats]{ccf}} for cross-correlation analysis
#'   \item \code{\link[stats]{loess.smooth}} for LOESS smoothing
#'   \item MACS2 documentation for using fragment length in peak calling
#' }
#' 
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom Rsamtools testPairedEndBam scanBam ScanBamParam scanBamFlag
#' scanBamHeader countBam 
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom stats loess.smooth ts ccf
#' @examples
#' 
#' # Example 1: Estimate fragment length for a single BAM file
#' \dontrun{
#' bamfile <- "sample.bam"
#' frag_length <- estFragmentLength(bamfile)
#' frag_length
#' }
#' 
#' # Example 2: Estimate fragment length for multiple BAM files
#' \dontrun{
#' bamfiles <- c("sample1.bam", "sample2.bam", "sample3.bam")
#' frag_lengths <- estFragmentLength(bamfiles)
#' frag_lengths
#' }
#' 
#' # Example 3: Without plotting (faster for batch processing)
#' \dontrun{
#' frag_length <- estFragmentLength("sample.bam", plot = FALSE)
#' }
#' 
#' # Example 4: Adjust parameters for single-end data
#' \dontrun{
#' # Increase lag.max to detect longer fragments
#' frag_length <- estFragmentLength("single_end.bam", 
#'                                  lag.max = 2000L,
#'                                  minFragmentSize = 150L)
#' }
#' 
#' # Example 5: Using with example data (if available)
#' if(interactive() || Sys.getenv("USER")=="jou"){
#'     path <- system.file("extdata", "reads", package="MMDiffBamSubset")
#'     if(file.exists(path)){
#'         WT.AB2 <- file.path(path, "WT_2.bam")
#'         Null.AB2 <- file.path(path, "Null_2.bam")
#'         Resc.AB2 <- file.path(path, "Resc_2.bam")
#'         estFragmentLength(c(WT.AB2, Null.AB2, Resc.AB2))
#'     }
#' }
#' 
estFragmentLength <- function(bamfiles,
                              index = bamfiles,
                              plot = TRUE,
                              lag.max = 1000,
                              minFragmentSize = 100,
                              ...) {
    #message("The fragment size is being calculated for DNA-seq.")
    res <- mapply(function(f, i) {
        if (suppressMessages(testPairedEndBam(f, index = i))) {
            isize <- sort(
                abs(scanBam(f,
                            index = i,
                            param = ScanBamParam(
                                flag = scanBamFlag(
                                    isPaired = TRUE,
                                    isProperPair = TRUE,
                                    isSecondaryAlignment = FALSE,
                                    isFirstMateRead = TRUE,
                                    isNotPassingQualityControls = FALSE
                                ),
                                what = "isize"
                            ))[[1]]$isize
                )
            )
            isize <- table(isize)
            isize.x <- as.numeric(names(isize))
            isize.y <- isize
            ## do smooth and select the max point
            isize.smooth <- loess.smooth(isize.x,
                                          isize.y,
                                          span = 1/3,
                                          evaluation = length(isize.x))
            isize.smooth$y[seq.int(minFragmentSize)] <- 0
            pos <- round(isize.smooth$x[which.max(isize.smooth$y)],
                         digits = 0)
            if (plot) {
                try({
                    plot(isize.x, isize.y,
                         xlab = "Lag",
                         ylab = "ACF",
                         main = f,
                         type = "l")
                    abline(v = pos, col = "red")
                })
            }
            pos
        } else {
            ## count reads for all chromosomes
            seqinfo <- lapply(scanBamHeader(f, index = i),
                               function(.ele) .ele$targets)
            seqn <- sort(table(unlist(lapply(seqinfo, names))))
            seqn <- names(seqn[which(seqn == max(seqn))])
            seqinfo <- unique(do.call(rbind, lapply(seqinfo, `[`, seqn)))[1, ]
            seqinfo.gr <- GRangesList(lapply(seq_along(seqinfo), function(.ele) {
                GRanges(names(seqinfo)[.ele], IRanges(1L, seqinfo[.ele]))
            }))
            names(seqinfo.gr) <- names(seqinfo)
            cnts <- countBam(f,
                             index = i,
                             param = ScanBamParam(
                                 flag = scanBamFlag(
                                     isPaired = FALSE,
                                     isUnmappedQuery = FALSE,
                                     isSecondaryAlignment = FALSE,
                                     isNotPassingQualityControls = FALSE
                                 ),
                                 which = seqinfo.gr
                             ))
            cnts <- cnts[cnts$nucleotides > 0L, , drop = FALSE]
            if (nrow(cnts) < 1L) {
                stop("no reads detected", call. = FALSE)
            }
            depth <- cnts$nucleotides / cnts$width
            ## select middle one
            select.chr <- as.character(cnts$space[which.min(abs(depth - median(depth)))])
            # read reads
            reads <- readGAlignments(f,
                                     index = i,
                                     param = ScanBamParam(
                                         flag = scanBamFlag(
                                             isPaired = FALSE,
                                             isUnmappedQuery = FALSE,
                                             isSecondaryAlignment = FALSE,
                                             isNotPassingQualityControls = FALSE
                                         ),
                                         which = seqinfo.gr[[select.chr]],
                                         what = scanBamWhat()
                                     ))
            #convert to GRanges
            reads.gr <- as(reads, "GRanges")
            mcols(reads.gr) <- NULL
            reads.gr <- split(reads.gr, strand(reads.gr))
            reads.gr <- sapply(reads.gr, coverage)
            reads.gr.pos <- as.integer(reads.gr[["+"]][[select.chr]])
            reads.gr.neg <- as.integer(reads.gr[["-"]][[select.chr]])
            ## remove 0 from both end
            pos0 <- max(which(reads.gr.pos != 0)[1], which(reads.gr.neg != 0)[1])
            pos1 <- length(reads.gr.pos) -
                max(which(rev(reads.gr.pos) != 0)[1],
                    which(rev(reads.gr.neg) != 0)[1])
            if (pos1 > pos0 + 5000L) {
                reads.gr.pos <- reads.gr.pos[pos0:pos1]
                reads.gr.neg <- reads.gr.neg[pos0:pos1]
            }
            if (length(reads.gr.pos) > 100000L) {
                block <- 50000L
                len <- ceiling(length(reads.gr.pos) / block)
                signal <- split(reads.gr.pos,
                                 rep(formatC(seq_len(len),
                                              width = nchar(as.character(len)),
                                              flag = "0"),
                                     each = block)[seq_along(reads.gr.pos)])
                signal <- vapply(signal, mean, FUN.VALUE = numeric(1))
                signal <- signal[-len]
                sid <- which.max(signal)
                pos0 <- (sid - 1L) * block + 1L
                pos1 <- min(sid * block, length(reads.gr.pos))
                reads.gr.pos <- reads.gr.pos[pos0:pos1]
                reads.gr.neg <- reads.gr.neg[pos0:pos1]
            }
            ## cross-correlation
            reads.gr.pos.ts <- ts(reads.gr.pos)
            reads.gr.neg.ts <- ts(reads.gr.neg)
            ccf <- ccf(reads.gr.neg.ts,
                       reads.gr.pos.ts,
                       type = "correlation",
                       plot = FALSE,
                       lag.max = lag.max)
            keep <- ccf$lag >= 0L
            lag <- ccf$lag[keep]
            acf1 <- acf <- ccf$acf[keep]
            acf[seq.int(minFragmentSize)] <- 0L
            pos <- lag[which.max(acf)]
            if (plot) {
                try({
                    plot(lag, acf1,
                         xlab = "Lag",
                         ylab = "ACF",
                         main = paste(f, "insertion size (not including reads length)"),
                         type = "l")
                    abline(v = pos, col = "red")
                }, silent = TRUE)
            }
            pos + mean(width(reads))
        }
    }, bamfiles, index, SIMPLIFY = TRUE)
    res
}
