#' Plot cumulative percentage tag allocation
#' 
#' @description 
#' Plots the cumulative percentage of tag allocation in paired samples to
#' visualize differences in signal distribution. This function is useful for
#' comparing ChIP-seq samples and identifying systematic biases or differences
#' in tag distribution across genomic regions.
#' 
#' The function performs the following steps:
#' \enumerate{
#'   \item Tiles the input genomic regions (\code{gr}) into bins of specified
#'         width using \code{\link{tileGRanges}}
#'   \item Counts reads from each BAM file overlapping each bin using
#'         \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}
#'   \item For each sample (excluding the input), calculates cumulative
#'         percentages of tags by sorting bins based on sample signal and input
#'         signal
#'   \item Plots cumulative percentage curves showing the percentage of bins
#'         (x-axis) versus the percentage of tags (y-axis) for both input and
#'         sample
#'   \item Identifies and marks key points: where input signal becomes
#'         significant (> 1/binWidth) and where the maximum distance between
#'         input and sample curves occurs
#' }
#' 
#' This analysis is particularly useful for:
#' \itemize{
#'   \item Detecting systematic biases between ChIP-seq samples and their
#'         controls
#'   \item Identifying regions with differential tag enrichment
#'   \item Quality control of ChIP-seq experiments
#'   \item Normalization assessment (see references)
#' }
#' 
#' @param bamfiles A character vector of BAM file paths. The files should be
#'        indexed (have corresponding .bai files). All BAM files will be
#'        analyzed, with one file designated as the input/control.
#' @param gr An object of \link[GenomicRanges:GRanges-class]{GRanges} specifying
#'        the genomic regions to analyze. The function will tile these regions
#'        into bins and count reads overlapping each bin.
#' @param input An integer specifying which BAM file (by position in
#'        \code{bamfiles}) should be used as the input/control sample. Default
#'        is \code{1L} (first file). The input sample is used as a reference for
#'        comparison with all other samples.
#' @param binWidth An integer specifying the width of each bin (in base pairs)
#'        for tiling the genomic regions. Default is \code{1000L} (1kb bins).
#'        Smaller bins provide higher resolution but require more computation.
#' @param ... Additional parameters to be passed to
#'        \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}},
#'        such as \code{ignore.strand}, \code{inter.feature}, \code{singleEnd},
#'        etc. See the \code{summarizeOverlaps} documentation for details.
#' 
#' @return Returns (invisibly) a list of data frames, one for each sample
#'        (excluding the input). The function primarily generates plots, but the
#'        data can be captured by assignment (e.g., \code{result <-
#'        cumulativePercentage(...)}). Each data frame contains:
#'        \itemize{
#'          \item \code{Rank}: The rank/position of bins after sorting by signal
#'                (ranging from 1 to number of bins)
#'          \item \code{input_name}: Cumulative percentage of tags from the input
#'                sample (ranging from 0 to 1)
#'          \item \code{sample_name}: Cumulative percentage of tags from the
#'                sample (ranging from 0 to 1)
#'        }
#'        The bins are sorted by sample signal (ascending), then by input signal
#'        (ascending) as a tie-breaker. This allows visualization of how tags are
#'        distributed across bins of different signal intensities.
#'        
#'        The function also generates plots (one per sample) showing:
#'        \itemize{
#'          \item Cumulative percentage curves for both input and sample
#'          \item Vertical lines (yellow-green, dashed) marking:
#'                \itemize{
#'                  \item The point where input signal becomes significant
#'                        (> 1/binWidth)
#'                  \item The point of maximum distance between input and sample
#'                        curves
#'                }
#'        }
#' @details
#' 
#' \strong{Interpreting the plots:}
#' 
#' The cumulative percentage plots show how tags are distributed across bins:
#' \itemize{
#'   \item \strong{X-axis}: Percentage of bins (0-100\%), sorted by signal
#'         intensity
#'   \item \strong{Y-axis}: Cumulative percentage of tags (0-100\%)
#'   \item \strong{Diagonal line}: Represents uniform distribution (if tags were
#'         evenly distributed, the curve would follow the diagonal)
#'   \item \strong{Curves above diagonal}: Indicate enrichment in high-signal bins
#'   \item \strong{Curves below diagonal}: Indicate enrichment in low-signal bins
#' }
#' 
#' The vertical lines mark:
#' \itemize{
#'   \item \strong{Green line (left)}: Point where input signal becomes
#'         significant (> 1 tag per bin on average)
#'   \item \strong{Green line (right)}: Point of maximum divergence between
#'         input and sample curves, indicating the region with the largest
#'         difference in tag distribution
#' }
#' 
#' \strong{Expected patterns:}
#' \itemize{
#'   \item Well-normalized samples should have curves that are close to each
#'         other and near the diagonal
#'   \item ChIP samples typically show curves above the diagonal, indicating
#'         enrichment in specific regions
#'   \item Large differences between input and sample curves suggest differential
#'         binding or normalization issues
#' }
#' 
#' @author Jianhong Ou
#' @references Normalization, bias correction, and peak calling for ChIP-seq.
#' Aaron Diaz, Kiyoub Park, Daniel A. Lim, Jun S. Song. Stat Appl Genet Mol
#' Biol. Author manuscript; available in PMC 2012 May 3. Published in final
#' edited form as: Stat Appl Genet Mol Biol. 2012 Mar 31; 11(3):
#' 10.1515/1544-6115.1750
#' /j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml. Published
#' online 2012 Mar 31. doi: 10.1515/1544-6115.1750 PMCID: PMC3342857
#' @importFrom SummarizedExperiment assays
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom graphics abline axis legend matlines par plot
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Example 1: Basic usage with multiple BAM files
#' path <- system.file("extdata", "reads", package="MMDiffBamSubset")
#' files <- dir(path, "bam$", full.names = TRUE)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' gr <- as(seqinfo(Hsapiens)["chr1"], "GRanges")
#' 
#' # First file is used as input by default
#' cumulativePercentage(files, gr)
#' 
#' ## Example 2: Specify which file is the input/control
#' # Use second file as input
#' cumulativePercentage(files, gr, input = 2L)
#' 
#' ## Example 3: Use smaller bins for higher resolution
#' cumulativePercentage(files, gr, binWidth = 500L)
#' 
#' ## Example 4: Analyze specific genomic regions
#' peaks <- GRanges("chr1", IRanges(start = 1000000, end = 2000000))
#' cumulativePercentage(files, peaks, binWidth = 2000L)
#' 
#' ## Example 5: Access the returned data
#' result <- cumulativePercentage(files, gr)
#' # View cumulative percentages for first sample
#' head(result[[1]])
#' } 
#'  
#' 
cumulativePercentage <- function(bamfiles, gr, input = 1L, binWidth = 1e3L,
                                 ...) {
    stopifnot(inherits(gr, "GRanges"))
    
    tileTargetRegions <- tileGRanges(gr, windowSize = binWidth, 
                                     step = binWidth)
    se <- summarizeOverlaps(features = tileTargetRegions, reads = bamfiles, ...)
    
    # Resample and calculate cumulative percentages
    sampleName <- colnames(assays(se)[[1L]])[-input]
    input_name <- colnames(assays(se)[[1L]])[input]
    
    sigBin <- lapply(sampleName, function(.n) {
        sig <- assays(se)[[1L]][, c(input_name, .n)]
        sig[order(sig[, .n], sig[, input_name]), ]
    })
    
    sigCumsum <- lapply(sigBin, function(.ele) {
        .ele <- apply(.ele, 2L, cumsum)
        .ele <- cbind(Rank = seq.int(nrow(.ele)), .ele)
        sweep(.ele, MARGIN = 2L, STATS = .ele[nrow(.ele), ], FUN = `/`)
    })
    
    pin <- par("pin")
    if (pin[2L] > 0L) {
        ratio <- 2^round(diff(log2(pin)))
        n <- length(sampleName)
        ncol <- ceiling(sqrt(n / ratio))
        nrow <- ceiling(n / ncol)
        op <- par(mfrow = c(nrow, ncol), pty = "s")
        on.exit(par(op))
        
        for (i in seq_len(n)) {
            # Plot cumulative percentage
            plot(c(0, 1), c(0, 1), type = "n", 
                 xlab = "% of bins", ylab = "% of tags", 
                 main = sampleName[i])
            matlines(x = sigCumsum[[i]][, 1L], 
                     y = sigCumsum[[i]][, -1L],
                     lty = 1L)
            
            zero <- which(sigCumsum[[i]][, 2L] > 1L / binWidth)
            maxDist <- which.max(abs(sigCumsum[[i]][, 2L] - 
                                    sigCumsum[[i]][, 3L]))
            
            if (length(zero) > 0L || length(maxDist) > 0L) {
                x <- numeric(0)
                if (length(zero) > 0L) {
                    x <- c(x, zero[1L])
                }
                if (length(maxDist) > 0L) {
                    x <- c(x, maxDist[1L])
                }
                x_tick <- sigCumsum[[i]][x, 1L]
                abline(v = x_tick, col = "yellowgreen", lty = 3L)
                axis(3, at = x_tick, labels = formatC(x_tick, digits = 2L), 
                     lwd = 0L, lwd.ticks = 1L)
            }
            legend("topleft", legend = colnames(sigCumsum[[i]])[-1L], 
                   col = seq_len(6L), lty = 1L, pch = NA, box.col = NA)
        }
    }
    
    invisible(sigCumsum)
}
