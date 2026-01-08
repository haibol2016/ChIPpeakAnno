#' Count reads in genome-wide windows
#' 
#' @description 
#' Splits the genome into fixed-size tiling windows and counts overlapping reads
#' in each window. This function extends \code{\link[GenomicAlignments]{summarizeOverlaps}}
#' by automatically generating genome-wide tiling windows from sequence length
#' information, eliminating the need to manually create windows.
#' 
#' This function is particularly useful for:
#' \itemize{
#'   \item Creating genome-wide coverage tracks
#'   \item Identifying regions of high read density
#'   \item Generating uniform bins for comparative analysis
#'   \item Creating input for visualization tools (e.g., genome browsers)
#'   \item Downsampling large datasets to fixed-size windows
#' }
#' 
#' @param reads A \code{\link[GenomicRanges]{GRanges}},
#'        \code{\link[GenomicRanges]{GRangesList}},
#'        \code{\link[GenomicAlignments]{GAlignments}},
#'        \code{\link[GenomicAlignments]{GAlignmentsList}},
#'        \code{\link[GenomicAlignments]{GAlignmentPairs}}, or
#'        \code{\link[Rsamtools]{BamFileList}} object representing the
#'        sequencing reads to be counted. Can be a single sample or multiple
#'        samples.
#' @param genome An object with sequence length information from which to
#'        extract chromosome/contig lengths. Can be:
#'        \itemize{
#'          \item A \code{BSgenome} object (e.g., \code{BSgenome.Hsapiens.UCSC.hg19})
#'          \item A \code{\link[GenomeInfoDb]{Seqinfo}} object
#'          \item A \code{\link[GenomicRanges]{GRanges}} object with
#'                \code{seqlengths} metadata
#'          \item Any object with a \code{seqlengths} method that returns a
#'                named numeric vector of sequence lengths
#'        }
#'        The function uses \code{seqlengths(genome)} to extract lengths and
#'        will stop with an error if any lengths are \code{NA}.
#' @param windowSize An integer specifying the size of each window in base
#'        pairs. Default is \code{1e6L} (1 Mb). Windows are created across all
#'        chromosomes/contigs in the genome. Smaller windows provide finer
#'        resolution but create more windows and require more computation.
#' @param step An integer specifying the step size for sliding windows in base
#'        pairs. Default is \code{1e6L} (1 Mb), which creates non-overlapping
#'        windows when equal to \code{windowSize}. When \code{step < windowSize},
#'        windows overlap by \code{windowSize - step} base pairs. Smaller step
#'        sizes create more windows and provide smoother coverage but require
#'        more computation.
#' @param keepPartialWindow A logical value. If \code{TRUE}, includes the last
#'        partial window at the end of each chromosome/contig (even if it's
#'        smaller than \code{windowSize}). If \code{FALSE} (default), partial
#'        windows are excluded. Partial windows occur when the chromosome length
#'        is not evenly divisible by \code{step}.
#' @param mode A function specifying the count method to use. Default is
#'        \code{countByOverlaps}, which counts overlaps between reads and windows
#'        using \code{countOverlaps(features, reads, ignore.strand=ignore.strand)}.
#'        See \code{\link[GenomicAlignments]{summarizeOverlaps}} for other
#'        options (e.g., \code{Union}, \code{IntersectionStrict},
#'        \code{IntersectionNotEmpty}).
#' @param ... Additional arguments passed to
#'        \code{\link[GenomicAlignments]{summarizeOverlaps}}, such as
#'        \code{ignore.strand}, \code{inter.feature}, etc.
#' 
#' @return Returns a \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#'        object with:
#'        \itemize{
#'          \item \code{assays}: A \code{SimpleList} containing a \code{counts}
#'                matrix with read counts per window. Rows correspond to tiling
#'                windows, columns correspond to samples/read files (from
#'                \code{reads}). Values are integer counts of overlapping reads.
#'          \item \code{rowRanges}: A \code{\link[GenomicRanges]{GRanges}}
#'                object containing the tiling windows. Each window has:
#'                \itemize{
#'                  \item Chromosome/contig name (seqnames)
#'                  \item Start and end positions
#'                  \item Strand (typically "*" for unstranded windows)
#'                  \item Width equal to \code{windowSize} (or smaller for
#'                        partial windows if \code{keepPartialWindow = TRUE})
#'                }
#'          \item \code{colData}: Sample metadata from the input \code{reads}
#'                object (if available).
#'        }
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Extracts sequence lengths from \code{genome} using
#'         \code{seqlengths(genome)}
#'   \item Validates that all sequence lengths are non-NA (stops with error if
#'         any are missing)
#'   \item Creates a GRanges object for each chromosome/contig spanning from
#'         position 1 to the sequence length
#'   \item Uses \code{\link{tileGRanges}} to create tiling windows:
#'         \itemize{
#'           \item Windows of size \code{windowSize}
#'           \item Sliding by \code{step} base pairs
#'           \item Optionally including partial windows at chromosome ends
#'         }
#'   \item Counts reads in each window using
#'         \code{\link[GenomicAlignments]{summarizeOverlaps}}
#'   \item Returns a SummarizedExperiment with counts and window coordinates
#' }
#' 
#' \strong{Window creation:}
#' Windows are created across all chromosomes/contigs in the genome:
#' \itemize{
#'   \item Each chromosome is divided into windows independently
#'   \item Windows are created by sliding a window of size \code{windowSize} by
#'         \code{step} base pairs
#'   \item When \code{step == windowSize}, windows are non-overlapping
#'   \item When \code{step < windowSize}, windows overlap by
#'         \code{windowSize - step} base pairs
#'   \item Partial windows at chromosome ends are included or excluded based on
#'         \code{keepPartialWindow}
#' }
#' 
#' Example: For a 5 Mb chromosome with \code{windowSize = 1e6} and
#' \code{step = 5e5}:
#' \itemize{
#'   \item Window 1: 1-1,000,000
#'   \item Window 2: 500,001-1,500,000 (overlaps with window 1)
#'   \item Window 3: 1,000,001-2,000,000
#'   \item ... (continues to end of chromosome)
#' }
#' 
#' \strong{Sequence length extraction:}
#' The function uses \code{seqlengths(genome)} to extract chromosome lengths.
#' This method works with:
#' \itemize{
#'   \item \code{BSgenome} objects: Automatically extracts sequence lengths
#'   \item \code{Seqinfo} objects: Uses stored sequence length information
#'   \item \code{GRanges} objects: Uses \code{seqlengths} metadata
#'   \item Any object with a \code{seqlengths} method
#' }
#' 
#' \strong{Use cases:}
#' \itemize{
#'   \item \strong{Genome-wide coverage}: Create uniform bins for visualizing
#'         read coverage across the entire genome
#'   \item \strong{Comparative analysis}: Compare read density across samples
#'         using identical windows
#'   \item \strong{Downsampling}: Reduce large datasets to fixed-size windows
#'         for faster analysis
#'   \item \strong{Visualization}: Generate data for genome browser tracks or
#'         coverage plots
#' }
#' 
#' \strong{Performance considerations:}
#' \itemize{
#'   \item Smaller \code{windowSize} and \code{step} values create more windows
#'         and increase computation time
#'   \item For large genomes, even 1 Mb windows can create thousands of windows
#'   \item The function processes all chromosomes and all samples, so large
#'         datasets may require significant memory and time
#'   \item Consider using larger windows for initial exploratory analysis, then
#'         refine with smaller windows for regions of interest
#' }
#' 
#' @note
#' \itemize{
#'   \item All sequence lengths must be non-NA; the function stops with an error
#'         if any are missing
#'   \item Windows are created independently for each chromosome/contig
#'   \item Partial windows are excluded by default; set
#'         \code{keepPartialWindow = TRUE} to include them
#'   \item The function uses \code{tileGRanges} internally for window creation
#'   \item When \code{step == windowSize}, windows are non-overlapping
#'   \item When \code{step < windowSize}, windows overlap
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link[GenomicAlignments]{summarizeOverlaps}} for the underlying
#'         read counting function
#'   \item \code{\link{tileGRanges}} for the window creation function used
#'         internally
#'   \item \code{\link[SummarizedExperiment]{SummarizedExperiment}} for the
#'         output object class
#'   \item \code{\link[GenomeInfoDb]{seqlengths}} for extracting sequence
#'         lengths
#' }
#' 
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom GenomicAlignments summarizeOverlaps
#' @examples
#' 
#' # Example 1: Basic usage with BAM files
#' fls <- list.files(system.file("extdata", package="GenomicAlignments"),
#'               recursive=TRUE, pattern="*bam$", full=TRUE)
#' names(fls) <- basename(fls)
#' genome <- GRanges(seqlengths = c(chr2L=7000, chr2R=10000))
#' se <- tileCount(fls, genome, windowSize=1000, step=500)
#' 
#' # Example 2: Using BSgenome object
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' se <- tileCount(fls, Hsapiens, windowSize=1e6, step=1e6)
#' }
#' 
#' # Example 3: Non-overlapping windows (step == windowSize)
#' \dontrun{
#' # Create 1 Mb non-overlapping windows
#' se <- tileCount(fls, genome, windowSize=1e6, step=1e6)
#' }
#' 
#' # Example 4: Overlapping windows for smoother coverage
#' \dontrun{
#' # 1 Mb windows with 500 kb step (50% overlap)
#' se <- tileCount(fls, genome, windowSize=1e6, step=5e5)
#' }
#' 
#' # Example 5: Smaller windows for finer resolution
#' \dontrun{
#' # 10 kb windows with 5 kb step
#' se <- tileCount(fls, genome, windowSize=10000, step=5000)
#' }
#' 
#' # Example 6: Include partial windows
#' \dontrun{
#' # Include last partial window at chromosome ends
#' se <- tileCount(fls, genome, windowSize=1e6, step=1e6,
#'                 keepPartialWindow = TRUE)
#' }
#' 
#' # Example 7: Using GRanges as reads
#' \dontrun{
#' reads <- GRanges("chr1", IRanges(1:1000000, width=50))
#' se <- tileCount(reads, genome, windowSize=10000, step=5000)
#' }
#' 
#' # Example 8: Accessing results
#' \dontrun{
#' se <- tileCount(fls, genome, windowSize=1000, step=500)
#' 
#' # Get count matrix
#' counts <- assay(se)
#' 
#' # Get window coordinates
#' windows <- rowRanges(se)
#' 
#' # Get sample metadata
#' sample_info <- colData(se)
#' }
#' 
tileCount <- function(reads,
                     genome,
                     windowSize = 1e6L,
                     step = 1e6L,
                     keepPartialWindow = FALSE,
                     mode = countByOverlaps,
                     ...) {
    targetRegions <- seqlengths(genome)
    if (any(is.na(targetRegions))) {
        stop("Cannot get sequence lengths from genome.", call. = FALSE)
    }
    
    targetRegions <- GRanges(
        names(targetRegions),
        IRanges(1L, width = targetRegions)
    )
    
    tileTargetRegions <- tileGRanges(targetRegions,
                                     windowSize,
                                     step,
                                     keepPartialWindow)
    se <- summarizeOverlaps(features = tileTargetRegions,
                            reads = reads,
                            mode = mode,
                            ...)
    se
}
