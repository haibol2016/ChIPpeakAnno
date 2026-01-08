#' Generate random peaks for permutation testing
#' 
#' @description 
#' Generates random genomic peaks for permutation testing. The function samples
#' positions from distance-binned genomic regions (provided by
#' \code{\link{preparePool}}) and assigns random widths based on the width
#' distribution of the input peaks. This ensures that random peaks match both
#' the distance distribution (from annotation features) and the width
#' distribution of observed peaks, making permutation tests more appropriate.
#' 
#' This is an internal function used by \code{\link{peakPermTest}} via the
#' \code{regioneR::permTest} framework. It is called once per permutation to
#' generate a set of random peaks that replace the observed peaks.
#' 
#' @param A A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        representing the original peaks. The width distribution of these peaks
#'        is used to generate random widths for the random peaks. The number of
#'        random peaks generated equals \code{length(A)}.
#' @param grs A \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}
#'        object containing genomic regions grouped by distance bins. Typically
#'        created by \code{\link{preparePool}}. Each element of the list
#'        corresponds to one distance bin (e.g., regions at -5000 to -4900 bp
#'        from TSS). Regions in each bin are used as sampling pools for random
#'        positions.
#' @param N An integer vector specifying the number of random peaks to generate
#'        from each corresponding element of \code{grs}. Must have the same
#'        length as \code{grs}. Values in \code{N} correspond to the number of
#'        observed peaks in each distance bin (from the binding distribution),
#'        ensuring that random peaks match the observed distance distribution.
#'        The sum of \code{N} should equal \code{length(A)}.
#' @param ... Additional arguments (currently unused, reserved for future
#'        extensions)
#' 
#' @return Returns a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        containing randomly generated peaks. The returned object has:
#'        \itemize{
#'          \item Same number of peaks as \code{length(A)}
#'          \item Positions sampled uniformly from regions in \code{grs},
#'                distributed across distance bins according to \code{N}
#'          \item Widths sampled from an exponential distribution matching the
#'                width distribution of \code{A}
#'          \item All ranges trimmed to valid genomic coordinates
#'        }
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates that all elements in \code{grs} have at least one region
#'         and that \code{length(N) == length(grs)}
#'   \item For each distance bin i:
#'         \itemize{
#'           \item Samples \code{N[i]} positions uniformly from \code{grs[[i]]}
#'                 using \code{\link{runifGR}}
#'           \item Positions are sampled uniformly across all positions in all
#'                 regions (weighted by region width)
#'         }
#'   \item Combines all sampled positions into a single GRanges object
#'   \item Generates random widths from an exponential distribution:
#'         \itemize{
#'           \item Rate parameter: \code{1 / (mean(width(A)) - min(width(A)))}
#'           \item Shift: \code{min(width(A))}
#'           \item This ensures the mean of random widths matches the mean of
#'                 observed widths
#'         }
#'   \item Centers each random peak around its sampled position:
#'         \itemize{
#'           \item Calculates half-width: \code{floor(width / 2)}
#'           \item Shifts start position upstream by half-width
#'           \item Sets width to \code{2 * half-width} (ensures even widths)
#'         }
#'   \item Trims all peaks to valid genomic coordinates
#' }
#' 
#' \strong{Position sampling:}
#' Positions are sampled uniformly across all positions in all regions within
#' each distance bin. This means:
#' \itemize{
#'   \item Larger regions contribute more positions to the sampling pool
#'   \item Each base pair in the regions has equal probability of being sampled
#'   \item Sampling is done with replacement (same position can be sampled
#'         multiple times)
#' }
#' 
#' \strong{Width generation:}
#' Random widths are generated from an exponential distribution with:
#' \itemize{
#'   \item \strong{Rate}: \code{1 / (mean(width(A)) - min(width(A)))}
#'   \item \strong{Shift}: \code{min(width(A))}
#'   \item This ensures:
#'         \itemize{
#'           \item Mean of random widths = mean of observed widths
#'           \item Minimum of random widths = minimum of observed widths
#'           \item Distribution shape approximates the observed width distribution
#'         }
#' }
#' 
#' \strong{Peak centering:}
#' After sampling positions and generating widths, peaks are centered around
#' their sampled positions:
#' \itemize{
#'   \item The sampled position becomes the center of the peak
#'   \item Half-width is calculated and used to extend the peak symmetrically
#'   \item Widths are forced to be even (2 * half-width) to ensure symmetric
#'         centering
#' }
#' 
#' \strong{Use in permutation testing:}
#' This function is called by \code{regioneR::permTest} as the
#' \code{randomize.function} parameter. For each permutation:
#' \itemize{
#'   \item \code{randPeaks} generates random peaks matching the binding and
#'         width distributions
#'   \item Overlaps between random peaks and the test set are counted
#'   \item The distribution of random overlap counts is compared to the
#'         observed overlap count
#' }
#' 
#' @author Internal utility function
#' @keywords internal
#' @importFrom stats rexp
#' @seealso \code{\link{peakPermTest}} for the permutation test framework,
#'          \code{\link{preparePool}} for creating the sampling pool,
#'          \code{\link{runifGR}} for uniform position sampling
randPeaks <- function(A, grs, N, ...) {
    len <- lengths(grs)
    if (any(len < 1L)) {
        stop("All elements of 'grs' must have length greater than 0", 
             call. = FALSE)
    }
    if (length(N) != length(grs)) {
        stop("Length of 'N' and 'grs' must be identical", call. = FALSE)
    }
    
    s <- mapply(runifGR, grs, N)
    s <- unlist(GRangesList(s))
    
    wid <- width(A)
    # Generate random widths from exponential distribution
    # Mean is adjusted to match the distribution of input widths
    wid <- rexp(length(s), 1L / (mean(wid) - min(wid))) + min(wid)
    halfwid <- floor(wid / 2L)
    suppressWarnings({
        start(s) <- start(s) - halfwid
        width(s) <- 2L * halfwid
    })
    s <- trim(s)
    s
}

#' Sample random positions from genomic regions
#' 
#' @description 
#' Samples n random positions uniformly from a set of genomic regions. Positions
#' are sampled uniformly across all base pairs in all regions, meaning larger
#' regions contribute more positions to the sampling pool. This is an internal
#' helper function used by \code{\link{randPeaks}}.
#' 
#' @param grNoN A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        containing genomic regions from which to sample positions. Can contain
#'        multiple regions (they will be treated as a continuous pool).
#' @param n An integer specifying the number of positions to sample. Must be >= 1.
#'        Sampling is done with replacement, so \code{n} can be larger than the
#'        total number of positions in \code{grNoN}.
#' 
#' @return Returns a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        with \code{n} ranges, each of width 1 bp. Each range represents one
#'        randomly sampled position from the input regions.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Randomly shuffles the order of regions in \code{grNoN} (for
#'         randomization)
#'   \item Calculates the total length of all regions combined
#'   \item Creates a continuous coordinate system by concatenating all regions
#'         end-to-end
#'   \item Samples \code{n} positions uniformly from this continuous coordinate
#'         system (with replacement)
#'   \item Maps sampled positions back to the original genomic coordinates:
#'         \itemize{
#'           \item Determines which region contains each sampled position
#'           \item Calculates the offset within that region
#'           \item Shifts the region by the offset to get the final position
#'         }
#'   \item Trims positions to valid genomic coordinates
#'   \item Sets all widths to 1 bp
#' }
#' 
#' \strong{Uniform sampling:}
#' Positions are sampled uniformly across all base pairs in all regions. This
#' means:
#' \itemize{
#'   \item A region of width 1000 bp has 10x more probability of contributing a
#'         position than a region of width 100 bp
#'   \item Each base pair in the regions has equal probability of being sampled
#'   \item Sampling is done with replacement, so the same position can be
#'         sampled multiple times
#' }
#' 
#' \strong{Coordinate mapping:}
#' The function treats all regions as if they were concatenated end-to-end in a
#' continuous coordinate system:
#' \itemize{
#'   \item Region 1: positions 1 to width1
#'   \item Region 2: positions width1+1 to width1+width2
#'   \item Region 3: positions width1+width2+1 to width1+width2+width3
#'   \item etc.
#' }
#' After sampling from this continuous system, positions are mapped back to
#' their original genomic coordinates by shifting the appropriate region.
#' 
#' @keywords internal
runifGR <- function(grNoN, n) {
    grNoN <- grNoN[sample.int(length(grNoN))]
    wid <- width(grNoN)
    totL <- sum(as.numeric(wid))
    ends <- cumsum(as.numeric(wid))
    starts <- c(1L, ends[-length(ends)] + 1L)
    idx <- sample.int(n = totL, size = n, replace = TRUE)
    pos <- cut(idx, breaks = c(0L, ends), labels = seq_along(grNoN))
    pos <- as.integer(as.character(pos))
    off <- idx - starts[pos] - 1L
    gr <- grNoN[pos]
    suppressWarnings(gr <- shift(gr, shift = off))
    gr <- trim(gr)
    width(gr) <- 1L
    gr
}