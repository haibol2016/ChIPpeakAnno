#' Plot signal distribution across features
#' 
#' @description 
#' Plots the average signal distribution across genomic features. The function
#' calculates the mean signal density across all features and visualizes it as
#' a line plot. This is useful for visualizing ChIP-seq or RNA-seq signal
#' patterns around features like TSS, gene bodies, or peaks, providing a
#' metagene-like visualization.
#' 
#' The function can work with either pre-computed signal matrices (from
#' \code{\link{featureAlignedSignal}}) or directly with RleList objects. When
#' RleList objects are provided, the function internally calls
#' \code{featureAlignedSignal} to compute the signal matrices.
#' 
#' @param cvglists Signal data in one of the following formats:
#'        \itemize{
#'          \item Output of \code{\link{featureAlignedSignal}}: A list of matrices
#'                where each matrix has \code{n.tile} columns (one per tile/bin)
#'                and rows corresponding to features
#'          \item A list of \code{\link[IRanges:AtomicList-class]{SimpleRleList}}
#'                or \code{\link[IRanges:AtomicList-class]{RleList}} objects
#'                containing coverage/signal data
#'          \item A single \code{SimpleRleList} or \code{RleList} object (will be
#'                automatically converted to a list)
#'        }
#'        Each element in the list represents a different sample or condition.
#'        List names will be used as labels in the legend.
#' @param feature.gr An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#'        with identical width for all ranges. The features around which signal
#'        will be plotted (e.g., TSS, gene centers, peak centers).
#' @param upstream An integer specifying the number of base pairs upstream from
#'        the feature to include in the plot. If provided, \code{downstream} must
#'        also be provided. When both are provided, \code{feature.gr} will be
#'        adjusted: ranges are centered at their midpoint, width is set to 1, and
#'        then extended by \code{upstream} and \code{downstream}. In this case,
#'        \code{zeroAt} is automatically calculated and any provided \code{zeroAt}
#'        value is ignored.
#' @param downstream An integer specifying the number of base pairs downstream
#'        from the feature to include in the plot. Must be provided together with
#'        \code{upstream}. See \code{upstream} for details.
#' @param zeroAt A numeric value specifying the zero point position within
#'        \code{feature.gr}. This determines where the x-axis zero is positioned.
#'        \itemize{
#'          \item If \code{0 <= zeroAt <= 1}: Treated as a fraction of the feature
#'                width (e.g., 0.5 means the center of the feature)
#'          \item If \code{zeroAt > 1}: Treated as an absolute position in base pairs
#'                from the start of the feature
#'        }
#'        Default is 0.5 (center of feature) when \code{upstream} and
#'        \code{downstream} are not provided. Ignored when \code{upstream} and
#'        \code{downstream} are provided.
#' @param n.tile An integer specifying the number of tiles/bins to divide each
#'        feature into. Default is 100. If \code{cvglists} contains pre-computed
#'        matrices from \code{featureAlignedSignal}, this must match the \code{n.tile}
#'        value used when creating those matrices.
#' @param ... Additional parameters passed to \code{\link[graphics]{matplot}} for
#'        customizing the plot appearance (e.g., \code{type}, \code{col}, \code{lty},
#'        \code{lwd}, \code{xlab}, \code{ylab}, \code{main}, etc.)
#' 
#' @return Returns invisibly a matrix with:
#'        \itemize{
#'          \item Rows: Tiles/bins (1 to \code{n.tile})
#'          \item Columns: Samples/conditions (one per element in \code{cvglists})
#'          \item Values: Average signal density across all features for each
#'                tile and sample
#'        }
#'        The matrix is also plotted as a line plot with appropriate x-axis labels.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates that all ranges in \code{feature.gr} have identical width
#'   \item If \code{upstream} and \code{downstream} are provided, adjusts
#'         \code{feature.gr} to be centered and extends it by the specified distances
#'   \item Calculates the zero point position for x-axis labeling
#'   \item If \code{cvglists} contains matrices, uses them directly; otherwise
#'         calls \code{featureAlignedSignal} to compute signal matrices
#'   \item Calculates mean signal density across all features for each tile
#'   \item Plots the results using \code{matplot} with automatic legend
#' }
#' 
#' \strong{X-axis labeling:}
#' The x-axis is labeled based on genomic coordinates relative to the zero point.
#' Labels are automatically generated using \code{grid.pretty} to create nice
#' round numbers. The zero point corresponds to the position specified by
#' \code{zeroAt} (or automatically calculated from \code{upstream} and
#' \code{downstream}).
#' 
#' \strong{NA handling:}
#' If the signal matrices contain NA values, they are omitted when calculating
#' column means (with a warning). This allows the function to work with sparse
#' or incomplete data.
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{featureAlignedSignal}} for computing signal matrices,
#'          \code{\link{featureAlignedHeatmap}} for heatmap visualization,
#'          \code{\link{featureAlignedExtendSignal}} for extending signal beyond
#'          feature boundaries
#' @keywords misc
#' @export
#' @importFrom BiocGenerics start end width strand
#' @importFrom graphics matplot axis legend
#' @importFrom grid grid.pretty
#' @examples
#' \dontrun{
#' ## Example 1: Using RleList objects
#' cvglists <- list(A = RleList(chr1 = Rle(sample.int(5000, 100), 
#'                                         sample.int(300, 100))), 
#'                  B = RleList(chr1 = Rle(sample.int(5000, 100), 
#'                                         sample.int(300, 100))))
#' feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width = 100))
#' featureAlignedDistribution(cvglists, feature.gr, zeroAt = 50, type = "l")
#' 
#' ## Example 2: Using upstream and downstream
#' feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width = 100))
#' featureAlignedDistribution(cvglists, feature.gr, 
#'                             upstream = 2000, downstream = 2000,
#'                             type = "l", col = c("red", "blue"))
#' 
#' ## Example 3: Using pre-computed matrices from featureAlignedSignal
#' signal_matrices <- featureAlignedSignal(cvglists, feature.gr, n.tile = 100)
#' featureAlignedDistribution(signal_matrices, feature.gr, 
#'                             zeroAt = 0.5, type = "l", lwd = 2)
#' }
#' 
featureAlignedDistribution <- function(cvglists, feature.gr, 
                                       upstream, downstream, 
                                       n.tile = 100L, zeroAt, ...) {
    stopifnot(inherits(feature.gr, "GRanges"))
    dots <- list(...)
    
    grWidr <- unique(width(feature.gr))
    if (missing(upstream) || missing(downstream)) {
        if (length(grWidr) != 1L) {
            stop("The width of feature.gr is not identical", call. = FALSE)
        }
        if (missing(zeroAt)) {
            zeroAt <- 0.5
            message("zero is set as the center of the feature.gr")
        }
    } else {
        if (!is.numeric(upstream) || !is.numeric(downstream)) {
            stop("'upstream' and 'downstream' must be numeric", call. = FALSE)
        }
        if (upstream < 0L || downstream < 0L) {
            stop("'upstream' and 'downstream' must be non-negative", 
                 call. = FALSE)
        }
        upstream <- as.integer(upstream)
        downstream <- as.integer(downstream)
        if (length(grWidr) != 1L || any(grWidr != 1L)) {
            start(feature.gr) <- start(feature.gr) + 
                floor(width(feature.gr) / 2L)
            width(feature.gr) <- 1L
            warning("feature.gr is set to the center of feature.gr", 
                   call. = FALSE)
        }
        if (!missing(zeroAt)) {
            warning("'zeroAt' will be ignored when 'upstream' and ",
                    "'downstream' are provided", call. = FALSE)
        }
        zeroAt <- upstream / (upstream + downstream)
        end(feature.gr) <- start(feature.gr) + downstream
        start(feature.gr) <- start(feature.gr) - upstream
        grWidr <- unique(width(feature.gr))
    }
    stopifnot(is.numeric(zeroAt))
    stopifnot(zeroAt >= 0L)
    if (zeroAt <= 1L) {
        zero <- round(grWidr * zeroAt)
    } else {
        zero <- round(zeroAt)
    }
    
    grWid <- c(0, grWidr) - zero
    grWidLab <- grid.pretty(grWid)
    grWidAt <- (grWidLab + zero) / grWidr * n.tile
    if (inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))) {
        cvglistsName <- substitute(deparse(cvglists))
        cvglists <- list(cvglists)
        names(cvglists) <- cvglistsName
    }
    if (!is.list(cvglists)) {
        stop("'cvglists' must be output of featureAlignedSignal or ",
             "a list of SimpleRleList or RleList", call. = FALSE)
    }
    
    cls <- vapply(cvglists, is.matrix, FUN.VALUE = logical(1))
    if (all(cls)) {
        cov <- cvglists
        if (ncol(cov[[1L]]) != n.tile) {
            stop("'n.tile' must match the value used in featureAlignedSignal", 
                 call. = FALSE)
        }
    } else {
        cls <- vapply(cvglists, inherits,
                      what = c("SimpleRleList", "RleList", "CompressedRleList"))
        if (any(!cls)) {
            stop("'cvglists' must be a list of SimpleRleList or RleList",
                 call. = FALSE)
        }
        cov <- featureAlignedSignal(cvglists, feature.gr, n.tile = n.tile)
    }
    
    ## normalized read density
    if (any(sapply(cov, function(.ele) any(is.na(.ele))))) {
        warning("cvglists contain NA values. ", 
                "NA value will be omit.")
    }
    density <- sapply(cov, colMeans, na.rm = TRUE)
    try({
        matplot(density, ..., xaxt = "n")
        axis(1, at = grWidAt, labels = grWidLab)
        lty <- if (!is.null(dots$lty)) dots$lty else 1:5
        lwd <- if (!is.null(dots$lwd)) dots$lwd else 1
        col <- if (!is.null(dots$col)) dots$col else 1:6
        legend("topright", legend = colnames(density), col = col,
               lty = lty, lwd = lwd)
    })
    
    return(invisible(density))
}
