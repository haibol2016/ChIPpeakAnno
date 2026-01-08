#' Plot heatmap of signals across features
#' 
#' @description 
#' Creates a heatmap visualization of signal intensity across genomic features.
#' Each row represents a feature, and columns represent bins/tiles along the
#' feature. This function is useful for visualizing ChIP-seq or RNA-seq signal
#' patterns across many features simultaneously, allowing identification of
#' common patterns and outliers.
#' 
#' The function supports multiple samples (displayed as adjacent heatmap columns),
#' annotation tracks (displayed on the right side), and flexible sorting options.
#' Features can be sorted by metadata columns or by signal intensity. The heatmap
#' uses a color scale to represent signal intensity, with customizable color
#' schemes and value ranges.
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
#'        List names will be used as labels below each heatmap column.
#' @param feature.gr An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#'        with identical width for all ranges. The features around which signal
#'        will be plotted (e.g., TSS, gene centers, peak centers). Each range
#'        becomes one row in the heatmap.
#' @param upstream An integer specifying the number of base pairs upstream from
#'        the feature to include in the plot. If provided, \code{downstream} must
#'        also be provided. When both are provided, \code{feature.gr} will be
#'        adjusted: ranges are centered at their midpoint, width is set to 1, and
#'        then extended using \code{promoters()} by \code{upstream} and
#'        \code{downstream + 1}. In this case, \code{zeroAt} is automatically
#'        calculated and any provided \code{zeroAt} value is ignored.
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
#' @param annoMcols A character vector specifying column names from
#'        \code{mcols(feature.gr)} to display as annotation tracks on the right
#'        side of the heatmap. Each specified column will be displayed as a
#'        vertical bar. For numeric columns, a color gradient is used. For
#'        character/factor columns, distinct colors are assigned to each unique
#'        value. If a column contains valid color names or numeric values 1-8
#'        (R color palette indices), those are used directly.
#' @param sortBy A character vector specifying how to sort features. Features are
#'        sorted by the specified columns in order. Can include:
#'        \itemize{
#'          \item Column names from \code{mcols(feature.gr)}: Sort by metadata
#'          \item Names from \code{cvglists}: Sort by total signal intensity
#'                (sum across all tiles) for that sample
#'        }
#'        If \code{annoMcols} is provided, sorting first uses \code{annoMcols}
#'        columns, then the \code{sortBy} columns. Default is the first sample
#'        name. Set to \code{NULL} or empty vector to disable sorting (features
#'        will be displayed in their original order).
#' @param color A character vector of colors used for the heatmap color scale.
#'        Colors should be ordered from low to high intensity. Default is
#'        \code{colorRampPalette(c("yellow", "red"))(50)}, creating a yellow-to-red
#'        gradient with 50 levels.
#' @param lower.extreme A numeric vector specifying the lower boundary value for
#'        each sample. If a single value is provided, it is used for all samples.
#'        If missing, the minimum value in each sample's data is used. Values
#'        below this threshold are mapped to the first color in \code{color}.
#' @param upper.extreme A numeric vector specifying the upper boundary value for
#'        each sample. If a single value is provided, it is used for all samples.
#'        If missing, the maximum value in each sample's data is used. Values
#'        above this threshold are mapped to the last color in \code{color}.
#' @param margin A numeric vector of length 4 specifying the plot margins:
#'        \code{c(bottom, left, top, right)}. Default is \code{c(0.1, 0.01, 0.15, 0.1)}.
#'        Values are in normalized device coordinates (0-1).
#' @param gap A numeric value specifying the gap between heatmap columns (for
#'        multiple samples) and between annotation tracks. Default is 0.01 in
#'        normalized device coordinates.
#' @param newpage A logical value. If \code{TRUE} (default), calls
#'        \code{grid.newpage()} before drawing, creating a new plot. If
#'        \code{FALSE}, draws on the current graphics device.
#' @param gp A \code{\link[grid]{gpar}} object specifying graphical parameters
#'        for text labels (font size, family, etc.). Default is
#'        \code{gpar(fontsize=10)}.
#' @param ... Additional parameters (currently not used)
#' 
#' @return Returns invisibly a \code{\link[grid]{gList}} object containing all
#'        the graphical elements of the heatmap. This can be saved and redrawn
#'        later using \code{grid.draw()}, or modified using grid editing
#'        functions.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates inputs and processes \code{feature.gr} (adjusts for
#'         \code{upstream}/\code{downstream} if provided)
#'   \item If \code{cvglists} contains RleList objects, calls
#'         \code{featureAlignedSignal} to compute signal matrices
#'   \item Sorts features if \code{sortBy} is specified
#'   \item Converts signal values to colors based on \code{color},
#'         \code{lower.extreme}, and \code{upper.extreme}
#'   \item Creates annotation tracks if \code{annoMcols} is specified
#'   \item Draws heatmap columns for each sample
#'   \item Adds x-axis labels showing genomic coordinates
#'   \item Adds color scale legend(s)
#'   \item Draws all elements using grid graphics
#' }
#' 
#' \strong{Sorting behavior:}
#' When \code{sortBy} is specified, features are sorted in descending order
#' (highest values first) by:
#' \enumerate{
#'   \item First, by \code{annoMcols} columns (if provided), in the order
#'         specified
#'   \item Then, by \code{sortBy} columns, in the order specified
#' }
#' For character/factor columns, levels are reversed so that the last unique
#' value appears first. For numeric columns, higher values appear first.
#' 
#' \strong{Color mapping:}
#' Signal values are mapped to colors using \code{cut()} with breaks determined
#' by \code{lower.extreme}, \code{upper.extreme}, and the length of \code{color}.
#' Values are clamped to the specified range before mapping. If all samples use
#' the same range (single \code{lower.extreme} and \code{upper.extreme}), a
#' single legend is displayed. Otherwise, each sample gets its own legend.
#' 
#' \strong{Annotation tracks:}
#' Annotation tracks are displayed as vertical bars on the right side of the
#' heatmap. Each track corresponds to one metadata column. The function
#' automatically:
#' \itemize{
#'   \item Detects if values are valid colors (color names or R palette indices)
#'   \item For numeric columns: Creates a color gradient using predefined color
#'         groups
#'   \item For character/factor columns: Assigns distinct colors to each unique
#'         value
#'   \item Adds a legend for each annotation track
#' }
#' 
#' \strong{X-axis labeling:}
#' The x-axis is labeled based on genomic coordinates relative to the zero point.
#' Labels are automatically generated using \code{grid.pretty} to create nice
#' round numbers. For multiple samples, labels are only shown on the first and
#' last columns to avoid clutter.
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{featureAlignedSignal}} for computing signal matrices,
#'          \code{\link{featureAlignedDistribution}} for average signal plots,
#'          \code{\link{featureAlignedExtendSignal}} for extending signal beyond
#'          feature boundaries
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom BiocGenerics start end width strand
#' @importFrom S4Vectors mcols
#' @importFrom grid grid.newpage viewport legendGrob gpar gList gEdit rasterGrob
#' gTree rasterGrob grid.pretty yaxisGrob grid.draw xaxisGrob textGrob
#' @importFrom grDevices colorRampPalette col2rgb rgb
#' @examples
#' \dontrun{
#' ## Example 1: Basic heatmap with annotation
#' cvglists <- list(A = RleList(chr1 = Rle(sample.int(5000, 100), 
#'                                         sample.int(300, 100))), 
#'                  B = RleList(chr1 = Rle(sample.int(5000, 100), 
#'                                         sample.int(300, 100))))
#' feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width = 100))
#' feature.gr$anno <- rep(c("type1", "type2"), c(25, 24))
#' featureAlignedHeatmap(cvglists, feature.gr, zeroAt = 50, 
#'                       annoMcols = "anno")
#' 
#' ## Example 2: Sort by signal intensity
#' featureAlignedHeatmap(cvglists, feature.gr, zeroAt = 50,
#'                       sortBy = "A")
#' 
#' ## Example 3: Custom colors and value range
#' featureAlignedHeatmap(cvglists, feature.gr, zeroAt = 50,
#'                       color = colorRampPalette(c("blue", "white", "red"))(100),
#'                       lower.extreme = 0, upper.extreme = 5000)
#' 
#' ## Example 4: Using upstream and downstream
#' feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width = 100))
#' featureAlignedHeatmap(cvglists, feature.gr, 
#'                       upstream = 2000, downstream = 2000,
#'                       annoMcols = "anno")
#' 
#' ## Example 5: Multiple annotation tracks
#' feature.gr$score <- runif(length(feature.gr), 0, 100)
#' featureAlignedHeatmap(cvglists, feature.gr, zeroAt = 50,
#'                       annoMcols = c("anno", "score"),
#'                       sortBy = c("anno", "A"))
#' }
#' 
featureAlignedHeatmap <- 
    function(cvglists, feature.gr, upstream, downstream, 
             zeroAt, n.tile = 100,
             annoMcols = c(), sortBy = names(cvglists)[1],
             color = colorRampPalette(c("yellow", "red"))(50),
             lower.extreme, upper.extreme,
             margin = c(0.1, 0.01, 0.15, 0.1), gap = 0.01, 
             newpage = TRUE, gp = gpar(fontsize = 10),
             ...) {
    if (!missing(lower.extreme)) {
        stopifnot(is.numeric(lower.extreme))
    }
    if (!missing(upper.extreme)) {
        stopifnot(is.numeric(upper.extreme))
    }
    stopifnot(inherits(gp, "gpar"))
    stopifnot(is.logical(newpage))
    stopifnot(inherits(feature.gr, "GRanges"))
    
    grWidr <- unique(width(feature.gr))
    if (missing(upstream) || missing(downstream)) {
        if (length(grWidr) != 1) {
            stop("The width of feature.gr is not identical.")
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
        if (length(grWidr) != 1) {
            start(feature.gr) <- start(feature.gr) + floor(grWidr / 2)
            width(feature.gr) <- 1
            warning("feature.gr is set to the center of feature.gr",
                    call. = FALSE)
        }
        zeroAt <- upstream / (upstream + downstream)
        feature.gr <- promoters(feature.gr, upstream = upstream,
                                downstream = downstream + 1)
        grWidr <- unique(width(feature.gr))
    }
    stopifnot(is.numeric(zeroAt))
    stopifnot(zeroAt >= 0)
    if (zeroAt <= 1) {
        zero <- round(grWidr * zeroAt)
    } else {
        zero <- round(zeroAt)
    }
    
    grWid <- c(0, grWidr) - zero
    grWidLab <- grid.pretty(grWid)
    grWidAt <- (grWidLab + zero) / grWidr
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
                      what = c("SimpleRleList", "RleList", "CompressedRleList"),
                      FUN.VALUE = logical(1))
        if (any(!cls)) {
            stop("'cvglists' must be a list of SimpleRleList or RleList",
                 call. = FALSE)
        }
        cov <- featureAlignedSignal(cvglists, feature.gr, n.tile = n.tile)
    }
    
    annoMcols <- annoMcols[annoMcols %in% colnames(mcols(feature.gr))]
    if (length(sortBy) > 0) {
        covTotalCoverage <- 
            do.call(cbind, lapply(cov, rowSums, na.rm = TRUE))
        colnames(covTotalCoverage) <- names(cvglists)
        if (length(annoMcols) > 0) {
            sortMatrix <- mcols(feature.gr)[, annoMcols, drop = FALSE]
            for (i in seq.int(ncol(sortMatrix))) {
                if (!is.numeric(sortMatrix[, i])) {
                    sortMatrix[, i] <- 
                        factor(as.character(sortMatrix[, i]), 
                               levels =
                                 rev(unique(as.character(sortMatrix[, i]))))
                }
            }
            if (all(sortBy %in% colnames(mcols(feature.gr)))) {
                sortMatrix <- 
                    cbind(sortMatrix, 
                          as.data.frame(mcols(feature.gr)[, sortBy, 
                                                          drop = FALSE]))
            } else {
                if (all(sortBy %in% colnames(covTotalCoverage))) {
                    sortMatrix <- 
                        cbind(sortMatrix,
                              data.frame(covTotalCoverage[, sortBy,
                                                          drop = FALSE]))
                } else {
                    stop("'sortBy' is incorrect", call. = FALSE)
                }
            }
        } else {
            if (all(sortBy %in% colnames(mcols(feature.gr)))) {
                sortMatrix <- 
                    as.data.frame(mcols(feature.gr)[, sortBy, drop = FALSE])
            } else {
                if (all(sortBy %in% colnames(covTotalCoverage))) {
                    sortMatrix <- data.frame(covTotalCoverage[, sortBy,
                                                              drop = FALSE])
                } else {
                    stop("'sortBy' is incorrect", call. = FALSE)
                }
            }
        }
        names(sortMatrix) <- paste("V", seq_len(ncol(sortMatrix)))
        soid <- do.call(order, c(as.list(sortMatrix), decreasing = TRUE))
        cov <- lapply(cov, function(.ele) .ele[soid, ])
        feature.gr <- feature.gr[soid]
    }
    
    #covert to color
    if (missing(lower.extreme)) {
        lower.extreme <- rep(0, length(cvglists))
    }
    if (missing(upper.extreme)) {
        lim <- list()
    } else {
        upper.extreme <- 
            rep(upper.extreme, length(cvglists))[seq_along(cvglists)]
        lim <- as.list(data.frame(rbind(lower.extreme, upper.extreme)))
        names(lim) <- NULL
    }
    
    if (length(lim) == 0) {
        lim <- sapply(cov, range, na.rm = TRUE, simplify = FALSE)
    }
    if (length(lim) == 1) {
        lim <- rep(lim, length(cov))
        legend <- 1
    } else {
        legend <- length(cov)
    }
    if (any(sapply(lim, length) < 2L)) {
        stop("each limit should be a numeric vector with length 2", 
             call. = FALSE)
    }
    cov <- mapply(function(.ele, .lim) {
        if (min(.lim) == max(.lim)) {
            .lim <- seq(min(.lim), max(.lim) + 1, length.out = length(color) + 1)
        } else {
            .lim <- seq(min(.lim), max(.lim), length.out = length(color) + 1)
        }
        .lim[1] <- min(.lim[1], 0)
        .lim[length(.lim)] <- max(max(.lim), .Machine$integer.max)
        t(apply(.ele, 1, 
                cut, breaks = .lim, labels = color, include.lowest = TRUE))
    }, cov, lim, SIMPLIFY = FALSE)
    # library(grid)
    # grid.raster
    if (newpage) grid.newpage()
    allGrob <- NULL
    height <- 1 - margin[3] - margin[1]
    width <- 1 - margin[4] - margin[2]
    ## draw y labels
    isFloat <- function(x) {
        if (is.numeric(x)) {
            return(x != floor(x))
        }
        return(FALSE)
    }
    areColors <- function(x) {
        sapply(x, function(X) {
            if (isFloat(X)) return(FALSE)
            if (is.numeric(X)) {
                if (X < 9 && X > 0) {
                    return(TRUE)
                } else {
                    return(FALSE)
                }
            }
            tryCatch(is.matrix(col2rgb(X)), 
                     error = function(e) FALSE)
        })
    }
    colorGroup <- 
        list(c("#295AA5", "#52BDBD", "#FFAD6B", "#F784A5", "#CEB55A", 
               "#31B5D6", "#BD2119", "#F79463", "#FFE608", "#D69419", 
               "#29B5CE", "#737BB5"), 
             c("#4A4221", "#A5BD42", "#8CC584", "#CEE6D6", "#C5D64A", 
               "#FFF79C", "#84AD31", "#DBDECE", "#D6AD4A", "#4A4229", 
               "#637B7B", "#525229", "#F7F7AD", "#EFAD63"), 
             c("#6B5252", "#DE6B73", "#000000", "#4A4221", "#3AA55A", 
               "#DEE6C5", "#102119", "#DEE6CE", "#738CC5", "#31317B", 
               "#F7F7AD", "#8C9CA5", "#732163", "#001000"), 
             c("#192919", "#8C946B", "#6B1921", "#F7F7AD", "#4A4229", 
               "#CEA54A", "#000000", "#6B736B", "#524A29", "#423119", 
               "#213129", "#7B4229", "#7B2121", "#C57B4A"), 
             c("#9D2932", "#EFDEB0", "#789262", "#494166", "#A88462", 
               "#D9B612", "#177CB0", "#F3F8F1", "#424B50", "#B36D61", 
               "#C3272B", "#A98175", "#D4F2E8", "#FF3300", "#76664D", 
               "#ED5736", "#A3E2C5", "#415065", "#D7ECF1"), 
             c("#009C8C", "#0094AD", "#BDDEC5", "#00ADCE", "#DEFF19", 
               "#D6EFB5", "#EFF7E6", "#42CE73", "#42A5CE", "#BDFFF7", 
               "#DEFFDE", "#4A6BCE"), 
             c("#D6B5BD", "#EFD6BD", "#FFE6D6", "#DEE6CE", "#F7F7C5", 
               "#D6E6B5", "#73BD8C", "#BDDEDE", "#DEBD8C", "#E6A5BD", 
               "#F7ADAD", "#FFE6D6", "#F7F7D6", "#E6BDBD"), 
             c("#E6E6BD", "#637342", "#425A52", "#298473", "#BDDEC5", 
               "#527384", "#73847B", "#CED6AD", "#94D6CE", "#527384", 
               "#4A6B5A", "#318442", "#C5DEDE"), 
             c("#F7846B", "#F79C9C", "#CEE6E6", "#FFEFBD", "#EFEFBD", 
               "#FFE6DE", "#FFBDA5", "#FFEF00", "#63BD84", "#D6E6AD", 
               "#FFF7D6", "#FFFFEF", "#FFE6DE", "#FFBD9C", "#EF9C94"), 
             c("#4A4229", "#212108", "#6B8484", "#EFE6BD", "#E68C3A", 
               "#A56352", "#6B2921", "#525229", "#844229", "#8C524A", 
               "#E68C3A", "#D6AD4A", "#212108", "#524A3A"), 
             c("#312184", "#AD2919", "#EF3129", "#D69C4A", "#FF6B00", 
               "#6B4A31", "#CE3121", "#AD2919", "#5A3A5A", "#292931", 
               "#3A3184", "#9C5229", "#EFCE19"), 
             c("#6373B5", "#FF5200", "#4AF54A", "#FFCE10", "#EF3A21", 
               "#EF4A94", "#5A63AD", "#FFEDA5", "#EF3121", "#FFC510", 
               "#FF5200"), 
             c("#FFB6C1", "#DC143C", "#DB7093", "#FF69B4", "#FF1493", 
               "#C71585", "#DA70D6", "#FF00FF", "#8B008B", "#9400D3", 
               "#4B0082", "#FFA07A", "#FF7F50", "#FF4500", "#E9967A", 
               "#FF6347", "#F08080", "#CD5C5C", "#FF0000", "#A52A2A", 
               "#B22222", "#8B0000"), 
             c("#7B68EE", "#6A5ACD", "#483D8B", "#E6E6FA", "#0000FF", 
               "#191970", "#4169E1", "#6495ED", "#B0C4DE", "#1E90FF", 
               "#87CEFA", "#00BFFF", "#00FFFF"), 
             c("#008B8B", "#48D1CC", "#7FFFD4", "#66CDAA", "#00FA9A", 
               "#00FF7F", "#3CB371", "#2E8B57", "#90EE90", "#32CD32", 
               "#00FF00", "#008000", "#006400", "#ADFF2F", "#556B2F", 
               "#9ACD32", "#6B8E23"), 
             c("#FAFAD2", "#FFFF00", "#808000", "#BDB76B", "#FFFACD", 
               "#EEE8AA", "#F0E68C", "#FFD700", "#DAA520", "#B8860B", 
               "#F5DEB3", "#FFA500", "#FFEBCD", "#FFDEAD", "#D2B48C", 
               "#FF8C00", "#CD853F", "#8B4513"),
             c("#F5F5F5", "#DCDCDC", "#D3D3D3", "#C0C0C0", "#A9A9A9", 
               "#808080", "#696969", "#5C5C5C", "#4D4D4D", "#333333", 
               "#1A1A1A", "#000000"))
    if(length(annoMcols)>0){
        ci <- 1
        ytop <- 0.75-margin[3]
        mcols <- mcols(feature.gr)
        vp.annoMcols <- viewport(x=.5, y=height/2+margin[1], 
                                 width=1, height=height,
                                 name="vp.annoMcols")
        for(i in annoMcols){
            mc <- as.character(mcols[, i])
            mc.name <- ifelse(is.numeric(i), colnames(mcols)[i], i)
            if(!all(areColors(mcols[, i]))){
                if(is.numeric(mcols[, i])){
                    mc.range <- range(mcols[, i])
                    mc.color <- range(colorGroup[[ci]])
                    mc.color <- colorRampPalette(mc.color)(100)
                    mc <- cut(mcols[, i], breaks = 100, labels = mc.color)
                    mc <- as.character(mc)
                    mc.label <- sort(grid.pretty(mc.range), decreasing = TRUE)
                    mc.label.color <-
                        mc.color[round(100*(mc.label-min(mc.range))/
                                     diff(mc.range))]
                    mc.label <- mc.label[!is.na(mc.label.color)]
                    mc.label.color <- mc.label.color[!is.na(mc.label.color)]
                }else{
                    mc.label <- unique(mc)
                    mc <- factor(mc, levels=mc.label)
                    mc.label.color <- rep(colorGroup[[ci]], 
                                          length(mc.label))[seq_along(mc.label)]
                    levels(mc) <- mc.label.color
                    mc <- as.character(mc)
                }
                
                ci <- ci+1
            }else{
                mc.label <- unique(mcols[, i])
                if(is.numeric(mcols[, i])){
                    mc <- col2rgb(mc);
                    mc <- rgb(red=mc[1, ],
                              green=mc[2, ],
                              blue=mc[3, ], maxColorValue = 255)
                }
                mc.label.color <- unique(mc)
            }
            ## mc is color
            ## mc.label is legends
            ## plot legned
            ht.annoMcols.legend <- length(mc.label) * 0.03
            wd.annoMcols.legend <- margin[4]
            vp.annoMcols.legend <- 
                viewport(x=1-wd.annoMcols.legend+gap, 
                         y=ytop-ht.annoMcols.legend/2, 
                         width=wd.annoMcols.legend, 
                         height=ht.annoMcols.legend,
                         just=0,
                         name="vp.annoMcols.legend")
            annoMcols.legend <- legendGrob(labels=mc.label, 
                                           ncol=1,
                                           pch=15,
                                           gp=gpar(col=mc.label.color))
            allGrob <- gList(allGrob, 
                             gTree(children=gList(annoMcols.legend),
                                   vp=vp.annoMcols.legend,
                                   name=paste("gTree.annoMcols.legend",
                                              mc.name,
                                              sep=".")))
            ytop <- ytop-ht.annoMcols.legend-gap
            
            ## add annoMcols
            vp.annoMcols.sub <- viewport(x=width - .0075,
                                         y=.5, 
                                         width=.015,
                                         height=1,
                                         name=paste("vp.annoMcols",
                                                    mc.name,
                                                    sep="."))
            raster.annoMcols.sub <- 
                rasterGrob(mc,
                           width=1, height=1,
                           interpolate=FALSE,
                           name=paste("raster.annoMcols",
                                      mc.name,
                                      sep="."),
                           vp=vp.annoMcols.sub)
            xaxis.annoMcols.sub <- 
                xaxisGrob(at=.5,
                          label=mc.name,
                          gp=gp,
                          edits=gEdit("labels", rot=90),
                          name=paste("xaxis.annoMcols",
                                     mc.name,
                                     sep="."),
                          vp=vp.annoMcols.sub)
            allGrob <- gList(allGrob, 
                             gTree(children=gList(raster.annoMcols.sub,
                                                  xaxis.annoMcols.sub),
                                   vp=vp.annoMcols,
                                   name=paste("gTree.annoMcols",
                                              mc.name,
                                              sep=".")))
            width <- width - .015
        }
        width <- width - gap
    }
    
    vp.heatmap <- viewport(x=width/2+margin[2], y=height/2+margin[1], 
                           width=width, height=height,
                           name="vp.heatmap")
    l <- length(cov)
    wid <- 1/l
    x <- -wid/2
    y <- 1
    for(i in seq_len(l)){
        x <- x+wid
        x.name <- names(cvglists)[i]
        vp.heatmap.sub <- viewport(x=x, width=wid-gap, height=y,
                                   name=paste("vp.heatmap",
                                              x.name,
                                              sep="."))
        raster.heatmap.sub <- rasterGrob(cov[[i]], 
                                         width=1, height=1, 
                                         interpolate=FALSE,
                                         name=paste("raster.heatmap",
                                                    x.name,
                                                    sep="."),
                                         vp=vp.heatmap.sub)
        
        grWidLab.this <- grWidLab
        if(i==1){
            if(i!=l){
                if(length(grWidAt)>1){
                    grWidLab.this[length(grWidLab)] <- " "
                }
            }
        }else{
            if(i==l){##i!=1
                if(length(grWidAt)>1){
                    grWidLab.this[1] <- " "
                }
            }else{## middle
                if(length(grWidAt)>2){
                    grWidLab.this[1] <- " "
                    grWidLab.this[length(grWidLab)] <- " "
                }
            }
        }
        xaxis.heatmap.sub <- xaxisGrob(at=grWidAt,
                                       label=grWidLab.this,
                                       gp=gp,
                                       edits=gEdit("labels", rot=90),
                                       name=paste("xaxis.heatmap",
                                                  x.name,
                                                  sep="."),
                                       vp=vp.heatmap.sub)
        allGrob <- gList(allGrob, 
                         gTree(children=gList(raster.heatmap.sub,
                                              xaxis.heatmap.sub),
                               vp=vp.heatmap,
                               name=paste("gTree.heatmap",
                                          x.name,
                                          sep=".")))
    }
    
    ## add heatmap legend
    if(legend==1){
        vp.legend <- viewport(x=1.01-margin[4]+gap, y=.875-margin[3], 
                              width=.02, height=.25,
                              name="vp.legend")
        raster.legend <- rasterGrob(matrix(rev(color), ncol=1), 
                                    width=1, height=1,
                                    name="raster.legend")
        .lim <- lim[[1]]
        label <- grid.pretty(.lim)
        ll <- length(label)
        if(ll>=2){
            if(ll%%2==0){
                label <- label[c(1, ll)]
            }else{
                label <- label[c(1, ceiling(ll/2), ll)]
            }
            at <- (label-min(.lim))/abs(diff(.lim))
            yaxis.legend <- yaxisGrob(at=at,
                                      label=label,
                                      main=FALSE,
                                      gp=gp,
                                      name="yaxis.legend")
        }
        allGrob <- gList(allGrob, 
                         gTree(children=gList(raster.legend,
                                              yaxis.legend),
                               vp=vp.legend,
                               name="gTree.legend"))
        
        ## draw x labels
        vp.xlabels <- viewport(x=width/2+margin[2], 
                               y=1-margin[3]/2+gap, 
                               width=width, height=margin[3]-gap,
                               name="vp.xlabels")
        x <- -wid/2
        for(i in seq_len(l)){
            x <- x+wid
            x.name <- names(cvglists)[i]
            text.xlabels.sub <- 
                textGrob(label = x.name,
                         gp=gp,
                         just=1, hjust=0,
                         vp=viewport(x=x, y=.5, width=wid, height=1, 
                                     name=paste("vp.xlabels", x.name, sep=".")),
                         name=paste("text.xlabels", x.name, sep="."))
            allGrob <- gList(allGrob, 
                             gTree(children=gList(text.xlabels.sub),
                                   vp=vp.xlabels,
                                   name=paste("gTree.xlabels",
                                              x.name,
                                              sep=".")))
        }
    }else{
        vp.legend <- viewport(x=width/2+margin[2], 
                              y=1-margin[3]/3+gap, 
                              width=width, height=margin[3]*2/3-gap,
                              name="vp.legend")
        
        vp.xlabels <- viewport(x=width/2+margin[2], 
                               y=1-margin[3]*5/6+gap, 
                               width=width, height=margin[3]/3-gap,
                               name="vp.xlabels")
        
        x <- -wid/2
        for(i in seq_len(l)){
            x <- x+wid
            x.name=names(cvglists)[i]
            vp.legend.sub <- 
                viewport(x=x, y=0.25, width=wid/2, height=0.2,
                         name=paste("vp.legend",
                                    x.name,
                                    sep="."))
            raster.legend.sub <- 
                rasterGrob(matrix(color, nrow=1), 
                           width=1, height=1,
                           name=paste("raster.legend", x.name, sep="."),
                           vp=vp.legend.sub)
            .lim <- lim[[i]]
            label <- grid.pretty(.lim)
            ll <- length(label)
            xaxis.legned.sub <- NULL
            if(ll>=2){
                if(ll%%2==0){
                    label <- label[c(1, ll)]
                }else{
                    label <- label[c(1, ceiling(ll/2), ll)]
                }
                at <- (label-min(.lim))/abs(diff(.lim))
                xaxis.legned.sub <- xaxisGrob(at=at,
                                              label=label,
                                              main=FALSE,
                                              gp=gp,
                                              name=paste("xaxis.legend",
                                                         x.name,
                                                         sep="."),
                                              vp=vp.legend.sub)
            }
            allGrob <- gList(allGrob, 
                             gTree(children=gList(raster.legend.sub,
                                                  xaxis.legned.sub),
                                   vp=vp.legend,
                                   name=paste("gTree.legend",
                                              x.name,
                                              sep=".")))
            
            ## draw x labels
            text.xlabels.sub <- 
                textGrob(label = x.name,
                         gp=gp,
                         vp=viewport(x=x, y=.5, width=wid, height=1, 
                                     name=paste("vp.xlabels", x.name, sep=".")),
                         name=paste("text.xlabels", x.name, sep="."))
            allGrob <- gList(allGrob, 
                             gTree(children=gList(text.xlabels.sub),
                                   vp=vp.xlabels,
                                   name=paste("gTree.xlabels",
                                              x.name,
                                              sep=".")))
            
        }
    }    
    
    ## consider how to return all and could be redraw by grid.draw.
    grid.draw(allGrob)
    return(invisible(allGrob))
}
