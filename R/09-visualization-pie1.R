#' Draw pie chart with percentages
#' 
#' @description 
#' Creates a pie chart with optional percentage or raw number labels for each
#' slice. This function extends the base \code{\link[graphics]{pie}} function
#' with additional features such as automatic percentage calculation, cutoff
#' thresholds for label display, flexible legend positioning, and customizable
#' label placement. Useful for visualizing categorical data distributions with
#' clear percentage or count annotations.
#' 
#' @param x A vector of non-negative numerical quantities. The values in \code{x}
#'        are displayed as the areas of pie slices. Values are normalized to
#'        percentages automatically. Must not contain NA or negative values.
#' @param labels One or more expressions or character strings giving names for
#'        the slices. Default is \code{names(x)}. Other objects are coerced by
#'        \code{as.graphicsAnnot}. For empty or NA (after coercion to character)
#'        labels, no label nor pointing line is drawn. If \code{legend = TRUE},
#'        labels are shown in the legend instead of on the pie chart.
#' @param edges An integer specifying the number of edges used to approximate
#'        the circular outline of the pie. Higher values create smoother circles
#'        but may be slower. Default is 200. The actual number of edges per slice
#'        is proportional to the slice size.
#' @param radius A numeric value between 0 and 1 specifying the radius of the
#'        pie. The pie is drawn centered in a square box whose sides range from
#'        -1 to 1. Default is 0.8. If the character strings labeling the slices
#'        are long, it may be necessary to use a smaller radius.
#' @param clockwise A logical value indicating if slices are drawn clockwise
#'        (\code{TRUE}) or counter-clockwise (\code{FALSE}, default). When
#'        \code{clockwise = TRUE}, \code{init.angle} defaults to 90 degrees.
#' @param init.angle A numeric value specifying the starting angle (in degrees)
#'        for the slices. Defaults to 0 (i.e., "3 o'clock") when
#'        \code{clockwise = FALSE}, or 90 (i.e., "12 o'clock") when
#'        \code{clockwise = TRUE}.
#' @param density The density of shading lines, in lines per inch. The default
#'        value of \code{NULL} means that no shading lines are drawn. Non-positive
#'        values also inhibit the drawing of shading lines. If specified, can be
#'        a vector with one value per slice.
#' @param angle The slope of shading lines, given as an angle in degrees
#'        (counter-clockwise). Default is 45. If specified, can be a vector with
#'        one value per slice.
#' @param col A vector of colors to be used in filling or shading the slices.
#'        If \code{NULL} (default), a set of 7 pastel colors is used:
#'        \code{c("white", "lightblue", "mistyrose", "lightcyan", "lavender",
#'        "cornsilk", "pink")}. If \code{density} is specified, \code{par("fg")}
#'        is used instead. Colors are recycled if the vector is shorter than the
#'        number of slices.
#' @param border A vector of colors for the border of each slice. If \code{NULL},
#'        no border is drawn. Can be a single value or a vector (recycled as needed).
#' @param lty A vector of line types for slice borders. If \code{NULL}, uses default.
#'        Can be a single value or a vector (recycled as needed).
#' @param main A character string giving an overall title for the plot. Can be
#'        \code{NULL} for no title.
#' @param percentage A logical value. If \code{TRUE} (default), percentage
#'        labels are displayed on each slice. Percentages are calculated as
#'        \code{100 * (slice_value / sum(x))}. Only slices with percentage >=
#'        \code{cutoff} are labeled.
#' @param rawNumber A logical value. If \code{TRUE}, raw numbers (from \code{x})
#'        are displayed on each slice instead of percentages. Default is
#'        \code{FALSE}. If both \code{percentage} and \code{rawNumber} are
#'        \code{TRUE}, only percentages are shown. Only slices with percentage >=
#'        \code{cutoff} are labeled.
#' @param digits An integer specifying the number of significant digits to use
#'        for percentage formatting. Default is 3. See \code{\link[base]{formatC}}
#'        for details. Only used when \code{percentage = TRUE}.
#' @param cutoff A numeric value between 0 and 1 specifying the minimum
#'        percentage threshold for displaying labels. Slices with percentage <
#'        \code{cutoff} will not show percentage or raw number labels. Default is
#'        0.01 (1\%). This prevents cluttering the chart with labels for very small
#'        slices.
#' @param legend A logical value. If \code{TRUE}, a legend is drawn instead of
#'        labels with pointing lines on the pie chart. Default is \code{FALSE}.
#'        When \code{legend = TRUE}, slice labels are shown in the legend, and
#'        no pointing lines are drawn.
#' @param legendpos A character string or numeric vector specifying the legend
#'        position. See \code{\link[graphics]{legend}} for valid values (e.g.,
#'        "topright", "bottomleft", c(x, y) coordinates). Default is "topright".
#'        Only used when \code{legend = TRUE}.
#' @param legendcol An integer specifying the number of columns for the legend.
#'        Default is 2. See \code{\link[graphics]{legend}} for details. Only used
#'        when \code{legend = TRUE}.
#' @param radius.innerlabel A numeric value between 0 and 1 specifying the
#'        position of percentage or raw number labels relative to the circle center.
#'        Values between 0 (center) and 1 (edge). Default equals \code{radius},
#'        placing labels at the edge of the pie. Smaller values place labels
#'        closer to the center.
#' @param ... Additional graphical parameters passed to \code{text()} for label
#'        formatting and \code{title()} for the main title. Common parameters
#'        include \code{cex}, \code{font}, \code{col}, etc.
#' 
#' @return Returns \code{invisible(NULL)}. The function is called for its side
#'        effect of drawing a pie chart.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates that \code{x} contains only non-negative numeric values
#'   \item Normalizes values to cumulative proportions (0 to 1)
#'   \item Calculates slice angles based on proportions
#'   \item Draws each slice as a polygon with appropriate colors and borders
#'   \item If \code{legend = FALSE}: Draws labels with pointing lines
#'   \item If \code{percentage = TRUE}: Adds percentage labels (filtered by
#'         \code{cutoff})
#'   \item If \code{rawNumber = TRUE} and \code{percentage = FALSE}: Adds raw
#'         number labels (filtered by \code{cutoff})
#'   \item If \code{legend = TRUE}: Draws a legend instead of on-chart labels
#'   \item Adds main title if provided
#' }
#' 
#' \strong{Label display logic:}
#' \itemize{
#'   \item If \code{percentage = TRUE}: Shows percentages (e.g., "25.0\%")
#'   \item Else if \code{rawNumber = TRUE}: Shows raw numbers from \code{x}
#'   \item Else: No numeric labels (only slice names if \code{legend = FALSE})
#'   \item Labels are only shown for slices with percentage >= \code{cutoff}
#' }
#' 
#' \strong{Label positioning:}
#' \itemize{
#'   \item Slice names (when \code{legend = FALSE}): Positioned at 1.1 * radius
#'         with pointing lines from 1.0 to 1.05 * radius
#'   \item Percentage/raw numbers: Positioned at \code{radius.innerlabel} *
#'         radius from center
#'   \item Label alignment: Slice names are left-aligned for right-side slices,
#'         right-aligned for left-side slices. Percentages/numbers are centered.
#' }
#' 
#' \strong{Color handling:}
#' \itemize{
#'   \item If \code{col = NULL}: Uses default pastel color palette (7 colors,
#'         recycled as needed)
#'   \item If \code{density} is specified: Uses \code{par("fg")} for shading
#'         instead of fill colors
#'   \item Colors are recycled if the vector is shorter than the number of slices
#' }
#' 
#' \strong{Differences from base \code{pie()}:}
#' \itemize{
#'   \item Adds percentage or raw number labels on slices
#'   \item Supports cutoff threshold for small slices
#'   \item Supports legend mode instead of on-chart labels
#'   \item Customizable label positioning (\code{radius.innerlabel})
#'   \item Uses \code{dev.hold()}/\code{dev.flush()} for smoother rendering
#' }
#' 
#' @author Jianhong Ou
#' @seealso \code{\link[graphics]{pie}} for the base pie chart function,
#'          \code{\link[graphics]{legend}} for legend customization,
#'          \code{\link[base]{formatC}} for number formatting
#' @keywords misc
#' @export
#' @importFrom grDevices as.graphicsAnnot dev.hold dev.flush 
#' @importFrom graphics plot.new plot.window polygon lines text title legend par
#' @examples
#' \dontrun{
#' ## Example 1: Basic pie chart with percentages
#' pie1(1:5, labels = LETTERS[1:5])
#' 
#' ## Example 2: Pie chart with raw numbers
#' pie1(c(10, 20, 30, 40), 
#'      labels = c("A", "B", "C", "D"),
#'      percentage = FALSE, rawNumber = TRUE)
#' 
#' ## Example 3: Pie chart with legend
#' pie1(c(25, 30, 20, 15, 10),
#'      labels = c("Type1", "Type2", "Type3", "Type4", "Type5"),
#'      legend = TRUE, legendpos = "bottomleft")
#' 
#' ## Example 4: Custom colors and higher cutoff
#' pie1(c(50, 30, 15, 4, 1),
#'      labels = c("Large", "Medium", "Small", "Tiny1", "Tiny2"),
#'      col = c("red", "blue", "green", "yellow", "purple"),
#'      cutoff = 0.05)  # Only show labels for slices >= 5%
#' 
#' ## Example 5: Clockwise with custom starting angle
#' pie1(1:4, labels = c("Q1", "Q2", "Q3", "Q4"),
#'      clockwise = TRUE, init.angle = 0)
#' 
#' ## Example 6: Labels closer to center
#' pie1(1:6, labels = LETTERS[1:6],
#'      radius.innerlabel = 0.5)  # Labels at 50% from center
#' 
#' ## Example 7: No percentage labels
#' pie1(c(10, 20, 30), labels = c("A", "B", "C"),
#'      percentage = FALSE, rawNumber = FALSE)
#' }
#' 
pie1 <- function (x, labels = names(x), edges = 200L, 
                  radius = 0.8, clockwise = FALSE, 
                  init.angle = if (clockwise) 90L else 0L, 
                  density = NULL, angle = 45L, 
                  col = NULL, border = NULL, lty = NULL, 
                  main = NULL, percentage = TRUE, rawNumber = FALSE, 
                  digits = 3L, cutoff = 0.01, 
                  legend = FALSE, legendpos = "topright", legendcol = 2L, 
                  radius.innerlabel = radius, ...) {
    if (!is.numeric(x) || any(is.na(x) | x < 0L)) {
        stop("'x' values must be non-negative", call. = FALSE)
    }
    if (is.null(labels)) {
        labels <- as.character(seq_along(x))
    } else {
        labels <- as.graphicsAnnot(labels)
    }
    rawX <- x
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    dev.hold()
    on.exit(dev.flush())
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "lightblue", "mistyrose", "lightcyan", 
              "lavender", "cornsilk", "pink")
    else par("fg")
    if (!is.null(col)) 
        col <- rep_len(col, nx)
    if (!is.null(border)) 
        border <- rep_len(border, nx)
    if (!is.null(lty)) 
        lty <- rep_len(lty, nx)
    angle <- rep(angle, nx)
    if (!is.null(density)) 
        density <- rep_len(density, nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
                border = border[i], col = col[i], lty = lty[i])
        if (!legend) {
            P <- t2xy(mean(x[i + 0:1]))
            lab <- as.character(labels[i])
            if (!is.na(lab) && nzchar(lab)) {
                lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
                text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
                     adj = ifelse(P$x < 0, 1, 0), ...)
            }
        }
    }
    if (percentage) {
        for (i in 1L:nx) {
            if (dx[i] > cutoff) {
                P <- t2xy(mean(x[i + 0:1]))
                text(radius.innerlabel * P$x, radius.innerlabel * P$y, 
                     paste(formatC(dx[i] * 100, digits = digits), "%", 
                           sep = ""), 
                     xpd = TRUE, 
                     adj = .5, ...)
            }
        }
    } else {
        if (rawNumber) {
            for (i in 1L:nx) {
                if (dx[i] > cutoff) {
                    P <- t2xy(mean(x[i + 0:1]))
                    text(radius.innerlabel * P$x, radius.innerlabel * P$y, 
                         rawX[i], xpd = TRUE, 
                         adj = .5, ...)
                }
            }
        }
    }
    if (legend) legend(legendpos, legend = labels, fill = col, 
                      border = "black", bty = "n", ncol = legendcol)
    title(main = main, ...)
    invisible(NULL)
}

#' Draw donut chart with ggplot2
#' 
#' @description 
#' Creates a donut chart (pie chart with a hole in the center) using ggplot2
#' with optional percentage or raw number labels for each slice. This function
#' provides a modern ggplot2-based alternative to \code{\link{pie1}} with similar
#' functionality for visualizing categorical data distributions.
#' 
#' @param x A namedvector of non-negative numerical quantities or a frequency table. 
#'        Must not contain NA or negative values. 
#' @param labels One or more expressions or character strings giving names for
#'        the slices. Default is \code{names(x)}. If \code{NULL}, labels are
#'        generated as \code{as.character(seq_along(x))}. 
#' @param percentage A logical value. If \code{TRUE} (default), percentage
#'        labels are displayed on each slice. Percentages are calculated as
#'        \code{100 * (slice_value / sum(x))}. Only slices with percentage >=
#'        \code{cutoff} are labeled. If \code{FALSE}, raw numbers (from \code{x})
#'        are displayed instead.
#' @param digits An integer specifying the number of significant digits to use
#'        for percentage formatting. Default is 3. See \code{\link[base]{formatC}}
#'        for details. Only used when \code{percentage = TRUE}.
#' @param cutoff A numeric value between 0 and 1 specifying the minimum
#'        percentage threshold for displaying labels. Slices with percentage <
#'        \code{cutoff} will not show percentage or raw number labels. Default is
#'        0.01 (1\%). This prevents cluttering the chart with labels for very small
#'        slices.
#' @param radii A numeric vector of length 2 specifying the inner and outer radii
#'        of the donut chart. The first value is the inner radius (hole size) and
#'        the second value is the outer radius (edge of donut). Default is
#'        \code{c(3, 4)}. In polar coordinates, these values represent the radial
#'        distance from the center. The first value must be less than the second.
#' 
#' @return Returns a ggplot2 object that can be further customized or printed.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates that \code{x} contains only non-negative numeric values
#'   \item Creates a data frame with values, labels, percentages, and positions
#'   \item Calculates cumulative positions for each slice
#'   \item Uses \code{ggplot2::geom_bar()} with \code{coord_polar()} to create
#'         the donut shape
#'   \item Sets \code{radii} to create the inner hole (donut effect)
#'   \item Adds percentage or raw number labels using \code{geom_text()}
#' }
#' 
#' \strong{Label display logic:}
#' \itemize{
#'   \item If \code{percentage = TRUE}: Shows percentages (e.g., "25.0\%")
#'   \item If \code{percentage = FALSE}: Shows raw numbers from \code{x}
#'   \item Labels are only shown for slices with percentage >= \code{cutoff}
#' }
#' 
#' \strong{Differences from \code{pie1()}:}
#' \itemize{
#'   \item Uses ggplot2 instead of base graphics
#'   \item Creates a donut chart (with hole) instead of a pie chart
#'   \item Returns a ggplot2 object that can be further customized
#'   \item More limited customization options (focuses on core functionality)
#' }
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{pie1}} for the base graphics pie chart function
#' @keywords misc
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_rect coord_polar geom_text theme_void
#' @importFrom ggplot2 scale_fill_manual xlim guides
#' @examples
#' \dontrun{
#' ## Example 1: Basic donut chart with percentages
#' donut(1:5, labels = LETTERS[1:5])
#' 
#' ## Example 2: Donut chart with raw numbers
#' donut(c(10, 20, 30, 40), 
#'        labels = c("A", "B", "C", "D"),
#'        percentage = FALSE)
#' 
#' ## Example 3: Higher cutoff
#' donut(c(50, 30, 15, 4, 1),
#'        labels = c("Large", "Medium", "Small", "Tiny1", "Tiny2"),
#'        cutoff = 0.05)  # Only show labels for slices >= 5%
#' 
#' ## Example 4: Thinner donut (larger hole)
#' donut(1:4, labels = c("Q1", "Q2", "Q3", "Q4"),
#'        radii = c(3.5, 4))
#' }
#' 
donut <- function(x, labels = names(x), 
                   percentage = TRUE, 
                   digits = 3L, 
                   cutoff = 0.01,
                   radii = c(3, 4)){
    # Validate input
    if (!is.numeric(x) || any(is.na(x) | any(x < 0L))) {
        stop("'x' values must be non-negative", call. = FALSE)
    }
    if (!is.vector(x)) {
        x_tmp <- as.vector(x)
        names(x_tmp) <- names(x)
        x <- x_tmp
    }
    if (radii[1] >= radii[2] || radii[1] < 0) {
        stop("'radii' must be a vector of two numbers with the first number less than the second number", call. = FALSE)
    }
    
    # Handle labels
    if (is.null(labels)) {
        labels <- as.character(seq_along(x))
    } else {
        labels <- as.character(labels)
    }
    
    # Calculate percentages
    total <- sum(x)
    percentages <- x / total
    
    # Create data frame
    df <- data.frame(
        label = labels,
        value = x,
        percentage = percentages,
        stringsAsFactors = FALSE
    )
    
    # Calculate positions for slices
    df$ymax <- cumsum(df$percentage)
    df$ymin <- c(0, head(df$ymax, n = -1))
    
    # Add xmin and xmax as columns for geom_rect (donut inner and outer radii)
    df$xmin <- radii[1]
    df$xmax <- radii[2]
    
    # Calculate label positions (middle of each slice)
    df$label_y <- (df$ymin + df$ymax) / 2
    
    # Create annotation text
    df$annotation <- ""
    if (percentage) {
        df$annotation[df$percentage >= cutoff] <- 
            paste0(formatC(df$percentage[df$percentage >= cutoff] * 100, 
                          digits = digits), "%")
    } else {
        df$annotation[df$percentage >= cutoff] <- 
            as.character(df$value[df$percentage >= cutoff])
    }
    
  # Determine if we need space for labels outside the donut
    has_labels <- any(nzchar(df$annotation))
    df$annotation <- ifelse(nzchar(df$annotation), df$annotation, "")
 
    # Create the plot
    # In coord_polar(theta = "y"):
    # - x (xmin/xmax) represents the RADIUS (distance from center)
    # - y (ymin/ymax) represents the ANGLE (position around circle, 0-1 maps to 0-2π)
    # For donut chart:
    # - xmin = radii[1] (inner radius, hole size)
    # - xmax = radii[2] (outer radius, edge of donut)
    # - ymin/ymax = angular positions of each slice (already calculated)
    
    # For donut chart to work properly:
    # - Rectangles are drawn from radii[1] (inner radius) to radii[2] (outer radius)
    # - xlim() must start from 0 or less than radii[1] to show the inner hole
    # - xlim() extends beyond radii[2] to show labels outside
    plot_xlim_min <- 0  # Start from center to show the hole
    plot_xlim_max <- if (has_labels) radii[2] + 1 else radii[2]
    p <- ggplot(df, aes_string(xmin = "xmin", xmax = "xmax", 
                               ymin = "ymin", ymax = "ymax", 
                               fill = "label")) +
        geom_rect() +
        coord_polar(theta = "y") +
        xlim(c(plot_xlim_min, plot_xlim_max)) +
        theme_void()+
        guides(fill = guide_legend(title = "Category"))
    
    # Add annotations if needed
    if (any(nzchar(df$annotation))) {
        df_label <- df[df$percentage >= cutoff & nzchar(df$annotation), ]
        
        if (nrow(df_label) > 0) {
            # Label position outside the donut
            df_label$x_pos <- radii[2] + 0.8
            df_label$y_pos <- df_label$label_y
            
            # Add text labels outside the donut
            p <- p + geom_text(data = df_label, 
                              aes_string(x = "x_pos", y = "y_pos",
                              label = "annotation"),
                              inherit.aes = FALSE,
                              size = 4, color = "black")
        }
    }
    
   p
}
