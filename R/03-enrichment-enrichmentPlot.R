#' Plot enrichment results
#' 
#' @description 
#' Creates visualization plots for GO, KEGG, or Reactome enrichment results.
#' The function automatically detects whether the input contains single-sample or
#' multi-sample data and generates appropriate visualizations:
#' \itemize{
#'   \item \strong{Single-sample}: Bar plots showing enriched terms with
#'         -log10(p-value) on the y-axis and gene counts as fill color
#'   \item \strong{Multi-sample}: Dot plots showing enriched terms across
#'         multiple samples, with terms on the y-axis, samples on the x-axis,
#'         colored by -log10(p-value) and sized by gene ratio
#' }
#' 
#' The function supports multiple enrichment categories (e.g., Biological Process,
#' Molecular Function, Cellular Component for GO) and displays them in separate
#' facets. Terms can be ordered by p-value, term ID, or left unordered.
#' 
#' @param res Output of \code{\link{getEnrichedGO}} or \code{\link{getEnrichedPATH}}.
#'        Can be:
#'        \itemize{
#'          \item A data frame: Single enrichment result (will be wrapped in a list)
#'          \item A list with one element: Single-sample enrichment result
#'          \item A list of lists: Multi-sample enrichment results (one list per sample)
#'        }
#'        The data frame(s) must contain columns: term ID (ending in ".id"),
#'        term description (ending in ".term"), "pvalue", "count.InDataset",
#'        "count.InGenome", and optionally "source" (for multi-sample data).
#' @param n An integer specifying the number of top terms to plot. Default is
#'        \code{20}. Terms are selected based on the ordering specified by
#'        \code{orderBy} parameter.
#' @param strlength An integer or \code{Inf} specifying the maximum length of
#'        term descriptions. Descriptions longer than this will be truncated
#'        with "..." appended. Default is \code{Inf} (no truncation).
#' @param style A character string specifying the plot orientation:
#'        \itemize{
#'          \item \code{"v"} (default): Vertical orientation (bar plots) or
#'                horizontal axis (dot plots)
#'          \item \code{"h"}: Horizontal orientation (bar plots flipped) or
#'                vertical axis (dot plots flipped)
#'        }
#' @param label_wrap An integer specifying the number of characters for soft
#'        wrapping of term labels. Default is \code{40}. Labels longer than
#'        this will be wrapped across multiple lines.
#' @param label_substring_to_remove A character string (or \code{NULL}) specifying
#'        a common substring to remove from all labels. Useful for removing
#'        redundant prefixes like organism names. Default is \code{NULL}.
#'        Special characters must be escaped. For example, to remove
#'        "Homo sapiens (human)" from labels, use: \code{"Homo sapiens \\(human\\)"}.
#' @param orderBy A character string specifying how to order the terms:
#'        \itemize{
#'          \item \code{"pvalue"} (default): Order by p-value (most significant first)
#'          \item \code{"termId"}: Order by term ID (alphabetical/numerical)
#'          \item \code{"none"}: No ordering (preserves original order)
#'        }
#' @author Jianhong Ou, Kai Hu
#' @return Returns an object of class \code{\link[ggplot2]{ggplot}} that can be
#'        further customized using ggplot2 syntax. The plot shows:
#'        \itemize{
#'          \item \strong{Single-sample bar plots}: Terms on x-axis, -log10(p-value)
#'                on y-axis, gene counts as fill color, with facets for each
#'                enrichment category
#'          \item \strong{Multi-sample dot plots}: Terms on y-axis, samples on
#'                x-axis, colored by -log10(p-value), sized by gene ratio
#'                (count.InDataset / count.InGenome), with facets for each
#'                enrichment category
#'        }
#'        If there are fewer than 2 terms to plot, a warning is issued and
#'        the plot data frame is returned instead of a plot.
#' @importFrom ggplot2 ggplot aes geom_bar geom_point scale_x_discrete 
#' scale_y_continuous geom_text xlab ylab theme_classic theme 
#' facet_grid expansion element_text coord_flip
#' @importFrom stats reorder
#' @importFrom scales label_wrap
#' @export
#' @seealso \code{\link{getEnrichedGO}}, \code{\link{getEnrichedPATH}}
#' @examples 
#' ## Example 1: Plot single-sample enrichment results
#' data(enrichedGO)
#' enrichmentPlot(enrichedGO)
#' 
#' ## Example 2: Customize plot appearance
#' enrichmentPlot(enrichedGO, n = 10, style = "h", 
#'                orderBy = "pvalue", strlength = 50)
#' 
#' ## Example 3: Remove common substring from labels
#' enrichmentPlot(enrichedGO, 
#'                label_substring_to_remove = "Homo sapiens \\(human\\)")
#' 
#' if (interactive() || Sys.getenv("USER") == "jou") {
#'     ## Example 4: Multi-sample comparison
#'     library(org.Hs.eg.db)
#'     library(GO.db)
#'     bed <- system.file("extdata", "MACS_output.bed", package = "ChIPpeakAnno")
#'     gr1 <- toGRanges(bed, format = "BED", header = FALSE)
#'     gff <- system.file("extdata", "GFF_peaks.gff", package = "ChIPpeakAnno")
#'     gr2 <- toGRanges(gff, format = "GFF", header = FALSE, skip = 3)
#'     library(EnsDb.Hsapiens.v75) ## (hg19)
#'     annoData <- toGRanges(EnsDb.Hsapiens.v75)
#'     gr1.anno <- annoPeaks(gr1, annoData)
#'     gr2.anno <- annoPeaks(gr2, annoData)
#'     over <- lapply(GRangesList(gr1 = gr1.anno, gr2 = gr2.anno), 
#'                    getEnrichedGO, orgAnn = "org.Hs.eg.db",
#'                    maxP = .05, minGOterm = 10, condense = TRUE)
#'     ## Single-sample plot
#'     enrichmentPlot(over$gr1)
#'     ## Horizontal orientation
#'     enrichmentPlot(over$gr2, style = "h")
#'     ## Multi-sample dot plot (when input is list of lists)
#'     enrichmentPlot(over)
#' }
enrichmentPlot <- function(res, n = 20, strlength = Inf,
                           style = c("v", "h"),
                           label_wrap = 40,
                           label_substring_to_remove = NULL,
                           orderBy = c("pvalue", "termId", "none")) {
    if (is.data.frame(res)) {
        res <- list(path = res)
    }
    stopifnot("n must be a numeric(1)" = is.numeric(n) && length(n) == 1)
    orderBy <- match.arg(orderBy)
    style <- match.arg(style)
    stopifnot(is.integer(as.integer(label_wrap)))
    if (is.list(res[[1]]) && is.data.frame(res[[1]][[1]])) {
        ## List of list, output of getEnrichedGO for multiple samples
        ## Dot plot
        res <- swapList(res)
        res <- lapply(res, function(.ele) {
            ## Unlist and add source
            .e <- do.call(rbind, .ele)
            .e$source <- rep(names(.ele), vapply(.ele, nrow, FUN.VALUE = 0))
            .e
        })
    }
    p <- lapply(res, function(.ele) {
        cn <- colnames(.ele)
        cn.id <- cn[grepl("\\.id$", cn)]
        cn.term <- cn[grepl("\\.term$", cn)]
        if (length(.ele$source) != nrow(.ele)) {
            .ele$source <- "undefined"
        }
        if (nrow(.ele) == 0) {
            return(data.frame())
        }
        if (!all(c(cn.id, cn.term, "pvalue", "count.InDataset", 
                  "count.InGenome", "source") %in% colnames(.ele))) {
            return(data.frame())
        }
        plotdata <- .ele[!is.na(.ele$pvalue), 
                         c(cn.id, cn.term, "pvalue", "count.InDataset", 
                           "count.InGenome", "source")]
        plotdata <- as.data.frame(plotdata)
        plotdata <- unique(plotdata)
        plotdata$qvalue <- -1 * log10(plotdata$pvalue)
        plotdata <- switch(orderBy,
                          "pvalue" = plotdata[order(plotdata$pvalue), ],
                          "termId" = plotdata[order(plotdata[, cn.term]), ],
                          plotdata)
        if (nrow(plotdata) > n) {
            plotdata <- plotdata[seq.int(n), ]
        }
        plotdata$Description <- shortStrs(plotdata[, cn.term], len = strlength)
        plotdata$Count <- plotdata$count.InDataset
        plotdata$GeneRatio <- plotdata$count.InDataset / plotdata$count.InGenome
        plotdata
    })
    plotdata <- do.call(rbind, p)
    plotdata$category <- rep(names(res), vapply(p, nrow, FUN.VALUE = 0))
    plotdata$Description <- removeLabelSubstring(plotdata, label_substring_to_remove)

    if (nrow(plotdata) < 2) {
        warning("two less data to plot")
        return(plotdata)
    }

    decreasing <- style == "h"
    if (all(plotdata$source == "undefined")) {
        p <- ggplot(plotdata, 
                   aes(x = switch(orderBy,
                                  "pvalue" = reorder(Description, pvalue, 
                                                     decreasing = decreasing),
                                  "termId" = reorder(Description, plotdata[, 2]),
                                  "Description"),
                       y = qvalue, fill = Count, label = Count))
        p <- p +
            geom_bar(stat = "identity") +
            scale_x_discrete(label = label_wrap(as.integer(label_wrap))) +
            scale_y_continuous(expand = expansion(mult = c(0, .1))) +
            xlab("") + ylab("-log10(p-value)") +
            theme_classic() + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
            facet_grid(~ category, scales = "free_x", space = "free_x")
        
        if (style == "v") {
            p <- p + geom_text(vjust = -.1)
        } else if (style == "h") {
            p <- p + geom_text(hjust = -.1) + coord_flip()
        }
    } else {
        ## Multiple samples dot plot
        p <- ggplot(plotdata,
                   aes(y = switch(orderBy, 
                                  "pvalue" = reorder(Description, pvalue, 
                                                     decreasing = decreasing),
                                  "termId" = reorder(Description, plotdata[, 2]),
                                  Description),
                       x = source, color = qvalue, size = GeneRatio))
        p <- p +
            geom_point() + theme_classic() +
            facet_grid(~ category, scales = "free_x", space = "free_x")
        
        if (style == "h") {
            p <- p + coord_flip()
        }
    }
    p
}

## Internal helper function to remove common substrings from labels
## @param plotdata Data frame with Description column
## @param label_substring_to_remove Character string to remove (or NULL)
## @return Character vector of modified descriptions
removeLabelSubstring <- function(plotdata, label_substring_to_remove = NULL) {
    if (!is.null(label_substring_to_remove)) {
        gsub(label_substring_to_remove, "", plotdata$Description)
    } else {
        plotdata$Description
    }
}

## Internal helper function to truncate strings intelligently
## Truncates strings to specified length, preserving the last word
## @param strs Character vector of strings to truncate
## @param len Maximum length (default 60)
## @return Character vector of truncated strings (made unique)
shortStrs <- function(strs, len = 60) {
    if (length(strs) == 0) {
        return(strs)
    }
    strs <- as.character(strs)
    shortStr <- function(str, len = 60) {
        stopifnot(length(str) == 1)
        stopifnot(is.character(str))
        if (nchar(str) <= len) {
            return(str)
        }
        strs <- strsplit(str, " ")[[1]]
        nc <- nchar(strs)
        nclast <- nc[length(nc)] + 3
        paste0(substring(str, first = 1, last = len - nclast), "...",
               strs[length(strs)])
    }
    strs <- sapply(strs, shortStr, len = len)
    make.unique(strs)
}
