#' Create Venn diagrams for peak overlap analysis
#' 
#' @description 
#' Generates Venn diagrams to visualize overlaps between two to five sets of
#' peaks. The function also calculates statistical significance (p-values) to
#' determine whether observed overlaps are significant using either hypergeometric
#' test or permutation testing. The function can accept either a list of GRanges
#' objects or an \code{overlappingPeaks} object (from \code{\link{findOverlapsOfPeaks}}).
#' 
#' @param Peaks Either a list of two to five \link[GenomicRanges:GRanges-class]{GRanges}
#'        objects containing peak sets, or an \code{overlappingPeaks} object (from
#'        \code{\link{findOverlapsOfPeaks}}) that contains a \code{venn_cnt} element.
#'        If a list is provided, each element must be a GRanges object. Missing peak
#'        names will be automatically generated.
#' @param NameOfPeaks A character vector specifying the names of peak sets,
#'        e.g., \code{c("TF1", "TF2")}. These will be used as labels in the Venn
#'        diagram. If missing, names will be extracted from \code{names(Peaks)} or
#'        auto-generated. If \code{Peaks} is an \code{overlappingPeaks} object,
#'        names will be extracted from the \code{venn_cnt} column names if not provided.
#' @param maxgap An integer specifying the maximum gap (in base pairs) allowed
#'        between ranges for them to be considered overlapping. Default is \code{-1L},
#'        which means ranges must actually overlap (no gap allowed). See
#'        \code{\link[IRanges:findOverlaps-methods]{findOverlaps}} for details.
#' @param minoverlap An integer or numeric value specifying the minimum overlap
#'        required. If an integer >= 1, it specifies the minimum number of base pairs
#'        that must overlap. If \code{0 < minoverlap < 1}, it specifies the minimum
#'        percentage of the interval that must be covered. Default is \code{0L}
#'        (any overlap is considered).
#' @param totalTest A numeric value specifying the total number of possible peaks
#'        in the testing space. This is required for the hypergeometric test and
#'        should be much larger than the number of peaks in the largest peak set.
#'        If missing and \code{method = "hyperG"}, the function will estimate
#'        \code{totalTest} based on average peak width using the formula:
#'        \code{5e+7 / averagePeakWidth} (assuming human genome size and coding
#'        regions). See details below.
#' @param by A character string specifying the level at which overlaps are calculated.
#'        Options:
#'        \itemize{
#'          \item \code{"region"} (default): Based on genomic coordinates (chromosome,
#'                start, end). This is the most common approach for ChIP-seq peak analysis.
#'          \item \code{"feature"}: Based on feature annotations stored in the metadata
#'                columns of the GRanges object. Useful when peaks are annotated with
#'                gene IDs or other feature identifiers.
#'          \item \code{"base"}: At nucleotide-level resolution, calculating the actual
#'                number of overlapping base pairs.
#'        }
#' @param ignore.strand A logical value. When \code{TRUE} (default), strand
#'        information is ignored in overlap calculations. When \code{FALSE}, only
#'        ranges on the same strand are considered for overlap.
#' @param connectedPeaks A character string specifying how to count connected/overlapping
#'        peak groups. Options:
#'        \itemize{
#'          \item \code{"min"} (default): Counts the minimal number of involved peaks
#'                in each group of connected/overlapped peaks
#'          \item \code{"merge"}: Counts each group of connected peaks as only 1,
#'                regardless of how many peaks are involved
#'          \item \code{"keepAll"}: Counts all involved peaks for each peak list,
#'                providing detailed counts per list while also outputting counts
#'                as if \code{connectedPeaks} were set to "min"
#'          \item \code{"keepFirstListConsistent"}: Keeps the counts consistent with
#'                the first list. This option requires that \code{findOverlapsOfPeaks}
#'                was called with \code{connectedPeaks = "keepAll"}. If not, a warning
#'                is issued and the default is used.
#'        }
#' @param method A character string specifying the method for p-value calculation.
#'        Options:
#'        \itemize{
#'          \item \code{"hyperG"} (default): Uses hypergeometric test to calculate
#'                p-values for pairwise overlaps. Requires \code{totalTest} (or will
#'                estimate it if missing). Works for any number of peak sets.
#'          \item \code{"permutation"}: Uses permutation testing via
#'                \code{\link{peakPermTest}} to calculate p-values. Only works for
#'                pairwise comparisons (2 peak sets) and requires \code{TxDb}. Cannot
#'                be used with \code{overlappingPeaks} objects.
#'        }
#' @param TxDb An object of \link[GenomicFeatures:TxDb-class]{TxDb}. Required
#'        when \code{method = "permutation"}. Used to define the genomic space for
#'        permutation testing.
#' @param plot A logical value. If \code{TRUE} (default), a Venn diagram or Upset
#'        plot is created using \code{\link{plotggVennDiagram}}. If \code{FALSE},
#'        only the statistical results are returned without plotting.
#' @param plot.type Character string specifying the plot type: \code{"auto"}
#'        (default) to automatically select based on the number of sets,
#'        \code{"Venn"} for traditional Venn diagrams, or \code{"Upset"} for
#'        Upset plots. When \code{"auto"}, Venn diagrams are used for 5 or fewer
#'        sets, and Upset plots are used for more than 5 sets.
#' @param \dots Additional arguments to be passed to
#'        \code{\link{plotggVennDiagram}} for customizing the plot appearance.
#'        For Venn diagrams: \code{fill_color}, \code{stroke_color},
#'        \code{show_percentage}, etc. For Upset plots: \code{nsets},
#'        \code{nintersects}, \code{sets.bar.color}, etc. See
#'        \code{\link{plotggVennDiagram}} documentation for all available options.
#' 
#' @details
#' 
#' \strong{Input handling:}
#' \itemize{
#'   \item If \code{Peaks} is a list of GRanges objects, the function calculates
#'         Venn counts using \code{\link{getVennCounts}}
#'   \item If \code{Peaks} is an \code{overlappingPeaks} object (from
#'         \code{\link{findOverlapsOfPeaks}}), the function uses the existing
#'         \code{venn_cnt} data, which is more efficient
#'   \item Missing peak names are automatically generated using zero-padded
#'         numeric indices
#' }
#' 
#' \strong{Statistical testing:}
#' \itemize{
#'   \item \code{method = "hyperG"}: Calculates p-values for all pairwise overlaps
#'         using the hypergeometric test. If \code{totalTest} is missing, it is
#'         estimated as \code{round(5e+7 / averagePeakWidth)}, where
#'         \code{averagePeakWidth} is the median peak width across all input peaks.
#'         This estimation assumes human genome size (3.3e9 bp) with 3% coding/regulatory
#'         regions divided by average peak width.
#'   \item \code{method = "permutation"}: Calculates p-values using permutation
#'         testing, which is more computationally intensive but may be more appropriate
#'         for certain analyses. Only works for pairwise comparisons and requires
#'         \code{TxDb} to define the genomic space.
#' }
#' 
#' \strong{Visualization:}
#' \itemize{
#'   \item The function uses \code{\link{plotggVennDiagram}} to create visualizations
#'   \item When \code{connectedPeaks = "keepAll"}, additional count information
#'         is displayed in parentheses for overlapping regions
#'   \item If \code{totalTest} is provided, the count of items not in any set
#'         ("others") is displayed as a caption in Venn diagrams only. For Upset plots,
#'         this information is not displayed on the plot but is available in the
#'         return value (\code{vennCounts})
#'   \item Plot type is automatically selected: Venn diagrams for ≤5 sets,
#'         Upset plots for >5 sets
#'   \item Venn diagrams return standard ggplot objects that can be customized
#'         with ggplot2 syntax. Upset plots return \code{aplot} objects which
#'         cannot be modified with the ggplot2 \code{+} operator
#'   \item For customized graph options, see \code{\link{plotggVennDiagram}}
#' }
#' 
#' @return Returns a list with the following components:
#'        \itemize{
#'          \item \code{p.value}: A data frame or matrix containing p-values for
#'                pairwise overlaps. For \code{method = "hyperG"}, this is a data frame
#'                with binary columns for each peak set and a \code{pval} column.
#'                For \code{method = "permutation"}, this is similar but calculated
#'                via permutation testing. If permutation testing cannot be performed
#'                (e.g., with \code{overlappingPeaks} objects or missing \code{TxDb}),
#'                \code{p.value} is set to \code{NA}.
#'          \item \code{vennCounts}: An object of class \code{VennCounts} containing
#'                the overlap counts used for Venn diagram generation. This is a matrix
#'                with binary columns (0/1) for each peak set and a \code{"Counts"}
#'                column with the number of peaks in each overlap category.
#'        }
#' @export
#' @author Lihua Julie Zhu, Jianhong Ou
#' @seealso \link{findOverlapsOfPeaks}, 
#' \link{plotggVennDiagram}, \link{peakPermTest}
#' @examples 
#' if (interactive()){
#' peaks1 <- GRanges(seqnames=c("1", "2", "3"),
#'                   IRanges(start=c(967654, 2010897, 2496704),
#'                           end=c(967754, 2010997, 2496804), 
#'                           names=c("Site1", "Site2", "Site3")),
#'                   strand="+",
#'                   feature=c("a","b","f"))
#' peaks2 = GRanges(seqnames=c("1", "2", "3", "1", "2"), 
#'                  IRanges(start = c(967659, 2010898,2496700,
#'                                    3075866,3123260),
#'                          end = c(967869, 2011108, 2496920, 
#'                                  3076166, 3123470),
#'                          names = c("t1", "t2", "t3", "t4", "t5")), 
#'                  strand = c("+", "+", "-", "-", "+"), 
#'                  feature=c("a","b","c","d","a"))
#' makeVennDiagram(list(peaks1, peaks2), NameOfPeaks=c("TF1", "TF2"),
#'                 totalTest=100)
#' 
#' ###### 4-way diagram using annotated feature instead of chromosome ranges
#' 
#' makeVennDiagram(list(peaks1, peaks2, peaks1, peaks2), 
#'                 NameOfPeaks=c("TF1", "TF2","TF3", "TF4"), 
#'                 totalTest=100, by="feature")
#' }
#' @keywords graph
#' 

makeVennDiagram <- function(Peaks, NameOfPeaks, maxgap = -1L, minoverlap = 0L,
                            totalTest, by = c("region", "feature", "base"), 
                            ignore.strand = TRUE, 
                            connectedPeaks = c("min", "merge", "keepAll", 
                                               "keepFirstListConsistent"), 
                            method = c("hyperG", "permutation"), 
                            TxDb, plot = TRUE, 
                            plot.type = c("auto", "Venn", "Upset"),
                            ...) {
    ### Functions to be used
    getCountsList <- function(counts) {
        CountsList <- list()
        cnt <- 1
        for (i in seq_along(counts)) {
            CountsList[[i]] <- seq(cnt, length.out = counts[i])
            cnt <- cnt + counts[i]
        }
        CountsList
    }
    getVennList <- function(a, NameOfPeaks, CountsList) {
        .x <- lapply(NameOfPeaks, function(.ele, cl, a) {
            .y <- c()
            for (i in seq_len(nrow(a))) {
                if (a[i, .ele] != 0) {
                    .y <- c(.y, cl[[i]])
                }
            }
            .y
        }, CountsList, a)
        names(.x) <- NameOfPeaks
        .x
    }

    getPval <- function(venn_cnt, totalTest) {
        n <- which(colnames(venn_cnt) == "Counts") - 1 ## ncol(venn_cnt) - 1
        s <- apply(venn_cnt[, 1:n], 1, sum)
        venn_cnt_s <- venn_cnt[s == 2, , drop = FALSE]
        for (i in seq_len(nrow(venn_cnt_s))) {
            venn_cnt_s[i, n + 1] <- 
                sum(venn_cnt[as.numeric(venn_cnt[, 1:n] %*% venn_cnt_s[i, 1:n]) == 2,
                             n + 1, drop = TRUE])
        }
        cnt <- venn_cnt[, n + 1]
        cnt_s <- apply(venn_cnt[, 1:n], 2, function(.ele) {
            sum(cnt[as.logical(.ele)])
        })
        p.value <- apply(venn_cnt_s, 1, function(.ele) {
            ab <- cnt_s[as.logical(.ele[1:n])]
            a <- ab[1]
            b <- ab[2]
            a.and.b <- .ele[(n + 1)]
            phyper <- phyper(a.and.b - 1, b, totalTest - b, a, 
                            lower.tail = FALSE, log.p = FALSE)
        })
        cbind(venn_cnt_s[, 1:n, drop = FALSE], pval = p.value)
    }
    ### Check inputs
    method <- match.arg(method)
    plot.type <- match.arg(plot.type)
    if (missing(totalTest) && method == "hyperG") {
        message("Missing totalTest! totalTest is required for HyperG test. 
If totalTest is missing, pvalue will be calculated by estimating 
the total binding sites of encoding region of human.
totalTest = humanGenomeSize * (2%(codingDNA) + 
             1%(regulationRegion)) / ( 2 * averagePeakWidth )
          = 3.3e+9 * 0.03 / ( 2 * averagePeakWidth)
          = 5e+7 /averagePeakWidth")
    }
    if (missing(Peaks)) {
        stop("Missing 'Peaks' which is a list of peaks in GRanges", call. = FALSE)
    }
    connectedPeaks <- match.arg(connectedPeaks)
    stopifnot(is.logical(plot))
    by <- match.arg(by)
    n1 <- length(Peaks)
    olout_flag <- FALSE
    if (inherits(Peaks, "overlappingPeaks")) {
        if (!is.null(Peaks$venn_cnt) && inherits(Peaks$venn_cnt, "VennCounts")) {
            venn_cnt <- Peaks$venn_cnt
            n1 <- which(colnames(venn_cnt) == "Counts") - 1 ## ncol(venn_cnt) - 1
            if (missing(NameOfPeaks)) {
                NameOfPeaks <- colnames(venn_cnt)[1:n1]
            }
            olout_flag <- TRUE
        } else {
            stop("Input is object of overlappingPeaks, ",
                 "but it does not have venn counts data", call. = FALSE)
        }
    }
    if (missing(NameOfPeaks) || mode(NameOfPeaks) != "character") {
        warning("Missing required character vector NameOfPeaks. 
            NameOfPeaks will be extract from the names of input Peaks.")
        NameOfPeaks <- names(Peaks)
        if (is.null(NameOfPeaks)) {
            NameOfPeaks <- paste("peaks", seq_along(Peaks), sep = "")
        }
    }
    n2 <- length(NameOfPeaks)
    if (n1 < 2L) {
        stop("At least 2 peak sets are required", call. = FALSE)
    }
    if (n1 > 5L) {
        stop("Maximum 5 peak sets are supported", call. = FALSE)
    }
    if (n1 > n2) {
        stop("The number of elements in 'NameOfPeaks' is less ",
             "than the number of elements in 'Peaks'. They must be equal", 
             call. = FALSE)
    }
    if (n1 < n2) {
        warning("The number of elements in 'NameOfPeaks' is larger than ",
                "the number of elements in 'Peaks'. ",
                "'NameOfPeaks' will be truncated to match 'Peaks'", 
                call. = FALSE)
        NameOfPeaks <- NameOfPeaks[seq_len(n1)]
    }
    NameOfPeaks <- make.names(NameOfPeaks, unique = TRUE, allow_ = TRUE)
    if (!olout_flag) {
        for (i in seq_len(n1)) {
            if (!inherits(Peaks[[i]], "GRanges")) {
                stop("Element ", i, " in 'Peaks' is not a valid GRanges object", 
                     call. = FALSE)
            }
            if (is.null(names(Peaks[[i]]))) {
                n_peaks <- length(Peaks[[i]])
                names(Peaks[[i]]) <- formatC(seq_len(n_peaks), 
                                             width = nchar(n_peaks), 
                                             flag = "0")
            }
        }
        if (!missing(totalTest)) {
            max_peaks <- max(vapply(Peaks, length, FUN.VALUE = 0L))
            if (totalTest < max_peaks) {
                stop("'totalTest' specifies the total number of possible peaks ",
                     "in the testing space. It should be larger than the largest ",
                     "peak count in the input sets. Please see more details ",
                     "at http://pgfe.umassmed.edu/ChIPpeakAnno/FAQ.html", 
                     call. = FALSE)
            }
        }
        names(Peaks) <- NameOfPeaks
        venn_cnt <- getVennCounts(Peaks, maxgap = maxgap, 
                                  minoverlap = minoverlap, 
                                  by = by, ignore.strand = ignore.strand, 
                                  connectedPeaks = 
                                      ifelse(connectedPeaks == "keepFirstListConsistent", 
                                             "keepAll", connectedPeaks))
    }
    colnames(venn_cnt)[1:n1] <- NameOfPeaks
    venn_cnt1 <- venn_cnt
    if (connectedPeaks == "keepFirstListConsistent") {
        ## keepAll to getVennCounts
        if (!grepl("^count\\.", colnames(venn_cnt)[ncol(venn_cnt)])) {
            warning("If connectedPeaks set to keepFirstListConsistent,",
                    "connectedPeaks of findOverlapsOfPeaks must be keepAll.",
                    "Setting connectedPeaks to default.")
        } else {
            venn_cnt1[venn_cnt[, 1] == 1, "Counts"] <- venn_cnt[venn_cnt[, 1] == 1, n1 + 2]
        }
        connectedPeaks <- "min"
    }
    Counts <- getCountsList(venn_cnt1[, "Counts"])
    vennx <- getVennList(venn_cnt1, NameOfPeaks, Counts)
    if (method == "hyperG") {
        if (!missing(totalTest)) {
            otherCount <- totalTest - sum(venn_cnt[, "Counts"])
            venn_cnt[1, "Counts"] <- otherCount
            p.value <- getPval(venn_cnt, totalTest)
        } else {
            otherCount <- venn_cnt[1, "Counts"]
            if (inherits(Peaks, "overlappingPeaks")) {
                averagePeakWidth <- median(unlist(lapply(Peaks$peaklist, width)))
            } else {
                averagePeakWidth <- median(unlist(lapply(Peaks, width))) 
            }
            p.value <- getPval(venn_cnt, round(5e+7 / averagePeakWidth))
        }
    } else {
        ## method == "permutation"
        p.value <- NA
        if (olout_flag) {
            warning("Input is an object of overlappingPeaks.
                  Can not do permutation test. Please try ?peakPermTest")
        } else {
            if (missing(TxDb)) {
                warning("TxDb is missing. Please try ?peakPermTest later.")
            } else {
                n <- which(colnames(venn_cnt) == "Counts") - 1 ## ncol(venn_cnt) - 1
                s <- apply(venn_cnt[, 1:n], 1, sum)
                venn_cnt_s <- venn_cnt[s == 2, , drop = FALSE]
                p.value <- apply(venn_cnt_s, 1, function(.ele) {
                    AB <- Peaks[as.logical(.ele)]
                    pt <- peakPermTest(AB[[1]], AB[[2]], TxDb = TxDb, ...)
                    pt$pval
                })
                p.value <- cbind(venn_cnt_s[, 1:n, drop = FALSE], pval = p.value)
            }
        }
        if (!missing(totalTest)) {
            otherCount <- totalTest - sum(venn_cnt[, "Counts"])
            venn_cnt[1, "Counts"] <- otherCount
        } else {
            otherCount <- venn_cnt[1, "Counts"]
        }
    }
    if (otherCount == 0) {
        otherCount <- NULL
    }
    if (plot) {
        p <- plotggVennDiagram(venn_cnt1, vennx, plot.type = plot.type, 
                               otherCounts = otherCount, ...)
        # Ensure plot is displayed - print() will handle device creation if needed
        # For aplot objects (Upset plots), print() should render to the active device
        print(p)
    }
    return(list(p.value = p.value, vennCounts = venn_cnt))
}

#' Create Venn diagrams or Upset plots using ggvenn and ggVennDiagram
#' 
#' @description 
#' Generates Venn diagrams or Upset plots to visualize overlaps between peak sets
#' using the \code{ggvenn} package for Venn diagrams and \code{ggVennDiagram} for
#' Upset plots. This function provides an improved interface for creating
#' visualizations. Venn diagrams return standard ggplot objects that can be
#' customized with ggplot2 syntax, while Upset plots return \code{aplot} objects
#' (composite plots) that cannot be modified with the ggplot2 \code{+} operator.
#' 
#' @param venn_cnt A data frame or matrix containing Venn counts, typically obtained from
#'   \code{\link{getVennCounts}} or as part of an \code{overlappingPeaks} object.
#'   Must have a "Counts" column and binary columns for each set.
#' @param vennx A named list where each element is a vector of indices
#'   representing items in that set. Typically obtained from \code{getVennList()}
#'   or as part of the Venn diagram calculation process.
#' @param plot.type Character string specifying the plot type: \code{"Venn"} for
#'   traditional Venn diagrams, \code{"Upset"} for Upset plots, or \code{"auto"}
#'   (default) to automatically select based on the number of sets. When
#'   \code{"auto"}, Venn diagrams are used for 5 or fewer sets, and Upset plots
#'   are used for more than 5 sets.
#' @param otherCounts Optional numeric value. Count of items not in any of the
#'   sets (for display purposes).
#' @param \dots Additional arguments to be passed to either
#'   \code{\link[ggvenn:ggvenn]{ggvenn}} (for Venn diagrams) or
#'   \code{\link[ggVennDiagram:plot_upset]{plot_upset}} (for Upset plots).
#'   For Venn diagrams: \code{fill_color}, \code{stroke_color}, \code{stroke_linetype},
#'   \code{show_elements}, \code{show_percentage}, \code{set_name_size}, etc.
#'   For Upset plots: \code{nintersects}, \code{order.intersect.by},
#'   \code{order.set.by}, \code{relative_height}, \code{relative_width}, etc.
#'   See \code{\link[ggVennDiagram:plot_upset]{plot_upset}} for all available options.
#' 
#' @return 
#'   \itemize{
#'     \item For Venn diagrams: A \code{ggplot} object that can be further customized
#'           using ggplot2 syntax (e.g., \code{+ theme()}, \code{+ labs()}, etc.)
#'     \item For Upset plots: An \code{aplot} object (inherits from \code{ggplot})
#'           that cannot be modified with the ggplot2 \code{+} operator. Attempting
#'           to add layers or modify the plot using \code{+} will return \code{NULL}
#'           and break the plot object.
#'   }
#' 
#' @details 
#' This function uses the
#' \code{ggvenn} package for Venn diagrams and \code{ggVennDiagram} for Upset plots
#' instead of \code{VennDiagram}. The main advantages are:
#' \itemize{
#'   \item Better integration with ggplot2 (Venn diagrams return standard ggplot objects)
#'   \item More flexible customization options with \code{ggvenn} for Venn diagrams
#'   \item Better visualization quality and aesthetics
#'   \item Support for both Venn diagrams and Upset plots in one function
#' }
#' 
#' \strong{Important limitations for Upset plots:}
#' \itemize{
#'   \item Upset plots return \code{aplot} objects which cannot be customized with
#'         the ggplot2 \code{+} operator (e.g., \code{+ labs()}, \code{+ theme()})
#'   \item Attempting to modify Upset plots with \code{+} will return \code{NULL}
#'         and break the plot object
#'   \item The \code{otherCounts} parameter is only displayed for Venn diagrams.
#'         For Upset plots, this information is not shown on the plot but can be
#'         accessed from the \code{vennCounts} in the return value of
#'         \code{\link{makeVennDiagram}}
#' }
#' 
#' \strong{Automatic selection} (default): When \code{plot.type = "auto"}, the
#' function automatically chooses the appropriate visualization:
#' \itemize{
#'   \item \code{<= 5 sets}: Uses Venn diagrams
#'   \item \code{> 5 sets}: Uses Upset plots
#' }
#' 
#' \strong{Venn diagrams} are best for visualizing overlaps between 2-5 sets,
#' providing an intuitive circular representation of set relationships.
#' 
#' \strong{Upset plots} are better for visualizing overlaps between many sets
#' (typically more than 5), where traditional Venn diagrams become difficult to
#' interpret. Upset plots provide:
#' \itemize{
#'   \item Set sizes in a bar chart on the left
#'   \item Intersection sizes in a bar chart on the top
#'   \item A binary matrix showing which sets participate in each intersection
#'   \item Customizable display of intersections via function arguments
#' }
#' 
#' \strong{Note:} Upset plots return \code{aplot} objects that cannot be modified
#' after creation. All customization must be done through function arguments passed
#' via \code{...}, not through ggplot2 syntax.
#' 
#' The function expects the same data format as \code{makeVennDiagram()}, making it
#' easy to switch between visualization methods.
#' 
#' @note 
#' \itemize{
#'   \item The \code{ggvenn} package must be installed for Venn diagrams. Install it with:
#'         \code{install.packages("ggvenn")}
#'   \item For Upset plots, the \code{ggVennDiagram} package is required:
#'         \code{install.packages("ggVennDiagram")}
#'   \item Venn diagrams return standard \code{ggplot} objects that can be customized
#'         with ggplot2 syntax
#'   \item Upset plots return \code{aplot} objects that \strong{cannot} be modified
#'         with the ggplot2 \code{+} operator. Attempting to do so will break the plot
#'   \item The \code{otherCounts} caption is only displayed for Venn diagrams, not
#'         Upset plots
#' }
#' 
#' @author Haibo Liu
#' @seealso \code{\link{makeVennDiagram}}, \code{\link{getVennCounts}},
#'   \code{\link[ggvenn:ggvenn]{ggvenn}},
#'   \code{\link[ggVennDiagram:plot_upset]{plot_upset}}
#' @keywords graph
#' @export
#' @importFrom ggvenn ggvenn
#' @importFrom ggVennDiagram Venn plot_upset
#' @importFrom ggplot2 annotate labs
#' @examples
#' if (interactive() && requireNamespace("ggvenn", quietly = TRUE)) {
#'   peaks1 <- GRanges(seqnames = c("1", "2", "3"),
#'                     IRanges(start = c(967654, 2010897, 2496704),
#'                             end = c(967754, 2010997, 2496804),
#'                             names = c("Site1", "Site2", "Site3")),
#'                     strand = "+")
#'   peaks2 <- GRanges(seqnames = c("1", "2", "3"),
#'                     IRanges(start = c(967659, 2010898, 2496700),
#'                             end = c(967869, 2011108, 2496920),
#'                             names = c("t1", "t2", "t3")),
#'                     strand = "+")
#'   peaks3 <- GRanges(seqnames = c("1", "2", "3"),
#'                     IRanges(start = c(967660, 2010900, 2496705),
#'                             end = c(967870, 2011110, 2496925),
#'                             names = c("p1", "p2", "p3")),
#'                     strand = "+")
#'   
#'   # Get Venn counts
#'   venn_cnt <- getVennCounts(list(peaks1, peaks2, peaks3))
#'   
#'   # Create vennx list: indices representing items in each set
#'   CountsList <- lapply(venn_cnt[, "Counts"], function(x) seq_len(x))
#'   n_sets <- ncol(venn_cnt) - 1L
#'   set_names <- colnames(venn_cnt)[seq_len(n_sets)]
#'   vennx <- lapply(set_names, function(name) {
#'     indices <- c()
#'     for (i in seq_len(nrow(venn_cnt))) {
#'       if (venn_cnt[i, name] != 0) {
#'         indices <- c(indices, CountsList[[i]])
#'       }
#'     }
#'     indices
#'   })
#'   names(vennx) <- set_names
#'   
#'   # Plot using Venn diagram (default)
#'   plotggVennDiagram(venn_cnt, vennx, plot.type = "Venn")
#'   
#'   # Plot using Upset plot (better for many sets)
#'   # Note: Requires ggVennDiagram package
#'   if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
#'     plotggVennDiagram(venn_cnt, vennx, plot.type = "Upset",
#'                       nintersects = 10)
#'   }
#' }
plotggVennDiagram <- function(venn_cnt, vennx, 
                              plot.type = c("auto", "Venn", "Upset"),
                              otherCounts = NULL, ...) {
  # Validate inputs (shared validation for both plot types)
  if (missing(venn_cnt)) {
    stop("'venn_cnt' is required", call. = FALSE)
  }
  # venn_cnt can be a data frame or matrix (VennCounts object)
  if (!is.data.frame(venn_cnt) && !is.matrix(venn_cnt)) {
    stop("'venn_cnt' must be a data frame or matrix with Venn counts", 
         call. = FALSE)
  }
  if (missing(vennx) || !is.list(vennx)) {
    stop("'vennx' must be a named list of vectors", call. = FALSE)
  }
  
  # Extract set names from venn_cnt
  # Work with matrix or data frame directly - both support colnames()
  n_sets <- which(colnames(venn_cnt) == "Counts") - 1L
  if (n_sets < 1L) {
    stop("'venn_cnt' must contain a 'Counts' column and set columns", 
         call. = FALSE)
  }
  
  set_names <- colnames(venn_cnt)[seq_len(n_sets)]
  if (is.null(names(vennx))) {
    if (length(vennx) == length(set_names)) {
      names(vennx) <- set_names
    } else {
      stop("Length of 'vennx' does not match number of sets in 'venn_cnt'",
         call. = FALSE)
    }
  }
  
  # Match plot.type argument and handle automatic selection
  plot.type <- match.arg(plot.type)
  if (plot.type == "auto") {
    # Automatically select: Venn for <= 5 sets, Upset for > 5 sets
    plot.type <- if (n_sets > 5L) "Upset" else "Venn"
  }
  
  # Prepare arguments
  dots <- list(...)
  
  # Branch based on plot type
  if (plot.type == "Venn") {
    # Check for ggvenn package
    if (!requireNamespace("ggvenn", quietly = TRUE)) {
      stop("The 'ggvenn' package is required for Venn diagrams but not installed. ",
           "Please install it with: install.packages('ggvenn')",
           call. = FALSE)
    }
    
    # ggvenn expects a named list where each element contains items in that set
    # vennx already has this format (indices representing items)
    # ggvenn takes the data as the first argument (named 'data')
    args <- c(list(data = vennx, show_percentage = FALSE), dots)
    p <- do.call(ggvenn::ggvenn, args)
    
    # Add annotation for otherCounts if provided
    # Use labs(caption = ...) which positions text at the bottom with proper spacing
    if (!is.null(otherCounts) && otherCounts > 0L) {
      p <- p + 
        ggplot2::labs(caption = paste("others:", otherCounts)) + 
        ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0.5), # Centers the caption
                       plot.caption.position = "plot" # Aligns caption relative to the whole plot area
        )
    }
   } else {
    # plot.type == "Upset"
    # Check for ggVennDiagram package
    if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
      stop("The 'ggVennDiagram' package is required for Upset plots but not installed. ",
           "Please install it with: install.packages('ggVennDiagram')",
           call. = FALSE)
    }
    
    # ggVennDiagram::plot_upset() expects a Venn object
    # Create a Venn object from vennx
    # ggVennDiagram::Venn expects a named list where each element contains items in that set
    # Ensure all sets have the same type (numeric) to avoid "All sets should have same classes" error
    vennx_typed <- lapply(vennx, function(x) as.numeric(x))
    venn_obj <- ggVennDiagram::Venn(vennx_typed)

    # Call plot_upset with the Venn object
    args <- c(list(venn = venn_obj), dots)
    p <- do.call(ggVennDiagram::plot_upset, args)
    
    # Note: ggVennDiagram::plot_upset() returns an aplot object (composite plot),
    # which cannot be customized with ggplot2 + operator.
  }
  p
}
