#' Plot genomic element distribution
#' 
#' @description 
#' Categorizes and visualizes the distribution of peaks across genomic elements
#' (promoters, gene bodies, exons, introns, UTRs, etc.) using pie charts or
#' bar plots. Supports multiple categorization levels (gene level, exon/intron,
#' and exon sub-types) with customizable precedence rules for overlapping
#' annotations. The function can work with a single peak set (GRanges) or
#' multiple peak sets (GRangesList) for comparative analysis.
#' 
#' @details 
#' The distribution is calculated at multiple levels:
#' \itemize{
#'   \item \code{geneLevel}: Categorizes peaks as promoter region, gene body,
#'         gene downstream, and distal intergenic region. Promoter and downstream
#'         regions are defined relative to transcript boundaries.
#'   \item \code{ExonIntron}: Categorizes peaks as exon, intron, or intergenic.
#'         Exons and introns are extracted from the TxDb annotation.
#'   \item \code{Exons}: Categorizes peaks as 5' UTR, 3' UTR, CDS, or other exon.
#'         This provides fine-grained classification of exonic regions.
#'   \item \code{promoterLevel} (optional): Subdivides the promoter region into
#'         multiple bins based on distance from TSS. Only used if
#'         \code{promoterLevel} parameter is provided.
#' }
#' 
#' \strong{Precedence rules:} When a peak overlaps multiple genomic elements,
#' the precedence follows the order of labels definition. For example, in
#' \code{ExonIntron}, if a peak overlaps both exon and intron, and "exon" is
#' specified before "intron" in the labels, then only "exon" will be counted.
#' The precedence is applied in the order: promoter > gene body > downstream >
#' intergenic (for geneLevel), and similarly for other levels.
#' 
#' \strong{Counting methods:} When \code{nucleotideLevel = FALSE} (default),
#' each peak is counted once and assigned to the highest-priority overlapping
#' element. When \code{nucleotideLevel = TRUE}, the function counts each
#' nucleotide within peaks, providing a nucleotide-centric view that accounts
#' for peak width.
#' 
#' \strong{Visualization:} For a single GRanges object, a donut/pie chart is
#' generated showing the distribution across all categories. For a GRangesList,
#' bar plots are generated with separate facets for each category level,
#' allowing comparison across multiple peak sets.
#' 
#' @param peaks A \link[GenomicRanges:GRanges-class]{GRanges} object or
#'        \link[GenomicRanges:GRangesList-class]{GRangesList} containing peak
#'        regions to categorize. If a GRangesList is provided, each element
#'        represents a separate peak set for comparative analysis.
#' @param TxDb An object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#'        containing transcript annotation data. Used to extract gene, transcript,
#'        exon, intron, UTR, and CDS regions.
#' @param seqlev A character vector specifying which sequence levels (chromosomes)
#'        should be included in the analysis. Default is all sequence levels
#'        present in both peaks and TxDb (intersection).
#' @param nucleotideLevel A logical value. When \code{FALSE} (default), each peak
#'        is counted once (peak-centric view). When \code{TRUE}, each nucleotide
#'        within peaks is counted (nucleotide-centric view), accounting for peak
#'        width. This is useful when peaks vary significantly in size.
#' @param ignore.strand A logical value. When \code{TRUE} (default), strand
#'        information is ignored when determining overlaps. When \code{FALSE},
#'        only ranges on the same strand are considered overlapping.
#' @param promoterRegion A named numeric vector with elements \code{upstream}
#'        and \code{downstream} specifying the promoter region boundaries relative
#'        to transcript start sites (TSS). Default is \code{c(upstream = 2000,
#'        downstream = 100)}, meaning 2kb upstream and 100bp downstream of TSS.
#' @param promoterLevel An optional named list with elements \code{breaks},
#'        \code{labels}, and \code{colors} to subdivide the promoter region into
#'        multiple bins. \code{breaks} must be a numeric vector in ascending
#'        order (from 5' to 3'), defining the boundaries of promoter sub-regions.
#'        \code{labels} and \code{colors} must have length equal to
#'        \code{length(breaks) - 1}. The precedence follows 3' -> 5' order
#'        (closer to TSS takes precedence). For example:
#'        \code{list(breaks = c(-2000, -1000, -500, 0, 100),
#'        labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", "upstream <500b",
#'        "TSS - 100b"), colors = c("#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32"))}
#' @param geneDownstream A named numeric vector with elements \code{upstream}
#'        and \code{downstream} specifying the gene downstream region boundaries
#'        relative to transcript end sites (TES). Default is
#'        \code{c(upstream = 0, downstream = 1000)}, meaning from TES to 1kb
#'        downstream.
#' @param labels A named list of named character vectors specifying custom labels
#'        for genomic elements. Each element corresponds to a category level:
#'        \itemize{
#'          \item \code{geneLevel}: Labels for promoter, geneDownstream,
#'                geneBody, distalIntergenic
#'          \item \code{ExonIntron}: Labels for exon, intron, intergenic
#'          \item \code{Exons}: Labels for utr5, utr3, CDS, otherExon
#'          \item \code{group}: Labels for category group names (geneLevel,
#'                promoterLevel, Exons, ExonIntron)
#'        }
#'        The order of elements in each vector determines precedence when peaks
#'        overlap multiple categories. Default labels are provided if not specified.
#' @param labelColors A named character vector specifying colors for each
#'        genomic element type. Names should match the element types (e.g.,
#'        "promoter", "geneBody", "exon", "intron", "utr5", "utr3", "CDS", etc.).
#'        Colors not specified will use default colors. Default colors follow
#'        a colorblind-friendly palette.
#' @param plot A logical value. When \code{TRUE} (default), the plot is displayed.
#'        When \code{FALSE}, only the categorized peaks and plot object are
#'        returned without displaying the plot.
#' @param keepExonsInGenesOnly A logical value. When \code{TRUE} (default), only
#'        exons and introns that are within annotated gene boundaries are kept.
#'        Exons/introns outside gene boundaries are filtered out. When
#'        \code{FALSE}, all exons and introns from the TxDb are used, which may
#'        include orphan exons/introns not associated with any gene.
#' @return Returns invisibly a list containing:
#'        \itemize{
#'          \item \code{peaks}: The input peaks with added metadata columns
#'                containing the assigned genomic element categories. For each
#'                category level (geneLevel, ExonIntron, Exons, and optionally
#'                promoterLevel), a column named after the level contains the
#'                assigned element type (e.g., "promoter", "exon", "utr5", etc.).
#'                Peaks that don't overlap any annotated element are labeled as
#'                "undefined".
#'          \item \code{plot}: A ggplot2 object containing the visualization.
#'                For a single GRanges object, this is a donut/pie chart. For a
#'                GRangesList, this is a faceted bar plot with separate panels
#'                for each category level.
#'        }
#'        The plot can be further customized using standard ggplot2 functions.
#' 
#' @export
#' @importFrom ggplot2 ggplot geom_rect xlim coord_polar aes_string geom_bar
#' coord_flip scale_fill_manual theme_void theme_bw facet_wrap geom_col 
#' geom_text guide_legend
#' @importFrom GenomicFeatures intronsByTranscript exons fiveUTRsByTranscript 
#' threeUTRsByTranscript genes transcripts cds
#' @importFrom S4Vectors DataFrame
#' @importFrom stats as.formula
#' @examples 
#' \dontrun{
#'   data(myPeakList)
#'   if (require(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
#'     seqinfo(myPeakList) <- 
#'       seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)[seqlevels(myPeakList)]
#'     myPeakList <- GenomicRanges::trim(myPeakList)
#'     myPeakList <- myPeakList[width(myPeakList) > 0]
#'     
#'     ## Basic usage: peak-centric counting
#'     result <- genomicElementDistribution(myPeakList, 
#'                                          TxDb.Hsapiens.UCSC.hg19.knownGene)
#'     ## Access categorized peaks
#'     peaks_with_annotation <- result$peaks
#'     
#'     ## Nucleotide-centric counting (accounts for peak width)
#'     genomicElementDistribution(myPeakList, 
#'                                TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                nucleotideLevel = TRUE)
#'     
#'     ## Custom promoter levels with multiple bins
#'     ## Breaks are from 5' -> 3', precedence follows 3' -> 5'
#'     genomicElementDistribution(myPeakList, 
#'                                TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                promoterLevel = list(
#'                                  breaks = c(-2000, -1000, -500, 0, 100),
#'                                  labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
#'                                             "upstream <500b", "TSS - 100b"),
#'                                  colors = c("#FFE5CC", "#FFCA99", 
#'                                             "#FFAD65", "#FF8E32")))
#'     
#'     ## Compare multiple peak sets using GRangesList
#'     peaks_list <- GRangesList(set1 = myPeakList, set2 = myPeakList)
#'     genomicElementDistribution(peaks_list, 
#'                                TxDb.Hsapiens.UCSC.hg19.knownGene)
#'   }
#' }
genomicElementDistribution <- 
  function(peaks, TxDb, seqlev, nucleotideLevel = FALSE, ignore.strand = TRUE,
           promoterRegion = c(upstream = 2000, downstream = 100),
           geneDownstream = c(upstream = 0, downstream = 1000),
           labels = list(geneLevel = c(promoter = "Promoter",
                                   geneDownstream = "Downstream",
                                   geneBody = "Gene body",
                                   distalIntergenic = "Distal Intergenic"),
                       ExonIntron = c(exon = "Exon",
                                    intron = "Intron",
                                    intergenic = "Intergenic"),
                       Exons = c(utr5 = "5' UTR",
                               utr3 = "3' UTR",
                               CDS = "CDS",
                               otherExon = "Other exon"),
                       group = c(geneLevel = "Transcript Level",
                               promoterLevel = "Promoter Level",
                               Exons = "Exon level",
                               ExonIntron = "Exon/Intron/Intergenic")),
           labelColors = c(promoter = "#E1F114",
                           geneBody = "#9EFF00",
                           geneDownstream = "#57CB1B",
                           distalIntergenic = "#066A4B",
                           exon = "#6600FF",
                           intron = "#8F00FF",
                           intergenic = "#DA00FF",
                           utr5 = "#00FFDB",
                           utr3 = "#00DFFF",
                           CDS = "#00A0FF",
                           otherExon = "#006FFF"),
           plot = TRUE,
           keepExonsInGenesOnly = TRUE,
           promoterLevel) {
    stopifnot("peaks must be an object of GRanges or GRangesList" =
                inherits(peaks, c("GRanges", "GRangesList")))
    if (inherits(peaks, "GRanges")) {
        n <- deparse(substitute(peaks))
        peaks <- GRangesList(peaks)
        names(peaks) <- n
        isGRanges <- TRUE
    } else {
        isGRanges <- FALSE
    }
    stopifnot("TxDb must be an object of TxDb" = inherits(TxDb, "TxDb"))
    stopifnot("nuleotideLevel is not logical" = is.logical(nucleotideLevel))
    stopifnot("promoterRegion should contain element upstream and downstream" =
                all(c("upstream", "downstream") %in% names(promoterRegion)))
    stopifnot("geneDownstream should contain element upstream and downstream" =
                all(c("upstream", "downstream") %in% names(geneDownstream)))
    stopifnot("Elements in promoterRegion should be numeric" =
                is.numeric(promoterRegion))
    stopifnot("Elements in geneDownstream should be numeric" =
                is.numeric(geneDownstream))
    labs <- list(geneLevel = c(promoter = "Promoter",
                             geneDownstream = "Downstream",
                             geneBody = "Gene body",
                             distalIntergenic = "Distal Intergenic"),
                 ExonIntron = c(exon = "Exon",
                              intron = "Intron",
                              intergenic = "Intergenic"),
                 Exons = c(utr5 = "5' UTR",
                         utr3 = "3' UTR",
                         CDS = "CDS",
                         otherExon = "Other exon"))
    groupLabels <- c(geneLevel = "Gene Level",
                     promoterLevel = "Promoter Level",
                     Exons = "Exon level",
                     ExonIntron = "Exon/Intron/Intergenic")
    labelCols <- c(promoter = "#D55E00",
                  geneDownstream = "#E69F00",
                  geneBody = "#51C6E6",
                  distalIntergenic = "#AAAAAA",
                  exon = "#009DDA",
                  intron = "#666666",
                  intergenic = "#DDDDDD",
                  utr5 = "#0072B2",
                  utr3 = "#56B4E9",
                  CDS = "#0033BF",
                  otherExon = "#009E73",
                  undefined = "#FFFFFF")
    labelCols[names(labelColors)] <- labelColors
    for (i in names(labs)) {
        if (i %in% names(labels)) {
            ## keep the orders in labels
            n <- intersect(names(labs[[i]]), names(labels[[i]]))
            labs[[i]] <- c(labels[[i]][n], labs[[i]][!names(labs[[i]]) %in% n])
        }
    }
    
    for (i in names(groupLabels)) {
        if ("group" %in% names(labels)) {
            groupLabels[names(labels[["group"]])] <- labels[["group"]]
        }
    }
    if (!missing(promoterLevel)) {
        # stopifnot("promoterLevel must be within promoterRegion" =
        #             all(abs(promoterLevel$breaks[promoterLevel$breaks < 0]) <=
        #                   promoterRegion["upstream"]) && 
        #             all(abs(promoterLevel$breaks[promoterLevel$breaks > 0]) <=
        #                   promoterRegion["downstream"]))
        promoterLevel$breaks <- sort(promoterLevel$breaks)
        stopifnot("breaks, labels and colors of promoterLevel are not paired" =
                  length(promoterLevel$breaks) ==
                  length(promoterLevel$labels) + 1 &&
                  length(promoterLevel$labels) ==
                  length(promoterLevel$colors))
        proK <- paste0("promoter", seq_along(promoterLevel$labels))
        promoterLevel$upstream <- ifelse(promoterLevel$breaks < 0,
                                         abs(promoterLevel$breaks),
                                         0)
        promoterLevel$upstream <- 
            promoterLevel$upstream[-length(promoterLevel$upstream)]
        promoterLevel$downstream <- ifelse(promoterLevel$breaks > 0,
                                          promoterLevel$breaks,
                                          0)
        promoterLevel$downstream <- promoterLevel$downstream[-1]
        proV <- promoterLevel$labels
        names(proV) <- proK
        labs <- c(list("promoterLevel" = proV), 
                  labs)
        proV <- promoterLevel$colors
        names(proV) <- proK
        labelCols <- c(proV, labelCols)
    } else {
        promoterLevel <- NULL
    }
    
    defaultW <- getOption("warn")
    options(warn = -1)
    on.exit(options(warn = defaultW))
    
    ## Set annotation
    anno <- GRangesList()
    for (i in names(labs)) {
        anno[[i]] <- switch (i,
        "promoterLevel" = {
            suppressMessages(g <- transcripts(TxDb))
            pro <- mapply(FUN = function(upstream, downstream) {
                unique(promoters(g, upstream = upstream, downstream = downstream))
            }, promoterLevel$upstream, 
            promoterLevel$downstream,
            SIMPLIFY = FALSE)
            names(pro) <- names(labs[["promoterLevel"]])
            pro <- rev(pro)
            current_anno <- GRanges()
            for (j in seq_along(pro)) {
                ca <- filterByOverlaps(pro[[j]], current_anno,
                                       ignore.strand = ignore.strand)
                mcols(ca) <- DataFrame(type = rep(names(pro)[j], length(ca)))
                current_anno <- c(current_anno, ca)
            }
            current_anno <- formatSeqnames(current_anno, peaks)
            current_anno
        },
        "geneLevel" = {
            suppressMessages(g <- transcripts(TxDb))
            pro <- unique(promoters(g, 
                             upstream = promoterRegion["upstream"],
                             downstream = promoterRegion["downstream"]))
            dws <- downstreams(g,
                               upstream = geneDownstream["upstream"],
                               downstream = geneDownstream["downstream"])
            pro <- GenomicRanges::trim(pro)
            dws <- GenomicRanges::trim(dws)
            intergenic <- gaps(reduce(c(pro, g, dws), ignore.strand = FALSE))
            intergenic <- intergenic[!strand(intergenic) %in% "*"]
            current_anno <- GRanges()
            ## set precedence
            for (j in names(labs[["geneLevel"]])) {
                s <- 
                    switch(j,
                           "promoter" = pro,
                           "geneDownstream" = dws,
                           "geneBody" = g,
                           "distalIntergenic" = intergenic)
                ca <- 
                    filterByOverlaps(s, current_anno,
                                   ignore.strand = ignore.strand)
                mcols(ca) <- DataFrame(type = rep(j, length(ca)))
                current_anno <- c(current_anno, ca)
            }
            current_anno <- formatSeqnames(current_anno, peaks)
            current_anno
        },
        "ExonIntron" = {
            exon <- exons(TxDb)
            intron <- unlist(intronsByTranscript(TxDb))
            if (keepExonsInGenesOnly) {
                suppressMessages(g <- transcripts(TxDb))
                ole <- findOverlaps(exon, g, type = "within")
                oli <- findOverlaps(intron, g, type = "within")
                ole <- !seq_along(exon) %in% queryHits(ole)
                oli <- !seq_along(intron) %in% queryHits(oli)
                if (sum(ole) > 0 || sum(oli) > 0) {
                    warning(paste(sum(ole), "exons were dropped because there is no",
                                "relative gene level annotations.",
                                sum(oli), "introns were dropped because there is no",
                                "relative gene level annotations."))
                    exon <- exon[!ole]
                    intron <- intron[!oli]
                }
            }
            intergenic <- gaps(reduce(c(exon, intron), ignore.strand = FALSE))
            intergenic <- intergenic[!strand(intergenic) %in% "*"] 
            current_anno <- GRanges()
            ## set precedence
            for (j in names(labs[["ExonIntron"]])) {
                s <- 
                    switch(j,
                           "exon" = exon,
                           "intron" = intron,
                           "intergenic" = intergenic)
                ca <- 
                    filterByOverlaps(s, current_anno,
                                   ignore.strand = ignore.strand)
                mcols(ca) <- DataFrame(type = rep(j, length(ca)))
                current_anno <- c(current_anno, ca)
            }
            current_anno <- formatSeqnames(current_anno, peaks)
            current_anno
        },
        "Exons" = {
            utr5 <- unlist(fiveUTRsByTranscript(TxDb))
            utr3 <- unlist(threeUTRsByTranscript(TxDb))
            CDS <- cds(TxDb)
            exon <- exons(TxDb)
            if (keepExonsInGenesOnly) {
                suppressMessages(g <- transcripts(TxDb))
                ole <- findOverlaps(exon, g, type = "within")
                olc <- findOverlaps(CDS, g, type = "within")
                ol5 <- findOverlaps(utr5, g, type = "within")
                ol3 <- findOverlaps(utr3, g, type = "within")
                ole <- !seq_along(exon) %in% queryHits(ole)
                olc <- !seq_along(CDS) %in% queryHits(olc)
                ol5 <- !seq_along(utr5) %in% queryHits(ol5)
                ol3 <- !seq_along(utr3) %in% queryHits(ol3)
                if (sum(ole) > 0 || sum(olc) > 0 || sum(ol5) > 0 || sum(ol3)) {
                    warning(paste(sum(ole), "exons were dropped because there is no",
                                "relative gene level annotations.",
                                sum(olc), "CDS were dropped because there is no",
                                "relative gene level annotations.",
                                sum(ol5), "utr5 were dropped because there is no",
                                "relative gene level annotations.",
                                sum(ol3), "utr3 were dropped because there is no",
                                "relative gene level annotations."))
                    exon <- exon[!ole]
                    CDS <- CDS[!olc]
                    utr5 <- utr5[!ol5]
                    utr3 <- utr3[!ol3]
                }
            }
            current_anno <- GRanges()
            ## set precedence
            for (j in names(labs[["Exons"]])) {
                s <- 
                    switch(j,
                           "utr5" = utr5,
                           "utr3" = utr3,
                           "CDS" = CDS,
                           "otherExon" = exon)
                ca <- 
                    filterByOverlaps(s, current_anno,
                                   ignore.strand = ignore.strand)
                mcols(ca) <- DataFrame(type = rep(j, length(ca)))
                current_anno <- c(current_anno, ca)
            }
            current_anno <- formatSeqnames(current_anno, peaks)
            current_anno
        }
        )
    }
    groupLabels <- groupLabels[names(anno)]
    groupLabels[is.na(groupLabels)] <- names(anno)[is.na(groupLabels)]
    
    ## filter peaks by seqlev
    if (!missing(seqlev)) {
        if (length(seqlev) > 0) {
            peaks <- lapply(peaks, 
                          function(.ele) .ele[seqnames(.ele) %in% seqlev])
        }
    }
    
    if (nucleotideLevel) {
        peaks <- lapply(peaks, function(.ele) {
            y <- disjoin(c(.ele, unlist(anno)), ignore.strand = ignore.strand)
            subsetByOverlaps(y, .ele, ignore.strand = ignore.strand)
        })
    }
    
    peaks <- lapply(peaks, FUN = function(.peaks) {
        pct <- lapply(anno, FUN = function(.ele) {
            y <- .peaks
            ol <- findOverlaps(y, .ele, ignore.strand = ignore.strand)
            ol <- as.data.frame(ol)
            ol <- ol[order(ol$queryHits, ol$subjectHits), ]
            ol <- ol[!duplicated((ol$queryHits)), ]
            mcols(y)[, 'anno'] <- rep("undefined", length(y))
            y$anno[ol$queryHits] <- .ele$type[ol$subjectHits]
            y$anno
        })
        
        pct <- do.call(cbind, pct)
        if (ncol(mcols(.peaks))) {
            mcols(.peaks) <- cbind(mcols(.peaks), pct)
        } else {
            mcols(.peaks) <- pct
        }
        
        .peaks
    })
    if (isGRanges) {
        peaks <- peaks[[1]]
    }
    
    melt <- function(m) {
        if (nucleotideLevel) {
            m <- m[rep(seq_along(m), width(m))]
        }
        m <- mcols(m)[, names(anno)]
        p <- lapply(names(anno), function(.ele) {
            tt <- table(m[, .ele])
            data.frame(category = factor(rep(groupLabels[.ele],
                                           length(tt)), 
                                       levels = groupLabels),
                       type = factor(names(tt), levels = rev(names(labelCols))), 
                       percentage = as.numeric(tt) / sum(tt))
        })
        do.call(rbind, p)
    }
    
    if (isGRanges) { ## donut-plot
        #dat <- reshape2::melt()
        dat <- melt(peaks)
        l <- unlist(unname(labs))
        l1 <- paste0(l, " (", 
                     round(dat$percentage[match(names(l), dat$type)] * 100,
                           digits = 1), "%)")
        names(l1) <- names(l)
        l1 <- c(l1, undefined = "")
        p <- ggplot(dat, 
                    aes_string(x = "category", y = "percentage", 
                             fill = "type")) +
            geom_col() +
            coord_polar("y") + 
            geom_text(data = subset(dat, !duplicated(dat$category)),
                      aes_string(x = "category", label = "category"),
                      y = 1) +
            theme_void()

    } else { ## bar-plot
        dat <- lapply(peaks, melt)
        dat1 <- do.call(rbind, dat)
        dat1$source <- rep(names(peaks), vapply(dat, nrow, FUN.VALUE = 0))
        l1 <- c(unlist(unname(labs)), undefined = "")
        p <- ggplot(dat1, 
                    aes_string(x = "source", y = "percentage", fill = "type")) +
            geom_bar(stat = "identity") + coord_flip() +
            facet_wrap(as.formula("~ category"), ncol = 1) + 
            theme_bw()
    }
    p <- p + 
        scale_fill_manual(values = labelCols, labels = l1, name = NULL,
                          guide = guide_legend(reverse = TRUE))
    
    if (plot) {
        print(p)
    }
    return(invisible(list(peaks = peaks, plot = p)))
  }
