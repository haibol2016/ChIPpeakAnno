#' Prepare data for UpSet plot of genomic element distribution
#' 
#' @description 
#' Prepares data for UpSet plots showing the distribution of peaks across
#' multiple genomic element categories. Unlike \code{\link{genomicElementDistribution}},
#' this function does not apply precedence rules, allowing peaks to be counted
#' in multiple categories simultaneously. A peak can overlap with multiple
#' genomic elements and will be counted in all relevant categories. The output
#' can be used with the \code{UpSetR} package to create intersection plots
#' showing how peaks are distributed across different genomic element combinations.
#' 
#' @details 
#' This function annotates peaks by their overlap with various genomic elements
#' defined in the \code{breaks} parameter. For each peak, it determines which
#' genomic element categories it overlaps with, creating a binary matrix where
#' each row represents a peak and each column represents a genomic element category.
#' A value of 1 indicates overlap, 0 indicates no overlap.
#' 
#' \strong{Key differences from \code{genomicElementDistribution}:}
#' \itemize{
#'   \item No precedence rules: peaks can be counted in multiple categories
#'   \item Returns data in UpSetR-compatible format
#'   \item Supports both single GRanges and GRangesList inputs
#' }
#' 
#' \strong{Genomic element definition via \code{breaks}:}
#' Each element in the \code{breaks} list can be either:
#' \itemize{
#'   \item A function that takes a \code{TxDb} object and returns a
#'         \code{GRanges} or \code{GRangesList} (e.g., \code{genes},
#'         \code{exons}, \code{cds}, \code{fiveUTRsByTranscript})
#'   \item A numeric vector of length 4: \code{c(upstream_point, downstream_point,
#'         direction, remove_gene_body)}
#'         \itemize{
#'           \item \code{upstream_point}: Upstream boundary (negative for upstream,
#'                 positive for downstream)
#'           \item \code{downstream_point}: Downstream boundary (negative for upstream,
#'                 positive for downstream)
#'           \item \code{direction}: \code{-1} for promoter/upstream regions,
#'                 \code{1} for downstream regions
#'           \item \code{remove_gene_body}: \code{1} to exclude gene body overlaps,
#'                 \code{0} to keep gene body overlaps
#'         }
#' }
#' 
#' @param peaks A \link[GenomicRanges:GRanges-class]{GRanges} object or
#'        \link[GenomicRanges:GRangesList-class]{GRangesList} containing peaks
#'        to be annotated. If a single \code{GRanges} is provided, it will be
#'        converted to a \code{GRangesList} internally.
#' @param TxDb An object of class \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#'        containing transcript annotation data.
#' @param seqlev A character vector specifying which sequence levels (chromosomes)
#'        should be included in the analysis. If missing, all sequence levels
#'        present in both \code{peaks} and \code{TxDb} will be used.
#' @param ignore.strand A logical value. When \code{TRUE} (default), strand
#'        information is ignored when determining overlaps between peaks and
#'        genomic elements. When \code{FALSE}, only ranges on the same strand
#'        are considered for overlap.
#' @param breaks A named list defining genomic element categories. Each element
#'        can be either:
#'        \itemize{
#'          \item A function that extracts genomic features from \code{TxDb}
#'                (e.g., \code{genes}, \code{exons}, \code{cds},
#'                \code{fiveUTRsByTranscript}, \code{threeUTRsByTranscript},
#'                \code{intronsByTranscript})
#'          \item A numeric vector of length 4 defining a region relative to genes
#'        }
#'        Default categories include:
#'        \itemize{
#'          \item \code{"distal_upstream"}: 100kb to 10kb upstream
#'          \item \code{"proximal_upstream"}: 10kb to 5kb upstream
#'          \item \code{"distal_promoter"}: 5kb to 2kb upstream
#'          \item \code{"proximal_promoter"}: 2kb upstream to 200bp downstream
#'          \item \code{"5'UTR"}: 5' untranslated regions
#'          \item \code{"3'UTR"}: 3' untranslated regions
#'          \item \code{"CDS"}: Coding sequences
#'          \item \code{"exon"}: Exons
#'          \item \code{"intron"}: Introns
#'          \item \code{"gene_body"}: Entire gene regions
#'          \item \code{"immediate_downstream"}: 0 to 2kb downstream
#'          \item \code{"proximal_downstream"}: 2kb to 5kb downstream
#'          \item \code{"distal_downstream"}: 5kb to 100kb downstream
#'        }
#' 
#' @return Returns a list with two components:
#'        \itemize{
#'          \item \code{peaks}: The input peaks with additional metadata columns:
#'                \itemize{
#'                  \item \code{id}: Original peak index
#'                  \item \code{annoType}: Genomic element type(s) the peak overlaps
#'                        with. Peaks overlapping multiple elements will have
#'                        multiple entries (one per overlap). Peaks with no
#'                        overlaps will have \code{annoType = "undefined"}
#'                }
#'                If input was a \code{GRangesList}, this will be a
#'                \code{GRangesList}. If input was a single \code{GRanges},
#'                this will be a \code{GRanges} object.
#'          \item \code{plotData}: A binary matrix (data frame) suitable for
#'                \code{UpSetR::upset()}. Each row represents a peak, each column
#'                represents a genomic element category. Values are 1 (overlap) or
#'                0 (no overlap). If input was a \code{GRangesList}, this will be
#'                a list of data frames (one per peak set). If input was a single
#'                \code{GRanges}, this will be a single data frame.
#'        }
#' 
#' @export
#' @importFrom GenomicFeatures intronsByTranscript exons fiveUTRsByTranscript 
#' threeUTRsByTranscript genes cds promoters
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @seealso \code{\link{genomicElementDistribution}} for precedence-based
#'          annotation, \code{\link[UpSetR:upset]{upset}} for creating UpSet plots
#' @keywords misc
#' @examples 
#' \dontrun{
#'   ## Example 1: Single peak set
#'   data(myPeakList)
#'   library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'   seqinfo(myPeakList) <- 
#'     seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)[seqlevels(myPeakList)]
#'   myPeakList <- GenomicRanges::trim(myPeakList)
#'   myPeakList <- myPeakList[width(myPeakList) > 0]
#'   
#'   x <- genomicElementUpSetR(myPeakList, 
#'                             TxDb.Hsapiens.UCSC.hg19.knownGene)
#'   
#'   ## Create UpSet plot
#'   library(UpSetR)
#'   upset(x$plotData, nsets = 13, nintersects = NA)
#'   
#'   ## Example 2: Multiple peak sets (GRangesList)
#'   peaks1 <- myPeakList[1:100]
#'   peaks2 <- myPeakList[101:200]
#'   peaks_list <- GRangesList(peaks1 = peaks1, peaks2 = peaks2)
#'   
#'   x2 <- genomicElementUpSetR(peaks_list, 
#'                              TxDb.Hsapiens.UCSC.hg19.knownGene)
#'   
#'   ## Plot for first peak set
#'   upset(x2$plotData[[1]], nsets = 13)
#'   
#'   ## Example 3: Custom breaks
#'   custom_breaks <- list(
#'     "promoter" = c(-2000, 200, -1, 0),
#'     "gene_body" = genes,
#'     "downstream" = c(0, 5000, 1, 1)
#'   )
#'   x3 <- genomicElementUpSetR(myPeakList, 
#'                                TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                breaks = custom_breaks)
#'   upset(x3$plotData, nsets = 3)
#' }
genomicElementUpSetR <- 
  function(peaks, TxDb, seqlev, ignore.strand = TRUE,
           breaks = 
             list("distal_upstream" = c(-100000, -10000, -1, 1), 
                  "proximal_upstream" = c(-10000, -5000, -1, 1),
                  "distal_promoter" = c(-5000, -2000, -1, 1),
                  "proximal_promoter" = c(-2000, 200, -1, 0),
                  "5'UTR" = fiveUTRsByTranscript,
                  "3'UTR" = threeUTRsByTranscript,
                  "CDS" = cds,
                  "exon" = exons,
                  "intron" = intronsByTranscript,
                  "gene_body" = genes,
                  "immediate_downstream" = c(0, 2000, 1, 1),
                  "proximal_downstream" = c(2000, 5000, 1, 1),
                  "distal_downstream" = c(5000, 100000, 1, 1))) {
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
    
    defaultW <- getOption("warn")
    options(warn = -1)
    on.exit(options(warn = defaultW))
    
    seql <- seqlevelsStyle(peaks)
    
    ## set annotation
    suppressMessages(g <- genes(TxDb, single.strand.genes.only = TRUE))
    
    anno <- lapply(breaks, function(.ele) {
        stopifnot("Elements of breaks must be function or numeric(4)" =
                  is.function(.ele) || is.numeric(.ele))
        if (is.numeric(.ele)) {
            stopifnot("Elements of breaks must be function or numeric(4)" =
                      length(.ele) == 4)
            ## upstream or downstream
            ups_dws <- function(n) {
                ifelse(n[3] == -1,
                       ifelse(all(n[1:2] < 0), "uu", 
                              ifelse(all(n[1:2] > 0), "pp", "ups")),
                       ifelse(all(n[1:2] > 0), "dd", 
                              ifelse(all(n[1:2] < 0), "ww", "dws")))
            }
            x <- ups_dws(.ele)
            fil <- .ele[4]
            .ele <- .ele[1:2]
            .ele <- switch(x,
                          "uu" = {
                              a <- promoters(g, upstream = abs(min(.ele)), 
                                           downstream = 0)
                              b <- promoters(g, upstream = abs(max(.ele)),
                                           downstream = 0)
                              filterByOverlaps(a, b,
                                               ignore.strand = ignore.strand)
                          },
                          "ups" = {
                              promoters(g, upstream = abs(min(.ele)), 
                                      downstream = abs(max(.ele)))
                          },
                          "pp" = {
                              a <- promoters(g, upstream = 0, 
                                           downstream = max(.ele))
                              b <- promoters(g, upstream = 0,
                                           downstream = min(.ele))
                              filterByOverlaps(a, b,
                                               ignore.strand = ignore.strand)
                          },
                          "dd" = {
                              a <- downstreams(g, upstream = 0, 
                                             downstream = max(.ele))
                              b <- downstreams(g, upstream = 0,
                                             downstream = min(.ele))
                              filterByOverlaps(a, b,
                                               ignore.strand = ignore.strand)
                          },
                          "dws" = {
                              downstreams(g, upstream = abs(min(.ele)),
                                        downstream = abs(max(.ele)))
                          },
                          "ww" = {
                              a <- downstreams(g, upstream = abs(min(.ele)), 
                                             downstream = 0)
                              b <- downstreams(g, upstream = abs(max(.ele)),
                                             downstream = 0)
                              filterByOverlaps(a, b,
                                               ignore.strand = ignore.strand)
                          })
            if (fil == 1) {
                .ele <- filterByOverlaps(.ele, g, ignore.strand = ignore.strand)
            }
            seqlevelsStyle(.ele) <- seql[1]
            return(.ele)
        }
        if (is.function(.ele)) {
            .ele <- .ele(TxDb)
            if (inherits(.ele, 'GRangesList')) .ele <- unlist(.ele)
            seqlevelsStyle(.ele) <- seql[1]
            return(.ele)
        }
    })
    l <- lengths(anno)
    n <- names(anno)
    anno <- unlist(GRangesList(anno))
    mcols(anno) <- DataFrame(type = rep(n, l))
    
    ## filter peaks by seqlev
    if (!missing(seqlev)) {
        if (length(seqlev) > 0) {
            peaks <- lapply(peaks, 
                          function(.ele) .ele[seqnames(.ele) %in% seqlev])
        }
    }
    
    peaks <- lapply(peaks, FUN = function(.peaks) {
        .peaks$id <- seq_along(.peaks)
        ol <- findOverlaps(.peaks, anno, ignore.strand = ignore.strand)
        .peaks$annoType <- "undefined"
        .p1 <- .peaks[-unique(queryHits(ol))]
        .p2 <- .peaks[queryHits(ol)]
        .p2$annoType <- anno$type[subjectHits(ol)]
        .peaks <- c(.p1, .p2)
        .peaks[order(.peaks$id)]
    })
    if (isGRanges) {
        peaks <- peaks[[1]]
    }
    melt <- function(.ele) {
        .i <- unique(.ele$id)
        .j <- c(names(breaks), "undefined")
        x <- matrix(0, nrow = length(.i), ncol = length(.j))
        colnames(x) <- .j
        rownames(x) <- .i
        .ele <- split(.ele$id, .ele$annoType)
        for (i in seq_along(.ele)) {
            x[.ele[[i]], names(.ele)[i]] <- 1
        }
        as.data.frame(x)
    }
    if (isGRanges) { ## for upset
        dat <- melt(peaks)
    } else { ## bar-plot
        dat <- lapply(peaks, melt)
    }
    
    return(list(peaks = peaks, plotData = dat))
  }
