#' Find possible enhancers using DNA interaction data
#' 
#' @description 
#' Identifies potential enhancer regions by integrating ChIP-seq peak data with
#' chromosome conformation capture (3C, 4C, 5C, Hi-C) interaction data. The
#' function uses DNA interaction information to shift gene annotation coordinates,
#' allowing identification of peaks that may be associated with distal regulatory
#' elements (enhancers) through chromatin looping, even when they are far from
#' genes in linear genomic coordinates.
#' 
#' This function is particularly useful for identifying long-range enhancer-gene
#' interactions where enhancers are located far from their target genes but are
#' brought into proximity through three-dimensional chromatin structure.
#' 
#' @param peaks A \code{\link[GenomicRanges]{GRanges}} object containing ChIP-seq
#'        peaks to be analyzed for enhancer associations. These are typically
#'        transcription factor binding sites or histone modification peaks that
#'        may function as enhancers.
#' @param annoData A \code{\link[GenomicRanges]{GRanges}} or \code{annoGR} object
#'        containing gene annotations (e.g., TSS, gene bodies, promoters). The
#'        annotations will be shifted based on DNA interaction regions to find
#'        potential enhancer-target gene associations.
#' @param DNAinteractiveData A \code{\link[GenomicRanges]{GRanges}} object with
#'        \code{blocks} metadata (for BED format), a \code{GInteractions} object
#'        (from InteractionSet package), or a file path to a BEDPE file
#'        containing DNA interaction data from Hi-C, 3C, 4C, 5C, or similar
#'        experiments. The interaction data defines regions that are in physical
#'        proximity in 3D space.
#' @param bindingType A character string specifying the criteria for associating
#'        peaks with shifted annotations. Must be one of:
#'        \itemize{
#'          \item \code{"nearestBiDirectionalPromoters"} (default): Finds the
#'                nearest enhancer regions from both directions relative to the
#'                shifted annotation positions. This is the most comprehensive
#'                option.
#'          \item \code{"startSite"}: Associates peaks with shifted TSS
#'                (transcription start site) positions. Useful for finding
#'                enhancers that regulate gene transcription initiation.
#'          \item \code{"endSite"}: Associates peaks with shifted gene/exon end
#'                positions. Useful for finding enhancers that may affect gene
#'                termination or 3' end processing.
#'        }
#' @param bindingRegion An integer vector of length 2 specifying the annotation
#'        range relative to shifted feature positions. The first value must be <=
#'        0 (upstream offset), and the second value must be >= 1 (downstream
#'        offset). Default is \code{c(-5000, 5000)}, meaning the function will
#'        search for peaks within 5 kb upstream and 5 kb downstream of the
#'        shifted annotation positions.
#' @param ignore.peak.strand A logical value. If \code{TRUE} (default), ignores
#'        peak strand information when finding overlaps. If \code{FALSE}, strand
#'        information is considered, which may be important for strand-specific
#'        analyses.
#' @param ... Additional arguments (currently not used).
#' 
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object containing
#'        annotated peaks with extensive metadata columns indicating:
#'        \itemize{
#'          \item \code{feature}: Index to the original annotation that was
#'                shifted
#'          \item \code{feature.shift.*}: Annotation information after shifting
#'                based on DNA interactions (e.g., \code{feature.shift.start},
#'                \code{feature.shift.end}, \code{feature.shift.strand})
#'          \item \code{feature.ranges}: Original annotation ranges (before
#'                shifting)
#'          \item \code{feature.strand}: Original annotation strand
#'          \item \code{distance}: Distance from peak to shifted annotation
#'          \item \code{DNAinteractive.*}: Information about the DNA interaction
#'                region (e.g., \code{DNAinteractive.ranges},
#'                \code{DNAinteractive.blocks})
#'          \item Standard peak annotation columns (peak, start, end, etc.)
#'        }
#'        The returned object contains only peaks that:
#'        \itemize{
#'          \item Overlap with DNA interaction regions
#'          \item Can be associated with annotations through shifted coordinates
#'          \item Pass the binding region and binding type criteria
#'        }
#'        If no peaks meet these criteria, an empty GRanges object is returned.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' The function implements a complex workflow to identify enhancers through
#' chromatin interaction data:
#' \enumerate{
#'   \item \strong{Process DNA interaction data}: Converts input to interaction
#'         regions with two anchor points (A-B and C-D, representing the two
#'         interacting regions)
#'   \item \strong{Define interaction regions}: Creates five regions around each
#'         interaction: A (upstream), A-B (first anchor), B-C (between anchors),
#'         C-D (second anchor), and D (downstream)
#'   \item \strong{Find peak overlaps}: Identifies which peaks overlap with each
#'         of the five interaction regions
#'   \item \strong{Find annotation overlaps}: Identifies which annotations
#'         overlap with each interaction region
#'   \item \strong{Shift annotations}: For each interaction, shifts annotation
#'         coordinates to the corresponding position in the interacting region
#'         (e.g., if a gene is in region A, shift it to region C or D)
#'   \item \strong{Annotate peaks}: Uses \code{\link{annoPeaks}} to find peaks
#'         near the shifted annotation positions
#'   \item \strong{Filter and return}: Returns peaks that are associated with
#'         annotations through the shifted coordinates
#' }
#' 
#' \strong{DNA interaction data formats:}
#' The function supports multiple input formats:
#' \itemize{
#'   \item \strong{GRanges with blocks}: BED format with \code{blocks} metadata
#'         containing two intervals per range (representing the two interacting
#'         regions)
#'   \item \strong{GInteractions}: InteractionSet package format for storing
#'         pairwise genomic interactions
#'   \item \strong{File path}: Path to a BEDPE file (will be converted to
#'         GInteractions format)
#' }
#' 
#' \strong{Coordinate shifting logic:}
#' The function shifts annotations based on interaction regions. For example:
#' \itemize{
#'   \item If a gene is in region A and a peak is in region C, the gene's
#'         coordinates are shifted to region C to check if the peak is near the
#'         shifted gene position
#'   \item Multiple shifting strategies are applied (AC, AD, BC, BD) to capture
#'         all possible enhancer-gene associations
#'   \item Reversed coordinates are also considered to handle different
#'         interaction orientations
#' }
#' 
#' \strong{Interaction region definitions:}
#' Each DNA interaction is divided into five regions:
#' \itemize{
#'   \item \code{A.ups}: Upstream of the first anchor (A)
#'   \item \code{AB}: First anchor region (A to B)
#'   \item \code{BC}: Between anchors (B to C)
#'   \item \code{CD}: Second anchor region (C to D)
#'   \item \code{D.dws}: Downstream of the second anchor (D)
#' }
#' Peaks and annotations are mapped to these regions, and coordinates are shifted
#' accordingly to find enhancer-gene associations.
#' 
#' \strong{Binding type options:}
#' \itemize{
#'   \item \code{"nearestBiDirectionalPromoters"}: Most comprehensive, finds
#'         nearest annotations from both directions
#'   \item \code{"startSite"}: Focuses on TSS/promoter regions
#'   \item \code{"endSite"}: Focuses on gene/exon end regions (strand-reversed
#'         for proper handling)
#' }
#' 
#' \strong{Filtering and deduplication:}
#' The function applies several filtering steps:
#' \itemize{
#'   \item Only considers interactions where both peaks and annotations overlap
#'   \item Removes duplicate peak-annotation pairs
#'   \item Sorts by peak ID and distance
#'   \item Returns only peaks that pass all criteria
#' }
#' 
#' @note
#' \itemize{
#'   \item The function requires that peaks, annotations, and interaction data
#'         have compatible seqlevels styles (e.g., all use "UCSC" or all use
#'         "NCBI")
#'   \item Large interaction datasets may require significant memory and
#'         computation time
#'   \item The function returns an empty GRanges object if no peaks can be
#'         associated with annotations through interactions
#'   \item Metadata columns are extensively modified to include interaction
#'         information
#'   \item The function uses \code{annoPeaks} internally for final peak
#'         annotation
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{annotatePeakInBatch}} for standard peak annotation
#'   \item \code{\link{annoPeaks}} for region-based peak annotation (used
#'         internally)
#'   \item \code{\link{toGRanges}} for converting various formats to GRanges
#'   \item InteractionSet package for \code{GInteractions} objects
#' }
#' 
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom InteractionSet GInteractions
#' @importMethodsFrom S4Vectors first second
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom BiocGenerics start end width strand
#' @importFrom S4Vectors elementNROWS mcols queryHits subjectHits
#' @examples
#' 
#' # Example 1: Basic usage with BED format interaction data
#' bed <- system.file("extdata", 
#'                    "wgEncodeUmassDekker5CGm12878PkV2.bed.gz",
#'                    package="ChIPpeakAnno")
#' DNAinteractiveData <- toGRanges(gzfile(bed))
#' library(EnsDb.Hsapiens.v75)
#' annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
#' data("myPeakList")
#' enhancers <- findEnhancers(myPeakList[500:1000], annoData, DNAinteractiveData)
#' 
#' # Example 2: Using different binding type (TSS only)
#' enhancers_tss <- findEnhancers(myPeakList[500:1000], annoData, 
#'                                 DNAinteractiveData,
#'                                 bindingType = "startSite")
#' 
#' # Example 3: Adjusting binding region (wider search window)
#' enhancers_wide <- findEnhancers(myPeakList[500:1000], annoData, 
#'                                  DNAinteractiveData,
#'                                  bindingRegion = c(-10000, 10000))
#' 
#' # Example 4: Using GInteractions object (if available)
#' \dontrun{
#' library(InteractionSet)
#' # Assuming you have a GInteractions object
#' # enhancers <- findEnhancers(peaks, annoData, ginteractions_obj)
#' }
#' 
#' # Example 5: Inspecting results
#' \dontrun{
#' enhancers <- findEnhancers(myPeakList[500:1000], annoData, DNAinteractiveData)
#' # View metadata columns
#' colnames(mcols(enhancers))
#' # Check distances to shifted annotations
#' enhancers$distance
#' # View DNA interaction information
#' enhancers$DNAinteractive.ranges
#' }
#'   
findEnhancers <- function(peaks, annoData, DNAinteractiveData,
                          bindingType = c("nearestBiDirectionalPromoters",
                                          "startSite", "endSite"),
                          bindingRegion = c(-5000, 5000),
                          ignore.peak.strand = TRUE, ...) {
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    stopifnot(length(intersect(seqlevelsStyle(peaks),
                               seqlevelsStyle(annoData))) > 0)
    bindingType <- match.arg(bindingType)
    stopifnot(length(bindingRegion) == 2)
    stopifnot(bindingRegion[1] <= 0 && bindingRegion[2] >= 1)
    
    if (inherits(annoData, "annoGR")) {
        annoData <- as(annoData, "GRanges")
    }
    
    peaks$peak.oid.to.be.deleted <- seq_along(peaks)
    
    if (inherits(DNAinteractiveData, "GRanges")) {
        stopifnot(length(DNAinteractiveData$blocks) > 0)
        stopifnot(all(elementNROWS(DNAinteractiveData$blocks) == 2))
        stopifnot(length(intersect(seqlevelsStyle(peaks),
                                   seqlevelsStyle(DNAinteractiveData))) > 0)
        HiC_FIRST <- lapply(DNAinteractiveData$blocks, `[`, 1)
        HiC_SECOND <- lapply(DNAinteractiveData$blocks, `[`, 2)
        HiC_FIRST <- unlist(IRangesList(HiC_FIRST))
        HiC_SECOND <- unlist(IRangesList(HiC_SECOND))
        ## BED file blocks are half open half close.
        HiC_FIRST <- shift(HiC_FIRST, start(DNAinteractiveData) - 1)
        HiC_SECOND <- shift(HiC_SECOND, start(DNAinteractiveData) - 1)
        HiC_FIRST_GR <- HiC_SECOND_GR <- DNAinteractiveData
        ranges(HiC_FIRST_GR) <- HiC_FIRST
        ranges(HiC_SECOND_GR) <- HiC_SECOND
    }
    ## TODO
    # if(is.character(DNAinteractiveData)){ 
    #     if(DNAinteractiveData %in% c("hg38", "hg19", "mm10", 
    #                                  "danRer10", "danRer11")){
    #         DNAinteractiveData <- 
    #             readRDS(system.file("extdata", 
    #                                 paste0(DNAinteractiveData, 
    #                                        "interactions.rds"),
    #                                 package = "ChIPpeakAnno",
    #                                 mustWork = TRUE))
    #     }else{
    #         DNAinteractiveData <- toGInteractions(DNAinteractiveData)
    #     }
    # }
    if (inherits(DNAinteractiveData, c("GInteractions", "Pairs"))) {
        HiC_FIRST_GR <- first(DNAinteractiveData)
        HiC_SECOND_GR <- second(DNAinteractiveData)
        stopifnot("seqlevels style of peaks and interaction data are different" =
                  length(intersect(seqlevelsStyle(peaks),
                                  seqlevelsStyle(HiC_FIRST_GR))) > 0)
        stopifnot("seqlevels style of peaks and interaction data are different" =
                  length(intersect(seqlevelsStyle(peaks),
                                  seqlevelsStyle(HiC_SECOND_GR))) > 0)
    }
    # peaks overlap with interaction region A_B, C_D
    # in upstream A, A_B, B_C, C_D, downstream D.
    HiC.A.ups <- shift(HiC_FIRST_GR, -max(abs(bindingRegion)))
    width(HiC.A.ups) <- max(abs(bindingRegion))
    HiC.D.dws <- shift(HiC_SECOND_GR, max(abs(bindingRegion)))
    start(HiC.D.dws) <- end(HiC_SECOND_GR) + 1
    HiC.BC <- HiC_SECOND_GR
    start(HiC.BC) <- end(HiC_FIRST_GR) + 1
    end(HiC.BC) <- start(HiC_SECOND_GR) - 1
    HiC.groups <- list(A.ups = HiC.A.ups,
                       AB = HiC_FIRST_GR,
                       BC = HiC.BC,
                       CD = HiC_SECOND_GR,
                       D.dws = HiC.D.dws)
    peaks.ol.HiCdata <- lapply(HiC.groups, function(hic) {
        ol <- findOverlaps(peaks, hic)
        this.peaks <- peaks[queryHits(ol)]
        this.peaks$HiC.idx <- subjectHits(ol)
        this.peaks
    })
    if (all(elementNROWS(peaks.ol.HiCdata) == 0L)) {
        # No peaks in interaction regions
        return(GRanges())
    }
    ## refine annoData by HiCdata
    anno.ol.HiCdata <- lapply(HiC.groups, function(hic) {
        ol <- findOverlaps(annoData, hic)
        annoData.ol.HiC <- annoData[queryHits(ol)]
        annoData.ol.HiC.pos <-
            switch(bindingType,
                   nearestBiDirectionalPromoters = {
                       promoters(annoData.ol.HiC, upstream = 0, downstream = 1)
                   },
                   startSite = {
                       promoters(annoData.ol.HiC, upstream = 0, downstream = 1)
                   },
                   endSite = {
                       tmp <- annoData.ol.HiC
                       strand(tmp) <- ifelse(strand(tmp) == "-", "+", "-")
                       promoters(tmp, upstream = 0, downstream = 1)
                   },
                   stop("Not supported binding type", bindingType))
        HiC_FIRST_GR.anno <- HiC_FIRST_GR[subjectHits(ol)]
        HiC_SECOND_GR.anno <- HiC_SECOND_GR[subjectHits(ol)]
        annoData.ol.HiC$point_A <- start(HiC_FIRST_GR.anno)
        annoData.ol.HiC$point_B <- end(HiC_FIRST_GR.anno)
        annoData.ol.HiC$point_C <- start(HiC_SECOND_GR.anno)
        annoData.ol.HiC$point_D <- end(HiC_SECOND_GR.anno)
        annoData.ol.HiC$point_X <- start(annoData.ol.HiC.pos)
        annoData.ol.HiC$HiC.idx <- subjectHits(ol)
        annoData.ol.HiC
    })
    if (all(elementNROWS(anno.ol.HiCdata) == 0L)) {
        return(GRanges())
    }
    HiC.idx.peaks <- unique(unlist(lapply(peaks.ol.HiCdata,
                                          function(.ele) .ele$HiC.idx),
                                   use.names = FALSE))
    HiC.idx.anno <- unique(unlist(lapply(anno.ol.HiCdata,
                                         function(.ele) .ele$HiC.idx),
                                  use.names = FALSE))
    HiC.idx <- intersect(HiC.idx.peaks, HiC.idx.anno)
    if (length(HiC.idx) == 0L) {
        return(GRanges())
    }
    anno.refined <- list()
    rotate.gr <- function(gr, anchor) {
        strand(gr) <- ifelse(strand(gr) == "-", "+", "-")
        off.pos <- end(gr) > anchor
        start(gr[off.pos]) <- anchor[off.pos] - width(gr[off.pos])
        end(gr[off.pos]) <- anchor[off.pos]
        end(gr[!off.pos]) <- anchor[!off.pos] + width(gr[!off.pos])
        start(gr[!off.pos]) <- anchor[!off.pos]
        gr
    }
    rev.gr <- function(gr, p1, p2, px) {
        tmp.shift <- p1 + p2 - 2 * px
        tmp <- shift(gr, shift = tmp.shift)
        rotate.gr(tmp, px + tmp.shift)
    }
    addListInfo <- function(l, info, infoname) {
        lapply(l, function(.ele) {
            if (length(.ele) > 0L) {
                mcols(.ele)[, infoname] <- rep(info, length(.ele))
            }
            .ele
        })
    }
    unList1level <- function(l) {
        l_offs <- unique(unlist(lapply(l, names)))
        sapply(l_offs, function(.ele) {
            .ele <- lapply(l, `[[`, .ele)
            .ele <- .ele[vapply(.ele, length, FUN.VALUE = integer(1)) > 0L]
            if (length(.ele) > 0L) {
                unlist(GRangesList(.ele), use.names = FALSE)
            } else {
                NULL
            }
        }, simplify = FALSE)
    }
    
    for (i in seq_along(anno.ol.HiCdata)) {
        this.name <- names(anno.ol.HiCdata)[i]
        this.data <- anno.ol.HiCdata[[i]]
        this.data <- this.data[this.data$HiC.idx %in% HiC.idx]
        if (length(this.data) > 0L) {
            AC.rev <- rev.gr(this.data, this.data$point_C,
                             this.data$point_A, this.data$point_X)
            AC <- shift(this.data, shift = this.data$point_C - this.data$point_A)
            CA <- shift(this.data, shift = this.data$point_A - this.data$point_C)
            AD.rev <- rev.gr(this.data, this.data$point_D,
                             this.data$point_A, this.data$point_X)
            AD <- shift(this.data, shift = this.data$point_D - this.data$point_A)
            DA <- shift(this.data, shift = this.data$point_A - this.data$point_D)
            BC.rev <- rev.gr(this.data, this.data$point_C,
                             this.data$point_B, this.data$point_X)
            BC <- shift(this.data, shift = this.data$point_C - this.data$point_B)
            CB <- shift(this.data, shift = this.data$point_B - this.data$point_C)
            BD.rev <- rev.gr(this.data, this.data$point_D,
                             this.data$point_B, this.data$point_X)
            BD <- shift(this.data, shift = this.data$point_D - this.data$point_B)
            DB <- shift(this.data, shift = this.data$point_B - this.data$point_D)
            shiftAnn <-
                switch(this.name,
                       A.ups = list(AC = list(A.ups = NULL,
                                              AB = AC.rev,
                                              BC = AC.rev,
                                              CD = AC,
                                              D.dws = AC),
                                    AD = list(A.ups = NULL,
                                              AB = AD.rev,
                                              BC = AD.rev,
                                              CD = AD.rev,
                                              D.dws = AD),
                                    BC = list(A.ups = NULL,
                                              AB = NULL,
                                              BC = BC.rev,
                                              CD = BC,
                                              D.dws = BC),
                                    BD = list(A.ups = NULL,
                                              AB = NULL,
                                              BC = BD.rev,
                                              CD = BD.rev,
                                              D.dws = BD)),
                       AB = list(AC = list(A.ups = AC.rev,
                                           AB = c(AC, CA),
                                           BC = AC,
                                           CD = AC.rev,
                                           D.dws = AC.rev),
                                 AD = list(A.ups = AD.rev,
                                           AB = c(AD, DA),
                                           BC = AD,
                                           CD = AD,
                                           D.dws = AD.rev),
                                 BC = list(A.ups = NULL,
                                           AB = NULL,
                                           BC = BC.rev,
                                           CD = BC,
                                           D.dws = BC),
                                 BD = list(A.ups = NULL,
                                           AB = NULL,
                                           BC = BD.rev,
                                           CD = BD.rev,
                                           D.dws = BD)),
                       BC = list(AC = list(A.ups = AC.rev,
                                           AB = CA,
                                           BC = c(AC, CA),
                                           CD = AC.rev,
                                           D.dws = AC.rev),
                                 AD = list(A.ups = AD.rev,
                                           AB = DA,
                                           BC = c(AD, DA),
                                           CD = AD,
                                           D.dws = AD.rev),
                                 BC = list(A.ups = BC.rev,
                                           AB = BC.rev,
                                           BC = c(BC, CB),
                                           CD = BC.rev,
                                           D.dws = BC.rev),
                                 BD = list(A.ups = BD.rev,
                                           AB = BD.rev,
                                           BC = c(BD, DB),
                                           CD = BD,
                                           D.dws = BD.rev)),
                       CD = list(AC = list(A.ups = CA,
                                           AB = AC.rev,
                                           BC = AC.rev,
                                           CD = NULL,
                                           D.dws = NULL),
                                 AD = list(A.ups = AD.rev,
                                           AB = DA,
                                           BC = DA,
                                           CD = c(AD, DA),
                                           D.dws = AD.rev),
                                 BC = list(A.ups = CB,
                                           AB = CB,
                                           BC = BC.rev,
                                           CD = NULL,
                                           D.dws = NULL),
                                 BD = list(A.ups = BD.rev,
                                           AB = BD.rev,
                                           BC = DB,
                                           CD = c(BD, DB),
                                           D.dws = BD.rev)),
                       D.dws = list(AC = list(A.ups = CA,
                                              AB = AC.rev,
                                              BC = AC.rev,
                                              CD = NULL,
                                              D.dws = NULL),
                                    AD = list(A.ups = DA,
                                              AB = AD.rev,
                                              BC = AD.rev,
                                              CD = AD.rev,
                                              D.dws = NULL),
                                    BC = list(A.ups = CB,
                                              AB = CB,
                                              BC = BC.rev,
                                              CD = NULL,
                                              D.dws = NULL),
                                    BD = list(A.ups = DB,
                                              AB = DB,
                                              BC = BD.rev,
                                              CD = BD.rev,
                                              D.dws = NULL)))
            shiftAnn <-
                mapply(addListInfo,
                       shiftAnn, names(shiftAnn),
                       infoname = "cross.link.region",
                       SIMPLIFY = FALSE)
            anno.refined[[this.name]] <- unList1level(shiftAnn)
        }
    }
    anno.refined <-
        mapply(addListInfo,
               anno.refined, names(anno.refined),
               infoname = "raw.annotation.region",
               SIMPLIFY = FALSE)
    anno.refined <- unList1level(anno.refined)
    anno.refined <- anno.refined[names(peaks.ol.HiCdata)]
    enhancer <- mapply(function(.anno, .peaks) {
        ol.HiC.idx <- intersect(.peaks$HiC.idx, .anno$HiC.idx)
        if (length(ol.HiC.idx) == 0) return(NULL)
        annot <- lapply(ol.HiC.idx, function(.HiC.id) {
            .a <- .anno[.anno$HiC.idx == .HiC.id]
            .p <- .peaks[.peaks$HiC.idx == .HiC.id]
            annoPeaks(.p, .a, bindingType = bindingType,
                      bindingRegion = bindingRegion,
                      ignore.peak.strand = ignore.peak.strand)
        })
        annot <- annot[vapply(annot, length, FUN.VALUE = integer(1)) > 0L]
        if (length(annot) == 0L) {
            return(NULL)
        }
        annot <- unlist(GRangesList(annot))
        colnames(mcols(annot)) <-
            gsub("feature\\.", "feature.shift.", colnames(mcols(annot)))
        annot
    }, anno.refined, peaks.ol.HiCdata, SIMPLIFY = FALSE)
    enhancer <- enhancer[vapply(enhancer, length, FUN.VALUE = integer(1)) > 0L]
    if (length(enhancer) == 0L) {
        return(GRanges())
    }
    enhancer <- mapply(function(.a, .n) {
        mcols(.a)[, "peak.annotation.region"] <- .n
        .a
    }, enhancer, names(enhancer), SIMPLIFY = FALSE)
    enhancer <- unlist(GRangesList(enhancer), use.names = FALSE)
    enhancer$feature.ranges <- ranges(annoData[enhancer$feature])
    enhancer$feature.strand <- strand(annoData[enhancer$feature])
    ncols <- ncol(mcols(enhancer))
    feature.col.id <- which(colnames(mcols(enhancer)) == "feature")
    mcols(enhancer) <- mcols(enhancer)[, c(1:(feature.col.id - 2),
                                           feature.col.id,
                                           ncols - 1, ncols,
                                           (feature.col.id + 1):(ncols - 2),
                                           feature.col.id - 1)]
    enhancer$point_X <- NULL
    enhancer$HiC.idx.1 <- NULL
    colnames(mcols(enhancer)) <-
        gsub("point_", "DNAinteractive_point_", colnames(mcols(enhancer)))
    colnames(mcols(enhancer)) <-
        gsub("HiC", "DNAinteractive", colnames(mcols(enhancer)))
    peak.gpid <- rle(enhancer$peak.oid.to.be.deleted)
    peak.gpid$values <- seq_along(peak.gpid$values)
    peak.gpid <- inverse.rle(peak.gpid)
    enhancer <- enhancer[order(peak.gpid, enhancer$distance)]
    enhancer <- enhancer[!duplicated(paste(enhancer$feature, enhancer$peak))]
    enhancer$peak.oid.to.be.deleted <- NULL
    enhancer$DNAinteractive.ranges <-
        ranges(DNAinteractiveData[enhancer$DNAinteractive.idx])
    enhancer$DNAinteractive.blocks <-
        DNAinteractiveData[enhancer$DNAinteractive.idx]$blocks
    enhancer$DNAinteractive_point_A <- NULL
    enhancer$DNAinteractive_point_B <- NULL
    enhancer$DNAinteractive_point_C <- NULL
    enhancer$DNAinteractive_point_D <- NULL
    enhancer$cross.link.region <- NULL
    enhancer$raw.annotation.region <- NULL
    enhancer$peak.annotation.region <- NULL
    enhancer$DNAinteractive.idx <- NULL
    enhancer
}
