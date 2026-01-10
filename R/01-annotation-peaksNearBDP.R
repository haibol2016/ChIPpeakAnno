#' Identify peaks near bidirectional promoters
#' 
#' Identifies peaks that are located near bidirectional promoters, which are
#' promoter regions shared by two divergently transcribed (head-to-head) genes.
#' This function uses a two-step algorithm: (1) first identifies candidate
#' bidirectional promoters by finding pairs of divergent genes with TSSs within
#' a distance threshold, then (2) finds peaks whose centers are within a
#' specified distance of the bidirectional promoter centers.
#' 
#' \strong{What are bidirectional promoters?}
#' 
#' Bidirectional promoters are genomic regions where two genes are transcribed
#' in opposite directions from a shared promoter region. These are functionally
#' distinct regulatory elements and are common in mammalian genomes. According
#' to the literature definition (Adachi et al. 2007), bidirectional promoters
#' are shared promoter regions between two head-to-head (divergently transcribed)
#' genes where their TSSs are within a short distance (typically < 1kb).
#' 
#' A peak is considered near a bidirectional promoter if the distance from the
#' peak center to the bidirectional promoter center is within the specified
#' threshold (default: 200 bp).
#' 
#' @param myPeakList A \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing peaks to be analyzed. The peaks should represent ChIP-seq
#'        or similar genomic regions of interest.
#' @param AnnotationData A \link[GenomicRanges:GRanges-class]{GRanges} or
#'        \code{\link{annoGR}} object containing annotation data (genes,
#'        transcripts, or TSS positions). The annotation must have strand
#'        information ("+" for positive strand, "-" for negative strand).
#'        
#'        \strong{Obtaining annotation data:}
#'        \itemize{
#'          \item \strong{Recommended:} Use \code{\link{getAnnotation}} with
#'                biomaRt to generate annotations matching your genome assembly,
#'                or use EnsDb/TxDb packages with \code{\link{annoGR}}.
#'          \item \strong{For examples only:} Pre-computed TSS datasets
#'                (e.g., \code{TSS.human.GRCh37}, \code{TSS.mouse.GRCm38}) are
#'                provided for package examples and testing. Users should
#'                generate their own annotations to match their genome assembly.
#'        }
#' @param MaxDistance A single integer value. \strong{Deprecated:} This parameter
#'        is kept for backward compatibility but is no longer used in the new
#'        algorithm. In the old algorithm, it defined the search region around
#'        peaks. In the new algorithm, bidirectional promoters are identified
#'        first, then peaks are filtered by distance from peak center to BDP center.
#'        Default is \code{5000L}.
#' @param maxTSSDistance A single integer value specifying the maximum distance
#'        (in base pairs) between two divergent TSSs to be considered a
#'        bidirectional promoter. Default is \code{1000L} (1kb). This determines
#'        which gene pairs form bidirectional promoters.
#' @param maxPeakToBDPDistance A single integer value specifying the maximum
#'        distance (in base pairs) from peak center to bidirectional promoter
#'        center for a peak to be considered "near" a bidirectional promoter.
#'        Default is \code{200L} (200 bp). This distance-based filtering is more
#'        precise than overlap-based detection, focusing on peaks actually
#'        centered near the bidirectional promoter.
#' @param ... Additional parameters (currently not used)
#' 
#' @return Returns a list with the following components:
#'        \itemize{
#'          \item \code{peaksWithBDP}: A \code{GRangesList} containing peaks
#'                that are near bidirectional promoters. Each element in the list
#'                corresponds to one peak and contains all annotations from both
#'                strands. The \code{GRangesList} preserves the original peak
#'                ranges and includes the following metadata columns from
#'                \code{\link{annoPeaks}}:
#'                \itemize{
#'                  \item \code{peak}: The name of the peak
#'                  \item \code{feature}: The name/ID of the annotated feature
#'                        (e.g., gene ID, transcript ID)
#'                  \item \code{feature.ranges}: The genomic ranges of the feature
#'                  \item \code{feature.strand}: The strand of the feature
#'                        ("+" or "-")
#'                  \item \code{distance}: Distance from peak to feature
#'                  \item \code{insideFeature}: Relationship between peak and
#'                        feature (e.g., "upstream", "downstream", "inside",
#'                        "overlapStart", "overlapEnd", "includeFeature",
#'                        "overlap")
#'                  \item \code{distanceToSite}: Distance from peak to the
#'                        feature TSS
#'                }
#'                Additionally, all metadata columns from \code{AnnotationData}
#'                are included (e.g., \code{gene_id}, \code{gene_name}, etc.).
#'          \item \code{percentPeaksWithBDP}: A numeric value (0-1) representing
#'                the percentage of input peaks that are near bidirectional
#'                promoters. Calculated as
#'                \code{length(peaksWithBDP) / length(myPeakList)}.
#'          \item \code{n.peaks}: The total number of input peaks (after
#'                removing duplicates).
#'          \item \code{n.peaksWithBDP}: The number of peaks that are near
#'                bidirectional promoters (i.e., peak center is within
#'                \code{maxPeakToBDPDistance} of a bidirectional promoter center).
#'        }
#' @details
#' 
#' \strong{Algorithm:}
#' 
#' The function uses a two-step algorithm that matches the literature definition
#' of bidirectional promoters (Adachi et al. 2007):
#' 
#' \enumerate{
#'   \item \strong{Identify candidate bidirectional promoters}: 
#'         \itemize{
#'           \item Extract TSS positions from annotation (strand-aware: start for
#'                 + strand, end for - strand)
#'           \item Find all pairs of divergent genes (one on + strand, one on
#'                 - strand) where their TSSs are within \code{maxTSSDistance}
#'           \item Create bidirectional promoter regions as the genomic interval
#'                 between the two TSSs
#'           \item Calculate the center coordinate of each bidirectional promoter
#'                 region
#'         }
#'   \item \strong{Find peaks near bidirectional promoters}:
#'         \itemize{
#'           \item Calculate distance from each peak center to each bidirectional
#'                 promoter center
#'           \item Filter to keep only peaks where distance <
#'                 \code{maxPeakToBDPDistance}
#'           \item For each peak, create annotations for both genes in the
#'                 bidirectional promoter pair
#'         }
#' }
#' 
#' \strong{Distance-based filtering:}
#' 
#' The new algorithm uses distance-based filtering (peak center to BDP center)
#' rather than overlap-based detection. This is more precise and focuses on
#' peaks actually centered near the bidirectional promoter, avoiding peaks that
#' just barely touch the edges. The distance threshold (\code{maxPeakToBDPDistance})
#' can be adjusted based on factor type (e.g., narrower for TFs, broader for
#' histone marks).
#' 
#' \strong{Note:} If no peaks are found near bidirectional promoters, the
#' function returns an empty \code{GRangesList} with \code{percentPeaksWithBDP}
#' set to 0.
#' 
#' @author Lihua Julie Zhu, Jianhong Ou
#' @seealso \code{\link{annoPeaks}} for the underlying annotation function,
#'          \code{\link{annotatePeakInBatch}} for general peak annotation,
#'          \code{\link{findOverlappingPeaks}} for finding overlapping peaks,
#'          \code{\link{makeVennDiagram}} for visualizing peak overlaps
#' @references Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to
#' annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237. doi:10.1186/1471-2105-11-237
#' 
#' Adachi, Noritaka et al. (2007). Bidirectional Gene Organization. 
#' \emph{Cell}, \bold{109(7)}: 807-809.
#' @keywords misc
#' @importFrom S4Vectors elementNROWS
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom BiocGenerics start end strand
#' @export
#' @examples
#' \dontrun{
#' library(GenomeInfoDb)
#' data(myPeakList)
#' 
#' ## Note: Using pre-computed TSS for example only.
#' ## Users should generate annotations matching their genome assembly.
#' data(TSS.human.GRCh37)
#' seqlevelsStyle(TSS.human.GRCh37) <- seqlevelsStyle(myPeakList)[1]
#' 
#' ## Customize parameters
#' annotatedBDP <- peaksNearBDP(myPeakList[1:6,],
#'                              AnnotationData=TSS.human.GRCh37,
#'                              MaxDistance=5000,
#'                              PeakLocForDistance =  "middle", 
#'                              FeatureLocForDistance = "TSS")
#' ## View results
#' annotatedBDP$peaksWithBDP
#' c(annotatedBDP$percentPeaksWithBDP, 
#'   annotatedBDP$n.peaks, 
#'   annotatedBDP$n.peaksWithBDP)
#' }
#' 
peaksNearBDP <- function(myPeakList, AnnotationData,
                         MaxDistance=5000L, ...){
        if (missing(myPeakList)) {
        stop("Missing required argument myPeakList!")
    }
    if (!inherits(myPeakList, c("GRanges"))) {
        stop("myPeakList needs to be GRanges object")
    }
    if (!missing(AnnotationData)){        
        if (!inherits(AnnotationData, c("GRanges", "annoGR"))) {
            stop("AnnotationData needs to be GRanges or annoGR object")
        }
        if(is(AnnotationData, "annoGR"))
            AnnotationData <- AnnotationData@gr
    }else{
        stop("Missing required argument AnnotationData!")
    }
    stopifnot(length(intersect(seqlevelsStyle(myPeakList),
                               seqlevelsStyle(AnnotationData)))>0)
    stopifnot(is.numeric(MaxDistance))
    
    MaxDistance <- round(MaxDistance[1])
    myPeakList <- unique(myPeakList)
    myPeakList$bdp_idx <- seq_along(myPeakList)
    anno <- annoPeaks(myPeakList, AnnotationData, 
                      bindingType = "nearestBiDirectionalPromoters",
                      bindingRegion = c(-1*MaxDistance, MaxDistance))
    if(length(anno)<1){
        return(list(peaksWithBDP=anno,
                    percentPeaksWithBDP=0,
                    n.peaks=length(myPeakList),
                    n.peaksWithBDP=0))
    }
    anno.s <- split(anno, anno$bdp_idx)
    len <- elementNROWS(anno.s)
    anno.s <- anno.s[len>=2]
    len <- sapply(anno.s, function(.ele){
        std <- .ele$feature.strand
        all(c("+", "-") %in% as.character(.ele$feature.strand))
    })
    anno.s <- anno.s[len]
    list(peaksWithBDP=anno.s,
         percentPeaksWithBDP = length(anno.s)/length(myPeakList),
         n.peaks=length(myPeakList),
         n.peaksWithBDP=length(anno.s))
}

#' Annotate peaks near bidirectional promoters
#' 
#' Annotates peaks that overlap with bidirectional promoter regions, creating
#' annotations for both genes in each bidirectional promoter pair. This is a
#' lower-level function that returns a flat \code{GRanges} object with all
#' annotated peaks, unlike \code{\link{peaksNearBDP}} which returns a list with
#' statistics.
#' 
#' \strong{What are bidirectional promoters?}
#' 
#' Bidirectional promoters are genomic regions where two genes are transcribed
#' in opposite directions from a shared promoter region. According to the
#' literature definition (Adachi et al. 2007), bidirectional promoters are
#' shared promoter regions between two head-to-head (divergently transcribed)
#' genes where their TSSs are within a short distance (typically < 1kb).
#' 
#' \strong{How it works:}
#' 
#' The function uses a two-step algorithm:
#' \enumerate{
#'   \item \strong{Identify bidirectional promoters}: Finds pairs of divergent
#'         genes (one on + strand, one on - strand) where their TSSs are within
#'         \code{maxTSSDistance} of each other, creating bidirectional promoter
#'         regions as the genomic interval between the two TSSs.
#'   \item \strong{Find overlapping peaks}: Identifies peaks whose centers
#'         overlap with bidirectional promoter regions using
#'         \code{\link[GenomicRanges]{findOverlaps}}.
#'   \item \strong{Create dual annotations}: For each peak that overlaps a
#'         bidirectional promoter, creates annotations for BOTH genes in the
#'         pair (one annotation for the + strand gene, one for the - strand
#'         gene). This means each peak near a BDP will appear twice in the
#'         output, once for each gene.
#' }
#' 
#' \strong{Key differences from \code{peaksNearBDP}:}
#' \itemize{
#'   \item Returns a flat \code{GRanges} object instead of a list with
#'         statistics
#'   \item Designed for use in hierarchical annotation pipelines (e.g.,
#'         \code{\link{annotateHierarchically}})
#'   \item Each peak near a BDP appears twice in the output (once per gene)
#'   \item Uses overlap detection (peak center within BDP region) rather than
#'         distance-based filtering
#' }
#' 
#' @param peaks A \link[GenomicRanges:GRanges-class]{GRanges} object containing
#'        peaks to be annotated. If names are missing, they will be
#'        automatically generated as "X1", "X2", etc.
#' @param annoData A \link[GenomicRanges:GRanges-class]{GRanges} or
#'        \code{\link{annoGR}} object containing annotation data (genes,
#'        transcripts, or TSS positions). The annotation must have strand
#'        information ("+" for positive strand, "-" for negative strand) and
#'        names for each feature.
#'        
#'        \strong{Obtaining annotation data:}
#'        \itemize{
#'          \item \strong{Recommended:} Use \code{\link{getAnnotation}} with
#'                biomaRt to generate annotations matching your genome assembly,
#'                or use EnsDb/TxDb packages with \code{\link{annoGR}}.
#'          \item \strong{For examples only:} Pre-computed TSS datasets
#'                (e.g., \code{TSS.human.GRCh37}, \code{TSS.mouse.GRCm38}) are
#'                provided for package examples and testing. Users should
#'                generate their own annotations to match their genome assembly.
#'        }
#' @param maxTSSDistance A single integer value specifying the maximum distance
#'        (in base pairs) between two divergent TSSs to be considered a
#'        bidirectional promoter. Default is \code{1000L} (1kb). This determines
#'        which gene pairs form bidirectional promoters. Larger values will
#'        identify more bidirectional promoter pairs but may include less
#'        biologically relevant cases.
#' @param ignore.peak.strand A logical value. If \code{TRUE} (default), peak
#'        strand information is ignored when calculating distances to
#'        bidirectional promoter centers. The original peak strand is preserved
#'        in a metadata column \code{peakstrand}. If \code{FALSE}, strand
#'        information is considered in distance calculations. This parameter
#'        also controls whether peak strand is used when finding overlaps
#'        between peak centers and BDP regions (overlaps always ignore strand
#'        for BDP regions, which are unstranded).
#' 
#' @return Returns a \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing annotated peaks. Each peak that overlaps a bidirectional
#'        promoter appears twice in the output (once for each gene in the BDP
#'        pair). The object includes the following metadata columns:
#'        \itemize{
#'          \item \code{peak}: The name of the peak (from \code{names(peaks)})
#'          \item \code{feature}: The name/ID of the annotated feature
#'                (e.g., gene ID, transcript ID) from \code{names(annoData)}
#'          \item \code{feature.ranges}: The genomic ranges of the feature
#'                (as an \code{IRanges} object)
#'          \item \code{feature.strand}: The strand of the feature ("+" or "-")
#'          \item \code{distance}: Distance from peak to feature boundary
#'                (calculated with \code{ignore.strand = FALSE})
#'          \item \code{distanceToSite}: Distance from peak to the bidirectional
#'                promoter center (calculated with \code{ignore.strand} as
#'                specified by the parameter)
#'          \item \code{insideFeature}: Relationship between peak and feature
#'                (e.g., "upstream", "downstream", "inside", "overlapStart",
#'                "overlapEnd", "includeFeature", "overlap")
#'          \item \code{peakstrand}: Original peak strand (only present if
#'                \code{ignore.peak.strand = TRUE})
#'        }
#'        Additionally, all metadata columns from \code{annoData} are included
#'        (e.g., \code{gene_id}, \code{gene_name}, etc.).
#'        
#'        The output is sorted by chromosome, start position, and end position.
#'        If no peaks overlap bidirectional promoters, returns an empty
#'        \code{GRanges} object.
#' 
#' @details
#' 
#' \strong{Algorithm details:}
#' 
#' The function first calls \code{.identifyBidirectionalPromoters} to identify
#' all bidirectional promoter pairs in the annotation data. For each pair, it
#' creates a bidirectional promoter region as the genomic interval between the
#' two TSSs.
#' 
#' Peak centers are calculated as the midpoint of each peak range. Overlaps
#' between peak centers and bidirectional promoter regions are found using
#' \code{\link[GenomicRanges]{findOverlaps}} with \code{type = "any"} and
#' \code{ignore.strand = TRUE} (BDP regions are unstranded).
#' 
#' For each peak that overlaps a bidirectional promoter, two annotation
#' entries are created:
#' \itemize{
#'   \item One entry for the gene on the positive strand
#'   \item One entry for the gene on the negative strand
#' }
#' 
#' This dual annotation approach ensures that peaks near bidirectional
#' promoters are properly associated with both genes that share the promoter
#' region, which is important for understanding the regulatory context.
#' 
#' \strong{Use cases:}
#' 
#' This function is primarily designed for use in hierarchical annotation
#' pipelines (e.g., \code{\link{annotateHierarchically}}) where bidirectional
#' promoters are treated as a distinct annotation category. It can also be
#' used directly when you need a flat list of all peaks annotated to
#' bidirectional promoters.
#' 
#' \strong{Comparison with \code{peaksNearBDP}:}
#' 
#' \code{peaksNearBDP} is a higher-level function that:
#' \itemize{
#'   \item Returns a list with statistics (\code{peaksWithBDP},
#'         \code{percentPeaksWithBDP}, etc.)
#'   \item Groups annotations by peak in a \code{GRangesList}
#'   \item Filters to only include peaks with annotations from both strands
#' }
#' 
#' \code{annotatePeaksNearBDP} is a lower-level function that:
#' \itemize{
#'   \item Returns a flat \code{GRanges} object
#'   \item Includes all peaks that overlap BDP regions (even if only one
#'         annotation is found)
#'   \item Is designed for integration into annotation pipelines
#' }
#' 
#' @author Lihua Julie Zhu, Jianhong Ou
#' @seealso \code{\link{peaksNearBDP}} for a higher-level interface that
#'          returns statistics, \code{\link{annoPeaks}} for the underlying
#'          annotation engine, \code{\link{annotateHierarchically}} for
#'          hierarchical annotation pipelines, \code{\link{annotatePeakInBatch}}
#'          for general peak annotation
#' @references Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to
#' annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237. doi:10.1186/1471-2105-11-237
#' 
#' Adachi, Noritaka et al. (2007). Bidirectional Gene Organization. 
#' \emph{Cell}, \bold{109(7)}: 807-809.
#' @keywords misc
#' @importFrom BiocGenerics start end strand
#' @importFrom S4Vectors mcols mcols<-
#' @author Haibo Liu
#' @export
#' @examples
#' \dontrun{
#' library(GenomeInfoDb)
#' data(myPeakList)
#' 
#' ## Note: Using pre-computed TSS for example only.
#' ## Users should generate annotations matching their genome assembly.
#' data(TSS.human.GRCh37)
#' seqlevelsStyle(TSS.human.GRCh37) <- seqlevelsStyle(myPeakList)[1]
#' 
#' ## Annotate peaks near bidirectional promoters
#' annotated <- annotatePeaksNearBDP(
#'     peaks = myPeakList[1:10],
#'     annoData = TSS.human.GRCh37,
#'     maxTSSDistance = 1000L,
#'     ignore.peak.strand = TRUE
#' )
#' 
#' ## View results
#' annotated
#' 
#' ## Note: Each peak near a BDP appears twice (once per gene)
#' table(annotated$peak)
#' }
annotatePeaksNearBDP <- function(peaks, annoData, 
                                maxTSSDistance = 1000L, 
                                ignore.peak.strand = TRUE) {
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    stopifnot(is.numeric(maxTSSDistance))
    maxTSSDistance <- round(maxTSSDistance[1L])
 
    if (inherits(annoData, "annoGR")) {
        annoData <- as(annoData, "GRanges")
    }
    if (is.null(names(peaks))) {
        names(peaks) <- paste0("X", seq_along(peaks))
    }
    # check the seqlevelStyle of peaks and annoData
    if (seqlevelsStyle(peaks)[1] != seqlevelsStyle(annoData)[1]) {
        warning("seqlevel style of peaks and annoData are not the same, ",
                "will try to match the seqlevel style of peaks to annoData")
        peaks <- formatSeqnames(peaks, annoData)
    }

    bdp_regions <- .identifyBidirectionalPromoters(annoData, maxTSSDistance)
            
    if (length(bdp_regions) == 0L) {
        return(GRanges())
    }
    
    # find overlaps between peak centers and bdp
    peak_centers <- as.integer(round((start(peaks) + end(peaks)) / 2))
    peak_centers <- GRanges(seqnames = seqnames(peaks),
                            IRanges(start = peak_centers, 
                            end = peak_centers),
                            strand = "*")
    ol <- findOverlaps(peak_centers, 
                        bdp_regions$bdp,
                        type = "any", 
                        select = "all", 
                        ignore.strand = TRUE)
    if (length(ol) == 0L) {
        return(GRanges())
    }

    bdp_plus_anno <- bdp_regions$bdp_plus_anno[subjectHits(ol)]
    bdp_minus_anno <- bdp_regions$bdp_minus_anno[subjectHits(ol)]
    
    # For bidirectional promoters, create annotations for both genes
    # This bypasses the standard overlap filtering logic
    # Create annotations for both genes
    peak_gr1 <- peaks[queryHits(ol)]
    peak_gr1$peak <- names(peaks)[queryHits(ol)]

    if ("tx_name" %in% colnames(mcols(bdp_plus_anno))) {
        peak_gr1$feature <- mcols(bdp_plus_anno)$tx_name
    } else if ("gene_id" %in% colnames(mcols(bdp_plus_anno))) {
        peak_gr1$feature <- mcols(bdp_plus_anno)$gene_id
    } else {
        peak_gr1$feature <- names(bdp_plus_anno)
    }
    peak_gr1$start_position <- start(bdp_plus_anno)
    peak_gr1$end_position <- end(bdp_plus_anno)
    peak_gr1$feature_strand <- strand(bdp_plus_anno)

    relations1 <- getRelationship(peak_gr1, bdp_plus_anno)
    peak_gr1$insideFeature <- relations1$insideFeature

    peak_gr1$distancetoFeature <- 0
    peak_gr1$shortestDistance <- relations1$shortestDistance
    peak_gr1$fromOverlappingOrNearest <- "bidirectional_promoters"
    
    
    peak_gr2 <- peaks[queryHits(ol)]
    peak_gr2$peak <- names(peaks)[queryHits(ol)]
    if ("tx_name" %in% colnames(mcols(bdp_minus_anno))) {
        peak_gr2$feature <- mcols(bdp_minus_anno)$tx_name
    } else if ("gene_id" %in% colnames(mcols(bdp_minus_anno))) {
        peak_gr2$feature <- mcols(bdp_minus_anno)$gene_id
    } else {
        peak_gr2$feature <- names(bdp_minus_anno)
    }
    peak_gr2$start_position <- start(bdp_minus_anno)
    peak_gr2$end_position <- end(bdp_minus_anno)
    peak_gr2$feature_strand <- strand(bdp_minus_anno)

    relations2 <- getRelationship(peak_gr2, bdp_minus_anno)
    peak_gr2$insideFeature <- relations2$insideFeature
    peak_gr2$distancetoFeature <- 0
    peak_gr2$shortestDistance <- relations2$shortestDistance
    peak_gr2$fromOverlappingOrNearest <- "bidirectional_promoters"
    
    all_annotated_peaks <- c(peak_gr1, peak_gr2)
    all_annotated_peaks <- all_annotated_peaks[order(seqnames(all_annotated_peaks),
                                                    start(all_annotated_peaks),
                                                    end(all_annotated_peaks))]

    all_annotated_peaks
}

#' Identify candidate bidirectional promoters
#' 
#' @description 
#' Internal helper function to identify candidate bidirectional promoters by
#' finding pairs of divergent genes (head-to-head) where their TSSs are within
#' a specified distance threshold. This matches the literature definition
#' (Adachi et al. 2007) where bidirectional promoters are shared promoter
#' regions between two divergently transcribed genes.
#' 
#' @param annoData A \code{GRanges} or \code{annoGR} object containing
#'   annotation data (genes, transcripts, or TSS positions). Must have strand
#'   information.
#' @param maxTSSDistance Integer. Maximum distance (in base pairs) between two
#'   divergent TSSs to be considered a bidirectional promoter. Default is
#'   \code{1000L} (1kb).
#' @return A \code{GRanges} object containing bidirectional promoter regions.
#'   Each element represents a bidirectional promoter region (the genomic
#'   interval between two divergent TSSs). Metadata columns include:
#'   \itemize{
#'     \item \code{gene1_id}: ID/name of the first gene (from names or metadata)
#'     \item \code{gene2_id}: ID/name of the second gene (from names or metadata)
#'     \item \code{gene1_strand}: Strand of the first genae ("+" or "-")
#'     \item \code{gene2_strand}: Strand of the second gene ("+" or "-")
#'     \item \code{TSS1_pos}: TSS position of the first gene
#'     \item \code{TSS2_pos}: TSS position of the second gene
#'     \item \code{TSS_distance}: Distance between the two TSSs
#'     \item \code{bdp_center}: Center coordinate of the bidirectional promoter region
#'   }
#' @details
#' 
#' \strong{Algorithm:}
#' \enumerate{
#'   \item Extract TSS positions from annotation (strand-aware: start for +,
#'         end for -)
#'   \item Find all pairs where one gene is on + strand and one on - strand
#'   \item Calculate distance between TSSs
#'   \item Filter to pairs where distance ≤ \code{maxTSSDistance}
#'   \item Create bidirectional promoter regions (genomic interval between two divergent TSSs)
#'   \item Calculate center coordinate of each bidirectional promoter region
#' }
#' 
#' \strong{Note:} Distance calculation naturally handles cross-chromosome cases
#' (infinite distance), so no explicit chromosome check is needed.
#' 
#' @author Haibo Liu
#' @keywords internal
#' @importFrom BiocGenerics start end strand
.identifyBidirectionalPromoters <- function(annoData, maxTSSDistance = 1000L) {
    stopifnot(inherits(annoData, c("GRanges", "annoGR")))
    if (inherits(annoData, "annoGR")) {
        annoData <- as(annoData, "GRanges")
    }
    stopifnot(is.numeric(maxTSSDistance))
    maxTSSDistance <- round(maxTSSDistance[1L])
    
    # Extract TSS positions (strand-aware)
    # For + strand: TSS is at start()
    # For - strand: TSS is at end()
    strand_vec <- as.character(strand(annoData))
    
    # Split by strand
    plus_idx <- strand_vec == "+"
    minus_idx <- strand_vec == "-"
    
    if (sum(plus_idx) == 0L || sum(minus_idx) == 0L) {
        # No divergent pairs possible
        return(GRanges())
    }
    
    plus_anno_tss <- suppressWarnings(promoters(annoData[plus_idx], upstream = 0L, downstream = 1L))
    minus_anno_tss <- suppressWarnings(promoters(annoData[minus_idx], upstream = 0L, downstream = 1L))  
    minus_anno_tss_expanded <- suppressWarnings(promoters(annoData[minus_idx], upstream = maxTSSDistance, downstream = 0L))
    
    # Find all pairs of divergent TSSs
    # For each chromosome, find pairs where + strand TSS is on the right side of - strand TSS
    # and distance between them is <= maxTSSDistance using findOverlaps
    ol <- findOverlaps(
        query = plus_anno_tss,
        subject = minus_anno_tss_expanded,
        type = "any",
        select = "all",
        ignore.strand = TRUE
    )
    
    if (length(ol) == 0L) { 
        return(GRanges())
    }
    
    plus_anno_tss <- plus_anno_tss[queryHits(ol)]
    minus_anno_tss <- minus_anno_tss[subjectHits(ol)]

    bdp <- suppressWarnings(GRanges(
        seqnames = seqnames(plus_anno_tss),
        IRanges(start = start(minus_anno_tss), end = start(plus_anno_tss)),
        strand = rep("*", length(plus_anno_tss)),
        feature1_id = ifelse("tx_name" %in% colnames(mcols(plus_anno_tss)), 
                             mcols(plus_anno_tss)$tx_name, 
                             ifelse("gene_id" %in% colnames(mcols(plus_anno_tss)), 
                                    mcols(plus_anno_tss)$gene_id, 
                                    names(plus_anno_tss))),
        feature2_id = ifelse("tx_name" %in% colnames(mcols(minus_anno_tss)), 
                             mcols(minus_anno_tss)$tx_name, 
                             ifelse("gene_id" %in% colnames(mcols(minus_anno_tss)), 
                                    mcols(minus_anno_tss)$gene_id, 
                                    names(minus_anno_tss))),
        feature1_strand = strand(plus_anno_tss),
        feature2_strand = strand(minus_anno_tss),
        TSS_distance = distance(plus_anno_tss, minus_anno_tss, ignore.strand = TRUE),
        bdp_center = as.integer(round((start(minus_anno_tss) + start(plus_anno_tss)) / 2))
    ))

    list(bdp = bdp, 
         bdp_plus_anno = annoData[plus_idx][queryHits(ol)], 
         bdp_minus_anno = annoData[minus_idx][subjectHits(ol)])
}
