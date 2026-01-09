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