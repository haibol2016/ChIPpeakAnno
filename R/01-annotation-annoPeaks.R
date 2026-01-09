#' Annotate peaks to genomic features using region expansion and overlap detection
#' 
#' This function annotates peaks to genomic features (genes, transcripts, etc.)
#' using a region-based algorithm that expands annotation regions (or peaks) by
#' a specified distance and finds overlaps. The algorithm works in three main
#' steps:
#' 
#' \enumerate{
#'   \item \strong{Region Expansion}: Based on \code{bindingType}, either annotation
#'         regions or peak regions are expanded by \code{bindingRegion} (e.g., 
#'         ±5kb from TSS). This creates expanded search regions around features.
#'   \item \strong{Overlap Detection}: Peaks that overlap with expanded annotation
#'         regions are identified using \code{findOverlaps()}. Only overlapping
#'         peaks are retained (non-overlapping peaks are excluded, not annotated
#'         with NA).
#'   \item \strong{Distance Calculation}: For overlapping peaks, distances are
#'         calculated from peaks to original (non-expanded) feature boundaries
#'         using GenomicRanges \code{distance()} function. If \code{select = "bestOne"},
#'         peaks are prioritized by shortest distance, then by highest overlap score
#'         (Jaccard index).
#' }
#' 
#' \strong{Key Characteristics:}
#' \itemize{
#'   \item \strong{Region-based}: Uses overlap detection with expanded regions,
#'         not point-to-point distance calculations
#'   \item \strong{Requires \code{bindingRegion}}: Must specify expansion distances
#'         (e.g., \code{c(-5000, 3000)} for 5kb upstream, 3kb downstream)
#'   \item \strong{Excludes non-overlapping peaks}: Peaks outside expanded regions
#'         are not returned (unlike \code{annotatePeakInBatch} which returns all
#'         peaks with NA annotations)
#'   \item \strong{Optimized for bidirectional promoters}: Special handling for
#'         finding promoters from both directions
#' }
#' 
#' Unlike internal logic of \code{\link{annotatePeakInBatch}}, which uses point-to-point 
#' distance calculations (e.g., distance from peak center to TSS), this function uses
#' region expansion and overlap detection, making it better suited for queries
#' like "find all peaks within 5kb of any TSS" or "find peaks near bidirectional
#' promoters".
#' 
#' @param peaks A \link[GenomicRanges:GRanges-class]{GRanges} object containing
#'        peaks to be annotated. If names are missing, they will be automatically
#'        generated as "X1", "X2", etc.
#' @param annoData A \link[GenomicRanges:GRanges-class]{GRanges} or
#'        \code{\link{annoGR}} object containing annotation data (genes,
#'        transcripts, etc.). The annotation object must have names for each
#'        feature.
#' @param bindingType A character specifying how to define the binding region
#'        relative to features. Default is \code{"nearestBiDirectionalPromoters"}.
#'        Available options:
#'        \itemize{
#'          \item \code{"startSite"}: Defines binding region relative to the
#'                feature start site (TSS for positive strand, gene end for
#'                negative strand). Annotation regions are expanded by
#'                \code{bindingRegion}, but constrained to not exceed original
#'                feature boundaries downstream. Use for promoter-proximal binding
#'                analysis.
#'          \item \code{"endSite"}: Defines binding region relative to the
#'                feature end site (gene end for positive strand, TSS for
#'                negative strand). Annotation regions are expanded by
#'                \code{bindingRegion}, but constrained to not exceed original
#'                feature boundaries upstream. Use for 3' end analysis.
#'          \item \code{"fullRange"}: Uses the entire feature range. Annotation
#'                regions are expanded by \code{bindingRegion} in both directions
#'                without constraints. Use for gene body or full transcript
#'                analysis.
#'          \item \code{"nearestBiDirectionalPromoters"}: Identifies peaks near
#'                bidirectional promoters using a two-step algorithm: (1) first
#'                identifies candidate bidirectional promoters by finding pairs
#'                of divergent genes (head-to-head) where their TSSs are within
#'                \code{maxTSSDistance}, then (2) finds peaks whose centers are
#'                within \code{maxPeakToBDPDistance} of the bidirectional promoter
#'                centers. This matches the literature definition (Adachi et al.
#'                2007) where bidirectional promoters are shared promoter regions
#'                between two divergently transcribed genes. Reports annotations
#'                for both genes in each bidirectional promoter pair. Note:
#'                \code{bindingRegion} parameter is ignored for this bindingType.
#'                \code{select = "bestOne"} is not supported and both genes are
#'                always kept.
#'          \item \code{"bothSidesNearest"} (deprecated): Similar to
#'                \code{"nearestBiDirectionalPromoters"} but expands peak regions
#'                instead of annotation regions. Kept for backward compatibility.
#'        }
#'        
#'        \strong{How bindingType affects region expansion:}
#'        \itemize{
#'          \item \code{"startSite"}, \code{"endSite"}, \code{"fullRange"}:
#'                Annotation regions are expanded by \code{bindingRegion}
#'          \item \code{"nearestBiDirectionalPromoters"}: Uses distance-based
#'                filtering (peak center to BDP center) instead of region expansion.
#'                \code{bindingRegion} is ignored; use \code{maxTSSDistance} and
#'                \code{maxPeakToBDPDistance} instead.
#'          \item \code{"bothSidesNearest"} (deprecated): Peak regions are
#'                expanded by \code{bindingRegion}
#'        }
#' 
#' @param bindingRegion A vector with two integer values specifying relative
#'        offsets (in base pairs) from the feature reference point defined by
#'        \code{bindingType}. Default is \code{c(-5000, 5000)}.
#'        \itemize{
#'          \item First value: Upstream offset (must be <= 0, typically negative)
#'          \item Second value: Downstream offset (must be >= 1, typically positive)
#'        }
#'        For example, \code{c(-5000, 3000)} means 5kb upstream and 3kb
#'        downstream of the feature reference point. Only peaks within this
#'        expanded region are considered for annotation.
#'        
#'        \strong{Note:} For \code{bindingType = "nearestBiDirectionalPromoters"},
#'        this parameter is ignored. Use \code{maxTSSDistance} and
#'        \code{maxPeakToBDPDistance} instead.
#' 
#' @param ignore.peak.strand A logical value indicating whether to ignore peak
#'        strand information when calculating distances. When \code{TRUE}
#'        (default), peak strand is temporarily set to "*" for distance
#'        calculations, then restored. This is appropriate for most ChIP-seq
#'        experiments where peaks are not strand-specific. Set to \code{FALSE}
#'        only if you have stranded peaks (e.g., from RNA-based methods) and
#'        want strand-aware distance calculations.
#' 
#' @param select A character specifying how to handle multiple overlapping
#'        features for a single peak. Options:
#'        \itemize{
#'          \item \code{"all"} (default): Returns all features meeting the
#'                association criteria for each peak. A peak may have multiple
#'                annotations if it overlaps with multiple features.
#'          \item \code{"bestOne"}: Returns the best feature for each peak based
#'                on: (1) shortest distance to the feature site
#'                (\code{distanceToSite}), then (2) highest overlapping score
#'                (Jaccard index of peak and feature ranges). Only one annotation
#'                per peak is returned.
#'        }
#'        Note: \code{select = "bestOne"} is automatically changed to "all"
#'        when \code{bindingType = "nearestBiDirectionalPromoters"}. For
#'        bidirectional promoters, both genes are always kept regardless of
#'        \code{select} value.
#' 
#' @param maxTSSDistance A single integer value specifying the maximum distance
#'        (in base pairs) between two divergent TSSs to be considered a
#'        bidirectional promoter. Default is \code{1000L} (1kb). This parameter
#'        is only used when \code{bindingType = "nearestBiDirectionalPromoters"}.
#'        This determines which gene pairs form bidirectional promoters.
#' 
#' @param maxPeakToBDPDistance A single integer value specifying the maximum
#'        distance (in base pairs) from peak center to bidirectional promoter
#'        center for a peak to be considered "near" a bidirectional promoter.
#'        Default is \code{200L} (200 bp). This parameter is only used when
#'        \code{bindingType = "nearestBiDirectionalPromoters"}. This
#'        distance-based filtering is more precise than overlap-based detection,
#'        focusing on peaks actually centered near the bidirectional promoter.
#' 
#' @param ... Additional parameters (currently not used)
#' 
#' @return Returns a \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing the annotated peaks. If no overlaps are found between peaks
#'        and (expanded) annotation regions, an empty \code{GRanges} object is
#'        returned.
#'        
#'        The returned object includes only peaks that overlap with (expanded)
#'        annotation regions and the following metadata columns:
#'        \itemize{
#'          \item \code{peak}: The name of the peak (from \code{names(peaks)}).
#'                If peaks had no names, auto-generated names ("X1", "X2", etc.)
#'                are used.
#'          \item \code{feature}: The name of the annotated feature (from
#'                \code{names(annoData)}). Only present if \code{annoData} has
#'                names.
#'          \item \code{feature.ranges}: The genomic ranges of the feature
#'                (as an \code{IRanges} object), representing the original
#'                feature boundaries before expansion.
#'          \item \code{feature.strand}: The strand of the feature (as a
#'                \code{Rle} object: "+", "-", or "*").
#'          \item \code{distance}: Distance from peak to feature using
#'                GenomicRanges \code{distance()} function. This is the distance
#'                between peak and feature ranges (strand-aware, but ignores peak
#'                strand if \code{ignore.peak.strand = TRUE}). Distance is 0 for
#'                overlapping ranges, positive for non-overlapping ranges.
#'          \item \code{insideFeature}: Relationship between peak and feature,
#'                determined by \code{getRelationship()} function. Possible values:
#'                \itemize{
#'                  \item \code{"upstream"}: Peak is upstream of the feature
#'                  \item \code{"downstream"}: Peak is downstream of the feature
#'                  \item \code{"inside"}: Peak is completely inside the feature
#'                  \item \code{"overlapStart"}: Peak overlaps with the start of
#'                        the feature
#'                  \item \code{"overlapEnd"}: Peak overlaps with the end of the
#'                        feature
#'                  \item \code{"includeFeature"}: Peak completely includes the
#'                        feature
#'                  \item \code{"overlap"}: Peak exactly overlaps with the feature
#'                }
#'          \item \code{distanceToSite}: Distance from peak to the feature
#'                reference site (TSS for \code{bindingType = "startSite"}, gene
#'                end for \code{bindingType = "endSite"}). This uses the original
#'                annotation ranges (before expansion) and respects
#'                \code{ignore.peak.strand} setting. This is the distance used
#'                for selecting "bestOne" when \code{select = "bestOne"}.
#'        }
#'        Additionally, all metadata columns from \code{annoData} are included
#'        (e.g., \code{gene_id}, \code{gene_name}, \code{tx_id}, \code{tx_name},
#'        etc., depending on the annotation source and feature type).
#'        
#'        \strong{Note:} For \code{bindingType = "nearestBiDirectionalPromoters"},
#'        a peak may be associated with up to two features (one from each
#'        direction) if bidirectional promoters are detected.
#'        
#'        \strong{Deprecated options:} \code{bindingType = "bothSidesNearest"}
#'        and \code{"bothSidesNSS"} are deprecated but kept for backward
#'        compatibility. \code{"bothSidesNSS"} is automatically converted to
#'        \code{"nearestBiDirectionalPromoters"}.
#' 
#' @details
#' 
#' \strong{Algorithm Steps (Detailed):}
#' 
#' The annotation process follows these steps:
#' 
#' \enumerate{
#'   \item \strong{Prepare annotation regions}: Based on \code{bindingType}, 
#'         annotation features are converted to reference points or kept as ranges:
#'         \itemize{
#'           \item \code{"startSite"}: Features are reduced to TSS points (strand-aware)
#'           \item \code{"endSite"}: Features are reduced to gene end points (strand-aware)
#'           \item \code{"fullRange"}: Features are kept as full ranges
#'           \item \code{"nearestBiDirectionalPromoters"}: Features are kept as ranges
#'         }
#'   
#'   \item \strong{Expand regions}: Based on \code{bindingType}, regions are expanded:
#'         \itemize{
#'           \item \code{"startSite"}, \code{"endSite"}, \code{"fullRange"}: 
#'                 Annotation regions are expanded by \code{bindingRegion} (e.g., 
#'                 ±5kb from TSS). For \code{"startSite"} and \code{"endSite"}, 
#'                 expansion is constrained to not exceed original feature boundaries.
#'           \item \code{"nearestBiDirectionalPromoters"}: Uses \code{promoters()} 
#'                 function to create promoter regions (upstream/downstream from TSS).
#'           \item \code{"bothSidesNearest"} (deprecated): Peak regions are expanded 
#'                 instead of annotation regions.
#'         }
#'   
#'   \item \strong{Find overlaps}: Uses \code{findOverlaps()} to identify peaks that 
#'         overlap with expanded annotation regions. Only overlapping peaks are 
#'         retained. If no overlaps are found, an empty GRanges is returned.
#'   
#'   \item \strong{Filter for bidirectional promoters} (if applicable): For 
#'         \code{bindingType = "nearestBiDirectionalPromoters"}, filters to keep 
#'         peaks associated with promoters from both directions (both strands) 
#'         or the nearest promoter in one direction.
#'   
#'   \item \strong{Calculate distances}: For each overlapping peak-feature pair:
#'         \itemize{
#'           \item \code{distance}: Distance between peak and original (non-expanded) 
#'                 feature range using \code{distance()} function (0 for overlapping, 
#'                 positive for non-overlapping)
#'           \item \code{distanceToSite}: Distance from peak to feature reference 
#'                 point (TSS for \code{"startSite"}, gene end for \code{"endSite"})
#'         }
#'   
#'   \item \strong{Select best match} (if \code{select = "bestOne"}): For each peak, 
#'         selects the best feature based on:
#'         \enumerate{
#'           \item Shortest \code{distanceToSite} (ascending order)
#'           \item Highest overlap score (Jaccard index, descending order)
#'         }
#'         The overlap score (Jaccard index) is calculated as: 
#'         \code{width(intersection) / width(union)} of peak and feature ranges.
#' }
#' 
#' \strong{Key Differences from \code{annotatePeakInBatch}:}
#' \itemize{
#'   \item \strong{Algorithm}: Uses region expansion + overlap detection vs. 
#'         point-to-point distance calculations
#'   \item \strong{Input requirement}: Always requires \code{bindingRegion} to be 
#'         specified (cannot be NULL)
#'   \item \strong{Output behavior}: Returns only overlapping peaks (excludes 
#'         non-overlapping peaks), whereas \code{annotatePeakInBatch} returns all 
#'         peaks with NA annotations for non-matching peaks
#'   \item \strong{Prioritization}: When \code{select = "bestOne"}, uses distance 
#'         + Jaccard index (overlap score) for prioritization
#'   \item \strong{Use case}: Optimized for queries like "find peaks within X kb 
#'         of TSS" or "find peaks near bidirectional promoters"
#' }
#' 
#' \strong{Region Expansion Examples:}
#' \itemize{
#'   \item \code{bindingType = "startSite"}, \code{bindingRegion = c(-2000, 500)}:
#'         Creates a region from 2kb upstream to 500bp downstream of TSS. If a gene 
#'         is 1kb long, the downstream expansion is constrained to the gene end.
#'   \item \code{bindingType = "fullRange"}, \code{bindingRegion = c(-5000, 5000)}:
#'         Expands the entire gene range by 5kb in both directions without constraints.
#'   \item \code{bindingType = "nearestBiDirectionalPromoters"}, 
#'         \code{bindingRegion = c(-5000, 3000)}: Creates promoter regions from 
#'         5kb upstream to 3kb downstream of TSS using \code{promoters()} function.
#' }
#' 
#' @importFrom GenomeInfoDb seqlevelsStyle seqlengths
#' @importFrom BiocGenerics strand start end width pos
#' @importFrom S4Vectors queryHits subjectHits
#' @author Jianhong Ou, Haibo Liu
#' @seealso \code{\link{annotatePeakInBatch}} for more flexible annotation with
#'          multiple output modes and point-based distance calculations,
#'          \code{\link{annoGR}} for creating annotation objects from TxDb or
#'          EnsDb packages
#' 
#' @references
#' Adachi, Noritaka et al. (2007). Bidirectional Gene Organization. 
#' \emph{Cell}, \bold{109(7)}: 807-809.
#' 
#' @keywords misc
#' 
#' @examples
#' \dontrun{
#' ## Example 1: Basic usage with EnsDb
#' library(EnsDb.Hsapiens.v75)
#' library(ChIPpeakAnno)
#' data("myPeakList")
#' annoGR <- annoGR(EnsDb.Hsapiens.v75, feature="gene")
#' seqlevelsStyle(myPeakList) <- seqlevelsStyle(annoGR)[1]
#' 
#' # Annotate peaks to genes within 5kb upstream and 3kb downstream of TSS
#' annotated_peaks <- annoPeaks(myPeakList, annoGR,
#'                              bindingType = "startSite",
#'                              bindingRegion = c(-5000, 3000))
#' head(annotated_peaks)
#' 
#' ## Example 2: Bidirectional promoter detection
#' annotated_peaks <- annoPeaks(myPeakList, annoGR,
#'                              bindingType = "nearestBiDirectionalPromoters",
#'                              bindingRegion = c(-5000, 3000))
#' 
#' ## Example 3: Full gene range annotation
#' annotated_peaks <- annoPeaks(myPeakList, annoGR,
#'                              bindingType = "fullRange",
#'                              bindingRegion = c(-5000, 5000))
#' 
#' ## Example 4: Select best match for each peak
#' annotated_peaks <- annoPeaks(myPeakList, annoGR,
#'                              bindingType = "startSite",
#'                              bindingRegion = c(-2000, 500),
#'                              select = "bestOne")
#' 
#' ## Example 5: 3' end analysis
#' annotated_peaks <- annoPeaks(myPeakList, annoGR,
#'                              bindingType = "endSite",
#'                              bindingRegion = c(-5000, 3000))
#' }
#' 

annoPeaks <- function(peaks,
                      annoData,
                      bindingType = c("nearestBiDirectionalPromoters",
                                     "startSite", "endSite", "fullRange"),
                      bindingRegion = c(-5000, 5000),
                      ignore.peak.strand = TRUE,
                      select = c("all", "bestOne"),
                      maxTSSDistance = 1000L,
                      ...) {
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    stopifnot(length(bindingRegion) == 2L)
    stopifnot(bindingRegion[1] <= 0L && bindingRegion[2] >= 1L)
    stopifnot(is.numeric(maxTSSDistance))
    select <- match.arg(select)
    if (is.null(names(annoData))) {
        stop("annoData must have names")
    }
    
    maxTSSDistance <- round(maxTSSDistance[1L])
    maxPeakToBDPDistance <- round(maxPeakToBDPDistance[1L])
    if (bindingType[1] %in%
        # bothSidesNearest and bothSidesNSS are deprecated, but kept for
        # backward compatibility
        c("bothSidesNearest", "nearestBiDirectionalPromoters", "bothSidesNSS")) {
        bindingType <- bindingType[1]
        if (bindingType == "bothSidesNSS") {
            bindingType <- "nearestBiDirectionalPromoters"
        }
        if (select != "all") {
            select <- "all"
            message("nearestBiDirectionalPromoters does not support ",
                    "select = 'bestOne'; using select = 'all'")
        }
    } else {
        bindingType <- match.arg(bindingType)
    }
    
    # check the seqlevelStyle of peaks and annoData
    check_seqlevel <- tryCatch({
        seqlevelsStyle(peaks)
        seqlevelsStyle(annoData)
    }, error = function(w) {
        warning("seqlevel style not recognized, you are probably ",
                "using custom GRanges objects: ", w$message)
        return(NA)
    })
    if (!any(is.na(check_seqlevel))) {
        stopifnot(length(intersect(seqlevelsStyle(peaks),
                                   seqlevelsStyle(annoData))) > 0L)
    }

    if (ignore.peak.strand) {
        peaks$peakstrand <- strand(peaks)
        strand(peaks) <- "*"
    }
    if (is.null(names(peaks))) {
        names(peaks) <- paste0("X", seq_along(peaks))
    }


    tmp <- annoData
    annotation <- switch(
        bindingType,
        startSite = {
            idx <- as.character(strand(tmp)) == "-"
            start(tmp)[idx] <- end(tmp)[idx]
            width(tmp) <- 1L
            tmp
        },
        endSite = {
            idx <- as.character(strand(tmp)) != "-"
            start(tmp)[idx] <- end(tmp)[idx]
            width(tmp) <- 1L
            tmp
        },
        fullRange = annoData,
        bothSidesNearest = annoData,
        nearestBiDirectionalPromoters = annoData,
        annoData
    )
    annotation.bck <- annotation
    rm(tmp)
    gc()

    # If bindingType is bothSidesNearest or nearestBiDirectionalPromoters,
    # find overlaps between peaks and annotation
    if (bindingType %in% c("bothSidesNearest", "nearestBiDirectionalPromoters")) {
        if (bindingType == "bothSidesNearest") { # deprecated
            # Expand peaks by bindingRegion upstream and downstream
            peaks.tmp <- .expandGRangesByBindingRegion(peaks, bindingRegion)
            ol <- findOverlaps(
                query = peaks.tmp,
                subject = annotation,
                type = "any",
                select = "all",
                ignore.strand = FALSE
            )
        } else {
            ## bindingType == "nearestBiDirectionalPromoters"
            # New algorithm: First identify bidirectional promoters, then find peaks
            # whose centers are within the bidirectional promoter regions
            all_annotated_peaks <- annotatePeaksNearBDP(peaks = peaks, annoData = annoData,
                                                        maxTSSDistance = maxTSSDistance,
                                                        ignore.peak.strand = ignore.peak.strand)
            return(all_annotated_peaks)
        }
    } else {
        # Expand annotation by bindingRegion upstream and downstream
        idx <- as.character(strand(annotation)) != "-"
        annotation <- .expandGRangesByBindingRegion(annotation, bindingRegion)
        
        # Constrain expanded annotation to not exceed original feature boundaries
        if (bindingType == "startSite") {
            ## Make sure the downstream is inside gene
            start(annotation)[!idx] <- pmax(
                start(annotation)[!idx],
                start(annoData)[!idx]
            )
            end(annotation)[idx] <- pmin(
                end(annotation)[idx],
                end(annoData)[idx]
            )
        } else if (bindingType == "endSite") {
            start(annotation)[idx] <- pmax(
                start(annotation)[idx],
                start(annoData)[idx]
            )
            end(annotation)[!idx] <- pmin(
                end(annotation)[!idx],
                end(annoData)[!idx]
            )
        }
        # Find overlaps between peaks and annotation
        ol <- findOverlaps(
            query = peaks,
            subject = annotation,
            type = "any",
            select = "all",
            ignore.strand = FALSE
        )
    }
    if (length(ol) < 1L) {
        return(GRanges())
    }
    peaks <- peaks[queryHits(ol)]
    anno <- annoData[subjectHits(ol)]
    annotation.bck.hits <- annotation.bck[subjectHits(ol)]

    # Filter overlaps results and save the nearest
    if (bindingType %in% c("bothSidesNearest", "nearestBiDirectionalPromoters")) {
        if (bindingType == "bothSidesNearest") {
            relations <- getRelationship(peaks, anno)
        } else {
            relations <- getRelationship(
                peaks,
                promoters(
                    unname(as(anno, "GRanges")),
                    upstream = 0L,
                    downstream = 1L
                )
            )
        }
        ## Filter results and save the nearest
        keep <- rep(FALSE, length(peaks))
        anno.strand <- as.character(strand(anno)) != "-"
        if (bindingType == "nearestBiDirectionalPromoters") {
            # Keep the peaks that are inside the feature or overlap the feature
            keep[relations$insideFeature %in%
                 c("includeFeature", "overlap")] <- TRUE

            ## 1. Keep the peaks that are upstream of the feature or overlap the
            ##    start of the feature on the left side (on the negative strand) of
            ##    the feature
            ## 2. Peak is completely inside the feature
            ## 3. Peak overlaps the end of the feature on the right side (on the
            ##    positive strand)
            keep.left <- (relations$insideFeature %in%
                          c("upstream", "overlapStart") & !anno.strand) |
                         (relations$insideFeature %in% "inside") |
                         ((relations$insideFeature %in% "overlapEnd") &
                          anno.strand)

            # Keep the peaks that are upstream of the feature or overlap the start
            ## 1. On the right side (on the positive strand) of the feature
            ## 2. Peak is completely inside the feature
            ## 3. Peak overlaps the start of the feature on the left side (on the
            ##    negative strand)
            keep.right <- (relations$insideFeature %in%
                           c("upstream", "overlapStart") & anno.strand) |
                          (relations$insideFeature %in% "inside") |
                          ((relations$insideFeature %in% "overlapEnd") &
                           !anno.strand)
            shortestDist <- relations$distanceToStart
        } else {
            keep[relations$insideFeature %in%
                 c("includeFeature", "inside", "overlap",
                   "overlapEnd", "overlapStart")] <- TRUE
            keep.left <- (relations$insideFeature == "downstream" &
                          anno.strand) |
                         (relations$insideFeature == "upstream" &
                          !anno.strand)
            keep.right <- (relations$insideFeature == "upstream" &
                           anno.strand) |
                           (relations$insideFeature == "downstream" &
                            !anno.strand)
            shortestDist <- relations$shortestDistance
        }
        names(shortestDist) <- seq_along(peaks)
        whichismin <- function(.ele) {
            as.numeric(names(.ele)[.ele == min(.ele)])
        }
        if (sum(keep.left) >= 1L) {
            nearest.left <- tapply(
                shortestDist[keep.left],
                queryHits(ol)[keep.left],
                whichismin,
                simplify = FALSE
            )
            keep[unlist(nearest.left)] <- TRUE
        }
        if (sum(keep.right) >= 1L) {
            nearest.right <- tapply(
                shortestDist[keep.right],
                queryHits(ol)[keep.right],
                whichismin,
                simplify = FALSE
            )
            keep[unlist(nearest.right)] <- TRUE
        }
        peaks <- peaks[keep]
        anno <- anno[keep]
        annotation.bck.hits <- annotation.bck.hits[keep]
    }
    peaks$peak <- names(peaks)
    if (!is.null(names(anno))) {
        peaks$feature <- names(anno)
    }
    peaks$feature.ranges <- unname(ranges(anno))
    peaks$feature.strand <- strand(anno)
    peaks$distance <- distance(peaks, anno, ignore.strand = FALSE)
    relations <- getRelationship(peaks, anno)
    peaks$insideFeature <- relations$insideFeature
    peaks$distanceToSite <- distance(
        peaks,
        annotation.bck.hits,
        ignore.strand = ignore.peak.strand
    )
    if (ignore.peak.strand) {
        strand(peaks) <- peaks$peakstrand
        peaks$peakstrand <- NULL
    }
    mcols(peaks) <- cbind(mcols(peaks), mcols(anno))
    if (select == "bestOne") {
        if (length(peaks) == 0L) {
            return(peaks)
        }
        annoscore <- -1L * annoScore(peaks, anno)
        peaks$ANNOPEAKS__peak.oid <- seq_along(peaks)
        # Order the peaks by the peak name, the distance to the site of the
        # feature in ascending order, and the overlapping score in descending
        # order, then keep the first one for each peak
        peaks <- peaks[order(peaks$peak, peaks$distanceToSite, annoscore)]
        peaks <- peaks[!duplicated(peaks$peak)]
        peaks <- peaks[order(peaks$ANNOPEAKS__peak.oid)]
        peaks$ANNOPEAKS__peak.oid <- NULL
    }
    peaks
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
        strand = "*",
        feature1_id = names(plus_anno_tss),
        feature2_id = names(minus_anno_tss),
        feature1_strand = strand(plus_anno_tss),
        feature2_strand = strand(minus_anno_tss),
        TSS_distance = distance(plus_anno_tss, minus_anno_tss, ignore.strand = FALSE),
        bdp_center = as.integer(round((start(minus_anno_tss) + start(plus_anno_tss)) / 2))
    ))

    list(bdp =bdp, 
         bdp_plus_anno = annoData[plus_idx][queryHits(ol)], 
         bdp_minus_anno = annoData[minus_idx][subjectHits(ol)])
}

## Helper function to expand GRanges by binding region
## Expands ranges by b[1] upstream and b[2] downstream, strand-aware
.expandGRangesByBindingRegion <- function(gr, bindingRegion) {
    str_pos <- as.character(strand(gr)) != "-"
    # For positive strand: extend upstream (subtract) and downstream (add)
    # For negative strand: extend upstream (add) and downstream (subtract)
    s1 <- ifelse(
        str_pos,
        start(gr) + bindingRegion[1],
        start(gr) - bindingRegion[2]
    )
    s1[s1 < 1L] <- 1L
    start(gr) <- s1
    e1 <- ifelse(
        str_pos,
        end(gr) + bindingRegion[2],
        end(gr) - bindingRegion[1]
    )
    # Bound by chromosome lengths
    seql <- seqlengths(gr)
    if (length(seql) > 0L) {
        e1_seql <- seql[as.character(seqnames(gr))]
        e1.idx <- which(e1 > e1_seql)
        if (length(e1.idx) > 0L) {
            e1[e1.idx] <- e1_seql[e1.idx]
        }
    }
    end(gr) <- e1
    gr
}


annotatePeaksNearBDP <- function(peaks, annoData, 
                                maxTSSDistance = 1000L, 
                                ignore.peak.strand = TRUE) {
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    stopifnot(is.numeric(maxTSSDistance))
    stopifnot(is.logical(ignore.peak.strand))
    maxTSSDistance <- round(maxTSSDistance[1L])
    ignore.peak.strand <- as.logical(ignore.peak.strand[1L])
    if (inherits(annoData, "annoGR")) {
        annoData <- as(annoData, "GRanges")
    }
    if (ignore.peak.strand) {
        peaks$peakstrand <- strand(peaks)
        strand(peaks) <- "*"
    }
    if (is.null(names(peaks))) {
        names(peaks) <- paste0("X", seq_along(peaks))
    }
    bdp_regions <- .identifyBidirectionalPromoters(annoData, maxTSSDistance)
            
    if (length(bdp_regions) == 0L) {
        return(peaks)
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
        return(peaks)
    }
    peaks <- peaks[queryHits(ol)]
    bdp_plus_anno <- bdp_regions$bdp_plus_anno[subjectHits(ol)]
    bdp_minus_anno <- bdp_regions$bdp_minus_anno[subjectHits(ol)]
    
    # For bidirectional promoters, create annotations for both genes
    # This bypasses the standard overlap filtering logic
    # Create annotations for both genes
    peak_gr1 <- peaks[queryHits(ol)]
    peak_gr1$peak <- names(peaks)[queryHits(ol)]
    peak_gr1$feature <- names(bdp_plus_anno)
    peak_gr1$feature.ranges <- ranges(bdp_plus_anno)
    peak_gr1$feature.strand <- strand(bdp_plus_anno)
    peak_gr1$distance <- distance(peak_gr1, bdp_plus_anno, 
                                    ignore.strand = FALSE)
    peak_gr1$distanceToSite <- distance(peak_gr1, bdp_regions$bdp, 
                                        ignore.strand = ignore.peak.strand)
    relations1 <- getRelationship(peak_gr1, bdp_plus_anno)
    peak_gr1$insideFeature <- relations1$insideFeature
    mcols(peak_gr1) <- cbind(mcols(peak_gr1), mcols(bdp_plus_anno))
    
    peak_gr2 <- peaks[queryHits(ol)]
    peak_gr2$peak <- names(peaks)[queryHits(ol)]
    peak_gr2$feature <- names(bdp_minus_anno)
    peak_gr2$feature.ranges <- ranges(bdp_minus_anno)
    peak_gr2$feature.strand <- strand(bdp_minus_anno)
    peak_gr2$distance <- distance(peak_gr2, bdp_minus_anno,
                                    ignore.strand = FALSE)
    peak_gr2$distanceToSite <- distance(peak_gr2, bdp_regions$bdp, 
                                        ignore.strand = ignore.peak.strand)
    relations2 <- getRelationship(peak_gr2, bdp_minus_anno)
    peak_gr2$insideFeature <- relations2$insideFeature
    mcols(peak_gr2) <- cbind(mcols(peak_gr2), mcols(bdp_minus_anno))
    all_annotated_peaks <- c(peak_gr1, peak_gr2)
    all_annotated_peaks <- all_annotated_peaks[order(seqnames(all_annotated_peaks),
                                                    start(all_annotated_peaks),
                                                    end(all_annotated_peaks))]

    all_annotated_peaks
}