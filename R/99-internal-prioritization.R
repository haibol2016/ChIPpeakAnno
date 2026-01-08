#' Factor-specific annotation prioritization functions
#' 
#' @description 
#' Internal helper functions for factor-specific annotation prioritization.
#' These functions implement the prioritization strategies described in
#' Peak_Annotation_Prioritization_Strategies.md.
#' 
#' @keywords internal
#' @importFrom BiocGenerics start end strand
#' @importFrom S4Vectors mcols

# Helper function: Default value operator
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

# Helper function: Calculate Jaccard index
.calculateJaccardIndex <- function(peak_ranges, feature_ranges) {
    # Ensure both are GRanges objects
    if (!inherits(peak_ranges, "GRanges") || !inherits(feature_ranges, "GRanges")) {
        stop("Both inputs must be GRanges objects", call. = FALSE)
    }
    if (length(peak_ranges) != length(feature_ranges)) {
        stop("peak_ranges and feature_ranges must have the same length", call. = FALSE)
    }
    
    # Calculate intersection and union using GenomicRanges functions
    intersection <- pintersect(peak_ranges, feature_ranges, ignore.strand = TRUE)
    union <- punion(peak_ranges, feature_ranges, 
                   fill.gap = TRUE, 
                   ignore.strand = TRUE)
    
    # Jaccard index = intersection / union
    jaccard <- width(intersection) / width(union)
    
    # Handle cases where union is 0 (no overlap)
    jaccard[is.na(jaccard)] <- 0
    
    return(jaccard)
}

# Helper function: Get Jaccard index from annotations
.getJaccardIndex <- function(annotated_peaks) {
    if ("annoScore" %in% colnames(mcols(annotated_peaks))) {
        # Use existing annoScore from annoPeaks() output (recommended)
        return(mcols(annotated_peaks)$annoScore)
    } else if ("jaccard" %in% colnames(mcols(annotated_peaks))) {
        return(mcols(annotated_peaks)$jaccard)
    } else {
        # Try to calculate from ranges if available
        if ("feature.ranges" %in% colnames(mcols(annotated_peaks))) {
            peak_ranges <- GRanges(
                seqnames = seqnames(annotated_peaks),
                ranges = IRanges(start = start(annotated_peaks), end = end(annotated_peaks))
            )
            feature_ranges <- mcols(annotated_peaks)$feature.ranges
            if (inherits(feature_ranges, "GRanges")) {
                return(.calculateJaccardIndex(peak_ranges, feature_ranges))
            }
        }
        # Cannot calculate - return 0
        return(rep(0, length(annotated_peaks)))
    }
}

# Helper function: Two-part exponential decay for distance scoring
# This provides continuity at the threshold and different decay rates
# for near vs far distances
.calculateDistanceScore <- function(distance, 
                                     max_distance,
                                     weight,
                                     near_decay_rate,
                                     far_decay_rate = NULL) {
    # Default far_decay_rate to 4x near_decay_rate if not specified
    if (is.null(far_decay_rate)) {
        far_decay_rate <- near_decay_rate * 4
    }
    
    # Two-part exponential decay:
    # - Within threshold: sharp decay (near_decay_rate)
    # - Beyond threshold: slower decay (far_decay_rate) with continuity at boundary
    ifelse(
        distance <= max_distance,
        weight * exp(-distance / near_decay_rate),
        weight * exp(-max_distance / near_decay_rate) * 
            exp(-(distance - max_distance) / far_decay_rate)
    )
}

# Helper function: Apply Tier 1 - Feature relationship priority
.applyFeaturePriority <- function(annotated_peaks, feature_priority) {
    annotated_peaks$priority_score <- 
        feature_priority[mcols(annotated_peaks)$insideFeature]
    annotated_peaks$priority_score[is.na(annotated_peaks$priority_score)] <- 0
    return(annotated_peaks)
}

# Helper function: Apply Tier 2 - Distance scoring
.applyDistanceScore <- function(annotated_peaks,
                                 max_distance,
                                 weight,
                                 near_decay_rate,
                                 far_decay_rate = NULL) {
    if ("distanceToSite" %in% colnames(mcols(annotated_peaks))) {
        distance <- abs(mcols(annotated_peaks)$distanceToSite)
        annotated_peaks$distance_score <- .calculateDistanceScore(
            distance, max_distance, weight, near_decay_rate, far_decay_rate
        )
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + annotated_peaks$distance_score
    }
    return(annotated_peaks)
}

# Helper function: Apply Tier 3 - Overlap quality (Jaccard index)
.applyOverlapScore <- function(annotated_peaks, weight) {
    jaccard <- .getJaccardIndex(annotated_peaks)
    overlap_score <- ifelse(jaccard > 0, weight * jaccard, 0)
    annotated_peaks$priority_score <- 
        annotated_peaks$priority_score + overlap_score
    return(annotated_peaks)
}

# Helper function: Apply biotype bonus
.applyBiotypeBonus <- function(annotated_peaks, biotype_priority = NULL, weight = NULL) {
    if ("gene_biotype" %in% colnames(mcols(annotated_peaks))) {
        if (is.null(biotype_priority)) {
            # Simple binary: protein_coding gets bonus
            if (is.null(weight)) {
                stop("Either 'biotype_priority' or 'weight' must be provided", call. = FALSE)
            }
            biotype_score <- ifelse(
                mcols(annotated_peaks)$gene_biotype == "protein_coding", weight, 0
            )
        } else {
            # Use priority mapping (weight is ignored if biotype_priority is provided)
            biotype_score <- ifelse(
                mcols(annotated_peaks)$gene_biotype %in% names(biotype_priority),
                biotype_priority[mcols(annotated_peaks)$gene_biotype], 0
            )
            biotype_score[is.na(biotype_score)] <- 0
        }
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + biotype_score
    }
    return(annotated_peaks)
}

#' Prioritize TF annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param max_promoter_distance Maximum distance for promoter scoring
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeTFAnnotations <- function(annotated_peaks, max_promoter_distance = 2000L) {
    # Tier 1: Regulatory relationship (weight: 10000)
    feature_priority <- c(
        "overlapStart" = 10000,
        "upstream" = 9500,
        "inside" = 8000,
        "overlapEnd" = 6000,
        "overlap" = 5000,
        "downstream" = 4000,
        "includeFeature" = 3000
    )
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Distance to TSS (weight: 1000, two-part exponential decay)
    annotated_peaks <- .applyDistanceScore(
        annotated_peaks, 
        max_distance = max_promoter_distance,
        weight = 1000,
        near_decay_rate = 500,
        far_decay_rate = 2000
    )
    
    # Tier 3: Overlap quality (weight: 100)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 100)
    
    # Tier 4: Feature type preference (weight: 10)
    annotated_peaks <- .applyBiotypeBonus(annotated_peaks, biotype_priority = NULL, weight = 10)
    
    return(annotated_peaks)
}

#' Prioritize promoter histone mark annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param mark_type "H3K4me3" or "H3K4me2"
#' @param max_promoter_distance Maximum distance for promoter scoring
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizePromoterHistoneAnnotations <- function(annotated_peaks,
                                                  mark_type = c("H3K4me3", "H3K4me2"),
                                                  max_promoter_distance = 2000L) {
    mark_type <- match.arg(mark_type)
    
    # Tier 1: Promoter relationship (weight: 10000)
    if (mark_type == "H3K4me3") {
        feature_priority <- c(
            "overlapStart" = 10000,
            "upstream" = 9500,
            "inside" = 7000,
            "overlapEnd" = 5000,
            "overlap" = 4000,
            "downstream" = 3000,
            "includeFeature" = 2000
        )
    } else {  # H3K4me2
        feature_priority <- c(
            "overlapStart" = 10000,
            "upstream" = 9500,
            "inside" = 8500,
            "overlapEnd" = 6000,
            "overlap" = 5000,
            "downstream" = 4000,
            "includeFeature" = 3000
        )
    }
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Distance to TSS (weight: 800, two-part exponential decay)
    annotated_peaks <- .applyDistanceScore(
        annotated_peaks,
        max_distance = max_promoter_distance,
        weight = 800,
        near_decay_rate = 500,
        far_decay_rate = 3000
    )
    
    # Tier 3: Overlap quality (weight: 200)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 200)
    
    return(annotated_peaks)
}

#' Prioritize enhancer histone mark annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param mark_type "H3K27ac" or "H3K4me1"
#' @param max_enhancer_distance Maximum distance for enhancer scoring
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeEnhancerHistoneAnnotations <- function(annotated_peaks,
                                                   mark_type = c("H3K27ac", "H3K4me1"),
                                                   max_enhancer_distance = 5000L) {
    mark_type <- match.arg(mark_type)
    
    # Tier 1: Enhancer relationship (weight: 10000)
    # H3K27ac: active enhancers/promoters (prioritizes promoters)
    # H3K4me1: poised/active enhancers (prioritizes intergenic enhancers)
    if (mark_type == "H3K27ac") {
        feature_priority <- c(
            "overlapStart" = 10000,   # Active promoters
            "upstream" = 9500,        # Enhancers upstream
            "inside" = 8000,          # Gene body (active genes)
            "overlap" = 7000,
            "downstream" = 6000,      # Enhancers downstream
            "overlapEnd" = 5000,
            "includeFeature" = 4000
        )
    } else {  # H3K4me1
        feature_priority <- c(
            "upstream" = 10000,       # Enhancers (primary)
            "downstream" = 9500,      # Enhancers downstream
            "overlapStart" = 9000,    # Promoter-proximal enhancers
            "inside" = 8000,          # Gene body enhancers
            "overlap" = 7000,
            "overlapEnd" = 6000,
            "includeFeature" = 4000
        )
    }
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Distance to regulatory site (weight: 500, two-part exponential decay)
    # For enhancers, distance is less critical than for promoters
    # Enhancers can be far from genes (>10kb), so we use a slower decay rate
    annotated_peaks <- .applyDistanceScore(
        annotated_peaks,
        max_distance = max_enhancer_distance,
        weight = 500,
        near_decay_rate = 1000,
        far_decay_rate = 10000
    )
    
    # Tier 3: Overlap quality (weight: 300, conditional)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 300)
    
    # Tier 4: Intergenic bonus for enhancer marks (weight: 50)
    # Both H3K27ac and H3K4me1 can be found at intergenic enhancers
    # H3K4me1 is more enhancer-specific, but H3K27ac also marks active enhancers
    # The bonus helps prioritize intergenic enhancer annotations for both marks
    intergenic_bonus <- ifelse(
        mcols(annotated_peaks)$insideFeature %in% c("upstream", "downstream"), 50, 0
    )
    annotated_peaks$priority_score <- 
        annotated_peaks$priority_score + intergenic_bonus
    
    return(annotated_peaks)
}

#' Prioritize gene body histone mark annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param mark_type "H3K36me3" or "H3K27me3"
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeGeneBodyHistoneAnnotations <- function(annotated_peaks,
                                                  mark_type = c("H3K36me3", "H3K27me3")) {
    mark_type <- match.arg(mark_type)
    
    # Tier 1: Gene body relationship (weight: 10000)
    feature_priority <- c(
        "inside" = 10000,
        "includeFeature" = 9500,
        "overlap" = 9000,
        "overlapStart" = 7000,
        "overlapEnd" = 7500,
        "upstream" = 5000,
        "downstream" = 4500
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Overlap quality (weight: 5000, PRIMARY for gene body marks)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 5000)
    
    # Tier 3: Feature type (weight: 200)
    if (mark_type == "H3K36me3") {
        biotype_priority <- c(
            "protein_coding" = 200,
            "known" = 150,
            "lncRNA" = 100,
            "pseudogene" = 50
        )
    } else {  # H3K27me3
        biotype_priority <- c(
            "protein_coding" = 200,
            "lncRNA" = 150,
            "pseudogene" = 100,
            "known" = 150
        )
    }
    annotated_peaks <- .applyBiotypeBonus(annotated_peaks, biotype_priority, weight = NULL)
    
    # Tier 4: Distance to feature center (weight: 50)
    # Calculate center score if feature information is available
    if ("feature.ranges" %in% colnames(mcols(annotated_peaks))) {
        feature_ranges <- mcols(annotated_peaks)$feature.ranges
        if (inherits(feature_ranges, "GRanges")) {
            peak_centers <- start(annotated_peaks) + width(annotated_peaks) / 2
            feature_centers <- start(feature_ranges) + width(feature_ranges) / 2
            feature_widths <- width(feature_ranges)
            
            distance_to_center <- abs(peak_centers - feature_centers)
            normalized_distance <- ifelse(
                feature_widths > 0,
                distance_to_center / feature_widths,
                Inf
            )
            
            center_score <- ifelse(
                normalized_distance <= 0.05,
                50,
                ifelse(
                    normalized_distance <= 1.0,
                    50 * exp(-normalized_distance * 2),
                    50 * exp(-2.0) * exp(-(normalized_distance - 1.0) * 0.5)
                )
            )
            center_score[is.na(center_score)] <- 0
            center_score[center_score < 0] <- 0
            annotated_peaks$priority_score <- 
                annotated_peaks$priority_score + center_score
        }
    }
    
    return(annotated_peaks)
}

#' Prioritize Pol II annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param sequencing_method "ChIP-seq", "GRO-seq", or "PRO-seq"
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizePolIIAnnotations <- function(annotated_peaks,
                                        sequencing_method = c("ChIP-seq", "GRO-seq", "PRO-seq")) {
    sequencing_method <- match.arg(sequencing_method)
    
    # Tier 1: Gene body relationship (weight: 10000)
    feature_priority <- c(
        "inside" = 10000,
        "overlapStart" = 9500,
        "includeFeature" = 9000,
        "overlap" = 8500,
        "overlapEnd" = 8000,
        "upstream" = 6000,
        "downstream" = 5000
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Strand matching (weight: 5000, CRITICAL for RNA-based methods)
    if (sequencing_method %in% c("GRO-seq", "PRO-seq")) {
        strand_match <- (strand(annotated_peaks) == mcols(annotated_peaks)$feature.strand) |
                       (strand(annotated_peaks) == "*") |
                       (mcols(annotated_peaks)$feature.strand == "*")
        annotated_peaks$strand_score <- ifelse(strand_match, 5000, 0)
    } else {
        annotated_peaks$strand_score <- 0
    }
    annotated_peaks$priority_score <- 
        annotated_peaks$priority_score + annotated_peaks$strand_score
    
    # Tier 3: Overlap quality (weight: 2000)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 2000)
    
    # Tier 4: Distance to TSS (weight: 500, two-part exponential decay)
    annotated_peaks <- .applyDistanceScore(
        annotated_peaks,
        max_distance = 1000,
        weight = 500,
        near_decay_rate = 500,
        far_decay_rate = 5000
    )
    
    # Tier 5: Feature type (weight: 100)
    annotated_peaks <- .applyBiotypeBonus(annotated_peaks, biotype_priority = NULL, weight = 100)
    
    return(annotated_peaks)
}

#' Prioritize RBP annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param sequencing_method "ChIP-seq", "CLIP-seq", "iCLIP", or "eCLIP"
#' @param use_transcript_level Logical, prefer transcript-level annotation
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeRBPAnnotations <- function(annotated_peaks,
                                       sequencing_method = c("ChIP-seq", "CLIP-seq", "iCLIP", "eCLIP"),
                                       use_transcript_level = TRUE) {
    sequencing_method <- match.arg(sequencing_method)
    
    # Tier 1: Gene body relationship (weight: 10000)
    feature_priority <- c(
        "inside" = 10000,
        "overlap" = 9500,
        "includeFeature" = 9000,
        "overlapStart" = 7000,
        "overlapEnd" = 7500,
        "upstream" = 5000,
        "downstream" = 4500
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Strand matching (weight: 5000, CRITICAL for CLIP-seq)
    if (sequencing_method %in% c("CLIP-seq", "iCLIP", "eCLIP")) {
        strand_match <- (strand(annotated_peaks) == mcols(annotated_peaks)$feature.strand) |
                       (strand(annotated_peaks) == "*") |
                       (mcols(annotated_peaks)$feature.strand == "*")
        annotated_peaks$strand_score <- ifelse(strand_match, 5000, 0)
    } else {
        annotated_peaks$strand_score <- 0
    }
    annotated_peaks$priority_score <- 
        annotated_peaks$priority_score + annotated_peaks$strand_score
    
    # Tier 3: Overlap quality (weight: 3000)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 3000)
    
    # Tier 4: Feature type preference (weight: 500)
    if (use_transcript_level && "feature_type" %in% colnames(mcols(annotated_peaks))) {
        transcript_bonus <- ifelse(
            mcols(annotated_peaks)$feature_type == "transcript", 500, 0
        )
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + transcript_bonus
    }
    
    # Tier 5: Exon preference (weight: 200)
    if ("feature_subtype" %in% colnames(mcols(annotated_peaks))) {
        exon_bonus <- ifelse(
            mcols(annotated_peaks)$feature_subtype == "exon", 200, 0
        )
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + exon_bonus
    }
    
    return(annotated_peaks)
}

#' Prioritize 3' end annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param max_tes_distance Maximum distance to TES
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeThreeEndAnnotations <- function(annotated_peaks, max_tes_distance = 3000L) {
    # Tier 1: 3' end relationship (weight: 10000)
    feature_priority <- c(
        "overlapEnd" = 10000,
        "downstream" = 9000,
        "inside" = 7000,
        "overlap" = 6000,
        "overlapStart" = 4000,
        "upstream" = 3000,
        "includeFeature" = 5000
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Distance to TES (weight: 1000, two-part exponential decay)
    annotated_peaks <- .applyDistanceScore(
        annotated_peaks,
        max_distance = max_tes_distance,
        weight = 1000,
        near_decay_rate = 500,
        far_decay_rate = 2000
    )
    
    # Tier 3: 3' UTR preference (weight: 500)
    if ("feature_subtype" %in% colnames(mcols(annotated_peaks))) {
        utr3_bonus <- ifelse(mcols(annotated_peaks)$feature_subtype == "utr3", 500, 0)
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + utr3_bonus
    }
    
    return(annotated_peaks)
}

#' Prioritize exon-specific annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeExonSpecificAnnotations <- function(annotated_peaks) {
    # Tier 1: Exon relationship (weight: 10000)
    feature_priority <- c(
        "inside" = 10000,
        "overlap" = 9500,
        "overlapStart" = 9000,
        "overlapEnd" = 9000,
        "includeFeature" = 8500,
        "upstream" = 5000,
        "downstream" = 4500
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Exon subtype preference (weight: 5000)
    if ("feature_subtype" %in% colnames(mcols(annotated_peaks))) {
        exon_bonus <- ifelse(mcols(annotated_peaks)$feature_subtype == "exon", 5000, 0)
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + exon_bonus
    }
    
    # Tier 3: Overlap quality (weight: 500)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 500)
    
    return(annotated_peaks)
}

#' Prioritize intron-specific annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeIntronSpecificAnnotations <- function(annotated_peaks) {
    # Tier 1: Intron relationship (weight: 10000)
    feature_priority <- c(
        "inside" = 10000,
        "overlap" = 9500,
        "includeFeature" = 9000,
        "overlapStart" = 7000,
        "overlapEnd" = 7000,
        "upstream" = 5000,
        "downstream" = 4500
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Intron subtype preference (weight: 5000)
    if ("feature_subtype" %in% colnames(mcols(annotated_peaks))) {
        intron_bonus <- ifelse(mcols(annotated_peaks)$feature_subtype == "intron", 5000, 0)
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + intron_bonus
    }
    
    # Tier 3: Overlap quality (weight: 500)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 500)
    
    return(annotated_peaks)
}

#' Prioritize architectural feature annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeArchitecturalAnnotations <- function(annotated_peaks) {
    # Tier 1: Relationship type (weight: 10000)
    feature_priority <- c(
        "overlap" = 10000,
        "includeFeature" = 9500,
        "inside" = 9000,
        "overlapStart" = 8000,
        "overlapEnd" = 8000,
        "upstream" = 6000,
        "downstream" = 5500
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Overlap quality (weight: 5000)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 5000)
    
    return(annotated_peaks)
}

#' Prioritize non-coding RNA annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeNonCodingRNAAnnotations <- function(annotated_peaks) {
    # Tier 1: Relationship type (weight: 10000)
    feature_priority <- c(
        "inside" = 10000,
        "overlap" = 9500,
        "includeFeature" = 9000,
        "overlapStart" = 8000,
        "overlapEnd" = 8000,
        "upstream" = 6000,
        "downstream" = 5500
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Non-coding RNA type preference (weight: 3000)
    if ("gene_biotype" %in% colnames(mcols(annotated_peaks))) {
        ncRNA_types <- c("lncRNA", "miRNA", "snoRNA", "snRNA", "rRNA", "tRNA")
        ncRNA_bonus <- ifelse(
            mcols(annotated_peaks)$gene_biotype %in% ncRNA_types, 3000, 0
        )
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + ncRNA_bonus
    }
    
    # Tier 3: Overlap quality (weight: 1000)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 1000)
    
    return(annotated_peaks)
}

#' Prioritize repetitive element annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeRepetitiveElementAnnotations <- function(annotated_peaks) {
    # Tier 1: Relationship type (weight: 10000)
    feature_priority <- c(
        "overlap" = 10000,
        "inside" = 9500,
        "includeFeature" = 9000,
        "overlapStart" = 8000,
        "overlapEnd" = 8000,
        "upstream" = 6000,
        "downstream" = 5500
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Repetitive element type preference (weight: 5000)
    if ("repeat_type" %in% colnames(mcols(annotated_peaks))) {
        repeat_bonus <- ifelse(
            !is.na(mcols(annotated_peaks)$repeat_type), 5000, 0
        )
        annotated_peaks$priority_score <- 
            annotated_peaks$priority_score + repeat_bonus
    }
    
    # Tier 3: Overlap quality (weight: 2000)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 2000)
    
    return(annotated_peaks)
}

#' Prioritize intergenic annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param max_intergenic_distance Maximum distance for intergenic scoring
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeIntergenicAnnotations <- function(annotated_peaks, max_intergenic_distance = 50000L) {
    # Tier 1: Intergenic relationship (weight: 10000)
    feature_priority <- c(
        "upstream" = 10000,
        "downstream" = 9500,
        "overlap" = 8000,
        "overlapStart" = 7000,
        "overlapEnd" = 7000,
        "inside" = 6000,
        "includeFeature" = 5000
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Distance (weight: 5000, important for intergenic, two-part exponential decay)
    annotated_peaks <- .applyDistanceScore(
        annotated_peaks,
        max_distance = max_intergenic_distance,
        weight = 5000,
        near_decay_rate = 5000,
        far_decay_rate = 10000
    )
    
    # Tier 3: Overlap quality (weight: 1000, conditional)
    # Most intergenic peaks don't overlap genes
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 1000)
    
    return(annotated_peaks)
}

#' Prioritize chromatin remodeler annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeChromatinRemodelerAnnotations <- function(annotated_peaks) {
    # Tier 1: Relationship type (weight: 5000, less important)
    feature_priority <- c(
        "includeFeature" = 5000,
        "inside" = 4500,
        "overlap" = 4000,
        "overlapStart" = 3500,
        "overlapEnd" = 3500,
        "upstream" = 3000,
        "downstream" = 3000
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Overlap quality (weight: 10000, PRIMARY)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 10000)
    
    # Tier 3: Feature size (weight: 1000)
    if ("feature_width" %in% colnames(mcols(annotated_peaks))) {
        # Prefer larger genes (more likely to be target)
        feature_width <- mcols(annotated_peaks)$feature_width
        max_width <- max(feature_width, na.rm = TRUE)
        if (max_width > 0) {
            size_score <- 1000 * (feature_width / max_width)
            size_score[is.na(size_score)] <- 0
            annotated_peaks$priority_score <- 
                annotated_peaks$priority_score + size_score
        }
    }
    
    return(annotated_peaks)
}

#' Prioritize multimodal annotations
#' 
#' @param annotated_peaks GRanges with annotations
#' @param ... Additional parameters
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeMultiModalAnnotations <- function(annotated_peaks, ...) {
    # For multimodal factors, use a balanced approach
    # Tier 1: Relationship type (weight: 10000)
    feature_priority <- c(
        "overlapStart" = 10000,
        "inside" = 9500,
        "overlap" = 9000,
        "includeFeature" = 8500,
        "overlapEnd" = 8000,
        "upstream" = 7000,
        "downstream" = 6500
    )
    
    annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
    
    # Tier 2: Overlap quality (weight: 3000)
    annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 3000)
    
    # Tier 3: Distance (weight: 2000, two-part exponential decay)
    annotated_peaks <- .applyDistanceScore(
        annotated_peaks,
        max_distance = 5000,
        weight = 2000,
        near_decay_rate = 1000,
        far_decay_rate = 5000
    )
    
    return(annotated_peaks)
}

#' Unified prioritization function
#' 
#' @param annotated_peaks GRanges with annotations
#' @param factor_type Factor type
#' @param sequencing_method Sequencing method (for PolII, RBP)
#' @param ... Additional parameters
#' @return GRanges with priority_score column
#' @keywords internal
.prioritizeAnnotations <- function(annotated_peaks,
                                  factor_type,
                                  sequencing_method = NULL,
                                  ...) {
    # Map factor type to prioritization function
    switch(factor_type,
        "TF" = .prioritizeTFAnnotations(annotated_peaks, ...),
        "H3K4me3" = .prioritizePromoterHistoneAnnotations(annotated_peaks, "H3K4me3", ...),
        "H3K4me2" = .prioritizePromoterHistoneAnnotations(annotated_peaks, "H3K4me2", ...),
        "H3K27ac" = .prioritizeEnhancerHistoneAnnotations(annotated_peaks, "H3K27ac", ...),
        "H3K4me1" = .prioritizeEnhancerHistoneAnnotations(annotated_peaks, "H3K4me1", ...),
        "H3K36me3" = .prioritizeGeneBodyHistoneAnnotations(annotated_peaks, "H3K36me3"),
        "H3K27me3" = .prioritizeGeneBodyHistoneAnnotations(annotated_peaks, "H3K27me3"),
        "PolII" = .prioritizePolIIAnnotations(annotated_peaks, sequencing_method %||% "ChIP-seq"),
        "RBP" = .prioritizeRBPAnnotations(annotated_peaks, sequencing_method %||% "ChIP-seq", ...),
        "3end" = .prioritizeThreeEndAnnotations(annotated_peaks, ...),
        "exon" = .prioritizeExonSpecificAnnotations(annotated_peaks),
        "intron" = .prioritizeIntronSpecificAnnotations(annotated_peaks),
        "architectural" = .prioritizeArchitecturalAnnotations(annotated_peaks),
        "ncRNA" = .prioritizeNonCodingRNAAnnotations(annotated_peaks),
        "repeat" = .prioritizeRepetitiveElementAnnotations(annotated_peaks),
        "intergenic" = .prioritizeIntergenicAnnotations(annotated_peaks, ...),
        "chromatin_remodeler" = .prioritizeChromatinRemodelerAnnotations(annotated_peaks),
        "multimodal" = .prioritizeMultiModalAnnotations(annotated_peaks, ...),
        # Default: generic prioritization
        {
            warning("Unknown factor type '", factor_type, 
                   "'. Using generic prioritization.", call. = FALSE)
            # Generic prioritization (simplified)
            feature_priority <- c(
                "overlapStart" = 10000,
                "upstream" = 9500,
                "inside" = 8000,
                "overlapEnd" = 6000,
                "overlap" = 5000,
                "downstream" = 4000,
                "includeFeature" = 3000
            )
            annotated_peaks <- .applyFeaturePriority(annotated_peaks, feature_priority)
            
            # Tier 2: Distance (weight: 1000, two-part exponential decay)
            annotated_peaks <- .applyDistanceScore(
                annotated_peaks,
                max_distance = 2000,
                weight = 1000,
                near_decay_rate = 500,
                far_decay_rate = 2000
            )
            
            # Tier 3: Overlap quality (weight: 100)
            annotated_peaks <- .applyOverlapScore(annotated_peaks, weight = 100)
            
            annotated_peaks
        }
    )
}

