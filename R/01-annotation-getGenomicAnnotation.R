# Helper function to get detailed genomic annotation
# This function determines detailed genomic features (promoter, UTR, exon, intron, etc.)
# and prioritizes them for each peak

#' Get detailed genomic annotation for peaks
#' 
#' Determines detailed genomic annotation features (promoter regions, UTRs,
#' exons, introns, downstream, intergenic) for peaks, similar to ChIPseeker's
#' \code{annotatePeak()} function. This function requires a TxDb or EnsDb
#' object to extract genomic features.
#' 
#' @param peaks A \code{GRanges} object containing peaks to be annotated.
#' @param TxDb An object of class \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#'        or \code{\link[ensembldb:EnsDb-class]{EnsDb}} containing genomic
#'        annotation. Required for extracting detailed genomic features. If the
#'        seqlevels style of the TxDb is not the same as the peaks, the function
#'        will try to match the seqlevels style of the TxDb to the peaks. 
#'        However, due to the different naming conventions of seqlevels, 
#'        the function may not be able to match the seqlevels style of
#'        the TxDb to the peaks perfectly. Therefore, it is strongly RECOMMENDED
#'        to use the same seqlevels style as the peaks when using this function.
#' @param tssRegion A numeric vector of length 2 specifying the TSS region
#'        boundaries. Default is \code{c(-3000, 3000)}. The first element is
#'        upstream distance, second is downstream distance. Used to define
#'        promoter regions. Promoter regions are dynamically binned based on
#'        this parameter (e.g., "Promoter (<=1kb)", "Promoter (1-2kb)", etc.).
#' @param immediateDownstreamLength Integer. Length in base pairs for immediate
#'        downstream regions from gene end. Default is \code{3000L}. Peaks within
#'        this distance downstream of a gene end are annotated as "Immediate Downstream".
#' @param level Character string. One of "transcript" or "gene". Default is
#'        "transcript". Determines whether to use transcript-level or gene-level
#'        features.
#' @param genomicAnnotationPriority A character vector specifying the priority
#'        order for genomic annotations. Default is:
#'        \code{c("promoter", "fiveUTR", "firstExon", "firstIntron", "threeUTR", 
#'        "otherExon", "otherIntron", "immediateDownstream", "distalIntergenic")}.
#'        Valid annotation types are: "promoter", "fiveUTR", "threeUTR", "firstExon",
#'        "otherExon", "firstIntron", "otherIntron", "immediateDownstream", 
#'        "distalIntergenic". Annotations are assigned based on this priority order
#'        when a peak overlaps multiple features. Higher priority features (earlier
#'        in the vector) take precedence over lower priority features.
#' @param ignore_strand Logical. If \code{TRUE}, strand information is ignored
#'        when finding overlaps. Default is \code{FALSE}. If all peaks have
#'        strand "*", this is automatically set to \code{TRUE}.
#' @param usePeakCenter Logical. If \code{TRUE} (default), uses peak centers
#'        (midpoints) as reference points for overlap checking. This is more
#'        precise, especially for large peaks that span multiple genomic regions,
#'        as it represents where the actual binding event (peak summit) is
#'        located. If \code{FALSE}, checks if any part of the peak overlaps with
#'        features, which is the traditional approach but may be less precise for
#'        large peaks.
#' 
#' @return A list with two elements:
#'        \itemize{
#'          \item \code{peaks_annotated_with_overlapping_features}: A
#'                \link[GenomicRanges:GRanges-class]{GRanges} object containing
#'                the input peaks with two additional metadata columns:
#'                \itemize{
#'                  \item \code{annotation}: Character vector containing all
#'                        overlapping annotations for each peak (comma-separated).
#'                        Possible values include Promoter bins (e.g., "Promoter (<=1kb)",
#'                        "Promoter (1-2kb)"), "5' UTR", "3' UTR", "1st Exon",
#'                        "Other Exon", "1st Intron", "Other Intron", "Immediate Downstream",
#'                        and "Distal Intergenic".
#'                  \item \code{priorityAnnotation}: Character vector containing
#'                        the highest priority annotation for each peak based on
#'                        \code{genomicAnnotationPriority}. Only one annotation
#'                        type per peak (the highest priority one).
#'                }
#'          \item \code{plot_priority_annotation}: A plot object (donut chart)
#'                visualizing the distribution of priority annotations across all peaks.
#'        }
#' 
#' @details
#' This function provides detailed genomic annotation similar to ChIPseeker's
#' \code{annotatePeak()} function. Key features:
#' \itemize{
#'   \item \strong{Flexible overlap checking}: By default, uses peak centers
#'         (midpoints) as reference points for overlap checking, which is more
#'         precise for large peaks. Can be changed to use full peak ranges via
#'         the \code{usePeakCenter} parameter.
#'   \item \strong{Dynamic promoter binning}: Promoter regions are automatically
#'         binned based on \code{tssRegion}. For example, if \code{tssRegion = c(-3000, 3000)},
#'         bins are created at 1kb, 2kb, and 3kb, resulting in "Promoter (<=1kb)",
#'         "Promoter (1-2kb)", and "Promoter (2-3kb)" annotations.
#'   \item \strong{Per-transcript handling}: First exons and introns are identified
#'         per transcript, ensuring that each transcript's first exon/intron is
#'         correctly recognized. This is important because the same genomic region
#'         may be the first exon of one transcript but a later exon of another.
#'   \item \strong{Priority-based annotation}: When a peak center overlaps multiple
#'         features, the annotation is assigned based on \code{genomicAnnotationPriority},
#'         with higher priority features taking precedence. Features are checked in
#'         reverse priority order (lowest to highest) for efficiency.
#'   \item \strong{Visualization}: Automatically generates a donut chart visualizing
#'         the distribution of priority annotations across all peaks using the
#'         \code{\link{donut}} function.
#' }
#' 
#' @importFrom GenomicFeatures exons intronsByTranscript fiveUTRsByTranscript
#'             threeUTRsByTranscript transcripts promoters exonsBy genes
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom S4Vectors queryHits
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' data(peaks1)
#' result <- getGenomicAnnotation(peaks1, txdb)
#' # Access annotated peaks
#' annotated_peaks <- result$peaks_annotated_with_overlapping_features
#' head(mcols(annotated_peaks)[, c("annotation", "priorityAnnotation")])
#' # View the plot
#' result$plot_priority_annotation
#' @export
getGenomicAnnotation <- function(peaks, TxDb, tssRegion = c(-3000, 3000),
                                 immediateDownstreamLength = 3000,
                                 level = c("transcript", "gene"),
                                 genomicAnnotationPriority = c("promoter", 
                                                            "fiveUTR", 
                                                            "firstExon",  
                                                            "firstIntron", 
                                                            "threeUTR", 
                                                            "otherExon", 
                                                            "otherIntron", 
                                                            "immediateDownstream",
                                                            "distalIntergenic"),
                                 ignore_strand = TRUE,
                                 usePeakCenter = TRUE) {   
    if (!inherits(TxDb, c("TxDb", "EnsDb"))) {
        stop("TxDb must be a TxDb or EnsDb object", call. = FALSE)
    }
    if (!inherits(peaks, "GRanges")) {
        stop("peaks must be a GRanges object", call. = FALSE)
    }

    level <- match.arg(level)
    
    # check genomic annotation priority
    if (!all(genomicAnnotationPriority %in% c("promoter", "fiveUTR", "firstExon", "firstIntron", "threeUTR", "otherExon", "otherIntron", "immediateDownstream", "distalIntergenic"))) {
        warning("Invalid priority annotation types ignored: ", 
               paste(genomicAnnotationPriority[!genomicAnnotationPriority %in% c("promoter", "fiveUTR", "firstExon", "firstIntron", "threeUTR", "otherExon", "otherIntron", "immediateDownstream", "distalIntergenic")], collapse = ", "), call. = FALSE)
        genomicAnnotationPriority <- genomicAnnotationPriority[genomicAnnotationPriority %in% c("promoter", "fiveUTR", "firstExon", "firstIntron", "threeUTR", "otherExon", "otherIntron", "immediateDownstream", "distalIntergenic")]
    }
    
    # match seqlevels styles between peaks and TxDb if needed
    GenomeInfoDb::seqlevelsStyle(TxDb) <- GenomeInfoDb::seqlevelsStyle(peaks)[1]
    
    # remove seqlevels that are not in peaks
    seqlevels(TxDb) <- intersect(seqlevels(TxDb), seqlevels(peaks))

    # Extract genomic features
    if (level == "transcript") {
        transcripts <- transcripts(TxDb, columns = NULL)
    } else {
        transcripts <- genes(TxDb, columns = NULL)
    }
    
    # Get UTRs
    five_utrs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
    three_utrs <- unique(unlist(threeUTRsByTranscript(TxDb)))
    
    # Bin promoter regions dynamically based on tssRegion
    max_flank_width <- max(abs(tssRegion))  # Maximum promoter distance (e.g., 3000 for c(-3000, 3000))
    
    # Create promoter bins dynamically
    promoter_regions <- list()
    promoter_labels <- character()
    
  
    kb_bins <- as.integer(seq(1, ceiling(max_flank_width / 1000)) * 1000)
    
    for (i in seq_along(kb_bins)) {
        upstream_dist <- kb_bins[i]
        if (i == 1) {
            # First bin: <= 1kb
            promoter_gr <- suppressWarnings(unique(promoters(transcripts, upstream = upstream_dist, 
                                            downstream = upstream_dist)))
            promoter_gr <- trim(promoter_gr)
            promoter_regions[[i]] <- promoter_gr
            promoter_labels[i] <- paste0("Promoter (<=1kb)")
        } else {
            # Subsequent bins: (i-1)kb - ikb
            # Note: These are nested (2kb includes 1kb, 3kb includes 2kb, etc.)
            # We'll process from smallest to largest and only assign to smallest matching bin
            promoter_gr <- suppressWarnings(unique(promoters(transcripts, upstream = upstream_dist, 
                                            downstream = upstream_dist)))
            promoter_gr <- trim(promoter_gr)
            promoter_regions[[i]] <- promoter_gr
            promoter_labels[i] <- paste0("Promoter (", (i - 1), "-", i, "kb)")
        }
    }
    
    # Define immediate downstream regions (within immediateDownstreamLength of gene end by default)
    immediate_downstream_gr <- suppressWarnings(unique(downstreams(transcripts, upstream = 0L, 
                                                  downstream = immediateDownstreamLength)))
    immediate_downstream_gr <- trim(immediate_downstream_gr)
    
    # Identify first exons and introns per transcript
    # For each transcript, get its exons and introns in order
    exons_by_tx <- exonsBy(TxDb, by = "tx", use.names = TRUE)
    introns_by_tx <- intronsByTranscript(TxDb, use.names = TRUE)
    
    # Get first exons and introns
    exons_by_tx_combined <- unlist(exons_by_tx)
    first_exons_index <- mcols(exons_by_tx_combined)$exon_rank == 1
    first_exons <- unique(exons_by_tx_combined[first_exons_index])
    all_exons <- unique(exons_by_tx_combined)

    # Get first introns
    introns_by_tx_combined <- unlist(introns_by_tx)
    first_introns_index <- mcols(introns_by_tx_combined)$intron_rank == 1
    first_introns <- unique(introns_by_tx_combined[first_introns_index])
    all_introns <- unique(introns_by_tx_combined)

    # Other exons and introns (excluding first)
    other_exons <- unique(exons_by_tx_combined[!first_exons_index])
    other_introns <- unique(introns_by_tx_combined[!first_introns_index])
    
    # Define mapping from priority names to feature data
    # This allows dynamic ordering based on genomicAnnotationPriority
    feature_mapping <- list(
        "promoter" = list(
            regions = promoter_regions,
            labels = promoter_labels,
            priority_label = "Promoter"
        ),
        "fiveUTR" = list(
            regions = list(five_utrs),
            labels = "5' UTR"
        ),
        "threeUTR" = list(
            regions = list(three_utrs),
            labels =  "3' UTR"
        ),
        "firstExon" = list(
            regions = list(first_exons),
            labels =  "1st Exon"
        ),
        "otherExon" = list(
            regions = list(other_exons),
            labels = "Other Exon"
        ),
        "firstIntron" = list(
            regions = list(first_introns),
            labels = "1st Intron"
        ),
        "otherIntron" = list(
            regions = list(other_introns),
            labels = "Other Intron"
        ),
        "immediateDownstream" = list(
            regions = list(immediate_downstream_gr),
            labels = "Immediate Downstream"
        ),
        "distalIntergenic" = list(
            regions = list(GRanges()),  # Will be calculated per chromosome
            labels = "Distal Intergenic"
        )
    )
    
    # Get feature mappings in the order specified by genomicAnnotationPriority
    # Filter out any NULL entries (invalid priority names)
    priority_order <- feature_mapping[genomicAnnotationPriority]
    priority_order <- priority_order[!sapply(priority_order, is.null)]
     
    # Reverse the order so we process from lowest to highest priority
    # (higher priority features will overwrite lower priority ones)
    priority_order <- rev(priority_order)
    
    # Determine reference points for overlap checking
    if (usePeakCenter) {
        # Use peak centers (midpoints) as reference points for overlap checking
        # This is more precise, especially for large peaks that span multiple
        # genomic regions. Peak center represents the actual binding event location.
        peaks_for_overlap <- GRanges(
            seqnames = GenomicRanges::seqnames(peaks),
            ranges = IRanges(start = start(peaks) + width(peaks) %/% 2L,
                            width = 1L),
            strand = strand(peaks)
        )
    } else {
        # Use full peak ranges for overlap checking (traditional approach)
        # Checks if any part of the peak overlaps with features
        peaks_for_overlap <- peaks
    }
    
    # Split peaks and annotation features into GRangesList by chromosome
    # This enables efficient chromosome-by-chromosome processing
    peaks_by_chr <- split(peaks_for_overlap, GenomicRanges::seqnames(peaks_for_overlap), drop =TRUE)
    
    # Split all annotation features by chromosome
    # Prepare feature GRangesList objects for efficient filtering
    promoter_regions_by_chr <- lapply(promoter_regions, function(pr) {
        split(pr, GenomicRanges::seqnames(pr))
    })
    five_utrs_by_chr <- split(five_utrs, GenomicRanges::seqnames(five_utrs))
    three_utrs_by_chr <- split(three_utrs, GenomicRanges::seqnames(three_utrs))
    immediate_downstream_by_chr <- split(immediate_downstream_gr,
                                          GenomicRanges::seqnames(immediate_downstream_gr))
    all_exons_by_chr <- split(all_exons, GenomicRanges::seqnames(all_exons))
    all_introns_by_chr <- split(all_introns, GenomicRanges::seqnames(all_introns))
    first_exons_by_chr <- split(first_exons, GenomicRanges::seqnames(first_exons))
    first_introns_by_chr <- split(first_introns, GenomicRanges::seqnames(first_introns))
    other_exons_by_chr <- split(other_exons, GenomicRanges::seqnames(other_exons))
    other_introns_by_chr <- split(other_introns, GenomicRanges::seqnames(other_introns))
    
    # Initialize annotation vectors
    # annotation: will accumulate ALL overlapping annotations (comma-separated)
    # priority_annotation: will contain only the highest priority annotation
    n_peaks <- length(peaks)
    annotation <- character(n_peaks)  # Start empty, will accumulate annotations
    priority_annotation <- character(n_peaks)
    
    # Store original peak indices for mapping back results
    peak_indices <- seq_len(n_peaks)
    peaks_idx_by_chr <- split(peak_indices, GenomicRanges::seqnames(peaks))
    
    # Process each chromosome separately using GRangesList structure
    chr_names <- names(peaks_by_chr)
    for (chr in chr_names) {
        chr_peaks <- peaks_by_chr[[chr]]
        if (length(chr_peaks) == 0L) next
        
        chr_peaks_idx <- peaks_idx_by_chr[[chr]]
        
        # Get chromosome-specific features from GRangesList
        # This reduces memory usage significantly for large genomes
        feature_mapping_chr <- lapply(feature_mapping, function(feat_info) {
            if (is.list(feat_info$regions)) {
                # For promoters (list of regions) or simple features (list with one element)
                regions_chr <- lapply(feat_info$regions, function(reg) {
                    if (chr %in% names(split(reg, GenomicRanges::seqnames(reg)))) {
                        split(reg, GenomicRanges::seqnames(reg))[[chr]]
                    } else {
                        GRanges()
                    }
                })
                list(regions = regions_chr,
                     labels = feat_info$labels,
                     priority_label = feat_info$labels)
            } else {
                # Should not happen, but handle gracefully
                feat_info
            }
        })
        
        # Get chromosome-specific annotation features
        largest_promoter_chr <- if (length(promoter_regions) > 0L) {
            if (chr %in% names(promoter_regions_by_chr[[length(promoter_regions)]])) {
                promoter_regions_by_chr[[length(promoter_regions)]][[chr]]
            } else {
                GRanges()
            }
        } else {
            GRanges()
        }
        all_exons_chr <- if (chr %in% names(all_exons_by_chr)) {
            all_exons_by_chr[[chr]]
        } else {
            GRanges()
        }
        all_introns_chr <- if (chr %in% names(all_introns_by_chr)) {
            all_introns_by_chr[[chr]]
        } else {
            GRanges()
        }
        first_exons_chr <- if (chr %in% names(first_exons_by_chr)) {
            first_exons_by_chr[[chr]]
        } else {
            GRanges()
        }
        first_introns_chr <- if (chr %in% names(first_introns_by_chr)) {
            first_introns_by_chr[[chr]]
        } else {
            GRanges()
        }
        other_exons_chr <- if (chr %in% names(other_exons_by_chr)) {
            other_exons_by_chr[[chr]]
        } else {
            GRanges()
        }
        other_introns_chr <- if (chr %in% names(other_introns_by_chr)) {
            other_introns_by_chr[[chr]]
        } else {
            GRanges()
        }
        five_utrs_chr <- if (chr %in% names(five_utrs_by_chr)) {
            five_utrs_by_chr[[chr]]
        } else {
            GRanges()
        }
        three_utrs_chr <- if (chr %in% names(three_utrs_by_chr)) {
            three_utrs_by_chr[[chr]]
        } else {
            GRanges()
        }
        immediate_downstream_chr <- if (chr %in% names(immediate_downstream_by_chr)) {
            immediate_downstream_by_chr[[chr]]
        } else {
            GRanges()
        }
        
        # Calculate intergenic regions for this chromosome
        all_annotated_chr <- reduce(c(all_exons_chr, all_introns_chr, five_utrs_chr, 
                                     three_utrs_chr, largest_promoter_chr, 
                                     immediate_downstream_chr))
        
        # Calculate intergenic for this chromosome
        # Get chromosome lengths from transcripts object (TxDb doesn't support seqlengths directly)
        chr_seqlength <- seqlengths(transcripts)[chr]
        if (is.na(chr_seqlength)) {
            # If chromosome length is unknown, skip intergenic calculation
            intergenic_chr <- GRanges()
        } else {
            if (ignore_strand) {
                all_annotated_strandless_chr <- all_annotated_chr
                strand(all_annotated_strandless_chr) <- "*"
                all_annotated_strandless_chr <- reduce(trim(all_annotated_strandless_chr))
                intergenic_chr <- gaps(all_annotated_strandless_chr, 
                                      start = 1L, 
                                      end = chr_seqlength)
                intergenic_chr <- intergenic_chr[strand(intergenic_chr) == "*"]
            } else {
                all_annotated_chr <- reduce(trim(all_annotated_chr))
                intergenic_chr <- gaps(all_annotated_chr, 
                                     start = 1L, 
                                     end = chr_seqlength)
                intergenic_chr <- intergenic_chr[strand(intergenic_chr) != "*"]
            }
        }
        
        # Update intergenic in feature mapping for this chromosome
        if ("distalIntergenic" %in% names(feature_mapping_chr)) {
            feature_mapping_chr[["distalIntergenic"]]$regions[[1]] <- intergenic_chr
        }
        
        # Get priority order for this chromosome
        priority_order_chr <- feature_mapping_chr[genomicAnnotationPriority]
        priority_order_chr <- priority_order_chr[!sapply(priority_order_chr, is.null)]
        priority_order_chr <- rev(priority_order_chr)
        
        # Process from lowest to highest priority (reversed order)
        # For annotation: accumulate all overlapping types
        # For priority_annotation: overwrite with highest priority only
        for (i in seq_along(priority_order_chr)) {
            priority_name <- names(priority_order_chr)[i]
            feature_info <- priority_order_chr[[i]]
            
            if (priority_name == "promoter") {
                # Process promoters from largest to smallest bin
                # if there is any overlap, the smallest bin will be the highest priority
                for (j in rev(seq_along(feature_info$regions))) {
                    ol_promoter <- findOverlaps(chr_peaks, 
                                              feature_info$regions[[j]], 
                                              ignore.strand = ignore_strand)
                    idx_promoter_chr <- queryHits(ol_promoter)
                    if (length(idx_promoter_chr) > 0L) {
                        # Map back to original peak indices
                        idx_promoter <- chr_peaks_idx[idx_promoter_chr]
                        
                        # Accumulate annotation (add to existing, comma-separated)
                        # Use vectorized operations for efficiency
                        has_existing <- nchar(annotation[idx_promoter]) > 0L
                        annotation[idx_promoter[has_existing]] <- 
                            paste(annotation[idx_promoter[has_existing]], 
                                 feature_info$labels[j], 
                                 sep = ", ")
                        annotation[idx_promoter[!has_existing]] <- feature_info$labels[j]
                        
                        # Overwrite priority_annotation (highest priority only)
                        priority_annotation[idx_promoter] <- feature_info$priority_label[j]
                    }
                }
                # if multiple promoters are found, the highest priority promoter will be the one with the smallest bin
                promoter_matches <- gregexpr("Promoter", annotation)
                promoter_matches <- lapply(seq_along(promoter_matches), function(i) {
                    match_starts <- unlist(promoter_matches[[i]])
                    if (length(match_starts) > 1L) { # multiple promoters found
                        annotation[i] <- paste0(substr(annotation[i], start = match_starts[1], stop = match_starts[1] - 1), 
                                              substr(annotation[i], start = match_starts[length(match_starts)], 
                                              stop = nchar(annotation[i])))
                    }
                    annotation[i]
                })
                annotation <- unlist(promoter_matches)
            } else {
                # Process simple features (single region type)
                ol_feature <- findOverlaps(chr_peaks, 
                                         feature_info$regions[[1]], 
                                         ignore.strand = ignore_strand)
                idx_feature_chr <- queryHits(ol_feature)
                if (length(idx_feature_chr) > 0L) {
                    # Map back to original peak indices
                    idx_feature <- chr_peaks_idx[idx_feature_chr]
                    
                    # Accumulate annotation (add to existing, comma-separated)
                    # Use vectorized operations for efficiency
                    has_existing <- nchar(annotation[idx_feature]) > 0L
                    annotation[idx_feature[has_existing]] <- 
                        paste(annotation[idx_feature[has_existing]], 
                             feature_info$labels, 
                             sep = ", ")
                    annotation[idx_feature[!has_existing]] <- feature_info$labels
                    
                    # Overwrite priority_annotation (highest priority only)
                    priority_annotation[idx_feature] <- feature_info$labels
                }
            }
        }
    }
    
    # Set default for peaks with no annotations
    annotation[annotation == ""] <- "Distal Intergenic"
    priority_annotation[priority_annotation == ""] <- "Distal Intergenic"
    
    # Return as data.frame
    annotation_df <- data.frame(
        annotation = annotation,
        priorityAnnotation = priority_annotation,
        stringsAsFactors = FALSE
    )
    # add peaks to the annotation_df
    mcols(peaks) <- cbind(mcols(peaks), annotation_df)

    # plot the priority annotation
    p <- donut(table(priority_annotation), labels = names(table(priority_annotation)))

    list(peaks_annotated_with_overlapping_features = peaks, plot_priority_annotation = p)
}

