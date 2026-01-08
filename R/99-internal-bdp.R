#' Obtain peaks near bi-directional promoters
#' 
#' @description 
#' Identifies peaks that are located near bi-directional promoters (BDPs), which
#' are genomic regions where two genes are transcribed in opposite directions from
#' a shared promoter region. A peak is considered near a BDP if it is within
#' \code{maxgap} distance of promoters on both strands.
#' 
#' 
#' @param peaks A \code{GRanges} object containing peak regions.
#' @param annoData An \code{annoGR} or \code{GRanges} object containing
#'   annotation data (e.g., TSS positions).
#' @param maxgap Integer. Maximum distance (in base pairs) between a peak and
#'   TSS to be considered "near". Default is 2000.
#' @param ... Additional arguments (currently not used).
#' @return A \code{GRangesList} object containing peaks near bi-directional
#'   promoters. Each element in the list corresponds to a peak that is within
#'   \code{maxgap} distance of promoters on both the plus and minus strands.
#'   Returns \code{NA} if no peaks meet the criteria.
#' @author Jianhong Ou
#' @seealso See Also as \code{\link{annoPeaks}}, \code{\link{annoGR}}
#' @keywords misc
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors elementNROWS
#' @examples
#' 
#'   if(interactive() || Sys.getenv("USER")=="jou"){
#'     library(ensembldb)
#'     library(EnsDb.Hsapiens.v75)
#'     data("myPeakList")
#'     annoGR <- annoGR(EnsDb.Hsapiens.v75)
#'     seqlevelsStyle(myPeakList) <- seqlevelsStyle(annoGR)
#'     ChIPpeakAnno:::bdp(myPeakList, annoGR)
#'   }
#' 
bdp <- function(peaks, annoData, maxgap = 2000L, ...) {
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    if (length(intersect(seqlevelsStyle(peaks), seqlevelsStyle(annoData))) == 0L) {
        stop("No matching seqlevelsStyle between 'peaks' and 'annoData'. ",
             "Please ensure they use compatible naming conventions.",
             call. = FALSE)
    }
    stopifnot(is.numeric(maxgap))
    maxgap <- round(maxgap[1L])
    peaks <- unique(peaks)
    peaks$bdp_idx <- seq_along(peaks)
    anno <- annoPeaks(peaks, annoData, 
                      bindingType = "nearestBiDirectionalPromoters",
                      bindingRegion = c(-1L * maxgap, maxgap))
    if (length(anno) < 1L) {
        return(NA)
    }
    anno.s <- split(anno, anno$bdp_idx)
    len <- elementNROWS(anno.s)
    anno.s <- anno.s[len >= 2L]
    len <- vapply(anno.s, function(.ele) {
        all(c("+", "-") %in% as.character(.ele$feature.strand))
    }, FUN.VALUE = logical(1))
    anno.s <- anno.s[len]
    anno.s
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
#'   \item Create bidirectional promoter regions (genomic interval between TSSs)
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
    tss_pos <- ifelse(strand_vec == "+", start(annoData), end(annoData))
    
    # Split by strand
    plus_idx <- strand_vec == "+"
    minus_idx <- strand_vec == "-"
    
    if (sum(plus_idx) == 0L || sum(minus_idx) == 0L) {
        # No divergent pairs possible
        return(GRanges())
    }
    
    plus_anno <- annoData[plus_idx]
    minus_anno <- annoData[minus_idx]
    plus_tss <- tss_pos[plus_idx]
    minus_tss <- tss_pos[minus_idx]
    
    # Get gene IDs (from names or metadata) and store original indices
    if (!is.null(names(plus_anno))) {
        plus_ids <- names(plus_anno)
    } else if ("gene_id" %in% colnames(mcols(plus_anno))) {
        plus_ids <- mcols(plus_anno)$gene_id
    } else {
        plus_ids <- paste0("gene_plus_", seq_along(plus_anno))
    }
    
    if (!is.null(names(minus_anno))) {
        minus_ids <- names(minus_anno)
    } else if ("gene_id" %in% colnames(mcols(minus_anno))) {
        minus_ids <- mcols(minus_anno)$gene_id
    } else {
        minus_ids <- paste0("gene_minus_", seq_along(minus_anno))
    }
    
    # Store mapping from plus/minus indices to original annotation indices
    plus_orig_idx <- which(plus_idx)
    minus_orig_idx <- which(minus_idx)
    
    # Find all pairs of divergent TSSs
    # For each chromosome, find pairs where + strand TSS is upstream of - strand TSS
    # and distance between them is <= maxTSSDistance
    bdp_list <- list()
    
    # Process by chromosome for efficiency
    all_chrs <- unique(c(as.character(seqnames(plus_anno)), as.character(seqnames(minus_anno))))
    
    for (chr in all_chrs) {
        plus_chr_idx <- as.character(seqnames(plus_anno)) == chr
        minus_chr_idx <- as.character(seqnames(minus_anno)) == chr
        
        if (sum(plus_chr_idx) == 0L || sum(minus_chr_idx) == 0L) {
            next
        }
        
        plus_chr_tss <- plus_tss[plus_chr_idx]
        minus_chr_tss <- minus_tss[minus_chr_idx]
        plus_chr_ids <- plus_ids[plus_chr_idx]
        minus_chr_ids <- minus_ids[minus_chr_idx]
        
        # For each + strand TSS, find - strand TSSs that are downstream
        # (head-to-head configuration: + strand TSS should be < - strand TSS)
        for (i in seq_along(plus_chr_tss)) {
            # Find - strand TSSs that are downstream (greater coordinate)
            # and within maxTSSDistance
            downstream_minus <- minus_chr_tss >= plus_chr_tss[i] &
                               minus_chr_tss <= plus_chr_tss[i] + maxTSSDistance
            
            if (any(downstream_minus)) {
                for (j in which(downstream_minus)) {
                    tss_dist <- minus_chr_tss[j] - plus_chr_tss[i]
                    
                    # Create bidirectional promoter region (between the two TSSs)
                    bdp_start <- plus_chr_tss[i]
                    bdp_end <- minus_chr_tss[j]
                    bdp_center <- as.integer(round((bdp_start + bdp_end) / 2))
                    
                    # Get original annotation indices
                    plus_chr_orig_idx <- plus_orig_idx[plus_chr_idx]
                    minus_chr_orig_idx <- minus_orig_idx[minus_chr_idx]
                    
                    bdp_list[[length(bdp_list) + 1L]] <- list(
                        seqnames = chr,
                        start = bdp_start,
                        end = bdp_end,
                        gene1_id = plus_chr_ids[i],
                        gene2_id = minus_chr_ids[j],
                        gene1_idx = plus_chr_orig_idx[i],  # Original index in annoData
                        gene2_idx = minus_chr_orig_idx[j],  # Original index in annoData
                        gene1_strand = "+",
                        gene2_strand = "-",
                        TSS1_pos = plus_chr_tss[i],
                        TSS2_pos = minus_chr_tss[j],
                        TSS_distance = tss_dist,
                        bdp_center = bdp_center
                    )
                }
            }
        }
    }
    
    if (length(bdp_list) == 0L) {
        return(GRanges())
    }
    
    # Convert to GRanges
    bdp_gr <- GRanges(
        seqnames = vapply(bdp_list, function(x) x$seqnames, character(1)),
        IRanges(
            start = vapply(bdp_list, function(x) x$start, integer(1)),
            end = vapply(bdp_list, function(x) x$end, integer(1))
        ),
        strand = "*",  # Bidirectional promoters are not strand-specific
        gene1_id = vapply(bdp_list, function(x) x$gene1_id, character(1)),
        gene2_id = vapply(bdp_list, function(x) x$gene2_id, character(1)),
        gene1_idx = vapply(bdp_list, function(x) x$gene1_idx, integer(1)),
        gene2_idx = vapply(bdp_list, function(x) x$gene2_idx, integer(1)),
        gene1_strand = vapply(bdp_list, function(x) x$gene1_strand, character(1)),
        gene2_strand = vapply(bdp_list, function(x) x$gene2_strand, character(1)),
        TSS1_pos = vapply(bdp_list, function(x) x$TSS1_pos, integer(1)),
        TSS2_pos = vapply(bdp_list, function(x) x$TSS2_pos, integer(1)),
        TSS_distance = vapply(bdp_list, function(x) x$TSS_distance, integer(1)),
        bdp_center = vapply(bdp_list, function(x) x$bdp_center, integer(1))
    )
    
    # Remove duplicates (same gene pair)
    bdp_gr <- unique(bdp_gr)
    
    bdp_gr
}
