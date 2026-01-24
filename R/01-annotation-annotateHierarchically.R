#' Hierarchical annotation of peaks with factor-specific prioritization
#' 
#' @description 
#' This function provides a comprehensive annotation system that applies multiple
#' annotation strategies in parallel, then ranks and prioritizes results based
#' on factor-specific binding patterns. Unlike single-method annotation, this
#' approach captures annotations from multiple perspectives (promoters, inside
#' and downstream regions, bidirectional promoters) using transcripts extracted
#' from an EnsDb object, and intelligently prioritizes them based on the
#' biological factor type.
#' 
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Parallel Annotation}: Applies multiple annotation strategies
#'         simultaneously (promoter-focused, upstream to downstream, bidirectional
#'         promoters) using transcripts from an EnsDb object
#'   \item \strong{Factor-Specific Prioritization}: Ranks results using
#'         factor-specific prioritization strategies (TF, histone marks, Pol II,
#'         etc.)
#'   \item \strong{Factor-Specific}: Requires explicit factor type specification
#'         for accurate prioritization
#'   \item \strong{Special Case Handling}: Automatically handles bidirectional
#'         promoters and broad peaks covering multiple features
#' }
#' 
#' @param peaks A \link[GenomicRanges:GRanges-class]{GRanges} object containing
#'        peaks to be annotated.
#' @param EnsDb A \link[ensembldb:EnsDb-class]{EnsDb} object containing 
#'        genome annotation data (genes, transcripts, exons, etc.) for an organism.
#'        If \code{NULL} or missing, the function will use the global EnsDb set by
#'        \code{\link{setChIPpeakAnnoEnsDb}}.
#' @param factor_type A character string specifying the biological factor type.
#'        \strong{Required if using default prioritization} (when
#'        \code{prioritization_function} is \code{NULL}). Options include: "TF", 
#'        "H3K4me3", "H3K4me2", "H3K27ac", "H3K36me3", "H3K27me3", "PolII", "RBP",
#'        (RNA-binding proteins),"ATAC" (accessibility), "3end" (3' end factors), 
#'        "exon" (exon-specific), "intron" (intron-specific), "architectural" 
#'        (architectural features),"intergenic" (intergenic regions), 
#'        "chromatin_remodeler" (chromatin remodelers), "multimodal" (multimodal factors).
#'        Must be specified based on biological knowledge of the factor being studied.
#'        Can be \code{NULL} if providing a custom \code{prioritization_function}.
#' @param is_promoter_binding Logical. If \code{TRUE} (default), treats the
#'        factor as promoter-binding, which enables special handling of
#'        bidirectional promoter annotations. When \code{TRUE}, peaks located in
#'        bidirectional promoter regions will have both divergent gene annotations
#'        preserved as special cases (excluded from priority score-based
#'        ranking). This is particularly relevant for transcription factors (TF),
#'        H3K4me3, H3K4me2, H3K27ac, and other factors that typically bind near
#'        promoters. Set to \code{FALSE} for factors that do not typically bind
#'        at promoters (e.g., H3K36me3, H3K27me3, intergenic factors) to disable
#'        bidirectional promoter special case handling.
#' @param annotation_strategies A list of custom annotation strategies to apply.
#'        If \code{NULL}, default strategies will be used. See Details.
#' @param prioritization_weights A list of custom prioritization weights. If
#'        \code{NULL}, factor-specific default weights will be used.
#' @param prioritization_function A custom prioritization function. If provided,
#'        this function will be used instead of the default factor-specific
#'        prioritization. The function should accept a \code{GRanges} object with
#'        annotations as the first argument and return a \code{GRanges} object
#'        with a \code{priority_score} column. Additional arguments can be passed
#'        via \code{...}. If \code{NULL} (default), uses the built-in
#'        factor-specific prioritization based on \code{factor_type}.
#' @param keep A character string specifying which annotations to return.
#'        Options are:
#'        \itemize{
#'          \item \code{"best"} (default): Returns only the best annotation per
#'                peak (plus special cases) based on priority score.
#'          \item \code{"all"}: Returns all annotations with priority scores.
#'        }
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'        object specifying the parallel backend to use for applying multiple
#'        annotation strategies. If \code{NULL} (default), uses
#'        \code{\link[BiocParallel]{bpparam}()} which automatically detects the
#'        best available backend. For sequential processing, use
#'        \code{BiocParallel::SerialParam()}.
#' @param sequencing_method Sequencing method. For PolII, it can be "ChIP-seq", 
#'        "GRO-seq", "PRO-seq". For RBP, it can be "ChIP-seq", "CLIP-seq", "iCLIP", "eCLIP".
#' @param ... Additional parameters passed to custom prioritization function.
#' 
#' @return Returns a \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing annotated peaks. The object includes:
#'        \itemize{
#'          \item All standard annotation columns from \code{annotatePeakInBatch}
#'                or \code{annotatePeaksNearBDP}
#'          \item \code{priority_score}: Numeric priority score (higher = better)
#'                for prioritized annotations
#'          \item \code{strategy_name}: Name of annotation strategy that produced
#'                this annotation (e.g., "promoter_focused", "upstream_to_downstream",
#'                "bidirectional_promoters")
#'        }
#' 
#' @details
#' 
#' \strong{Algorithm:}
#' 
#' \enumerate{
#'   \item \strong{Factor Type Validation}: Validates that \code{factor_type} is
#'         specified (required for accurate prioritization).
#'   \item \strong{Parallel Annotation}: Applies multiple annotation strategies
#'         simultaneously using transcripts extracted from the \code{EnsDb} object:
#'         \itemize{
#'           \item Promoter-focused: Uses \code{annotatePeakInBatch} with
#'                 \code{output = "both"}, \code{maxgap = 3000},
#'                 \code{FeatureLocForDistance = "TSS"}, \code{PeakLocForDistance = "middle"}
#'           \item Upstream to downstream: Uses \code{annotatePeakInBatch} with
#'                 \code{output = "upstream2downstream"}, \code{maxgap = 5000},
#'                 \code{FeatureLocForDistance = "TSS"}, \code{PeakLocForDistance = "middle"}
#'           \item Bidirectional promoters: Uses \code{annotatePeaksNearBDP} with
#'                 \code{maxTSSDistance = 1000}
#'         }
#'   \item \strong{Special Case Handling}:
#'         \itemize{
#'           \item \strong{Bidirectional Promoters}: For promoter-binding factors
#'                 (TF, H3K4me3, H3K27ac, H3K4me2), annotations from the
#'                 "bidirectional_promoters" strategy are kept without ranking.
#'           \item \strong{Broad Peaks}: If a peak completely covers multiple
#'                 features (\code{insideFeature == "includeFeature"}), all
#'                 completely covered features are kept without ranking.
#'         }
#'   \item \strong{Factor-Specific Prioritization}: For annotations not in special
#'         cases, applies factor-specific prioritization function to calculate
#'         priority scores.
#'   \item \strong{Result Selection}: If \code{keep = "best"}, returns annotation 
#'         with highest priority score per peak (plus special cases). 
#'         If \code{keep = "all"}, returns all annotations with priority scores.
#' }
#' 
#' \strong{When to Use Hierarchical Annotation:}
#' 
#' - When you need comprehensive annotation from multiple perspectives
#' - When factor-specific prioritization is important
#' - When dealing with complex binding patterns (bidirectional promoters, broad
#'   domains)
#' 
#' \strong{When to Use Single-Method Annotation:}
#' 
#' - When you have a specific annotation question (e.g., "find peaks near TSS")
#' - When computational speed is critical
#' - When you want full control over annotation parameters
#' 
#' @author Haibo Liu
#' @seealso \code{\link{annotatePeakInBatch}} for single-method annotation,
#'          \code{\link{annoPeaks}} for region-based annotation,
#'          \code{\link{annotatePeaksNearBDP}} for bidirectional promoter detection,
#'          \code{\link[ensembldb]{EnsDb-class}} for EnsDb objects,
#'          \code{\link{setChIPpeakAnnoEnsDb}} for setting global EnsDb option
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics start end
#' @importFrom ensembldb transcripts
#' @export
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v75)
#' library(GenomeInfoDb)
#' data(myPeakList)
#' 
#' # Create EnsDb object (or use pre-loaded EnsDb package)
#' EnsDb <- EnsDb.Hsapiens.v75
#' seqlevelsStyle(myPeakList) <- seqlevelsStyle(EnsDb)[1]
#' 
#' # Specify factor type explicitly (required)
#' anno_tf <- annotateHierarchically(
#'     myPeakList[1:100],
#'     EnsDb = EnsDb,
#'     factor_type = "TF",
#'     is_promoter_binding = TRUE
#' )
#' 
#' # Return all annotations with priority scores
#' anno_all <- annotateHierarchically(
#'     myPeakList[1:100],
#'     EnsDb = EnsDb,
#'     factor_type = "H3K27me3",
#'     is_promoter_binding = FALSE,
#'     keep = "all"
#' )
#' 
#' # Example 3: Using a custom prioritization function
#' customPrioritize <- function(annotated_peaks, weight_distance = 1000) {
#'     # Simple custom prioritization: prioritize by distance only
#'     if ("distanceToSite" %in% colnames(mcols(annotated_peaks))) {
#'         distance <- abs(mcols(annotated_peaks)$distanceToSite)
#'         annotated_peaks$priority_score <- weight_distance * exp(-distance / 500)
#'     } else {
#'         annotated_peaks$priority_score <- 0
#'     }
#'     return(annotated_peaks)
#' }
#' 
#' anno_custom <- annotateHierarchically(
#'     myPeakList[1:100],
#'     EnsDb = EnsDb,
#'     is_promoter_binding = FALSE,
#'     factor_type = "H3K27me3",
#'     prioritization_function = customPrioritize,
#'     weight_distance = 2000  # Pass additional parameters via ...
#' )
#' }
annotateHierarchically <- function(peaks,
                                   EnsDb = NULL,
                                   factor_type = c("TF", "H3K4me3", "H3K4me2", 
                                                   "H3K27ac", "H3K36me3", "H3K27me3", 
                                                   "PolII", "RBP", "ATAC", "3end",
                                                   "exon", "intron", "architectural",
                                                   "ncRNA", "repeat", "intergenic", 
                                                   "chromatin_remodeler", "multimodal"),
                                    is_promoter_binding = TRUE,
                                    annotation_strategies = NULL,
                                    prioritization_weights = NULL,
                                    prioritization_function = NULL,
                                    keep = c("best", "all"),
                                    BPPARAM = NULL,
                                    sequencing_method = c("ChIP-seq", "GRO-seq", "PRO-seq",
                                                          "CLIP-seq", "iCLIP", "eCLIP"),
                                    ...
                                ) {
    # Input validation
    if (missing(peaks)) {
        stop("Missing required argument 'peaks'!", call. = FALSE) 
    }
    if (!inherits(peaks, "GRanges")) {
        stop("'peaks' must be a GRanges object", call. = FALSE)
    }

    # Check for EnsDb: use provided, then global option, then error
    if (missing(EnsDb) || is.null(EnsDb)) {
        EnsDb <- getChIPpeakAnnoEnsDb()
        if (is.null(EnsDb)) {
            stop("Missing required argument 'EnsDb'! ",
                 "Either provide EnsDb explicitly or set it globally using ",
                 "setChIPpeakAnnoEnsDb().", call. = FALSE)
        }
        message("Using global EnsDb for annotation")
    }
    if (!inherits(EnsDb, "EnsDb")) {
        stop("'EnsDb' must be an EnsDb object", call. = FALSE)
    }
    
    # Ensure peaks have names
    if (is.null(names(peaks))) {
        names(peaks) <- paste0("peak_", seq_along(peaks))
    }
    
    # Step 1: Parameter validation
    keep <- match.arg(keep)
    
    # Factor type validation (only required if using default prioritization)
    if (is.null(prioritization_function)) {
        if (missing(factor_type) || is.null(factor_type)) {
            stop("'factor_type' must be specified when using default prioritization. ",
                 "Please specify the biological factor type ",
                 "(e.g., 'TF', 'H3K4me3', 'H3K36me3', 'H3K27me3', 'PolII', 'RBP', etc.). ",
                 "Alternatively, provide a custom 'prioritization_function'. ",
                 "See ?annotateHierarchically for available factor types.",
                 call. = FALSE)
        }
    }

    if (is.null(BPPARAM)) {
        if (requireNamespace("BiocParallel", quietly = TRUE)) {
            BPPARAM <- BiocParallel::bpparam()
        } else {
            BPPARAM <- NULL  # Will use sequential processing below
        }
    }

    # Step 2: Apply multiple annotation strategies in parallel
    message("Applying multiple annotation strategies...")
    strategy_results <- applyMultipleStrategies(peaks, EnsDb, 
                                                 annotation_strategies, 
                                                 BPPARAM = BPPARAM)
    
    # Step 3: Combine all annotations
    # Extract non-empty annotations more efficiently
    all_annotations <- lapply(strategy_results, function(x) x$annotations)
    all_annotations <- all_annotations[lengths(all_annotations) > 0L]
    
    if (length(all_annotations) == 0L) {
        warning("No annotations found from any strategy.", call. = FALSE)
        return(GRanges())
    }
    
    # Combine all annotations at once
    combined_anno <- do.call(c, all_annotations)
    
    # Ensure peak names are present and preserve original peak names
    # Annotation functions (annotatePeakInBatch, annoPeaks) preserve original
    # peak names in the "peak" column, but we verify consistency
    if (!"peak" %in% colnames(mcols(combined_anno))) {
        # Peak column missing - this should not happen
        stop("Peak names not found in annotation results. ",
             "Annotation functions should preserve peak names in 'peak' column.",
             call. = FALSE)
    }
    
    # Verify that all peak names in annotations match original peak names
    # This ensures we're using the original peak names from the input
    unique_anno_peak_names <- unique(combined_anno$peak)
    if (!all(unique_anno_peak_names %in% unique(names(peaks)))) {
        stop("Some peak names in annotations don't match original peak names. ",
               "This may indicate an issue with annotation functions. ",
               "Original peak names should be preserved.",
               call. = FALSE)
    }
    
    # Step 4: Identify and separate special cases (hierarchical annotation)
    # Peaks with bidirectional promoters or completely covered features should
    # NOT be considered for priority score-based annotation
    
    # Special Case 1: Bidirectional promoters for promoter-binding factors
    bdp_annotations <- GRanges()
    bdp_peaks <- character(0)
    if (is_promoter_binding) {
        bdp_mask <- mcols(combined_anno)$strategy_name == "bidirectional_promoters"
        bdp_annotations <- combined_anno[bdp_mask]
        if (length(bdp_annotations) > 0L) {
            bdp_peaks <- unique(bdp_annotations$peak)
            message("Found ", length(bdp_peaks), 
                   " peaks in bidirectional promoters. ",
                   "Excluding from priority score-based annotation.")
        }
    }
    
    # Special Case 2: Broad peaks completely covering multiple features
    # Exclude peaks already annotated by bidirectional promoters (special case 1)
    completely_covered_annotations <- GRanges()
    covered_peaks <- character(0)
    if ("insideFeature" %in% colnames(mcols(combined_anno))) {
        # Exclude bidirectional promoter peaks from completely covered analysis
        covered_mask <- (mcols(combined_anno)$insideFeature == "includeFeature" &
                        !combined_anno$peak %in% bdp_peaks)
     
        completely_covered_annotations <- combined_anno[covered_mask]
        if (length(completely_covered_annotations) > 0L) {
            covered_peaks <- unique(completely_covered_annotations$peak)
            message("Found ", length(covered_peaks), 
                   " peaks completely covering features",
                   if (length(bdp_peaks) > 0L) {
                       paste0(" (excluding ", length(bdp_peaks), 
                             " already in bidirectional promoters)")
                   } else {
                       ""
                   }, ". ",
                   "Excluding from priority score-based annotation.")
        }
    }
    
    # Identified peaks that should be excluded from prioritization: bidirectional 
    # promoters and completely covered peaks
    excluded_peaks <- unique(c(bdp_peaks, covered_peaks))
    
    # Step 5: Apply factor-specific prioritization only to non-excluded peaks
    # Filter out excluded peaks for prioritization
    prioritization_mask <- !combined_anno$peak %in% excluded_peaks
    anno_for_prioritization <- combined_anno[prioritization_mask]
    
    # Apply prioritization only to peaks not in special cases
    if (length(anno_for_prioritization) > 0L) {
        if (!is.null(prioritization_function)) {
            # Use user-provided prioritization function
            message("Applying custom prioritization function (", 
                   length(unique(anno_for_prioritization$peak)), " peaks)")
            anno_prioritized <- prioritization_function(
                anno_for_prioritization, 
                factor_type = factor_type,
                sequencing_method = sequencing_method,
                ...
            )
            # Validate that priority_score column exists
            if (!"priority_score" %in% colnames(mcols(anno_prioritized))) {
                stop("Custom prioritization function must return a GRanges object ",
                     "with a 'priority_score' column in metadata.", call. = FALSE)
            }
        } else {
            # Use default factor-specific prioritization
            message("Applying factor-specific prioritization for: ", factor_type,
                   " (", length(unique(anno_for_prioritization$peak)), " peaks)")
            anno_prioritized <- prioritizeAnnotations(
                anno_for_prioritization, 
                factor_type = factor_type,
                sequencing_method = sequencing_method,
                ...
            )
        }
    } else {
        message("No peaks remaining for prioritization (all in special cases).")
    }
    
    # Step 6: Result selection and combination
    result_parts <- list()
    if (keep == "all") {
        # Return all annotations: special cases + prioritized annotations
        
        if (length(bdp_annotations) > 0L) {
            result_parts[["bidirectional_promoters"]] <- bdp_annotations
        }
        if (length(completely_covered_annotations) > 0L) {
            result_parts[["completely_covered"]] <- completely_covered_annotations
        }
        if (length(anno_prioritized) > 0L) {
            result_parts[["prioritized"]] <- anno_prioritized
        }
    } else {
        # Return best annotation per peak for prioritized annotations,
        # plus all special case annotations
        # Add bidirectional promoter annotations (all of them)
        if (length(bdp_annotations) > 0L) {
            result_parts[["bidirectional_promoters"]] <- bdp_annotations
        }
        
        # Add completely covered annotations (all of them)
        if (length(completely_covered_annotations) > 0L) {
            result_parts[["completely_covered"]] <- completely_covered_annotations
        }
        
        # Get best annotation per peak for prioritized annotations
        if (length(anno_for_prioritization) > 0L) {
            # Sort by peak name, then by priority score (descending)
            anno_sorted <- anno_prioritized[order(
                anno_prioritized$peak, 
                -anno_prioritized$priority_score, 
                na.last = TRUE
            )]
            
            # Get first (best) annotation for each peak
            best_annotations <- anno_sorted[!duplicated(anno_sorted$peak)]
            result_parts[["prioritized"]] <- best_annotations
        }
    }
        if (length(result_parts) == 0L) {
        warning("No annotations found from any strategy.", call. = FALSE)
            return(GRanges())
        }
        
    # Combine all annotations
        result <- do.call(c, result_parts)
        return(result)
}


#' Apply multiple annotation strategies in parallel
#' 
#' @description 
#' Internal helper function to apply multiple annotation strategies
#' simultaneously to the same set of peaks. Each strategy uses different
#' parameters optimized for different genomic contexts (promoters, gene bodies,
#' intergenic regions, bidirectional promoters).
#' 
#' This function uses \code{\link[BiocParallel]{bplapply}} for parallel
#' execution, allowing multiple strategies to run simultaneously and
#' potentially providing significant speedup (up to N-fold where N is the
#' number of strategies, limited by available CPU cores). See Details for
#' setup instructions.
#' 
#' @param peaks A \code{GRanges} object containing peaks to be annotated.
#' @param EnsDb A \code{EnsDb} object containing
#'   genome annotation data (genes, transcripts, exons, etc.) for an organism.
#' @param strategies A list of strategy definitions. Each strategy is a list
#'   with:
#'   \itemize{
#'     \item \code{name}: Character string identifying the strategy
#'     \item \code{method}: "annotatePeakInBatch" or "annoPeaks"
#'     \item \code{params}: List of parameters to pass to the annotation method
#'   }
#'   If \code{NULL}, default strategies will be used.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'        object specifying the parallel backend to use. If \code{NULL}
#'        (default), uses \code{\link[BiocParallel]{bpparam}()} which
#'        automatically detects the best available backend. For sequential
#'        processing, use \code{BiocParallel::SerialParam()}.
#' @return A list of annotation results, one per strategy. Each result is a
#'   list with:
#'   \itemize{
#'     \item \code{annotations}: \code{GRanges} object with annotations
#'     \item \code{metadata}: Additional information (overlap counts, etc.)
#'   }
#' @details
#' 
#' \strong{Parallel Execution:}
#' 
#' This function uses \code{\link[BiocParallel]{bplapply}} for parallel
#' execution of annotation strategies. Parallel processing is automatically
#' enabled if \code{BiocParallel} is available. To control parallel processing:
#' \itemize{
#'   \item Install the \code{BiocParallel} package (included in Bioconductor)
#'   \item Use \code{BPPARAM} parameter to specify the parallel backend:
#'         \code{BiocParallel::MulticoreParam(workers = 4)} for local
#'         multicore, or \code{BiocParallel::SnowParam(workers = 4)} for
#'         cluster computing
#'   \item If \code{BPPARAM = NULL}, the function automatically uses
#'         \code{BiocParallel::bpparam()} which detects the best available
#'         backend
#' }
#' 
#' \strong{Performance and Speedup:}
#' \itemize{
#'   \item \strong{Parallel execution}: Multiple strategies run simultaneously,
#'         potentially reducing total computation time by up to N-fold where N
#'         is the number of strategies (limited by available CPU cores)
#'   \item \strong{Expected speedup}: With 3 default strategies and sufficient
#'         cores, expect 2-3x speedup for large peak sets (>1000 peaks)
#'   \item \strong{Best performance when}:
#'         \itemize{
#'           \item Multiple strategies are applied (default: 3 strategies)
#'           \item Large peak sets are being annotated (>1000 peaks)
#'           \item Annotation operations are computationally intensive
#'           \item Multiple CPU cores are available
#'         }
#'   \item \strong{Overhead considerations}: For small peak sets (< 100 peaks),
#'         sequential execution may be faster due to parallelization overhead.
#'         The function automatically uses parallel execution if
#'         \code{BiocParallel} is available.
#' }
#' 
#' \strong{Default Strategies:}
#' 
#' The default strategies are optimized for transcription factor (TF) binding
#' factors. Transcripts are extracted from the \code{EnsDb} object using
#' \code{ensembldb::transcripts()}. The default strategies include:
#' 
#' \enumerate{
#'   \item \strong{Promoter-focused}: Uses \code{annotatePeakInBatch} with
#'         \code{output = "both"}, \code{maxgap = 3000},
#'         \code{FeatureLocForDistance = "TSS"}, \code{PeakLocForDistance = "middle"}.
#'         Optimized for finding promoter-proximal binding and nearest features.
#'   \item \strong{Upstream to downstream}: Uses \code{annotatePeakInBatch} with
#'         \code{output = "upstream2downstream"}, \code{maxgap = 5000},
#'         \code{FeatureLocForDistance = "TSS"}, \code{PeakLocForDistance = "middle"}.
#'         Captures peaks upstream of TSS, inside gene bodies, and downstream of gene end
#'         (within maxgap distance).
#'   \item \strong{Bidirectional promoters}: Uses \code{annotatePeaksNearBDP} with
#'         \code{maxTSSDistance = 1000}. Detects peaks near bidirectional promoters.
#' }
#' 
#' Note: Strategy parameters (including \code{myPeakList} and \code{AnnotationData})
#' must be included in the \code{params} list when defining custom strategies.
#' 
#' @author Haibo Liu
#' @export
#' @importFrom BiocGenerics start end
#' @importFrom BiocParallel bplapply bpparam SerialParam
applyMultipleStrategies <- function(peaks, EnsDb, 
                                    strategies = NULL, 
                                    BPPARAM = NULL) {
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(EnsDb, "EnsDb"))
    
    # Default strategies for TF binding factors if none provided
    if (is.null(strategies)) {
        transcripts <- ensembldb::transcripts(EnsDb)
        strategies <- list(
            list(
                name = "promoter_focused",
                method = "annotatePeakInBatch",
                params = list(
                    myPeakList = peaks,
                    AnnotationData = transcripts,
                    output = "both",
                    maxgap = 3000L,
                    PeakLocForDistance = "middle",
                    FeatureLocForDistance = "TSS",
                    select = "all",
                    ignore.strand = TRUE
                )
            ),
            list(
                name = "upstream_to_downstream",
                method = "annotatePeakInBatch",
                params = list(
                    myPeakList = peaks,
                    AnnotationData = transcripts,
                    output = "upstream2downstream",
                    maxgap = 5000L,
                    FeatureLocForDistance = "TSS",
                    PeakLocForDistance = "middle",
                    select = "all",
                    ignore.strand = TRUE
                )
            ),
            list(
                name = "bidirectional_promoters",
                method = "annotatePeaksNearBDP",
                params = list(
                    peaks = peaks,
                    annoData = transcripts,
                    maxTSSDistance = 1000L,
                    ignore.peak.strand = TRUE
                )
            )
        )
    }
    
    # Helper function to apply a single strategy
    applySingleStrategy <- function(strategy) {
        strategy_name <- strategy$name
        method <- strategy$method
        params <- strategy$params
        
        tryCatch({
            if (method == "annotatePeakInBatch") {
                # Call annotatePeakInBatch with params
                anno_result <- do.call(annotatePeakInBatch, params)
                
                # Add strategy metadata
                if (length(anno_result) > 0L) {
                    anno_result$strategy_name <- strategy_name
                }
                
                return(list(
                    annotations = anno_result,
                    metadata = list(
                        n_annotations = length(anno_result),
                        n_unique_peaks = if (length(anno_result) > 0L && "peak" %in% colnames(mcols(anno_result))) {
                            length(unique(anno_result$peak))
                        } else {
                            0L
                        }
                    )
                ))
            } else if (method == "annotatePeaksNearBDP"){
                # Call annoPeaks with params
                anno_result <- do.call(annotatePeaksNearBDP, params)
                
                # Add strategy metadata
                if (length(anno_result) > 0L) {
                    anno_result$strategy_name <- strategy_name
                }
                
                return(list(
                    annotations = anno_result,
                    metadata = list(
                        n_annotations = length(anno_result),
                        n_unique_peaks = if (length(anno_result) > 0L && "peak" %in% colnames(mcols(anno_result))) {
                            length(unique(anno_result$peak))
                        } else {
                            0L
                        }
                    )
                ))
            } else if (method == "annoPeaks") {
                anno_result <- do.call(annoPeaks, params)
                
                # Add strategy metadata
                if (length(anno_result) > 0L) {
                    anno_result$strategy_name <- strategy_name
                }
            
                return(list(
                    annotations = anno_result,
                    metadata = list(
                        n_annotations = length(anno_result),
                        n_unique_peaks = if (length(anno_result) > 0L && "peak" %in% colnames(mcols(anno_result))) {
                            length(unique(anno_result$peak))
                        } else {
                            0L
                        }
                    )
                ))
            } else {
                stop("Unknown annotation method: ", method, 
                       " for strategy: ", strategy_name, call. = FALSE)
            }}, error = function(e) {
                stop("Error applying strategy '", strategy_name, "': ", 
                   conditionMessage(e), call. = FALSE)
            }
        )
    }

    # Setup BiocParallel backend
    if (is.null(BPPARAM)) {
        if (requireNamespace("BiocParallel", quietly = TRUE)) {
            BPPARAM <- BiocParallel::bpparam()
        } else {
            BPPARAM <- NULL  # Will use sequential processing below
        }
    }
    
    # Apply strategies in parallel using BiocParallel
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
        results_list <- BiocParallel::bplapply(strategies, applySingleStrategy,
                                              BPPARAM = BPPARAM)
    } else {
        # Fallback to sequential if BiocParallel not available
        results_list <- lapply(strategies, applySingleStrategy)
    }
    
    # Convert to named list using strategy names
    results <- setNames(results_list, sapply(strategies, function(s) s$name))
    
    return(results)
}
