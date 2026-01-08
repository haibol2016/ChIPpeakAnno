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
#' @param annoData A \code{GRanges} or \code{annoGR} object containing
#'   annotation data.
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
#'     \item \code{strategy_name}: Name of the strategy
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
#'   \item \strong{Expected speedup}: With 4-5 default strategies and sufficient
#'         cores, expect 2-4x speedup for large peak sets (>1000 peaks)
#'   \item \strong{Best performance when}:
#'         \itemize{
#'           \item Multiple strategies are applied (default: 5 strategies)
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
#' \enumerate{
#'   \item \strong{Promoter-focused}: Uses \code{annotatePeakInBatch} with
#'         \code{output = "overlapping"}, \code{maxgap = 2000},
#'         \code{FeatureLocForDistance = "TSS"}. Optimized for finding
#'         promoter-proximal binding.
#'   \item \strong{Gene body}: Uses \code{annoPeaks} with
#'         \code{bindingType = "fullRange"}, \code{bindingRegion = c(-1000, 1000)}.
#'         Optimized for finding gene body binding.
#'   \item \strong{Nearest (comprehensive)}: Uses \code{annotatePeakInBatch}
#'         with \code{output = "nearestLocation"}, \code{maxgap = 10000}.
#'         Finds nearest features regardless of relationship type.
#'   \item \strong{Bidirectional promoters}: Uses fixed bidirectional promoter
#'         detection (from \code{.identifyBidirectionalPromoters}).
#'   \item \strong{Intergenic/distal}: Uses \code{annotatePeakInBatch} with
#'         \code{output = "nearestLocation"}, \code{maxgap = 50000}. Finds
#'         distal/intergenic features.
#' }
#' 
#' @author Haibo Liu
#' @keywords internal
#' @importFrom BiocGenerics start end
#' @importFrom BiocParallel bplapply bpparam SerialParam
.applyMultipleStrategies <- function(peaks, annoData, strategies = NULL, BPPARAM = NULL) {
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("GRanges", "annoGR")))
    
    # Default strategies if none provided
    if (is.null(strategies)) {
        strategies <- list(
            list(
                name = "promoter",
                method = "annotatePeakInBatch",
                params = list(
                    output = "overlapping",
                    maxgap = 2000L,
                    FeatureLocForDistance = "TSS",
                    PeakLocForDistance = "middle",
                    select = "all"
                )
            ),
            list(
                name = "gene_body",
                method = "annoPeaks",
                params = list(
                    bindingType = "fullRange",
                    bindingRegion = c(-1000L, 1000L),
                    select = "all"
                )
            ),
            list(
                name = "nearest",
                method = "annotatePeakInBatch",
                params = list(
                    output = "nearestLocation",
                    maxgap = 10000L,
                    FeatureLocForDistance = "TSS",
                    PeakLocForDistance = "middle",
                    select = "all"
                )
            ),
            list(
                name = "bidirectional",
                method = "bidirectional_promoters",
                params = list(
                    maxTSSDistance = 1000L,
                    maxPeakToBDPDistance = 200L
                )
            ),
            list(
                name = "intergenic",
                method = "annotatePeakInBatch",
                params = list(
                    output = "nearestLocation",
                    maxgap = 50000L,
                    FeatureLocForDistance = "TSS",
                    PeakLocForDistance = "middle",
                    select = "all"
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
                anno_result <- do.call(
                    annotatePeakInBatch,
                    c(list(myPeakList = peaks, AnnotationData = annoData), params)
                )
                
                # Add strategy metadata
                if (length(anno_result) > 0L) {
                    anno_result$strategy_name <- strategy_name
                }
                
                return(list(
                    annotations = anno_result,
                    strategy_name = strategy_name,
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
                # Call annoPeaks with params
                anno_result <- do.call(
                    annoPeaks,
                    c(list(peaks = peaks, annoData = annoData), params)
                )
                
                # Add strategy metadata
                if (length(anno_result) > 0L) {
                    anno_result$strategy_name <- strategy_name
                }
                
                return(list(
                    annotations = anno_result,
                    strategy_name = strategy_name,
                    metadata = list(
                        n_annotations = length(anno_result),
                        n_unique_peaks = if (length(anno_result) > 0L && "peak" %in% colnames(mcols(anno_result))) {
                            length(unique(anno_result$peak))
                        } else {
                            0L
                        }
                    )
                ))
            } else if (method == "bidirectional_promoters") {
                # Use bidirectional promoter detection
                maxTSSDist <- if (!is.null(params$maxTSSDistance)) params$maxTSSDistance else 1000L
                bdp_regions <- .identifyBidirectionalPromoters(
                    annoData, 
                    maxTSSDistance = maxTSSDist
                )
                
                if (length(bdp_regions) > 0L) {
                    # Calculate distance from peak center to BDP center
                    peak_centers <- as.integer(round((start(peaks) + end(peaks)) / 2))
                    bdp_centers <- mcols(bdp_regions)$bdp_center
                    maxPeakToBDPDist <- if (!is.null(params$maxPeakToBDPDistance)) params$maxPeakToBDPDistance else 200L
                    
                    # Find peaks near bidirectional promoters
                    peak_bdp_pairs <- list()
                    for (j in seq_along(peaks)) {
                        peak_chr <- seqnames(peaks)[j]
                        peak_center <- peak_centers[j]
                        
                        bdp_chr_idx <- seqnames(bdp_regions) == peak_chr
                        if (sum(bdp_chr_idx) == 0L) {
                            next
                        }
                        
                        bdp_chr_centers <- bdp_centers[bdp_chr_idx]
                        distances <- abs(bdp_chr_centers - peak_center)
                        near_bdp_idx <- distances <= maxPeakToBDPDist
                        
                        if (any(near_bdp_idx)) {
                            bdp_global_idx <- which(bdp_chr_idx)
                            for (k in which(near_bdp_idx)) {
                                peak_bdp_pairs[[length(peak_bdp_pairs) + 1L]] <- list(
                                    peak_idx = j,
                                    bdp_idx = bdp_global_idx[k],
                                    distance = distances[k]
                                )
                            }
                        }
                    }
                    
                    # Create annotations for both genes in each BDP
                    if (length(peak_bdp_pairs) > 0L) {
                        all_annotations <- list()
                        for (pair in peak_bdp_pairs) {
                            bdp <- bdp_regions[pair$bdp_idx]
                            gene1_idx <- mcols(bdp)$gene1_idx
                            gene2_idx <- mcols(bdp)$gene2_idx
                            
                            if (gene1_idx > 0L && gene1_idx <= length(annoData) &&
                                gene2_idx > 0L && gene2_idx <= length(annoData)) {
                                gene1_anno <- annoData[gene1_idx]
                                gene2_anno <- annoData[gene2_idx]
                                peak_gr <- peaks[pair$peak_idx]
                                
                                # Create annotation entries for both genes
                                # (simplified - full implementation would include all metadata)
                                anno1 <- peak_gr
                                anno1$feature <- if (!is.null(names(gene1_anno))) names(gene1_anno) else mcols(bdp)$gene1_id
                                anno1$strategy_name <- strategy_name
                                anno1$is_bidirectional_promoter <- TRUE
                                
                                anno2 <- peak_gr
                                anno2$feature <- if (!is.null(names(gene2_anno))) names(gene2_anno) else mcols(bdp)$gene2_id
                                anno2$strategy_name <- strategy_name
                                anno2$is_bidirectional_promoter <- TRUE
                                
                                all_annotations[[length(all_annotations) + 1L]] <- anno1
                                all_annotations[[length(all_annotations) + 1L]] <- anno2
                            }
                        }
                        
                        if (length(all_annotations) > 0L) {
                            anno_result <- do.call(c, all_annotations)
                            return(list(
                                annotations = anno_result,
                                strategy_name = strategy_name,
                                metadata = list(
                                    n_annotations = length(anno_result),
                                    n_unique_peaks = if ("peak" %in% colnames(mcols(anno_result))) {
                                        length(unique(anno_result$peak))
                                    } else {
                                        0L
                                    },
                                    n_bdp_regions = length(unique(sapply(peak_bdp_pairs, function(p) p$bdp_idx)))
                                )
                            ))
                        } else {
                            return(list(
                                annotations = GRanges(),
                                strategy_name = strategy_name,
                                metadata = list(n_annotations = 0L, n_unique_peaks = 0L)
                            ))
                        }
                    } else {
                        return(list(
                            annotations = GRanges(),
                            strategy_name = strategy_name,
                            metadata = list(n_annotations = 0L, n_unique_peaks = 0L)
                        ))
                    }
                } else {
                    return(list(
                        annotations = GRanges(),
                        strategy_name = strategy_name,
                        metadata = list(n_annotations = 0L, n_unique_peaks = 0L)
                    ))
                }
            } else {
                warning("Unknown annotation method: ", method, 
                       " for strategy: ", strategy_name, call. = FALSE)
                return(list(
                    annotations = GRanges(),
                    strategy_name = strategy_name,
                    metadata = list(n_annotations = 0L, n_unique_peaks = 0L)
                ))
            }
        }, error = function(e) {
            warning("Error applying strategy '", strategy_name, "': ", 
                   conditionMessage(e), call. = FALSE)
            return(list(
                annotations = GRanges(),
                strategy_name = strategy_name,
                metadata = list(n_annotations = 0L, n_unique_peaks = 0L, error = conditionMessage(e))
            ))
        })
    }

    # Setup BiocParallel backend
    if (is.null(BPPARAM)) {
        if (requireNamespace("BiocParallel", quietly = TRUE)) {
            BPPARAM <- BiocParallel::bpparam()
        } else {
            BPPARAM <- BiocParallel::SerialParam()
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

