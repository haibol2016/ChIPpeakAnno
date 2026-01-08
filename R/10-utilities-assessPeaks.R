#' Assess peak widths and scores
#' 
#' @description 
#' Provides comprehensive assessment of peak widths and scores (if present) in
#' a GRanges object or a list of GRanges objects (e.g., from multiple replicates).
#' This function is particularly useful for quality control after loading peaks
#' using \code{\link{toGRanges}}, as it summarizes the distribution of peak
#' widths and any score-related metadata columns.
#' 
#' The function automatically detects score-related columns in the metadata
#' (e.g., "score", "signalValue", "pValue", "qValue", "fold_enrichment",
#' "-log10(pvalue)", "-log10(qvalue)") and provides summary statistics for
#' each. When multiple peak sets are provided, the function also generates
#' comparative statistics to facilitate replicate comparison.
#' 
#' @param peaks A \link[GenomicRanges:GRanges-class]{GRanges} object or a list
#'        of GRanges objects containing peaks to be assessed. Typically loaded
#'        using \code{\link{toGRanges}} from BED, narrowPeak, broadPeak, MACS,
#'        or MACS2 output files. If a list is provided, each element should be
#'        a GRanges object representing one replicate or sample. The list can be
#'        named (names will be used for labeling) or unnamed (names will be
#'        auto-generated as "Replicate1", "Replicate2", etc.).
#' @param scoreColumns A character vector specifying which metadata columns to
#'        assess as scores. If \code{NULL} (default), the function automatically
#'        detects common score-related columns: "score", "signalValue", "pValue",
#'        "qValue", "fold_enrichment", "-log10(pvalue)", "-log10(qvalue)",
#'        "pileup", "tags", "-10*log10(pvalue)", "FDR". Users can specify
#'        custom column names if needed. When multiple peak sets are provided,
#'        only columns present in all sets will be assessed for comparison.
#' @param verbose Logical. If \code{TRUE} (default), prints a brief summary to
#'        the console. If \code{FALSE}, only returns the results without
#'        printing.
#' 
#' @return Returns a list with the following components:
#'        \itemize{
#'          \item \code{width_plot}: A ggplot object visualizing peak width
#'                distribution. For single GRanges: histogram or density plot.
#'                For multiple GRanges: boxplot or violin plot comparing widths
#'                across replicates.
#'          \item \code{score_plots}: A list of ggplot objects, one for each
#'                detected score column, visualizing score distributions. For
#'                single GRanges: histogram or density plot. For multiple
#'                GRanges: boxplot or violin plot comparing scores across
#'                replicates.
#'          \item \code{width_data}: A data frame containing the raw width data
#'                used for plotting (includes replicate names for multi-replicate
#'                analysis)
#'          \item \code{score_data}: A list of data frames, one for each score
#'                column, containing the raw score data used for plotting
#'          \item \code{detected_score_columns}: A character vector of score
#'                column names that were found and assessed (common across all
#'                replicates when multiple sets are provided)
#'          \item \code{all_metadata_columns}: A list of character vectors, one
#'                for each replicate, containing all metadata column names
#'          \item \code{is_multi_replicate}: Logical indicating whether multiple
#'                peak sets were assessed
#'          \item \code{replicate_names}: Character vector of replicate names
#'        }
#' 
#' @details
#' 
#' \strong{Peak Width Assessment:}
#' 
#' Peak widths are calculated as \code{width(peaks)} (which equals
#' \code{end - start + 1}). The function provides standard summary statistics
#' (min, quartiles, median, mean, max) to help users understand the size
#' distribution of their peaks. This is particularly useful for:
#' - Identifying unusually narrow or broad peaks
#' - Understanding peak size distribution (e.g., narrow peaks for TFs vs broad
#'   peaks for histone marks)
#' - Quality control (detecting potential issues with peak calling)
#' 
#' \strong{Score Assessment:}
#' 
#' The function automatically searches for common score-related columns in the
#' metadata. For each detected score column, it provides:
#' - Standard summary statistics (min, quartiles, median, mean, max)
#' - Counts of NA, zero, negative, and positive values
#' - This helps users understand the distribution and quality of scores
#' 
#' Common score columns from different peak calling tools:
#' - \strong{BED}: "score" (integer 0-1000)
#' - \strong{narrowPeak/broadPeak}: "score", "signalValue", "pValue", "qValue"
#' - \strong{MACS}: "tags", "-10*log10(pvalue)", "fold_enrichment", "FDR"
#' - \strong{MACS2}: "pileup", "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)"
#' 
#' @author Haibo Liu
#' @seealso \code{\link{toGRanges}} for loading peaks from various formats,
#'          \code{\link[GenomicRanges]{width}} for calculating peak widths
#' @importFrom S4Vectors mcols
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin geom_histogram
#' @importFrom ggplot2 geom_density theme_bw theme_classic xlab ylab
#' @importFrom ggplot2 scale_fill_manual scale_color_manual facet_wrap
#' @importFrom ggplot2 element_text element_blank labs theme .data
#' @export
#' @examples
#' \dontrun{
#' # Single peak set
#' peaks <- toGRanges("peaks.bed", format = "BED")
#' assessment <- assessPeaks(peaks)
#' 
#' # Multiple replicates (as a list)
#' peaks_rep1 <- toGRanges("rep1_peaks.bed", format = "BED")
#' peaks_rep2 <- toGRanges("rep2_peaks.bed", format = "BED")
#' peaks_rep3 <- toGRanges("rep3_peaks.bed", format = "BED")
#' 
#' # Named list for better labeling
#' peak_list <- list(Rep1 = peaks_rep1, Rep2 = peaks_rep2, Rep3 = peaks_rep3)
#' multi_assessment <- assessPeaks(peak_list)
#' 
#' # Access plots
#' multi_assessment$width_plot
#' multi_assessment$score_plots$signalValue
#' 
#' # Customize plots
#' library(ggplot2)
#' multi_assessment$width_plot + labs(title = "Peak Width Comparison")
#' 
#' # Unnamed list (auto-generated names)
#' peak_list2 <- list(peaks_rep1, peaks_rep2, peaks_rep3)
#' multi_assessment2 <- assessPeaks(peak_list2)
#' 
#' # Silent mode (no console output)
#' assessment <- assessPeaks(peaks, verbose = FALSE)
#' }
assessPeaks <- function(peaks, scoreColumns = NULL, verbose = TRUE) {
    # Input validation
    if (missing(peaks)) {
        stop("Missing required argument 'peaks'!", call. = FALSE)
    }
    
    # Check if input is a list (multiple replicates) or single GRanges
    # A GRanges object inherits from many classes, so we check if it's
    # specifically a GRanges first, then check if it's a list
    if (inherits(peaks, "GRanges")) {
        # Single GRanges object
        return(.assessSinglePeaks(peaks, scoreColumns, verbose))
    } else if (is.list(peaks) && length(peaks) > 0L) {
        # List of peak sets (multiple replicates)
        return(.assessMultiplePeaks(peaks, scoreColumns, verbose))
    } else {
        stop("'peaks' must be a GRanges object or a list of GRanges objects", 
             call. = FALSE)
    }
}

# Internal function to assess a single GRanges object
.assessSinglePeaks <- function(peaks, scoreColumns = NULL, verbose = TRUE) {
    if (!inherits(peaks, "GRanges")) {
        stop("'peaks' must be a GRanges object or a list of GRanges objects", 
             call. = FALSE)
    }
    
    if (length(peaks) == 0L) {
        warning("Input GRanges object is empty.", call. = FALSE)
        return(list(
            width_plot = NULL,
            score_plots = list(),
            width_data = data.frame(),
            score_data = list(),
            detected_score_columns = character(0),
            all_metadata_columns = character(0),
            is_multi_replicate = FALSE,
            replicate_names = character(0)
        ))
    }
    
    # Get all metadata columns
    all_metadata_cols <- colnames(mcols(peaks))
    
    # Get peak widths
    peak_widths <- width(peaks)
    width_data <- data.frame(width = peak_widths)
    
    # Create width plot
    width_plot <- .plotWidths(width_data, is_multi = FALSE)
    
    # Detect score columns if not specified
    scoreColumns <- .detectScoreColumns(all_metadata_cols, scoreColumns)
    
    # Get score data and create plots
    score_data <- list()
    score_plots <- list()
    if (length(scoreColumns) > 0L) {
        for (col in scoreColumns) {
            score_values <- mcols(peaks)[[col]]
            
            # Convert to numeric if possible
            if (!is.numeric(score_values)) {
                score_values <- tryCatch(
                    as.numeric(score_values),
                    warning = function(w) {
                        warning("Column '", col, "' could not be converted to ",
                               "numeric. Skipping.", call. = FALSE)
                        return(NULL)
                    }
                )
                if (is.null(score_values)) {
                    next
                }
            }
            
            score_data[[col]] <- data.frame(score = score_values)
            score_plots[[col]] <- .plotScores(score_data[[col]], 
                                             col_name = col, 
                                             is_multi = FALSE)
        }
    }
    
    # Prepare results
    results <- list(
        width_plot = width_plot,
        score_plots = score_plots,
        width_data = width_data,
        score_data = score_data,
        detected_score_columns = scoreColumns,
        all_metadata_columns = all_metadata_cols,
        is_multi_replicate = FALSE,
        replicate_names = character(0)
    )
    
    # Print brief summary if verbose
    if (verbose) {
        cat("Assessed", length(peaks), "peaks.\n")
        cat("Detected", length(scoreColumns), "score column(s):", 
           paste(scoreColumns, collapse = ", "), "\n")
    }
    
    return(results)
}

# Internal function to assess multiple GRanges objects
.assessMultiplePeaks <- function(peaks_list, scoreColumns = NULL, verbose = TRUE) {
    # Validate list
    if (length(peaks_list) == 0L) {
        stop("Peak list is empty", call. = FALSE)
    }
    
    # Check all elements are GRanges
    for (i in seq_along(peaks_list)) {
        if (!inherits(peaks_list[[i]], "GRanges")) {
            stop("All elements in 'peaks' list must be GRanges objects. ",
                 "Element ", i, " is not a GRanges object.", call. = FALSE)
        }
    }
    
    # Get replicate names
    replicate_names <- names(peaks_list)
    if (is.null(replicate_names) || any(replicate_names == "")) {
        replicate_names <- paste0("Replicate", seq_along(peaks_list))
    }
    
    # Assess each replicate individually
    individual_results <- lapply(seq_along(peaks_list), function(i) {
        .assessSinglePeaks(peaks_list[[i]], scoreColumns = NULL, verbose = FALSE)
    })
    names(individual_results) <- replicate_names
    
    # Find common score columns across all replicates
    all_score_cols <- lapply(individual_results, function(x) {
        x$detected_score_columns
    })
    common_score_cols <- Reduce(intersect, all_score_cols)
    
    # If user specified scoreColumns, use intersection with common columns
    if (!is.null(scoreColumns)) {
        common_score_cols <- intersect(scoreColumns, common_score_cols)
        if (length(common_score_cols) == 0L) {
            warning("None of the specified score columns are present in all ",
                   "replicates.", call. = FALSE)
        }
    }
    
    # Re-assess with common score columns only
    if (length(common_score_cols) > 0L) {
        individual_results <- lapply(seq_along(peaks_list), function(i) {
            .assessSinglePeaks(peaks_list[[i]], 
                             scoreColumns = common_score_cols, 
                             verbose = FALSE)
        })
        names(individual_results) <- replicate_names
    }
    
    # Combine width data from all replicates
    width_data_list <- lapply(seq_along(individual_results), function(i) {
        data <- individual_results[[i]]$width_data
        data$replicate <- replicate_names[i]
        return(data)
    })
    width_data_combined <- do.call(rbind, width_data_list)
    width_data_combined$replicate <- factor(width_data_combined$replicate,
                                            levels = replicate_names)
    
    # Create comparative width plot
    width_plot <- .plotWidths(width_data_combined, is_multi = TRUE)
    
    # Combine score data from all replicates
    score_data_combined <- list()
    score_plots_combined <- list()
    if (length(common_score_cols) > 0L) {
        for (col in common_score_cols) {
            score_data_list <- lapply(seq_along(individual_results), function(i) {
                if (col %in% names(individual_results[[i]]$score_data)) {
                    data <- individual_results[[i]]$score_data[[col]]
                    data$replicate <- replicate_names[i]
                    return(data)
                } else {
                    return(NULL)
                }
            })
            score_data_list <- score_data_list[!sapply(score_data_list, is.null)]
            if (length(score_data_list) > 0L) {
                score_combined <- do.call(rbind, score_data_list)
                score_combined$replicate <- factor(score_combined$replicate,
                                                   levels = replicate_names)
                score_data_combined[[col]] <- score_combined
                score_plots_combined[[col]] <- .plotScores(score_combined,
                                                           col_name = col,
                                                           is_multi = TRUE)
            }
        }
    }
    
    # Get all metadata columns (as a list, one per replicate)
    all_metadata_cols_list <- lapply(individual_results, function(x) {
        x$all_metadata_columns
    })
    names(all_metadata_cols_list) <- replicate_names
    
    # Prepare results
    results <- list(
        width_plot = width_plot,
        score_plots = score_plots_combined,
        width_data = width_data_combined,
        score_data = score_data_combined,
        detected_score_columns = common_score_cols,
        all_metadata_columns = all_metadata_cols_list,
        is_multi_replicate = TRUE,
        replicate_names = replicate_names,
        individual_results = individual_results
    )
    
    # Print brief summary if verbose
    if (verbose) {
        cat("Assessed", length(peaks_list), "replicate(s):", 
           paste(replicate_names, collapse = ", "), "\n")
        total_peaks <- sum(sapply(peaks_list, length))
        cat("Total peaks across all replicates:", total_peaks, "\n")
        cat("Common score columns:", 
           if (length(common_score_cols) > 0L) {
               paste(common_score_cols, collapse = ", ")
           } else {
               "none"
           }, "\n")
    }
    
    return(results)
}

# Internal function to detect score columns
.detectScoreColumns <- function(all_metadata_cols, scoreColumns = NULL) {
    if (!is.null(scoreColumns)) {
        # Validate user-specified columns
        missing_cols <- scoreColumns[!scoreColumns %in% all_metadata_cols]
        if (length(missing_cols) > 0L) {
            warning("The following score columns were not found in metadata: ",
                   paste(missing_cols, collapse = ", "), call. = FALSE)
        }
        return(intersect(scoreColumns, all_metadata_cols))
    }
    
    # Common score-related column names from various formats
    common_score_names <- c(
        "score", "signalValue", "pValue", "qValue", 
        "fold_enrichment", "-log10(pvalue)", "-log10(qvalue)",
        "pileup", "tags", "-10*log10(pvalue)", "FDR",
        "signal.value", "p.value", "q.value", "fold.enrichment"
    )
    
    # Find matching columns (case-insensitive partial matching)
    detected_cols <- character(0)
    for (pattern in common_score_names) {
        # Exact match
        if (pattern %in% all_metadata_cols) {
            detected_cols <- c(detected_cols, pattern)
        } else {
            # Case-insensitive partial match
            pattern_lower <- tolower(pattern)
            pattern_escaped <- gsub("([()])", "\\\\\\1", pattern_lower)
            matches <- grep(pattern_escaped, all_metadata_cols, 
                           ignore.case = TRUE, value = TRUE)
            if (length(matches) > 0L) {
                detected_cols <- c(detected_cols, matches)
            }
        }
    }
    
    # Remove duplicates while preserving order
    return(unique(detected_cols))
}

# Internal function to plot peak widths
.plotWidths <- function(width_data, is_multi = FALSE) {
    if (nrow(width_data) == 0L) {
        return(NULL)
    }
    
    if (is_multi) {
        # Violin plot for multiple replicates
        p <- ggplot(width_data, aes(x = replicate, y = width, fill = replicate)) +
            geom_violin(alpha = 0.7, trim = FALSE) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, hjust = 1)) +
            xlab("Replicate") +
            ylab("Peak Width (bp)") +
            labs(title = "Peak Width Distribution Across Replicates")
    } else {
        # Violin plot for single peak set (using a dummy grouping variable)
        width_data_with_group <- width_data
        width_data_with_group$dummy <- "All Peaks"
        p <- ggplot(width_data_with_group, 
                   aes(x = .data$dummy, y = .data$width, fill = .data$dummy)) +
            geom_violin(alpha = 0.7, trim = FALSE, fill = "steelblue") +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()) +
            xlab("") +
            ylab("Peak Width (bp)") +
            labs(title = "Peak Width Distribution")
    }
    
    return(p)
}

# Internal function to plot scores
.plotScores <- function(score_data, col_name, is_multi = FALSE) {
    if (nrow(score_data) == 0L || length(score_data) == 0L) {
        return(NULL)
    }
    
    # Remove NA values for plotting
    score_data_clean <- score_data[!is.na(score_data$score), , drop = FALSE]
    if (nrow(score_data_clean) == 0L || length(score_data_clean) == 0L) {
        return(NULL)
    }
    
    if (is_multi) {
        # Violin plot for multiple replicates
        p <- ggplot(score_data_clean, 
                   aes(x = replicate, y = score, fill = replicate)) +
            geom_violin(alpha = 0.7, trim = FALSE) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, hjust = 1)) +
            xlab("Replicate") +
            ylab(col_name) +
            labs(title = paste("Score Distribution:", col_name))
    } else {
        # Violin plot for single peak set (using a dummy grouping variable)
        score_data_with_group <- score_data_clean
        score_data_with_group$dummy <- "All Peaks"
        p <- ggplot(score_data_with_group, 
                   aes(x = .data$dummy, y = .data$score, fill = .data$dummy)) +
            geom_violin(alpha = 0.7, trim = FALSE, fill = "steelblue") +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()) +
            xlab("") +
            ylab(col_name) +
            labs(title = paste("Score Distribution:", col_name))
    }
    
    p
}

