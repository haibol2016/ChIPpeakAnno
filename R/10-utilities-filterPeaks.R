#' Filter peaks based on width and/or score thresholds
#' 
#' @description 
#' Filters a GRanges object containing peaks based on peak width and/or score
#' thresholds. This function is useful for quality control, removing outliers,
#' or selecting peaks within specific size or signal strength ranges.
#' 
#' The function supports filtering by:
#' \itemize{
#'   \item \strong{Peak width}: Minimum and/or maximum width in base pairs
#'   \item \strong{Score columns}: Minimum and/or maximum values for any
#'         metadata column containing numeric scores (e.g., "score",
#'         "signalValue", "pValue", "qValue", "fold_enrichment")
#' }
#' 
#' All thresholds are optional (default \code{NULL}), allowing flexible
#' filtering strategies. Peaks that pass all specified thresholds are retained.
#' 
#' @param peaks A \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing peaks to be filtered.
#' @param min_width An integer or \code{NULL} (default). Minimum peak width in
#'        base pairs. Peaks with width < \code{min_width} will be removed. If
#'        \code{NULL}, no width-based filtering is applied.
#' @param max_width An integer or \code{NULL} (default). Maximum peak width in
#'        base pairs. Peaks with width > \code{max_width} will be removed. If
#'        \code{NULL}, no width-based filtering is applied.
#' @param score_thresholds A named list or \code{NULL} (default). Each element
#'        should be a numeric vector of length 2: \code{c(min, max)}. The
#'        names should correspond to metadata column names in the GRanges
#'        object. For example:
#'        \code{list(score = c(100, 1000), signalValue = c(5, Inf))}
#'        filters peaks where:
#'        \itemize{
#'          \item \code{score} is between 100 and 1000 (inclusive)
#'          \item \code{signalValue} is >= 5 (no upper limit)
#'        }
#'        If a threshold vector contains \code{NULL} or \code{NA}, that limit
#'        is ignored. For example, \code{list(score = c(100, NULL))} filters
#'        peaks with \code{score >= 100} (no upper limit).
#' @param verbose Logical. If \code{TRUE} (default), prints filtering summary
#'        to the console. If \code{FALSE}, only returns the filtered peaks
#'        without printing.
#' 
#' @return Returns a \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing only peaks that pass all specified thresholds. The
#'        returned object:
#'        \itemize{
#'          \item Preserves all original metadata columns
#'          \item Maintains peak names and other attributes
#'          \item Returns an empty GRanges object if no peaks pass the filters
#'        }
#' 
#' @details
#' 
#' \strong{Filtering Logic:}
#' 
#' The function applies filters in the following order:
#' \enumerate{
#'   \item \strong{Width filtering}: If \code{min_width} or \code{max_width} is
#'         specified, peaks are filtered by width (calculated as
#'         \code{width(peaks)} = \code{end - start + 1})
#'   \item \strong{Score filtering}: For each score column specified in
#'         \code{score_thresholds}, peaks are filtered by the corresponding
#'         min/max values
#' }
#' 
#' Peaks must pass \strong{all} specified thresholds to be retained. If a peak
#' fails any threshold, it is removed.
#' 
#' \strong{Handling Missing Values:}
#' 
#' - Peaks with \code{NA} or \code{NaN} values in score columns are removed if
#'   that score column is used for filtering
#' - Peaks with invalid widths (e.g., \code{NA}, negative, or zero width) are
#'   removed if width filtering is applied
#' 
#' \strong{Common Use Cases:}
#' 
#' \itemize{
#'   \item \strong{Remove narrow peaks}: \code{min_width = 50} removes peaks
#'         narrower than 50 bp (potential artifacts)
#'   \item \strong{Remove broad peaks}: \code{max_width = 5000} removes peaks
#'         wider than 5 kb (may represent large domains or artifacts)
#'   \item \strong{Filter by significance}: \code{score_thresholds = list(pValue
#'         = c(0.05, NULL))} keeps only peaks with p-value < 0.05
#'   \item \strong{Filter by signal strength}: \code{score_thresholds =
#'         list(signalValue = c(10, Inf))} keeps only peaks with signalValue >=
#'         10
#'   \item \strong{Combined filtering}: Filter by both width and score to
#'         obtain high-confidence peaks
#' }
#' 
#' @author Haibo Liu
#' @seealso \code{\link{assessPeaks}} for assessing peak widths and scores
#'          before filtering, \code{\link{IDRfilter}} for IDR-based filtering
#'          of replicate peaks
#' @importFrom S4Vectors mcols
#' @export
#' @examples
#' \dontrun{
#' # Load example peaks
#' data(myPeakList)
#' peaks <- myPeakList[1:100]
#' 
#' # Add some score columns for demonstration
#' mcols(peaks)$score <- sample(0:1000, length(peaks), replace = TRUE)
#' mcols(peaks)$signalValue <- rnorm(length(peaks), mean = 10, sd = 3)
#' mcols(peaks)$pValue <- runif(length(peaks), min = 0, max = 1)
#' 
#' # Example 1: Filter by width only
#' filtered1 <- filterPeaks(peaks, min_width = 100, max_width = 1000)
#' 
#' # Example 2: Filter by score only
#' filtered2 <- filterPeaks(peaks, 
#'                         score_thresholds = list(score = c(100, 500)))
#' 
#' # Example 3: Filter by multiple score columns
#' filtered3 <- filterPeaks(peaks,
#'                         score_thresholds = list(
#'                             score = c(100, 500),
#'                             signalValue = c(5, Inf),
#'                             pValue = c(NULL, 0.05)  # pValue < 0.05
#'                         ))
#' 
#' # Example 4: Combined width and score filtering
#' filtered4 <- filterPeaks(peaks,
#'                         min_width = 100,
#'                         max_width = 1000,
#'                         score_thresholds = list(
#'                             score = c(200, NULL),
#'                             signalValue = c(8, Inf)
#'                         ))
#' 
#' # Example 5: Silent mode (no console output)
#' filtered5 <- filterPeaks(peaks, min_width = 100, verbose = FALSE)
#' }
filterPeaks <- function(peaks, min_width = NULL, max_width = NULL,
                        score_thresholds = NULL, verbose = TRUE) {
    # Input validation
    if (missing(peaks)) {
        stop("Missing required argument 'peaks'!", call. = FALSE)
    }
    if (!inherits(peaks, "GRanges")) {
        stop("'peaks' must be a GRanges object", call. = FALSE)
    }
    
    if (length(peaks) == 0L) {
        warning("Input GRanges object is empty.", call. = FALSE)
        return(peaks)
    }
    
    original_count <- length(peaks)
    filtered_peaks <- peaks
    filter_reasons <- character(0)
    
    # Filter by width
    if (!is.null(min_width) || !is.null(max_width)) {
        peak_widths <- width(filtered_peaks)
        
        # Create width filter mask
        width_mask <- rep(TRUE, length(filtered_peaks))
        
        if (!is.null(min_width)) {
            if (!is.numeric(min_width) || length(min_width) != 1L) {
                stop("'min_width' must be a single numeric value", call. = FALSE)
            }
            width_mask <- width_mask & (peak_widths >= min_width)
            filter_reasons <- c(filter_reasons, 
                               paste0("width >= ", min_width, " bp"))
        }
        
        if (!is.null(max_width)) {
            if (!is.numeric(max_width) || length(max_width) != 1L) {
                stop("'max_width' must be a single numeric value", call. = FALSE)
            }
            width_mask <- width_mask & (peak_widths <= max_width)
            filter_reasons <- c(filter_reasons, 
                               paste0("width <= ", max_width, " bp"))
        }
        
        # Also remove peaks with invalid widths (NA, NaN, negative, zero)
        width_mask <- width_mask & !is.na(peak_widths) & 
                     !is.nan(peak_widths) & peak_widths > 0
        
        filtered_peaks <- filtered_peaks[width_mask]
        
        if (verbose) {
            removed <- original_count - length(filtered_peaks)
            cat("Width filtering: ", length(filtered_peaks), " of ", 
               original_count, " peaks retained", 
               if (removed > 0L) paste0(" (", removed, " removed)") else "",
               "\n", sep = "")
        }
    }
    
    # Filter by score thresholds
    if (!is.null(score_thresholds)) {
        if (!is.list(score_thresholds)) {
            stop("'score_thresholds' must be a named list", call. = FALSE)
        }
        
        # Get available metadata columns
        available_cols <- colnames(mcols(filtered_peaks))
        
        for (col_name in names(score_thresholds)) {
            if (!col_name %in% available_cols) {
                warning("Score column '", col_name, "' not found in metadata. ",
                       "Skipping.", call. = FALSE)
                next
            }
            
            threshold <- score_thresholds[[col_name]]
            
            # Validate threshold format
            if (!is.numeric(threshold) && !is.null(threshold)) {
                stop("Threshold for '", col_name, "' must be numeric or NULL",
                     call. = FALSE)
            }
            
            # Handle threshold vector
            if (length(threshold) == 0L || all(is.na(threshold))) {
                next  # Skip if threshold is empty or all NA
            }
            
            # Extract min and max
            # Handle NULL values: c(NULL, 5) becomes c(NA, 5) in R
            # Handle Inf values: c(5, Inf) means min=5, no max
            min_score <- NULL
            max_score <- NULL
            
            if (length(threshold) >= 1L) {
                val1 <- threshold[1L]
                if (!is.na(val1) && !is.null(val1) && is.finite(val1)) {
                    min_score <- val1
                } else if (is.infinite(val1) && val1 > 0) {
                    # Inf as min means no lower limit
                    min_score <- NULL
                }
            }
            if (length(threshold) >= 2L) {
                val2 <- threshold[2L]
                if (!is.na(val2) && !is.null(val2) && is.finite(val2)) {
                    max_score <- val2
                } else if (is.infinite(val2) && val2 > 0) {
                    # Inf as max means no upper limit
                    max_score <- NULL
                }
            }
            
            # Get score values
            score_values <- mcols(filtered_peaks)[[col_name]]
            
            # Convert to numeric if needed
            if (!is.numeric(score_values)) {
                score_values <- tryCatch(
                    as.numeric(score_values),
                    warning = function(w) {
                        warning("Column '", col_name, "' could not be ",
                               "converted to numeric. Skipping.", 
                               call. = FALSE)
                        return(NULL)
                    }
                )
                if (is.null(score_values)) {
                    next
                }
            }
            
            # Create score filter mask
            score_mask <- rep(TRUE, length(filtered_peaks))
            
            # Apply min threshold
            if (!is.null(min_score)) {
                if (is.infinite(min_score)) {
                    # Handle Inf case
                    if (min_score > 0) {
                        # min = Inf means no lower limit
                        min_score <- NULL
                    } else {
                        # min = -Inf means all values pass
                        score_mask <- rep(TRUE, length(filtered_peaks))
                    }
                } else {
                    score_mask <- score_mask & (score_values >= min_score)
                }
                filter_reasons <- c(filter_reasons, 
                                   paste0(col_name, " >= ", min_score))
            }
            
            # Apply max threshold
            if (!is.null(max_score)) {
                if (is.infinite(max_score)) {
                    # Handle Inf case
                    if (max_score < 0) {
                        # max = -Inf means no upper limit
                        max_score <- NULL
                    } else {
                        # max = Inf means all values pass
                        # (already handled by score_mask initialization)
                    }
                } else {
                    score_mask <- score_mask & (score_values <= max_score)
                }
                filter_reasons <- c(filter_reasons, 
                                   paste0(col_name, " <= ", max_score))
            }
            
            # Remove peaks with NA/NaN scores
            score_mask <- score_mask & !is.na(score_values) & 
                         !is.nan(score_values)
            
            # Apply filter
            before_count <- length(filtered_peaks)
            filtered_peaks <- filtered_peaks[score_mask]
            
            if (verbose) {
                removed <- before_count - length(filtered_peaks)
                cat("Score filtering (", col_name, "): ", length(filtered_peaks),
                   " of ", before_count, " peaks retained",
                   if (removed > 0L) paste0(" (", removed, " removed)") else "",
                   "\n", sep = "")
            }
        }
    }
    
    # Print summary
    if (verbose) {
        final_count <- length(filtered_peaks)
        total_removed <- original_count - final_count
        cat("\n=== Filtering Summary ===\n")
        cat("Original peaks:", original_count, "\n")
        cat("Filtered peaks:", final_count, "\n")
        cat("Removed:", total_removed, 
           if (total_removed > 0L) {
               paste0(" (", round(100 * total_removed / original_count, 1), "%)")
           } else {
               ""
           }, "\n")
        if (length(filter_reasons) > 0L) {
            cat("Applied filters:", paste(filter_reasons, collapse = ", "), "\n")
        }
        cat("========================\n")
    }
    
    return(filtered_peaks)
}

