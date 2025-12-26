#' Add metadata columns to overlapping peaks
#' 
#' Adds metadata columns from the original GRanges objects to the overlapping
#' peaks result object returned by \code{\link{findOverlapsOfPeaks}}. This
#' function aggregates metadata values from overlapping peaks using a specified
#' function (default: \code{c} for concatenation).
#' 
#' @param ol An object of class \code{overlappingPeaks}, which is the output
#'   of \code{\link{findOverlapsOfPeaks}}.
#' @param colNames Character vector of metadata column names to be added. If
#'   \code{NULL} (default), the function will automatically detect columns that
#'   are present in all input peak sets and add those.
#' @param FUN A function used to aggregate metadata values when multiple peaks
#'   overlap. Default is \code{c} (concatenation). Other useful functions
#'   include \code{mean}, \code{sum}, \code{max}, \code{min}, or custom
#'   functions. The function should accept a vector and return a single value.
#' @param ... Additional arguments passed to \code{FUN}.
#' 
#' @return An object of class \code{overlappingPeaks} with metadata columns
#'   added to the \code{peaklist} element. The structure is the same as the
#'   input, but each element in \code{peaklist} now contains the aggregated
#'   metadata columns from the overlapping peaks.
#' 
#' @details
#' This function is useful when you want to preserve metadata (e.g., scores,
#' p-values, fold changes) from the original peak sets in the overlap analysis
#' results. When multiple peaks overlap, the metadata values are aggregated
#' using the specified function.
#' 
#' \strong{Automatic column detection:} If \code{colNames} is \code{NULL}, the
#' function identifies columns that exist in all input peak sets and adds only
#' those columns. This ensures consistency across peak sets.
#' 
#' \strong{Metadata aggregation:} When peaks from different sets overlap, their
#' metadata values are combined using \code{FUN}. For example:
#' \itemize{
#'   \item \code{FUN = c}: Concatenates values (useful for IDs, names)
#'   \item \code{FUN = mean}: Computes mean (useful for scores, signals)
#'   \item \code{FUN = max}: Takes maximum (useful for p-values, confidence)
#' }
#' 
#' @note
#' \itemize{
#'   \item All specified columns must exist in all input peak sets.
#'   \item The classes of metadata columns must be identical across all peak
#'     sets.
#'   \item The function modifies the input object in place (returns the same
#'     object with added metadata).
#' }
#' 
#' @export
#' @importFrom S4Vectors mcols aggregate
#' @importFrom GenomicRanges GRangesList
#' @author Jianhong Ou
#' @seealso \code{\link{findOverlapsOfPeaks}} for creating the overlap object
#' @keywords misc
#' @examples
#' 
#' # Create example peak sets with metadata
#' peaks1 <- GRanges(
#'     seqnames = c(6, 6, 6, 6, 5),
#'     IRanges(
#'         start = c(1543200, 1557200, 1563000, 1569800, 167889600),
#'         end = c(1555199, 1560599, 1565199, 1573799, 167893599),
#'         names = c("p1", "p2", "p3", "p4", "p5")
#'     ),
#'     strand = "+",
#'     score = 1:5,
#'     id = letters[1:5],
#'     pvalue = c(0.001, 0.01, 0.05, 0.1, 0.2)
#' )
#' 
#' peaks2 <- GRanges(
#'     seqnames = c(6, 6, 6, 6, 5),
#'     IRanges(
#'         start = c(1549800, 1554400, 1565000, 1569400, 167888600),
#'         end = c(1550599, 1560799, 1565399, 1571199, 167888999),
#'         names = c("f1", "f2", "f3", "f4", "f5")
#'     ),
#'     strand = "+",
#'     score = 6:10,
#'     id = LETTERS[1:5],
#'     pvalue = c(0.002, 0.02, 0.03, 0.15, 0.25)
#' )
#' 
#' # Find overlaps
#' ol <- findOverlapsOfPeaks(peaks1, peaks2)
#' 
#' # Add all common metadata columns (automatic detection)
#' ol_with_metadata <- addMetadata(ol)
#' 
#' # Add specific columns with custom aggregation function
#' ol_mean_scores <- addMetadata(ol, colNames = "score", FUN = mean)
#' ol_max_pvalues <- addMetadata(ol, colNames = "pvalue", FUN = min)
#' 
addMetadata <- function(ol, colNames = NULL, FUN = c, ...) {
    # Input validation
    if (!inherits(ol, "overlappingPeaks")) {
        stop("'ol' must be an object of class 'overlappingPeaks' ",
             "(output from findOverlapsOfPeaks)")
    }
    
    if (length(ol$all.peaks) == 0L) {
        stop("'ol$all.peaks' is empty")
    }
    
    if (!is.function(FUN)) {
        stop("'FUN' must be a function")
    }
    
    peaks_list <- ol$all.peaks
    
    # Auto-detect common columns if colNames is NULL
    if (is.null(colNames)) {
        all_colnames <- lapply(peaks_list, function(.ele) {
            colnames(mcols(.ele))
        })
        colname_counts <- table(unlist(all_colnames))
        # Only include columns present in all peak sets
        colNames <- names(colname_counts)[colname_counts == length(peaks_list)]
        
        if (length(colNames) == 0L) {
            stop("No common metadata columns found across all peak sets. ",
                 "Please specify 'colNames' explicitly.")
        }
    }
    
    # Validate colNames
    if (!is.character(colNames) || length(colNames) == 0L) {
        stop("'colNames' must be a non-empty character vector")
    }
    
    # Check that all specified columns exist in all peak sets
    for (i in seq_along(peaks_list)) {
        missing_cols <- colNames[!colNames %in% colnames(mcols(peaks_list[[i]]))]
        if (length(missing_cols) > 0L) {
            stop("Column(s) '", paste(missing_cols, collapse = "', '"),
                 "' not found in peak set ", i)
        }
    }
    
    # Extract only specified columns from each peak set
    peaks_list_subset <- lapply(peaks_list, function(.ele) {
        .ele[, colNames]
    })
    
    # Verify that column classes are identical across all peak sets
    col_classes <- lapply(peaks_list_subset, function(.ele) {
        vapply(mcols(.ele), class, character(1))
    })
    
    for (i in seq.int(2L, length(col_classes))) {
        if (!identical(col_classes[[1L]], col_classes[[i]])) {
            stop("Metadata column classes are not identical across peak sets. ",
                 "Column classes must match for aggregation.")
        }
    }
    
    # Combine all peaks into a single GRanges object
    all_peaks <- unlist(GRangesList(peaks_list_subset), use.names = FALSE)
    
    # Add aggregated metadata to each overlap group
    ol$peaklist <- lapply(ol$peaklist, function(.ele) {
        # Get peak indices for this overlap group
        peak_indices <- unlist(.ele$peakNames)
        
        # Aggregate metadata for overlapping peaks
        group_ids <- rep(seq_along(.ele), lengths(.ele$peakNames))
        aggregated_mcols <- aggregate(
            mcols(all_peaks[peak_indices]),
            by = list(addMetadata_group = group_ids),
            FUN = FUN,
            ...
        )
        
        # Sort by group ID to maintain order
        aggregated_mcols <- aggregated_mcols[order(aggregated_mcols$addMetadata_group), ]
        
        # Add aggregated columns to the overlap result
        n_existing_cols <- ncol(mcols(.ele))
        mcols(.ele) <- cbind(mcols(.ele), aggregated_mcols[, colNames, drop = FALSE])
        
        # Update column names
        if (n_existing_cols > 0L) {
            colnames(mcols(.ele))[-seq.int(n_existing_cols)] <- colNames
        } else {
            colnames(mcols(.ele)) <- colNames
        }
        
        .ele
    })
    
    ol
}
