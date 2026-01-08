#' Calculate z-scores for log2 ratio between two numeric vectors
#' 
#' Calculate z-scores for the log2-transformed ratio between two numeric
#' vectors. This is useful for normalizing fold-change values (e.g., ChIP-seq
#' signal ratios) to a standard normal distribution.
#' 
#' @param A,B Numeric vectors of equal length. Typically representing signal
#'        values (e.g., ChIP-seq read counts) in two conditions.
#' @param background Numeric. Background value to add before log2
#'        transformation to avoid log(0). Default is \code{1}.
#' 
#' @return Numeric vector of z-scores with the same length as input vectors.
#'        Z-scores are calculated as: \code{(log2_ratio - mean) / sd}, where
#'        \code{log2_ratio = log2(A + background) - log2(B + background)}.
#' 
#' @details
#' This function:
#' \enumerate{
#'   \item Calculates log2 ratios: \code{log2(A + background) - log2(B + background)}
#'   \item Computes population standard deviation (using \code{n-1} correction)
#'   \item Standardizes to z-scores: \code{(value - mean) / sd}
#' }
#' 
#' The resulting z-scores follow a standard normal distribution (mean = 0, sd = 1),
#' which is useful for identifying outliers or significant changes.
#' 
#' @note
#' This function assumes the log2 ratios follow a normal distribution. For
#' non-normal distributions, consider using rank-based methods or other
#' transformations.
#' 
#' @author Internal utility function
#' @keywords internal
#' 
#' @examples
#' # Example: ChIP-seq signal in two conditions
#' condition1 <- c(100, 200, 150, 300, 250)
#' condition2 <- c(50, 100, 75, 150, 125)
#' 
#' z_scores <- ratio.zScore(condition1, condition2)
#' z_scores
ratio.zScore <- function(A, B, background = 1) {
    stopifnot(
        is.numeric(A),
        is.numeric(B),
        length(A) == length(B),
        is.numeric(background),
        length(background) == 1L,
        background >= 0
    )
    
    # Calculate log2 ratio
    r <- log2(A + background) - log2(B + background)
    
    # Calculate population statistics
    n <- length(r)
    pop_mean <- mean(r)
    pop_sd <- sd(r) * sqrt((n - 1L) / n)
    
    # Calculate z-scores
    z <- (r - pop_mean) / pop_sd
    
    z
}