#' Global variables for package
#' 
#' @description 
#' Declares global variables used in package functions to avoid R CMD check
#' warnings. These variables are typically column names in data frames created
#' dynamically within functions (e.g., in \code{enrichmentPlot}).
#' 
#' @name globalVariables
#' @keywords internal
NULL

utils::globalVariables(
    c("Description", "pvalue", "qvalue", "Count", "GeneRatio", "base")
)