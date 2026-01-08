#' Count occurrences of unique GO IDs
#' 
#' @description 
#' Counts the number of occurrences for each unique GO ID in a character vector.
#' This function efficiently aggregates GO term counts by sorting the input vector
#' and counting consecutive occurrences of each unique GO term. It is used internally
#' by enrichment analysis functions (e.g., \code{\link{getEnrichedGO}}) to aggregate
#' GO term counts before performing statistical tests.
#' 
#' The function handles edge cases gracefully: if the input is not a character
#' vector or is empty, it returns empty vectors for both \code{GOterm} and
#' \code{GOcount}.
#' 
#' @param goList A character vector containing GO IDs (e.g., "GO:0000075",
#'        "GO:0000082"). The vector may contain duplicates. GO IDs are typically
#'        in the format "GO:XXXXXXXX" but the function accepts any character
#'        strings. If \code{goList} is not a character vector or has length 0,
#'        the function returns empty results.
#' 
#' @return Returns a list with two elements:
#'        \itemize{
#'          \item \code{GOterm}: A character vector of unique GO terms, sorted
#'                alphabetically. Each GO term appears only once, even if it
#'                occurred multiple times in the input.
#'          \item \code{GOcount}: A numeric vector of the same length as
#'                \code{GOterm}, containing the count of occurrences for each
#'                corresponding GO term in the input. The counts are in the same
#'                order as the GO terms in \code{GOterm}.
#'        }
#'        If the input is empty or invalid, both vectors will be empty (length 0).
#' 
#' @details
#' The function uses an efficient algorithm:
#' \enumerate{
#'   \item Sorts the input GO list alphabetically
#'   \item Identifies unique GO terms
#'   \item Counts consecutive occurrences of each GO term in the sorted list
#'   \item Returns unique GO terms and their counts in matching order
#' }
#' 
#' This approach is more efficient than using \code{table()} for large vectors
#' because it processes the sorted list in a single pass.
#' 
#' @note This is an internal function primarily used by enrichment analysis
#'       functions. While it is exported and can be called directly, users
#'       typically interact with higher-level functions like
#'       \code{\link{getEnrichedGO}}.
#' 
#' @author Lihua Julie Zhu
#' @seealso \code{\link{getEnrichedGO}} for GO enrichment analysis,
#'          \code{\link{getEnrichedPATH}} for pathway enrichment analysis
#' @keywords internal
#' @export
#' @examples
#' ## Example 1: Basic usage with duplicate GO terms
#' goList <- c("GO:0000075", "GO:0000082", "GO:0000082", "GO:0000122", 
#'             "GO:0000122", "GO:0000075", "GO:0000082", "GO:0000082", 
#'             "GO:0000122", "GO:0000122", "GO:0000122", "GO:0000122", 
#'             "GO:0000075", "GO:0000082", "GO:000012")
#' result <- getUniqueGOidCount(goList)
#' result
#' ## $GOterm
#' ## [1] "GO:0000075" "GO:0000082" "GO:000012"  "GO:0000122"
#' ## 
#' ## $GOcount
#' ## [1] 3 5 1 6
#' 
#' ## Example 2: Empty input
#' getUniqueGOidCount(character(0))
#' ## $GOterm
#' ## character(0)
#' ## 
#' ## $GOcount
#' ## numeric(0)
#' 
#' ## Example 3: All unique GO terms (no duplicates)
#' getUniqueGOidCount(c("GO:0000001", "GO:0000002", "GO:0000003"))
#' ## $GOterm
#' ## [1] "GO:0000001" "GO:0000002" "GO:0000003"
#' ## 
#' ## $GOcount
#' ## [1] 1 1 1
#' 
getUniqueGOidCount <- function(goList) {
    if (!is.character(goList) || length(goList) == 0L) {
        return(list(GOterm = character(0), GOcount = numeric(0)))
    }
    
    x <- goList[order(goList)]
    duplicated_go <- duplicated(x)
    unique_go <- unique(x)
    n_unique <- length(unique_go)
    go_count <- numeric(n_unique)
    
    count <- 1L
    j <- 1L
    
    for (i in seq.int(2L, length(duplicated_go))) {
        if (!duplicated_go[i]) {
            go_count[j] <- count  # previous GO
            j <- j + 1L
            count <- 1L
        } else {
            count <- count + 1L
        }
    }
    go_count[j] <- count  # last GO
    
    list(GOterm = unique_go, GOcount = go_count)
}

