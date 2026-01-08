#' Condense matrix rows by grouping on a specified column
#' 
#' @description 
#' Condenses a matrix by grouping rows based on values in a specified column.
#' For each unique value in the grouping column, values in other columns are
#' concatenated with automatic deduplication (only unique values are kept) and
#' optionally counted. This is useful for aggregating data where multiple rows
#' share the same key value but have different associated values in other
#' columns.
#' 
#' @param mx A matrix to be condensed. Must have column names. The matrix will
#'        be grouped by one column, and all other columns will be condensed.
#' @param iname A character string or integer specifying the column to use for
#'        grouping rows. Can be specified as:
#'        \itemize{
#'          \item A character string: Column name (e.g., \code{"gene_id"})
#'          \item An integer: Column index (e.g., \code{1L} for the first
#'                column)
#'        }
#'        If a numeric index is provided, it is converted to the corresponding
#'        column name. The specified column must exist in the matrix.
#' @param sep A character string used as a separator when concatenating values
#'        within each group. Default is \code{";"}. This separator is used both
#'        for concatenation and for splitting when deduplicating values.
#' @param cnt A logical value. If \code{TRUE}, adds a count column for each
#'        condensed column showing the number of unique values that were
#'        aggregated. Count columns are named as \code{<column_name>.count}.
#'        Default is \code{FALSE}.
#' 
#' @return Returns a data frame with condensed rows. The structure is:
#'        \itemize{
#'          \item One row per unique value in the grouping column
#'          \item The grouping column (\code{iname}) is placed first
#'          \item All other columns are condensed (concatenated and deduplicated)
#'          \item If \code{cnt = TRUE}, each condensed column is followed by a
#'                corresponding count column (e.g., \code{column.count})
#'        }
#'        Column order: \code{iname}, then each condensed column (and its count
#'        column if \code{cnt = TRUE}).
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates that \code{mx} is a matrix with column names
#'   \item Converts \code{iname} to a column name if provided as an index
#'   \item Splits the matrix into groups based on unique values in the grouping
#'         column
#'   \item For each group and each non-grouping column:
#'         \itemize{
#'           \item Collects all values from that column for rows in the group
#'           \item Splits concatenated values by \code{sep} (if values were
#'                 already concatenated)
#'           \item Removes duplicates (keeps only unique values)
#'           \item Concatenates unique values back using \code{sep}
#'           \item If \code{cnt = TRUE}, counts the number of unique values
#'         }
#'   \item Combines all groups into a data frame with one row per group
#'   \item Reorders columns so the grouping column is first
#' }
#' 
#' \strong{Value deduplication:}
#' The function automatically deduplicates values before concatenation. This
#' means:
#' \itemize{
#'   \item If multiple rows in a group have the same value, it appears only
#'         once in the condensed result
#'   \item If values are already concatenated (contain \code{sep}), they are
#'         split, deduplicated, and re-concatenated
#'   \item This ensures each unique value appears exactly once per group
#' }
#' 
#' \strong{Count columns:}
#' When \code{cnt = TRUE}, for each condensed column, an additional count
#' column is added showing the number of unique values that were aggregated.
#' For example:
#' \itemize{
#'   \item If column \code{GO_term} is condensed, a column \code{GO_term.count}
#'         is added
#'   \item The count represents the number of unique GO terms for that group
#'   \item Count columns are placed immediately after their corresponding
#'         condensed columns
#' }
#' 
#' \strong{Use cases:}
#' This function is particularly useful for:
#' \itemize{
#'   \item Gene annotation data: Multiple rows per gene with different GO
#'         terms, pathways, or other annotations
#'   \item Peak annotation data: Multiple annotations per peak region
#'   \item Any tabular data where one-to-many relationships need to be
#'         collapsed into one row per key
#' }
#' 
#' \strong{Example transformation:}
#' Input matrix:
#' \preformatted{
#'   gene_id  GO_term
#'   G1       GO:0001
#'   G1       GO:0002
#'   G1       GO:0001
#'   G2       GO:0003
#' }
#' 
#' Output (with \code{sep = ";"}, \code{cnt = TRUE}):
#' \preformatted{
#'   gene_id  GO_term      GO_term.count
#'   G1       GO:0001;GO:0002  2
#'   G2       GO:0003          1
#' }
#' 
#' Note: Duplicate \code{GO:0001} for G1 is removed, and the count shows 2
#' unique GO terms.
#' 
#' \strong{Validation:}
#' The function performs the following checks:
#' \itemize{
#'   \item \code{mx} must be a matrix
#'   \item \code{mx} must have column names
#'   \item \code{iname} must be a single value (not a vector)
#'   \item The specified column must exist in the matrix
#' }
#' 
#' @note
#' \itemize{
#'   \item Values are always deduplicated, even if \code{cnt = FALSE}
#'   \item The separator \code{sep} is used for both concatenation and
#'         splitting, so choose a separator that doesn't appear in your data
#'         values
#'   \item If a column contains values that already include the separator, they
#'         will be split during deduplication, which may not be desired
#'   \item The function converts the matrix to a data frame for the output
#'   \item Row names of the output are the unique values from the grouping
#'         column (before reordering)
#' }
#' 
#' @author Jianhong Ou, Lihua Julie Zhu
#' @keywords misc
#' @export
#' 
#' @examples
#' 
#' # Example 1: Basic condensing by column name
#' a <- matrix(c(rep(rep(1:5, 2), 2), rep(1:10, 2)), ncol = 4)
#' colnames(a) <- c("con.1", "con.2", "index.1", "index.2")
#' condenseMatrixByColnames(a, "con.1")
#' 
#' # Example 2: Condense by column index
#' condenseMatrixByColnames(a, 2)
#' 
#' # Example 3: With count columns showing number of unique values
#' condenseMatrixByColnames(a, "con.1", cnt = TRUE)
#' 
#' # Example 4: Gene annotation example (more realistic)
#' gene_anno <- matrix(
#'     c(rep(c("G1", "G2", "G1"), 3),
#'       c("GO:0001", "GO:0002", "GO:0001", "GO:0003", "GO:0004", "GO:0001",
#'         "pathway1", "pathway2", "pathway1", "pathway3", "pathway4", "pathway1")),
#'     ncol = 3,
#'     dimnames = list(NULL, c("gene_id", "GO_term", "pathway"))
#' )
#' # Condense: one row per gene with all GO terms and pathways
#' condenseMatrixByColnames(gene_anno, "gene_id")
#' 
#' # Example 5: With counts to see how many unique annotations per gene
#' condenseMatrixByColnames(gene_anno, "gene_id", cnt = TRUE)
#' 
#' # Example 6: Custom separator
#' condenseMatrixByColnames(gene_anno, "gene_id", sep = "|")
condenseMatrixByColnames <- function(mx, iname, sep = ";", cnt = FALSE) {
    if (!is.matrix(mx)) {
        stop("'mx' must be a matrix", call. = FALSE)
    }
    
    if (length(iname) != 1L) {
        stop("'iname' must be a single column name or index", call. = FALSE)
    }
    
    m_cname <- colnames(mx)
    if (is.null(m_cname)) {
        stop("'mx' must have column names", call. = FALSE)
    }
    
    # Convert numeric index to column name
    if (is.numeric(iname) && iname <= length(m_cname)) {
        iname <- m_cname[as.integer(iname)]
    }
    
    cnames <- m_cname[m_cname != iname]
    
    if (length(m_cname) == length(cnames)) {
        stop("the column name specified for condense does not exist", 
             call. = FALSE)
    }
    
    # Split matrix by grouping column
    m_split <- split(mx[, cnames, drop = FALSE], mx[, iname])
    col_n <- length(cnames)
    
    m_list <- lapply(m_split, function(.ele) {
        x <- apply(matrix(.ele, nrow = col_n, byrow = TRUE),
                   1L,
                   paste,
                   collapse = sep)
        if (cnt) {
            unlist(lapply(x, function(w) {
                tmp <- unique(as.character(unlist(strsplit(w,
                                                           sep,
                                                           fixed = TRUE))))
                c(paste(tmp, collapse = sep), length(tmp))
            }))
        } else {
            unlist(lapply(x, function(w) {
                tmp <- unique(as.character(unlist(strsplit(w,
                                                           sep,
                                                           fixed = TRUE))))
                paste(tmp, collapse = sep)
            }))
        }
    })
    
    m_dat <- as.data.frame(do.call(rbind, m_list))
    m_dat$index <- rownames(m_dat)
    
    # Build column names
    cnames_cnt <- character(0)
    for (i in cnames) {
        if (cnt) {
            cnames_cnt <- c(cnames_cnt, i, paste0(i, ".count"))
        } else {
            cnames_cnt <- c(cnames_cnt, i)
        }
    }
    
    colnames(m_dat) <- c(cnames_cnt, iname)
    m_dat <- m_dat[, c(iname, cnames_cnt), drop = FALSE]
    
    m_dat
}
