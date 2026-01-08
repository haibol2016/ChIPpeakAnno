## funtion not needed

#' Add GO IDs of the ancestors for a given vector of GO IDs
#' 
#' @description 
#' Adds GO IDs of the ancestors for a given vector of GO IDs leveraging the
#' \code{GO.db} package. This function expands a set of GO terms by including
#' their ancestor terms from the Gene Ontology hierarchy. For each input GO ID,
#' all ancestor terms (parent terms up to the root) are retrieved and associated
#' with the same Entrez IDs as the original GO ID.
#' 
#' The function uses the GO ancestor mappings from \code{GO.db}:
#' \code{GOBPANCESTOR} for biological process, \code{GOCCANCESTOR} for cellular
#' component, and \code{GOMFANCESTOR} for molecular function. The root term
#' "all" is automatically excluded from the results.
#' 
#' @param go.ids A character matrix with at least 4 columns where:
#'   \itemize{
#'     \item Column 1: GO IDs (character) - the GO term identifiers to expand
#'     \item Column 2: Evidence codes (not used in this function, but required
#'           for input format compatibility)
#'     \item Column 3: Ontology type (not used in this function, but required
#'           for input format compatibility)
#'     \item Column 4: Entrez IDs (character) - gene identifiers associated
#'           with each GO ID
#'   }
#'   The matrix should contain unique GO ID and Entrez ID pairs. This format
#'   is typically obtained from \code{\link{getEnrichedGO}} or similar
#'   enrichment analysis functions.
#' @param ontology A character string specifying the Gene Ontology type.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"bp"} (default): Biological process ontology
#'     \item \code{"cc"}: Cellular component ontology
#'     \item \code{"mf"}: Molecular function ontology
#'   }
#'   The ontology type determines which ancestor mapping is used from
#'   \code{GO.db}.
#' 
#' @return Returns a character matrix with 2 columns:
#'   \itemize{
#'     \item \code{GO_ID}: GO term identifiers, including both the original
#'           input GO IDs and all their ancestor terms from the GO hierarchy
#'     \item \code{EntrezID}: Entrez gene identifiers, associated with each
#'           GO ID (both original and ancestor terms)
#'   }
#'   The matrix contains unique rows (duplicates are removed) and excludes
#'   any rows with empty or NA GO IDs. Each ancestor GO term is associated
#'   with the same Entrez IDs as its descendant term(s).
#' 
#'   \strong{Special cases:}
#'   \itemize{
#'     \item If no ancestors are found for any input GO IDs, the function
#'           returns the original input GO IDs with their Entrez IDs
#'     \item If the result matrix has fewer than 3 rows, the function returns
#'           only the original input GO IDs (without ancestors)
#'     \item The root GO term "all" is always excluded from the results
#'   }
#' 
#' @details
#' This function is useful for expanding GO enrichment results to include
#' broader (more general) GO terms in the hierarchy. For example, if a gene
#' is annotated to "DNA repair" (GO:0006281), adding ancestors would also
#' include broader terms like "DNA metabolic process" (GO:0006259) and
#' "metabolic process" (GO:0008152), all associated with the same gene.
#' 
#' The function requires the \code{GO.db} package to be installed. If it's
#' not available, the function will stop with an error message.
#' 
#' @export
#' @author Lihua Julie Zhu, Haibo Liu
#' @keywords misc
#' @importFrom GO.db GOBPANCESTOR GOCCANCESTOR GOMFANCESTOR
#' @seealso \code{\link{getEnrichedGO}} for obtaining GO enrichment results
#'          in the required format, \code{\link[GO.db]{GOBPANCESTOR}} for
#'          biological process ancestor mappings
#' @examples 
#' \dontrun{
#' ## Example 1: Add ancestors for biological process ontology
#' ## First, get enriched GO terms (example format)
#' ## go.ids <- getEnrichedGO(peaks, annoData, ...)
#' 
#' ## Add ancestors to expand the GO term set
#' result <- addAncestors(go.ids, ontology = "bp")
#' head(result)
#' 
#' ## Example 2: Add ancestors for molecular function ontology
#' result <- addAncestors(go.ids, ontology = "mf")
#' 
#' ## Example 3: Add ancestors for cellular component ontology
#' result <- addAncestors(go.ids, ontology = "cc")
#' }
#' 

addAncestors <- function(go.ids, ontology = c("bp", "cc", "mf")) {
    # Check for required package
    stopifnot("The 'GO.db' package is required" = 
              requireNamespace("GO.db", quietly = TRUE))
    
    # Validate and match ontology argument
    ontology <- match.arg(ontology)
    
    # Validate go.ids parameter
    if (missing(go.ids)) {
        stop("Missing required parameter 'go.ids'")
    }
    
    if (!is.matrix(go.ids) || ncol(go.ids) < 4L) {
        stop("'go.ids' must be a matrix with at least 4 columns.\n",
             "  Column 1: GO IDs\n",
             "  Column 2: Evidence codes\n",
             "  Column 3: Ontology type\n",
             "  Column 4: Entrez IDs", call. = FALSE)
    }
    
    # Get GO ancestor mappings based on ontology
    ancestor_maps <- switch(
        ontology,
        "bp" = as.list(GO.db::GOBPANCESTOR),
        "cc" = as.list(GO.db::GOCCANCESTOR),
        "mf" = as.list(GO.db::GOMFANCESTOR)
    )
    
    # Remove NA entries
    ancestor_maps <- ancestor_maps[!is.na(ancestor_maps)]
    
    if (length(ancestor_maps) == 0) {
        stop("No GO ancestors found for ontology '", ontology, 
             "'. This should not happen - please check GO.db installation.")
    }
    
    # Extract GO IDs and Entrez IDs from input
    input_go_ids <- as.character(go.ids[, 1])
    input_entrez_ids <- as.character(go.ids[, 4])
    
    # Find GO IDs that have ancestors in the mapping
    go_ids_with_ancestors <- intersect(input_go_ids, names(ancestor_maps))
    
    # If no GO IDs have ancestors, return original input
    if (length(go_ids_with_ancestors) == 0) {
        result <- unique(cbind(
            GO_ID = input_go_ids,
            EntrezID = input_entrez_ids
        ))
        return(result)
    }
    
    # Extract ancestors for each GO ID
    ancestors_list <- lapply(go_ids_with_ancestors, function(go_id) {
        ancestor_set <- ancestor_maps[[go_id]]
        if (length(ancestor_set) > 0) {
            # Remove "all" term (root of GO hierarchy)
            ancestor_set <- ancestor_set[ancestor_set != "all"]
            if (length(ancestor_set) > 0) {
                return(cbind(
                    child = rep(go_id, length(ancestor_set)),
                    ancestor = ancestor_set
                ))
            }
        }
        return(NULL)
    })
    
    # Combine all ancestors into a matrix
    ancestors_matrix <- do.call(rbind, ancestors_list[!sapply(ancestors_list, is.null)])
    
    if (is.null(ancestors_matrix) || nrow(ancestors_matrix) == 0) {
        # No ancestors found, return original input
        result <- unique(cbind(
            GO_ID = input_go_ids,
            EntrezID = input_entrez_ids
        ))
        return(result)
    }
    
    # Create data frames for merging
    children_df <- data.frame(
        child = input_go_ids,
        entrezID = input_entrez_ids,
        stringsAsFactors = FALSE
    )
    
    ancestors_df <- data.frame(
        child = ancestors_matrix[, "child"],
        ancestor = ancestors_matrix[, "ancestor"],
        stringsAsFactors = FALSE
    )
    
    # Merge children with their ancestors
    go_all <- merge(children_df, ancestors_df, by = "child", all.x = TRUE)
    
    # Combine original GO IDs and ancestor GO IDs with their Entrez IDs
    result <- rbind(
        cbind(GO_ID = go_all$child, EntrezID = go_all$entrezID),
        cbind(GO_ID = go_all$ancestor, EntrezID = go_all$entrezID)
    )
    
    # Remove duplicates and empty entries
    result <- unique(result)
    result <- result[!is.na(result[, "GO_ID"]) & result[, "GO_ID"] != "", , drop = FALSE]
    
    # If result is too small, return just the children
    if (nrow(result) < 3) {
        result <- unique(cbind(
            GO_ID = input_go_ids,
            EntrezID = input_entrez_ids
        ))
    }
    
    return(result)
}

