#' Add GO IDs of the ancestors for a given vector of GO IDs
#' 
#' @description 
#' Adds GO IDs of the ancestors for a given vector of GO IDs leveraging GO.db.
#' This function expands a set of GO terms by including their ancestor terms
#' from the Gene Ontology hierarchy.
#' 
#' @param go.ids A matrix with at least 4 columns where:
#'   \itemize{
#'     \item Column 1: GO IDs (character)
#'     \item Column 2: Evidence codes (not used in this function)
#'     \item Column 3: Ontology type (not used in this function)
#'     \item Column 4: Entrez IDs (character)
#'   }
#' @param ontology Character string specifying the ontology type:
#'   \itemize{
#'     \item \code{"bp"}: Biological process
#'     \item \code{"cc"}: Cellular component
#'     \item \code{"mf"}: Molecular function
#'   }
#' 
#' @return A character matrix with 2 columns:
#'   \itemize{
#'     \item Column 1: GO IDs (including ancestors)
#'     \item Column 2: Entrez IDs
#'   }
#'   The matrix contains unique rows with non-empty GO IDs.
#' 
#' @export
#' @author Lihua Julie Zhu, Haibo Liu
#' @keywords misc
#' @importFrom methods is
#' @importFrom GO.db GOBPANCESTOR GOCCANCESTOR GOMFANCESTOR
#' @examples 
#' # Add ancestors for biological process ontology
#' if (requireNamespace("GO.db", quietly = TRUE)) {
#'     result <- addAncestors(go.ids, ontology = "bp")
#'     head(result)
#' }
#' 
#' addAncestors(go.ids, ontology="bp")
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
    
    if (!is(go.ids, "matrix") || ncol(go.ids) < 4) {
        stop("'go.ids' must be a matrix with at least 4 columns.\n",
             "  Column 1: GO IDs\n",
             "  Column 2: Evidence codes\n",
             "  Column 3: Ontology type\n",
             "  Column 4: Entrez IDs")
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

