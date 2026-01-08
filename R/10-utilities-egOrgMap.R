#' Map between organism annotation package names and species names
#' 
#' @description 
#' Converts between the name of a Bioconductor organism annotation package
#' (OrgDb) and its corresponding species name. The function provides
#' bidirectional mapping: it can convert from package names to species names,
#' or from species names to package names. If an exact match is not found, it
#' uses fuzzy matching based on string distance to suggest the closest match.
#' 
#' This function is useful when you need to look up the correct package name
#' for a species, or when you have a package name and need the full species
#' name for display or further processing.
#' 
#' @param name A character string (length 1) specifying either:
#'   \itemize{
#'     \item An organism annotation package name (e.g., \code{"org.Hs.eg.db"},
#'           \code{"org.Mm.eg.db"})
#'     \item A species name (e.g., \code{"Homo sapiens"}, \code{"Mus musculus"},
#'           \code{"Drosophila melanogaster"})
#'   }
#'   The function is case-sensitive and requires exact matches (or close
#'   matches for fuzzy suggestions).
#' 
#' @return Returns a character string containing:
#'   \itemize{
#'     \item The corresponding species name (if a package name was provided)
#'     \item The corresponding package name (if a species name was provided)
#'   }
#'   If no exact match is found, the function stops with an error message that
#'   suggests the closest match based on string distance.
#' 
#' @details
#' 
#' \strong{Supported organisms:}
#' The function supports 20 common model organisms with their corresponding
#' Bioconductor annotation packages:
#' \itemize{
#'   \item \code{org.Ag.eg.db} ↔ \code{Anopheles gambiae}
#'   \item \code{org.At.eg.db} ↔ \code{Arabidopsis thaliana}
#'   \item \code{org.Bt.eg.db} ↔ \code{Bos taurus}
#'   \item \code{org.Ce.eg.db} ↔ \code{Caenorhabditis elegans}
#'   \item \code{org.Cf.eg.db} ↔ \code{Canis lupus familiaris}
#'   \item \code{org.Dm.eg.db} ↔ \code{Drosophila melanogaster}
#'   \item \code{org.Dr.eg.db} ↔ \code{Danio rerio}
#'   \item \code{org.EcK12.eg.db} ↔ \code{Escherichia coli str. K12a}
#'   \item \code{org.EcSakai.eg.db} ↔ \code{Escherichia coli O157:H7 str. Sakai}
#'   \item \code{org.Gg.eg.db} ↔ \code{Gallus gallus}
#'   \item \code{org.Hs.eg.db} ↔ \code{Homo sapiens}
#'   \item \code{org.Mm.eg.db} ↔ \code{Mus musculus}
#'   \item \code{org.Mmu.eg.db} ↔ \code{Macaca mulatta}
#'   \item \code{org.Pf.plasmo.db} ↔ \code{Plasmodium falciparum}
#'   \item \code{org.Pt.eg.db} ↔ \code{Pan troglodytes}
#'   \item \code{org.Rn.eg.db} ↔ \code{Rattus norvegicus}
#'   \item \code{org.Sc.sgd.db} ↔ \code{Saccharomyces cerevisiae}
#'   \item \code{org.Sco.eg.db} ↔ \code{Streptomyces coelicolor}
#'   \item \code{org.Ss.eg.db} ↔ \code{Sus scrofa}
#'   \item \code{org.Tgondii.eg.db} ↔ \code{Toxoplasma gondii}
#'   \item \code{org.Xl.eg.db} ↔ \code{Xenopus laevis}
#' }
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates that \code{name} is a single character string
#'   \item Checks for exact match in package names (forward mapping)
#'   \item Checks for exact match in species names (reverse mapping)
#'   \item If no exact match, calculates string distance using
#'         \code{\link[utils]{adist}} between the input and all known package
#'         and species names
#'   \item Finds the closest match (minimum string distance)
#'   \item Stops with an error message suggesting the closest match
#' }
#' 
#' \strong{Bidirectional mapping:}
#' The function works in both directions:
#' \itemize{
#'   \item \strong{Package → Species}: If you provide a package name like
#'         \code{"org.Hs.eg.db"}, it returns \code{"Homo sapiens"}
#'   \item \strong{Species → Package}: If you provide a species name like
#'         \code{"Mus musculus"}, it returns \code{"org.Mm.eg.db"}
#' }
#' 
#' \strong{Fuzzy matching:}
#' If an exact match is not found, the function uses the Levenshtein distance
#' (edit distance) to find the closest match. This is useful for:
#' \itemize{
#'   \item Handling typos or misspellings
#'   \item Providing helpful error messages
#'   \item Suggesting the correct name when close matches exist
#' }
#' The function calculates the string distance between the input and all known
#' package and species names, then suggests the one with the minimum distance.
#' 
#' \strong{Error handling:}
#' The function performs validation and provides informative error messages:
#' \itemize{
#'   \item If \code{name} is not a character string or has length != 1, stops
#'         with a validation error
#'   \item If no exact match is found, stops with an error that includes a
#'         suggestion for the closest match
#' }
#' 
#' @note
#' \itemize{
#'   \item The function is case-sensitive; \code{"Homo sapiens"} is different
#'         from \code{"homo sapiens"}
#'   \item Exact matches are required; partial matches are not supported
#'         (use fuzzy matching suggestions as a guide)
#'   \item The function only supports organisms that have Bioconductor
#'         annotation packages
#'   \item The mapping is hardcoded and includes 20 common model organisms
#' }
#' 
#' @seealso
#' \itemize{
#'   \item Organism annotation packages:
#'         \code{\link[AnnotationDbi]{AnnotationDbi}}
#'   \item String distance calculation: \code{\link[utils]{adist}}
#'   \item For a complete list of available organism packages:
#'         \url{https://bioconductor.org/packages/release/data/annotation/}
#' }
#' 
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom utils adist
#' 
#' @examples
#' 
#' # Example 1: Convert package name to species name (human)
#' egOrgMap("org.Hs.eg.db")
#' 
#' # Example 2: Convert package name to species name (mouse)
#' egOrgMap("org.Mm.eg.db")
#' 
#' # Example 3: Convert species name to package name (mouse)
#' egOrgMap("Mus musculus")
#' 
#' # Example 4: Convert species name to package name (human)
#' egOrgMap("Homo sapiens")
#' 
#' # Example 5: Convert other organisms
#' egOrgMap("org.Dm.eg.db")  # Drosophila melanogaster
#' egOrgMap("org.Ce.eg.db")  # Caenorhabditis elegans
#' egOrgMap("Danio rerio")   # org.Dr.eg.db
#' 
#' # Example 6: Fuzzy matching for misspellings (shows error with suggestion)
#' \dontrun{
#' egOrgMap("Homo sapien")      # Suggests "Homo sapiens"
#' egOrgMap("org.Hs.db")       # Suggests "org.Hs.eg.db"
#' egOrgMap("Mus muscullus")   # Suggests "Mus musculus"
#' }
#' 
#' # Example 7: Using in a workflow to get package name from species
#' species <- "Homo sapiens"
#' pkg_name <- egOrgMap(species)
#' # Then use pkg_name to load the package: library(pkg_name, character.only = TRUE)
egOrgMap <- function(name) {
    if (!is.character(name) || length(name) != 1L) {
        stop("'name' must be a single character string", call. = FALSE)
    }
    
    organism <- c(
        "org.Ag.eg.db" = "Anopheles gambiae",
        "org.At.eg.db" = "Arabidopsis thaliana",
        "org.Bt.eg.db" = "Bos taurus",
        "org.Ce.eg.db" = "Caenorhabditis elegans",
        "org.Cf.eg.db" = "Canis lupus familiaris",
        "org.Dm.eg.db" = "Drosophila melanogaster",
        "org.Dr.eg.db" = "Danio rerio",
        "org.EcK12.eg.db" = "Escherichia coli str. K12a",
        "org.EcSakai.eg.db" = "Escherichia coli O157:H7 str. Sakai",
        "org.Gg.eg.db" = "Gallus gallus",
        "org.Hs.eg.db" = "Homo sapiens",
        "org.Mm.eg.db" = "Mus musculus",
        "org.Mmu.eg.db" = "Macaca mulatta",
        "org.Pf.plasmo.db" = "Plasmodium falciparum",
        "org.Pt.eg.db" = "Pan troglodytes",
        "org.Rn.eg.db" = "Rattus norvegicus",
        "org.Sc.sgd.db" = "Saccharomyces cerevisiae",
        "org.Sco.eg.db" = "Streptomyces coelicolor",
        "org.Ss.eg.db" = "Sus scrofa",
        "org.Tgondii.eg.db" = "Toxoplasma gondii",
        "org.Xl.eg.db" = "Xenopus laevis"
    )
    
    # Exact match: package name -> species name
    if (name %in% names(organism)) {
        return(organism[name])
    }
    
    # Exact match: species name -> package name
    if (name %in% organism) {
        return(names(organism)[organism == name])
    }
    
    # Fuzzy matching: calculate string distance and suggest closest match
    org <- as.character(c(names(organism), organism))
    dis <- adist(name, org)
    closest_match <- org[dis == min(dis)]
    
    stop(
        "No exact match found for \"", name, "\". ",
        "Did you mean \"", closest_match, "\"?",
        call. = FALSE
    )
}
