#' Test enrichment among peak-associated features using MSigDB gene sets
#'
#' @description
#' Performs enrichment analysis on peak-associated features using gene sets from
#' MSigDB. The function uses the \code{fgsea} package's \code{fora} (Fisher's
#' Over-Representation Analysis) function to test for enrichment of gene sets
#' in the annotated peaks.
#'
#' The test is performed for each collection or subcollection of terms from
#' MSigDB. Gene sets with adjusted p-value not greater than \code{max_adjP}
#' are returned in the results.
#'
#' @param annotated_peak A \code{GRanges} object containing annotated peaks
#'   with gene identifiers in one of the metadata columns. The peaks should
#'   have been annotated using functions like \code{annotatePeakInBatch()} or
#'   \code{annoPeaks()}.
#' @param id_column Character string specifying the name of the column in
#'   \code{annotated_peak} that contains gene identifiers. Default is
#'   \code{"gene_symbol"}.
#' @param id_type Character string specifying the type of gene identifier in
#'   \code{annotated_peak[[id_column]]}. Options:
#'   \itemize{
#'     \item \code{"gene_symbol"} (default): Gene symbols (e.g., "TP53", "BRCA1")
#'     \item \code{"ensembl_gene"}: Ensembl gene IDs (version numbers are
#'           automatically removed)
#'     \item \code{"ncbi_gene"}: NCBI Gene IDs (Entrez IDs)
#'   }
#'   This must match the \code{id_type} used when loading terms with
#'   \code{load_terms_from_msigdbr()}.
#' @param collections A list of gene set collections as returned by
#'   \code{load_terms_from_msigdbr()}. Each element should be a named list
#'   where names are gene set names and values are character vectors of gene IDs.
#'   Required parameter.
#' @param universe Optional character vector of all gene IDs to use as the
#'   background/universe for enrichment testing. If \code{NULL}, all unique
#'   genes in \code{collections} will be used as the universe (with a warning).
#'   All IDs in \code{annotated_peak[[id_column]]} must be present in the
#'   universe. For testing enrichment among all peak-associated features,
#'   the universe should be the set of all genes in the genome. For testing
#'   enrichment among features associated with differential peaks,
#'   the universe should be all features associated with all peaks in the 
#'   dataset.
#' @param max_adjP Numeric value specifying the maximum adjusted p-value
#'   (FDR) threshold for significance. Gene sets with adjusted p-value greater
#'   than this threshold are filtered out. Default is \code{0.1}.
#' @param min_term_size Integer specifying the minimum number of genes a gene
#'   set must have to be included in the analysis. Gene sets smaller than this
#'   are excluded. Default is \code{5}.
#' @param max_term_size Integer or \code{Inf} specifying the maximum number of
#'   genes a gene set can have to be included in the analysis. Gene sets larger
#'   than this are excluded. Default is \code{Inf} (no upper limit).
#'
#' @return A list of data frames, one for each collection/subcollection in
#'   \code{collections}. Each data frame contains enrichment results with
#'   columns:
#'   \itemize{
#'     \item \code{pathway}: Gene set name
#'     \item \code{pval}: Raw p-value from Fisher's exact test
#'     \item \code{padj}: Adjusted p-value (FDR)
#'     \item \code{log2err}: Log2 of the standard error of the log2 fold change
#'     \item \code{ES}: Enrichment score
#'     \item \code{NES}: Normalized enrichment score
#'     \item \code{size}: Number of genes in the gene set
#'     \item \code{leadingEdge}: Character vector of leading edge genes
#'   }
#'   Only gene sets with \code{padj <= max_adjP} are included in each data frame.
#'
#' @details
#' \strong{Statistical Method:}
#' The function uses Fisher's Over-Representation Analysis (FORA) from the
#' \code{fgsea} package, which is a hypergeometric test comparing the overlap
#' between the query gene set and the gene set of interest against the
#' background universe.
#'
#' \strong{Input Requirements:}
#' \itemize{
#'   \item \code{annotated_peak} must be a \code{GRanges} object
#'   \item \code{id_column} must exist in \code{colnames(annotated_peak)}
#'   \item All IDs in \code{annotated_peak[[id_column]]} must be present in
#'         \code{universe}
#'   \item \code{id_type} must match the type used in \code{collections}
#' }
#'
#' \strong{ID Type Matching:}
#' \itemize{
#'   \item If \code{id_type == "ensembl_gene"}, version numbers (e.g., ".1",
#'         ".2") are automatically removed from Ensembl IDs
#'   \item Ensure the \code{id_type} parameter matches the \code{id_type} used
#'         when calling \code{load_terms_from_msigdbr()}
#' }
#'
#' @importFrom fgsea fora
#' @export
#'
#' @examples
#' \dontrun{
#' # Load MSigDB gene sets
#' collections <- load_terms_from_msigdbr(
#'   species = "Homo sapiens",
#'   collection = "H",
#'   id_type = "gene_symbol"
#' )
#' library(org.Hs.eg.db)
#' all_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")  # For gene symbols
#' # Perform enrichment analysis
#' # Assuming 'annotated_peaks' is a GRanges object with a 'gene_symbol' column
#' results <- test_enrichment(
#'   annotated_peak = annotated_peaks,
#'   id_column = "gene_symbol",
#'   id_type = "gene_symbol",
#'   universe = all_genes,  # character vector of all gene symbols
#'   collections = collections,
#'   max_adjP = 0.05,
#'   min_term_size = 10
#' )
#'
#' # Access results for a specific collection
#' hallmark_results <- results[[1]]
#' }
test_enrichment <- function(annotated_peak,
                            id_column = "gene_symbol",
                            id_type = c("gene_symbol", "ensembl_gene", "ncbi_gene"),
                            collections,
                            universe = NULL,
                            max_adjP = 0.1,
                            min_term_size = 5,
                            max_term_size = Inf) {
    id_type <- match.arg(id_type)
    
    if (missing(collections)) {
        stop("Missing required argument 'collections'! ",
             "Use load_terms_from_msigdbr() to create collections.",
             call. = FALSE)
    }
    if (missing(annotated_peak)) {
        stop("Missing required argument 'annotated_peak'!", call. = FALSE)
    }
    if (!inherits(annotated_peak, "GRanges")) {
        stop("'annotated_peak' must be a GRanges object", call. = FALSE)
    }
    
    if (!id_column %in% colnames(annotated_peak)) {
        stop("'id_column' ('", id_column, "') must be a column in 'annotated_peak' metadata. ",
             "Available columns: ", paste(colnames(mcols(annotated_peak)), collapse = ", "),
             call. = FALSE)
    }
    
    ids <- unique(annotated_peak[[id_column]])
    ids <- ids[!is.na(ids)]
    
    if (length(ids) == 0) {
        stop("No valid gene IDs found in 'annotated_peak[[\"", id_column, "\"]]'", call. = FALSE)
    }
    
    if (id_type == "ensembl_gene") {
        ids <- sub("\\.[0-9]+$", "", ids)
    }
    
    if (is.null(universe)) {
        stop("Missing 'universe' argument. ",
             "Please provide a character vector of all gene IDs to use as the ",
             "background/universe for enrichment testing.",
             call. = FALSE)
    }
    
    if (!all(ids %in% universe)) {
        missing_ids <- setdiff(ids, universe)
        stop("Not all IDs in 'annotated_peak[[id_column]]' are present in 'universe'. ",
             "Missing IDs (first 10): ", paste(head(missing_ids, 10), collapse = ", "),
             call. = FALSE)
    }
    
    # Test enrichment for each collection
    results <- lapply(collections, function(collection) {
        collection <- lapply(collection, function(x) tolower(x))
        res <- fgsea::fora(pathways = collection,
                           genes = tolower(ids),
                           universe = tolower(universe),
                           minSize = min_term_size,
                           maxSize = max_term_size)
        res <- res[res$padj <= max_adjP, ]
        res
    })
    
    return(results)
}


#' List available msigdbr species
#'
#' @description
#' Retrieves a data frame listing all species supported by the \code{msigdbr}
#' package. The \code{msigdbr} package provides access to the Molecular
#' Signatures Database (MSigDB) gene sets for multiple species.
#'
#' @details
#' MSigDB contains curated gene sets including:
#' \itemize{
#'   \item Hallmark gene sets
#'   \item Positional gene sets
#'   \item Curated gene sets
#'   \item Motif gene sets
#'   \item Computational gene sets
#'   \item GO gene sets
#'   \item Oncogenic signatures
#'   \item Immunologic signatures
#' }
#'
#' Use this function to check which species are available before calling
#' \code{load_terms_from_msigdbr()}.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{species_name}: Latin name of the species (e.g., "Homo sapiens")
#'   \item \code{species_common_name}: Common name of the species (e.g., "human")
#'   \item Additional metadata columns as provided by \code{msigdbr}
#' }
#'
#' @importFrom msigdbr msigdbr_species
#' @export
#'
#' @examples
#' \dontrun{
#' # List all available species
#' species_list <- list_available_msigdbr_species()
#' print(species_list)
#' }
list_available_msigdbr_species <- function() {
    msigdbr::msigdbr_species()
}

#' Load gene set terms from MSigDB via msigdbr
#'
#' @description
#' Loads gene set terms from the Molecular Signatures Database (MSigDB) using
#' the \code{msigdbr} package. MSigDB provides curated gene sets for over 20
#' species, including hallmark gene sets, positional gene sets, curated gene
#' sets, motif gene sets, computational gene sets, GO gene sets, oncogenic
#' signatures, and immunologic signatures.
#'
#' The function retrieves gene sets from MSigDB and organizes them into a
#' nested list structure where outer keys are collection/subcollection names
#' and inner lists map gene IDs to vectors of gene set names.
#'
#' @param species Character string specifying the species Latin name.
#'   Default is \code{"Homo sapiens"}. Use \code{list_available_msigdbr_species()}
#'   to see all available species.
#' @param collection Optional character string specifying the MSigDB collection
#'   to retrieve. If \code{NULL} (default), all collections are retrieved.
#'   Common collections include: \code{"H"}, \code{"C1"}, \code{"C2"}, \code{"C3"},
#'   \code{"C4"}, \code{"C5"}, \code{"C6"}, \code{"C7"}, \code{"C8"}.
#' @param subcollection Optional character string specifying the MSigDB
#'   subcollection to retrieve. If \code{NULL} (default), all subcollections
#'   are retrieved. Subcollections vary by collection (e.g., \code{"CP:REACTOME"},
#'   \code{"CP:BIOCARTA"} for C2).
#' @param id_type Character string specifying the type of gene identifier to
#'   use. Options:
#'   \itemize{
#'     \item \code{"ensembl_gene"} (default): Ensembl gene IDs
#'     \item \code{"ncbi_gene"}: NCBI Gene IDs (Entrez IDs)
#'     \item \code{"gene_symbol"}: Gene symbols (e.g., "TP53", "BRCA1")
#'   }
#'
#' @return A nested list where:
#' \itemize{
#'   \item Outer keys are collection/subcollection names (e.g., "H", "C2_CP:REACTOME")
#'   \item Inner lists are named vectors where:
#'     \itemize{
#'       \item Names are gene set names (e.g., "HALLMARK_APOPTOSIS")
#'       \item Values are character vectors of gene IDs of the specified type
#'     }
#' }
#'
#' @details
#' \strong{Collection Structure:}
#' \itemize{
#'   \item Gene sets are organized by collection and subcollection
#'   \item If a subcollection is empty, the collection name is used as the key
#'   \item If a subcollection exists, the key format is "COLLECTION_SUBCOLLECTION"
#'   \item Special characters (e.g., ":") in names are replaced with "-"
#' }
#'
#' \strong{Species Support:}
#' The \code{msigdbr} package supports over 20 species. Common species include:
#' \itemize{
#'   \item \code{"Homo sapiens"} (human)
#'   \item \code{"Mus musculus"} (mouse)
#'   \item \code{"Rattus norvegicus"} (rat)
#'   \item \code{"Drosophila melanogaster"} (fly)
#'   \item \code{"Caenorhabditis elegans"} (worm)
#'   \item \code{"Saccharomyces cerevisiae"} (yeast)
#'   \item \code{"Danio rerio"} (zebrafish)
#' }
#'
#' @importFrom msigdbr msigdbr msigdbr_species
#' @export
#'
#' @examples
#' \dontrun{
#' # Load all gene sets for human using Ensembl gene IDs
#' human_terms <- load_terms_from_msigdbr(species = "Homo sapiens")
#'
#' # Load only Hallmark gene sets
#' hallmark_terms <- load_terms_from_msigdbr(
#'   species = "Homo sapiens",
#'   collection = "H"
#' )
#'
#' # Load Reactome pathways using gene symbols
#' reactome_terms <- load_terms_from_msigdbr(
#'   species = "Homo sapiens",
#'   collection = "C2",
#'   subcollection = "CP:REACTOME",
#'   id_type = "gene_symbol"
#' )
#' }
load_msigdbr_gene_sets <- function(species = "Homo sapiens",
                                    collection = NULL,
                                    subcollection = NULL,
                                    id_type = c("ensembl_gene", "ncbi_gene", "gene_symbol")) {
    id_type <- match.arg(id_type)
    
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("msigdbr package is required. Please install it using: install.packages('msigdbr')")
    }

    supported_species <- msigdbr::msigdbr_species()
    if (!species %in% supported_species$species_name) {
        message("Available species in msigdbr:")
        print(supported_species)
        stop(paste0("Species '", species, "' not found in msigdbr. ",
                    "Use list_available_msigdbr_species() to see available species."),
             call. = FALSE)
    }
    
    all_gene_sets <- msigdbr::msigdbr(species = species,
                                      collection = collection,
                                      subcollection = subcollection)

    msigdbr_list <- split(x = all_gene_sets,
                          f = ifelse(all_gene_sets$gs_subcollection == "",
                                     all_gene_sets$gs_collection,
                                     paste(all_gene_sets$gs_collection,
                                           all_gene_sets$gs_subcollection,
                                           sep = "_")))

    terms <- lapply(msigdbr_list, function(.x) {
        split(.x[[id_type]], f = .x$gs_name)
    })
    names(terms) <- gsub(":", "-", names(terms))
    terms
}

#' Load gene set libraries from EnrichR database
#'
#' @description
#' Loads gene set libraries from the EnrichR database
#' (\url{https://maayanlab.cloud/Enrichr/}). EnrichR is a comprehensive gene
#' set enrichment analysis web server that provides access to over 100 gene
#' set libraries covering various biological processes, pathways, and
#' ontologies.
#'
#' This function downloads gene set libraries from EnrichR and organizes them
#' into a nested list structure where outer keys are library names and inner
#' lists map gene symbols to vectors of term names.
#'
#' @param species Character string specifying the species. Options:
#'   \itemize{
#'     \item \code{"Homo sapiens"} (default): Human
#'     \item \code{"Mus musculus"}: Mouse
#'     \item \code{"Drosophila melanogaster"}: Fruit fly
#'     \item \code{"Caenorhabditis elegans"}: Worm
#'     \item \code{"Saccharomyces cerevisiae"}: Yeast
#'     \item \code{"Danio rerio"}: Zebrafish
#'   }
#' @param sleep_time Numeric value specifying the sleep time (in seconds)
#'   between API requests to EnrichR. This helps avoid rate limiting.
#'   Default is \code{2} seconds.
#' @param geneset_library Optional character vector specifying the names of
#'   gene set libraries to download. Use \code{get_enrichr_libraries_info()}
#'   to see available libraries. If \code{NULL} (default), a warning is issued
#'   and all libraries will be downloaded (this can take a long time).
#'
#' @return A nested list where:
#' \itemize{
#'   \item Outer keys are library names (e.g., "GO_Biological_Process_2021")
#'   \item Inner lists are named vectors where:
#'     \itemize{
#'       \item Names are gene symbols
#'       \item Values are character vectors of term names that contain that gene
#'     }
#' }
#'
#' @details
#' \strong{EnrichR Database:}
#' EnrichR provides access to over 100 gene set libraries including:
#' \itemize{
#'   \item Gene Ontology (GO) terms (biological process, molecular function,
#'         cellular component)
#'   \item Pathway databases (KEGG, Reactome, WikiPathways)
#'   \item Transcription factor targets
#'   \item Disease associations
#'   \item Drug signatures
#'   \item Cell type signatures
#'   \item And many more
#' }
#'
#' \strong{Species Support:}
#' EnrichR supports different species through different web servers:
#' \itemize{
#'   \item \code{"Enrichr"}: Human and mouse
#'   \item \code{"FlyEnrichr"}: Drosophila melanogaster
#'   \item \code{"WormEnrichr"}: Caenorhabditis elegans
#'   \item \code{"YeastEnrichr"}: Saccharomyces cerevisiae
#'   \item \code{"FishEnrichr"}: Danio rerio
#' }
#'
#' \strong{Performance Notes:}
#' \itemize{
#'   \item Downloading all libraries can take a very long time (hours)
#'   \item It is recommended to specify only the libraries of interest using
#'         \code{geneset_library}
#'   \item Use \code{get_enrichr_libraries_info()} to see available libraries
#'         before downloading
#'   \item The function includes a progress bar when downloading multiple
#'         libraries
#' }
#'
#' @importFrom enrichR listEnrichrSites setEnrichrSite listEnrichrDbs
#' @importFrom progress progress_bar
#' @export
#'
#' @examples
#' \dontrun{
#' # Get information about available libraries first
#' library_info <- get_enrichr_libraries_info(species = "Homo sapiens")
#' print(head(library_info))
#'
#' # Download specific libraries
#' libraries <- get_enrichr_libraries(
#'   species = "Homo sapiens",
#'   geneset_library = c("GO_Biological_Process_2021", "KEGG_2021_Human"),
#'   sleep_time = 2
#' )
#'
#' # Access a specific library
#' go_terms <- libraries[["GO_Biological_Process_2021"]]
#' }
get_enrichr_libraries <- function(species = c("Homo sapiens",
                                              "Mus musculus",
                                              "Drosophila melanogaster",
                                              "Caenorhabditis elegans",
                                              "Saccharomyces cerevisiae",
                                              "Danio rerio"),
                                  sleep_time = 2,
                                  geneset_library = NULL) {
    species <- match.arg(species)
    enrichr_libraries <- get_enrichr_libraries_info(species = species)
    
    if (is.null(geneset_library)) {
        warning("Missing 'geneset_library' argument. ",
                "All EnrichR libraries will be downloaded. This may take a very long time. ",
                "Consider using get_enrichr_libraries_info() to select specific libraries.",
                call. = FALSE)
        print(enrichr_libraries)
        geneset_library <- enrichr_libraries$libraryName
    } else if (!all(geneset_library %in% enrichr_libraries$libraryName)) {
        missing <- setdiff(geneset_library, enrichr_libraries$libraryName)
        stop("Some geneset libraries are not available: ",
             paste(missing, collapse = ", "), ". ",
             "Use get_enrichr_libraries_info() to see available libraries.",
             call. = FALSE)
    } 
    pb <- progress::progress_bar$new(
        format = "[:bar] :percent ETA: :eta",
        total = length(geneset_library),
        clear = FALSE,
        width = 60,
        show_after = 0
    )
    libraries <- lapply(geneset_library, function(.library) {
        Sys.sleep(sleep_time)
        terms <- enrichR:::.read_gmt(.library)
        if (!is.null(pb)) pb$tick()
        terms
    })
    names(libraries) <- geneset_library
    if (!is.null(pb)) pb$terminate()
    message("Downloaded ", length(libraries), " libraries from EnrichR")
    return(libraries)
}

#' Get information about available EnrichR libraries for a species
#'
#' @description
#' Retrieves a data frame listing all available gene set libraries in EnrichR
#' for a given species. This function is useful for exploring available
#' libraries before downloading them with \code{get_enrichr_libraries()}.
#'
#' The function automatically sets the appropriate EnrichR web server based on
#' the species and retrieves the list of available databases.
#'
#' @param species Character string specifying the species. Options:
#'   \itemize{
#'     \item \code{"Homo sapiens"} (default): Human
#'     \item \code{"Mus musculus"}: Mouse
#'     \item \code{"Drosophila melanogaster"}: Fruit fly
#'     \item \code{"Caenorhabditis elegans"}: Worm
#'     \item \code{"Saccharomyces cerevisiae"}: Yeast
#'     \item \code{"Danio rerio"}: Zebrafish
#'   }
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{libraryName}: Name of the gene set library (e.g., "GO_Biological_Process_2021")
#'   \item \code{description}: Description of the library
#'   \item Additional metadata columns as provided by EnrichR
#' }
#'
#' @details
#' \strong{EnrichR Web Servers:}
#' EnrichR uses different web servers for different species:
#' \itemize{
#'   \item \code{"Enrichr"}: Default server for human and mouse
#'   \item \code{"FlyEnrichr"}: For Drosophila melanogaster
#'   \item \code{"WormEnrichr"}: For Caenorhabditis elegans
#'   \item \code{"YeastEnrichr"}: For Saccharomyces cerevisiae
#'   \item \code{"FishEnrichr"}: For Danio rerio
#' }
#'
#' The function automatically maps the species to the appropriate server and
#' sets it before retrieving the library list.
#'
#' \strong{Internet Connection:}
#' This function requires an active internet connection and access to the
#' EnrichR web servers. If the website is not accessible, the function will
#' stop with an error message.
#'
#' @importFrom enrichR listEnrichrSites setEnrichrSite listEnrichrDbs
#' @export
#'
#' @examples
#' \dontrun{
#' # Get available libraries for human
#' human_libs <- get_enrichr_libraries_info(species = "Homo sapiens")
#' print(head(human_libs))
#'
#' # Get available libraries for mouse
#' mouse_libs <- get_enrichr_libraries_info(species = "Mus musculus")
#' print(head(mouse_libs))
#'
#' # Search for specific libraries
#' go_libs <- human_libs[grepl("GO", human_libs$libraryName), ]
#' print(go_libs)
#' }
get_enrichr_libraries_info <- function(species = c("Homo sapiens",
                                                   "Mus musculus",
                                                   "Drosophila melanogaster",
                                                   "Caenorhabditis elegans",
                                                   "Saccharomyces cerevisiae",
                                                   "Danio rerio")) {
    if (!requireNamespace("enrichR", quietly = TRUE)) {
        stop("enrichR package is required. Please install it using: ",
             "install.packages('enrichR')", call. = FALSE)
    }
    
    library(enrichR)
    website_live <- getOption("enrichR.live")
    if (!website_live) {
        stop("EnrichR website is not accessible. ",
             "Please check your internet connection and try again.",
             call. = FALSE)
    }

    species <- match.arg(species)

    # Map species to EnrichR site names
    # EnrichR uses different sites for different species:
    # - "Enrichr" (default) for human and mouse
    # - "FlyEnrichr" for Drosophila
    # - "WormEnrichr" for C. elegans
    # - "YeastEnrichr" for S. cerevisiae
    # - "FishEnrichr" for zebrafish
    site_name <- switch(species,
        "Homo sapiens" = "Enrichr",
        "Mus musculus" = "Enrichr",
        "Drosophila melanogaster" = "FlyEnrichr",
        "Caenorhabditis elegans" = "WormEnrichr",
        "Saccharomyces cerevisiae" = "YeastEnrichr",
        "Danio rerio" = "FishEnrichr",
        stop(paste0("Species '", species, "' not supported by EnrichR"),
             call. = FALSE)
    )

    # Set the EnrichR site for the species
    enrichR::setEnrichrSite(site_name)

    # Get available databases for this species
    dbs <- enrichR::listEnrichrDbs()
    dbs
}