#' Add common gene identifiers to annotated peaks
#' 
#' @description 
#' Adds common gene identifiers (such as gene symbol, Entrez ID, Ensembl gene ID,
#' and RefSeq ID) to annotated peaks by leveraging organism annotation databases
#' or biomaRt. This function is useful for enriching peak annotations with
#' additional gene identifiers that may be needed for downstream analysis,
#' visualization, or integration with other datasets.
#' 
#' The function supports two annotation sources:
#' \itemize{
#'   \item \strong{Organism annotation databases} (via \code{orgAnn}): Uses
#'         AnnotationDbi-based packages like \code{org.Hs.eg.db} for human or
#'         \code{org.Mm.eg.db} for mouse. These provide fast, local annotation
#'         lookups.
#'   \item \strong{biomaRt} (via \code{mart}): Uses online Ensembl or other
#'         biomaRt databases for annotation. Requires internet connection but
#'         provides access to the latest annotations.
#' }
#' 
#' @param annotatedPeak A \link[GenomicRanges:GRanges-class]{GRanges} object
#'        with a \code{feature} metadata column containing feature IDs, or a
#'        character vector of feature IDs. If a GRanges object, the function
#'        extracts unique feature IDs from \code{annotatedPeak$feature} and
#'        adds the new ID columns to the metadata. If a character vector, the
#'        function returns a data frame mapping feature IDs to the requested
#'        identifiers.
#' @param orgAnn A character string specifying an organism annotation database
#'        package name (e.g., \code{"org.Hs.eg.db"} for human,
#'        \code{"org.Mm.eg.db"} for mouse). The package must be installed and
#'        will be loaded automatically. Alternatively, an \code{OrgDb} object
#'        can be provided, which will be converted to its package name.
#'        Required if \code{mart} is not provided.
#' @param IDs2Add A character vector specifying which annotation identifiers
#'        to add. Default is \code{c("symbol")}. The available options depend
#'        on whether \code{orgAnn} or \code{mart} is used (see Details).
#' @param feature_id_type A character string specifying the type of ID in the
#'        input feature column. Default is \code{"ensembl_gene_id"}. When using
#'        \code{orgAnn}, must be one of: \code{"ensembl_gene_id"},
#'        \code{"entrez_id"}, \code{"gene_symbol"}, \code{"gene_alias"}, or
#'        \code{"refseq_id"}. When using \code{mart}, must be a valid filter
#'        name from \code{listFilters(mart)}.
#' @param silence A logical value. If \code{TRUE} (default), suppresses
#'        messages about unmapped IDs and progress updates. If \code{FALSE},
#'        displays messages when feature IDs cannot be mapped to Entrez IDs
#'        and when adding each identifier.
#' @param mart A \code{Mart} object from the \code{biomaRt} package. Can be
#'        created using \code{\link[biomaRt:useMart]{useMart}}. Required if
#'        \code{orgAnn} is not provided.
#' 
#' @details 
#' \strong{One of \code{orgAnn} or \code{mart} must be provided.}
#' 
#' \strong{Using \code{orgAnn} (AnnotationDbi-based databases):} 
#' \itemize{
#'   \item \code{feature_id_type} must be one of: \code{"ensembl_gene_id"},
#'         \code{"entrez_id"}, \code{"gene_symbol"}, \code{"gene_alias"}, or
#'         \code{"refseq_id"}. The function first maps these IDs to Entrez IDs
#'         (if not already Entrez IDs), then uses Entrez IDs as the key to
#'         retrieve additional identifiers.
#'   \item \code{IDs2Add} can be any combination of the following identifiers:
#'         \itemize{
#'           \item \code{"accnum"}: GenBank accession numbers
#'           \item \code{"ensembl"}: Ensembl gene accession numbers
#'           \item \code{"ensemblprot"}: Ensembl protein accession numbers
#'           \item \code{"ensembltrans"}: Ensembl transcript accession numbers
#'           \item \code{"entrez_id"}: Entrez gene identifiers
#'           \item \code{"enzyme"}: EC (Enzyme Commission) numbers
#'           \item \code{"genename"}: Full gene names
#'           \item \code{"pfam"}: Pfam protein domain identifiers
#'           \item \code{"pmid"}: PubMed identifiers
#'           \item \code{"prosite"}: PROSITE protein domain identifiers
#'           \item \code{"refseq"}: RefSeq identifiers
#'           \item \code{"symbol"}: Gene symbols (abbreviations)
#'           \item \code{"unigene"}: UniGene cluster identifiers
#'           \item \code{"uniprot"}: UniProt accession numbers
#'           \item \code{"omim"}: OMIM (Online Mendelian Inheritance in Man)
#'                 identifiers (human only, for \code{org.Hs.eg.db})
#'           \item \code{"mgi"}: Jackson Laboratory MGI gene accession numbers
#'                 (mouse only, for \code{org.Mm.eg.db})
#'         }
#'   \item The function handles one-to-many mappings (e.g., one Ensembl ID
#'         mapping to multiple Entrez IDs) by condensing multiple values into
#'         semicolon-separated strings.
#'   \item If a feature ID cannot be mapped to an Entrez ID, it will be
#'         excluded from the output (unless \code{silence = FALSE}, in which
#'         case a message is displayed).
#' }
#' 
#' \strong{Using \code{mart} (biomaRt):}
#' \itemize{
#'   \item \code{feature_id_type} must be a valid filter name from
#'         \code{\link[biomaRt:listFilters]{listFilters(mart)}}, such as
#'         \code{"ensembl_gene_id"}.
#'   \item \code{IDs2Add} must be one or more valid attribute names from
#'         \code{\link[biomaRt:listAttributes]{listAttributes(mart)}}, such as
#'         \code{"hgnc_symbol"}, \code{"entrezgene"}, \code{"wikigene_name"},
#'         or \code{"mirbase_transcript_name"}.
#'   \item The function directly queries biomaRt using
#'         \code{\link[biomaRt:getBM]{getBM}} with the specified filters and
#'         attributes.
#'   \item Requires an active internet connection.
#' }
#' 
#' \strong{Input processing:}
#' \itemize{
#'   \item Feature IDs are extracted from \code{annotatedPeak$feature} if input
#'         is a GRanges object, or from the character vector itself.
#'   \item Empty strings and NA values are automatically removed.
#'   \item Only unique feature IDs are used for annotation lookup (duplicates
#'         are handled during merging back to the original input).
#' }
#' 
#' \strong{Output:}
#' \itemize{
#'   \item If input is a GRanges object: Returns the same GRanges object with
#'         additional metadata columns for each identifier in \code{IDs2Add}.
#'         The order of peaks is preserved.
#'   \item If input is a character vector: Returns a data frame with columns
#'         \code{feature_id_type} and all identifiers in \code{IDs2Add}.
#'         Rows are ordered by the input feature IDs.
#' }
#' 
#' @return 
#' \itemize{
#'   \item If \code{annotatedPeak} is a GRanges object: Returns a GRanges
#'         object with the same structure, but with additional metadata columns
#'         (one for each identifier in \code{IDs2Add}) added to
#'         \code{mcols(annotatedPeak)}. Peaks that could not be mapped will
#'         have \code{NA} values in the new columns.
#'   \item If \code{annotatedPeak} is a character vector: Returns a data frame
#'         with columns \code{feature_id_type} and all identifiers in
#'         \code{IDs2Add}. Each row represents a unique feature ID and its
#'         mapped identifiers. If multiple mappings exist (e.g., one Ensembl
#'         ID mapping to multiple Entrez IDs), values are condensed into
#'         semicolon-separated strings.
#' }
#' @references http://www.bioconductor.org/packages/release/data/annotation/
#' @author Jianhong Ou, Lihua Julie Zhu
#' @seealso \link[biomaRt:getBM]{getBM}, AnnotationDb
#' @export
#' @importFrom AnnotationDbi mget
#' @importFrom biomaRt getBM
#' @keywords misc
#' @examples
#' \dontrun{
#' ## Example 1: Using organism annotation database (org.Hs.eg.db)
#' data(annotatedPeak)
#' library(org.Hs.eg.db)
#' 
#' # Add gene symbols and OMIM IDs to annotated peaks
#' annotated_with_ids <- addGeneIDs(annotatedPeak[1:6, ],
#'                                   orgAnn = "org.Hs.eg.db",
#'                                   IDs2Add = c("symbol", "omim"))
#' 
#' # Check the added columns
#' mcols(annotated_with_ids)[, c("feature", "symbol", "omim")]
#' 
#' ## Example 2: Using character vector input (returns data frame)
#' feature_ids <- annotatedPeak$feature[1:6]
#' id_mapping <- addGeneIDs(feature_ids,
#'                         orgAnn = "org.Hs.eg.db",
#'                         IDs2Add = c("symbol", "genename"))
#' head(id_mapping)
#' 
#' ## Example 3: Using biomaRt (requires internet connection)
#' if (interactive()) {
#'     library(biomaRt)
#'     mart <- useMart("ENSEMBL_MART_ENSEMBL",
#'                     host = "www.ensembl.org",
#'                     dataset = "hsapiens_gene_ensembl")
#'     
#'     annotated_with_ids <- addGeneIDs(annotatedPeak[1:6, ],
#'                                      mart = mart,
#'                                      feature_id_type = "ensembl_gene_id",
#'                                      IDs2Add = c("hgnc_symbol", "entrezgene"))
#' }
#' 
#' ## Example 4: Adding multiple identifiers
#' annotated_with_ids <- addGeneIDs(annotatedPeak[1:10, ],
#'                                  orgAnn = "org.Hs.eg.db",
#'                                  IDs2Add = c("symbol", "genename", "refseq",
#'                                             "uniprot", "pmid"))
#' }


addGeneIDs <- function(annotatedPeak, orgAnn, IDs2Add = c("symbol"), 
                       feature_id_type = "ensembl_gene_id",
                       silence = TRUE, 
                       mart) {
    # Input validation
    if (missing(annotatedPeak)) {
        stop("Missing required argument annotatedPeak!", call. = FALSE)
    }
    
    if (missing(orgAnn) && missing(mart)) {
        stop("no annotation database selected", call. = FALSE)
    }
    
    # Extract feature IDs from input
    if (inherits(annotatedPeak, "GRanges")) {
        feature_ids <- unique(annotatedPeak$feature)
    } else if (is.character(annotatedPeak)) {
        feature_ids <- unique(annotatedPeak)
    } else {
        stop("annotatedPeak needs to be GRanges type with ",
             "feature variable holding the feature id or a ",
             "character vector holding the IDs of the features ",
             "used to annotate the peaks!", call. = FALSE)
    }
    
    # Clean feature IDs
    feature_ids <- feature_ids[!is.na(feature_ids)]
    feature_ids <- feature_ids[feature_ids != ""]
    
    if (length(feature_ids) == 0) {
        stop("There is no feature column in annotatedPeak or ",
             "annotatedPeak has size 0!", call. = FALSE)
    }
    # Process orgAnn path
    if (!missing(orgAnn)) {
        if (inherits(orgAnn, "OrgDb")) {
            orgAnn <- deparse(substitute(orgAnn))
        }
        if (!is.character(orgAnn)) {
            stop("orgAnn must be a character.", call. = FALSE)
        }
        if (!grepl(".eg.db", orgAnn, ignore.case = TRUE)) {
            stop("Annotation database must be *.eg.db", call. = FALSE)
        }
        
        # Check if package is installed and available
        if (!requireNamespace(orgAnn, quietly = TRUE)) {
            stop("Package '", orgAnn, "' is not installed. ",
                 "Please install it with: BiocManager::install('", orgAnn, "')",
                 call. = FALSE)
        }
        orgAnn <- sub("\\.db$", "", orgAnn, ignore.case = TRUE)
        # Get Entrez IDs
        if (feature_id_type == "entrez_id") {
            m_ent <- as.data.frame(feature_ids, stringsAsFactors = FALSE)
            colnames(m_ent) <- c("entrez_id")
        } else {
            prefix <- switch(feature_id_type,
                             gene_alias = "ALIAS",
                             gene_symbol = "SYMBOL",
                             ensembl_gene_id = "ENSEMBL",
                             refseq_id = "REFSEQ",
                             "UNKNOWN"
            )
            if (prefix == "UNKNOWN") {
                stop("Currently only the following type of IDs are supported: ",
                     "entrez_id, gene_alias, ensembl_gene_id, refseq_id and ",
                     "gene_symbol!", call. = FALSE)
            }
            
            # Get mapping environment
            env_name <- paste0(orgAnn, prefix, "2EG")
            tryCatch({
                env <- get(env_name)
            }, error = function(e) {
                stop("Annotation database ", env_name, " does not exist!\n",
                     "\tPlease try to load annotation database by library(",
                     orgAnn, ".db)", call. = FALSE)
            })
            
            # Map feature IDs to Entrez IDs
            entrez <- AnnotationDbi::mget(feature_ids, env, ifnotfound = NA)
            gene_ids <- names(entrez)
            m_ent <- do.call(rbind, lapply(gene_ids, function(.ele) {
                r <- entrez[[.ele]]
                if (!is.na(r[1])) {
                    cbind(rep(.ele, length(r)), r)
                } else {
                    if (!silence) {
                        message("entrez id for '", .ele, "' not found\n")
                    }
                    c(.ele, NA)
                }
            }))
            m_ent <- as.data.frame(m_ent, stringsAsFactors = FALSE)
            m_ent <- m_ent[!is.na(m_ent[, 1]), , drop = FALSE]
            colnames(m_ent) <- c(feature_id_type, "entrez_id")
        }
        
        entrezIDs <- as.character(m_ent$entrez_id)
        entrezIDs <- unique(entrezIDs)
        entrezIDs <- entrezIDs[!is.na(entrezIDs)]
        
        if (length(entrezIDs) == 0L) {
            stop("No entrez identifier can be mapped by input data based on ",
                 "the feature_id_type. Please consider to use correct ",
                 "feature_id_type, orgAnn or annotatedPeak", call. = FALSE)
        }
        # Add additional IDs
        IDs2Add <- unique(IDs2Add)
        IDs2Add <- IDs2Add[IDs2Add != feature_id_type]
        IDs <- unique(entrezIDs[!is.na(entrezIDs)])
        
        for (IDtoAdd in IDs2Add) {
            x <- NULL
            if (!silence) {
                message("Adding ", IDtoAdd, " ... ")
            }
            
            if (IDtoAdd != "entrez_id") {
                orgDB <- NULL
                db_name <- paste0(orgAnn, toupper(IDtoAdd))
                tryCatch({
                    orgDB <- get(db_name)
                }, error = function(e) {
                    if (!silence) {
                        message("The IDs2Add you input, \"", IDtoAdd, 
                                "\", is not supported!\n")
                    }
                })
                
                if (is.null(orgDB)) {
                    IDs2Add <- IDs2Add[IDs2Add != IDtoAdd]
                    next
                }
                
                if (!inherits(orgDB, "AnnDbBimap") && !inherits(orgDB, "IpiAnnDbMap")) {
                    if (!silence) {
                        message("The IDs2Add you input, \"", IDtoAdd, 
                                "\", is not supported!\n")
                    }
                    IDs2Add <- IDs2Add[IDs2Add != IDtoAdd]
                    next
                }
                
                x <- AnnotationDbi::mget(IDs, orgDB, ifnotfound = NA)
                x <- sapply(x, base::paste, collapse = ";")
                x <- as.data.frame(x, stringsAsFactors = FALSE)
                m_ent <- merge(m_ent, x, 
                               by.x = "entrez_id",
                               by.y = "row.names",
                               all.x = TRUE)
                colnames(m_ent)[length(colnames(m_ent))] <- IDtoAdd
            }
            if (!silence) {
                message("done\n")
            }
        }
        m_ent <- m_ent[, c(feature_id_type, IDs2Add), drop = FALSE]
    } else {
        # Use biomaRt
        if (missing(mart) || !inherits(mart, "Mart")) {
            stop("No valid mart object is passed in!", call. = FALSE)
        }
        
        IDs2Add <- unique(IDs2Add)
        IDs2Add <- IDs2Add[IDs2Add != feature_id_type]
        
        tryCatch({
            m_ent <- getBM(attributes = c(feature_id_type, IDs2Add),
                           filters = feature_id_type, 
                           values = feature_ids, 
                           mart = mart)
        }, error = function(e) {
            stop("Get error when calling getBM: ", e, call. = FALSE)
        })
        
        if (any(colnames(m_ent) != c(feature_id_type, IDs2Add))) {
            colnames(m_ent) <- c(feature_id_type, IDs2Add)
        }
    }
    
    # Prepare output
    if (!silence) {
        message("prepare output ... ")
    }
    
    # Handle multiple entrez_id for single feature_id
    if (ncol(m_ent) == 1L) {
        stop("None of IDs could be appended. Please double check IDs2Add.",
             call. = FALSE)
    }
    
    duplicated_ids <- m_ent[duplicated(m_ent[, feature_id_type]), 
                            feature_id_type]
    if (length(duplicated_ids) > 0L) {
        m_ent.duplicated <- m_ent[m_ent[, feature_id_type] %in% duplicated_ids, ]
        m_ent.duplicated <- condenseMatrixByColnames(as.matrix(m_ent.duplicated),
                                                      feature_id_type)
        m_ent <- m_ent[!(m_ent[, feature_id_type] %in% duplicated_ids), ]
        m_ent <- rbind(m_ent, m_ent.duplicated)
    }
    
    # Merge back to original input
    if (inherits(annotatedPeak, "GRanges")) {
        # Rearrange m_ent by annotatedPeak$feature
        # data.frame is very important for order...
        orderlist <- data.frame(annotatedPeak$feature, stringsAsFactors = FALSE)
        orderlist <- cbind(seq_along(annotatedPeak), orderlist)
        colnames(orderlist) <- c("orderid___", feature_id_type)
        m_ent <- merge(orderlist, m_ent, by = feature_id_type, all.x = TRUE)
        m_ent <- m_ent[order(m_ent[, "orderid___"]), 
                       c(feature_id_type, IDs2Add)]
        for (IDtoAdd in IDs2Add) {
            mcols(annotatedPeak)[, IDtoAdd] <- m_ent[, IDtoAdd]
        }
    } else {
        annotatedPeak <- m_ent
    }
    
    if (!silence) {
        message("done\n")
    }
    
    return(annotatedPeak)
}
