#' Add common IDs to annotated peaks such as gene symbol, entrez ID, 
#' ensembl gene id and refseq id.
#' 
#' @description Add common IDs to annotated peaks such as gene symbol, 
#' entrez ID, ensembl gene id and refseq id leveraging organism annotation 
#' dataset. For example, org.Hs.eg.db is the dataset from orgs.Hs.eg.db 
#' package for human, while org.Mm.eg.db is the dataset from the org.Mm.eg.db
#' package for mouse.
#' 
#' @param annotatedPeak GRanges or a vector of feature IDs.
#' @param orgAnn organism annotation dataset such as org.Hs.eg.db. Can be a 
#' character string (e.g., "org.Hs.eg.db") or an OrgDb object (which will be 
#' converted to character).
#' @param IDs2Add a vector of annotation identifiers to be added
#' @param feature_id_type type of ID to be annotated, default is 
#' ensembl_gene_id
#' @param silence TRUE or FALSE. If TRUE, will not show unmapped entrez id
#' for feature ids.
#' @param mart mart object, see \link[biomaRt:useMart]{useMart} of biomaRt
#' package for details
#' @details One of orgAnn and mart should be assigned. 
#' \itemize{    
#'   \item If orgAnn is given, parameter feature_id_type should be 
#'   ensembl_gene_id, entrez_id, gene_symbol, gene_alias or refseq_id.
#'   And parameter IDs2Add can be set to any combination of identifiers 
#'   such as "accnum", "ensembl", "ensemblprot", "ensembltrans", "entrez_id",
#'   "enzyme", "genename", "pfam", "pmid", "prosite", "refseq", "symbol", 
#'   "unigene" and "uniprot". Some IDs are unique to an organism, 
#'   such as "omim" for org.Hs.eg.db and "mgi" for org.Mm.eg.db.
#'   
#'   Here is the definition of different IDs : 
#'     \itemize{
#'       \item accnum: GenBank accession numbers
#'       \item ensembl: Ensembl gene accession numbers
#'       \item ensemblprot: Ensembl protein accession numbers
#'       \item ensembltrans: Ensembl transcript accession numbers
#'       \item entrez_id: entrez gene identifiers
#'       \item enzyme: EC numbers
#'       \item genename: gene name
#'       \item pfam: Pfam identifiers
#'       \item pmid: PubMed identifiers
#'       \item prosite: PROSITE identifiers
#'       \item refseq: RefSeq identifiers
#'       \item symbol: gene abbreviations
#'       \item unigene: UniGene cluster identifiers
#'       \item uniprot: Uniprot accession numbers
#'       \item omim: OMIM(Mendelian Inheritance in Man) identifiers
#'       \item mgi: Jackson Laboratory MGI gene accession numbers
#'     }
#'   
#'   \item If mart is used instead of orgAnn, for valid parameter 
#'   feature_id_type and IDs2Add parameters, please refer to 
#'   \link[biomaRt:getBM]{getBM} in bioMart package. 
#'   Parameter feature_id_type should be one valid filter name listed by 
#'   \link[biomaRt:listFilters]{listFilters(mart)} such as ensembl_gene_id.
#'   And parameter IDs2Add should be one or more valid attributes name listed 
#'   by \link[biomaRt:listAttributes]{listAttributes(mart)} such as 
#'   external_gene_id, entrezgene, wikigene_name, or mirbase_transcript_name.
#'   
#' }
#' @return GRanges if the input is a GRanges (with added ID columns in mcols), 
#' or a data.frame if input is a character vector.
#' @references http://www.bioconductor.org/packages/release/data/annotation/
#' @author Jianhong Ou, Lihua Julie Zhu
#' @seealso \link[biomaRt:getBM]{getBM}, AnnotationDb
#' @export
#' @importFrom AnnotationDbi mget
#' @importFrom biomaRt getBM
#' @keywords misc
#' @examples
#' data(annotatedPeak)
#' library(org.Hs.eg.db)
#' addGeneIDs(annotatedPeak[1:6,],orgAnn="org.Hs.eg.db",
#'            IDs2Add=c("symbol","omim"))
#' ##addGeneIDs(annotatedPeak$feature[1:6],orgAnn="org.Hs.eg.db",
#' ##           IDs2Add=c("symbol","genename"))
#' if(interactive()){
#'   mart <- useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",
#'                   dataset="hsapiens_gene_ensembl")
#'   ##mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#'   addGeneIDs(annotatedPeak[1:6,], mart=mart,
#'              IDs2Add=c("hgnc_symbol","entrezgene"))
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
    if (is(annotatedPeak, "GRanges")) {
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
        if (is(orgAnn, "OrgDb")) {
            orgAnn <- deparse(substitute(orgAnn))
        }
        if (!is(orgAnn, "character")) {
            stop("orgAnn must be a character.")
        }
        if (!grepl(".eg.db", orgAnn, ignore.case = TRUE)) {
            stop("Annotation database must be *.eg.db", call. = FALSE)
        }
        
        # Check if package is installed and available
        if (!requireNamespace(orgAnn, quietly = TRUE)) {
            stop("Package '", orgAnn, "' is not installed.",
                 call. = FALSE)
        }
        
        # Attach package namespace (required for accessing annotation objects)
        if (!requireNamespace(orgAnn, quiet = TRUE)) {
            stop("package '", orgAnn, "' is not installed.",
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
        
        if (length(entrezIDs) == 0) {
            stop("No entrez identifier can be mapped by input data based on ",
                 "the feature_id_type.\nPlease consider to use correct ",
                 "feature_id_type, orgAnn or annotatedPeak\n", call. = FALSE)
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
                
                if (!is(orgDB, "AnnDbBimap") && !is(orgDB, "IpiAnnDbMap")) {
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
        if (missing(mart) || !is(mart, "Mart")) {
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
    if (ncol(m_ent) == 1) {
        stop("None of IDs could be appended. Please double check IDs2Add.")
    }
    
    duplicated_ids <- m_ent[duplicated(m_ent[, feature_id_type]), 
                            feature_id_type]
    if (length(duplicated_ids) > 0) {
        m_ent.duplicated <- m_ent[m_ent[, feature_id_type] %in% duplicated_ids, ]
        m_ent.duplicated <- condenseMatrixByColnames(as.matrix(m_ent.duplicated),
                                                      feature_id_type)
        m_ent <- m_ent[!(m_ent[, feature_id_type] %in% duplicated_ids), ]
        m_ent <- rbind(m_ent, m_ent.duplicated)
    }
    
    # Merge back to original input
    if (is(annotatedPeak, "GRanges")) {
        # Rearrange m_ent by annotatedPeak$feature
        # data.frame is very important for order...
        orderlist <- data.frame(annotatedPeak$feature)
        orderlist <- cbind(seq_len(nrow(orderlist)), orderlist)
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
