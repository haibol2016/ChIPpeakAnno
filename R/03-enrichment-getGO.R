#' Obtain Gene Ontology (GO) terms for given genes
#' 
#' @description 
#' Retrieves Gene Ontology (GO) terms for a set of genes using organism-specific
#' annotation packages (e.g., \code{org.Hs.eg.db} for human). The function maps
#' gene identifiers (Ensembl, RefSeq, or gene symbols) to Entrez IDs and then
#' retrieves associated GO terms with their definitions, evidence codes, and
#' ontology types.
#' 
#' This function is useful for obtaining all GO annotations for a set of genes
#' without performing enrichment analysis. For enrichment testing, use
#' \code{\link{getEnrichedGO}} instead.
#' 
#' @param all.genes A character vector of gene identifiers. The type of identifier
#'        is specified by \code{ID_type}. Duplicate and NA values are automatically
#'        removed. At least 2 unique genes must be provided after conversion to
#'        Entrez IDs.
#' @param orgAnn Character string specifying the organism annotation package name.
#'        Common options include:
#'        \itemize{
#'          \item \code{"org.Hs.eg.db"}: Human (Homo sapiens)
#'          \item \code{"org.Mm.eg.db"}: Mouse (Mus musculus)
#'          \item \code{"org.Dm.eg.db"}: Fly (Drosophila melanogaster)
#'          \item \code{"org.Rn.eg.db"}: Rat (Rattus norvegicus)
#'          \item \code{"org.Sc.eg.db"}: Yeast (Saccharomyces cerevisiae)
#'          \item \code{"org.Dr.eg.db"}: Zebrafish (Danio rerio)
#'        }
#'        The package must be installed and loaded before calling this function.
#'        Default is \code{"org.Hs.eg.db"}.
#' @param ID_type Character string specifying the type of gene identifier in
#'        \code{all.genes}. Options:
#'        \itemize{
#'          \item \code{"ensembl_gene_id"} (default): Ensembl gene IDs (e.g., 
#'                "ENSG00000139618")
#'          \item \code{"refseq_id"}: RefSeq IDs (e.g., "NM_000492")
#'          \item \code{"gene_symbol"}: Gene symbols (e.g., "BRCA2")
#'        }
#'        Default is \code{"gene_symbol"}.
#' @param writeTo Optional character string specifying a file path to write 
#'        the results table. If provided, results are written as a tab-separated
#'        file. If missing, results are not written to file.
#' 
#' @return Returns an invisible data frame with one row per gene-GO term
#'        association. The data frame contains the following columns:
#'        \itemize{
#'          \item \code{Alternate.ID}: Original gene identifier from \code{all.genes}
#'          \item \code{go.id}: GO term identifier (e.g., "GO:0003674")
#'          \item \code{go.term}: GO term name (e.g., "molecular_function")
#'          \item \code{Definition}: Full definition of the GO term
#'          \item \code{Evidence}: Evidence code indicating how the annotation
#'                was determined (e.g., "IDA", "IEA", "TAS", "IMP", etc.)
#'          \item \code{Ontology}: Ontology type: \code{"BP"} (Biological Process),
#'                \code{"CC"} (Cellular Component), or \code{"MF"} (Molecular Function)
#'          \item \code{EntrezID}: Entrez gene identifier (numeric)
#'        }
#'        
#'        Note: The function returns the result invisibly (using \code{invisible()}),
#'        so it can be assigned to a variable or printed directly. If a gene has
#'        multiple GO terms, it will appear in multiple rows.
#' 
#' @details
#' 
#' \strong{Algorithm:}
#' \enumerate{
#'   \item Converts input gene identifiers to Entrez IDs using the appropriate
#'         mapping (ENSEMBL2EG, SYMBOL2EG, or REFSEQ2EG)
#'   \item Retrieves all GO annotations for the mapped Entrez IDs from the
#'         organism-specific GO annotation package
#'   \item Annotates GO terms with their names, definitions, and ontology types
#'   \item Maps results back to original input identifiers
#'   \item Optionally writes results to a file
#' }
#' 
#' \strong{Requirements:}
#' \itemize{
#'   \item The specified organism annotation package (e.g., \code{org.Hs.eg.db})
#'         must be installed and loaded
#'   \item At least 2 unique genes must be successfully mapped to Entrez IDs
#'   \item Input genes that cannot be mapped to Entrez IDs are excluded from
#'         the results
#' }
#' 
#' \strong{Note:} If multiple Entrez IDs map to the same input identifier, only
#' the first mapping is used. Genes with no GO annotations will not appear in
#' the results.
#' @author Lihua Julie Zhu
#' @seealso \code{\link{getEnrichedGO}} for GO enrichment analysis,
#'          \code{\link{addGeneIDs}} for adding gene identifiers to annotated peaks
#' @keywords misc
#' @export
#' @importFrom AnnotationDbi mappedkeys Definition Ontology Term
#' @importFrom utils write.table
#' @importFrom S4Vectors elementNROWS
#' @examples
#' \dontrun{
#'   ## Example 1: Get GO terms for Ensembl gene IDs
#'   library(org.Hs.eg.db)
#'   data(annotatedPeak)
#'   go_terms <- getGO(annotatedPeak[1:6]$feature, 
#'                     orgAnn = "org.Hs.eg.db", 
#'                     ID_type = "ensembl_gene_id")
#'   head(go_terms)
#'   
#'   ## Example 2: Get GO terms for gene symbols
#'   gene_symbols <- c("BRCA1", "BRCA2", "TP53", "ATM", "CHEK2")
#'   go_terms <- getGO(gene_symbols, 
#'                     orgAnn = "org.Hs.eg.db", 
#'                     ID_type = "gene_symbol")
#'   
#'   ## Example 3: Write results to file
#'   getGO(gene_symbols, 
#'         orgAnn = "org.Hs.eg.db", 
#'         ID_type = "gene_symbol",
#'         writeTo = "go_terms.txt")
#' }
#' 
getGO <- function(all.genes, orgAnn = "org.Hs.eg.db", 
                  writeTo, ID_type = "gene_symbol") {
    
    GOgenome <- sub(".db", "", orgAnn)
    
    cov2EntrezID <- function(IDs, orgAnn, ID_type = "ensembl_gene_id") {
        GOgenome <- sub(".db", "", orgAnn)
        
        orgAnn <- switch(ID_type, 
                         ensembl_gene_id = "ENSEMBL2EG", 
                         gene_symbol = "SYMBOL2EG", 
                         refseq_id = "REFSEQ2EG", 
                         "BAD_NAME")
        
        if (orgAnn == "BAD_NAME") {
            stop("Currently only the following type of IDs are supported: ",
                 "ensembl_gene_id, refseq_id and gene_symbol!", call. = FALSE)
        }
        
        orgAnn <- get(paste(GOgenome, orgAnn, sep = ""))
        
        if (!inherits(orgAnn, "AnnDbBimap")) {
            stop("orgAnn is not a valid annotation dataset! ",
                 "For example, org.Hs.eg.db package for human and ",
                 "the org.Mm.eg.db package for mouse.", call. = FALSE)
        }
        
        IDs <- unique(IDs[!is.na(IDs)])
        
        ids <- mget(IDs, orgAnn, ifnotfound = NA)
        
        ids <- lapply(ids, `[`, 1)
        
        ids <- unlist(ids)
        
        ids <- unique(ids[!is.na(ids)])
        
        ids
    }
    
    entrezIDs <- cov2EntrezID(as.character(all.genes), orgAnn, ID_type)
    
    if (length(entrezIDs) < 2L) {
        stop("The number of genes is less than 2. ",
             "Please double check your feature_id_type.", call. = FALSE)
    }
    
    goAnn <- get(paste(GOgenome, "GO", sep = ""))
    
    mapped_genes <- mappedkeys(goAnn)
    
    totalN.genes <- length(unique(mapped_genes))
    
    thisN.genes <- length(unique(entrezIDs))
    
    xx <- mget(mapped_genes, goAnn, ifnotfound = NA)
    
    all.GO <- cbind(matrix(unlist(unlist(xx)), ncol = 3, byrow = TRUE), 
                    rep(names(xx), elementNROWS(xx)))
    
    all.GO <- unique(all.GO)
    
    this.GO <- all.GO[all.GO[, 4] %in% entrezIDs, , drop = FALSE]
    
    annoTerms <- function(goids) {
        if (length(goids) < 1) {
            goterm <- matrix(ncol = 4)
        } else {
            goids <- as.character(goids)
            
            terms <- Term(goids)
            
            definition <- Definition(goids)
            
            ontology <- Ontology(goids)
            
            goterm <- cbind(goids, terms[match(goids, names(terms))], 
                           definition[match(goids, names(definition))]) 
            
            ## ontology[match(goids, names(ontology))] - not currently used
        }
        
        colnames(goterm) <- c("go.id", "go.term", "Definition")
        
        rownames(goterm) <- NULL
        
        goterm
    }
    
    this.GO.withTerms <- unique(annoTerms(this.GO[, 1]))
    
    colnames(this.GO) <- c("go.id", "Evidence", "Ontology", "EntrezID")
    
    IDs <- as.character(all.genes)
    
    GOgenome <- sub(".db", "", orgAnn)
    
    orgAnn2 <- switch(ID_type, 
                      ensembl_gene_id = "ENSEMBL2EG", 
                      gene_symbol = "SYMBOL2EG", 
                      refseq_id = "REFSEQ2EG", 
                      "BAD_NAME")
    
    if (orgAnn2 == "BAD_NAME") {
        stop("Currently only the following type of IDs are supported: ",
             "ensembl_gene_id, refseq_id and gene_symbol!", call. = FALSE)
    }
    
    orgAnn2 <- get(paste(GOgenome, orgAnn2, sep = ""))
    
    if (!inherits(orgAnn2, "AnnDbBimap")) {
        stop("orgAnn is not a valid annotation dataset! ",
             "For example, org.Hs.eg.db package for human and ",
             "the org.Mm.eg.db package for mouse.", call. = FALSE)
    }
    
    ids <- mget(IDs, orgAnn2, ifnotfound = NA)
    
    ids.withInputIDs <- cbind(names(unlist(ids)), unlist(ids))
    
    ids.withInputIDs <- subset(ids.withInputIDs, !is.na(ids.withInputIDs[, 2]))
    
    this.GO.withTermsEntrezIDs <- merge(this.GO, this.GO.withTerms)
    
    this.GO.withTermsEntrezIDs <- 
        cbind(ids.withInputIDs[match(this.GO.withTermsEntrezIDs[, 4], 
                                     ids.withInputIDs[, 2]), 1],
              this.GO.withTermsEntrezIDs)
    
    colnames(this.GO.withTermsEntrezIDs)[1] <- "Alternate.ID"
    
    if (!missing(writeTo)) {
        write.table(this.GO.withTermsEntrezIDs,
                   file = writeTo, sep = "\t",
                   row.names = FALSE)
    }
    
    return(invisible(this.GO.withTermsEntrezIDs))
}
