#' Obtain enriched Gene Ontology (GO) terms for peaks
#' 
#' @description 
#' Identifies enriched Gene Ontology (GO) terms based on features associated
#' with peaks. The function uses GO.db and organism-specific annotation packages
#' (e.g., \code{org.Hs.eg.db}) to retrieve GO annotations and performs
#' hypergeometric tests to identify significantly enriched terms. Multiple
#' testing correction is available via the \code{multtest} package.
#' 
#' The function automatically adds ancestor GO terms to the analysis, meaning
#' that if a gene is annotated to a specific GO term, it is also considered
#' annotated to all ancestor terms in the GO hierarchy. This provides a more
#' comprehensive enrichment analysis.
#' 
#' @param annotatedPeak A \code{GRanges} object with a \code{feature} column
#'        containing feature IDs (e.g., gene IDs), a character vector of feature
#'        IDs, or a \code{GRangesList} of peak sets. For \code{GRangesList}, the
#'        function processes each element separately and returns a list of results.
#'        Note: \code{subGroupComparison} is not supported for \code{GRangesList}.
#' @param orgAnn A character string specifying the organism annotation package
#'        name (e.g., \code{"org.Hs.eg.db"} for human, \code{"org.Mm.eg.db"} for
#'        mouse, \code{"org.Dm.eg.db"} for fly, \code{"org.Rn.eg.db"} for rat,
#'        \code{"org.Sc.eg.db"} for yeast, \code{"org.Dr.eg.db"} for zebrafish).
#'        The package must be installed and loaded. The function extracts the
#'        organism name by removing ".db" from the package name.
#' @param feature_id_type A character string specifying the type of feature IDs
#'        in \code{annotatedPeak}. Supported types:
#'        \itemize{
#'          \item \code{"ensembl_gene_id"} (default): Ensembl gene IDs
#'          \item \code{"refseq_id"}: RefSeq IDs
#'          \item \code{"gene_symbol"}: Gene symbols
#'          \item \code{"entrez_id"}: Entrez gene IDs (no conversion needed)
#'        }
#'        The function converts these IDs to Entrez gene IDs internally for GO
#'        annotation lookup.
#' @param maxP A numeric value specifying the maximum p-value (or adjusted
#'        p-value if \code{multiAdjMethod} is specified) to be considered
#'        significant. Default is \code{0.01}. GO terms with p-values greater
#'        than this threshold are filtered out.
#' @param minGOterm An integer specifying the minimum count in the genome for a
#'        GO term to be included in the results. This filters out very specific
#'        GO terms with few annotated genes. Default is \code{10}.
#' @param multiAdjMethod A character string specifying the multiple testing
#'        correction method. If \code{NULL} (default), no correction is applied.
#'        Available methods (from \code{multtest} package): \code{"Bonferroni"},
#'        \code{"Holm"}, \code{"Hochberg"}, \code{"SidakSS"}, \code{"SidakSD"},
#'        \code{"BH"} (Benjamini-Hochberg), \code{"BY"} (Benjamini-Yekutieli),
#'        \code{"ABH"}, \code{"TSBH"}. When specified, an additional column
#'        \code{<method>.adjusted.p.value} is added to the results, and filtering
#'        by \code{maxP} uses the adjusted p-value instead of the raw p-value.
#' @param condense A logical value. If \code{TRUE}, condenses the GO term to
#'        Entrez ID mapping matrix using \code{condenseMatrixByColnames}. This
#'        can reduce memory usage for large datasets. Default is \code{FALSE}.
#' @param removeAncestorByPval A numeric value or \code{NULL} (default). If
#'        provided, removes parent GO terms when all their children terms have
#'        significantly more genes (based on Fisher's exact test p-value <
#'        \code{removeAncestorByPval}). This helps reduce redundancy in the GO
#'        hierarchy. Uses the \code{removeAncestor} function internally.
#' @param keepByLevel An integer or \code{NULL} (default). If provided, filters
#'        GO terms based on their level in the GO hierarchy. Terms with a
#'        shortest path to 'all' greater than the specified level are removed.
#'        Uses the \code{filterByLevel} function internally. This helps focus on
#'        terms at a specific level of specificity.
#' @param subGroupComparison A logical vector of the same length as
#'        \code{annotatedPeak} (when it's a \code{GRanges} or character vector),
#'        or \code{NULL} (default). If provided, splits the peaks into two groups
#'        based on \code{TRUE}/\code{FALSE} values and performs comparative
#'        enrichment analysis:
#'        \itemize{
#'          \item \code{TRUE} group: Enrichment analysis using hypergeometric test
#'          \item \code{FALSE} group: Used as background for comparison
#'        }
#'        The analysis adds two additional columns to the results:
#'        \itemize{
#'          \item \code{count.InBackgroundDataset}: Count of genes annotated to
#'                the GO term in the FALSE group
#'          \item \code{dataset.vs.background.pval}: P-value from Fisher's exact
#'                test comparing TRUE group vs FALSE group. Only GO terms with
#'                this p-value < \code{maxP} are retained.
#'        }
#'        Note: To compare FALSE vs TRUE, repeat the analysis with inverted
#'        \code{subGroupComparison} values.
#' @return Returns a list with 3 elements, each containing a data frame of
#'        enriched GO terms:
#'        \itemize{
#'          \item \code{bp}: Enriched biological process (BP) terms
#'          \item \code{mf}: Enriched molecular function (MF) terms
#'          \item \code{cc}: Enriched cellular component (CC) terms
#'        }
#'        
#'        Each data frame contains the following columns (minimum 9 columns):
#'        \itemize{
#'          \item \code{go.id}: GO term identifier (e.g., "GO:0008150")
#'          \item \code{go.term}: GO term name (e.g., "biological_process")
#'          \item \code{Definition}: Detailed description of the GO term
#'          \item \code{Ontology}: Ontology branch ("BP", "MF", or "CC")
#'          \item \code{count.InDataset}: Number of genes in the input dataset
#'                annotated to this GO term (including ancestor terms)
#'          \item \code{count.InGenome}: Number of genes in the genome annotated
#'                to this GO term (including ancestor terms)
#'          \item \code{pvalue}: P-value from the hypergeometric test
#'          \item \code{totaltermInDataset}: Total number of GO term annotations
#'                in the input dataset (sum across all terms)
#'          \item \code{totaltermInGenome}: Total number of GO term annotations
#'                in the genome (sum across all terms)
#'        }
#'        
#'        \strong{Additional columns when \code{multiAdjMethod} is specified:}
#'        \itemize{
#'          \item \code{<method>.adjusted.p.value}: Adjusted p-value using the
#'                specified multiple testing correction method
#'        }
#'        
#'        \strong{Additional columns when \code{subGroupComparison} is provided:}
#'        \itemize{
#'          \item \code{count.InBackgroundDataset}: Number of genes in the FALSE
#'                group annotated to this GO term
#'          \item \code{dataset.vs.background.pval}: P-value from Fisher's exact
#'                test comparing TRUE group vs FALSE group
#'        }
#'        
#'        \strong{Additional columns when \code{condense = TRUE}:}
#'        \itemize{
#'          \item \code{EntrezID}: CharacterList of Entrez gene IDs annotated to
#'                this GO term
#'        }
#'        
#'        \strong{Note:} If \code{annotatedPeak} is a \code{GRangesList}, the
#'        function returns a list of results, where each element corresponds to
#'        one element of the input \code{GRangesList}.
#' @details
#' 
#' \strong{Input Processing:}
#' \itemize{
#'   \item If \code{annotatedPeak} is a \code{GRanges} object, the function
#'         extracts unique feature IDs from the \code{feature} metadata column
#'   \item If \code{annotatedPeak} is a character vector, it uses the values
#'         directly as feature IDs
#'   \item If \code{annotatedPeak} is a \code{GRangesList}, each element is
#'         processed separately and results are returned as a list
#' }
#' 
#' \strong{ID Conversion:}
#' The function converts feature IDs to Entrez gene IDs using organism-specific
#' annotation packages. Supported conversions:
#' \itemize{
#'   \item Ensembl gene ID → Entrez ID (via \code{ENSEMBL2EG})
#'   \item Gene symbol → Entrez ID (via \code{SYMBOL2EG})
#'   \item RefSeq ID → Entrez ID (via \code{REFSEQ2EG})
#'   \item Entrez ID → Entrez ID (no conversion needed)
#' }
#' 
#' \strong{GO Annotation and Ancestor Addition:}
#' The function retrieves GO annotations for all genes in the genome and the
#' input dataset. It then automatically adds ancestor GO terms, meaning that if
#' a gene is annotated to a specific GO term, it is also considered annotated to
#' all parent terms in the GO hierarchy. This provides more comprehensive
#' enrichment results.
#' 
#' \strong{Statistical Testing:}
#' \itemize{
#'   \item \strong{Hypergeometric test}: Used to test for enrichment of GO terms
#'         in the input dataset compared to the genome background. The test
#'         parameters are:
#'         \itemize{
#'           \item \code{q}: Number of genes in dataset annotated to the GO term
#'           \item \code{m}: Number of genes in genome annotated to the GO term
#'           \item \code{n}: Total number of GO annotations in genome
#'           \item \code{k}: Total number of GO annotations in dataset
#'         }
#'   \item \strong{Fisher's exact test}: When \code{subGroupComparison} is used,
#'         this test compares the TRUE group vs FALSE group for each enriched
#'         GO term identified by the hypergeometric test
#' }
#' 
#' \strong{Filtering and Post-processing:}
#' \itemize{
#'   \item GO terms are filtered by \code{maxP} (raw or adjusted p-value)
#'   \item GO terms are filtered by \code{minGOterm} (minimum genome count)
#'   \item If \code{removeAncestorByPval} is specified, parent terms are removed
#'         when all children have significantly more genes
#'   \item If \code{keepByLevel} is specified, terms beyond the specified level
#'         are removed
#' }
#' 
#' @author Lihua Julie Zhu. Jianhong Ou for subGroupComparison
#' @seealso \code{\link[stats]{phyper}}, \code{\link{hyperGtest}},
#'          \code{\link{removeAncestor}}, \code{\link{filterByLevel}},
#'          \code{\link{condenseMatrixByColnames}}
#' @references Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) Univariate
#' Discrete Distributions, Second Edition. New York: Wiley
#' @keywords misc
#' @export
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom multtest mt.rawp2adjp
#' @importFrom AnnotationDbi Definition Ontology Term
#' @importFrom stats phyper fisher.test
#' @examples
#' 
#'   data(enrichedGO)
#'   enrichedGO$mf[1:10,]
#'   enrichedGO$bp[1:10,]
#'   enrichedGO$cc
#'   if (interactive()) {
#'      data(annotatedPeak)
#'      library(org.Hs.eg.db)
#'      library(GO.db)
#'      enriched.GO = getEnrichedGO(annotatedPeak[1:6,], 
#'                                  orgAnn="org.Hs.eg.db", 
#'                                  maxP=0.01,
#'                                  minGOterm=10,
#'                                  multiAdjMethod= NULL)
#'      dim(enriched.GO$mf)
#'      colnames(enriched.GO$mf)
#'      dim(enriched.GO$bp)
#'      enriched.GO$cc
#' }
#' 
getEnrichedGO <- function(annotatedPeak, orgAnn, 
                          feature_id_type = "ensembl_gene_id", 
                          maxP = 0.01,
                          minGOterm = 10, multiAdjMethod = NULL,
                          condense = FALSE,
                          removeAncestorByPval = NULL,
                          keepByLevel = NULL,
                          subGroupComparison = NULL) {
    if (inherits(annotatedPeak, "GRangesList")) {
        if (length(subGroupComparison) > 0L) {
            stop("subGroupComparison parameter is not supported for list of peaks.", 
                 call. = FALSE)
        }
        args <- as.list(match.call())
        res <- lapply(annotatedPeak, function(.ele) {
            args$annotatedPeak <- .ele
            do.call(getEnrichedGO, args = args)
        })
    }
    stopifnot("The 'GO.db' package is required" =
                requireNamespace("GO.db", quietly = TRUE)) 
    if (missing(annotatedPeak)) {
        stop("Missing required argument 'annotatedPeak'!", call. = FALSE)    
    }
    if (length(multiAdjMethod) > 0) {
        multiAdjMethod <- match.arg(multiAdjMethod, 
                                    c("Bonferroni", "Holm", "Hochberg", 
                                      "SidakSS", "SidakSD", "BH", "BY",
                                      "ABH", "TSBH"))
    }
    if (missing(orgAnn)) {
        stop("Missing required argument 'orgAnn'. ",
             "Please refer to ",
             "http://www.bioconductor.org/packages/release/data/annotation/ ",
             "for available org.xx.eg.db packages", call. = FALSE)
    }
    GOgenome <- sub(".db", "", orgAnn)
    if (nchar(GOgenome) < 1L) {
        stop("Invalid 'orgAnn' parameter. ",
             "Please refer to ",
             "http://www.bioconductor.org/packages/release/data/annotation/ ",
             "for available org.xx.eg.db packages", call. = FALSE)
    }
    groupFALSE <- NULL
    feature_ids_FALSE <- NULL
    entrezIDs_FALSE <- NULL
    if (length(subGroupComparison) > 0L) {
        if (length(subGroupComparison) != length(annotatedPeak)) {
            stop("Length of 'subGroupComparison' must match ",
                 "length of 'annotatedPeak'", call. = FALSE)
        }
        stopifnot(is.logical(subGroupComparison))
        groupFALSE <- annotatedPeak[!subGroupComparison]
        annotatedPeak <- annotatedPeak[subGroupComparison]
    }
    if (inherits(annotatedPeak, "GRanges")) {
        feature_ids <- unique(as.character(annotatedPeak$feature))
        if (length(groupFALSE)) {
            feature_ids_FALSE <- unique(as.character(groupFALSE$feature))
        }
    } else if (is.character(annotatedPeak)) {
        feature_ids <- unique(annotatedPeak)
        if (length(groupFALSE)) {
            feature_ids_FALSE <- unique(groupFALSE)
        }
    } else {
        stop("'annotatedPeak' must be a GRanges object with a 'feature' ",
             "column or a character vector of feature IDs", call. = FALSE)
    }
    if (feature_id_type == "entrez_id") {
        entrezIDs <- feature_ids
        if (length(groupFALSE)) {
            entrezIDs_FALSE <- feature_ids_FALSE
        }
    } else {
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
        entrezIDs <- cov2EntrezID(feature_ids, orgAnn, feature_id_type)
        if (length(groupFALSE)) {
            entrezIDs_FALSE <- cov2EntrezID(feature_ids_FALSE,
                                             orgAnn, feature_id_type)
        }
    }
    if (length(entrezIDs) < 2L) {
        stop("The number of genes is less than 2. ",
             "Please double check your feature_id_type.", call. = FALSE)
    }
    goAnn <- get(paste(GOgenome, "GO", sep = ""))
    mapped_genes <- mappedkeys(goAnn)
    totalN.genes <- length(unique(mapped_genes))
    thisN.genes <- length(unique(entrezIDs))
    #xx <- as.list(goAnn[mapped_genes])
    
    xx <- mget(mapped_genes, goAnn, ifnotfound = NA)
    all.GO <- cbind(matrix(unlist(unlist(xx)), ncol = 3, byrow = TRUE),
                    rep(names(xx), elementNROWS(xx)))
    all.GO <- unique(all.GO)  ## incase the database is not unique
    this.GO <- all.GO[all.GO[, 4] %in% entrezIDs, , drop = FALSE]
    FALSE.GO <- all.GO[all.GO[, 4] %in% entrezIDs_FALSE, , drop = FALSE]
    
    addAnc <- function(go.ids,
                       onto = c("bp", "cc", "mf")) {
        ##replace the function addAncestors
        GOIDs <- unique(as.character(go.ids[, 1]))
        empty <- matrix(nrow = 0, ncol = 2)
        colnames(empty) <- c("go.id", "EntrezID")
        if (length(GOIDs) < 1) {
            return(empty)
        }
        onto <- match.arg(onto)
        Ancestors <- switch(onto, 
                            bp = mget(GOIDs, GO.db::GOBPANCESTOR, ifnotfound = NA),
                            cc = mget(GOIDs, GO.db::GOCCANCESTOR, ifnotfound = NA),
                            mf = mget(GOIDs, GO.db::GOMFANCESTOR, ifnotfound = NA))
        Ancestors <- unique(unlist(Ancestors))
        Ancestors <- Ancestors[Ancestors != "all" & !is.na(Ancestors)]
        if (length(Ancestors) > 0) {
            children <- go.ids[, c(1, 4), drop = FALSE]
            children <- unique(children)
            children.s <- split(children[, 2], children[, 1])
            Ancestors <- switch(onto,
                               bp = mget(Ancestors, 
                                       GO.db::GOBPOFFSPRING, 
                                       ifnotfound = NA),
                               cc = mget(Ancestors, 
                                       GO.db::GOCCOFFSPRING, 
                                       ifnotfound = NA),
                               mf = mget(Ancestors, 
                                       GO.db::GOMFOFFSPRING, 
                                       ifnotfound = NA))
            Ancestors <- cbind(Ancestor = rep(names(Ancestors),
                                           elementNROWS(Ancestors)),
                              child = unlist(Ancestors))
            Ancestors <- 
                Ancestors[Ancestors[, "child"] %in% names(children.s), 
                          , drop = FALSE]
            temp <- cbind(Ancestor = children[, 1], child = children[, 1])
            Ancestors <- rbind(Ancestors, temp)
            Ancestors <- unique(Ancestors)
            Ancestors.ezid <- children.s[Ancestors[, "child"]]
            temp <- cbind(rep(Ancestors[, "Ancestor"], 
                                    elementNROWS(Ancestors.ezid)),
                          unlist(Ancestors.ezid))
            temp <- unique(temp)  #time....
            if (length(temp) < 3) {
                re <- unique(children)
            } else {
                re <- temp[!is.na(temp[, 1]) & temp[, 1] != "", ]
            }
        } else {
            re <- unique(cbind(as.character(go.ids[, 1]), go.ids[, 4]))
        }
        colnames(re) <- c("go.id", "EntrezID")
        re
    }
    
    bp.go.this.withEntrez <- addAnc(this.GO[this.GO[, 3] == "BP", ], "bp")
    cc.go.this.withEntrez <- addAnc(this.GO[this.GO[, 3] == "CC", ], "cc")
    mf.go.this.withEntrez <- addAnc(this.GO[this.GO[, 3] == "MF", ], "mf")
    
    bp.go.all.withEntrez <- addAnc(all.GO[all.GO[, 3] == "BP", ], "bp")
    cc.go.all.withEntrez <- addAnc(all.GO[all.GO[, 3] == "CC", ], "cc")
    mf.go.all.withEntrez <- addAnc(all.GO[all.GO[, 3] == "MF", ], "mf")
    
    bp.go.this <- bp.go.this.withEntrez[, 1]
    cc.go.this <- cc.go.this.withEntrez[, 1]
    mf.go.this <- mf.go.this.withEntrez[, 1]
    
    bp.go.all <- bp.go.all.withEntrez[, 1]
    cc.go.all <- cc.go.all.withEntrez[, 1]
    mf.go.all <- mf.go.all.withEntrez[, 1]
    
    total.mf <- length(mf.go.all)
    total.cc <- length(cc.go.all)
    total.bp <- length(bp.go.all)
    this.mf <- length(mf.go.this)
    this.cc <- length(cc.go.this)
    this.bp <- length(bp.go.this)
    
    this.bp.count <- table(as.character(bp.go.this[bp.go.this!=""]))
    this.mf.count <- table(as.character(mf.go.this[mf.go.this!=""]))
    this.cc.count <- table(as.character(cc.go.this[cc.go.this!=""]))
    
    all.bp.count <- table(as.character(bp.go.all[bp.go.all!=""]))
    all.mf.count <- table(as.character(mf.go.all[mf.go.all!=""]))
    all.cc.count <- table(as.character(cc.go.all[cc.go.all!=""]))
    
    if (length(groupFALSE)) {
        bp.go.FALSE.withEntrez <- addAnc(FALSE.GO[FALSE.GO[, 3] == "BP", ], "bp")
        cc.go.FALSE.withEntrez <- addAnc(FALSE.GO[FALSE.GO[, 3] == "CC", ], "cc")
        mf.go.FALSE.withEntrez <- addAnc(FALSE.GO[FALSE.GO[, 3] == "MF", ], "mf")
        
        bp.go.FALSE <- bp.go.FALSE.withEntrez[, 1]
        cc.go.FALSE <- cc.go.FALSE.withEntrez[, 1]
        mf.go.FALSE <- mf.go.FALSE.withEntrez[, 1]
        
        FALSE.mf <- length(mf.go.FALSE)
        FALSE.cc <- length(cc.go.FALSE)
        FALSE.bp <- length(bp.go.FALSE)
        
        FALSE.bp.count <- table(as.character(bp.go.FALSE[bp.go.FALSE != ""]))
        FALSE.mf.count <- table(as.character(mf.go.FALSE[mf.go.FALSE != ""]))
        FALSE.cc.count <- table(as.character(cc.go.FALSE[cc.go.FALSE != ""]))
    }
    
    hyperGT <- function(alltermcount, thistermcount, 
                        totaltermInGenome, totaltermInPeakList) {
        m <- as.numeric(alltermcount[names(thistermcount)])
        q <- as.numeric(thistermcount)
        n <- as.numeric(totaltermInGenome)
        k <- as.numeric(totaltermInPeakList)
        pvalue <- phyper(q - 1, m, n - m, k, lower.tail = FALSE, log.p = FALSE)
        data.frame(go.id = names(thistermcount), 
                   count.InDataset = q,
                   count.InGenome = m, 
                   pvalue = pvalue,
                   totaltermInDataset = k, 
                   totaltermInGenome = n)
    }
    
    bp.selected <- hyperGT(all.bp.count,
                             this.bp.count, 
                             total.bp, 
                             this.bp)
    mf.selected <- hyperGT(all.mf.count,
                             this.mf.count, 
                             total.mf, 
                             this.mf)
    cc.selected <- hyperGT(all.cc.count,
                             this.cc.count, 
                             total.cc, 
                             this.cc)
    
    if (length(multiAdjMethod) < 1) {
        bp.s <- 
            bp.selected[
                as.numeric(as.character(bp.selected[, 4])) < maxP & 
                    as.numeric(as.character(bp.selected[, 3])) >= minGOterm, ]
        mf.s <- 
            mf.selected[
                as.numeric(as.character(mf.selected[, 4])) < maxP & 
                    as.numeric(as.character(mf.selected[, 3])) >= minGOterm, ]
        cc.s <- 
            cc.selected[
                as.numeric(as.character(cc.selected[, 4])) < maxP & 
                    as.numeric(as.character(cc.selected[, 3])) >= minGOterm, ]
    } else {
        procs <- c(multiAdjMethod)
        res <- mt.rawp2adjp(as.numeric(as.character(bp.selected[, 4])), procs)
        adjp <- unique(res$adjp)
        colnames(adjp)[1] <- colnames(bp.selected)[4]
        colnames(adjp)[2] <- paste(multiAdjMethod, "adjusted.p.value", sep = ".")
        bp.selected[, 4] <- as.numeric(as.character(bp.selected[, 4]))
        bp1 <- merge(bp.selected, adjp, all.x = TRUE)
        
        res <- mt.rawp2adjp(as.numeric(as.character(mf.selected[, 4])), procs)
        adjp <- unique(res$adjp)
        colnames(adjp)[1] <- colnames(mf.selected)[4]
        colnames(adjp)[2] <- paste(multiAdjMethod, "adjusted.p.value", sep = ".")
        mf.selected[, 4] <- as.numeric(as.character(mf.selected[, 4]))
        mf1 <- merge(mf.selected, adjp, all.x = TRUE)
        
        res <- mt.rawp2adjp(as.numeric(as.character(cc.selected[, 4])), procs)
        adjp <- unique(res$adjp)
        colnames(adjp)[1] <- colnames(cc.selected)[4]
        colnames(adjp)[2] <- paste(multiAdjMethod, "adjusted.p.value", sep = ".")
        cc.selected[, 4] <- as.numeric(as.character(cc.selected[, 4]))
        cc1 <- merge(cc.selected, adjp, all.x = TRUE)
        
        bp.s <- bp1[as.numeric(as.character(bp1[, dim(bp1)[2]])) < maxP &  
                       !is.na(bp1[, dim(bp1)[2]]) & 
                       as.numeric(as.character(bp1[, 4])) >= minGOterm, ]
        mf.s <- mf1[as.numeric(as.character(mf1[, dim(mf1)[2]])) < maxP & 
                       !is.na(mf1[, dim(mf1)[2]]) & 
                       as.numeric(as.character(mf1[, 4])) >= minGOterm, ]
        cc.s <- cc1[as.numeric(as.character(cc1[, dim(cc1)[2]])) < maxP & 
                       !is.na(cc1[, dim(cc1)[2]]) & 
                       as.numeric(as.character(cc1[, 4])) >= minGOterm, ]
    }
    
    annoTerms <- function(goids) {
        if (length(goids) < 1) {
            goterm <- matrix(ncol = 4)
        } else {
            goids <- as.character(goids)
            terms <- Term(goids)
            definition <- Definition(goids)
            ontology <- Ontology(goids)
            goterm <- cbind(goids, 
                            terms[match(goids, names(terms))],
                            definition[match(goids, names(definition))],
                            ontology[match(goids, names(ontology))])
        }
        colnames(goterm) <- c("go.id", "go.term", "Definition", "Ontology")
        rownames(goterm) <- NULL
        goterm
    }
    goterm.bp <- annoTerms(bp.s$go.id)
    goterm.mf <- annoTerms(mf.s$go.id)
    goterm.cc <- annoTerms(cc.s$go.id)
    
    bp.selected1 <- merge(goterm.bp, bp.s, by = "go.id")
    mf.selected1 <- merge(goterm.mf, mf.s, by = "go.id")
    cc.selected1 <- merge(goterm.cc, cc.s, by = "go.id")
    
    if (length(groupFALSE)) {  ## add counts for FALSE group and do Fisher's exact test
        bp.selected1$count.InBackgroundDataset <- FALSE.bp.count[bp.selected1$go.id]
        mf.selected1$count.InBackgroundDataset <- FALSE.mf.count[mf.selected1$go.id]
        cc.selected1$count.InBackgroundDataset <- FALSE.cc.count[cc.selected1$go.id]
        
        multiFisher <- function(.data) {
            rowFisher <- function(x, ...) {
                return(fisher.test(matrix(x, nrow = 2, byrow = TRUE), ...)$p.value)
            }
            .data0 <- .data[, c("count.InDataset", "count.InBackgroundDataset")]
            .data0 <- as.matrix(.data0)
            .data0[is.na(.data0)] <- 0
            .data0 <- cbind(.data0, .data[, "count.InGenome"] - .data0)
            apply(.data0, 1, rowFisher, alternative = "greater")
        }
        bp.selected1$dataset.vs.background.pval <- multiFisher(bp.selected1)
        mf.selected1$dataset.vs.background.pval <- multiFisher(mf.selected1)
        cc.selected1$dataset.vs.background.pval <- multiFisher(cc.selected1)
        
        bp.selected1 <-
            bp.selected1[bp.selected1$dataset.vs.background.pval < maxP, , drop = FALSE]
        mf.selected1 <-
            mf.selected1[mf.selected1$dataset.vs.background.pval < maxP, , drop = FALSE]
        cc.selected1 <-
            cc.selected1[cc.selected1$dataset.vs.background.pval < maxP, , drop = FALSE]
    }
    
    if (condense) {
        bp.go.this.withEntrez <- 
            condenseMatrixByColnames(bp.go.this.withEntrez, iname = "go.id")
        mf.go.this.withEntrez <- 
            condenseMatrixByColnames(mf.go.this.withEntrez, iname = "go.id")
        cc.go.this.withEntrez <- 
            condenseMatrixByColnames(cc.go.this.withEntrez, iname = "go.id")
    }
    
    bp.selected <- merge(bp.selected1, bp.go.this.withEntrez)
    mf.selected <- merge(mf.selected1, mf.go.this.withEntrez)
    cc.selected <- merge(cc.selected1, cc.go.this.withEntrez)
    
    if (length(removeAncestorByPval) > 0) {
        if (is.numeric(removeAncestorByPval[1]) || 
            is.integer(removeAncestorByPval[1])) {
            bp.selected <- removeAncestor(bp.selected, onto = "BP", 
                                          cutoffPvalue = removeAncestorByPval)
            mf.selected <- removeAncestor(mf.selected, onto = "MF", 
                                          cutoffPvalue = removeAncestorByPval)
            cc.selected <- removeAncestor(cc.selected, onto = "CC", 
                                          cutoffPvalue = removeAncestorByPval)
        }
    }
    
    if (length(keepByLevel) > 0) {
        if (is.numeric(keepByLevel[1]) || is.integer(keepByLevel[1])) {
            bp.selected <- filterByLevel(bp.selected, onto = "BP", level = keepByLevel)
            mf.selected <- filterByLevel(mf.selected, onto = "MF", level = keepByLevel)
            cc.selected <- filterByLevel(cc.selected, onto = "CC", level = keepByLevel)
        }
    }
    list(bp = bp.selected, mf = mf.selected, cc = cc.selected)
}
