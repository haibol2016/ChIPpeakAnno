#' Obtain enriched pathway terms for peaks
#'
#' @description 
#' Identifies enriched pathway terms (KEGG, Reactome, etc.) based on features
#' associated with peaks. The function uses pathway annotation packages
#' (e.g., \code{reactome.db}, \code{KEGGREST}) and organism-specific annotation
#' packages (e.g., \code{org.Hs.eg.db}) to retrieve pathway annotations and
#' performs hypergeometric tests to identify significantly enriched pathways.
#' Multiple testing correction is available via the \code{multtest} package.
#'
#' The function extracts feature IDs from the input (either from a \code{GRanges}
#' object's \code{feature} column or directly from a character vector), converts
#' them to Entrez Gene IDs if necessary, maps them to pathways, and performs
#' statistical enrichment analysis. Results are filtered by p-value threshold
#' and minimum pathway size.
#'
#' @param annotatedPeak A \code{GRanges} object with a \code{feature} column
#'   containing feature IDs (e.g., from \code{annotatePeakInBatch()}), or a
#'   character vector of feature IDs. If a \code{GRanges} object, the function
#'   extracts unique values from the \code{feature} column.
#' @param orgAnn Character string specifying the organism annotation package
#'   (e.g., \code{"org.Hs.eg.db"} for human, \code{"org.Mm.eg.db"} for mouse,
#'   \code{"org.Dm.eg.db"} for fly, \code{"org.Rn.eg.db"} for rat,
#'   \code{"org.Sc.sgd.db"} for yeast, \code{"org.Dr.eg.db"} for zebrafish).
#'   The package must be loaded before calling this function. See
#'   \url{http://www.bioconductor.org/packages/release/data/annotation/} for
#'   available packages.
#' @param pathAnn Character string specifying the pathway annotation source:
#'   \itemize{
#'     \item \code{"reactome.db"}: Reactome pathways from the reactome.db package
#'     \item \code{"KEGGREST"}: KEGG pathways via KEGGREST API (recommended)
#'   }
#'   Note: \code{"KEGG.db"} is deprecated and will produce an error. The
#'   specified package must be loaded before calling this function.
#' @param feature_id_type Character string specifying the type of gene identifier
#'   in \code{annotatedPeak}. Options:
#'   \itemize{
#'     \item \code{"ensembl_gene_id"} (default): Ensembl gene IDs without version
#'           numbers (e.g., "ENSG00000139618", not "ENSG00000139618.2")
#'     \item \code{"refseq_id"}: RefSeq IDs
#'     \item \code{"gene_symbol"}: Gene symbols (e.g., "TP53", "BRCA1")
#'     \item \code{"entrez_id"}: Entrez Gene IDs (no conversion needed)
#'   }
#'   All identifiers are converted to Entrez Gene IDs internally for pathway
#'   mapping. If using Ensembl IDs, ensure version numbers are removed (use
#'   \code{sub('\\.[^.]+$', '', ensembl_gene_id)} if needed).
#' @param maxP Numeric value specifying the maximum p-value (or adjusted p-value
#'   if \code{multiAdjMethod} is specified) to be considered significant. Pathways
#'   with p-values >= \code{maxP} are filtered out. Default is \code{0.01}.
#' @param minPATHterm Integer specifying the minimum number of genes in the genome
#'   that must be annotated to a pathway for it to be included in the results.
#'   Pathways with fewer genes are filtered out. Default is \code{10}.
#' @param multiAdjMethod Character string specifying the multiple testing correction
#'   method. Options: \code{"Bonferroni"}, \code{"Holm"}, \code{"Hochberg"},
#'   \code{"SidakSS"}, \code{"SidakSD"}, \code{"BH"} (Benjamini-Hochberg),
#'   \code{"BY"} (Benjamini-Yekutieli), \code{"ABH"}, \code{"TSBH"}. See
#'   \code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} in the \code{multtest}
#'   package for details. Default is \code{NULL} (no correction, uses raw p-values).
#'   When specified, the adjusted p-value is used for filtering with \code{maxP}.
#' @param subGroupComparison Optional logical vector of the same length as
#'   \code{annotatedPeak} to split peaks into two groups for comparative enrichment
#'   analysis. When provided:
#'   \itemize{
#'     \item Peaks where \code{subGroupComparison == TRUE} are used as the test
#'           group (dataset)
#'     \item Peaks where \code{subGroupComparison == FALSE} are used as the
#'           background group
#'     \item Enrichment analysis is performed for the TRUE group using hypergeometric
#'           test (comparing against all genes in the genome)
#'     \item Comparative analysis is performed using Fisher's Exact test comparing
#'           TRUE vs FALSE groups
#'     \item Only pathways with \code{dataset.vs.background.pval < maxP} are
#'           returned
#'   }
#'   To compare FALSE vs TRUE, repeat the analysis with inverted logical values.
#'   Default is \code{NULL} (no subgroup comparison).
#' @return Returns a data frame of enriched pathways with the following columns:
#'   \itemize{
#'     \item \code{path.id}: Pathway identifier (KEGG pathway ID like "hsa00010"
#'           or Reactome pathway ID like "R-HSA-109581")
#'     \item \code{path.term}: Pathway name/description
#'     \item \code{entrez_id}: Entrez Gene ID(s) associated with this pathway in
#'           the dataset (may be a list/vector if multiple genes)
#'     \item \code{count.InDataset}: Number of genes in the dataset that are
#'           annotated to this pathway
#'     \item \code{count.InGenome}: Number of genes in the genome that are
#'           annotated to this pathway
#'     \item \code{pvalue}: P-value from the hypergeometric test (raw p-value,
#'           or adjusted if \code{multiAdjMethod} is specified)
#'     \item \code{totaltermInDataset}: Total number of pathway annotations in
#'           the dataset (sum of all pathway-gene associations)
#'     \item \code{totaltermInGenome}: Total number of pathway annotations in
#'           the genome (sum of all pathway-gene associations)
#'     \item \code{<multiAdjMethod>.adjusted.p.value}: (if \code{multiAdjMethod}
#'           is specified) Adjusted p-value from the specified multiple testing
#'           correction method
#'     \item \code{count.InBackgroundDataset}: (if \code{subGroupComparison} is
#'           used) Number of genes in the background group (FALSE) that are
#'           annotated to this pathway
#'     \item \code{dataset.vs.background.pval}: (if \code{subGroupComparison} is
#'           used) P-value from Fisher's Exact test comparing the test group
#'           (TRUE) vs background group (FALSE)
#'   }
#'   
#'   If no enriched pathways are found (or all are filtered out), returns a
#'   data frame with a single row containing \code{log = "No enriched pathway found!"}.
#' @details
#' 
#' \strong{Input Requirements:}
#' \itemize{
#'   \item The organism annotation package (e.g., \code{org.Hs.eg.db}) must be
#'         loaded before calling this function
#'   \item The pathway annotation package (e.g., \code{reactome.db}) or
#'         \code{KEGGREST} must be loaded before calling this function
#'   \item For \code{GRanges} input, the object must have a \code{feature} column
#'         containing feature IDs
#'   \item At least 2 unique genes must be present after conversion to Entrez IDs
#' }
#' 
#' \strong{Pathway Annotation Sources:}
#' \itemize{
#'   \item \code{"reactome.db"}: Uses local Reactome pathway annotations. Requires
#'         the \code{reactome.db} package to be installed and loaded. Provides
#'         Reactome pathway IDs (e.g., "R-HSA-109581") and names.
#'   \item \code{"KEGGREST"}: Uses KEGG REST API via the \code{KEGGREST} package.
#'         Requires internet connection. Provides KEGG pathway IDs (e.g., "hsa00010")
#'         and names. Automatically maps organism annotation packages to KEGG
#'         organism codes (e.g., "org.Hs.eg.db" -> "hsa").
#' }
#' 
#' \strong{Statistical Analysis:}
#' \itemize{
#'   \item Hypergeometric test: Compares the number of genes in a pathway in the
#'         dataset vs the expected number based on the genome-wide distribution
#'   \item Multiple testing correction: When \code{multiAdjMethod} is specified,
#'         p-values are adjusted using the selected method, and the adjusted
#'         p-value is used for filtering
#'   \item Fisher's Exact test: When \code{subGroupComparison} is used, compares
#'         pathway enrichment between the test group (TRUE) and background group
#'         (FALSE) using a 2x2 contingency table
#' }
#' 
#' \strong{Filtering:}
#' \itemize{
#'   \item Pathways are filtered by: \code{pvalue < maxP} (or adjusted p-value
#'         if \code{multiAdjMethod} is specified)
#'   \item Pathways are filtered by: \code{count.InGenome >= minPATHterm}
#'   \item If \code{subGroupComparison} is used, additional filtering by:
#'         \code{dataset.vs.background.pval < maxP}
#' }
#' 
#' @author Jianhong Ou, Kai Hu
#' @seealso \code{\link{hyperGtest}}, \code{\link{convert2EntrezID}},
#'          \code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}},
#'          \code{\link[KEGGREST:keggGet]{keggGet}},
#'          \code{\link[KEGGREST:keggLink]{keggLink}}
#' @references Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) Univariate
#' Discrete Distributions, Second Edition. New York: Wiley
#' @keywords misc
#' @export
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom multtest mt.rawp2adjp
#' @importFrom KEGGREST keggGet keggLink
#' @examples
#' \dontrun{
#' ## Example 1: Reactome pathway enrichment
#' data(annotatedPeak)
#' library(org.Hs.eg.db)
#' library(reactome.db)
#' enriched.PATH <- getEnrichedPATH(annotatedPeak, 
#'                                  orgAnn = "org.Hs.eg.db",
#'                                  feature_id_type = "ensembl_gene_id",
#'                                  pathAnn = "reactome.db", 
#'                                  maxP = 0.01,
#'                                  minPATHterm = 10, 
#'                                  multiAdjMethod = NULL)
#' head(enriched.PATH)
#' 
#' ## Example 2: KEGG pathway enrichment
#' library(KEGGREST)
#' enrichedKEGG <- getEnrichedPATH(annotatedPeak, 
#'                                 orgAnn = "org.Hs.eg.db",
#'                                 feature_id_type = "ensembl_gene_id",
#'                                 pathAnn = "KEGGREST", 
#'                                 maxP = 0.01,
#'                                 minPATHterm = 10, 
#'                                 multiAdjMethod = "BH")
#' enrichmentPlot(enrichedKEGG)
#' 
#' ## Example 3: With multiple testing correction
#' enriched.PATH.adj <- getEnrichedPATH(annotatedPeak, 
#'                                      orgAnn = "org.Hs.eg.db",
#'                                      feature_id_type = "ensembl_gene_id",
#'                                      pathAnn = "reactome.db", 
#'                                      maxP = 0.05,
#'                                      minPATHterm = 10, 
#'                                      multiAdjMethod = "BH")
#' 
#' ## Example 4: Subgroup comparison
#' ## Split peaks into two groups based on some criteria
#' group1 <- c(rep(TRUE, 100), rep(FALSE, 50))
#' enriched.PATH.comp <- getEnrichedPATH(annotatedPeak, 
#'                                       orgAnn = "org.Hs.eg.db",
#'                                       feature_id_type = "ensembl_gene_id",
#'                                       pathAnn = "reactome.db", 
#'                                       maxP = 0.01,
#'                                       minPATHterm = 10, 
#'                                       subGroupComparison = group1)
#' }
#'
getEnrichedPATH <- function(annotatedPeak, orgAnn, pathAnn,
                            feature_id_type = "ensembl_gene_id",
                            maxP = 0.01, minPATHterm = 10, multiAdjMethod = NULL,
                            subGroupComparison = NULL) {
    if (missing(annotatedPeak)) {
        stop("Missing required argument 'annotatedPeak'!", call. = FALSE)
    }
    # if (missing(feature_id_type)) {
    #     stop("Missing required argument feature_id_type!")
    # }
    if (length(multiAdjMethod) > 0L) {
        multiAdjMethod <- match.arg(multiAdjMethod,
                                    c("Bonferroni", "Holm", "Hochberg",
                                      "SidakSS", "SidakSD", "BH", "BY",
                                      "ABH", "TSBH"))
    }
    if (!grepl("^org\\...\\.eg\\.db", orgAnn)) {
        stop("Invalid 'orgAnn' parameter. ",
             "Please refer to ",
             "http://www.bioconductor.org/packages/release/data/annotation/ ",
             "for available org.xx.eg.db packages", call. = FALSE)
    }
    if (!isNamespaceLoaded(orgAnn)) {
        stop("Package '", orgAnn, "' must be loaded. ",
             "Try: library(", orgAnn, ")", call. = FALSE)
    }
    if (missing(pathAnn)) {
        stop("Missing required argument 'pathAnn'. ",
             "pathAnn should be a pathway annotation package (e.g., ",
             "'reactome.db' or 'KEGGREST') with objects that map ",
             "Entrez Gene IDs to pathway identifiers (xxxxxEXTID2PATHID) ",
             "and pathway identifiers to pathway names (xxxxxPATHID2NAME).",
             call. = FALSE)
    }
    if (!isNamespaceLoaded(pathAnn)) {
        stop("Package '", pathAnn, "' must be loaded. ",
             "Try: library(", pathAnn, ")", call. = FALSE)
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
        if (length(groupFALSE) > 0L) {
            feature_ids_FALSE <- unique(as.character(groupFALSE$feature))
        }
    } else if (is.character(annotatedPeak)) {
        feature_ids <- unique(annotatedPeak)
        if (length(groupFALSE) > 0L) {
            feature_ids_FALSE <- unique(groupFALSE)
        }
    } else {
        stop("'annotatedPeak' must be a GRanges object with a 'feature' ",
             "column or a character vector of feature IDs", call. = FALSE)
    }
    if (feature_id_type == "entrez_id") {
        entrezIDs <- feature_ids
        if (length(groupFALSE) > 0L) {
            entrezIDs_FALSE <- feature_ids_FALSE
        }
    } else {
        entrezIDs <- convert2EntrezID(feature_ids, orgAnn, feature_id_type)
        if (length(groupFALSE) > 0L) {
            entrezIDs_FALSE <- convert2EntrezID(feature_ids_FALSE,
                                                orgAnn, feature_id_type)
        }
    }
    if (length(entrezIDs) < 2L) {
        stop("The number of genes is less than 2. ",
             "Please double check your 'feature_id_type'. ",
             "If using 'ensembl_gene_id', do not include the version number. ",
             "E.g., 'ENSG00000139618' is okay while 'ENSG00000139618.2' is not. ",
             "Use `ensembl_gene_id <- sub('\\\\.[^.]+$', '', ensembl_gene_id)` ",
             "to batch remove the version number if needed.", call. = FALSE)
    }
    if (pathAnn == "KEGG.db") {
        stop("'pathAnn = \"KEGG.db\"' has been deprecated. ",
             "Use 'pathAnn = \"KEGGREST\"' instead.", call. = FALSE)
    }
    if (pathAnn %in% c("reactome.db", "KEGG.db")) {
        extid2path <- paste(gsub(".db$", "", pathAnn), "EXTID2PATHID", sep = "")
        path2name <- paste(gsub(".db$", "", pathAnn), "PATHID2NAME", sep = "")
        pkg_ns <- paste("package", pathAnn, sep = ":")
        if (length(objects(pkg_ns, pattern = extid2path)) != 1L ||
            length(objects(pkg_ns, pattern = path2name)) != 1L) {
            stop("'pathAnn' is not a valid annotation package with objects ",
                 "named as xxxxxEXTID2PATHID and/or xxxxxPATHID2NAME", 
                 call. = FALSE)
        }
        extid2path <- get(extid2path)
        mapped_genes <- mappedkeys(extid2path)
        # Get all the entrez_ids in the species
        org_symbol <- get(paste(gsub(".db", "", orgAnn), "SYMBOL", sep = ""))
        mapped_genes <- mapped_genes[mapped_genes %in% mappedkeys(org_symbol)]
        xx <- as.list(extid2path[mapped_genes])
    } else if (pathAnn == "KEGGREST") {
        # for KEGGREST db
        organismKEGGREST <- .findKEGGRESTOrganismName(orgAnn)
        EGID2PATHID <- keggLink("pathway", organismKEGGREST)
        # get rid of the leading organismKEGGREST in front of the EntrezID and the "path" in front of the PATHID
        names(EGID2PATHID) <- unlist(lapply(strsplit(names(EGID2PATHID), ":"), "[", 2))
        EGID2PATHID <- unlist(lapply(strsplit(EGID2PATHID, ":"), "[", 2))

        mapped_genes2 <- names(keggLink("pathway", organismKEGGREST))
        mapped_genes2 <- unique(unlist(lapply(strsplit(mapped_genes2, ":"), "[", 2)))

        org_symbol <- get(paste(gsub(".db", "", orgAnn), "SYMBOL", sep = ""))
        mapped_genes <- mapped_genes2[mapped_genes2 %in% mappedkeys(org_symbol)]

        xx <- with(stack(as.list(EGID2PATHID)), split(values, ind))
    }

    all.PATH <- do.call(rbind, lapply(mapped_genes, function(x1) {
        temp <- unlist(xx[names(xx) == x1])
        if (length(temp) > 0L) {
            temp1 <- matrix(temp, ncol = 1L, byrow = TRUE)
            cbind(temp1, rep(x1, nrow(temp1)))
        }
    }))
    this.PATH <- do.call(rbind, lapply(entrezIDs, function(x1) {
        temp <- unlist(xx[names(xx) == x1])
        if (length(temp) > 0L) {
            temp1 <- matrix(temp, ncol = 1L, byrow = TRUE)
            cbind(temp1, rep(x1, nrow(temp1)))
        }
    }))


    all.PATH <- unique(all.PATH)  # In case the database is not unique
    this.PATH <- unique(this.PATH)
    if (is.null(all.PATH) || is.null(this.PATH)) {
        return(data.frame(log = "No enriched pathway found!"))
    }
    colnames(all.PATH) <- c("path.id", "entrez_id")
    colnames(this.PATH) <- c("path.id", "entrez_id")
    path.all <- as.character(all.PATH[, "path.id"])
    path.this <- as.character(this.PATH[, "path.id"])

    total <- length(path.all)
    this <- length(path.this)

    all.count <- getUniqueGOidCount(as.character(path.all[path.all != ""]))
    this.count <- getUniqueGOidCount(as.character(path.this[path.this != ""]))

    if (length(groupFALSE) > 0L) {
        FALSE.PATH <- do.call(rbind, lapply(entrezIDs_FALSE, function(x1) {
            temp <- unlist(xx[names(xx) == x1])
            if (length(temp) > 0L) {
                temp1 <- matrix(temp, ncol = 1L, byrow = TRUE)
                cbind(temp1, rep(x1, nrow(temp1)))
            }
        }))
        FALSE.PATH <- unique(FALSE.PATH)
        colnames(FALSE.PATH) <- c("path.id", "entrez_id")
        path.FALSE <- as.character(FALSE.PATH[, "path.id"])
        FALSE.count <- getUniqueGOidCount(as.character(path.FALSE[path.FALSE != ""]))
        names(FALSE.count[[2L]]) <- FALSE.count[[1L]]
        FALSE.count <- FALSE.count[[2L]]
    }

    selected <- hyperGtest(all.count, this.count, total, this)

    selected <- data.frame(selected)

    colnames(selected) <- c("path.id", "count.InDataset", "count.InGenome",
                            "pvalue", "totaltermInDataset", "totaltermInGenome")
    
    annoTerms <- function(termids) {
        if (length(termids) < 1L) {
            goterm <- matrix(ncol = 2L)
        } else {
            if (pathAnn == "reactome.db") {
                termids <- as.character(termids)
                path2name_obj <- get(sub(".db", "PATHID2NAME", pathAnn))
                terms <- xget(termids, path2name_obj)
                goterm <- cbind(termids,
                                terms[match(termids, names(terms))])
            } else if (pathAnn %in% c("KEGG.db", "KEGGREST")) {
                # Must get rid of the leading organismKEGGREST of path.id
                # if using KEGG.db or KEGGREST
                organismKEGGREST <- .findKEGGRESTOrganismName(orgAnn)
                termidsPreRemoved <- sub(organismKEGGREST, "", termids)

                if (pathAnn == "KEGG.db") {
                    path2name_obj <- get(sub(".db", "PATHID2NAME", pathAnn))
                    terms <- xget(termidsPreRemoved, path2name_obj)
                    goterm <- cbind(termids,
                                    terms[match(termidsPreRemoved, names(terms))])
                } else if (pathAnn == "KEGGREST") {
                    getPathName <- function(pathid) {
                        namePath <- tryCatch(
                            {
                                keggGet(pathid)[[1L]]$NAME
                            },
                            error = function(c) {
                                "NA"
                            }
                        )
                        names(namePath) <- pathid
                        namePath
                    }
                    terms <- unlist(lapply(termids, getPathName))
                    goterm <- cbind(termids,
                                    terms[match(termids, names(terms))])
                }
            }
        }
        colnames(goterm) <- c("path.id", "path.term")
        rownames(goterm) <- NULL
        goterm
    }
    annoterm <- annoTerms(selected$path.id)
    selected <- merge(annoterm, selected, by = "path.id")

    if (is.null(multiAdjMethod)) {
        s <- selected[as.numeric(as.character(selected[, "pvalue"])) < maxP &
                      as.numeric(as.character(selected[, "count.InGenome"])) >= minPATHterm, ]
    } else {
        procs <- c(multiAdjMethod)
        res <- mt.rawp2adjp(as.numeric(as.character(selected[, "pvalue"])), procs)
        adjp <- unique(res$adjp)
        colnames(adjp)[1L] <- "pvalue"
        colnames(adjp)[2L] <- paste(multiAdjMethod, "adjusted.p.value", sep = ".")
        selected[, "pvalue"] <- as.numeric(as.character(selected[, "pvalue"]))
        bp1 <- merge(selected, adjp, all.x = TRUE)

        s <- bp1[as.numeric(as.character(bp1[, ncol(bp1)])) < maxP &
                 !is.na(bp1[, ncol(bp1)]) &&
                 as.numeric(as.character(bp1[, "count.InGenome"])) >= minPATHterm, ]
    }

    selected <- merge(this.PATH, s)

    if (length(groupFALSE) > 0L) {  # Add counts for FALSE group and do Fisher's exact test
        selected$count.InBackgroundDataset <- FALSE.count[selected$path.id]

        multiFisher <- function(.data) {
            rowFisher <- function(x, ...) {
                return(fisher.test(matrix(x, nrow = 2L, byrow = TRUE), ...)$p.value)
            }
            .data0 <- .data[, c("count.InDataset", "count.InBackgroundDataset")]
            .data0 <- as.matrix(.data0)
            .data0[is.na(.data0)] <- 0L
            .data0 <- cbind(.data0, .data[, "count.InGenome"] - .data0)
            apply(.data0, 1L, rowFisher, alternative = "greater")
        }
        selected$dataset.vs.background.pval <- multiFisher(selected)

        selected <- selected[selected$dataset.vs.background.pval < maxP, , drop = FALSE]
    }

    selected
}


.findKEGGRESTOrganismName <- function(name) {
    # Find the organism names used in the KEGGREST db, tested for Hs, Mm, Dm, Rn, Sc, and Dr
    # acroName <- egOrgMap(orgName)
    # KEGGOrganismList <- keggList("organism")
    # ad <- adist(acroName, KEGGOrganismList[, "species"])[1, ]
    # KEGGOrganism <- KEGGOrganismList[which.min(ad), "organism"][[1]]
    if (!is.character(name)) {
        stop("Input organism name must be a character string", call. = FALSE)
    }
    organism <- c(
        "org.Ag.eg.db" = "aga",
        "org.At.eg.db" = "ath",
        "org.Bt.eg.db" = "bta",
        "org.Ce.eg.db" = "cel",
        "org.Cf.eg.db" = "cfa",
        "org.Dm.eg.db" = "dme",
        "org.Dr.eg.db" = "dre",
        "org.EcK12.eg.db" = "eco",
        "org.EcSakai.eg.db" = "ecs",
        "org.Gg.eg.db" = "gga",
        "org.Hs.eg.db" = "hsa",
        "org.Mm.eg.db" = "mmu",
        "org.Mmu.eg.db" = "mcc",
        "org.Pf.plasmo.db" = "pfa",
        "org.Pt.eg.db" = "ptr",
        "org.Rn.eg.db" = "rno",
        "org.Sc.sgd.db" = "sce",
        "org.Sco.eg.db" = "sco",
        "org.Ss.eg.db" = "ssc",
        "org.Tgondii.eg.db" = "tgo",
        "org.Xl.eg.db" = "xla")
    if (name %in% names(organism)) {
        return(organism[name])
    } else {
        if (name %in% organism) {
            return(names(organism)[organism == name])
        }
    }
    # Calculate the string distance, need utils package
    org <- as.character(c(names(organism), organism))
    dis <- adist(name, org)
    org <- org[dis == min(dis)]
    stop("You input \"", name, "\". Do you mean \"", org, "\"?", call. = FALSE)
}
