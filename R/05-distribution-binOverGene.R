#' Calculate coverage of gene body per bin
#' 
#' @description 
#' Calculates the average coverage of gene bodies, upstream, and downstream
#' regions per bin across multiple samples. The function bins genes into
#' equal-sized windows and aggregates coverage from RleList objects (e.g.,
#' from BigWig files). Useful for visualizing signal distribution around
#' genes in ChIP-seq or RNA-seq data.
#' 
#' The function:
#' \itemize{
#'   \item Extracts gene features from a TxDb object (transcripts or exons)
#'   \item Filters genes by length (minGeneLen to maxGeneLen)
#'   \item Optionally adds upstream and downstream regions
#'   \item Bins each region into equal-sized windows
#'   \item Calculates average coverage per bin across all genes
#'   \item Returns coverage matrices organized by feature type (upstream, gene, downstream)
#' }
#' 
#' @param cvglists A list of \link[IRanges:AtomicList-class]{SimpleRleList} or
#'        \link[IRanges:AtomicList-class]{RleList} objects, where each element
#'        represents coverage data for one sample. Coverage data can be obtained
#'        from BigWig files using \code{rtracklayer::import()}. If a single
#'        RleList is provided, it will be wrapped in a list. The names of the list
#'        will be used as sample names in the output.
#' @param TxDb An object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}} used
#'        for extracting gene annotations. The function extracts either transcripts
#'        or exons depending on the \code{includeIntron} parameter.
#' @param upstream.cutoff An integer specifying the length (in base pairs) of the
#'        upstream region to include for each gene. Default is \code{0L}, which
#'        means no upstream region is included. When > 0, creates upstream regions
#'        using \code{promoters()} function. Only genes with upstream regions
#'        >= \code{nbinsUpstream} are retained.
#' @param downstream.cutoff An integer specifying the length (in base pairs) of the
#'        downstream region to include for each gene. Default is \code{upstream.cutoff}.
#'        When > 0, creates downstream regions by reversing strand and using
#'        \code{promoters()}. Only genes with downstream regions >=
#'        \code{nbinsDownstream} are retained.
#' @param nbinsGene An integer specifying the number of bins to divide each gene
#'        body into. Default is \code{100L}. Each bin will have equal width
#'        (gene_length / nbinsGene). Genes shorter than \code{nbinsGene} are
#'        filtered out.
#' @param nbinsUpstream An integer specifying the number of bins for upstream
#'        regions. Default is \code{20L}. Only used when \code{upstream.cutoff > 0}.
#' @param nbinsDownstream An integer specifying the number of bins for downstream
#'        regions. Default is \code{nbinsUpstream}. Only used when
#'        \code{downstream.cutoff > 0}.
#' @param includeIntron A logical value. When \code{FALSE} (default), uses only
#'        exons (introns are excluded). When \code{TRUE}, uses full transcripts
#'        including introns. This affects how gene length is calculated and which
#'        regions are binned.
#' @param minGeneLen An integer specifying the minimum gene length (in base pairs)
#'        to include. Default is \code{nbinsGene}. Genes shorter than this are
#'        filtered out. Must be >= \code{nbinsGene}.
#' @param maxGeneLen An integer or \code{Inf} specifying the maximum gene length
#'        (in base pairs) to include. Default is \code{Inf} (no upper limit). Genes
#'        longer than this are filtered out. Must be > \code{minGeneLen}.
#' 
#' @return Returns a named list with up to three elements, each containing a
#'        matrix of average coverage values:
#'        \itemize{
#'          \item \code{upstream}: A matrix with \code{nbinsUpstream} rows and
#'                one column per sample. Each value is the average coverage across
#'                all genes for that bin. Only present if \code{upstream.cutoff > 0}.
#'          \item \code{gene}: A matrix with \code{nbinsGene} rows and one column
#'                per sample. Each value is the average coverage across all genes
#'                for that bin. Rows are ordered from 5' to 3' end (strand-aware).
#'          \item \code{downstream}: A matrix with \code{nbinsDownstream} rows and
#'                one column per sample. Each value is the average coverage across
#'                all genes for that bin. Only present if \code{downstream.cutoff > 0}.
#'        }
#'        
#'        The matrices can be directly used with \code{\link{plotBinOverRegions}}
#'        to visualize metagene plots. Column names correspond to sample names
#'        from \code{cvglists}.
#' 
#' @details
#' 
#' \strong{Gene filtering:}
#' \itemize{
#'   \item Genes are filtered by length: \code{minGeneLen <= gene_length < maxGeneLen}
#'   \item Gene length is calculated as the sum of all exon/transcript widths
#'         for that gene (after reducing overlapping regions)
#'   \item At least 2 genes must remain after filtering, otherwise an error is raised
#' }
#' 
#' \strong{Upstream/Downstream regions:}
#' \itemize{
#'   \item Upstream regions are created using \code{promoters(genes, upstream = upstream.cutoff, downstream = 0)}
#'   \item Downstream regions are created by reversing strand, using \code{promoters()},
#'         then reversing strand back
#'   \item Only genes with upstream/downstream regions >= the respective bin counts
#'         are retained
#'   \item Regions are trimmed to chromosome boundaries using \code{trim()}
#' }
#' 
#' \strong{Binning process:}
#' \itemize{
#'   \item Each region (upstream, gene, downstream) is divided into equal-sized bins
#'   \item Coverage is calculated as the mean coverage within each bin
#'   \item For genes on negative strand, bins are reversed to maintain 5'->3' order
#'   \item Final values are averaged across all genes for each bin
#' }
#' 
#' \strong{Chromosome matching:}
#' \itemize{
#'   \item Only chromosomes present in both \code{cvglists} and \code{TxDb} are used
#'   \item If no common chromosomes are found, an error is raised
#'   \item Genes with start positions < 1 are removed
#' }
#' @author Jianhong Ou
#' @seealso \link{binOverRegions}, \link{plotBinOverRegions}
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics start end width strand do.call
#' @importFrom GenomeInfoDb seqinfo
#' @examples
#' 
#' if(Sys.getenv("USER")=="jou"){
#' path <- system.file("extdata", package="ChIPpeakAnno")
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(rtracklayer)
#' files <- dir(path, "bigWig")
#' if(.Platform$OS.type != "windows"){
#' cvglists <- lapply(file.path(path, files), import,
#'                    format="BigWig", as="RleList")
#' names(cvglists) <- sub(".bigWig", "", files)
#' d <- binOverGene(cvglists, TxDb.Hsapiens.UCSC.hg19.knownGene)
#' plotBinOverRegions(d)
#' }
#' }
#' 
binOverGene <- function(cvglists, TxDb, 
                        upstream.cutoff = 0L, 
                        downstream.cutoff = upstream.cutoff, 
                        nbinsGene = 100L,
                        nbinsUpstream = 20L,
                        nbinsDownstream = nbinsUpstream,
                        includeIntron = FALSE,
                        minGeneLen = nbinsGene,
                        maxGeneLen = Inf) {
    if (inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))) {
        cvglistsName <- substitute(deparse(cvglists))
        cvglists <- list(cvglists)
        names(cvglists) <- cvglistsName
    }
    if (!is.list(cvglists)) {
        stop("cvglists must be a list of SimpleRleList or RleList", call. = FALSE)
    }
    cls <- sapply(cvglists, inherits, 
                  what = c("SimpleRleList", "RleList", "CompressedRleList"))
    if (any(!cls)) {
        stop("cvglists must be a list of SimpleRleList or RleList")
    }
    stopifnot(inherits(TxDb, "TxDb"))
    stopifnot(maxGeneLen > minGeneLen)
    stopifnot(minGeneLen >= nbinsGene)
    if (includeIntron) {
        features <- toGRanges(TxDb, feature = "transcript")
    } else {
        features <- toGRanges(TxDb, feature = "exon")
    }
    ## Reduce by gene
    GRapply <- function(X, FUN, by, ...) {
        stopifnot(by != "by_id")
        stopifnot(inherits(X, "GRanges"))
        x <- GRanges(seqnames = mcols(X)[, by], 
                     ranges = ranges(X),
                     strand = strand(X))
        x <- FUN(x, ...)
        x$by_id <- as.character(seqnames(x))
        x <- GRanges(seqnames = 
                        as.character(seqnames(X[match(x$by_id, 
                                                      mcols(X)[, by])])),
                    ranges = ranges(x),
                    strand = strand(x),
                    by_id = x$by_id)
        colnames(mcols(x)) <- by
        seqinfo(x) <- seqinfo(X)[seqlevels(x)]
        x
    }
    features <- features[lengths(features$gene_id) > 0]
    features$gene_id <- sapply(features$gene_id, `[`, 1)
    features <- GRapply(features, reduce, "gene_id")
    
    ## Filter by minGeneLen, maxGeneLen
    features <- trim(features)
    f.width <- rowsum(width(features), features$gene_id, reorder = FALSE)
    features <- features[features$gene_id %in% 
                        rownames(f.width)[f.width[, 1] >= minGeneLen & 
                                         f.width[, 1] < maxGeneLen]]
    
    if (length(features) < 2L) {
        stop("Less than 2 genes remaining after filtering", call. = FALSE)
    }
    ## Add upstream and downstream
    if (upstream.cutoff > 0L && downstream.cutoff > 0L) {
        genes <- GRapply(features, range, "gene_id")
        suppressWarnings({
            upstream <- promoters(genes,
                                  upstream = upstream.cutoff,
                                  downstream = 0)
        })
        upstream <- trim(upstream)
        upstream <- upstream[width(upstream) >= nbinsUpstream]
        revert_strand <- function(gr) {
            l <- levels(strand(gr))
            l_plus <- which(l == "+")
            l_minus <- which(l == "-")
            if (length(l_plus) == 1L) {
                l[l_plus] <- "-"
            }
            if (length(l_minus) == 1L) {
                l[l_minus] <- "+"
            }
            levels(strand(gr)) <- l
            gr
        }
        genes_rev_strand <- revert_strand(genes)
        suppressWarnings({
            downstream <- promoters(genes_rev_strand,
                                   upstream = downstream.cutoff,
                                   downstream = 0L)
        })
        downstream <- revert_strand(downstream)
        downstream <- trim(downstream)
        downstream <- downstream[width(downstream) >= nbinsDownstream]
        upstream$feature_type <- "upstream"
        downstream$feature_type <- "downstream"
        features$feature_type <- "gene"
        gene_id <- unique(features$gene_id)
        features <- c(features, upstream, downstream)
        features <- features[order(as.numeric(factor(features$gene_id, 
                                                     levels = gene_id)),
                                   as.numeric(factor(features$feature_type,
                                                     levels = c("upstream", 
                                                               "gene", 
                                                               "downstream"))),
                                   start(features))]
    } else {
        ## Make sure features are sorted by position
        gene_id <- unique(features$gene_id)
        features$feature_type <- "gene"
        features <- features[order(as.numeric(factor(features$gene_id, 
                                                     levels = gene_id)),
                                   start(features))]
    }
    ## Split the features by seqnames and generate Views for cvglist
    seqn <- Reduce(intersect, lapply(cvglists, names))
    seqn <- intersect(seqn, seqlevels(features))
    if (length(seqn) < 1L) {
        stop("Please check the names of cvglist. ",
             "None of them in the seqlevels of TxDb.", call. = FALSE)
    }
    features <- features[seqnames(features) %in% seqn]
    features <- features[!features$gene_id %in% 
                        features$gene_id[start(features) < 1]] 
    ## Remove the start < 1
    len <- length(unique(features$gene_id))
    features.s <- split(features, seqnames(features))
    features.s <- features.s[seqn]
    cvglists <- lapply(cvglists, function(.ele) .ele[seqn])
    features.l <- as(features.s, "IntegerRangesList")
    
    bins <- c("upstream" = nbinsUpstream,
              "gene" = nbinsGene,
              "downstream" = nbinsDownstream)
    
    features.view <- lapply(cvglists, function(.ele) {
        vw <- Views(.ele, features.l)
        cntByChr <- lapply(vw, function(.vw) {
            .df <- elementMetadata(.vw)
            if (length(.df$feature_type) == 0) {
                .df$feature_type <- "gene"
            }
            .gr <- GRanges(paste(.df$feature_type, .df$gene_id), ranges(.vw), 
                          feature_type = .df$feature_type, 
                          gene_id = .df$gene_id)
            .gr.s <- split(.gr, seqnames(.gr))
            .cvg <- subject(.vw)
            .cvgs <- rep(list(.cvg), length(.gr.s))
            .gr.l <- as(.gr.s, "IntegerRangesList")
            names(.cvgs) <- names(.gr.l)
            .cvgs <- as(.cvgs, "SimpleRleList")
            .cvg.sub <- .cvgs[.gr.l]
            ### Split the RleList into bins
            .ir <- IRanges(1, lengths(.cvg.sub))
            .ir <- IRanges::tile(.ir, 
                                 n = bins[sub("^(.*?) .*$", "\\1", names(.cvg.sub))])
            names(.ir) <- names(.cvg.sub)
            .cnt <- viewMeans(Views(.cvg.sub, .ir))
            .cnt <- split(.cnt, sub("^(.*?) .*$", "\\1", names(.cnt)))
            ## Make sure 5'->3'
            .cnt <- lapply(.cnt, function(x) {
                strd <- as.character(strand(features))[match(sub("^(.*?) (.*$)", 
                                                                "\\2", names(x)), 
                                                            features$gene_id)]
                x <- split(x, strd)
                x <- lapply(x, BiocGenerics::do.call, what = rbind)
                if ("-" %in% names(x)) {
                    x[["-"]] <- x[["-"]][, rev(seq_len(ncol(x[["-"]]))), drop = FALSE]
                }
                do.call(rbind, x)
            })
            lapply(.cnt, colSums)
        })
        cntByFeature <- swapList(cntByChr)
        cntByFeature <- lapply(cntByFeature, function(.ele) {
            colSums(do.call(rbind, .ele)) / len
        })
    })
    d <- swapList(features.view)
    feature_type <- c("upstream", "gene", "downstream")
    d <- d[feature_type[feature_type %in% names(d)]]
    d <- lapply(d, do.call, what = cbind)
    return(d)
}
