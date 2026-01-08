#' Calculate coverage of transcript regions per bin
#' 
#' @description 
#' Calculates the average coverage of transcript regions (5' UTR, CDS, 3' UTR,
#' and optionally upstream/downstream regions) per bin across multiple samples.
#' The function bins transcript regions into equal-sized windows and aggregates
#' coverage from RleList objects (e.g., from BigWig files). This is useful for
#' visualizing signal distribution across transcript structures in ChIP-seq or
#' RNA-seq data, creating metagene plots that show how signal varies across
#' different transcript regions.
#' 
#' The function processes transcripts that contain all three core regions (5' UTR,
#' CDS, and 3' UTR), filters them by length criteria, and optionally includes
#' intronic regions or extends analysis to upstream/downstream regions. Coverage
#' is averaged across all qualifying transcripts for each bin within each region
#' type.
#' 
#' @param cvglists A list of \link[IRanges:AtomicList-class]{SimpleRleList} or
#'        \link[IRanges:AtomicList-class]{RleList} objects, where each element
#'        represents coverage data for one sample. Coverage data can be imported
#'        from BigWig files using \code{rtracklayer::import()}. If a single
#'        RleList is provided (not in a list), it will be automatically wrapped
#'        in a list. All RleList objects must have the same chromosome names
#'        (seqlevels) that match the TxDb object.
#' @param TxDb An object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#'        containing transcript annotations. The function extracts gene models
#'        using \code{toGRanges(TxDb, feature="geneModel")} to obtain 5' UTR,
#'        CDS, and 3' UTR regions.
#' @param upstream.cutoff An integer specifying the length (in base pairs) of
#'        the upstream region to include relative to the transcript start site.
#'        Default is \code{1000L}. If set to 0 or negative, upstream regions are
#'        not included. The upstream region is defined as the region immediately
#'        upstream of the transcript start site (TSS).
#' @param downstream.cutoff An integer specifying the length (in base pairs) of
#'        the downstream region to include relative to the transcript end site.
#'        Default is \code{upstream.cutoff}. If set to 0 or negative, downstream
#'        regions are not included. The downstream region is defined as the region
#'        immediately downstream of the transcript end site.
#' @param nbinsCDS An integer specifying the number of bins to divide CDS regions
#'        into. Default is \code{100L}. Each CDS region will be divided into this
#'        many equal-sized bins, and coverage will be averaged within each bin.
#' @param nbinsUTR An integer specifying the number of bins to divide UTR regions
#'        (both 5' and 3' UTR) into. Default is \code{20L}. Each UTR region will
#'        be divided into this many equal-sized bins.
#' @param nbinsUpstream An integer specifying the number of bins to divide
#'        upstream regions into. Default is \code{20L}. Only used if
#'        \code{upstream.cutoff > 0}.
#' @param nbinsDownstream An integer specifying the number of bins to divide
#'        downstream regions into. Default is \code{nbinsUpstream}. Only used if
#'        \code{downstream.cutoff > 0}.
#' @param includeIntron A logical value. When \code{TRUE}, intronic regions
#'        within CDS are included in the analysis. When \code{FALSE} (default),
#'        only exonic CDS regions are used. If \code{TRUE}, the function merges
#'        CDS regions from the same transcript on the same chromosome, effectively
#'        including introns between CDS exons.
#' @param minCDSLen An integer specifying the minimum total CDS length (in base
#'        pairs) required for a transcript to be included. Default is
#'        \code{nbinsCDS}. Transcripts with CDS shorter than this value are
#'        filtered out. Must be >= \code{nbinsCDS}.
#' @param minUTRLen An integer specifying the minimum total UTR length (in base
#'        pairs) required for both 5' and 3' UTR regions. Default is
#'        \code{nbinsUTR}. Transcripts with either UTR shorter than this value are
#'        filtered out. Must be >= \code{nbinsUTR}.
#' @param maxCDSLen A numeric value specifying the maximum total CDS length (in
#'        base pairs) allowed for a transcript to be included. Default is
#'        \code{Inf} (no maximum). Transcripts with CDS longer than this value are
#'        filtered out. Must be > \code{minCDSLen}.
#' @param maxUTRLen A numeric value specifying the maximum total UTR length (in
#'        base pairs) allowed for both 5' and 3' UTR regions. Default is
#'        \code{Inf} (no maximum). Transcripts with UTRs longer than this value
#'        are filtered out. Must be > \code{minUTRLen}.
#' 
#' @details
#' 
#' \strong{Transcript filtering:}
#' \itemize{
#'   \item Only transcripts containing all three core regions (5' UTR, CDS, and
#'         3' UTR) are included in the analysis
#'   \item Transcripts are filtered by length criteria (\code{minCDSLen},
#'         \code{maxCDSLen}, \code{minUTRLen}, \code{maxUTRLen})
#'   \item Upstream and downstream regions are filtered to ensure they meet
#'         minimum bin requirements (\code{nbinsUpstream}, \code{nbinsDownstream})
#'   \item At least 2 transcripts must remain after filtering, otherwise an error
#'         is raised
#' }
#' 
#' \strong{Region processing:}
#' \itemize{
#'   \item Transcript regions are disjoined to remove overlaps between different
#'         feature types
#'   \item For each feature type, regions are divided into equal-sized bins
#'   \item Coverage is calculated per bin using \code{viewMeans()} on Views objects
#'   \item Coverage is averaged across all transcripts for each bin
#'   \item For negative strand transcripts, bins are reversed to maintain 5'->3'
#'         orientation
#' }
#' 
#' \strong{Upstream and downstream regions:}
#' \itemize{
#'   \item Upstream regions are defined using \code{promoters()} function with
#'         \code{upstream = upstream.cutoff} and \code{downstream = 0}
#'   \item Downstream regions are defined by reversing strand, applying
#'         \code{promoters()}, then reversing strand back
#'   \item Overlapping parts between upstream/downstream and gene body regions
#'         are removed to avoid double-counting
#' }
#' 
#' @return Returns a named list of matrices, where each element corresponds to
#'        a feature type (e.g., "upstream", "5UTR", "CDS", "3UTR", "downstream").
#'        Each matrix has:
#'        \itemize{
#'          \item Rows: Bins within the feature type (e.g., 100 rows for CDS if
#'                \code{nbinsCDS = 100})
#'          \item Columns: Samples (one column per element in \code{cvglists})
#'          \item Values: Average coverage per bin, averaged across all qualifying
#'                transcripts
#'        }
#'        The list is ordered as: upstream (if included), 5UTR, CDS, 3UTR,
#'        downstream (if included). This structure is suitable for plotting with
#'        \code{\link{plotBinOverRegions}}.
#' @author Jianhong Ou
#' @seealso \code{\link{binOverGene}} for gene-level binning (simpler approach),
#'          \code{\link{plotBinOverRegions}} for visualizing the output,
#'          \code{\link[rtracklayer]{import}} for importing BigWig files as RleList
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom GenomicFeatures transcripts
#' @importFrom S4Vectors mcols elementMetadata
#' @importFrom BiocGenerics strand start
#' @importFrom GenomeInfoDb seqinfo seqlevels
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
#' d <- binOverRegions(cvglists, TxDb.Hsapiens.UCSC.hg19.knownGene)
#' plotBinOverRegions(d)
#' }
#' }
#' 
binOverRegions <- function(cvglists, TxDb, 
                           upstream.cutoff = 1000L, 
                           downstream.cutoff = upstream.cutoff, 
                           nbinsCDS = 100L, nbinsUTR = 20L, 
                           nbinsUpstream = 20L,
                           nbinsDownstream = nbinsUpstream,
                           includeIntron = FALSE,
                           minCDSLen = nbinsCDS,
                           minUTRLen = nbinsUTR,
                           maxCDSLen = Inf,
                           maxUTRLen = Inf) {
    if (inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))) {
        cvglistsName <- substitute(deparse(cvglists))
        cvglists <- list(cvglists)
        names(cvglists) <- cvglistsName
    }
    if (!is.list(cvglists)) {
        stop("cvglists must be a list of SimpleRleList or RleList", 
             call. = FALSE)
    }
    cls <- sapply(cvglists, inherits, 
                  what = c("SimpleRleList", "RleList", "CompressedRleList"))
    if (any(!cls))
        stop("cvglists must be a list of SimpleRleList or RleList")
    stopifnot(inherits(TxDb, "TxDb"))
    stopifnot(maxCDSLen > minCDSLen)
    stopifnot(minCDSLen >= nbinsCDS)
    stopifnot(maxUTRLen > minUTRLen)
    stopifnot(minUTRLen >= nbinsUTR)
    features <- toGRanges(TxDb, feature = "geneModel")
    features.reduced <- reduce(features)
    feature_type <- c("upstream", "5UTR", "CDS", "3UTR", "downstream")
    features <- features[features$feature_type %in% c("5UTR", "CDS", "3UTR")]
    transcripts <- transcripts(TxDb, columns = c("tx_name", "gene_id"))
    trx <- unique(mcols(transcripts))
    features$gene_id <- trx[match(features$tx_name, trx$tx_name), "gene_id"]
    features <- features[lengths(features$gene_id) > 0]
    features$tx_name <- sapply(features$gene_id, function(.ele) .ele[1])
    features$gene_id <- NULL
    
    features.disjoin <- disjoin(features, with.revmap = TRUE)
    features.disjoin.1 <- 
        features.disjoin[rep(seq_along(features.disjoin),
                             lengths(features.disjoin$revmap))]
    features.disjoin.1$revmap <- unlist(features.disjoin$revmap)
    features.disjoin.1$feature_type <- 
        features$feature_type[features.disjoin.1$revmap]
    features.disjoin.1$tx_name <- features$tx_name[features.disjoin.1$revmap]
    features.disjoin.1$revmap2 <- 
        rep(seq_along(features.disjoin), lengths(features.disjoin$revmap))
    feature_type2 <- 
        split(features.disjoin.1$feature_type, features.disjoin.1$revmap2)
    feature_type2 <- lapply(feature_type2, unique)
    feature_type2 <- feature_type2[lengths(feature_type2) == 1]
    features.disjoin.1 <- 
      features.disjoin.1[match(as.numeric(names(feature_type2)),
                               features.disjoin.1$revmap2)]
    features.disjoin.1$seqn <- paste(features.disjoin.1$tx_name, 
                                     features.disjoin.1$feature_type, 
                                     as.character(seqnames(features.disjoin.1)))
    features1 <- GRanges(seqnames = features.disjoin.1$seqn, 
                         ranges = ranges(features.disjoin.1),
                         strand = strand(features.disjoin.1))
    features1 <- reduce(features1)
    mcols(features1) <- 
        data.frame(do.call(rbind, 
                           strsplit(as.character(seqnames(features1)), " ")), 
                   stringsAsFactors = FALSE)
    colnames(mcols(features1)) <- c("tx_name", "feature_type", "seqn")
    features1 <- GRanges(seqnames = features1$seqn, 
                         ranges = ranges(features1),
                         strand = strand(features1),
                         tx_name = features1$tx_name,
                         feature_type = features1$feature_type)
    seqinfo(features1) <- seqinfo(features)[seqlevels(features1)]
    features <- features1
    rm(list = c("features.disjoin", "features.disjoin.1", 
                "feature_type2", "features1"))
    
    txs <- unique(features$tx_name)
    txs <- txs[!is.na(txs)]
    if (includeIntron) {
        cds <- features[features$feature_type == "CDS"]
        cds <- GRanges(seqnames = paste(as.character(seqnames(cds)), 
                                        cds$tx_name),
                       ranges = ranges(cds), strand = strand(cds))
        cds <- reduce(cds, min.gapwidth = 1e9)
        mcols(cds) <- 
            data.frame(do.call(rbind, 
                               strsplit(as.character(seqnames(cds)), " ")), 
                       stringsAsFactors = FALSE)
        colnames(mcols(cds)) <- c("seqn", "tx_name")
        cds <- GRanges(seqnames = cds$seqn, 
                       ranges = ranges(cds), strand = strand(cds),
                       tx_name = cds$tx_name, feature_type = "CDS")
        seqinfo(cds) <- seqinfo(features)[seqlevels(cds)]
        features <- features[features$feature_type %in% c("5UTR", "3UTR")]
        features <- c(features, cds)
    }
    ## resort by txs
    features <- features[order(as.numeric(factor(features$tx_name, 
                                                 levels = txs)))]
    ## features must contain 5UTR, CDS and 3UTR
    txs <- split(features$feature_type, features$tx_name)
    txs <- lapply(txs, unique)
    txs <- names(txs)[lengths(txs) == 3]
    features <- features[features$tx_name %in% txs]
    ## add upstream and downstream
    if (upstream.cutoff > 0 && downstream.cutoff > 0) {
        genes <- GRanges(seqnames = paste(as.character(seqnames(features)), 
                                          features$tx_name),
                         ranges = ranges(features), strand = strand(features))
        genes <- reduce(genes, min.gapwidth = 1e9)
        mcols(genes) <- 
          data.frame(do.call(rbind, 
                             strsplit(as.character(seqnames(genes)), " ")),
                     stringsAsFactors = FALSE)
        colnames(mcols(genes)) <- c("seqn", "tx_name")
        genes <- 
          GRanges(seqnames = genes$seqn, 
                  ranges = ranges(genes), strand = strand(genes),
                  tx_name = genes$tx_name)
        seqinfo(genes) <- seqinfo(features)[seqlevels(genes)]
        suppressWarnings({
            upstream <- promoters(genes,
                                  upstream = upstream.cutoff,
                                  downstream = 0)})
        upstream <- trim(upstream)
        upstream <- upstream[width(upstream) >= nbinsUpstream]
        revert_strand <- function(gr) {
            l <- as.character(strand(gr))
            l_plus <- which(l == "+")
            l_minus <- which(l == "-")
            if (length(l_plus) >= 1L) {
                l[l_plus] <- "-"
            }
            if (length(l_minus) >= 1L) {
                l[l_minus] <- "+"
            }
            strand(gr) <- l
            gr
        }
        genes_rev_strand <- revert_strand(genes)
        suppressWarnings({
            downstream <- 
                promoters(genes_rev_strand,
                         upstream = downstream.cutoff,
                         downstream = 0L)
        })
        downstream <- revert_strand(downstream)
        downstream <- trim(downstream)
        downstream <- downstream[width(downstream) >= nbinsDownstream]
        upstream$feature_type <- "upstream"
        downstream$feature_type <- "downstream"
        #remove the overlapping parts
        rmOverlaps <- function(gr, tobeRemoved = features.reduced) {
            gr.intersect <- intersect(gr, tobeRemoved, ignore.strand = TRUE)
            gr.intersect$tx_name <- NA
            gr.intersect$feature_type <- NA
            gr.disjoin <- disjoin(c(gr, gr.intersect), 
                                  with.revmap = TRUE, ignore.strand = TRUE)
            gr.disjoin <- gr.disjoin[lengths(gr.disjoin$revmap) == 1]
            gr.disjoin$revmap <- unlist(gr.disjoin$revmap)
            strand(gr.disjoin) <- strand(gr[gr.disjoin$revmap])
            mcols(gr.disjoin) <- mcols(gr[gr.disjoin$revmap])
            gr.disjoin
        }
        upstream <- rmOverlaps(upstream)
        downstream <- rmOverlaps(downstream)
        features <- c(features, upstream, downstream)
    }
    ## split the features by seqnames and generate Views for cvglist
    seqn <- Reduce(intersect, lapply(cvglists, names))
    seqn <- intersect(seqn, seqlevels(features))
    if (length(seqn) < 1L) {
        stop("Please check the names of cvglist. ",
             "None of them in the seqlevels of TxDb.", call. = FALSE)
    }
    features <- features[seqnames(features) %in% seqn]
    
    ## make sure no overlaps in all the region 
    features.disjoin <- disjoin(features, with.revmap = TRUE)
    features.disjoin <- features.disjoin[lengths(features.disjoin$revmap) == 1]
    features.disjoin$revmap <- unlist(features.disjoin$revmap)
    features.disjoin$feature_type <- 
      features$feature_type[features.disjoin$revmap]
    features.disjoin$tx_name <- features$tx_name[features.disjoin$revmap]
    features.disjoin$revmap <- NULL
    features.disjoin.tx_name <- 
      split(features.disjoin$feature_type, features.disjoin$tx_name)
    features.disjoin.tx_name <- lapply(features.disjoin.tx_name, unique)
    features.disjoin.tx_name <- 
        names(features.disjoin.tx_name)[lengths(features.disjoin.tx_name) ==
                                        max(lengths(features.disjoin.tx_name))]
    features.disjoin <- features.disjoin[features.disjoin$tx_name %in% 
                                           features.disjoin.tx_name]
    seqinfo(features.disjoin) <- seqinfo(features)[seqlevels(features.disjoin)]
    features <- features.disjoin
    rm(features.disjoin)
    ## make sure features are sorted by pos
    features <- 
        features[order(
            as.numeric(factor(features$tx_name, 
                            levels = txs)),
            as.numeric(factor(features$feature_type,
                            levels = 
                                feature_type[feature_type %in% 
                                           unique(
                                               features$feature_type)])),
            start(features))]
    
    ## filter by minCDSLen, min5UTR, min3UTR
    filterByLen <- function(type, minlen, maxLen) {
        x <- features[features$feature_type %in% type]
        if (length(x) < 1) return(x)
        x.width <- rowsum(width(x), x$tx_name, reorder = FALSE)
        rownames(x.width)[x.width[, 1] >= minlen & x.width[, 1] < maxLen]
    }
    txs <- Reduce(intersect, mapply(filterByLen, 
                                    c("CDS", "5UTR", "3UTR", 
                                      "upstream", "downstream"), 
                                    c(minCDSLen, minUTRLen, minUTRLen, 
                                      nbinsUpstream, nbinsDownstream),
                                    c(maxCDSLen, maxUTRLen, maxUTRLen,
                                      Inf, Inf)))
    if (length(txs) < 2L) {
        stop("Less than 2 transcripts remaining after filtering", 
             call. = FALSE)
    }
    features <- features[features$tx_name %in% txs]
    ## features must contain 5UTR, CDS, 3UTR, and/or upstream, downstream.
    txs <- split(features$feature_type, features$tx_name)
    txs <- lapply(txs, unique)
    txs <- names(txs)[lengths(txs) == length(unique(features$feature_type))]
    if (length(txs) < 2L) {
        stop("Less than 2 transcripts remaining after filtering", 
             call. = FALSE)
    }
    features <- features[features$tx_name %in% txs]
    ## calculate fators to balance the region be count multiple times
    len <- length(unique(features$tx_name))
    features.s <- split(features, seqnames(features))
    features.s <- features.s[seqn]
    cvglists <- lapply(cvglists, function(.ele) .ele[seqn])
    features.l <- as(features.s, "IntegerRangesList")
    
    bins <- c("upstream" = nbinsUpstream,
              "5UTR" = nbinsUTR,
              "CDS" = nbinsCDS,
              "3UTR" = nbinsUTR,
              "downstream" = nbinsDownstream)
    
    features.view <- lapply(cvglists, function(.ele) {
        vw <- Views(.ele, features.l)
        cntByChr <- lapply(vw, function(.vw) {
            .df <- elementMetadata(.vw)
            .gr <- GRanges(paste(.df$feature_type, .df$tx_name), ranges(.vw), 
                           feature_type = .df$feature_type, tx_name = .df$tx_name)
            .gr.s <- split(.gr, seqnames(.gr))
            .cvg <- subject(.vw)
            .cvgs <- rep(list(.cvg), length(.gr.s))
            .gr.l <- as(.gr.s, "IntegerRangesList")
            names(.cvgs) <- names(.gr.l)
            .cvgs <- as(.cvgs, "SimpleRleList")
            .cvg.sub <- .cvgs[.gr.l]
            ### split the RleList into bins
            .ir <- IRanges(rep(1, length(.cvg.sub)), lengths(.cvg.sub))
            .ir <- IRanges::tile(.ir, n = bins[sub("^(.*?) .*$", "\\1", 
                                                   names(.cvg.sub))])
            names(.ir) <- names(.cvg.sub)
            .cnt <- viewMeans(Views(.cvg.sub, .ir))
            .cnt <- split(.cnt, sub("^(.*?) .*$", "\\1", names(.cnt)))
            ## make sure 5'->3'
            .cnt <- lapply(.cnt, function(x) {
                strd <- 
                    as.character(strand(features))[match(sub("^(.*?) (.*$)", "\\2", 
                                                           names(x)), 
                                                       features$tx_name)]
                x <- split(x, strd)
                x <- lapply(x, do.call, what = rbind)
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
    d <- d[feature_type[feature_type %in% names(d)]]
    d <- lapply(d, do.call, what = cbind)
    return(d)
}




#' Plot coverage of transcript regions per bin
#' 
#' @description 
#' Creates a line plot visualizing the average coverage across transcript regions
#' (upstream, 5' UTR, CDS, 3' UTR, downstream) per bin. This function is designed
#' to visualize the output from \code{\link{binOverRegions}} or \code{\link{binOverGene}},
#' creating metagene plots that show signal distribution across transcript structures.
#' 
#' The plot displays:
#' \itemize{
#'   \item Coverage values on the y-axis (averaged across transcripts)
#'   \item Bin positions on the x-axis, with vertical lines separating different
#'         region types
#'   \item One line per sample (from the columns of input matrices)
#'   \item Region type labels centered in each region section
#'   \item A legend showing sample names
#' }
#' 
#' @param dat A named list of matrices, typically the output from
#'        \code{\link{binOverRegions}} or \code{\link{binOverGene}}. Each element
#'        of the list should be a matrix where:
#'        \itemize{
#'          \item Rows represent bins within a region type
#'          \item Columns represent different samples
#'          \item Values are average coverage per bin
#'        }
#'        The list names must be a subset of: "upstream", "5UTR", "CDS", "gene",
#'        "3UTR", "downstream". The function will plot regions in this order.
#' @param ... Additional parameters passed to \code{\link[graphics]{matplot}} and
#'        \code{\link[graphics]{legend}}, such as:
#'        \itemize{
#'          \item \code{col}: Colors for each sample line
#'          \item \code{lty}: Line types for each sample
#'          \item \code{lwd}: Line widths
#'          \item \code{xlab}, \code{ylab}: Axis labels
#'          \item \code{main}: Plot title
#'          \item \code{cex}: Character expansion for legend text
#'          \item Other standard graphics parameters
#'        }
#' 
#' @details
#' The function automatically:
#' \itemize{
#'   \item Orders regions in the standard order: upstream, 5UTR, CDS, gene (if
#'         present), 3UTR, downstream
#'   \item Combines all region matrices into a single matrix for plotting
#'   \item Adds vertical dashed lines at region boundaries
#'   \item Centers region type labels within each region section
#'   \item Creates a legend in the top-right corner with sample names
#'   \item Uses default colors, line types, and line widths from \code{matplot}
#'         unless overridden via \code{...}
#' }
#' 
#' @return Returns \code{invisible(NULL)}. The function is called for its side
#'        effect of creating a plot.
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{binOverRegions}} for calculating coverage per bin,
#'          \code{\link{binOverGene}} for gene-level binning,
#'          \code{\link[graphics]{matplot}} for the underlying plotting function
#' @export
#' @importFrom S4Vectors elementNROWS
#' @importFrom graphics abline axis legend matplot
#' @examples
#' 
#' if(interactive()){
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
plotBinOverRegions <- function(dat, ...) {
    feature_type <- c("upstream", "5UTR", "CDS", "gene", "3UTR", "downstream")
    if (!all(names(dat) %in% feature_type)) {
        stop("names of dat must be upstream, 5UTR, CDS, 3UTR or downstream")
    }
    feature_type <- feature_type[feature_type %in% names(dat)]
    dat <- dat[feature_type]
    bins <- elementNROWS(dat)
    dat <- do.call(rbind, dat)
    matplot(dat, type = "l", xaxt = "n", ...)
    bins.sum <- cumsum(bins) + .5
    bins.sum <- bins.sum[-length(bins.sum)]
    abline(v = bins.sum, lty = 2, col = "gray")
    axis(side = 1, at = bins.sum, 
         labels = rep("", length(bins.sum)), 
         tick = TRUE, lwd = -1, lwd.ticks = 1, ...)
    bins.sum2 <- c(0, bins.sum) + bins / 2
    axis(side = 1, at = bins.sum2, 
         labels = feature_type, lwd = -1, ...)
    legend.arg <- formals("legend")
    legend.arg.names <- names(legend.arg)
    legend.arg$x <- "topright"
    legend.arg$box.col <- NA
    legend.arg$merge <- FALSE
    legend.arg$legend <- colnames(dat)
    matplot.arg <- formals("matplot")
    for (arg in c("col", "lty", "lwd", "bg")) {
        legend.arg[[arg]] <- matplot.arg[[arg]]
    }
    legend.arg <- c(list(...), legend.arg)
    legend.arg <- legend.arg[!duplicated(names(legend.arg))]
    legend.arg <- legend.arg[legend.arg.names]
    legend.arg <- legend.arg[lengths(legend.arg) > 0]
    if (is.null(legend.arg[["cex"]])) {
        legend.arg[["cex"]] <- 1
    }
    do.call(legend, legend.arg)
}

