#' Extract signals in given genomic ranges
#' 
#' @description 
#' Extracts and aggregates signal coverage (e.g., from BigWig files or ChIP-seq
#' data) across specified genomic features. The function tiles each feature into
#' bins and calculates average signal per bin using \code{viewMeans()}, returning
#' matrices suitable for plotting heatmaps (via \code{\link{featureAlignedHeatmap}})
#' or line plots (via \code{\link{featureAlignedDistribution}}).
#' 
#' This is a core function for metagene analysis, allowing visualization of signal
#' patterns across many features simultaneously. The function handles strand
#' information (reverses tiles for negative strand features), filters out-of-bound
#' features, and manages missing seqlevels gracefully.
#' 
#' @param cvglists A list of \code{\link[IRanges:AtomicList-class]{SimpleRleList}}
#'        or \code{\link[IRanges:AtomicList-class]{RleList}} objects containing
#'        coverage/signal data. Each element represents a different sample or
#'        condition. List names will be preserved in the output. Alternatively, a
#'        single \code{SimpleRleList} or \code{RleList} can be provided (will be
#'        automatically converted to a list).
#'        
#'        The RleList objects should have seqlevels matching those in
#'        \code{feature.gr}. If a seqlevel is missing, a zero-filled Rle will be
#'        created for that chromosome (with a warning).
#' @param feature.gr An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#'        with identical width for all ranges. The features from which signal will
#'        be extracted. Each range will be divided into \code{n.tile} bins. Must
#'        contain at least 2 features.
#' @param upstream An integer specifying the number of base pairs upstream from
#'        the center of each feature to include. If provided, \code{downstream} must
#'        also be provided. When both are provided, \code{feature.gr} will be
#'        adjusted: ranges are centered at their midpoint (width set to 1), then
#'        extended using \code{promoters()} by \code{upstream} and
#'        \code{downstream + 1} base pairs.
#' @param downstream An integer specifying the number of base pairs downstream
#'        from the center of each feature to include. Must be provided together with
#'        \code{upstream}. See \code{upstream} for details.
#' @param n.tile An integer specifying the number of tiles/bins to divide each
#'        feature into. Default is 100. Each tile will have approximately equal
#'        width (using \code{tile()} function). Signal values are averaged within
#'        each tile using \code{viewMeans()}.
#' @param ... Additional parameters (currently not used)
#' 
#' @return Returns a list of matrices, one per element in \code{cvglists}. Each
#'        matrix has:
#'        \itemize{
#'          \item Rows: Features (one per range in \code{feature.gr}), with
#'                rownames in the format \code{"seqnames:start-end"}
#'          \item Columns: Tiles/bins (1 to \code{n.tile}), ordered from start
#'                to end of the feature (or end to start for negative strand
#'                features)
#'          \item Values: Average signal intensity within each tile, calculated
#'                using \code{viewMeans()}
#'        }
#'        The list names match the names of \code{cvglists} (or are auto-generated
#'        if a single RleList was provided).
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates inputs and processes \code{feature.gr} (adjusts for
#'         \code{upstream}/\code{downstream} if provided)
#'   \item Filters features with start positions < 1 or end positions beyond
#'         chromosome lengths (with warnings)
#'   \item Tiles each feature into \code{n.tile} bins using \code{tile()}
#'   \item For negative strand features, reverses the tile order to maintain
#'         consistent orientation (5' to 3' relative to the feature)
#'   \item For each sample in \code{cvglists}:
#'         \itemize{
#'           \item Extracts signal using \code{Views()} on the RleList
#'           \item Calculates mean signal per tile using \code{viewMeans()}
#'           \item Organizes results into a matrix (rows = features, columns = tiles)
#'         }
#' }
#' 
#' \strong{Strand handling:}
#' For features on the negative strand, tiles are reversed so that the first
#' column represents the 5' end of the feature and the last column represents
#' the 3' end. This ensures consistent visualization regardless of strand.
#' 
#' \strong{Edge case handling:}
#' \itemize{
#'   \item Features with start < 1: Filtered out (with warning)
#'   \item Features extending beyond chromosome length: Filtered out (with warning)
#'   \item Missing seqlevels in \code{cvglists}: Zero-filled Rle created (with warning)
#'   \item NA values in signal: Preserved in output (may cause issues in downstream
#'         plotting functions)
#'   \item Infinite values: Converted to NA (with warning)
#' }
#' 
#' \strong{Performance considerations:}
#' The function uses efficient \code{Views()} and \code{viewMeans()} operations
#' from the IRanges package. For large numbers of features or high-resolution
#' tiling, computation time scales with the number of features and \code{n.tile}.
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{featureAlignedHeatmap}} for heatmap visualization,
#'          \code{\link{featureAlignedDistribution}} for average signal plots,
#'          \code{\link{featureAlignedExtendSignal}} for extending signal beyond
#'          feature boundaries, \code{\link[IRanges]{Views}} and
#'          \code{\link[IRanges]{viewMeans}} for the underlying signal extraction
#' @keywords misc
#' @export
#' @import GenomicRanges
#' @importFrom BiocGenerics width start end strand `%in%`
#' @importFrom S4Vectors runLength runValue Rle
#' @importFrom GenomeInfoDb seqlengths
#' @examples
#' \dontrun{
#' ## Example 1: Basic usage
#' cvglists <- list(A = RleList(chr1 = Rle(sample.int(5000, 100), 
#'                                         sample.int(300, 100))), 
#'                  B = RleList(chr1 = Rle(sample.int(5000, 100), 
#'                                         sample.int(300, 100))))
#' feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width = 100))
#' signal_matrices <- featureAlignedSignal(cvglists, feature.gr)
#' 
#' ## Example 2: Using upstream and downstream
#' feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width = 100))
#' signal_matrices <- featureAlignedSignal(cvglists, feature.gr,
#'                                         upstream = 2000, downstream = 2000,
#'                                         n.tile = 200)
#' 
#' ## Example 3: Single RleList (auto-converted to list)
#' cvg <- RleList(chr1 = Rle(sample.int(5000, 100), sample.int(300, 100)))
#' signal_matrix <- featureAlignedSignal(cvg, feature.gr)
#' 
#' ## Example 4: Use with featureAlignedHeatmap
#' signal_matrices <- featureAlignedSignal(cvglists, feature.gr)
#' featureAlignedHeatmap(signal_matrices, feature.gr, zeroAt = 0.5)
#' }
#' 
featureAlignedSignal <- function(cvglists, feature.gr, 
                                 upstream, downstream, 
                                 n.tile = 100L, ...) {
    stopifnot(inherits(feature.gr, "GRanges"))
    
    grWidr <- unique(width(feature.gr))
    if (missing(upstream) || missing(downstream)) {
        if (length(grWidr) != 1L) {
            stop("The width of feature.gr is not identical", call. = FALSE)
        }
    } else {
        if (!is.numeric(upstream) || !is.numeric(downstream)) {
            stop("'upstream' and 'downstream' must be numeric", call. = FALSE)
        }
        if (upstream < 0L || downstream < 0L) {
            stop("'upstream' and 'downstream' must be non-negative", 
                 call. = FALSE)
        }
        upstream <- as.integer(upstream)
        downstream <- as.integer(downstream)
        if (length(grWidr) != 1L || any(grWidr != 1L)) {
            start(feature.gr) <- start(feature.gr) + 
                floor(width(feature.gr) / 2L)
            width(feature.gr) <- 1L
            warning("feature.gr is set to the center of feature.gr", 
                   call. = FALSE)
        }
        feature.gr <- promoters(feature.gr, upstream = upstream,
                                downstream = downstream + 1L)
        grWidr <- unique(width(feature.gr))
    }
    if (any(start(feature.gr) < 1L)) {
        warning("Some start positions of the peaks are less than 1. ",
                "They will be filtered.", call. = FALSE)
        feature.gr <- feature.gr[start(feature.gr) > 0L]
    }
    if (length(feature.gr) < 2L) {
        stop("Length of feature.gr must be at least 2", call. = FALSE)
    }
    if (inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))) {
        cvglistsName <- substitute(deparse(cvglists))
        cvglists <- list(cvglists)
        names(cvglists) <- cvglistsName
    }
    if (!is.list(cvglists)) {
        stop("'cvglists' must be a list of SimpleRleList or RleList", 
             call. = FALSE)
    }
    cls <- vapply(cvglists, inherits, 
                  what = c("SimpleRleList", "RleList", "CompressedRleList"),
                  FUN.VALUE = logical(1))
    if (any(!cls)) {
        stop("'cvglists' must be a list of SimpleRleList or RleList", 
             call. = FALSE)
    }
    seqLen <- lapply(cvglists, function(.ele) 
        sapply(.ele, function(.e) sum(runLength(.e))))
    seqLen_keep <- table(unlist(sapply(seqLen, names))) == length(cvglists)
    seqLen <- seqLen[[1L]][seqLen_keep]
    seqLen <- seqLen[!is.na(seqLen)]
    if (length(seqLen) > 0L) {
        feature_gr_subset <- 
            feature.gr[seqnames(feature.gr) %in% names(seqLen)]
        if (any(end(feature_gr_subset) >
                seqLen[as.character(seqnames(feature_gr_subset))])) {
            warning("Some end positions of the peaks are out of bound. ",
                    "They will be filtered.", call. = FALSE)
            feature.gr <- 
                feature_gr_subset[end(feature_gr_subset) <=
                                  seqLen[as.character(seqnames(feature_gr_subset))]]
        }
    }
    #feature.gr.bck <- feature.gr
    stopifnot(is.numeric(n.tile))
    n.tile <- round(n.tile)
    grL <- tile(feature.gr, n = n.tile)
    idx <- as.character(strand(feature.gr)) == "-"
    if (sum(idx) > 0) {
        grL.rev <- grL[idx]
        grL.rev.len <- lengths(grL.rev)
        grL.rev <- unlist(grL.rev, use.names = FALSE)
        grL.rev$oid <- rep(seq.int(length(grL[idx])), grL.rev.len)
        grL.rev <- rev(grL.rev)
        grL.rev.oid <- grL.rev$oid
        grL.rev$oid <- NULL
        grL.rev <- split(grL.rev, grL.rev.oid)
        grL[idx] <- as(grL.rev, "CompressedGRangesList")
        rm(grL.rev, grL.rev.len, grL.rev.oid)
    }
    grL.len <- lengths(grL)
    grL <- unlist(grL)
    grL$oid <- rep(seq_along(feature.gr), grL.len)
    grL$nid <- unlist(lapply(grL.len, seq.int))
    grL.len <- length(grL)
    grL.s <- split(grL, as.character(seqnames(grL)))
    grL <- unlist(grL.s) ## make sure grL and grL.s keep the same order
    seqn <- names(grL.s)
    seql <- seqlengths(feature.gr)
    seql <- seql[seqn]
    seql.f <- range(feature.gr)
    seql.f <- seql.f[match(seqn, seqnames(seql.f))]
    seql[is.na(seql)] <- width(seql.f)[is.na(seql)] + 1
    trimChar <- function(x, width) {
        len <- nchar(x)
        if (len <= width) return(x)
        return(paste0(strtrim(x, width = width), "..."))
    }
    rowname.feature.gr <- paste0(as.character(seqnames(feature.gr)), ":",
                                 start(feature.gr), "-",
                                 end(feature.gr))
    cov <- lapply(cvglists, function(.dat) {
        .dat <- .dat[seqn[seqn %in% names(.dat)]]
        if (length(.dat) != length(seqn)) {
            warning(paste(seqn[!seqn %in% names(.dat)], collapse = ", "), 
                    ifelse(length(seqn[!seqn %in% names(.dat)]) > 1, " are", " is"), 
                    " not in cvglists. seqlevels of cvglist are ", 
                    trimChar(paste(names(.dat), collapse = ", "), width = 60))
            for (i in seqn[!seqn %in% names(.dat)]) {
                .dat[[i]] <- Rle(0, seql[i])
            }
        }
        warn <- sapply(.dat, function(.ele) {
            any(is.na(runValue(.ele)))
        })
        if (any(warn)) {
            warning("cvglists contain NA values.")
        }
        .dat <- sapply(.dat, function(.ele) {
            if (any(is.infinite(runValue(.ele)))) {
                warning("cvglists contain infinite values. ", 
                        "infinite value will be converted to NA.")
                runValue(.ele)[is.infinite(runValue(.ele))] <- NA
            }
            .ele
        })
        .dat <- .dat[seqn]
        vw <- Views(as(.dat, "RleList"), grL.s)
        vm <- viewMeans(vw)
        vm <- unlist(vm)
        stopifnot(length(vm) == grL.len)
        mm <- matrix(0, nrow = length(feature.gr), ncol = n.tile)
        mm[grL$oid + length(feature.gr) * (grL$nid - 1)] <- vm
        rownames(mm) <- rowname.feature.gr
        mm
    })
    
    cov
}
