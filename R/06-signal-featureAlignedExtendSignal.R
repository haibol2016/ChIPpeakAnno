#' Extract and normalize signals from BAM files aligned to genomic features
#' 
#' @description 
#' Extracts and aggregates read coverage from BAM files across specified
#' genomic features. Reads are extended to estimated fragment length before
#' counting, which is important for ChIP-seq and ATAC-seq data where the
#' actual binding site is between the two ends of a sequenced fragment. The
#' function tiles each feature region (defined by upstream and downstream
#' distances) into bins and calculates normalized signal per bin, returning
#' matrices suitable for plotting heatmaps or line plots.
#' 
#' This function is designed for DNA-seq data (ChIP-seq, ATAC-seq, DNase-seq,
#' etc.) and requires features to be single base pair positions (width=1),
#' typically representing TSS, peak centers, or other point features.
#' 
#' @param bamfiles A character vector of BAM file paths to be processed. If
#'        provided, BAM files will be read. Alternatively, use \code{gal} to
#'        provide pre-loaded alignment data.
#' @param index A character vector of BAM index file names (without the '.bai'
#'        extension). Must have the same length as \code{bamfiles}. Defaults to
#'        \code{bamfiles} (assumes index files have the same name as BAM files
#'        with '.bai' extension).
#' @param feature.gr An object of \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#'        where all ranges must have width=1 (single base pair positions). These
#'        represent the center points around which signal will be extracted
#'        (e.g., TSS, peak centers, binding sites). The function will extract
#'        signal from \code{upstream} bp upstream to \code{downstream} bp
#'        downstream of each feature.
#' @param upstream An integer specifying the number of base pairs upstream from
#'        each feature to include in signal extraction. Required parameter.
#' @param downstream An integer specifying the number of base pairs downstream
#'        from each feature to include in signal extraction. Required parameter.
#' @param n.tile An integer specifying the number of tiles/bins to divide the
#'        region (upstream + downstream) into. Default is 100. Each tile will
#'        represent (upstream + downstream) / n.tile base pairs.
#' @param fragmentLength An integer or numeric vector specifying the estimated
#'        fragment length for each BAM file. Required parameter. For single-end
#'        data, reads are extended to this length. For paired-end data, this
#'        should match the average insert size. If a vector, must have length
#'        equal to the number of BAM files or GAlignments objects.
#' @param librarySize An integer or numeric vector specifying the estimated
#'        library size (total number of mapped reads) for each BAM file.
#'        Required parameter. Used for normalization. If a vector, must have
#'        length equal to the number of BAM files or GAlignments objects.
#' @param pe A character string specifying whether the data is paired-end.
#'        Options: \code{"auto"} (default, automatically detect from BAM files
#'        or GAlignments), \code{"PE"} (paired-end), or \code{"SE"} (single-end).
#'        When \code{"auto"} and using \code{gal}, detection is based on
#'        duplicate query names in GAlignments objects.
#' @param adjustFragmentLength An optional numeric value (length 1) to adjust
#'        all fragments/reads to a fixed length after reading. If provided,
#'        all reads will be re-centered and extended/trimmed to this length,
#'        overriding the \code{fragmentLength} parameter for the final signal
#'        calculation.
#' @param gal A \code{GAlignmentsList} object or a list of \code{GAlignments}
#'        or \code{GAlignmentPairs} objects. If provided, BAM files are not
#'        read and this pre-loaded alignment data is used instead. Required if
#'        \code{bamfiles} is missing.
#' @param \dots Additional parameters (currently not used)
#' 
#' @return Returns a list of matrices, one per BAM file or GAlignments object.
#'        Each matrix has:
#'        \itemize{
#'          \item Rows: Features (one row per element in \code{feature.gr})
#'          \item Columns: Tiles/bins (1 to \code{n.tile})
#'          \item Values: Normalized signal density (RPKM-like normalization:
#'                reads per 100 million reads per fragment length per base pair)
#'        }
#'        The normalization formula is: \code{count * 1e8 / librarySize * 100 /
#'        fragmentLength / totalBPinBin}, where \code{totalBPinBin} is the number
#'        of base pairs per bin.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates that all ranges in \code{feature.gr} have width=1
#'   \item Expands feature regions by \code{upstream} and \code{downstream}
#'   \item Tiles each expanded region into \code{n.tile} bins
#'   \item For negative strand features, reverses the tile order to maintain
#'         consistent orientation (5' to 3')
#'   \item Reads BAM files or uses provided GAlignments objects
#'   \item Extends single-end reads to \code{fragmentLength} (for negative strand,
#'         extends from the 3' end)
#'   \item For paired-end data, uses the fragment span directly
#'   \item Optionally adjusts all fragments to \code{adjustFragmentLength}
#'   \item Counts overlaps between extended reads and each tile
#'   \item Normalizes counts by library size, fragment length, and bin size
#' }
#' 
#' \strong{Read extension:}
#' For single-end reads, the function extends reads to \code{fragmentLength}:
#' \itemize{
#'   \item Positive strand: extends from the 5' end (start position)
#'   \item Negative strand: extends from the 3' end (end position), adjusting
#'         the start position accordingly
#' }
#' For paired-end reads, the fragment span is used directly (no extension needed).
#' 
#' \strong{Normalization:}
#' The signal is normalized using a RPKM-like approach:
#' \deqn{signal = \frac{count \times 10^8}{librarySize} \times \frac{100}{fragmentLength} \times \frac{1}{totalBPinBin}}
#' This accounts for library size differences, fragment length, and bin size.
#' 
#' \strong{Strand handling:}
#' The function ignores strand information when counting overlaps
#' (\code{ignore.strand=TRUE}), which is appropriate for most ChIP-seq and
#' ATAC-seq data. However, tiles are ordered consistently (5' to 3') for both
#' positive and negative strand features.
#' 
#' \strong{Filtering:}
#' When reading BAM files, the function filters out:
#' \itemize{
#'   \item Secondary alignments
#'   \item Reads failing quality controls
#'   \item For paired-end: only properly paired reads are used
#' }
#' @author Jianhong Ou
#' @seealso \code{\link{featureAlignedSignal}} for signal extraction from
#'          RleList objects, \code{\link{estLibSize}} for estimating library
#'          sizes, \code{\link{estFragmentLength}} for estimating fragment
#'          lengths, \code{\link{featureAlignedHeatmap}} for visualizing the
#'          results, \code{\link{reCenterPeaks}} for converting peaks to
#'          single-base positions
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors elementNROWS mcols
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs granges
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamWhat scanBamFlag
#' @importFrom BiocGenerics start end width strand
#' @examples
#' \dontrun{
#' ## Example 1: Extract signal from BAM files around peak centers
#' path <- system.file("extdata", package="MMDiffBamSubset")
#' if(file.exists(path)){
#'     WT.AB2 <- file.path(path, "reads", "WT_2.bam")
#'     Null.AB2 <- file.path(path, "reads", "Null_2.bam")
#'     Resc.AB2 <- file.path(path, "reads", "Resc_2.bam")
#'     peaks <- file.path(path, "peaks", "WT_2_Macs_peaks.xls")
#'     
#'     ## Estimate library sizes
#'     libSizes <- estLibSize(c(WT.AB2, Null.AB2, Resc.AB2))
#'     
#'     ## Load peaks and convert to single-base positions (peak centers)
#'     feature.gr <- toGRanges(peaks, format="MACS")
#'     feature.gr <- feature.gr[seqnames(feature.gr)=="chr1" & 
#'                              start(feature.gr)>3000000 & 
#'                              end(feature.gr)<75000000]
#'     feature.gr <- reCenterPeaks(feature.gr, width=1)
#'     
#'     ## Extract signal 505bp upstream and downstream of peak centers
#'     sig <- featureAlignedExtendSignal(c(WT.AB2, Null.AB2, Resc.AB2), 
#'                                       feature.gr = feature.gr, 
#'                                       upstream = 505,
#'                                       downstream = 505,
#'                                       n.tile = 101, 
#'                                       fragmentLength = 250,
#'                                       librarySize = libSizes)
#'     
#'     ## Visualize as heatmap
#'     featureAlignedHeatmap(sig, 
#'                          reCenterPeaks(feature.gr, width=1010), 
#'                          zeroAt = 0.5, n.tile = 101)
#' }
#' 
#' ## Example 2: Using GAlignments objects (pre-loaded alignments)
#' ## gal <- readGAlignments("sample.bam")
#' ## sig <- featureAlignedExtendSignal(gal = list(gal), 
#' ##                                  feature.gr = tss_gr,
#' ##                                  upstream = 2000, downstream = 2000,
#' ##                                  n.tile = 100,
#' ##                                  fragmentLength = 200,
#' ##                                  librarySize = 1e6)
#' }
#' 
featureAlignedExtendSignal <- function(bamfiles, index = bamfiles, 
                                       feature.gr, 
                                       upstream, downstream, 
                                       n.tile = 100, 
                                       fragmentLength,
                                       librarySize,
                                       pe = c("auto", "PE", "SE"),
                                       adjustFragmentLength,
                                       gal, ...) {
    #message("The signal is being calculated for DNA-seq.")
    if (missing(fragmentLength)) {
        stop("'fragmentLength' is required", call. = FALSE)
    }
    if (!missing(adjustFragmentLength)) {
        stopifnot(inherits(adjustFragmentLength, c("numeric", "integer")))
        stopifnot(length(adjustFragmentLength) == 1)
    }
    if (!missing(bamfiles)) {
        stopifnot(is.character(bamfiles))
        stopifnot(length(bamfiles) == length(index))
        galInput <- FALSE
    } else {
        if (missing(gal)) {
            stop("gal is required if missing bamfiles")
        } else {
            if (!inherits(gal, "GAlignmentsList")) {
                galInput <- sapply(gal, function(.ele) {
                    inherits(.ele, c("GAlignments", "GAlignmentPairs"))
                })
                if (any(!galInput)) {
                    stop("'gal' must be a GAlignmentsList object or ",
                         "a list of GAlignmentPairs", call. = FALSE)
                }
            }
        }
        galInput <- TRUE
    }
    stopifnot(inherits(feature.gr, "GRanges"))
    if (missing(upstream) || missing(downstream)) {
        stop("'upstream' and 'downstream' are required", call. = FALSE)
    }
    stopifnot(inherits(upstream, c("numeric", "integer")))
    stopifnot(inherits(downstream, c("numeric", "integer")))
    stopifnot(inherits(n.tile, c("numeric", "integer")))
    stopifnot(inherits(fragmentLength, c("numeric", "integer")))
    stopifnot(inherits(librarySize, c("numeric", "integer")))
    pe <- match.arg(pe)
    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    n.tile <- as.integer(n.tile)
    fragmentLength <- as.integer(fragmentLength)
    librarySize <- as.integer(librarySize)
    stopifnot(all(width(feature.gr) == 1))
    
    totalBPinBin <- floor((upstream + downstream) / n.tile) # * length(feature.gr)
    feature.gr$oid <- seq_along(feature.gr)
    feature.gr.expand <- 
        suppressWarnings(promoters(feature.gr, 
                                   upstream = upstream + max(fragmentLength), 
                                   downstream = downstream + max(fragmentLength)))
    strand(feature.gr.expand) <- "*"
    feature.gr <- suppressWarnings(promoters(feature.gr,
                                             upstream = upstream, 
                                             downstream = downstream))
    
    grL <- tile(feature.gr, n = n.tile)
    ## reorder the tiles for negative strand and set group id for them.
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
    grL.eleLen <- elementNROWS(grL)
    grs <- unlist(grL)
    grs$gpid <- sequence(grL.eleLen) 
    grs$oid <- rep(feature.gr$oid, grL.eleLen)
    if (galInput) {
        bams.gr <- mapply(function(ga, .fLen) {
            if (inherits(ga, "GAlignmentPairs")) {
                granges(ga)
            } else {
                qname <- mcols(ga)$qname
                if (pe == "auto") {
                    ## check qname
                    if (length(qname) == 0) {
                        pe <- FALSE
                    } else {
                        ## check duplicated qname
                        if (all(table(qname) == 1)) {
                            pe <- FALSE
                        } else {
                            if (all(table(qname) < 3)) {
                                pe <- TRUE
                            } else {
                                pe <- FALSE
                            }
                        }
                    }
                } else {
                    pe <- pe == "PE"
                }
                
                if (pe) {
                    ga_split <- split(ga, qname)
                    granges(ga_split, ignore.strand = TRUE)
                } else {
                    ga_gr <- granges(ga)
                    start(ga_gr[strand(ga_gr) == "-"]) <- 
                        end(ga_gr[strand(ga_gr) == "-"]) - .fLen + 1
                    width(ga_gr) <- .fLen
                    ga_gr
                }
            }
        }, gal, fragmentLength, SIMPLIFY = FALSE)
    } else {
        if (pe == "auto") {
            pe <- mapply(function(file, id) 
                suppressMessages(testPairedEndBam(file, index = id)), 
                bamfiles, index)
        } else {
            pe <- pe == "PE"
        }
        param <- ScanBamParam(which = reduce(feature.gr.expand), 
                              flag = scanBamFlag(isSecondaryAlignment = FALSE,
                                                 isNotPassingQualityControls = FALSE),
                              what = scanBamWhat())
        paramp <- ScanBamParam(which = reduce(feature.gr.expand), 
                               flag = scanBamFlag(
                                   isProperPair = TRUE,
                                   isSecondaryAlignment = FALSE,
                                   isNotPassingQualityControls = FALSE),
                               what = scanBamWhat())
        bams.gr <- mapply(function(f, i, p, .fLen) {
            if (!p) {
                ga <- readGAlignments(f, index = i, param = param)
                ga_gr <- granges(ga)
                start(ga_gr[strand(ga_gr) == "-"]) <- 
                    end(ga_gr[strand(ga_gr) == "-"]) - .fLen + 1
                width(ga_gr) <- .fLen
                ga_gr
            } else {
                gap <- readGAlignmentPairs(f, index = i, param = paramp)
                granges(gap)
            }
        }, bamfiles, index, pe, fragmentLength, SIMPLIFY = FALSE)
    }
    
    if (!missing(adjustFragmentLength)) {
        bams.gr <- lapply(bams.gr, reCenterPeaks, 
                          width = adjustFragmentLength)
        fragmentLength <- adjustFragmentLength
    }
    ## count overlaps
    co <- lapply(bams.gr, countOverlaps, 
                 query = grs, ignore.strand = TRUE)
    
    countTable <- do.call(cbind, co)
    colnames(countTable) <- names(bams.gr)
    stopifnot(nrow(countTable) == length(grs))
    #    sumByGpid <- rowsum(countTable, grs$gpid, reorder = FALSE)
    #    signal <- t(t(sumByGpid)*1e8/librarySize*100/fragmentLength)/totalBPinBin
    countTable.list <- as.list(as.data.frame(countTable))
    signal <- mapply(function(.ele, .libsize, .fLen) {
        tbl <- matrix(nrow = n.tile, ncol = length(feature.gr))
        tbl[(grs$oid - 1) * n.tile + grs$gpid] <- .ele
        t(tbl) * 1e8 / .libsize * 100 / .fLen / totalBPinBin
    }, countTable.list, librarySize, fragmentLength, SIMPLIFY = FALSE)
    return(signal)
}
