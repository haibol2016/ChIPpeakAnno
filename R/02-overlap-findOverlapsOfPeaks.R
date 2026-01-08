#' Find overlapping peaks among two or more peak sets
#' 
#' @description 
#' Identifies overlapping genomic regions between two to five sets of peaks using
#' an efficient interval tree algorithm implemented in GenomicRanges. This function
#' is the recommended approach for finding overlaps between multiple peak sets,
#' replacing the deprecated \code{\link{findOverlappingPeaks}} function.
#' 
#' The function efficiently performs overlap queries and provides detailed
#' information about spatial relationships between overlapping peaks, including
#' distance calculations and overlap classifications. It returns an object of class
#' \code{overlappingPeaks} containing Venn counts, merged peaks, unique peaks, and
#' detailed overlap annotations.
#' 
#' @aliases findOverlapsOfPeaks overlappingPeaks overlappingPeaks-class
#' @param \dots Two to five objects of \link[GenomicRanges:GRanges-class]{GRanges}
#'        containing peak sets to compare. Alternatively, a single list of GRanges
#'        objects can be provided. Peak sets must have unique names (auto-generated
#'        if missing). See examples below.
#' @param maxgap An integer specifying the maximum gap (in base pairs) allowed
#'        between ranges for them to be considered overlapping. Default is \code{-1L},
#'        which means ranges must actually overlap (no gap allowed). See
#'        \code{\link[IRanges:findOverlaps-methods]{findOverlaps}} for details.
#' @param minoverlap An integer or numeric value specifying the minimum overlap
#'        required. If an integer >= 1, it specifies the minimum number of base pairs
#'        that must overlap. If \code{0 < minoverlap < 1}, it specifies the minimum
#'        percentage of the interval that must be covered, and the filter condition
#'        will be set to the maximum covered percentage of overlapping peaks.
#'        Default is \code{0L} (any overlap is considered).
#' @param ignore.strand A logical value. When \code{TRUE} (default), strand
#'        information is ignored in overlap calculations. When \code{FALSE}, only
#'        ranges on the same strand are considered for overlap.
#' @param connectedPeaks A character string specifying how to count connected/overlapping
#'        peak groups. Options:
#'        \itemize{
#'          \item \code{"keepAll"} (default): Adds the number of involved peaks for
#'                each peak list to the corresponding overlapping counts. Also outputs
#'                counts as if \code{connectedPeaks} were set to "min". For example,
#'                if 5 peaks in group1 overlap with 2 peaks in group2, this will add
#'                5 peaks to count.group1, 2 to count.group2, and 2 to counts.
#'          \item \code{"min"}: Adds the minimal number of involved peaks in each group
#'                of connected/overlapped peaks to the overlapping counts. In the
#'                example above, this would add 2 to the overlapping counts.
#'          \item \code{"merge"}: Counts each group of connected peaks as only 1,
#'                regardless of how many peaks are involved. In the example above, this
#'                would add 1 to the overlapping counts.
#'        }
#'        See \url{https://support.bioconductor.org/p/133486/#133603} for more examples.
#' @return Returns an object of class \code{overlappingPeaks} containing:
#'        \itemize{
#'          \item \code{venn_cnt}: An object of class \code{VennCounts} containing
#'                the overlap counts for Venn diagram generation
#'          \item \code{peaklist}: A named list of \code{GRanges} objects, where each
#'                element represents peaks in a specific overlap category (e.g., peaks
#'                unique to one set, peaks shared by two sets, etc.). Names are
#'                constructed from the peak set names joined by "///"
#'          \item \code{uniquePeaks}: A \code{GRanges} object containing all peaks
#'                that are unique to a single peak set (not overlapping with any
#'                other set)
#'          \item \code{mergedPeaks}: A \code{GRanges} object containing all merged
#'                overlapping peaks (peaks that overlap between two or more sets)
#'          \item \code{peaksInMergedPeaks}: A \code{GRanges} object containing all
#'                original peaks from each sample that are involved in overlapping
#'                peaks
#'          \item \code{overlappingPeaks}: A named list of data frames, where each
#'                data frame contains detailed annotation of overlapping peaks between
#'                two specific peak sets. Each data frame includes:
#'                \itemize{
#'                  \item \code{peaks1}, \code{peaks2}: Peak names from each set
#'                  \item Genomic coordinates (seqnames, start, end, strand, width)
#'                  \item \code{overlapFeature}: Spatial relationship between peaks
#'                        ("upstream", "downstream", "inside", "overlapStart",
#'                        "overlapEnd", "includeFeature", "overlap")
#'                  \item \code{shortestDistance}: Shortest distance between the
#'                        overlapping peaks
#'                }
#'          \item \code{all.peaks}: A list of \code{GRanges} objects containing the
#'                input peaks with formatted names (peak set name + "__" + original
#'                peak name)
#'        }
#' @author Jianhong Ou
#' @seealso \link{annotatePeakInBatch}, \link{makeVennDiagram},
#' \link{getVennCounts}, \link{findOverlappingPeaks}
#' @references 1.Interval tree algorithm from: Cormen, Thomas H.; Leiserson,
#' Charles E.; Rivest, Ronald L.; Stein, Clifford. Introduction to Algorithms,
#' second edition, MIT Press and McGraw-Hill. ISBN 0-262-53196-8
#' 
#' 2.Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate
#' ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237doi:10.1186/1471-2105-11-237
#' 
#' 3. Zhu L (2013). "Integrative analysis of ChIP-chip and ChIP-seq dataset."
#' In Lee T and Luk ACS (eds.), Tilling Arrays, volume 1067, chapter 4, pp.
#' -19. Humana Press. http://dx.doi.org/10.1007/978-1-62703-607-8_8,
#' http://link.springer.com/protocol/10.1007\%2F978-1-62703-607-8_8
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @examples
#' 
#' peaks1 <- GRanges(seqnames=c(6,6,6,6,5),
#'                  IRanges(start=c(1543200,1557200,1563000,1569800,167889600),
#'                          end=c(1555199,1560599,1565199,1573799,167893599),
#'                          names=c("p1","p2","p3","p4","p5")),
#'                  strand="+")
#' peaks2 <- GRanges(seqnames=c(6,6,6,6,5),
#'                   IRanges(start=c(1549800,1554400,1565000,1569400,167888600),
#'                           end=c(1550599,1560799,1565399,1571199,167888999),
#'                           names=c("f1","f2","f3","f4","f5")),
#'                   strand="+")
#' t1 <- findOverlapsOfPeaks(peaks1, peaks2, maxgap=1000)
#' makeVennDiagram(t1)
#' t1$venn_cnt
#' t1$peaklist
#' t2 <- findOverlapsOfPeaks(peaks1, peaks2, minoverlap = .5)
#' makeVennDiagram(t2)
#' 
#' t3 <- findOverlapsOfPeaks(peaks1, peaks2, minoverlap = .90)
#' makeVennDiagram(t3)
#' 
findOverlapsOfPeaks <- function(..., maxgap = -1L, minoverlap = 0L,
                                ignore.strand = TRUE, 
                                connectedPeaks = c("keepAll", "min", "merge")) {
    ### Check inputs
    NAME_conn_string <- "___conn___"
    NAME_short_string <- "__"
    NAME_long_string <- "///"
    PeaksList <- list(...)
    n <- length(PeaksList)
    if (n == 1) {
        PeaksList <- PeaksList[[1]]
        n <- length(PeaksList)
        names <- names(PeaksList)
        if (is.null(names)) {
            names <- paste("peaks", seq_len(n), sep = "")
        }
    } else {
        ## Save dots arguments names
        dots <- substitute(list(...))[-1]
        names <- make.names(unlist(sapply(dots, deparse)))
        names(PeaksList) <- names
    }
    if (any(grepl(NAME_short_string, names))) {
        stop("The name of peaks could not contain '", NAME_short_string, "'", 
             call. = FALSE)
    }
    if (n < 2L) {
        stop("At least 2 peak sets are required", call. = FALSE)
    }
    if (n > 5L) {
        stop("Maximum 5 peak sets are supported", call. = FALSE)
    }
    connectedPeaks <- match.arg(connectedPeaks)
    if (any(duplicated(names))) {
        stop("Duplicate peak set names detected", call. = FALSE)
    }
    PeaksList <- lapply(PeaksList, trimPeakList, by = "region", 
                        ignore.strand = ignore.strand, keepMetadata = TRUE)
    venn_cnt <- vennCounts(PeaksList, n = n, names = names, 
                           maxgap = maxgap, minoverlap = minoverlap, by = "region",
                           ignore.strand = ignore.strand, 
                           connectedPeaks = connectedPeaks)
    
    outcomes <- venn_cnt$venn_cnt[, 1:n]
    xlist <- do.call(rbind, venn_cnt$xlist)
    xlist <- xlist - 1
    xlist <- xlist[rev(seq_len(nrow(xlist))), , drop = FALSE] 
    ## Reverse xlist to match the order of names
    xlist <- apply(xlist, 2, base::paste, collapse = "")
    if (length(venn_cnt$all) != length(xlist)) {
        stop("Length of 'xlist' and 'all' should be identical", call. = FALSE)
    }
    all <- cbind(.id = rep(seq_along(venn_cnt$all), sapply(venn_cnt$all, length)), 
                 .ele = unlist(venn_cnt$all), 
                 .gp = rep(xlist, sapply(venn_cnt$all, length)))
    all.peaks <- venn_cnt$Peaks[all[, 2]]
    names(all.peaks) <- gsub(NAME_conn_string,
                             NAME_short_string,
                             names(all.peaks))
    if (!is.null(all.peaks$old_strand_HH)) {
        strand(all.peaks) <- all.peaks$old_strand_HH
        all.peaks$old_strand_HH <- NULL
    }
    all.peaks$gpForFindOverlapsOfPeaks <- all[, 1]
    all.peaks$gpType <- all[, 3]
    all.peaks.split <- split(all.peaks, all.peaks$gpType)
    listname <- apply(outcomes, 1, 
                      function(id) 
                          paste(names[as.logical(id)], 
                                collapse = NAME_long_string))
    listcode <- apply(outcomes, 1, base::paste, collapse = "")
    listname <- listname[-1]
    listcode <- listcode[-1]
    names(listname) <- listcode
    names(all.peaks.split) <- listname[names(all.peaks.split)]
    peaklist <- lapply(all.peaks.split, reduce, 
                       min.gapwidth = maxgap + 1L, with.revmap = TRUE, 
                       ignore.strand = ignore.strand)
    peaklist <- mapply(function(peaks, info) {
        revmap <- peaks$revmap
        peakNames <- cbind(id = unlist(revmap),
                           gp = rep(seq_along(revmap), sapply(revmap, length)))
        peaks$peakNames <- CharacterList(split(names(info)[peakNames[, 1]], 
                                               peakNames[, 2]), 
                                         compress = TRUE)
        if (ignore.strand) {
            strand <- split(as.character(strand(info))[peakNames[, 1]], 
                            peakNames[, 2])
            strand <- lapply(strand, unique)
            l <- sapply(strand, length)
            strand[l >= 2] <- "*"
            strand(peaks) <- unlist(strand)
        }
        peaks$revmap <- NULL
        peaks
    }, peaklist, all.peaks.split, SIMPLIFY = FALSE)
    
    listcode <- strsplit(listcode, "")
    names(listcode) <- listname
    listcode <- sapply(listcode, function(.ele) sum(as.numeric(.ele)) == 2)
    overlappingPeaks <- sapply(names(listcode)[listcode], function(.ele) {
        peakListName <- strsplit(.ele, NAME_long_string, fixed = TRUE)[[1]]
        ps <- PeaksList[peakListName]
        ol <- findOverlaps(query = ps[[1]], subject = ps[[2]], 
                          maxgap = maxgap, minoverlap = minoverlap, 
                          ignore.strand = ignore.strand)
        q <- ps[[1]][queryHits(ol)]
        s <- ps[[2]][subjectHits(ol)]
        cl <- getRelationship(q, s)[, c("insideFeature", "shortestDistance")]
        colnames(cl)[grepl("insideFeature", colnames(cl))] <- "overlapFeature"
        correlation <- cbind(peaks1 = names(q), as.data.frame(unname(q)), 
                            peaks2 = names(s), as.data.frame(unname(s)), 
                            cl)
        rownames(correlation) <- make.names(paste(names(q), names(s), sep = "_"))
        correlation <- correlation[correlation[, "shortestDistance"] < maxgap | 
                                   correlation[, "overlapFeature"] %in% 
                                   c("includeFeature", "inside",
                                     "overlapEnd", "overlapStart",
                                     "overlap"), ]
    }, simplify = FALSE)
    PeaksList <- sapply(PeaksList, trimPeakList, by = "region",
                        ignore.strand = ignore.strand,
                        keepMetadata = TRUE)
    for (i in seq_len(n)) {
        names(PeaksList[[i]]) <- 
            paste(names[i], names(PeaksList[[i]]), sep = NAME_short_string)
    }
    sharedColnames <- Reduce(function(a, b) {
        shared <- intersect(colnames(mcols(a)), colnames(mcols(b)))
        if (length(shared) > 0) {
            shared <- shared[sapply(shared, function(.ele) {
                class(mcols(a)[, .ele]) == class(mcols(b)[, .ele])
            })]
        }
        if (length(shared) > 0) {
            mcols(a) <- mcols(a)[, shared, drop = FALSE]
        } else {
            mcols(a) <- NULL
        }
        a
    }, PeaksList)
    sharedColnames <- colnames(mcols(sharedColnames))
    uniquePeaks <- peaklist[names(peaklist) %in% 
                            listname[!grepl(NAME_long_string, 
                                            listname, fixed = TRUE)]]
    uniquePeaks <- mapply(function(a, b) {
        b <- b[unlist(a$peakNames)]
        if (length(sharedColnames) > 0) {
            mcols(b) <- mcols(b)[, sharedColnames, drop = FALSE]
        } else {
            mcols(b) <- NULL
        }
        b
    }, uniquePeaks, PeaksList[names(uniquePeaks)])
    if (!inherits(uniquePeaks, "GRangesList")) {
        uniquePeaks <- GRangesList(uniquePeaks)
    }
    uniquePeaks <- unlist(uniquePeaks, use.names = FALSE)
    mergedPeaks <- peaklist[names(peaklist) %in% 
                            listname[grepl(NAME_long_string,
                                           listname, fixed = TRUE)]]
    if (!inherits(mergedPeaks, "GRangesList")) {
        mergedPeaks <- GRangesList(mergedPeaks)
    }
    mergedPeaks <- unlist(mergedPeaks, use.names = FALSE)
    peaksInMergedPeaks <- lapply(PeaksList, function(.ele) {
        if (length(sharedColnames) > 0) {
            mcols(.ele) <- mcols(.ele)[, sharedColnames, drop = FALSE]
        } else {
            mcols(.ele) <- NULL
        }
        .ele
    })
    peaksInMergedPeaks <- unlist(GRangesList(peaksInMergedPeaks), 
                                 use.names = FALSE)
    peaksInMergedPeaks <- peaksInMergedPeaks[unlist(mergedPeaks$peakNames)]
    structure(list(venn_cnt = venn_cnt$venn_cnt, 
                   peaklist = peaklist, 
                   uniquePeaks = uniquePeaks,
                   mergedPeaks = mergedPeaks,
                   peaksInMergedPeaks = peaksInMergedPeaks,
                   overlappingPeaks = overlappingPeaks,
                   all.peaks = PeaksList), 
              class = "overlappingPeaks")
}
