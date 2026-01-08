#' Find overlapping peaks between two peak sets (deprecated)
#' 
#' \strong{This function is deprecated.} Use \code{\link{findOverlapsOfPeaks}}
#' instead, which supports multiple peak sets and provides better performance.
#' 
#' This function identifies overlapping genomic regions between two to five sets
#' of peaks using an interval tree algorithm implemented in IRanges and graph
#' algorithms to find connected components. It returns detailed information about
#' the spatial relationships between overlapping peaks, including Venn counts and
#' merged peak regions.
#' 
#' @aliases findOverlappingPeaks findOverlappingPeaks-deprecated
#' @param Peaks1 An optional \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing the first set of peaks. If provided along with \code{Peaks2},
#'        the function will return additional fields specific to the two-peak-set
#'        comparison.
#' @param Peaks2 An optional \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing the second set of peaks. If provided along with \code{Peaks1},
#'        the function will return additional fields specific to the two-peak-set
#'        comparison.
#' @param maxgap An integer specifying the maximum gap (in base pairs) allowed
#'        between ranges for them to be considered overlapping. Default is
#'        \code{-1L}, which means ranges must actually overlap (no gap allowed).
#'        See \code{\link[IRanges:findOverlaps-methods]{findOverlaps}} for details.
#' @param minoverlap An integer or numeric value specifying the minimum overlap
#'        required. If an integer >= 1, it specifies the minimum number of base
#'        pairs that must overlap. If \code{0 < minoverlap < 1}, it specifies the
#'        minimum percentage of the interval that must be covered. Default is
#'        \code{0L} (any overlap is considered).
#' @param multiple A logical value (deprecated). This parameter is kept for
#'        backward compatibility but is not used in the current implementation.
#'        Use \code{select} parameter instead (though \code{select} is also not
#'        currently used in the implementation).
#' @param NameOfPeaks1 A character string specifying the name for \code{Peaks1}.
#'        Used for generating column names and identifying peak sets in the output.
#'        Default is \code{"TF1"}. If \code{Peaks1} is missing, defaults to
#'        \code{"Peaks1"}.
#' @param NameOfPeaks2 A character string specifying the name for \code{Peaks2}.
#'        Used for generating column names and identifying peak sets in the output.
#'        Default is \code{"TF2"}. If \code{Peaks2} is missing, defaults to
#'        \code{"Peaks2"}.
#' @param select A character string (deprecated). Options: \code{"all"} (return
#'        multiple overlapping peaks), \code{"first"} (return the first overlapping
#'        peak), \code{"last"} (return the last overlapping peak), or
#'        \code{"arbitrary"} (return one of the overlapping peaks). This parameter
#'        is kept for backward compatibility but is not currently used in the
#'        implementation.
#' @param annotate An integer (deprecated). \code{1} means include
#'        \code{overlapFeature} and \code{shortestDistance} in the output,
#'        \code{0} means do not include them. Default is \code{0}. This parameter
#'        is kept for backward compatibility but is not currently used in the
#'        implementation. The function always includes overlap annotations when
#'        both \code{Peaks1} and \code{Peaks2} are provided.
#' @param ignore.strand A logical value. When \code{TRUE} (default), strand
#'        information is ignored in overlap calculations. When \code{FALSE}, only
#'        ranges on the same strand are considered for overlap.
#' @param connectedPeaks A character string specifying how to count
#'        connected/overlapping peak groups. Options:
#'        \itemize{
#'          \item \code{"min"} (default): Counts the minimal number of involved
#'                peaks in each group of connected/overlapped peaks
#'          \item \code{"merge"}: Counts each group of connected peaks as only 1,
#'                regardless of how many peaks are involved
#'        }
#' @param \dots Additional \link[GenomicRanges:GRanges-class]{GRanges} objects
#'        containing peak sets. The function supports 2-5 peak sets total (including
#'        \code{Peaks1}, \code{Peaks2}, and any additional sets provided via
#'        \code{...}). Peak sets can also be provided as a single list of GRanges
#'        objects. See \code{\link{findOverlapsOfPeaks}} for the recommended
#'        approach.
#' 
#' @return Returns an object of class \code{overlappingPeaks} containing:
#'        \itemize{
#'          \item \code{venn_cnt}: An object of class \code{VennCounts} containing
#'                the overlap counts for Venn diagram generation
#'          \item \code{peaklist}: A named list of \code{GRanges} objects, where
#'                each element represents peaks in a specific overlap category
#'                (e.g., peaks unique to one set, peaks shared by two sets, etc.).
#'                Names are constructed from the peak set names joined by "///"
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
#'        }
#'        
#'        \strong{Additional fields when both \code{Peaks1} and \code{Peaks2} are provided:}
#'        \itemize{
#'          \item \code{OverlappingPeaks}: A data frame containing the overlap
#'                annotations between \code{Peaks1} and \code{Peaks2} (same as
#'                \code{overlappingPeaks[[sampleName]]})
#'          \item \code{MergedPeaks}: A \code{GRanges} object containing merged
#'                overlapping peaks between \code{Peaks1} and \code{Peaks2}
#'          \item \code{Peaks1withOverlaps}: A \code{GRanges} object containing
#'                all peaks from \code{Peaks1} that overlap with \code{Peaks2}
#'          \item \code{Peaks2withOverlaps}: A \code{GRanges} object containing
#'                all peaks from \code{Peaks2} that overlap with \code{Peaks1}
#'        }
#'        
#'        \strong{Note:} The function uses graph algorithms (from the \code{graph}
#'        and \code{RBGL} packages) to identify connected components of overlapping
#'        peaks. Peaks with duplicated or missing names are automatically renamed
#'        using zero-padded numeric indices.
#' @details
#' 
#' \strong{Algorithm Overview:}
#' 
#' The function uses the following approach to identify overlapping peaks:
#' \enumerate{
#'   \item Combines all input peak sets into a single \code{GRangesList}
#'   \item Uses \code{findOverlaps()} to identify all pairwise overlaps
#'   \item Constructs a graph where nodes are peaks and edges represent overlaps
#'   \item Uses graph algorithms to find connected components (groups of peaks
#'         that are connected through overlaps)
#'   \item Calculates Venn counts based on which peak sets are represented in
#'         each connected component
#'   \item Merges overlapping peaks within each component using \code{reduce()}
#'   \item Annotates spatial relationships between overlapping peaks using
#'         \code{getRelationship()}
#' }
#' 
#' \strong{Metadata Handling:}
#' 
#' The function automatically identifies shared metadata columns across all peak
#' sets and keeps only those columns where the data types are compatible (same
#' class). This ensures that merged peaks can retain meaningful metadata.
#' 
#' \strong{Deprecated Parameters:}
#' 
#' The parameters \code{multiple}, \code{select}, and \code{annotate} are
#' defined in the function signature for backward compatibility but are not
#' currently used in the implementation. The function always returns all
#' overlapping peaks with full annotations when both \code{Peaks1} and
#' \code{Peaks2} are provided.
#' 
#' @author Lihua Julie Zhu
#' @seealso \code{\link{findOverlapsOfPeaks}} (recommended replacement),
#'          \code{\link{annotatePeakInBatch}}, \code{\link{makeVennDiagram}}
#' @references 1.Interval tree algorithm from: Cormen, Thomas H.; Leiserson,
#' Charles E.; Rivest, Ronald L.; Stein, Clifford. Introduction to Algorithms,
#' second edition, MIT Press and McGraw-Hill. ISBN 0-262-53196-8
#' 
#' 2.Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate
#' ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 11:237
#' doi:10.1186/1471-2105-11-237
#' 
#' 3. Zhu L (2013). Integrative analysis of ChIP-chip and ChIP-seq dataset.  In
#' Lee T and Luk ACS (eds.), Tilling Arrays, volume 1067, chapter 4, pp. -19.
#' Humana Press. http://dx.doi.org/10.1007/978-1-62703-607-8_8
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors mcols DataFrame
#' @importClassesFrom graph graphNEL
#' @importFrom graph ugraph
#' @importFrom RBGL connectedComp
#' @importFrom utils data
#' @examples
#' 
#'     if (interactive())
#'     {    
#'     peaks1 = 
#'         GRanges(seqnames=c(6,6,6,6,5), 
#'                 IRanges(start=c(1543200,1557200,1563000,1569800,167889600),
#'                         end=c(1555199,1560599,1565199,1573799,167893599),
#'                         names=c("p1","p2","p3","p4","p5")),
#'                 strand=as.integer(1))
#'     peaks2 = 
#'         GRanges(seqnames=c(6,6,6,6,5), 
#'                 IRanges(start=c(1549800,1554400,1565000,1569400,167888600),
#'                         end=c(1550599,1560799,1565399,1571199,167888999),
#'                         names=c("f1","f2","f3","f4","f5")),
#'                 strand=as.integer(1))
#'     t1 =findOverlappingPeaks(peaks1, peaks2, maxgap=1000, 
#'           NameOfPeaks1="TF1", NameOfPeaks2="TF2", select="all", annotate=1) 
#'     r = t1$OverlappingPeaks
#'     pie(table(r$overlapFeature))
#'     as.data.frame(t1$MergedPeaks)
#'     }
#' 
findOverlappingPeaks <- function(Peaks1, Peaks2, maxgap = -1L, minoverlap = 0L,
                                 multiple = c(TRUE, FALSE),
                                 NameOfPeaks1 = "TF1", NameOfPeaks2 = "TF2", 
                                 select = c("all", "first", "last", "arbitrary"), 
                                 annotate = 0, ignore.strand = TRUE, 
                                 connectedPeaks = c("min", "merge"), ...) {
        .Deprecated("findOverlapsOfPeaks")
        ### Check inputs
        NAME_conn_string <- "___conn___"
        NAME_short_string <- "__"
        NAME_long_string <- "///"
        PeaksList <- list(...)
        PeaksList <- lapply(PeaksList, function(Peaks) {
            if (!inherits(Peaks, "GRanges")) {
                stop("No valid Peaks passed in. It needs to be GRanges object")
            }
            if (any(is.na(names(Peaks))) || any(duplicated(names(Peaks)))) {
                message("duplicated or NA names found. 
                        Rename all the names by numbers.")
                n_peaks <- length(Peaks)
                names(Peaks) <- formatC(seq_len(n_peaks), 
                                        width = nchar(n_peaks), 
                                        flag = '0')
            }
            Peaks
        })
        n <- length(PeaksList)
        if (n > 0) {
            if (n == 1) {
                PeaksList <- PeaksList[[1]]
                n <- length(PeaksList)
                names <- names(PeaksList)
                if (is.null(names)) {
                    names <- paste("peaks", seq.int(n), sep = "")
                }
            } else {
                ## Save dots arguments names
                dots <- substitute(list(...))[-1]
                names <- unlist(sapply(dots, deparse))
            }
        } else {
            names <- NULL
        }
        if ((!missing(Peaks1)) || (!missing(Peaks2))) {
            if (!missing(Peaks2)) {
                if (!inherits(Peaks2, "GRanges")) {
                    stop("No valid Peaks passed in. 
                         It needs to be GRanges object")
                }
                if (n != 0) {
                    PeaksList <- c(list(Peaks2), PeaksList)
                } else {
                    PeaksList <- list(Peaks2)
                }
                n <- n + 1
                if (missing(NameOfPeaks2)) {
                    NameOfPeaks2 <- "Peaks2"
                }
                names <- c(NameOfPeaks2, names)
            }
            if (!missing(Peaks1)) {
                if (!inherits(Peaks1, "GRanges")) {
                    stop("No valid Peaks passed in. 
                         It needs to be GRanges object")
                }
                if (n != 0) {
                    PeaksList <- c(list(Peaks1), PeaksList)
                } else {
                    PeaksList <- list(Peaks1)
                }
                n <- n + 1
                if (missing(NameOfPeaks1)) {
                    NameOfPeaks1 <- "Peaks1"
                }
                names <- c(NameOfPeaks1, names)
            }
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
        ## Handle colnames of metadata
        metacolnames <- lapply(PeaksList, function(Peaks)
            colnames(mcols(Peaks)))
        metacolnames <- Reduce(intersect, metacolnames)
        metacolclass <- do.call(rbind, lapply(PeaksList, function(Peaks)
            sapply(mcols(Peaks)[, metacolnames, drop = FALSE], 
                   function(.ele) class(.ele)[1])))
        metacolclass <- apply(metacolclass, 2, 
                              function(.ele) length(unique(.ele)) == 1)
        metacolnames <- metacolnames[metacolclass]
        PeaksList <- lapply(PeaksList, function(Peaks) {
            mcols(Peaks) <- mcols(Peaks)[, metacolnames]
            Peaks
        })
        ## Get all merged peaks
        for (i in seq_len(n)) {
            names(PeaksList[[i]]) <- 
                paste(names[i], names(PeaksList[[i]]), 
                      sep = NAME_conn_string)
        }
        Peaks <- unlist(GRangesList(PeaksList))
        if (ignore.strand) {
            strand(Peaks) <- "*"
        }
        
        ol <- as.data.frame(findOverlaps(Peaks, maxgap = maxgap, 
                                         minoverlap = minoverlap, 
                                         select = "all",
                                         drop.self = TRUE, 
                                         drop.redundant = TRUE))
        olm <- cbind(names(Peaks[ol[, 1]]), names(Peaks[ol[, 2]]))
        edgeL <- c(split(olm[, 2], olm[, 1]), split(olm[, 1], olm[, 2]))
        nodes <- unique(as.character(olm))
        ## Use graph to extract all the connected peaks
        gR <- new("graphNEL", nodes = nodes, edgeL = edgeL)
        Merged <- connectedComp(ugraph(gR))        
        Left <- as.list(names(Peaks)[!names(Peaks) %in% nodes])
        all <- c(Merged, Left)
        
        ##venn count
        ncontrasts <- n
        noutcomes <- 2^ncontrasts
        outcomes <- matrix(0,noutcomes,ncontrasts)
        colnames(outcomes) <- names
        for (j in seq_len(ncontrasts))
            outcomes[,j] <- rep(0:1,times=2^(j-1), 
                                each=2^(ncontrasts-j))
        xlist <- list()
        xlist1 <- list()
        for (i in seq_len(ncontrasts)){
            xlist[[i]] <- factor(as.numeric(unlist(lapply(all, function(.ele) 
                any(grepl(paste("^",
                                names[ncontrasts-i+1], 
                                NAME_conn_string, sep=""), .ele))))),
                                 levels=c(0,1))
            if(connectedPeaks=="merge"){
                xlist1[[i]] <- xlist[[i]]
            }else{
                xlist1[[i]] <- 
                    factor(as.numeric(unlist(lapply(all, function(.ele) {
                    ##count involved nodes in each group
                    if(length(.ele)>2){
                        .ele <- gsub(paste(NAME_conn_string, ".*?$", sep=""), 
                                     "", .ele)
                        .ele <- table(.ele)
                        rep(names[ncontrasts-i+1] %in% names(.ele), 
                            min(.ele))
                    }else{
                        any(grepl(paste("^",names[ncontrasts-i+1], 
                                        NAME_conn_string, sep=""), 
                                  .ele))
                    }
                }))), levels=c(0,1))
            }
        }
        counts <- as.vector(table(xlist1))
        venn_cnt <- structure(cbind(outcomes, Counts=counts), 
                              class="VennCounts")
        xlist <- do.call(rbind, xlist)
        xlist <- xlist - 1
        xlist <- xlist[rev(seq_len(nrow(xlist))),,drop=FALSE] 
        ## reverse xlist to match the order of names
        xlist <- apply(xlist, 2, base::paste, collapse="")
        all <- do.call(rbind, mapply(function(.ele, .id) cbind(.id, .ele), 
                                     all, seq_along(all), SIMPLIFY=FALSE))
        all.peaks <- Peaks[all[,2]]
        all.peaks$gpForFindOverlapsOfPeaks <- all[, 1]
        all.peaks.rd <- reduce(all.peaks, min.gapwidth=maxgap+1L,
                               with.revmap=TRUE)
        mapping <- all.peaks.rd$revmap
        m <- sapply(mapping, length)
        mIndex <- rep(seq_along(mapping), m)
        mLists <- unlist(mapping)
        mcols <- mcols(all.peaks[mLists])
        mcols$peakNames <- gsub(NAME_conn_string, 
                                NAME_short_string, 
                                names(all.peaks[mLists]))
        mcolsn <- sapply(mcols[1, ], function(.ele) class(.ele)[1])
        mapping <- DataFrame(HHH_row___H=seq_along(all.peaks.rd))
        for(.name in names(mcolsn)){
            .dat <- split(mcols[, .name], mIndex)
            mapping[, .name] <- switch(mcolsn[.name],
                                       logical=LogicalList(.dat),
                                       integer=IntegerList(.dat),
                                       numeric=NumericList(.dat),
                                       character=CharacterList(.dat),
                                       rle=RleList(.dat),
                                       ComplexList(.dat))
        }
        mapping$HHH_row___H <- NULL
        mcols(all.peaks.rd) <- mapping
        names(all.peaks.rd) <- sapply(all.peaks.rd$peakNames, 
                                      base::paste,
                                      collapse=NAME_short_string)
        all.peaks.rd$gpForFindOverlapsOfPeaks <- 
            unlist(lapply(all.peaks.rd$gpForFindOverlapsOfPeaks, unique)) 
        all <- split(all.peaks.rd, all.peaks.rd$gpForFindOverlapsOfPeaks)
        all <- all[order(as.numeric(names(all)))] 
        ##important, and length(all)==length(xlist)
        if(length(all) != length(xlist)) {
            stop("Length of 'all' should be equal to length of 'xlist'. ",
                 "Please report the bug. Thanks.", call. = FALSE)
        }
        listname <- apply(outcomes, 1, 
                          function(id) paste(names[as.logical(id)], 
                                             collapse=NAME_long_string))
        listcode <- apply(outcomes, 1, base::paste, collapse="")
        listname <- listname[-1]
        listcode <- listcode[-1]
        peaklist <- list()
        for(i in seq_along(listcode)){
            sublist <- all[xlist==listcode[i]]
            if(length(sublist)>0) 
                peaklist[[listname[i]]]<-unlist(sublist, use.names=FALSE)
        }
        correlation <- list()
        names(Peaks) <- gsub(NAME_conn_string, NAME_short_string, names(Peaks))
        npl <- names(peaklist)
        for(i in seq_along(npl)){
            npln <- unlist(strsplit(npl[i], NAME_long_string))
            if(length(npln)==2){
                pl <- peaklist[[i]]$peakNames
                pl <- pl[sapply(pl, length)>1]
                pl1 <- lapply(pl, 
                              function(.ele) 
                                  .ele[grepl(paste("^", npln[1], 
                                                   NAME_short_string, 
                                                   sep=""), 
                                             .ele)][1])
                pl2 <- lapply(pl, 
                              function(.ele) 
                                  .ele[grepl(paste("^", npln[2], 
                                                   NAME_short_string, 
                                                   sep=""), .ele)][1])
                idsel <- mapply(function(.p1, .p2) is.na(.p1)+is.na(.p2)==0, 
                                pl1, pl2)
                pl1 <- Peaks[unlist(pl1[idsel])]
                pl2 <- Peaks[unlist(pl2[idsel])]
                cl <- 
                    getRelationship(pl1, pl2)[,c("insideFeature", 
                                                 "shortestDistance")]
                colnames(cl)[grepl("insideFeature", colnames(cl))] <- 
                    "overlapFeature"
                correlation[[npl[i]]] <- 
                    cbind(peaks1=names(pl1), as.data.frame(pl1), 
                          peaks2=names(pl2), as.data.frame(pl2), 
                          cl)
                rownames(correlation[[npl[i]]]) <- 
                    paste(names(pl1), names(pl2), sep="_")
            }
        }
        for(i in seq_along(peaklist)){
            peaklist[[i]]$gpForFindOverlapsOfPeaks <- NULL
        }
        
        if((!missing(Peaks1)) && (!missing(Peaks2))){
            sampleName <- 
                npl==paste(NameOfPeaks1, 
                           NAME_long_string, 
                           NameOfPeaks2, sep="")
            if(!any(sampleName)) 
                sampleName <- npl==paste(NameOfPeaks2, 
                                         NAME_long_string, 
                                         NameOfPeaks1, sep="")
            sampleName <- npl[sampleName]
            Peaks1withOverlaps <- 
                Peaks[names(Peaks) %in% correlation[[sampleName]]$peaks1]
            names(Peaks1withOverlaps) <- 
                gsub(paste(NameOfPeaks1, NAME_short_string, sep=""), 
                     "", 
                     names(Peaks1withOverlaps))
            Peaks2withOverlaps <- 
                Peaks[names(Peaks) %in% correlation[[sampleName]]$peaks2]
            names(Peaks2withOverlaps) <- 
                gsub(paste(NameOfPeaks2, NAME_short_string, sep=""), 
                     "", 
                     names(Peaks2withOverlaps))
            mergedPeaks <- peaklist[[sampleName]]##To fit old version
            mergedPeaks$peakNames <- NULL
            structure(list(venn_cnt=venn_cnt, 
                           peaklist=peaklist, 
                           overlappingPeaks=correlation,
                           OverlappingPeaks=correlation[[sampleName]],
                           MergedPeaks=mergedPeaks,##
                           Peaks1withOverlaps=Peaks1withOverlaps,
                           Peaks2withOverlaps=Peaks2withOverlaps), 
                      class="overlappingPeaks")
        }else{
            structure(list(venn_cnt=venn_cnt, 
                           peaklist=peaklist, 
                           overlappingPeaks=correlation), 
                      class="overlappingPeaks")
        }
    }
