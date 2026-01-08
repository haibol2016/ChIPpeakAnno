#' Obtain genomic sequences around peaks
#' 
#' @description 
#' Extracts genomic DNA sequences around peak regions using either a BSgenome
#' object (for direct sequence retrieval) or a biomaRt Mart object (for
#' annotation-based sequence retrieval). The function extends peaks by
#' specified upstream and downstream offsets and retrieves the corresponding
#' sequences from the genome.
#' 
#' \strong{Two modes of operation:}
#' \itemize{
#'   \item \code{BSgenome} mode: Directly retrieves sequences from a BSgenome
#'         object. Extends peaks by \code{upstream} and \code{downstream} offsets,
#'         with strand-aware handling (for negative strand peaks, upstream and
#'         downstream are swapped relative to the genomic coordinates).
#'   \item \code{Mart} mode: Uses biomaRt to retrieve sequences based on
#'         annotation. Requires \code{AnnotationData} (or queries biomaRt if
#'         missing). Annotates peaks first, then retrieves sequences relative to
#'         annotated features (e.g., TSS).
#' }
#' 
#' @param myPeakList A \code{GRanges} object containing peak regions. If peaks
#'        have no names, they will be automatically generated as "peak_1", "peak_2", etc.
#' @param upstream An integer specifying the upstream offset from peak start (in
#'        base pairs). Default is \code{200L}. For negative strand peaks, this
#'        offset is applied from the peak end (which is the 5' end on the negative
#'        strand).
#' @param downstream An integer specifying the downstream offset from peak end (in
#'        base pairs). Default is equal to \code{upstream}. For negative strand
#'        peaks, this offset is applied from the peak start.
#' @param genome Either a \code{BSgenome} object (for direct sequence
#'        retrieval) or a \code{Mart} object from biomaRt (for annotation-based
#'        retrieval). See \code{\link[BSgenome]{available.genomes}} and
#'        \code{\link[biomaRt]{useMart}} for details.
#' @param AnnotationData Optional \code{GRanges} object with annotation
#'        information. Required when \code{genome} is a Mart object. If missing
#'        and \code{genome} is a Mart object, the function will query biomaRt
#'        using \code{\link{getAnnotation}}, but it is recommended to call
#'        \code{getAnnotation} beforehand for better performance.
#' 
#' @return Returns a \code{GRanges} object with the same structure as the input
#'        peaks, with additional metadata columns:
#'        \itemize{
#'          \item \code{upstream}: Upstream offset used (integer)
#'          \item \code{downstream}: Downstream offset used (integer)
#'          \item \code{sequence}: Retrieved DNA sequence as a character string.
#'                Peaks that could not be retrieved (e.g., due to chromosome
#'                mismatches or boundary issues) will have \code{NA} for sequence.
#'        }
#'        
#'        \strong{BSgenome mode:}
#'        \itemize{
#'          \item Returns all input peaks with sequence metadata
#'          \item Peaks that extend beyond chromosome boundaries are trimmed to
#'                fit within chromosome limits
#'          \item Peaks that start beyond chromosome boundaries (start > chromosome
#'                length) are removed
#'          \item Strand information is preserved, but sequences are always
#'                returned in the forward strand orientation (5' to 3')
#'          \item If chromosome names don't match between peaks and genome, the
#'                function attempts to format them using \code{formatSeqnames}
#'        }
#'        
#'        \strong{Mart mode:}
#'        \itemize{
#'          \item Returns only peaks that were successfully annotated and had
#'                sequences retrieved
#'          \item Strand is set to "+" for all returned peaks
#'          \item Additional metadata from annotation may be included
#'        }
#' 
#' @details
#' 
#' \strong{Strand-aware extension (BSgenome mode):}
#' 
#' For positive strand peaks:
#' \itemize{
#'   \item Extension: [peak_start - upstream, peak_end + downstream]
#' }
#' 
#' For negative strand peaks:
#' \itemize{
#'   \item Extension: [peak_start - downstream, peak_end + upstream]
#'   \item This ensures that "upstream" refers to the 5' direction relative to
#'         the gene/feature orientation
#' }
#' 
#' \strong{Boundary handling:}
#' 
#' Peaks that extend beyond chromosome boundaries are handled as follows:
#' \itemize{
#'   \item If the extended start position is <= 0, it is set to 1
#'   \item If the extended end position exceeds the chromosome length, it is
#'         trimmed to the chromosome length
#'   \item Peaks where start > end after trimming are removed
#' }
#' 
#' \strong{Performance considerations:}
#' \itemize{
#'   \item For Mart mode, it is recommended to call \code{\link{getAnnotation}}
#'         beforehand and pass the result as \code{AnnotationData} to avoid
#'         repeated database queries
#'   \item BSgenome mode is generally faster as it directly accesses local
#'         genome sequences
#' }
#' @author Lihua Julie Zhu, Jianhong Ou
#' @references Durinck S. et al. (2005) BioMart and Bioconductor: a powerful
#' link between biological biomarts and microarray data analysis.
#' Bioinformatics, 21, 3439-3440.
#' @keywords misc
#' @export
#' @importFrom BiocGenerics start end width strand
#' @importFrom GenomeInfoDb seqlengths seqlevels `seqlengths<-`
#' @importFrom Biostrings getSeq
#' @examples
#' 
#' #### use Annotation data from BSgenome
#' peaks <- GRanges(seqnames=c("NC_008253", "NC_010468"),
#'                  IRanges(start=c(100, 500), end=c(300, 600), 
#'                          names=c("peak1", "peak2")))
#' library(BSgenome.Ecoli.NCBI.20080805)
#' seq <- getAllPeakSequence(peaks, upstream=20, downstream=20, genome=Ecoli)
#' write2FASTA(seq, file="test.fa")
#' 
getAllPeakSequence <- function(myPeakList, 
                               upstream = 200L, downstream = upstream, 
                               genome, AnnotationData) {
    if (!inherits(myPeakList, "GRanges")) {
        stop("'myPeakList' must be a GRanges object", call. = FALSE)
    }
    if (missing(genome)) {
        stop("'genome' is required. Please pass in either a BSgenome object ",
             "or a Mart object", call. = FALSE)
    }
    
    old_name <- names(myPeakList)
    if (length(old_name) != length(myPeakList)) {
        names(myPeakList) <- paste0("peak_", seq_along(myPeakList))
    }
    myPeakList_bk <- myPeakList
    if (inherits(genome, "BSgenome")) {
        strand_info <- strand(myPeakList)
        width_info <- width(myPeakList)
        
        start(myPeakList) <- ifelse(strand_info == "-", 
                                    start(myPeakList) - as.numeric(downstream),
                                    start(myPeakList) - as.numeric(upstream))
        lt0 <- start(myPeakList) <= 0L
        start_adj <- start(myPeakList) - 1L
        start_adj[!lt0] <- 0L
        start(myPeakList)[lt0] <- 1L
        width(myPeakList) <- width_info + as.numeric(upstream) + 
            as.numeric(downstream) + start_adj
        strand_info[strand_info != "-"] <- "+"
        strand(myPeakList) <- strand_info
        
        # Format seqnames to match genome
        if (!all(seqlevels(myPeakList) %in% seqnames(genome))) {
            genome <- formatSeqnames(genome, myPeakList)
        }
        
        ends <- vapply(seq_along(myPeakList), function(i) {
            chr <- as.character(seqnames(myPeakList)[i])
            min(end(myPeakList)[i], 
                seqlengths(genome)[chr], na.rm = TRUE)
        }, FUN.VALUE = integer(1))
        ends <- unname(ifelse(is.na(ends), end(myPeakList), ends))
        keep <- start(myPeakList) <= ends
        end(myPeakList)[keep] <- ends[keep]
        myPeakList <- myPeakList[keep]
        seq <- getSeq(genome, myPeakList, as.character = TRUE)
        
        myPeakList <- myPeakList_bk
        myPeakList$upstream <- rep(upstream, length(myPeakList_bk))
        myPeakList$downstream <- rep(downstream, length(myPeakList_bk))
        myPeakList$sequence <- NA_character_
        if (!all(names(myPeakList) %in% names(seq))) {
            warning("The genome assembly may not be identical to your peaks!",
                    call. = FALSE, immediate. = TRUE)
        }
        myPeakList[names(seq)]$sequence <- seq
        if (all(seqlevels(myPeakList) %in% seqlevels(genome))) {
            seqlengths(myPeakList) <- seqlengths(genome)[seqlevels(myPeakList)]
        }   
        names(myPeakList) <- old_name
        myPeakList
    } else if (inherits(genome, "Mart")) {
        if (missing(AnnotationData)) {
            message("No AnnotationData as GRanges is passed in, ",
                    "so now querying biomart database for AnnotationData ....")
            AnnotationData <- getAnnotation(genome)
            message("Done querying biomart database, start annotating .... ",
                    "Better way would be calling getAnnotation before ",
                    "querying for sequence")
        }
        if (!inherits(AnnotationData, "GRanges")) {
            stop("AnnotationData needs to be a GRanges object. ",
                 "Better way would be calling getAnnotation.", call. = FALSE)
        }
        
        downstream.bk <- downstream
        plusAnno <- AnnotationData[strand(AnnotationData) == "+"]
        temp <- annotatePeakInBatch(myPeakList, AnnotationData = plusAnno)
        TSSlength <- temp$end_position - temp$start_position
        downstream <- end(temp) - start(temp) + downstream
        temp$downstream <- downstream
        temp$TSSlength <- TSSlength
        
        myList3 <- as.data.frame(temp)
        rm(temp)
        
        if (nrow(myList3) > 0L) {
            l3 <- cbind(as.character(myList3$feature), 
                        as.numeric(as.character(myList3$distancetoFeature)), 
                        as.numeric(rep(upstream, nrow(myList3))), 
                        as.numeric(as.character(myList3$downstream)), 
                        as.numeric(as.character(myList3$start_position)), 
                        as.numeric(as.character(myList3$end_position)))
            r3 <- apply(l3, 1, getGeneSeq, genome)
        } else {
            r3 <- 0
        }
        
        if (is.list(r3)) {
            r <- as.data.frame(do.call("rbind", r3))
        } else {
            stop("No sequence found error!")
        }
        colnames(r) <- c("feature", "distancetoFeature", "upstream", 
                       "downstream", "seq")
        r4 <- merge(r, myList3)
        GRanges(seqnames = as.character(r4$space),
                ranges = IRanges(start = r4$start, 
                               end = r4$end, 
                               names = as.character(r4$name)),
                strand = "+",
                upstream = rep(upstream, nrow(r4)), 
                downstream = rep(downstream.bk, nrow(r4)),
                sequence = unlist(r4$seq))
    } else {
        stop("'genome' must be either a BSgenome object or Mart object!", 
             call. = FALSE)
    }
}
