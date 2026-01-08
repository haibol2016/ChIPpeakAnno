#' Write sequences to a FASTA file
#' 
#' @description 
#' Writes genomic sequences from a \code{GRanges} object (typically obtained
#' from \code{\link{getAllPeakSequence}}) to a file in FASTA format. FASTA is
#' a standard file format for biological sequence data where each sequence is
#' preceded by a header line starting with \code{>}. The function leverages
#' \code{writeXStringSet} from the Biostrings package to handle the actual
#' writing operation.
#' 
#' @param mySeq A \code{GRanges} object that must contain a \code{sequence}
#'        metadata column with DNA sequences (typically obtained from
#'        \code{\link{getAllPeakSequence}}). The sequences should be character
#'        strings or \code{DNAString} objects. If the object has names (from
#'        \code{names(mySeq)}), they will be used as sequence identifiers in
#'        the FASTA headers. If names are missing, identifiers will be
#'        automatically generated from genomic coordinates in the format:
#'        \code{"X001_chr1:1000-2000:+"} (where the number is zero-padded,
#'        followed by chromosome, start, end, and strand).
#' @param file A character string naming the output file path, or a connection
#'        object. If \code{""} (default), writes to standard output (console).
#'        The file will be created if it doesn't exist, or overwritten if it
#'        does. Can also be a connection opened for writing (e.g., from
#'        \code{file()} or \code{url()}).
#' @param width An integer specifying the maximum number of characters per line
#'        in the output FASTA file. Sequences longer than this will be wrapped
#'        across multiple lines. Default is \code{80}, which is the standard
#'        FASTA format line width. Set to a larger value (e.g., \code{1000}) for
#'        single-line sequences, or a smaller value for more compact output.
#' @return Invisibly returns the file path (as a character string) or connection
#'        object that was used for writing. The function primarily writes
#'        sequences in FASTA format to the specified destination. Each sequence
#'        is written with a header line starting with \code{>} followed by the
#'        sequence identifier, then the sequence itself (possibly wrapped
#'        according to \code{width}).
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that \code{mySeq} is a \code{GRanges} object
#'   \item Checks that the \code{sequence} metadata column exists
#'   \item Determines sequence identifiers: uses \code{names(mySeq)} if available,
#'         otherwise generates identifiers from genomic coordinates
#'   \item Converts sequences to an \code{XStringSet} object
#'   \item Writes to the specified file or connection using
#'         \code{writeXStringSet}
#' }
#' 
#' The auto-generated identifiers have the format:
#' \code{"X<zero-padded-index>_<seqname>:<start>-<end>:<strand>"}
#' 
#' For example: \code{"X001_chr1:1000-2000:+"} or \code{"X042_chr2:5000-5100:-"}
#' 
#' @author Lihua Julie Zhu
#' @seealso \code{\link{getAllPeakSequence}} for obtaining sequences from
#'          genomic ranges, \code{\link[Biostrings:writeXStringSet]{writeXStringSet}}
#'          for the underlying writing function
#' @keywords misc
#' @export
#' @importFrom Biostrings writeXStringSet
#' @examples
#' ## Example 1: Write sequences with names
#' peaksWithSequences <- GRanges(seqnames = c("1", "2"),
#'                                IRanges(start = c(1000, 2000), 
#'                                        end = c(1010, 2010), 
#'                                        names = c("peak1", "peak2")), 
#'                                sequence = c("CCCCCCCCGGGGG", "TTTTTTTAAAAAA"))
#' 
#' write2FASTA(peaksWithSequences, file = "testseq.fasta", width = 50)
#' 
#' ## Example 2: Write sequences without names (auto-generated identifiers)
#' peaksNoNames <- GRanges(seqnames = c("chr1", "chr2"),
#'                          IRanges(start = c(1000, 2000), 
#'                                  end = c(1010, 2010)), 
#'                          strand = c("+", "-"),
#'                          sequence = c("ATGCATGCATGC", "GCATGCATGCAT"))
#' 
#' write2FASTA(peaksNoNames, file = "testseq2.fasta")
#' 
#' ## Example 3: Write to standard output (console)
#' write2FASTA(peaksWithSequences, file = "", width = 60)
#' 
write2FASTA <- function(mySeq, file = "", width = 80) {
    if (!inherits(mySeq, "GRanges")) {
        stop("'mySeq' must be a GRanges object", call. = FALSE)
    }
    
    if (is.null(mySeq$sequence)) {
        stop("metadata must contain 'sequence' column", call. = FALSE)
    }
    
    if (!is.null(names(mySeq))) {
        descriptions <- names(mySeq)
    } else {
        n_seq <- length(mySeq)
        descriptions <- paste0(
            "X", 
            formatC(seq_len(n_seq), 
                    width = nchar(as.character(n_seq)),
                    flag = "0"), 
            "_", 
            as.character(seqnames(mySeq)), ":", 
            start(mySeq), "-",
            end(mySeq), ":",
            as.character(strand(mySeq))
        )
    }
    
    sequences <- mySeq$sequence
    names(sequences) <- descriptions
    writeXStringSet(as(sequences, "XStringSet"), file, width = width)
}
