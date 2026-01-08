#' Count sequences containing a DNA pattern
#' 
#' @description 
#' Counts the number of sequences that contain at least one occurrence of the
#' specified DNA pattern. The function searches for the pattern on both the
#' forward strand and the reverse complement strand, counting a sequence as
#' positive if the pattern is found on either strand.
#' 
#' IUPAC nucleotide ambiguity codes in the pattern are automatically converted
#' to regular expressions using \code{\link{translatePattern}}. Supported codes
#' include: Y (C or T), R (A or G), S (G or C), W (A or T), K (T or G), M (A or C),
#' B (C, G, or T), D (A, G, or T), H (A, C, or T), V (A, C, or G), and N (any
#' nucleotide).
#' 
#' @param pattern A \code{DNAStringSet} object containing the DNA pattern to
#'        search for. If the object contains multiple patterns, only the first
#'        pattern (element) is used. IUPAC nucleotide ambiguity codes are
#'        supported and will be converted to regular expressions for matching.
#' @param sequences A character vector of DNA sequences to search in. Each
#'        element should be a DNA sequence string (A, T, C, G, and IUPAC codes).
#' @return Returns an integer representing the total number of sequences that
#'        contain at least one occurrence of the pattern. A sequence is counted
#'        if the pattern (or its reverse complement) is found anywhere within
#'        that sequence. The function counts sequences, not the total number of
#'        pattern occurrences (i.e., a sequence with multiple matches is still
#'        counted only once).
#' @details
#' 
#' \strong{Algorithm:}
#' \enumerate{
#'   \item The pattern is converted to a character string and translated to a
#'         regular expression using \code{\link{translatePattern}} to handle
#'         IUPAC ambiguity codes
#'   \item The reverse complement of the pattern is computed and also translated
#'         to a regular expression
#'   \item For each sequence, the function searches for both the forward pattern
#'         and the reverse complement pattern using \code{regexpr()}
#'   \item A sequence is counted if either pattern is found (position > 0)
#' }
#' 
#' \strong{Note:} This function counts sequences with matches, not the total
#' number of matches. If you need to count total occurrences or get match
#' positions, consider using \code{\link{summarizePatternInPeaks}} or
#' \code{\link{getAllPeakSequence}} with pattern matching functions.
#' 
#' @author Lihua Julie Zhu
#' @seealso \code{\link{summarizePatternInPeaks}} for detailed pattern
#'          occurrence analysis in peak regions, \code{\link{translatePattern}}
#'          for IUPAC code translation, \code{\link{getAllPeakSequence}} for
#'          extracting sequences from peak regions
#' @keywords misc
#' @export
#' @importFrom Biostrings reverseComplement
#' @examples
#' \dontrun{
#'   library(Biostrings)
#'   
#'   ## Example 1: Using patterns from a FASTA file
#'   filepath <- system.file("extdata", "examplePattern.fa", 
#'                          package = "ChIPpeakAnno")
#'   dict <- readDNAStringSet(filepath = filepath, format = "fasta", 
#'                           use.names = TRUE)
#'   sequences <- c("ACTGGGGGGGGCCTGGGCCCCCAAAT", 
#'                  "AAAAAACCCCTTTTGGCCATCCCGGGACGGGCCCAT", 
#'                  "ATCGAAAATTTCC")
#'   countPatternInSeqs(pattern = dict[1], sequences = sequences)
#'   countPatternInSeqs(pattern = dict[2], sequences = sequences)
#'   
#'   ## Example 2: Using IUPAC ambiguity codes
#'   pattern <- DNAStringSet("ATNGMAA")  # N = any nucleotide, M = A or C
#'   countPatternInSeqs(pattern = pattern, sequences = sequences)
#'   
#'   ## Example 3: Pattern with multiple ambiguity codes
#'   pattern <- DNAStringSet("ATRYW")  # R = A or G, Y = C or T, W = A or T
#'   countPatternInSeqs(pattern = pattern, sequences = sequences)
#' }
#' 
countPatternInSeqs <- function(pattern, sequences) {
    if (missing(pattern) || !inherits(pattern, "DNAStringSet")) {
        stop("'pattern' is required as a DNAStringSet object", call. = FALSE)
    }
    
    if (missing(sequences) || length(sequences) == 0L) {
        stop("'sequences' must be provided", call. = FALSE)
    }
    
    revcomp_pattern <- reverseComplement(pattern)
    pattern_char <- as.character(pattern)[[1L]]
    pattern_regex <- translatePattern(pattern_char)
    revcomp_pattern_char <- as.character(revcomp_pattern)[[1L]]
    revcomp_pattern_regex <- translatePattern(revcomp_pattern_char)
    
    total <- 0L
    for (i in seq_along(sequences)) {
        pos_plus <- regexpr(pattern_regex, sequences[i], perl = TRUE)[1L]
        pos_minus <- regexpr(revcomp_pattern_regex, sequences[i], perl = TRUE)[1L]
        if (pos_plus > 0L || pos_minus > 0L) {
            total <- total + 1L
        }
    }
    
    total
}
