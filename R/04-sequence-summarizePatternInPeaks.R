#' Summarize pattern occurrence and enrichment in peak sequences
#' 
#' @description 
#' Analyzes the occurrence and statistical enrichment of DNA patterns (motifs)
#' in genomic peak sequences. The function extracts sequences from peak regions,
#' searches for patterns (with optional reverse complement search), and performs
#' statistical enrichment analysis. Two enrichment test methods are supported:
#' binomial test (with expected frequency calculated using naive or Markov chain
#' models) and permutation test (with chromosome-based or shuffled background).
#' 
#' Patterns can be specified using IUPAC nucleotide codes (A, C, G, T, R, Y, S,
#' W, K, M, B, D, H, V, N) and are automatically expanded to all possible
#' sequences. The function returns both enrichment statistics and detailed
#' occurrence information for each pattern match.
#'
#' @param patternFilePath A character string specifying the path to the file
#'        containing DNA patterns. The file must be in FASTA or FASTQ format.
#'        Each sequence in the file represents a pattern to search for. Pattern
#'        names are taken from the sequence names in the file.
#' 
#' @param format A character string specifying the file format. Must be either
#'        \code{"fasta"} (default) or \code{"fastq"}.
#' 
#' @param BSgenomeName A \code{BSgenome} object containing the reference genome
#'        sequences. Required for extracting peak sequences and for permutation
#'        testing with chromosome-based background. See
#'        \code{\link[BSgenome]{available.genomes}} for available genomes.
#' 
#' @param peaks A \link[GenomicRanges:GRanges-class]{GRanges} object containing
#'        the genomic regions (peaks) to analyze. Sequences are extracted from
#'        these regions using the specified \code{BSgenomeName}. Also accepts
#'        \code{RangedData} objects (automatically converted to GRanges).
#' 
#' @param revcomp A logical value. If \code{TRUE} (default), searches for both
#'        the forward pattern and its reverse complement. If \code{FALSE}, only
#'        searches for the forward pattern. Reverse complement patterns are
#'        automatically generated using \code{reverseComplement()}.
#' 
#' @param method A character string specifying the enrichment test method.
#'        Options:
#'        \itemize{
#'          \item \code{"binom.test"} (default): Uses binomial test to assess
#'                enrichment. Requires expected frequency calculation (see
#'                \code{expectFrequencyMethod}).
#'          \item \code{"permutation.test"}: Uses permutation testing to assess
#'                enrichment. Generates background sequences and compares
#'                observed pattern frequency to the distribution from permutations.
#'        }
#' 
#' @param expectFrequencyMethod A character string specifying the method for
#'        calculating expected pattern frequency (only used when
#'        \code{method = "binom.test"}). Options:
#'        \itemize{
#'          \item \code{"Markov"} (default): Uses Markov chain model to calculate
#'                expected frequency. More accurate for sequences with dependencies
#'                between nucleotides. Requires pattern length > 3 and < 13.
#'          \item \code{"Naive"}: Uses naive independence model (product of
#'                individual nucleotide frequencies). Simpler but less accurate
#'                for sequences with dependencies.
#'        }
#' 
#' @param MarkovOrder An integer specifying the order of the Markov chain model
#'        (only used when \code{expectFrequencyMethod = "Markov"}). Must be > 0,
#'        < 6, and < (pattern length - 2). Default is \code{3L}. Higher orders
#'        capture more complex dependencies but require more data.
#' 
#' @param bgdForPerm A character string specifying the method for generating
#'        background sequences for permutation testing (only used when
#'        \code{method = "permutation.test"}). Options:
#'        \itemize{
#'          \item \code{"chromosome"} (default): Selects random regions from
#'                chromosomes as background. Chromosome selection is controlled by
#'                the \code{chromosome} parameter.
#'          \item \code{"shuffle"}: Shuffles the peak sequences to create
#'                background while preserving k-mer composition. Uses
#'                \code{\link[universalmotif]{shuffle_sequences}} function.
#'                Additional parameters can be passed via \code{...}.
#'        }
#'
#' @param chromosome A character string specifying how to select chromosomes for
#'        background generation (only used when \code{bgdForPerm = "chromosome"}).
#'        Options:
#'        \itemize{
#'          \item \code{"asPeak"} (default): Uses the same chromosomes as the
#'                input peaks. Background regions are randomly selected from these
#'                chromosomes.
#'          \item \code{"random"}: Randomly selects chromosomes from all
#'                available chromosomes in the genome. Background regions are
#'                then selected from these chromosomes.
#'        }
#' 
#' @param nperm An integer specifying the number of permutations to perform
#'        (only used when \code{method = "permutation.test"}). Default is
#'        \code{1000}. Higher values provide more accurate p-values but take
#'        longer to compute. Must satisfy \code{round(nperm * alpha, 0) != 0}.
#' 
#' @param alpha A numeric value specifying the significance level for permutation
#'        testing (only used when \code{method = "permutation.test"}). Default
#'        is \code{0.05}. The cutoff is calculated as the value at the
#'        \code{round(nperm * alpha, 0)}-th position in the sorted permutation
#'        distribution. Must satisfy \code{round(nperm * alpha, 0) != 0}.
#'
#' @param ... Additional parameters passed to
#'        \code{\link[universalmotif]{shuffle_sequences}} when
#'        \code{bgdForPerm = "shuffle"}. See the documentation of that function
#'        for available options.
#'
#' @return Returns a list containing two data frames:
#'        \itemize{
#'          \item \code{motif_enrichment}: A data frame with enrichment statistics
#'                for each pattern. Columns:
#'                \itemize{
#'                  \item \code{patternNum}: Number of pattern matches found in
#'                        peak sequences (forward + reverse complement if
#'                        \code{revcomp = TRUE})
#'                  \item \code{totalNumPatternWithSameLen}: Total number of
#'                        possible positions for patterns of the same length
#'                        (i.e., total sequence length - pattern length + 1)
#'                  \item \code{expectedRate}: Expected frequency of pattern
#'                        occurrence (only for \code{method = "binom.test"}).
#'                        Calculated using \code{expectFrequencyMethod}.
#'                  \item \code{patternRate}: Observed pattern frequency
#'                        (\code{patternNum / totalNumPatternWithSameLen}).
#'                        Only for \code{method = "permutation.test"}.
#'                  \item \code{pValueBinomTest}: P-value from binomial test
#'                        (only for \code{method = "binom.test"}). Tests whether
#'                        observed frequency is significantly greater than expected.
#'                  \item \code{cutOffPermutationTest}: Cutoff value from
#'                        permutation test (only for \code{method =
#'                        "permutation.test"}). Pattern is considered enriched if
#'                        \code{patternRate > cutOffPermutationTest}.
#'                }
#'          \item \code{motif_occurrence}: A data frame with detailed information
#'                about each pattern match. One row per match. Columns:
#'                \itemize{
#'                  \item \code{motifChr}: Chromosome where the motif is located
#'                        (same as \code{peakChr})
#'                  \item \code{motifStartInChr}: Start position of the motif in
#'                        chromosome coordinates
#'                  \item \code{motifEndInChr}: End position of the motif in
#'                        chromosome coordinates
#'                  \item \code{motifName}: Name of the pattern (from FASTA/FASTQ
#'                        file)
#'                  \item \code{motifPattern}: The pattern sequence that was
#'                        matched
#'                  \item \code{motifStartInPeak}: Start position of the motif
#'                        relative to the peak start (1-based)
#'                  \item \code{motifEndInPeak}: End position of the motif
#'                        relative to the peak start (1-based)
#'                  \item \code{motifFound}: The actual DNA sequence found in the
#'                        peak (may differ from pattern if IUPAC codes were used)
#'                  \item \code{motifFoundStrand}: Strand where the motif was
#'                        found. \code{"+"} indicates forward strand match,
#'                        \code{"-"} indicates reverse complement match
#'                  \item \code{peakChr}: Chromosome of the peak
#'                  \item \code{peakStart}: Start position of the peak
#'                  \item \code{peakEnd}: End position of the peak
#'                  \item \code{peakWidth}: Width of the peak (in base pairs)
#'                  \item \code{peakStrand}: Strand of the peak
#'                }
#'        }
#' 
#' @details 
#' \strong{Pattern Matching:}
#' \itemize{
#'   \item Patterns are matched using IUPAC nucleotide codes, which are
#'         automatically expanded to all possible sequences (e.g., "R" matches
#'         "A" or "G")
#'   \item Pattern matching uses Perl regular expressions via \code{translatePattern()}
#'   \item If \code{revcomp = TRUE}, both forward and reverse complement patterns
#'         are searched, and matches are labeled accordingly
#' }
#' 
#' \strong{Enrichment Testing:}
#' \itemize{
#'   \item \code{method = "binom.test"}: Compares observed pattern frequency to
#'         expected frequency using binomial test. Expected frequency is
#'         calculated from peak sequences using either naive or Markov chain
#'         models.
#'   \item \code{method = "permutation.test"}: Generates background sequences
#'         and compares observed frequency to the distribution from permutations.
#'         More computationally intensive but does not require assumptions about
#'         sequence composition.
#' }
#' 
#' \strong{Background Generation for Permutation Test:}
#' \itemize{
#'   \item \code{bgdForPerm = "chromosome"}: Randomly selects regions from
#'         chromosomes. Preserves chromosome-specific composition but not local
#'         sequence structure.
#'   \item \code{bgdForPerm = "shuffle"}: Shuffles peak sequences while
#'         preserving k-mer composition. Preserves local sequence structure but
#'         not chromosome-specific composition. See
#'         \code{\link[universalmotif]{shuffle_sequences}} for details.
#' }
#' 
#' \strong{Markov Chain Model:}
#' \itemize{
#'   \item Only applicable when \code{expectFrequencyMethod = "Markov"}
#'   \item Pattern length must be > 3 and < 13
#'   \item Markov order must be < (pattern length - 2)
#'   \item Uses \code{oligoFrequency()} to calculate k-mer frequencies from peak
#'         sequences
#' }
#' 
#' @note 
#' \itemize{
#'   \item For permutation testing, ensure \code{round(nperm * alpha, 0) != 0}
#'         to get meaningful results
#'   \item Peak sequences must be shorter than chromosome lengths when using
#'         \code{bgdForPerm = "chromosome"}
#'   \item The function automatically handles IUPAC codes and expands patterns
#'         to all possible sequences
#'   \item If no patterns are found, \code{motif_occurrence} will be empty or
#'         contain a message
#' }
#' 
#' @author Lihua Julie Zhu, Junhui Li, Kai Hu
#' 
#' @export
#' 
#' @importFrom Biostrings readDNAStringSet reverseComplement oligonucleotideFrequency
#' @importFrom stats binom.test
#' @importFrom universalmotif shuffle_sequences 
#' @importFrom data.table data.table
#' @importFrom stringr str_split
#' @importFrom dplyr %>% pull left_join
#' @importFrom tidyr separate_rows
#' @importFrom tibble tibble
#' @keywords misc
#' @examples
#'                             
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' filepath <- system.file("extdata", "examplePattern.fa", 
#'                         package = "ChIPpeakAnno")
#' peaks <- GRanges(seqnames = c("chr17", "chr3", "chr12", "chr8"),
#'                  IRanges(start = c(41275784, 10076141, 4654135, 31024288),
#'                          end = c(41276382, 10076732, 4654728, 31024996),
#'                          names = paste0("peak", 1:4)))
#' result <- summarizePatternInPeaks(patternFilePath = filepath, peaks = peaks,
#'                                   BSgenomeName = Hsapiens)
#' 
summarizePatternInPeaks <- function(patternFilePath, 
                                    format = "fasta", 
                                    BSgenomeName,
                                    peaks,
                                    revcomp = TRUE,
                                    method = c("binom.test", "permutation.test"),
                                    expectFrequencyMethod = c("Markov", "Naive"),
                                    MarkovOrder = 3L,
                                    bgdForPerm = c("shuffle", "chromosome"),
                                    chromosome = c("asPeak", "random"),
                                    nperm = 1000,
                                    alpha = 0.05,
                                    ...) {
    method <- match.arg(method)
    expectFrequencyMethod <- match.arg(expectFrequencyMethod)
    bgdForPerm <- match.arg(bgdForPerm)
    chromosome <- match.arg(chromosome)
    
    if (missing(patternFilePath)) {
        stop("Missing required parameter 'patternFilePath'!", call. = FALSE)
    }
    if (!file.exists(patternFilePath)) {
        stop("patternFilePath specified as '", patternFilePath,
             "' does not exist!", call. = FALSE)
    }
    if (format != "fasta" && format != "fastq") {
        stop("'format' must be either 'fasta' or 'fastq'!", call. = FALSE)
    }
    if (missing(BSgenomeName) || !inherits(BSgenomeName, "BSgenome")) {
        stop("BSgenomeName is required as BSgenome object!", call. = FALSE)
    }
    if (missing(peaks)) {
        stop("Missing required parameter 'peaks'!", call. = FALSE)
    }
    if (!inherits(peaks, c("RangedData", "GRanges"))) {
        stop("'peaks' must be a RangedData or GRanges object", call. = FALSE)
    }
    if (inherits(peaks, "RangedData")) {
        peaks <- toGRanges(peaks, format = "RangedData")
    }

    seqs <- getAllPeakSequence(peaks, upstream = 0, downstream = 0, 
                               genome = BSgenomeName)
    
    # seqs <- getAllPeakSequence(peaks, BSgenome.Hsapiens.UCSC.hg19, upstream = 0, downstream = 0)
    # patternVec <- DNAStringSet("AACCCA")
    # names(patternVec) <- "testRun"
    # revcomp <- FALSE
    
    patternVec <- readDNAStringSet(patternFilePath, format, use.names = TRUE)
    temp <- do.call(rbind, lapply(seq_along(patternVec), function(i) {
        getPosInSeqs(thispattern = patternVec[i], seq = seqs, revcomp = revcomp)
    }))
    temp2 <- patternEnrichmentInPeak(patternVec,
                                     revcomp = revcomp,
                                     BSgenomeName = BSgenomeName,
                                     seqs = seqs, 
                                     method = method, 
                                     expectFrequencyMethod = expectFrequencyMethod,
                                     bgdForPerm = bgdForPerm,
                                     chromosome = chromosome,
                                     MarkovOrder = MarkovOrder,
                                     nperm = nperm,
                                     alpha = alpha,
                                     ...)
    
    output <- list()
    output["motif_enrichment"] <- list(temp2)
    output["motif_occurrence"] <- list(temp)
    return(output)
}

patternEnrichmentInPeak <- function(patternVec,
                                    revcomp = TRUE,
                                    seqs,
                                    method = c("binom.test", "permutation.test"),
                                    expectFrequencyMethod = c("Markov", "Naive"),
                                    bgdForPerm = c("chromosome", "shuffle"),
                                    chromosome = c("asPeak", "random"),
                                    BSgenomeName,
                                    MarkovOrder = 3L,
                                    nperm = 1000,
                                    alpha = 0.05,
                                    ...) {
    method <- match.arg(method)
    expectFrequencyMethod <- match.arg(expectFrequencyMethod)
    bgdForPerm <- match.arg(bgdForPerm)
    stopifnot("The 'nperm' or 'alpha' parameter should be increased to a 
            sufficient extent." = round(nperm * alpha, 0) != 0)
    
    pattern <- as.character(patternVec)
    if (method == "binom.test") {
        peakPvalue <- binomEnrichment(pattern, 
                                      expectFrequencyMethod = expectFrequencyMethod,
                                      revcomp = revcomp, 
                                      seqs = seqs, 
                                      MarkovOrder = MarkovOrder)
    } else if (method == "permutation.test") {
        peakPvalue <- permutationEnrichment(pattern, 
                                            revcomp = revcomp, 
                                            seqs = seqs, 
                                            BSgenomeName = BSgenomeName, 
                                            bgdForPerm = bgdForPerm,
                                            chromosome = chromosome,
                                            nperm = nperm, 
                                            alpha = alpha,
                                            ...)
    }
    return(peakPvalue)
}

binomEnrichment <- function(pattern, 
                            expectFrequencyMethod = c("Naive", "Markov"), 
                            revcomp = TRUE, 
                            seqs, 
                            MarkovOrder = 3L) {
    ## get frequency of single nucleotide in seqs
    expFreqSet <- expectMotifFrequency(pattern, 
                                       revcomp = revcomp, 
                                       seqs = seqs, 
                                       expectFrequencyMethod = expectFrequencyMethod,
                                       MarkovOrder = MarkovOrder)
    ## get real frequency of input seq in peak region
    peakFreq <- matrix(0, length(pattern), 3)
    rownames(peakFreq) <- pattern
    colnames(peakFreq) <- c("patternNum", "totalNumPatternWithSameLen",
                            "expectedRate")
    peakFreq[, 3] <- expFreqSet[, 1]
    
    ## get expected frequency of fasta
    pattern.t <- translatePattern(pattern)
    backwardpattern.t <- NULL
    backwardpattern <- NULL
    if (revcomp == TRUE) {
        backwardpattern <- reverseComplement(DNAStringSet(pattern))
        backwardpattern <- as.character(backwardpattern)
        backwardpattern.t <- translatePattern(backwardpattern)
    }
    peakFreq[, 1:2] <- t(sapply(seq.int(pattern), function(i) {
        oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),
                                                      width = nchar(pattern[i])))
        posPlus <- gregexpr(pattern.t[i], names(oligoNVec), perl = TRUE)
        posPlusFreq <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
        posMinusFreq <- 0
        if (revcomp == TRUE) {
            posMinus <- gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
            posMinusFreq <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
        }
        c(posPlusFreq + posMinusFreq, sum(oligoNVec))
    }))
    
    pValueBinomTest <- sapply(seq.int(nrow(peakFreq)), function(i) {
        stats::binom.test(peakFreq[i, 1], peakFreq[i, 2], peakFreq[i, 3],
                          alternative = "greater")$p.value
    })
    peakFreq <- data.frame(peakFreq, pValueBinomTest)
    return(peakFreq)
}

expectMotifFrequency <- function(pattern, 
                                 revcomp = TRUE, 
                                 seqs,
                                 expectFrequencyMethod = c("Naive", "Markov"),
                                 MarkovOrder = 3L) {
    
    if (expectFrequencyMethod == "Naive") {
        ACGTcount <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),
                                                      width = 1))
        ACGTfreq <- ACGTcount / sum(ACGTcount)
        iuapc <- data.table(code = c("A", "C", "G", "T", "R", "Y", "S", "W", "K", 
                                     "M", "B", "D", "H", "V", "N"),
                            base = c("A", "C", "G", "T", "AG", "CT", "GC", "AT", 
                                    "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT"))
        allBaseFreq <- vapply(seq_len(nrow(iuapc)), function(i) {
            sum(ACGTfreq[seqinr::s2c(as.character(iuapc[i, 2]))])
        }, FUN.VALUE = numeric(1))
        names(allBaseFreq) <- as.matrix(iuapc)[, 1]
        
        ## get expected frequency of forward and backward fasta 
        pattern.t <- translatePattern(pattern)
        backwardpattern.t <- NULL
        backwardpattern <- NULL
        if (revcomp == TRUE) {
            backwardpattern <- reverseComplement(DNAStringSet(pattern))
            backwardpattern <- as.character(backwardpattern)
            backwardpattern.t <- translatePattern(backwardpattern)
        }
        expFreqSet <- sapply(c(pattern, backwardpattern), function(i) {
            patternTab <- table(strsplit(i, "")[[1]])
            expFre <- prod(allBaseFreq[names(patternTab)]^patternTab)
            expFre
        })
        expFreqSet <- colSums(matrix(expFreqSet, ncol = length(pattern), byrow = TRUE))
        names(expFreqSet) <- pattern
        expectedFreq <- data.frame(expFreqSet)
        colnames(expectedFreq) <- "expected_frequency"
    } else if (expectFrequencyMethod == "Markov") {
        
        oligoLength <- nchar(pattern)
        stopifnot("motif length should > 3 and < 13 with Markov method" =
                    all(oligoLength > 3 & oligoLength < 13))
        MarkovOrder <- as.integer(MarkovOrder)
        stopifnot("MarkovOrder should > 0 and < 6" = MarkovOrder > 0 & MarkovOrder < 6)
        stopifnot("MarkovOrder should be less than motif length - 2" =
                    all(MarkovOrder < oligoLength - 2))
        # stopifnot("The seqinr package is required." = 
        #             requireNamespace("seqinr", quietly = TRUE))
        if (inherits(seqs, "GRanges")) {
            seqs_content <- seqs$sequence
        } else {
            if (inherits(seqs, "DNAStringSet")) {
                seqs_content <- as.character(seqs)
            } else {
                if (!is.character(seqs)) {
                    stop("seqs must be an object of DNAStringSet or DNAString ",
                         "or output of getAllPeakSequence", call. = FALSE)
                }
            }
        }
        sequence <- tolower(seqs_content)
        len <- length(sequence)
        #i=6
        expectedFreq <- sapply(unique(oligoLength), function(i) {
            oligoWords <- seqinr::words(i)
            dict <- PDict(DNAStringSet(oligoWords))
            sequence.tbl <- 1
            maxsize <- 100000
            if (len > maxsize) {
                sequence.tbl <- table(sequence)
                if (length(sequence.tbl) > maxsize) {
                    stop("Can not handle such huge dataset.")
                }
                cnt <- vcountPDict(dict, subject = DNAStringSet(names(sequence.tbl)), 
                                   max.mismatch = 0, min.mismatch = 0, 
                                   with.indels = FALSE, fixed = TRUE)
            } else {
                cnt <- vcountPDict(dict, subject = DNAStringSet(sequence), 
                                   max.mismatch = 0, min.mismatch = 0, 
                                   with.indels = FALSE, fixed = TRUE)
            }
            
            rownames(cnt) <- oligoWords
            cnt <- t(cnt)
            mergeRevcomp <- function(mat) {
                coln <- colnames(mat)
                index <- nrow(mat)
                map <- c(a = "t", c = "g", g = "c", t = "a")
                revComp <- function(.ele) {
                    paste(map[rev(seqinr::s2c(.ele))], collapse = "")
                }
                coln.rev <- sapply(coln, revComp)
                coln.id <- seq_along(coln)
                coln.rev.id <- match(coln.rev, coln)
                coln.ids <- apply(cbind(coln.id, coln.rev.id), 1, sort)
                coln.ids <- unique(t(coln.ids))
                mat1 <- mat[, coln.ids[, 1]] > 0
                mat2 <- mat[, coln.ids[, 2]] > 0
                if (index == 1) {
                    mat1 <- t(as.data.frame(mat1))
                    mat2 <- t(as.data.frame(mat2))
                }
                mat1 <- colSums(mat1)
                mat2 <- colSums(mat2)
                coln <- ifelse(mat1 > mat2, names(mat1), names(mat2))
                mat <- mat[, coln.ids[, 1]] + mat[, coln.ids[, 2]]
                if (index == 1) {
                    mat <- t(as.data.frame(mat))
                }
                colnames(mat) <- coln
                mat
            }
            
            if (revcomp) cnt <- mergeRevcomp(cnt)
            
            #cnt <- cnt>0
            mode(cnt) <- "logical"
            mode(cnt) <- "numeric"
            if (len > maxsize) cnt <- cnt * as.numeric(sequence.tbl)
            cntSum <- colSums(cnt)
            freqs <- oligoFrequency(sequence, MarkovOrder = MarkovOrder)
            namesFreqs <- unique(c(seqinr::words(1), 
                                   seqinr::words(MarkovOrder),
                                   seqinr::words(MarkovOrder + 1)))
            
            names(freqs) <- tolower(names(freqs))
            stopifnot(all(namesFreqs %in% names(freqs)))
            f <- sapply(names(cntSum), function(.ele) {
                m1 <- substring(.ele, 1:(i - MarkovOrder),
                            (MarkovOrder + 1):i)
                m0 <- substring(.ele, 2:(i - MarkovOrder),
                            (MarkovOrder + 1):(i - 1))
                prod(freqs[m1]) / prod(freqs[m0])
            })
            #i=6
            allPatternList <- expandPattern(pattern[oligoLength %in% i])
            
            expFreqSubset <- sapply(allPatternList, function(seq) {
                sum(f[names(f) %in% tolower(seq)])
            })
            expFreqSubset
        })
        expectedFreq <- as.data.frame(expectedFreq)
        colnames(expectedFreq) <- "expected_frequency"
    }
    return(expectedFreq)
}


expandPattern <- function(pattern) {
    pattern <- toupper(pattern)
    iuapc <- data.table(code = c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M",
                                 "B", "D", "H", "V", "N"),
                        base = c("A", "C", "G", "T", "AG", "CT", "GC", "AT", "GT",
                                 "AC", "CGT", "AGT", "ACT", "ACG", "ACGT"))
    
    allExpPattern <- lapply(pattern, function(seq) {
        allSeqMat <- tibble(seq) %>%
            separate_rows(seq, sep = '(?<=.)(?=.)') %>%
            left_join(iuapc, by = c("seq" = "code")) %>%
            pull(base) %>%
            str_split("") %>%
            expand.grid(stringsAsFactors = FALSE)
        apply(allSeqMat, 1, seqinr::c2s)
    })
    names(allExpPattern) <- pattern
    allExpPattern
}

permutationEnrichment <- function(pattern, 
                                  revcomp = TRUE, 
                                  seqs, 
                                  BSgenomeName, 
                                  bgdForPerm = c("chromosome", "shuffle"), 
                                  chromosome = c("asPeak", "random"),
                                  nperm = 1000, 
                                  alpha = 0.05,
                                  ...) {
    pattern.t <- translatePattern(pattern)
    backwardpattern.t <- NULL
    backwardpattern <- NULL
    if (revcomp == TRUE) {
        backwardpattern <- reverseComplement(DNAStringSet(pattern))
        backwardpattern <- as.character(backwardpattern)
        backwardpattern.t <- translatePattern(backwardpattern)
    }
    allseqLen <- seqlengths(BSgenomeName)
    seqWidth <- width(seqs)
    nPeak <- length(seqs)
    colnames(mcols(seqs))[4] <- "orig.sequence"
    names(seqs$orig.sequence) <- do.call(paste, c(as.data.frame(seqs)[, 1:3], 
                                                  sep = "_"))
    inputSeq <- DNAStringSet(seqs$orig.sequence)
    candiChrom <- names(allseqLen >= max(seqWidth))
    stopifnot("The length of peak sequence should be less than the length of 
  chromosomes" = all(seqWidth < allseqLen[as.character(seqnames(seqs)@values)]))
    
    peakFreq <- matrix(0, length(pattern), 4)
    rownames(peakFreq) <- pattern
    colnames(peakFreq) <- c("patternNum", "totalNumPatternWithSameLen",
                            "patternRate", "cutOffPermutationTest")
    peakFreq[, 1:2] <- t(sapply(seq.int(pattern), function(i) {
        oligoNVec <- colSums(oligonucleotideFrequency(
            DNAStringSet(seqs$orig.sequence), width = nchar(pattern[i])))
        posPlus <- gregexpr(pattern.t[i], names(oligoNVec), perl = TRUE)
        posPlusFreq <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
        posMinusFreq <- 0
        if (revcomp == TRUE) {
            posMinus <- gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
            posMinusFreq <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
        }
        c(posPlusFreq + posMinusFreq, sum(oligoNVec))
    }))
    peakFreq[, 3] <- peakFreq[, 1] / peakFreq[, 2]
    
    ep <- list(...)
    permutationFreq <- do.call(rbind, lapply(seq.int(nperm), 
                                             function(n, ep) {
                                                 if (bgdForPerm == "chromosome") {
                                                     if (chromosome == "random") {
                                                         chrs <- sample(candiChrom, nPeak, replace = TRUE)
                                                         givenSeqLen <- allseqLen[chrs]
                                                     } else if (chromosome == "asPeak") {
                                                         chrs <- as.character(seqnames(seqs))
                                                         givenSeqLen <- allseqLen[chrs]
                                                     }
                                                     startPos <- vapply(givenSeqLen - seqWidth, function(x) {sample(seq.int(x), 1)},
                                                                      numeric(1))
                                                     endPos <- startPos + seqWidth - 1
                                                     backgroudPeak <- GRanges(seqnames = chrs,
                                                                              IRanges(start = startPos,
                                                                                      end = endPos,
                                                                                      names = paste0("peak", seq.int(nPeak))))
                                                     backgroudPeakseq <- getAllPeakSequence(backgroudPeak, upstream = 0,
                                                                                           downstream = 0, 
                                                                                           genome = BSgenomeName)
                                                     
                                                 } 
                                                 else if (bgdForPerm == "shuffle") {
                                                     sequences.shuffled <- do.call(universalmotif::shuffle_sequences, c(list(inputSeq), ep))
                                                     seqs$sequence <- as.vector(sequences.shuffled)
                                                     backgroudPeakseq <- seqs
                                                 }
                                                 ## get frequency from permutation
                                                 sapply(seq.int(pattern), function(i) {
                                                     oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(
                                                         backgroudPeakseq$sequence), width = nchar(pattern[i])))
                                                     posPlus <- gregexpr(pattern.t[i], names(oligoNVec), perl = TRUE)
                                                     expFre <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
                                                     if (revcomp == TRUE) {
                                                         posMinus <- gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
                                                         expFreMinus <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
                                                         expFre <- expFre + expFreMinus
                                                     }
                                                     c(expFre / sum(oligoNVec))
                                                 })
                                             }, ep = ep))
    colnames(permutationFreq) <- names(pattern)
    permutationFreqSorted <- apply(permutationFreq, 2, sort, decreasing = TRUE)
    peakFreq[, 4] <- as.vector(permutationFreqSorted[round(nperm * alpha, 0), ])
    return(peakFreq)
}

getPosInSeqs <- function(thispattern, seq, revcomp = TRUE) {
    if (missing(thispattern) || !inherits(thispattern, "DNAStringSet")) {
        stop("thispattern is required as a DNAStringSet object!", call. = FALSE)
    }
    if (missing(seq)) {
        stop("No valid sequences passed in!")
    }
    patternName <- names(thispattern)
    thispattern.t <- as.character(thispattern)[[1]]
    thispattern.t <- translatePattern(thispattern.t)
    if (revcomp == TRUE) {
        revcomp.pattern <- reverseComplement(thispattern)
        revcomp.pattern.t <- as.character(revcomp.pattern)[[1]]
        revcomp.pattern.t <- translatePattern(revcomp.pattern.t)
    }
    sequences <- seq$sequence
    #seqInfo <- as.data.frame(seqGRanges)
    seqInfo <- as.data.frame(seq)[, 1:5]
    
    total <- do.call(rbind, lapply(seq_along(sequences), function(i) {
        pos.plus <- gregexpr(thispattern.t, sequences[i], perl = TRUE)[[1]]
        if (pos.plus[1] > 0) {
            pattern.start <- as.numeric(pos.plus)
            pattern.end <- pattern.start + attr(pos.plus, "match.length") - 1
            patternSeq.found.plus <- 
                do.call(rbind, lapply(seq_along(pattern.start), function(j) {
                    c(patternName, as.character(thispattern), 
                      pattern.start[j], pattern.end[j], 
                      substr(sequences[i], pattern.start[j], 
                             pattern.end[j]), "+", seqInfo[i, ])
                }))
        }
        if (revcomp == TRUE) {
            pos.minus <- gregexpr(revcomp.pattern.t, sequences[i], perl = TRUE)[[1]]
            if (pos.minus[1] > 0) {
                pattern.start <- as.numeric(pos.minus)
                pattern.end <- pattern.start + attr(pos.minus, "match.length") - 1
                patternSeq.found.minus <- 
                    do.call(rbind, lapply(seq_along(pattern.start), function(j) {
                        c(patternName, as.character(thispattern), 
                          pattern.start[j], pattern.end[j], 
                          substr(sequences[i], pattern.start[j],
                                 pattern.end[j]), "-", seqInfo[i, ])
                    }))
                if (pos.plus[1] > 0) {
                    rbind(patternSeq.found.plus, patternSeq.found.minus)
                } else {
                    patternSeq.found.minus
                }
            } else {
                if (pos.plus[1] > 0) {
                    patternSeq.found.plus
                }
            }
        } else {
            if (pos.plus[1] > 0) {
                patternSeq.found.plus
            }
        }
    }))
    if (length(total) == 0) {
        cat(c(patternName, as.character(thispattern), 
              "not found in the input sequences!\n\n"))
    } else {
        total <- as.data.frame(total, stringsAsFactors = FALSE)
        colnames(total)[1:6] <- c("motifName", "motifPattern", "motifStartInPeak",
                                  "motifEndInPeak", "motifFound", 
                                  "motifFoundStrand")
        motifstartInChromosome <- as.numeric(total$start) + as.numeric(total[, 3]) - 1
        motifendInChromosome <- as.numeric(total$start) + as.numeric(total[, 4]) - 1
        total <- cbind(motifChr = motifstartInChromosome,
                       motifStartInChr = motifstartInChromosome, 
                       motifEndInChr = motifendInChromosome, 
                       total)
        colnames(total)[(ncol(total) - 4):ncol(total)] <- c("peakChr",
                                                          "peakStart",
                                                          "peakEnd",
                                                          "peakWidth",
                                                          "peakStrand")
        for (i in seq_len(ncol(total))) {
            total[, i] <- unlist(total[, i])
        }
        total$motifChr <- total$peakChr
    }
    total
}


