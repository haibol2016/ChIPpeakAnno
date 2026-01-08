#' Calculate oligonucleotide frequencies
#' 
#' @description 
#' Calculates oligonucleotide frequencies for DNA sequences using Markov chain
#' models. This function computes frequencies for oligonucleotides of length 1,
#' \code{MarkovOrder}, and \code{MarkovOrder + 1}, which are used to estimate
#' expected frequencies for higher-order oligonucleotides in z-score
#' calculations. This function is primarily used internally by
#' \code{\link{oligoSummary}} but can be called directly to pre-compute
#' frequencies for multiple analyses.
#' 
#' @param sequence The input sequences. Can be:
#'        \itemize{
#'          \item A \code{\link[Biostrings]{DNAStringSet}} object
#'          \item A \code{\link[Biostrings]{DNAString}} object
#'          \item A \code{\link[GenomicRanges]{GRanges}} object with a
#'                \code{sequence} metadata column (output of
#'                \code{\link{getAllPeakSequence}})
#'          \item A character vector of DNA sequences
#'        }
#' @param MarkovOrder An integer specifying the order of the Markov chain model.
#'        Must be > 0. Default is \code{3L}. The function calculates frequencies
#'        for oligonucleotides of length 1, \code{MarkovOrder}, and
#'        \code{MarkovOrder + 1} to enable probability calculations for
#'        higher-order models.
#' 
#' @return Returns a named numeric vector containing frequencies (probabilities)
#'        for all oligonucleotides of length 1, \code{MarkovOrder}, and
#'        \code{MarkovOrder + 1}. Frequencies are normalized by the total number
#'        of nucleotides in all sequences, so they sum to 1 for each length
#'        category. The vector is named with the oligonucleotide sequences
#'        (e.g., "a", "c", "g", "t", "aaa", "aac", ...).
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates input and converts to appropriate format (handles GRanges
#'         with sequence metadata, character vectors, etc.)
#'   \item Calculates oligonucleotide frequencies for three lengths:
#'         \itemize{
#'           \item Length 1: Single nucleotides (A, C, G, T)
#'           \item Length \code{MarkovOrder}: k-mers of Markov order
#'           \item Length \code{MarkovOrder + 1}: (k+1)-mers for transition
#'                 probabilities
#'         }
#'   \item Normalizes frequencies by total sequence length to get probabilities
#'   \item Returns a flattened vector with all frequencies
#' }
#' 
#' \strong{Markov chain model:}
#' The function uses a Markov chain model of order \code{MarkovOrder} to
#' estimate expected frequencies. This model assumes that the probability of a
#' nucleotide depends on the previous \code{MarkovOrder} nucleotides. The
#' frequencies calculated here are used to compute transition probabilities in
#' \code{oligoSummary}.
#' 
#' \strong{Use cases:}
#' \itemize{
#'   \item Pre-computing frequencies for multiple \code{oligoSummary} calls
#'   \item Analyzing sequence composition independently
#'   \item Validating sequence data quality
#' }
#' 
#' @note
#' \itemize{
#'   \item The \code{seqinr} package must be installed and available
#'   \item Frequencies are normalized by total sequence length (not per sequence)
#'   \item The function handles sequences of different lengths
#'   \item Case is not preserved (sequences are converted to lowercase
#'         internally if needed)
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{oligoSummary}} for z-score calculations using these
#'         frequencies
#'   \item \code{\link{getAllPeakSequence}} for extracting sequences from peaks
#'   \item \code{\link[Biostrings]{oligonucleotideFrequency}} for the underlying
#'         frequency calculation
#' }
#' 
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet
#' @examples
#' 
#' # Example 1: Basic usage with DNAString
#' library(seqinr)
#' library(Biostrings)
#' seq <- DNAString("AATTCGACGTACAGATGACTAGACT")
#' freqs <- oligoFrequency(seq, MarkovOrder = 3L)
#' freqs
#' 
#' # Example 2: Multiple sequences with DNAStringSet
#' seqs <- DNAStringSet(c("AATTCGACGTACAGATGACTAGACT",
#'                        "GCTAGCTAGCTAGCTAGCTAGCTAG"))
#' freqs <- oligoFrequency(seqs, MarkovOrder = 3L)
#' 
#' # Example 3: Character vector
#' seqs_char <- c("AATTCGACGTACAGATGACTAGACT",
#'                "GCTAGCTAGCTAGCTAGCTAGCTAG")
#' freqs <- oligoFrequency(seqs_char, MarkovOrder = 3L)
#' 
#' # Example 4: Pre-compute frequencies for oligoSummary
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' seq <- getAllPeakSequence(peaks, genome = Hsapiens)
#' freqs <- oligoFrequency(seq, MarkovOrder = 3L)
#' # Use pre-computed frequencies in multiple calls
#' result1 <- oligoSummary(seq, freqs = freqs)
#' result2 <- oligoSummary(seq, freqs = freqs, oligoLength = 8L)
#' }
#' 
oligoFrequency <- function(sequence, MarkovOrder = 3L) {
    stopifnot(is.numeric(MarkovOrder))
    stopifnot(MarkovOrder > 0)
    stopifnot("The seqinr package is required." =
                  requireNamespace("seqinr", quietly = TRUE))
    if (inherits(sequence, "GRanges")) {
        sequence <- sequence$sequence
    }
    if (is.character(sequence)) {
        sequence <- DNAStringSet(sequence)
    }
    if (!inherits(sequence, c("DNAStringSet", "DNAString"))) {
        stop("sequence must be an object of DNAStringSet or DNAString ",
             "or output of getAllPeakSequence", call. = FALSE)
    }
    MarkovOrder <- as.integer(MarkovOrder)
    total <- sum(nchar(sequence))
    freqs <- lapply(unique(c(1, MarkovOrder, MarkovOrder + 1)), function(m) {
        of <- oligonucleotideFrequency(sequence,
                                       width = m,
                                       as.prob = FALSE)
        if (length(dim(of)) == 2) of <- colSums(of)
        of / total
    })
    unlist(freqs, recursive = FALSE)
}




#' Calculate z-scores for oligonucleotides using Markov chain models
#' 
#' @description 
#' Calculates z-scores for all oligonucleotides of a given length based on
#' Markov chain models. Z-scores measure how significantly enriched or depleted
#' each oligonucleotide is compared to the expected frequency under a Markov
#' chain model. This is useful for identifying overrepresented or
#' underrepresented sequence motifs in ChIP-seq peaks, regulatory regions, or
#' other DNA sequences.
#' 
#' Optionally, the function can generate position weight matrices (PWMs) from
#' top-scoring oligonucleotides, which can be used for motif discovery and
#' visualization.
#' 
#' @param sequence Input DNA sequences to analyze. Can be:
#'        \itemize{
#'          \item A \code{\link[Biostrings]{DNAStringSet}} object
#'          \item A \code{\link[Biostrings]{DNAString}} object
#'          \item A \code{\link[GenomicRanges]{GRanges}} object with a
#'                \code{sequence} metadata column (output of
#'                \code{\link{getAllPeakSequence}})
#'          \item A character vector of DNA sequences
#'        }
#'        Sequences are converted to lowercase internally.
#' @param oligoLength An integer specifying the length of oligonucleotides to
#'        analyze. Must be between 4 and 12 (inclusive). Default is \code{6L}.
#'        Longer oligonucleotides provide more specific motifs but require more
#'        computation and may have lower statistical power due to the large
#'        number of possible sequences (4^oligoLength).
#' @param freqs A numeric vector of pre-computed oligonucleotide frequencies
#'        from \code{\link{oligoFrequency}}. If \code{NULL} (default), frequencies
#'        are calculated automatically using the specified \code{MarkovOrder}.
#'        Pre-computing frequencies can save time when analyzing the same
#'        sequences with different parameters.
#' @param MarkovOrder An integer specifying the order of the Markov chain model.
#'        Must be between 1 and 5 (inclusive), and must be less than
#'        \code{oligoLength - 2}. Default is \code{3L}. Higher orders capture
#'        more complex sequence dependencies but require more data.
#' @param quickMotif A logical value. If \code{TRUE}, generates position weight
#'        matrices (PWMs) from top-scoring oligonucleotides. The function:
#'        \itemize{
#'          \item Selects top oligonucleotides based on z-scores
#'          \item Groups similar oligonucleotides using hierarchical clustering
#'          \item Creates consensus motifs and PWMs for each group
#'        }
#'        Default is \code{FALSE} (only returns z-scores and counts).
#' @param revcomp A logical value. If \code{TRUE}, considers both forward and
#'        reverse complement strands when counting oligonucleotides. This is
#'        useful for double-stranded DNA where motifs can appear on either
#'        strand. When \code{TRUE}, counts from complementary oligonucleotides
#'        are merged. Default is \code{FALSE}.
#' @param maxsize An integer specifying the maximum number of unique sequences
#'        to process before switching to table-based counting for efficiency.
#'        Default is \code{100000}. For datasets with many duplicate sequences,
#'        using table-based counting can significantly improve performance.
#' 
#' @return Returns a list with the following elements:
#'        \itemize{
#'          \item \code{zscore}: A named numeric vector of z-scores for each
#'                oligonucleotide. Positive z-scores indicate enrichment,
#'                negative z-scores indicate depletion. Z-scores are calculated
#'                as: \code{(observed - expected) / std_dev}. Names are the
#'                oligonucleotide sequences.
#'          \item \code{counts}: A named numeric vector of observed counts for
#'                each oligonucleotide across all sequences.
#'          \item \code{expCnt}: A named numeric vector of expected counts
#'                under the Markov chain model (only returned if
#'                \code{quickMotif = FALSE}). Expected counts are calculated
#'                based on the Markov chain transition probabilities.
#'          \item \code{motifs}: If \code{quickMotif = TRUE}, a list of position
#'                weight matrices (PWMs) representing discovered motifs. Each PWM
#'                is a 4 x motif_length matrix with rows A, C, G, T and columns
#'                representing positions. If \code{quickMotif = FALSE}, this is
#'                \code{NA}.
#'        }
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Validates input and converts sequences to appropriate format
#'   \item Counts all oligonucleotides of length \code{oligoLength} in the
#'         sequences (using efficient pattern matching with \code{PDict})
#'   \item If \code{revcomp = TRUE}, merges counts from complementary
#'         oligonucleotides
#'   \item Calculates expected frequencies using Markov chain model:
#'         \itemize{
#'           \item Uses transition probabilities from \code{MarkovOrder}-mer to
#'                 (\code{MarkovOrder}+1)-mer
#'           \item Computes expected probability for each oligonucleotide
#'           \item Converts to expected counts based on sequence lengths
#'         }
#'   \item Calculates z-scores: \code{(observed - expected) / std_dev}
#'   \item If \code{quickMotif = TRUE}, generates motif matrices from
#'         top-scoring oligonucleotides
#' }
#' 
#' \strong{Markov chain model:}
#' The function uses a Markov chain of order \code{MarkovOrder} to estimate
#' expected frequencies. The probability of an oligonucleotide is calculated as:
#' \preformatted{
#' P(oligo) = P(n1|n2...nk) * P(n2|n3...nk+1) * ... * P(nL-k|nL-k+1...nL)
#' }
#' where k = \code{MarkovOrder} and L = \code{oligoLength}. This accounts for
#' sequence context dependencies.
#' 
#' \strong{Z-score calculation:}
#' Z-scores are calculated as:
#' \preformatted{
#' z = (observed - expected) / sqrt(expected * (1 - expected) / N)
#' }
#' where N is the total number of possible positions. High positive z-scores
#' indicate significant enrichment, high negative z-scores indicate
#' significant depletion.
#' 
#' \strong{Motif generation (quickMotif = TRUE):}
#' When motif generation is enabled:
#' \itemize{
#'   \item Top 50 oligonucleotides by z-score are selected
#'   \item Oligonucleotides within half of the maximum z-score are kept
#'   \item Z-score cutoff filtering is applied (varies by \code{oligoLength})
#'   \item Similar oligonucleotides are grouped using hierarchical clustering
#'         based on edit distance
#'   \item Consensus sequences are generated using pairwise alignment
#'   \item Position weight matrices (PWMs) are created for each group
#' }
#' 
#' \strong{Performance optimization:}
#' For large datasets:
#' \itemize{
#'   \item If number of sequences > \code{maxsize}, uses table-based counting
#'         (counts unique sequences first, then counts oligonucleotides)
#'   \item This reduces memory usage and computation time for datasets with many
#'         duplicate sequences
#'   \item The function will stop with an error if the number of unique
#'         sequences still exceeds \code{maxsize} after table-based counting
#' }
#' 
#' \strong{Use cases:}
#' \itemize{
#'   \item Identifying enriched sequence motifs in ChIP-seq peaks
#'   \item Finding transcription factor binding sites
#'   \item Analyzing sequence composition of regulatory regions
#'   \item Discovering overrepresented k-mers in genomic regions
#' }
#' 
#' @note
#' \itemize{
#'   \item The \code{seqinr} package must be installed and available
#'   \item Sequences are converted to lowercase internally
#'   \item Z-scores with zero standard deviation are set to \code{NA}
#'   \item For motif generation, oligonucleotides must pass z-score cutoffs
#'         (varies by length: 2.155 for length 4, 2.66 for length 5, etc.)
#'   \item The function uses efficient pattern matching but can be slow for
#'         very large datasets or long oligonucleotides
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{oligoFrequency}} for pre-computing frequencies
#'   \item \code{\link{getAllPeakSequence}} for extracting sequences from peaks
#'   \item \code{\link[Biostrings]{oligonucleotideFrequency}} for basic
#'         frequency counting
#' }
#' 
#' @references van Helden, Jacques, Marcel li del Olmo, and Jose E.
#' Perez-Ortin. "Statistical analysis of yeast genomic downstream sequences
#' reveals putative polyadenylation signals." Nucleic Acids Research 28.4
#' (2000): 1000-1010.
#' 
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom Biostrings DNAStringSet PDict vcountPDict DNAString 
#' @importFrom pwalign pairwiseAlignment aligned
#' @importFrom utils adist combn
#' @importFrom stats hclust kmeans as.dendrogram nobs
#' @examples
#' 
#' # Example 1: Basic usage with peak sequences
#' if(interactive() || Sys.getenv("USER")=="jou"){
#'     data(annotatedPeak)
#'     library(BSgenome.Hsapiens.UCSC.hg19)
#'     library(seqinr)
#'     seq <- getAllPeakSequence(annotatedPeak[1:100], 
#'                  upstream=20, 
#'                  downstream=20, 
#'                  genome=Hsapiens)
#'     result <- oligoSummary(seq)
#'     # View top enriched oligonucleotides
#'     head(sort(result$zscore, decreasing=TRUE))
#' }
#' 
#' # Example 2: Different oligonucleotide length
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(seqinr)
#' seq <- getAllPeakSequence(peaks, genome = Hsapiens)
#' # Analyze 8-mers instead of 6-mers
#' result <- oligoSummary(seq, oligoLength = 8L)
#' }
#' 
#' # Example 3: With reverse complement consideration
#' \dontrun{
#' # Count oligonucleotides on both strands
#' result <- oligoSummary(seq, revcomp = TRUE)
#' }
#' 
#' # Example 4: Generate motif matrices
#' \dontrun{
#' # Generate PWMs from top-scoring oligonucleotides
#' result <- oligoSummary(seq, quickMotif = TRUE)
#' # Access motif matrices
#' result$motifs
#' }
#' 
#' # Example 5: Pre-compute frequencies for multiple analyses
#' \dontrun{
#' library(seqinr)
#' # Pre-compute frequencies once
#' freqs <- oligoFrequency(seq, MarkovOrder = 3L)
#' # Use in multiple analyses
#' result6 <- oligoSummary(seq, oligoLength = 6L, freqs = freqs)
#' result8 <- oligoSummary(seq, oligoLength = 8L, freqs = freqs)
#' }
#' 
#' # Example 6: Character vector input
#' \dontrun{
#' library(seqinr)
#' sequences <- c("AATTCGACGTACAGATGACTAGACT",
#'                "GCTAGCTAGCTAGCTAGCTAGCTAG",
#'                "AATTCGACGTACAGATGACTAGACT")
#' result <- oligoSummary(sequences, oligoLength = 6L)
#' }
#' 
oligoSummary <- function(sequence,
                         oligoLength = 6L,
                         freqs = NULL,
                         MarkovOrder = 3L,
                         quickMotif = FALSE,
                         revcomp = FALSE,
                         maxsize = 100000) {
    oligoLength <- as.integer(oligoLength)
    stopifnot(oligoLength > 3 & oligoLength < 13)
    MarkovOrder <- as.integer(MarkovOrder)
    stopifnot(MarkovOrder > 0 & MarkovOrder < 6)
    stopifnot(MarkovOrder < oligoLength - 2)
    stopifnot("The seqinr package is required." =
                  requireNamespace("seqinr", quietly = TRUE))
    if (inherits(sequence, "GRanges")) {
        sequence <- sequence$sequence
    } else if (inherits(sequence, "DNAStringSet")) {
        sequence <- as.character(sequence)
    } else if (!is.character(sequence)) {
        stop("'sequence' must be a DNAStringSet, DNAString, or output of ",
             "getAllPeakSequence", call. = FALSE)
    }
    
    sequence <- tolower(sequence)
    oligoWords <- seqinr::words(oligoLength)
    dict <- PDict(DNAStringSet(oligoWords))
    len <- length(sequence)
    sequence.tbl <- 1
    if (len > maxsize) {
        sequence.tbl <- table(sequence)
        if (length(sequence.tbl) > maxsize) {
            stop("oligoSummary cannot handle such huge dataset. ",
                 "Consider increasing 'maxsize' or reducing input size",
                 call. = FALSE)
        }
        cnt <- vcountPDict(dict,
                           subject = DNAStringSet(names(sequence.tbl)),
                           max.mismatch = 0,
                           min.mismatch = 0,
                           with.indels = FALSE,
                           fixed = TRUE)
    } else {
        cnt <- vcountPDict(dict,
                           subject = DNAStringSet(sequence),
                           max.mismatch = 0,
                           min.mismatch = 0,
                           with.indels = FALSE,
                           fixed = TRUE)
    }
    rownames(cnt) <- oligoWords
    cnt <- t(cnt)
    mergeRevcomp <- function(mat) {
        coln <- colnames(mat)
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
        mat1 <- colSums(mat1)
        mat2 <- colSums(mat2)
        coln <- ifelse(mat1 > mat2, names(mat1), names(mat2))
        mat <- mat[, coln.ids[, 1]] + mat[, coln.ids[, 2]]
        colnames(mat) <- coln
        mat
    }
    if (revcomp) {
        cnt <- mergeRevcomp(cnt)
    }
    # cnt <- cnt > 0
    mode(cnt) <- "logical"
    mode(cnt) <- "numeric"
    if (len > maxsize) {
        cnt <- cnt * as.numeric(sequence.tbl)
    }
    cntSum <- colSums(cnt)
    ## calculate Z score
    ## shuffle sequence based on 1 level markovmodel
    if (is.null(freqs)) {
        freqs <- oligoFrequency(sequence, MarkovOrder = MarkovOrder)
    }
    namesFreqs <- unique(c("a", "c", "g", "t",
                           seqinr::words(MarkovOrder),
                           seqinr::words(MarkovOrder + 1)))
    names(freqs) <- tolower(names(freqs))
    stopifnot(all(namesFreqs %in% names(freqs)))
    f <- sapply(names(cntSum), function(.ele) {
        m1 <- substring(.ele,
                         1:(oligoLength - MarkovOrder),
                         (MarkovOrder + 1):oligoLength)
        m0 <- substring(.ele,
                         2:(oligoLength - MarkovOrder),
                         (MarkovOrder + 1):(oligoLength - 1))
        prod(freqs[m1]) / prod(freqs[m0])
    })
    if (len > maxsize) {
        seqLen <- vapply(names(sequence.tbl),
                        nchar,
                        FUN.VALUE = integer(1),
                        USE.NAMES = FALSE) - oligoLength + 1L
    } else {
        seqLen <- vapply(sequence,
                        nchar,
                        FUN.VALUE = integer(1),
                        USE.NAMES = FALSE) - oligoLength + 1L
    }
    seqLen[seqLen < 0L] <- 0L
    f.m <- seqLen %*% t(f)
    f.m[f.m > 1] <- 1
    if (len > maxsize) {
        f.m <- f.m * as.numeric(sequence.tbl)
    }
    mu <- colSums(f.m)
    names(mu) <- names(cntSum)
    #     Kov <- sapply(names(cntSum), function(.ele){
    #         arr <- seqinr::s2c(.ele)
    #         k <- lapply(seq_len(length(arr)), function(.e){
    #             if(.e==1) return(TRUE)
    #             if(.e %% 2)
    #                 return(all(arr[1:ceiling(.e/2)]==arr[ceiling(.e/2):.e]))
    #             return(all(arr[1:(.e/2)]==arr[(.e/2+1):.e]))
    #         })
    #         k <- as.numeric(k)
    #         j <- sapply(seq_len(length(arr)), function(.e){
    #             (1/freqs[arr[.e]])^.e
    #         })
    #         sum(k*j)
    #     })
    # std <- sqrt(mu * (2*Kov - 1 - (2*oligoLength - 1) * mu/length(sequence)))
    # std <- sqrt(length(sequence) *
    #             (Kov/2-1+(2*oligoLength-1)/4^oligoLength) / 4^oligoLength)
    N <- len  # length(sequence)
    std <- sqrt(mu / N * (1 - mu / N) / N)
    zscore <- (cntSum - mu) / N / std
    # zscore <- (cntSum - mu)/std
    zscore[std == 0] <- NA
    
    if (!quickMotif) {
        return(list(zscore = zscore,
                   counts = cntSum,
                   expCnt = mu,
                   motifs = NA))
    }
    ## motif search
    seeds <- zscore[!is.na(zscore)]
    seeds <- sort(seeds, decreasing = TRUE)
    seeds <- seeds[seq_len(min(50L, length(seeds)))]
    zscore.max <- max(zscore)
    seeds <- seeds[zscore.max - seeds <= zscore.max / 2]
    zscore.cutoff <- c(NA, NA, 2.155, 2.66, 3.095, 3.49, 3.83, 4.1, 5, 5, 5, 5)
    seeds <- names(seeds)[seeds > zscore.cutoff[oligoLength]]  # top50 only
    
    str2motif <- function(s) {
        tmp <- matrix(c(1, 0, 0, 0,
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, 1),
                     ncol = 4,
                     nrow = 4,
                     dimnames = list(c("A", "C", "G", "T"),
                                     c("a", "c", "g", "t")))
        if (length(s) == 1) {
            return(tmp[, seqinr::s2c(s)])
        }
        ss <- table(s)
        if (length(ss) == 1) {
            return(tmp[, seqinr::s2c(names(ss))])
        }
        sss <- lapply(names(ss), function(.ele)
            as.numeric(tmp[, seqinr::s2c(.ele)]))
        sss <- mapply(function(mat, times, coln) mat * times,
                      sss, ss, SIMPLIFY = FALSE)
        consensusStr <- DNAString(names(ss)[1])
        for (i in seq.int(2L, length(ss))) {
            tempStr <- DNAString(names(ss)[i])
            consensusStr <- aligned(
                pairwiseAlignment(consensusStr, tempStr,
                                 type = "local"),
                degap = TRUE
            )
        }
        consensusStr <- tolower(as.character(consensusStr))
        idx <- regexpr(consensusStr, names(ss))
        tl <- max(idx)
        nc <- nchar(names(ss)[1])
        tmp0 <- rep(0, 4 * (tl + nc - 1))
        pcm <- mapply(function(.ele, .idx) {
            left <- max(idx) - .idx
            this.tmp <- tmp0
            this.tmp[(left * 4 + 1L):length(.ele)] <- .ele
            this.tmp
        }, sss, idx, SIMPLIFY = FALSE)
        pcm <- do.call(rbind, pcm)
        pcm <- colSums(pcm)
        pcm <- matrix(pcm,
                     nrow = 4,
                     dimnames = list(c("A", "C", "G", "T")))
        pcm.colsums <- colSums(pcm)
        pcm.colsums.max <- max(pcm.colsums)
        pcm <- t(t(pcm) + (pcm.colsums.max - pcm.colsums) / 4)
        pcm / colSums(pcm)
    }
    subgroupMotif <- function(.ele) {
        .ele <- .ele[!is.na(.ele)]
        if (length(.ele) == 0) {
            return(NA)
        }
        if (length(.ele) == 1) {
            return(str2motif(.ele))
        }
        .cnt <- cnt[, .ele]
        .ord <- colSums(.cnt)
        .cnt <- .cnt[, order(-.ord)]
        if (length(.ele) == 2) {
            .ele <- colnames(.cnt)
            .ele <- c(rep(.ele[1], sum(.cnt[, 1])),
                      rep(.ele[2], sum(.cnt[.cnt[, 2] > 0 &
                                          (!(.cnt[, 1] > 0 &
                                                .cnt[, 2] > 0)), 2])))
            return(str2motif(.ele))
        }
        ## get max hits of any combinations
        ad <- adist(.ele, partial = TRUE)
        rownames(ad) <- colnames(ad) <- .ele
        all.comb <- sapply(seq_len(min(length(.ele), 5L)), function(m) {
            .comb <- combn(.ele, m, simplify = FALSE)
            .comb.ad <- sapply(.comb, function(.cmb) {
                if (length(.cmb) > 1) {
                    .cmb <- combn(.cmb, 2)
                    all(ad[.cmb[1, ], .cmb[2, ]] <= 2)
                } else {
                    TRUE
                }
            })
            .comb[.comb.ad]
        })
        all.comb <- unlist(all.comb, recursive = FALSE, use.names = FALSE)
        prob.mul.evt <- function(x) {
            if (length(x) == 1) return(x)
            p <- x[1] + x[2] - x[1] * x[2]
            Recall(c(p, x[-(1:2)]))
        }
        comb.f <- sapply(all.comb,
                         function(.e) prob.mul.evt(f[.e]),
                         simplify = TRUE,
                         USE.NAMES = FALSE)
        comb.f.m <- seqLen %*% t(comb.f)
        comb.f.m[comb.f.m > 1] <- 1
        if (len > maxsize) {
            comb.f.m <- comb.f.m * as.numeric(sequence.tbl)
        }
        comb.mu <- colSums(comb.f.m)
        rm(comb.f.m)
        gc(reset = TRUE)
        comb.std <- sqrt(comb.mu / N * (1 - comb.mu / N) / N)
        comb.score <- sapply(all.comb,
                            function(.e) sum(.cnt[, .e]),
                            simplify = TRUE,
                            USE.NAMES = FALSE)
        score <- (comb.score - comb.mu) / N / comb.std
        best.combs <- all.comb[score == max(score)]
        best.combs.score <- sapply(best.combs, function(.ele) sum(.ord[.ele]))
        best.comb <- best.combs[best.combs.score == max(best.combs.score)][[1]]
        .cnt <- .cnt[, colnames(.cnt) %in% best.comb, drop = FALSE]
        ## resort the column order
        .cnt <- .cnt[, order(.ord[colnames(.cnt)], decreasing = TRUE), drop = FALSE]
        .cnt.gt0 <- .cnt > 0
        .ele <- rep(colnames(.cnt)[1], sum(.cnt[, 1]))
        if (ncol(.cnt) > 1) {
            for (i in 2:ncol(.cnt)) {
                .cnt.gt0[, i] <- .cnt.gt0[, i] &
                    !(apply(.cnt.gt0[, 1:(i - 1), drop = FALSE], 1, any) &
                         .cnt.gt0[, i])
                .ele <- c(.ele, rep(colnames(.cnt)[i],
                                    sum(.cnt[.cnt.gt0[, i], i])))
            }
        }
        str2motif(.ele)
    }
    if (length(seeds) > 2) {
        ### split the seeds by distance
        dist <- adist(seeds, partial = TRUE)
        hc <- hclust(dist(dist), method = "average")
        len <- length(hc$height)
        if (len > 2) {
            hc.height.diff <- diff(hc$height)
            km <- kmeans(hc.height.diff, centers = 2)
            idx <- km$cluster == km$cluster[length(km$cluster)]
            idx <- rev(seq_along(idx))[which(rev(!idx))[1L]] + 1L
            d <- cut(as.dendrogram(hc), hc$height[idx])
            lowers <- d$lower[sapply(d$lower, nobs) >= 1]
            if (length(lowers) > 0) {
                subgroup <- lapply(lowers, function(.ele)
                    seeds[as.numeric(labels(.ele))])
                subgroup <- subgroup[order(sapply(subgroup, function(.ele)
                    max(zscore[.ele])), decreasing = TRUE)]
                motifs <- lapply(subgroup, subgroupMotif)
            } else {
                motifs <- NA  ## this is impossible
            }
        } else {
            d <- as.dendrogram(hc)
            subgroup <- list()
            subgroup[[1L]] <- seeds[as.numeric(labels(d[[1L]]))]
            subgroup[[2L]] <- seeds[as.numeric(labels(d[[2L]]))]
            motifs <- lapply(subgroup, subgroupMotif)
        }
    } else {
        if (length(seeds) == 2) {
            consensus <- tolower(as.character(aligned(
                pairwiseAlignment(DNAString(seeds[1]),
                                 DNAString(seeds[2]),
                                 type = "local"),
                degap = TRUE)))
            if (nchar(consensus) >= oligoLength - 2) {
                motifs <- list(subgroupMotif(seeds))
            } else {
                motifs <- lapply(seeds[order(-zscore[seeds])], str2motif)
            }
        } else {
            motifs <- lapply(seeds[order(-zscore[seeds])], str2motif)
        }
    }
    
    return(list(zscore = zscore,
               counts = cntSum,
               motifs = motifs))
}
