#' ENCODE transcription factor binding site clusters (V3) for human (hg19)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing transcription factor
#' binding site clusters from ENCODE Project (version 3) for human genome
#' assembly hg19. This dataset has been processed to remove HOT (High Occupancy
#' Target) spots, making it suitable as a background pool for permutation testing
#' and statistical analysis.
#' 
#' @name wgEncodeTfbsV3
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with approximately
#' 617,916 genomic regions representing transcription factor binding site
#' clusters. The object contains:
#' \describe{
#'   \item{seqnames}{Chromosome names (hg19 format: "chr1", "chr2", etc.)}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information (typically "*" for unstranded data)}
#' }
#' 
#' @details
#' This dataset is derived from the ENCODE Project's comprehensive collection
#' of transcription factor binding sites identified across multiple cell types
#' and conditions. The data has been:
#' \itemize{
#'   \item Clustered to merge overlapping binding sites from different
#'         transcription factors
#'   \item Filtered to remove HOT spots (genomic regions with unusually high
#'         binding across many factors, which may represent artifacts)
#'   \item Reduced to merge overlapping regions
#' }
#' 
#' \strong{Use cases:}
#' \itemize{
#'   \item Background pool for permutation testing with \code{\link{peakPermTest}}
#'   \item Reference set for assessing binding site enrichment
#'   \item Comparison dataset for ChIP-seq peak analysis
#' }
#' 
#' \strong{Data generation:}
#' The original data was downloaded from the ENCODE Project and processed as
#' follows (see examples section for code):
#' \enumerate{
#'   \item Downloaded wgEncodeRegTfbsClusteredV3.bed.gz from UCSC
#'   \item Converted to GRanges format
#'   \item Removed regions overlapping with HOT spots (using \code{\link{HOT.spots}})
#'   \item Reduced overlapping regions
#' }
#' 
#' @source
#' ENCODE Project Consortium. The ENCODE (ENCyclopedia Of DNA Elements) Project.
#' Data downloaded from:
#' \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz}
#' 
#' @references
#' ENCODE Project Consortium (2012) An integrated encyclopedia of DNA elements
#' in the human genome. \emph{Nature} 489:57-74.
#' \doi{10.1038/nature11247}
#' 
#' @seealso
#' \code{\link{HOT.spots}}, \code{\link{peakPermTest}}, \code{\link{preparePool}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(wgEncodeTfbsV3)
#' 
#' # Inspect the structure
#' wgEncodeTfbsV3
#' length(wgEncodeTfbsV3)
#' 
#' # View first few regions
#' head(wgEncodeTfbsV3)
#' 
#' # Check chromosome distribution
#' table(seqnames(wgEncodeTfbsV3))
#' 
#' # Use as background pool for permutation testing
#' \dontrun{
#' # Example: Compare your peaks to ENCODE binding sites
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' perm_test <- peakPermTest(
#'     peaks1 = your_peaks,
#'     peaks2 = wgEncodeTfbsV3,
#'     TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene
#' )
#' }
#' 
#' \dontrun{
#' # How to regenerate this dataset:
#' temp <- tempfile()
#' download.file(
#'     file.path(
#'         "http://hgdownload.cse.ucsc.edu", "goldenPath",
#'         "hg19", "encodeDCC",
#'         "wgEncodeRegTfbsClustered",
#'         "wgEncodeRegTfbsClusteredV3.bed.gz"
#'     ),
#'     temp
#' )
#' data <- read.delim(gzfile(temp, "r"), header = FALSE)
#' unlink(temp)
#' 
#' colnames(data)[1:4] <- c("seqnames", "start", "end", "TF")
#' wgEncodeRegTfbsClusteredV3 <- GRanges(
#'     seqnames = as.character(data$seqnames),
#'     ranges = IRanges(data$start, data$end),
#'     TF = data$TF
#' )
#' 
#' data(HOT.spots)
#' hot <- reduce(unlist(HOT.spots))
#' ol <- findOverlaps(wgEncodeRegTfbsClusteredV3, hot)
#' wgEncodeTfbsV3 <- wgEncodeRegTfbsClusteredV3[-unique(queryHits(ol))]
#' wgEncodeTfbsV3 <- reduce(wgEncodeTfbsV3)
#' 
#' save(
#'     list = "wgEncodeTfbsV3",
#'     file = "data/wgEncodeTfbsV3.rda",
#'     compress = "xz",
#'     compression_level = 9
#' )
#' }
"wgEncodeTfbsV3"
