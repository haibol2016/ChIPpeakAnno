#' Transcription Start Site (TSS) annotation for human (GRCh38/hg38)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Homo sapiens} based on the GRCh38 (also
#' known as hg38) genome assembly. This is the current reference genome
#' assembly for human and is recommended for new analyses. The dataset was
#' obtained from Ensembl via biomaRt and can be used directly with
#' \code{\link{annotatePeakInBatch}} for peak annotation.
#' 
#' @name TSS.human.GRCh38
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names (e.g., "1", "2", "X", "Y", "MT")}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with TSS coordinates}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Ensembl gene IDs as character vector}
#'   \item{metadata}{Additional metadata columns from Ensembl}
#' }
#' 
#' @details
#' This dataset contains TSS coordinates for all annotated genes in the human
#' genome (GRCh38/hg38 assembly). The TSS is defined as the 5' end of the
#' transcript for plus-strand genes and the 3' end for minus-strand genes.
#' 
#' \strong{Genome assembly:} GRCh38 (also known as hg38) - **Current reference**  
#' \strong{Source:} Ensembl via biomaRt  
#' \strong{Use case:} Direct annotation of ChIP-seq peaks to nearest TSS
#' 
#' \strong{Data generation:}
#' The dataset was obtained using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ensembl",
#'     dataset = "hsapiens_gene_ensembl"
#' )
#' TSS.human.GRCh38 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @note
#' GRCh38 (hg38) is the current reference genome assembly for human. This
#' dataset is recommended for new analyses. Ensure your peak coordinates are
#' aligned to the same genome assembly (GRCh38/hg38).
#' 
#' @seealso
#' \code{\link{TSS.human.GRCh37}}, \code{\link{TSS.human.NCBI36}},
#' \code{\link{annotatePeakInBatch}}, \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.human.GRCh38)
#' 
#' # Inspect the structure
#' TSS.human.GRCh38
#' length(TSS.human.GRCh38)
#' 
#' # View first few TSS
#' head(TSS.human.GRCh38)
#' 
#' # Use for peak annotation
#' \dontrun{
#' library(ChIPpeakAnno)
#' annotated <- annotatePeakInBatch(
#'     peaks,
#'     AnnotationData = TSS.human.GRCh38,
#'     output = "nearestLocation",
#'     FeatureLocForDistance = "TSS"
#' )
#' }
"TSS.human.GRCh38"
