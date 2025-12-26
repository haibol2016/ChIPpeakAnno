#' Transcription Start Site (TSS) annotation for human (NCBI36/hg18)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Homo sapiens} based on the NCBI36 (also
#' known as hg18) genome assembly. This is a legacy genome assembly and is
#' provided for compatibility with older datasets.
#' 
#' @name TSS.human.NCBI36
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names (e.g., "1", "2", "X", "Y", "MT")}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with TSS coordinates}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Ensembl gene IDs as character vector}
#'   \item{description}{Gene description from Ensembl (in metadata columns)}
#' }
#' 
#' @details
#' This dataset contains TSS coordinates for all annotated genes in the human
#' genome (NCBI36/hg18 assembly). The TSS is defined as the 5' end of the
#' transcript for plus-strand genes and the 3' end for minus-strand genes.
#' 
#' \strong{Genome assembly:} NCBI36 (also known as hg18) - **Legacy assembly**  
#' \strong{Source:} Ensembl archive via biomaRt  
#' \strong{Use case:} Annotation of legacy datasets aligned to hg18
#' 
#' \strong{Data generation:}
#' The dataset was obtained from Ensembl archive using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ensembl_mart_47",
#'     dataset = "hsapiens_gene_ensembl",
#'     archive = TRUE
#' )
#' TSS.human.NCBI36 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @note
#' NCBI36 (hg18) is a legacy genome assembly that has been superseded by
#' GRCh37 (hg19) and GRCh38 (hg38). For new analyses, use
#' \code{\link{TSS.human.GRCh38}}. This dataset is maintained for compatibility
#' with older datasets that were aligned to hg18.
#' 
#' @seealso
#' \code{\link{TSS.human.GRCh37}}, \code{\link{TSS.human.GRCh38}},
#' \code{\link{annotatePeakInBatch}}, \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.human.NCBI36)
#' 
#' # Inspect the structure
#' TSS.human.NCBI36
#' length(TSS.human.NCBI36)
#' 
#' # View first few TSS
#' head(TSS.human.NCBI36)
"TSS.human.NCBI36"
