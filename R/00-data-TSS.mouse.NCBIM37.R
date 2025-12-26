#' Transcription Start Site (TSS) annotation for mouse (NCBIM37/mm9)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Mus musculus} based on the NCBIM37 (also
#' known as mm9) genome assembly. This is a legacy genome assembly and is
#' provided for compatibility with older datasets.
#' 
#' @name TSS.mouse.NCBIM37
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
#' This dataset contains TSS coordinates for all annotated genes in the mouse
#' genome (NCBIM37/mm9 assembly). The TSS is defined as the 5' end of the
#' transcript for plus-strand genes and the 3' end for minus-strand genes.
#' 
#' \strong{Genome assembly:} NCBIM37 (also known as mm9) - **Legacy assembly**  
#' \strong{Source:} Ensembl via biomaRt  
#' \strong{Use case:} Annotation of legacy datasets aligned to mm9
#' 
#' \strong{Data generation:}
#' The dataset was obtained using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ensembl",
#'     dataset = "mmusculus_gene_ensembl"
#' )
#' TSS.mouse.NCBIM37 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @note
#' NCBIM37 (mm9) is a legacy genome assembly that has been superseded by
#' GRCm38 (mm10). For new analyses, use \code{\link{TSS.mouse.GRCm38}}. This
#' dataset is maintained for compatibility with older datasets that were aligned
#' to mm9.
#' 
#' @seealso
#' \code{\link{TSS.mouse.GRCm38}}, \code{\link{annotatePeakInBatch}},
#' \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.mouse.NCBIM37)
#' 
#' # Inspect the structure
#' TSS.mouse.NCBIM37
#' length(TSS.mouse.NCBIM37)
#' 
#' # View first few TSS
#' head(TSS.mouse.NCBIM37)
"TSS.mouse.NCBIM37"
