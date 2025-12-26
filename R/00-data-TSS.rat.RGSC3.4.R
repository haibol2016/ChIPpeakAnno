#' Transcription Start Site (TSS) annotation for rat (RGSC3.4/rn4)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Rattus norvegicus} based on the RGSC3.4
#' (also known as rn4) genome assembly. This is a legacy genome assembly and
#' is provided for compatibility with older datasets.
#' 
#' @name TSS.rat.RGSC3.4
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
#' This dataset contains TSS coordinates for all annotated genes in the rat
#' genome (RGSC3.4/rn4 assembly). The TSS is defined as the 5' end of the
#' transcript for plus-strand genes and the 3' end for minus-strand genes.
#' 
#' \strong{Genome assembly:} RGSC3.4 (also known as rn4) - **Legacy assembly**  
#' \strong{Source:} Ensembl via biomaRt  
#' \strong{Use case:} Annotation of legacy datasets aligned to rn4
#' 
#' \strong{Data generation:}
#' The dataset was obtained using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ensembl",
#'     dataset = "rnorvegicus_gene_ensembl"
#' )
#' TSS.rat.RGSC3.4 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @note
#' RGSC3.4 (rn4) is a legacy genome assembly that has been superseded by
#' Rnor_5.0 (rn5) and Rnor_6.0 (rn6). For new analyses, use
#' \code{\link{TSS.rat.Rnor_5.0}} or consider using the latest assembly. This
#' dataset is maintained for compatibility with older datasets that were aligned
#' to rn4.
#' 
#' @seealso
#' \code{\link{TSS.rat.Rnor_5.0}}, \code{\link{annotatePeakInBatch}},
#' \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.rat.RGSC3.4)
#' 
#' # Inspect the structure
#' TSS.rat.RGSC3.4
#' length(TSS.rat.RGSC3.4)
#' 
#' # View first few TSS
#' head(TSS.rat.RGSC3.4)
"TSS.rat.RGSC3.4"
