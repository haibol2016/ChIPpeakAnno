#' Transcription Start Site (TSS) annotation for rat (Rnor_5.0/rn5)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Rattus norvegicus} based on the Rnor_5.0
#' (also known as rn5) genome assembly. This dataset was obtained from Ensembl
#' via \code{\link{biomaRt}} and can be used directly with \code{\link{annotatePeakInBatch}}
#' for peak annotation.
#' 
#' @name TSS.rat.Rnor_5.0
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
#' genome (Rnor_5.0/rn5 assembly). The TSS is defined as the 5' end of the
#' transcript for plus-strand genes and the 3' end for minus-strand genes.
#' 
#' \strong{Genome assembly:} Rnor_5.0 (also known as rn5)  
#' \strong{Source:} Ensembl via \code{\link{biomaRt}}  
#' \strong{Use case:} Direct annotation of ChIP-seq peaks to nearest TSS
#' 
#' \strong{Data generation:}
#' The dataset was obtained from Ensembl using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ensembl",
#'     dataset = "rnorvegicus_gene_ensembl"
#' )
#' TSS.rat.Rnor_5.0 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @seealso
#' \code{\link{TSS.rat.RGSC3.4}}, \code{\link{annotatePeakInBatch}},
#' \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.rat.Rnor_5.0)
#' 
#' # Inspect the structure
#' TSS.rat.Rnor_5.0
#' length(TSS.rat.Rnor_5.0)
#' 
#' # View first few TSS
#' head(TSS.rat.Rnor_5.0)
#' 
#' # Use for peak annotation
#' \dontrun{
#' library(ChIPpeakAnno)
#' annotated <- annotatePeakInBatch(
#'     peaks,
#'     AnnotationData = TSS.rat.Rnor_5.0,
#'     output = "nearestLocation",
#'     FeatureLocForDistance = "TSS"
#' )
#' }
"TSS.rat.Rnor_5.0"
