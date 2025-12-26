#' Transcription Start Site (TSS) annotation for mouse (GRCm38/mm10)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Mus musculus} based on the GRCm38 (also
#' known as mm10) genome assembly. This is the current reference genome
#' assembly for mouse and is recommended for new analyses. The dataset was
#' obtained from Ensembl via biomaRt and can be used directly with
#' \code{\link{annotatePeakInBatch}} for peak annotation.
#' 
#' @name TSS.mouse.GRCm38
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
#' genome (GRCm38/mm10 assembly). The TSS is defined as the 5' end of the
#' transcript for plus-strand genes and the 3' end for minus-strand genes.
#' 
#' \strong{Genome assembly:} GRCm38 (also known as mm10) - **Current reference**  
#' \strong{Source:} Ensembl via biomaRt  
#' \strong{Use case:} Direct annotation of ChIP-seq peaks to nearest TSS
#' 
#' \strong{Data generation:}
#' The dataset was obtained using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ensembl",
#'     dataset = "mmusculus_gene_ensembl"
#' )
#' TSS.mouse.GRCm38 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @note
#' GRCm38 (mm10) is the current reference genome assembly for mouse. Ensure
#' your peak coordinates are aligned to the same genome assembly (GRCm38/mm10).
#' 
#' @seealso
#' \code{\link{TSS.mouse.NCBIM37}}, \code{\link{annotatePeakInBatch}},
#' \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.mouse.GRCm38)
#' 
#' # Inspect the structure
#' TSS.mouse.GRCm38
#' length(TSS.mouse.GRCm38)
#' 
#' # View first few TSS
#' head(TSS.mouse.GRCm38)
#' 
#' # Use for peak annotation
#' \dontrun{
#' library(ChIPpeakAnno)
#' annotated <- annotatePeakInBatch(
#'     peaks,
#'     AnnotationData = TSS.mouse.GRCm38,
#'     output = "nearestLocation",
#'     FeatureLocForDistance = "TSS"
#' )
#' }
"TSS.mouse.GRCm38"
