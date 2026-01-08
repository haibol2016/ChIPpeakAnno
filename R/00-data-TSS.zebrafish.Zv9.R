#' Transcription Start Site (TSS) annotation for zebrafish (Zv9/danRer10)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Danio rerio} (zebrafish) based on the Zv9
#' (also known as danRer10) genome assembly. This dataset was obtained from
#' Ensembl archive via \code{\link{biomaRt}} and can be used directly with
#' \code{\link{annotatePeakInBatch}} for peak annotation.
#' 
#' @name TSS.zebrafish.Zv9
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names (e.g., "1", "2", "MT")}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with TSS coordinates}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Ensembl gene IDs as character vector}
#'   \item{description}{Gene description from Ensembl (in metadata columns)}
#' }
#' 
#' @details
#' This dataset contains TSS coordinates for all annotated genes in the
#' zebrafish genome (Zv9/danRer10 assembly). The TSS is defined as the 5' end
#' of the transcript for plus-strand genes and the 3' end for minus-strand
#' genes.
#' 
#' \strong{Genome assembly:} Zv9 (also known as danRer10)  
#' \strong{Source:} Ensembl archive via \code{\link{biomaRt}}  
#' \strong{Use case:} Direct annotation of ChIP-seq peaks to nearest TSS
#' 
#' \strong{Data generation:}
#' The dataset was obtained from Ensembl archive using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ENSEMBL_MART_ENSEMBL",
#'     host = "mar2015.archive.ensembl.org",
#'     path = "/biomart/martservice",
#'     dataset = "drerio_gene_ensembl"
#' )
#' TSS.zebrafish.Zv9 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @note
#' Zv9 (danRer10) is a commonly used zebrafish genome assembly. Ensure your
#' peak coordinates are aligned to the same genome assembly (Zv9/danRer10).
#' 
#' @seealso
#' \code{\link{TSS.zebrafish.Zv8}}, \code{\link{annotatePeakInBatch}},
#' \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.zebrafish.Zv9)
#' 
#' # Inspect the structure
#' TSS.zebrafish.Zv9
#' length(TSS.zebrafish.Zv9)
#' 
#' # View first few TSS
#' head(TSS.zebrafish.Zv9)
#' 
#' # Use for peak annotation
#' \dontrun{
#' library(ChIPpeakAnno)
#' annotated <- annotatePeakInBatch(
#'     peaks,
#'     AnnotationData = TSS.zebrafish.Zv9,
#'     output = "nearestLocation",
#'     FeatureLocForDistance = "TSS"
#' )
#' }
"TSS.zebrafish.Zv9"
