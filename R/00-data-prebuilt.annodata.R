#' Transcription Start Site (TSS) annotation for human (GRCh37/hg19)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Homo sapiens} based on the GRCh37 (also
#' known as hg19) genome assembly. This dataset was obtained from Ensembl via
#' biomaRt and can be used directly with \code{\link{annotatePeakInBatch}} for
#' peak annotation.
#' 
#' @name TSS.human.GRCh37
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
#' \strong{Important:} This dataset is provided for package examples and unit
#' testing only. Users should generate their own TSS annotations to match their
#' genome assembly using \code{\link{getAnnotation}} or EnsDb/TxDb packages.
#' 
#' This dataset contains TSS coordinates for all annotated genes in the human
#' genome (GRCh37/hg19 assembly). The TSS is defined as the 5' end of the
#' transcript for plus-strand genes and the 3' end for minus-strand genes.
#' 
#' \strong{Genome assembly:} GRCh37 (also known as hg19)  
#' \strong{Source:} Ensembl via biomaRt  
#' \strong{Intended use:} Package examples and unit testing only
#' 
#' \strong{Data generation:}
#' The dataset was obtained using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ENSEMBL_MART_ENSEMBL",
#'     host = "grch37.ensembl.org",
#'     path = "/biomart/martservice",
#'     dataset = "hsapiens_gene_ensembl"
#' )
#' TSS.human.GRCh37 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @note
#' \strong{For users:} Do not use this dataset for your own peak annotation.
#' Instead, generate annotations matching your genome assembly:
#' \preformatted{
#' # Recommended approach:
#' library(biomaRt)
#' mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#' TSS <- getAnnotation(mart, featureType = "TSS")
#' 
#' # Or use EnsDb packages:
#' library(EnsDb.Hsapiens.v86)
#' annoData <- annoGR(EnsDb.Hsapiens.v86)
#' }
#' 
#' GRCh37 (hg19) is a legacy genome assembly. For new analyses, consider using
#' GRCh38 (hg38) and generate annotations accordingly.
#' 
#' @seealso
#' \code{\link{annotatePeakInBatch}}, \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.human.GRCh37)
#' 
#' # Inspect the structure
#' TSS.human.GRCh37
#' length(TSS.human.GRCh37)
#' 
#' # View first few TSS
#' head(TSS.human.GRCh37)
#' 
#' # Use for peak annotation
#' \dontrun{
#' library(ChIPpeakAnno)
#' annotated <- annotatePeakInBatch(
#'     peaks,
#'     AnnotationData = TSS.human.GRCh37,
#'     output = "nearestLocation",
#'     FeatureLocForDistance = "TSS"
#' )
#' }
"TSS.human.GRCh37"


#' Transcription Start Site (TSS) annotation for mouse (GRCm38/mm10)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Mus musculus} based on the GRCm38 (also
#' known as mm10) genome assembly. This dataset was obtained from Ensembl via 
#' \code{\link{biomaRt}} and can be used directly with\code{\link{annotatePeakInBatch}} for
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
#' \strong{Important:} This dataset is provided for package examples and unit
#' testing only. Users should generate their own TSS annotations to match their
#' genome assembly using \code{\link{getAnnotation}} or EnsDb/TxDb packages.
#' 
#' This dataset contains TSS coordinates for all annotated genes in the mouse
#' genome (GRCm38/mm10 assembly). The TSS is defined as the 5' end of the
#' transcript for plus-strand genes and the 3' end for minus-strand genes.
#' 
#' \strong{Genome assembly:} GRCm38 (also known as mm10) 
#' \strong{Source:} Ensembl via \code{\link{biomaRt}}  
#' \strong{Intended use:} Package examples and unit testing only
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
#' \strong{For users:} Do not use this dataset for your own peak annotation.
#' Instead, generate annotations matching your genome assembly:
#' \preformatted{
#' # Recommended approach:
#' library(biomaRt)
#' mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#' TSS <- getAnnotation(mart, featureType = "TSS")
#' 
#' # Or use EnsDb packages:
#' library(EnsDb.Mmusculus.v79)
#' annoData <- annoGR(EnsDb.Mmusculus.v79)
#' }
#' 
#' 
#' @seealso
#' \code{\link{annotatePeakInBatch}},
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
