#' Gene model with exon, 5' UTR and 3' UTR information for the human genome (GRCh37/hg19)
#' 
#' A pre-computed annotation dataset containing exon, 5' UTR, and 3' UTR
#' information for *Homo sapiens* based on the GRCh37/hg19 genome assembly.
#' This dataset was obtained from Ensembl via biomaRt and can be used directly
#' for peak annotation without requiring an active biomaRt connection.
#' 
#' @name ExonPlusUtr.human.GRCh37
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following
#'   structure:
#'   \describe{
#'     \item{\code{seqnames}}{Chromosome names (e.g., "chr1", "chr2")}
#'     \item{\code{ranges}}{IRanges object with start and end positions of exons}
#'     \item{\code{strand}}{Strand information: "+" for positive strand, "-" for negative strand}
#'     \item{\code{names}}{Ensembl transcript IDs}
#'     \item{\code{description}}{Transcript description from Ensembl}
#'     \item{\code{ensembl_gene_id}}{Ensembl gene ID associated with the transcript}
#'     \item{\code{utr5start}}{Start position of the 5' UTR (untranslated region)}
#'     \item{\code{utr5end}}{End position of the 5' UTR}
#'     \item{\code{utr3start}}{Start position of the 3' UTR}
#'     \item{\code{utr3end}}{End position of the 3' UTR}
#'   }
#' 
#' @details
#' This dataset was generated using the following code:
#' 
#' \preformatted{
#' library(biomaRt)
#' mart <- useMart(
#'     biomart = "ENSEMBL_MART_ENSEMBL",
#'     host = "grch37.ensembl.org",
#'     path = "/biomart/martservice",
#'     dataset = "hsapiens_gene_ensembl"
#' )
#' ExonPlusUtr.human.GRCh37 <- getAnnotation(mart = mart, 
#'                                           featureType = "ExonPlusUtr")
#' }
#' 
#' **Note**: GRCh37 (also known as hg19) is an older genome assembly. For
#' current analyses, consider using GRCh38/hg38 annotations when available.
#' 
#' @source
#' Data obtained from Ensembl (GRCh37 release) via biomaRt:
#' \url{https://grch37.ensembl.org/}
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' # Load the dataset
#' data(ExonPlusUtr.human.GRCh37)
#' 
#' # Inspect the structure
#' ExonPlusUtr.human.GRCh37
#' 
#' # Check the number of exons
#' length(ExonPlusUtr.human.GRCh37)
#' 
#' # View metadata columns
#' mcols(ExonPlusUtr.human.GRCh37)
#' 
#' # Example: Find exons for a specific gene
#' # (Replace "ENSG00000139618" with an actual Ensembl gene ID)
#' # gene_exons <- ExonPlusUtr.human.GRCh37[
#' #     mcols(ExonPlusUtr.human.GRCh37)$ensembl_gene_id == "ENSG00000139618"
#' # ]
#' 
"ExonPlusUtr.human.GRCh37"

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
#' \code{\link{TSS.human.GRCh38}}, \code{\link{TSS.human.NCBI36}},
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


#' Transcription Start Site (TSS) annotation for human (GRCh38/hg38)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Homo sapiens} based on the GRCh38 (also
#' known as hg38) genome assembly. This dataset was obtained from Ensembl via 
#' \code{\link{biomaRt}} and can be used directly with\code{\link{annotatePeakInBatch}} for
#' peak annotation.
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
#' \strong{Genome assembly:} GRCh38 (also known as hg38) 
#' \strong{Source:} Ensembl via \code{\link{biomaRt}}  
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
#' \strong{Source:} Ensembl archive via \code{\link{biomaRt}}  
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


#' Transcription Start Site (TSS) annotation for mouse (NCBIM37/mm9)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Mus musculus} based on the NCBIM37 (also
#' known as mm9) genome assembly. This dataset was obtained from Ensembl via 
#' \code{\link{biomaRt}} and can be used directly with \code{\link{annotatePeakInBatch}} for
#' peak annotation. It is provided for compatibility with older datasets.
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
#' \strong{Source:} Ensembl via \code{\link{biomaRt}}  
#' \strong{Use case:} Annotation of legacy datasets aligned to mm9
#' 
#' \strong{Data generation:}
#' The dataset was obtained from Ensembl using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ensembl",
#'     dataset = "mmusculus_gene_ensembl"
#' )
#' TSS.mouse.NCBIM37 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
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
#' \strong{Source:} Ensembl via \code{\link{biomaRt}}  
#' \strong{Use case:} Annotation of legacy datasets aligned to rn4
#' 
#' \strong{Data generation:}
#' The dataset was obtained from Ensembl using:
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


#' Transcription Start Site (TSS) annotation for zebrafish (Zv8/danRer7)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Transcription Start
#' Site (TSS) coordinates for \emph{Danio rerio} (zebrafish) based on the Zv8
#' (also known as danRer7) genome assembly. This is a legacy genome assembly
#' and is provided for compatibility with older datasets.
#' 
#' @name TSS.zebrafish.Zv8
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
#' zebrafish genome (Zv8/danRer7 assembly). The TSS is defined as the 5' end
#' of the transcript for plus-strand genes and the 3' end for minus-strand
#' genes.
#' 
#' \strong{Genome assembly:} Zv8 (also known as danRer7) - **Legacy assembly**  
#' \strong{Source:} Ensembl archive via \code{\link{biomaRt}}  
#' \strong{Use case:} Annotation of legacy datasets aligned to danRer7
#' 
#' \strong{Data generation:}
#' The dataset was obtained from Ensembl archive using:
#' \preformatted{
#' mart <- useMart(
#'     biomart = "ENSEMBL_MART_ENSEMBL",
#'     host = "may2009.archive.ensembl.org",
#'     path = "/biomart/martservice",
#'     dataset = "drerio_gene_ensembl"
#' )
#' TSS.zebrafish.Zv8 <- getAnnotation(mart, featureType = "TSS")
#' }
#' 
#' @note
#' Zv8 (danRer7) is a legacy genome assembly that has been superseded by Zv9
#' (danRer10) and later assemblies. For new analyses, use
#' \code{\link{TSS.zebrafish.Zv9}} or consider using the latest assembly. This
#' dataset is maintained for compatibility with older datasets that were aligned
#' to danRer7.
#' 
#' @seealso
#' \code{\link{TSS.zebrafish.Zv9}}, \code{\link{annotatePeakInBatch}},
#' \code{\link{getAnnotation}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the TSS annotation
#' data(TSS.zebrafish.Zv8)
#' 
#' # Inspect the structure
#' TSS.zebrafish.Zv8
#' length(TSS.zebrafish.Zv8)
#' 
#' # View first few TSS
#' head(TSS.zebrafish.Zv8)
"TSS.zebrafish.Zv8"


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
