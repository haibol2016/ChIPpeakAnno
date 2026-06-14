#' Example Annotated Peaks Dataset
#' 
#' A \code{GRanges} object containing ChIP-seq peaks that have been annotated
#' with genomic features (transcription start sites). This dataset serves as
#' an example of the output format from \code{\link{annotatePeakInBatch}}.
#' 
#' @details
#' This dataset contains putative STAT1-binding regions identified in
#' un-stimulated cells using ChIP-seq technology (Robertson et al., 2007).
#' The peaks were annotated with transcription start sites (TSS) using the
#' \code{annotatePeakInBatch} function.
#' 
#' **How this dataset was created:**
#' \preformatted{
#' data(TSS.human.GRCh37)
#' data(myPeakList)
#' annotatedPeak <- annotatePeakInBatch(
#'     myPeakList, 
#'     AnnotationData = TSS.human.GRCh37,
#'     output = "both",
#'     multiple = FALSE)
#' }
#' 
#' @name annotatedPeak
#' @docType data
#' 
#' @format
#' A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' 
#' **Standard GRanges slots:**
#' \describe{
#'   \item{\code{seqnames}}{Chromosome/contig names (e.g., "chr1", "chr2")}
#'   \item{\code{ranges}}{IRanges object containing start and end positions}
#'   \item{\code{strand}}{Strand information ("+", "-", or "*")}
#'   \item{\code{names}}{Peak identifiers}
#' }
#' 
#' **Metadata columns (mcols):**
#' \describe{
#'   \item{\code{feature}}{Feature identifier (e.g., Ensembl gene ID, 
#'   transcript ID) that the peak is associated with}
#'   \item{\code{insideFeature}}{Spatial relationship between peak and feature:
#'   \itemize{
#'     \item \code{"upstream"}: Peak resides upstream of the feature
#'     \item \code{"downstream"}: Peak resides downstream of the feature
#'     \item \code{"inside"}: Peak is completely inside the feature
#'     \item \code{"overlapStart"}: Peak overlaps with the start of the feature
#'     \item \code{"overlapEnd"}: Peak overlaps with the end of the feature
#'     \item \code{"includeFeature"}: Peak completely includes/contains the feature
#'   }}
#'   \item{\code{distancetoFeature}}{Distance (in base pairs) from the peak's
#'   reference point to the feature's reference point (typically TSS). 
#'   Negative values indicate upstream, positive values indicate downstream.
#'   \code{NA} if the peak overlaps the feature.}
#'   \item{\code{start_position}}{Start position of the associated feature 
#'   (e.g., gene start, TSS position)}
#'   \item{\code{end_position}}{End position of the associated feature 
#'   (e.g., gene end, TES position)}
#' }
#' 
#' @source
#' Original ChIP-seq data from Robertson et al. (2007). Genome-wide profiles
#' of STAT1 DNA association using chromatin immunoprecipitation and massively
#' parallel sequencing. Nature Methods, 4(8), 651-657.
#' 
#' @references
#' Robertson, G., Hirst, M., Bainbridge, M., Bilenky, M., Zhao, Y., 
#' Zeng, T., ... & Jones, S. (2007). Genome-wide profiles of STAT1 DNA 
#' association using chromatin immunoprecipitation and massively parallel 
#' sequencing. \emph{Nature Methods}, 4(8), 651-657.
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' # Load the dataset
#' data(annotatedPeak)
#' 
#' # Inspect the structure
#' annotatedPeak
#' 
#' # View first few peaks
#' head(annotatedPeak, 4)
#' 
#' # Access metadata columns
#' mcols(annotatedPeak)
#' 
#' # Check distribution of peak-feature relationships
#' table(annotatedPeak$insideFeature)
#' 
#' # Plot distance distribution (interactive example)
#' if (interactive()) {
#'     distances <- annotatedPeak$distancetoFeature
#'     distances <- distances[!is.na(distances)]
#'     distances <- as.numeric(as.character(distances))
#'     
#'     hist(distances,
#'          xlab = "Distance To Nearest TSS (bp)",
#'          main = "Distribution of Peak-to-TSS Distances",
#'          breaks = 100,
#'          col = "steelblue",
#'          border = "white")
#' }
#' 
"annotatedPeak"


#' Example Enriched Gene Ontology (GO) Terms
#' 
#' A dataset containing example results from Gene Ontology enrichment analysis
#' performed using \code{\link{getEnrichedGO}}. This data object demonstrates
#' the structure and format of GO enrichment results returned by the
#' \code{getEnrichedGO} function.
#' 
#' @name enrichedGO
#' @docType data
#' 
#' @format A named list with 3 data frames, one for each GO ontology branch:
#' \describe{
#'   \item{\code{bp}}{Data frame containing enriched Biological Process (BP)
#'   terms with 9 columns (see Details below)}
#'   \item{\code{mf}}{Data frame containing enriched Molecular Function (MF)
#'   terms with 9 columns (see Details below)}
#'   \item{\code{cc}}{Data frame containing enriched Cellular Component (CC)
#'   terms with 9 columns (see Details below)}
#' }
#' 
#' @details
#' Each data frame in the list contains the following 9 columns:
#' \describe{
#'   \item{\code{go.id}}{Character. GO term identifier (e.g., "GO:0000001")}
#'   \item{\code{go.term}}{Character. GO term name (e.g., "mitochondrion inheritance")}
#'   \item{\code{go.Definition}}{Character. Detailed description of the GO term}
#'   \item{\code{Ontology}}{Character. Ontology branch: "BP" (Biological Process),
#'   "MF" (Molecular Function), or "CC" (Cellular Component)}
#'   \item{\code{count.InDataset}}{Integer. Number of genes in the input dataset
#'   annotated with this GO term}
#'   \item{\code{count.InGenome}}{Integer. Total number of genes in the genome
#'   annotated with this GO term}
#'   \item{\code{pvalue}}{Numeric. P-value from the hypergeometric test for
#'   enrichment}
#'   \item{\code{totaltermInDataset}}{Integer. Total number of GO terms
#'   associated with genes in the input dataset}
#'   \item{\code{totaltermInGenome}}{Integer. Total number of GO terms in the
#'   genome background}
#' }
#' 
#' @source
#' Generated using \code{\link{getEnrichedGO}} on example peak annotation data.
#' This is a demonstration dataset showing the expected output format.
#' 
#' @seealso
#' \code{\link{getEnrichedGO}} for performing GO enrichment analysis,
#' \code{\link{enrichmentPlot}} for visualizing enrichment results
#' 
#' @author Lihua Julie Zhu
#' @keywords datasets
#' 
#' @examples
#' # Load the example data
#' data(enrichedGO)
#' 
#' # Check the structure
#' str(enrichedGO)
#' 
#' # View dimensions of each ontology branch
#' dim(enrichedGO$bp)   # Biological Process
#' dim(enrichedGO$mf)   # Molecular Function
#' dim(enrichedGO$cc)   # Cellular Component
#' 
#' # View first few enriched terms in each category
#' head(enrichedGO$bp, n = 5)
#' head(enrichedGO$mf, n = 5)
#' head(enrichedGO$cc, n = 5)
#' 
#' # Access specific columns
#' enrichedGO$bp$go.term[1:5]
#' enrichedGO$bp$pvalue[1:5]
#' 
#' # Visualize enrichment results
#' if (requireNamespace("ChIPpeakAnno", quietly = TRUE)) {
#'     enrichmentPlot(enrichedGO)
#' }
"enrichedGO"


