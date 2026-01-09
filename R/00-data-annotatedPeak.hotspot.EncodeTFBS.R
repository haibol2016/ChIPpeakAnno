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


#' High Occupancy Transcription Factor Binding Regions (HOT Spots)
#' 
#' A dataset containing genomic regions with high occupancy of transcription
#' factor binding sites from the human genome (hg19/GRCh37). HOT (High
#' Occupancy Target) spots are genomic regions that are bound by many
#' transcription factors in ChIP-seq experiments, representing regions with
#' high regulatory potential.
#' 
#' @details
#' **What are HOT Spots?**
#' 
#' HOT spots are genomic regions that show high occupancy by multiple
#' transcription-related factors. These regions are identified by analyzing
#' binding sites from more than 100 transcription factors across multiple
#' ENCODE ChIP-seq experiments. HOT spots are important because they:
#' 
#' - Represent regions with high regulatory potential
#' - Are often associated with active promoters and enhancers
#' - Can cause false positives in permutation tests if not excluded
#' - Should be removed from peak pools in statistical testing to prevent
#'   overestimation of associations
#' 
#' **Data Structure:**
#' 
#' The dataset is a \code{GRangesList} object where each element contains
#' genomic regions (GRanges) for different HOT spot categories. The names of
#' the list elements correspond to different HOT spot classifications (e.g.,
#' "HOT_All", "HOT_intergenic").
#' 
#' **Usage in ChIPpeakAnno:**
#' 
#' HOT spots are primarily used in permutation testing with
#' \code{\link{peakPermTest}} to create unbiased peak pools. When creating a
#' custom peak pool using \code{\link{preparePool}}, HOT spots should be
#' removed to prevent overestimation of the association between peak sets.
#' 
#' **How the data was generated:**
#' 
#' The data was generated by downloading HOT spot regions from the Gerstein
#' Lab ENCODE MetaTracks database and converting BED files to GRanges objects.
#' The original data processing code:
#' 
#' \preformatted{
#' temp <- tempfile()
#' url <- "http://metatracks.encodenets.gersteinlab.org"
#' 
#' # Download HOT spot data
#' download.file(file.path(url, "HOT_All_merged.tar.gz"), temp)
#' temp2 <- tempfile()
#' download.file(file.path(url, "HOT_intergenic_All_merged.tar.gz"), temp2)
#' 
#' # Extract and process
#' untar(temp, exdir = dirname(temp))
#' untar(temp2, exdir = dirname(temp))
#' 
#' # Convert BED files to GRanges
#' f <- dir(dirname(temp), "bed$", full.names = TRUE)
#' HOT.spots <- lapply(f, function(x) {
#'     toGRanges(x, format = "BED")
#' })
#' 
#' # Set names and create GRangesList
#' names(HOT.spots) <- gsub("_merged.bed", "", basename(f))
#' HOT.spots <- GRangesList(HOT.spots)
#' 
#' # Save
#' save(HOT.spots, file = "data/HOT.spots.rda", 
#'      compress = "xz", compression_level = 9)
#' }
#' 
#' @name HOT.spots
#' @docType data
#' @format A \code{\link[GenomicRanges]{GRangesList}} object containing HOT
#'   spot regions. Each element in the list is a \code{GRanges} object
#'   representing a different category of HOT spots (e.g., all HOT spots,
#'   intergenic HOT spots).
#' 
#' @references
#' Yip KY, Cheng C, Bhardwaj N, Brown JB, Leng J, Kundaje A, Rozowsky J,
#' Birney E, Bickel P, Snyder M, Gerstein M. Classification of human genomic
#' regions based on experimentally determined binding sites of more than 100
#' transcription-related factors. \emph{Genome Biology} 2012 Sep
#' 26;13(9):R48. doi: 10.1186/gb-2012-13-9-r48. PubMed PMID: 22950945; PubMed
#' Central PMCID: PMC3491392.
#' 
#' @source
#' Gerstein Lab ENCODE MetaTracks Database:
#' \url{http://metatracks.encodenets.gersteinlab.org/}
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' # Load the data
#' data(HOT.spots)
#' 
#' # Check the structure
#' class(HOT.spots)
#' length(HOT.spots)
#' names(HOT.spots)
#' 
#' # Count regions in each category
#' elementNROWS(HOT.spots)
#' 
#' # Access a specific category
#' if (length(HOT.spots) > 0) {
#'     hot_all <- HOT.spots[[1]]
#'     head(hot_all)
#'     length(hot_all)
#' }
#' 
#' # Example: Remove HOT spots from a peak pool for permutation testing
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' # Create a peak pool excluding HOT spots
#' pool <- preparePool(
#'     TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'     template = myPeaks
#' )
#' 
#' # Remove HOT spots from the pool
#' if ("HOT_All" %in% names(HOT.spots)) {
#'     hot_regions <- unlist(HOT.spots[["HOT_All"]])
#'     pool$grs <- pool$grs[!overlapsAny(pool$grs, hot_regions)]
#' }
#' 
#' # Use the filtered pool in permutation test
#' perm_test <- peakPermTest(peaks1, peaks2, pool = pool)
#' }
#' 
"HOT.spots"

#' ENCODE transcription factor binding site clusters (V3) for human (hg19)
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing transcription factor
#' binding site clusters from ENCODE Project (version 3) for human genome
#' assembly hg19. This dataset has been processed to remove HOT (High Occupancy
#' Target) spots, making it suitable as a background pool for permutation testing
#' and statistical analysis.
#' 
#' @name wgEncodeTfbsV3
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with approximately
#' 617,916 genomic regions representing transcription factor binding site
#' clusters. The object contains:
#' \describe{
#'   \item{seqnames}{Chromosome names (hg19 format: "chr1", "chr2", etc.)}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information (typically "*" for unstranded data)}
#' }
#' 
#' @details
#' This dataset is derived from the ENCODE Project's comprehensive collection
#' of transcription factor binding sites identified across multiple cell types
#' and conditions. The data has been:
#' \itemize{
#'   \item Clustered to merge overlapping binding sites from different
#'         transcription factors
#'   \item Filtered to remove HOT spots (genomic regions with unusually high
#'         binding across many factors, which may represent artifacts)
#'   \item Reduced to merge overlapping regions
#' }
#' 
#' \strong{Use cases:}
#' \itemize{
#'   \item Background pool for permutation testing with \code{\link{peakPermTest}}
#'   \item Reference set for assessing binding site enrichment
#'   \item Comparison dataset for ChIP-seq peak analysis
#' }
#' 
#' \strong{Data generation:}
#' The original data was downloaded from the ENCODE Project and processed as
#' follows (see examples section for code):
#' \enumerate{
#'   \item Downloaded wgEncodeRegTfbsClusteredV3.bed.gz from UCSC
#'   \item Converted to GRanges format
#'   \item Removed regions overlapping with HOT spots (using \code{\link{HOT.spots}})
#'   \item Reduced overlapping regions
#' }
#' 
#' @source
#' ENCODE Project Consortium. The ENCODE (ENCyclopedia Of DNA Elements) Project.
#' Data downloaded from:
#' \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz}
#' 
#' @references
#' ENCODE Project Consortium (2012) An integrated encyclopedia of DNA elements
#' in the human genome. \emph{Nature} 489:57-74.
#' \doi{10.1038/nature11247}
#' 
#' @seealso
#' \code{\link{HOT.spots}}, \code{\link{peakPermTest}}, \code{\link{preparePool}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(wgEncodeTfbsV3)
#' 
#' # Inspect the structure
#' wgEncodeTfbsV3
#' length(wgEncodeTfbsV3)
#' 
#' # View first few regions
#' head(wgEncodeTfbsV3)
#' 
#' # Check chromosome distribution
#' table(seqnames(wgEncodeTfbsV3))
#' 
#' # Use as background pool for permutation testing
#' \dontrun{
#' # Example: Compare your peaks to ENCODE binding sites
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' perm_test <- peakPermTest(
#'     peaks1 = your_peaks,
#'     peaks2 = wgEncodeTfbsV3,
#'     TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene
#' )
#' }
#' 
#' \dontrun{
#' # How to regenerate this dataset:
#' temp <- tempfile()
#' download.file(
#'     file.path(
#'         "http://hgdownload.cse.ucsc.edu", "goldenPath",
#'         "hg19", "encodeDCC",
#'         "wgEncodeRegTfbsClustered",
#'         "wgEncodeRegTfbsClusteredV3.bed.gz"
#'     ),
#'     temp
#' )
#' data <- read.delim(gzfile(temp, "r"), header = FALSE)
#' unlink(temp)
#' 
#' colnames(data)[1:4] <- c("seqnames", "start", "end", "TF")
#' wgEncodeRegTfbsClusteredV3 <- GRanges(
#'     seqnames = as.character(data$seqnames),
#'     ranges = IRanges(data$start, data$end),
#'     TF = data$TF
#' )
#' 
#' data(HOT.spots)
#' hot <- reduce(unlist(HOT.spots))
#' ol <- findOverlaps(wgEncodeRegTfbsClusteredV3, hot)
#' wgEncodeTfbsV3 <- wgEncodeRegTfbsClusteredV3[-unique(queryHits(ol))]
#' wgEncodeTfbsV3 <- reduce(wgEncodeTfbsV3)
#' 
#' save(
#'     list = "wgEncodeTfbsV3",
#'     file = "data/wgEncodeTfbsV3.rda",
#'     compress = "xz",
#'     compression_level = 9
#' )
#' }
"wgEncodeTfbsV3"
