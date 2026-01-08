#' ChIPpeakAnno: Batch annotation and analysis of genomic peaks
#' 
#' @description
#' \code{ChIPpeakAnno} is a Bioconductor package for batch annotation and analysis
#' of peaks from ChIP-seq, ATAC-seq, CUT&RUN, CUT&Tag, CLIP-seq, and other
#' genomic interval experiments. It annotates peaks to genomic features (genes,
#' transcripts, exons, introns, UTRs, TSSs), performs enrichment analysis (GO, pathways),
#' extracts sequences for motif discovery, visualizes binding patterns, identifies peaks 
#' associated with bi-directional promoters, and tests statistical significance of peak 
#' overlaps. It also provides functionalities to identify peaks associated with genomic features (feature-centered), and perform
#' statistical testing on overlaps between peak sets.
#' 
#' @details
#' \strong{Core Capabilities:}
#' \itemize{
#'   \item \strong{Peak Annotation:} Associate peaks with nearest genes, exons,
#'         transcription start sites (TSS), or custom genomic features
#'   \item \strong{Multiple Annotation Sources:} Support for TxDb, EnsDb,
#'         biomaRt, and custom GRanges-based annotations
#'   \item \strong{Flexible Binding Types:} Annotate to TSS, gene ends, full
#'         gene ranges, or bi-directional promoters
#'   \item \strong{Enrichment Analysis:} Gene Ontology (GO) and pathway
#'         enrichment with multiple testing correction
#'   \item \strong{Sequence Analysis:} Retrieve sequences of peak regions for
#'         motif discovery and pattern matching
#'   \item \strong{Statistical Testing:} Hypergeometric tests, permutation
#'         tests, and overlap significance testing
#'   \item \strong{Visualization:} Venn diagrams, metagene plots, feature-aligned
#'         heatmaps, and signal profile plots
#'   \item \strong{Peak Comparison:} Find overlaps between multiple peak sets,
#'         assess reproducibility, and generate comparison statistics
#' }
#' 
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Bi-directional Promoter Detection:} Unique function
#'         (\code{peaksNearBDP}) to identify peaks near bi-directional promoters
#'   \item \strong{Enhancer Detection:} Use DNA interaction data (3C/HiC) to
#'         identify enhancers associated with peaks
#'   \item \strong{Signal Analysis:} Extract and visualize read coverage or
#'         signals around genomic features
#'   \item \strong{Pattern/Motif Enrichment:} Identify enriched DNA patterns
#'         with statistical testing
#'   \item \strong{Comprehensive Integration:} Works seamlessly with the
#'         Bioconductor ecosystem (TxDb, EnsDb, BSgenome, biomaRt, etc.)
#' }
#' 
#' \strong{Obtaining Annotation Data:}
#' \itemize{
#'   \item \strong{Recommended:} Use \code{\link{getAnnotation}} with 
#'         \code{\link{biomaRt}}, \code{\link{EnsDb}}, or \code{\link{TxDb}} 
#'         packages with \code{\link{annoGR}} to generate annotations
#' matching your genome assembly
#'   \item \strong{For examples only:} Pre-computed TSS datasets (e.g.,
#'         \code{TSS.human.GRCh37}, \code{TSS.mouse.GRCm38}) are provided for
#'         package examples and testing. Users should generate their own 
#'         annotations matching their genome assembly using the recommended
#'         methods.
#' }
#' 
#' \strong{Package Information:}
#' \tabular{ll}{
#'   Package: \tab ChIPpeakAnno\cr
#'   Type: \tab Package\cr
#'   Version: \tab See \code{DESCRIPTION} file\cr
#'   License: \tab GPL (>= 2)\cr
#'   LazyLoad: \tab yes\cr
#' }
#' 
#' @name ChIPpeakAnno-package
#' @aliases ChIPpeakAnno-package ChIPpeakAnno
#' @docType package
#' 
#' @author
#' Lihua Julie Zhu, Jianhong Ou, Jun Yu, Kai Hu, Haibo Liu, Junhui Li,
#' Hervé Pagès, Claude Gazin, Nathan Lawson, Ryan Thompson, Simon Lin,
#' David Lapointe, Michael Green
#' 
#' @references
#' \strong{Main Package Reference:}
#' \itemize{
#'   \item Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to
#'         annotate ChIP-seq and ChIP-chip data. \emph{BMC Bioinformatics}
#'         11:237. \doi{10.1186/1471-2105-11-237}
#' }
#' 
#' \strong{Statistical Methods:}
#' \itemize{
#'   \item Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate:
#'         a practical and powerful approach to multiple testing. \emph{J. R.
#'         Statist. Soc. B} 57:289-300.
#'   \item Benjamini Y, Yekutieli D (2001) The control of the false discovery
#'         rate in multiple hypothesis testing under dependency. \emph{Annals of
#'         Statistics} 29:1165-1188.
#'   \item Dudoit S, Shaffer JP, Boldrick JC (2003) Multiple hypothesis testing
#'         in microarray experiments. \emph{Statistical Science} 18:71-103.
#'   \item Hochberg Y (1988) A sharper Bonferroni procedure for multiple tests
#'         of significance. \emph{Biometrika} 75:800-802.
#'   \item Holm S (1979) A simple sequentially rejective multiple test procedure.
#'         \emph{Scand. J. Statist.} 6:65-70.
#' }
#' 
#' \strong{Integrated Packages:}
#' \itemize{
#'   \item Durinck S et al. (2005) BioMart and Bioconductor: a powerful link
#'         between biological biomarts and microarray data analysis.
#'         \emph{Bioinformatics} 21:3439-3440. \doi{10.1093/bioinformatics/bti525}
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{annotatePeakInBatch}} - Main function for peak annotation
#'   \item \code{\link{annoPeaks}} - Alternative annotation function with flexible
#'         binding types
#'   \item \code{\link{getAnnotation}} - Generate annotations from biomaRt
#'   \item \code{\link{annoGR}} - Convert TxDb/EnsDb to annotation GRanges
#'   \item \code{\link{getEnrichedGO}} - GO enrichment analysis
#'   \item \code{\link{getEnrichedPATH}} - Pathway enrichment analysis
#'   \item \code{\link{findOverlapsOfPeaks}} - Compare multiple peak sets
#'   \item \code{\link{peaksNearBDP}} - Find peaks near bi-directional promoters
#' }
#' 
#' For comprehensive examples and tutorials, see the package vignette:
#' \code{vignette("ChIPpeakAnno", package = "ChIPpeakAnno")}
#' 
#' @keywords package
"_PACKAGE"
