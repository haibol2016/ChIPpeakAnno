#' ChIPpeakAnno: Batch annotation and analysis of genomic peaks
#' 
#' @description
#' \code{ChIPpeakAnno} is a Bioconductor package for batch annotation and analysis
#' of peaks from ChIP-seq, ATAC-seq, CUT&RUN, CUT&Tag, CLIP-seq, and other
#' genomic interval experiments. It provides a comprehensive workflow for peak
#' annotation including quality assessment, consensus peak creation, factor-specific
#' hierarchical annotation, enhancer identification, and downstream enrichment analysis.
#' The package annotates peaks to genomic features (genes, transcripts, exons, introns,
#' UTRs, TSSs), integrates with epigenomic databases (ENCODE cCREs, Roadmap Epigenomics),
#' performs enrichment analysis (GO, pathways, MSigDB), extracts sequences for motif
#' discovery, visualizes binding patterns, identifies peaks associated with bi-directional
#' promoters and enhancers, and tests statistical significance of peak overlaps.
#' 
#' @details
#' \strong{Core Capabilities:}
#' \itemize{
#'   \item \strong{Quality Assessment:} Comprehensive peak quality evaluation
#'         per replicate including width distribution, score metrics, and visualization
#'         (\code{\link{assessPeaks}})
#'   \item \strong{Peak Filtering:} Filter peaks based on width, scores, and
#'         quality metrics (\code{\link{filterPeaks}})
#'   \item \strong{Consensus Peak Creation:} Combine replicates using IDR
#'         filtering or majority voting (\code{\link{IDRfilter}},
#'         \code{\link{findOverlapsOfPeaks}})
#'   \item \strong{Preliminary Annotation:} Quick overview of peak distribution
#'         across genomic features (promoters, UTRs, exons, introns, intergenic)
#'         (\code{\link{getGenomicAnnotation}})
#'   \item \strong{Hierarchical Annotation:} Factor-specific annotation with
#'         multiple parallel strategies and intelligent prioritization for TFs,
#'         histone marks, Pol II, and other factors (\code{\link{annotateHierarchically}})
#'   \item \strong{Multiple Annotation Sources:} Support for TxDb, EnsDb,
#'         biomaRt, and custom GRanges-based annotations
#'   \item \strong{Epigenomic Integration:} Annotate with ENCODE cCREs
#'         (\code{\link{annotatePeaksWithcCRE}}) and Roadmap Epigenomics data
#'         (\code{\link{annotatePeaksWithRoadmap}}) for enhancer and chromatin
#'         state annotations
#'   \item \strong{Enhancer Identification:} Link peaks to enhancers and target
#'         genes using Hi-C data or ENCODE CRE annotations (\code{\link{findEnhancers}})
#'   \item \strong{Flexible Binding Types:} Annotate to TSS, gene ends, full
#'         gene ranges, or bi-directional promoters
#'   \item \strong{Enrichment Analysis:} Gene Ontology (GO), pathway, MSigDB,
#'         and EnrichR enrichment with multiple testing correction
#'         (\code{\link{getEnrichedGO}}, \code{\link{getEnrichedPATH}},
#'         \code{\link{test_enrichment}})
#'   \item \strong{Sequence Analysis:} Retrieve sequences of peak regions for
#'         motif discovery and pattern matching (\code{\link{getAllPeakSequence}},
#'         \code{\link{write2FASTA}})
#'   \item \strong{Statistical Testing:} Hypergeometric tests, permutation
#'         tests, and overlap significance testing
#'   \item \strong{Visualization:} Venn diagrams, metagene plots, feature-aligned
#'         heatmaps, and signal profile plots
#'   \item \strong{Peak Comparison:} Find overlaps between multiple peak sets,
#'         assess reproducibility, and generate comparison statistics
#'   \item \strong{Global Options:} Convenient global settings for genome,
#'         species, TxDb, and EnsDb objects (\code{\link{setChIPpeakAnnoGlobals}})
#' }
#' 
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Workflow-Based Analysis:} Systematic pipeline from quality
#'         assessment through annotation to downstream enrichment analysis
#'   \item \strong{Factor-Specific Annotation:} Intelligent prioritization
#'         strategies tailored to transcription factors, histone marks, Pol II,
#'         and other factor types
#'   \item \strong{Bi-directional Promoter Detection:} Unique function
#'         (\code{\link{peaksNearBDP}}) to identify peaks near bi-directional promoters
#'   \item \strong{Enhancer Detection:} Use DNA interaction data (3C/HiC) or
#'         ENCODE CRE annotations to identify enhancers associated with peaks
#'         (\code{\link{findEnhancers}})
#'   \item \strong{Epigenomic Database Integration:} Direct access to ENCODE
#'         cCRE annotations (human hg38, mouse mm10) and Roadmap Epigenomics
#'         data (human hg19) with automatic liftOver support
#'   \item \strong{Signal Analysis:} Extract and visualize read coverage or
#'         signals around genomic features
#'   \item \strong{Pattern/Motif Enrichment:} Identify enriched DNA patterns
#'         with statistical testing, compatible with HOMER, STREME, FIMO,
#'         and R packages (memes, monaLisa, rGADEM)
#'   \item \strong{Comprehensive Integration:} Works seamlessly with the
#'         Bioconductor ecosystem (TxDb, EnsDb, BSgenome, biomaRt, etc.) and
#'         external tools (HOMER, MEME Suite)
#'   \item \strong{Automatic Genome Assembly Handling:} Detects genome assembly
#'         from peaks or global options, performs automatic liftOver when needed
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
#' \strong{Quality Control and Filtering:}
#' \itemize{
#'   \item \code{\link{assessPeaks}} - Quality assessment per replicate
#'   \item \code{\link{filterPeaks}} - Filter peaks based on quality metrics
#'   \item \code{\link{IDRfilter}} - IDR-based consensus peak creation (2 replicates)
#'   \item \code{\link{findOverlapsOfPeaks}} - Consensus peak creation (3+ replicates)
#' }
#' 
#' \strong{Annotation Functions:}
#' \itemize{
#'   \item \code{\link{getGenomicAnnotation}} - Preliminary annotation with
#'         genomic feature distribution
#'   \item \code{\link{annotateHierarchically}} - Factor-specific hierarchical
#'         annotation with multiple strategies
#'   \item \code{\link{annotatePeakInBatch}} - Main function for peak annotation
#'   \item \code{\link{annoPeaks}} - Alternative annotation function with flexible
#'         binding types
#'   \item \code{\link{annotatePeaksWithcCRE}} - Annotate with ENCODE cCRE data
#'   \item \code{\link{annotatePeaksWithRoadmap}} - Annotate with Roadmap
#'         Epigenomics data
#'   \item \code{\link{findEnhancers}} - Identify enhancer-associated peaks
#'   \item \code{\link{peaksNearBDP}} - Find peaks near bi-directional promoters
#' }
#' 
#' \strong{Annotation Data Preparation:}
#' \itemize{
#'   \item \code{\link{getAnnotation}} - Generate annotations from biomaRt
#'   \item \code{\link{annoGR}} - Convert TxDb/EnsDb to annotation GRanges
#'   \item \code{\link{listAvailablecCREs}} - List available ENCODE cCRE datasets
#'   \item \code{\link{listAvailableRoadmapEpigenomes}} - List available Roadmap
#'         epigenomes
#'   \item \code{\link{listAvailableRoadmapTissues}} - List available Roadmap
#'         tissue types
#' }
#' 
#' \strong{Enrichment Analysis:}
#' \itemize{
#'   \item \code{\link{getEnrichedGO}} - GO enrichment analysis
#'   \item \code{\link{getEnrichedPATH}} - Pathway enrichment analysis
#'   \item \code{\link{test_enrichment}} - MSigDB and EnrichR enrichment
#'   \item \code{\link{runLOLA}} - Locus Overlap Analysis (LOLA)
#' }
#' 
#' \strong{Utilities:}
#' \itemize{
#'   \item \code{\link{setChIPpeakAnnoGlobals}} - Set global options (genome,
#'         species, TxDb, EnsDb)
#'   \item \code{\link{getAllPeakSequence}} - Extract sequences from peaks
#'   \item \code{\link{write2FASTA}} - Write sequences to FASTA format
#' }
#' 
#' For comprehensive examples and tutorials, see the package vignette:
#' \code{vignette("ChIPpeakAnno", package = "ChIPpeakAnno")}
#' 
#' @keywords package
"_PACKAGE"
