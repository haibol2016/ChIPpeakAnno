#' Example ChIP-seq peak dataset: STAT1 binding sites
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing putative STAT1-binding
#' regions identified in un-stimulated cells using ChIP-seq technology. This
#' dataset is commonly used in package examples and vignettes to demonstrate
#' peak annotation workflows.
#' 
#' @name myPeakList
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names (e.g., "chr1", "chr2")}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Peak identifiers as character vector}
#'   \item{metadata}{Additional metadata columns if present}
#' }
#' 
#' @details
#' This dataset contains 11,004 STAT1 binding peaks identified from ChIP-seq
#' experiments. The data is provided in the \code{\link[GenomicRanges]{GRanges}} 
#' format, which is the standard format for genomic interval data in Bioconductor
#' packages. 
#' 
#' @source
#' Robertson G, Hirst M, Bainbridge M, Bilenky M, Zhao Y, et al. (2007)
#' Genome-wide profiles of STAT1 DNA association using chromatin
#' immunoprecipitation and massively parallel sequencing.
#' \emph{Nature Methods} 4:651-657.
#' \doi{10.1038/nmeth1068}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(myPeakList)
#' 
#' # Inspect the structure
#' myPeakList
#' 
#' # Check the number of peaks
#' length(myPeakList)
#' 
#' # View first few peaks
#' head(myPeakList)
#' 
#' # Check chromosome distribution
#' table(seqnames(myPeakList))
"myPeakList"

#' Ste12-binding sites from biological replicate 1 in yeasts
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Ste12 transcription
#' factor binding sites identified from biological replicate 1 in yeast
#' (\emph{Saccharomyces cerevisiae}) ChIP-seq experiments. This dataset is part
#' of a three-replicate series used to demonstrate overlap analysis and
#' reproducibility assessment.
#' 
#' @name Peaks.Ste12.Replicate1
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names (yeast chromosomes: "chrI", "chrII", etc.)}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Peak identifiers as character vector}
#' }
#' 
#' @details
#' Ste12 is a transcription factor involved in mating and filamentous growth in
#' yeast. These datasets contain binding sites identified from three biological
#' replicates, which can be used to:
#' \itemize{
#'   \item Demonstrate overlap analysis between replicates
#'   \item Assess reproducibility using \code{\link{findOverlapsOfPeaks}}
#'   \item Generate Venn diagrams with \code{\link{makeVennDiagram}}
#'   \item Perform permutation testing with \code{\link{peakPermTest}}
#' }
#' 
#' @references
#' Lefranois P, Euskirchen GM, Auerbach RK, Rozowsky J, Gibson T, Yellman CM,
#' Gerstein M, Snyder M (2009) Efficient yeast ChIP-Seq using multiplex
#' short-read DNA sequencing. \emph{BMC Genomics} 10:37.
#' \doi{10.1186/1471-2164-10-37}
#' 
#' @seealso
#' \code{\link{Peaks.Ste12.Replicate2}}, \code{\link{Peaks.Ste12.Replicate3}},
#' \code{\link{findOverlapsOfPeaks}}, \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load replicate 1
#' data(Peaks.Ste12.Replicate1)
#' 
#' # Inspect the data
#' Peaks.Ste12.Replicate1
#' length(Peaks.Ste12.Replicate1)
#' 
#' # Compare with other replicates
#' data(Peaks.Ste12.Replicate2)
#' data(Peaks.Ste12.Replicate3)
#' 
#' # Find overlaps between replicates
#' overlaps <- findOverlapsOfPeaks(
#'     Peaks.Ste12.Replicate1,
#'     Peaks.Ste12.Replicate2,
#'     Peaks.Ste12.Replicate3
#' )
"Peaks.Ste12.Replicate1"

#' Ste12-binding sites from biological replicate 2 in yeasts
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Ste12 transcription
#' factor binding sites identified from biological replicate 2 in yeast
#' (\emph{Saccharomyces cerevisiae}) ChIP-seq experiments. This dataset is part
#' of a three-replicate series used to demonstrate overlap analysis and
#' reproducibility assessment.
#' 
#' @name Peaks.Ste12.Replicate2
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names (yeast chromosomes: "chrI", "chrII", etc.)}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Peak identifiers as character vector}
#' }
#' 
#' @details
#' Ste12 is a transcription factor involved in mating and filamentous growth in
#' yeast. These datasets contain binding sites identified from three biological
#' replicates, which can be used to:
#' \itemize{
#'   \item Demonstrate overlap analysis between replicates
#'   \item Assess reproducibility using \code{\link{findOverlapsOfPeaks}}
#'   \item Generate Venn diagrams with \code{\link{makeVennDiagram}}
#'   \item Perform permutation testing with \code{\link{peakPermTest}}
#' }
#' 
#' @references
#' Lefranois P, Euskirchen GM, Auerbach RK, Rozowsky J, Gibson T, Yellman CM,
#' Gerstein M, Snyder M (2009) Efficient yeast ChIP-Seq using multiplex
#' short-read DNA sequencing. \emph{BMC Genomics} 10:37.
#' \doi{10.1186/1471-2164-10-37}
#' 
#' @seealso
#' \code{\link{Peaks.Ste12.Replicate1}}, \code{\link{Peaks.Ste12.Replicate3}},
#' \code{\link{findOverlapsOfPeaks}}, \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load replicate 2
#' data(Peaks.Ste12.Replicate2)
#' 
#' # Inspect the data
#' Peaks.Ste12.Replicate2
#' length(Peaks.Ste12.Replicate2)
"Peaks.Ste12.Replicate2"

#' Ste12-binding sites from biological replicate 3 in yeasts
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing Ste12 transcription
#' factor binding sites identified from biological replicate 3 in yeast
#' (\emph{Saccharomyces cerevisiae}) ChIP-seq experiments. This dataset is part
#' of a three-replicate series used to demonstrate overlap analysis and
#' reproducibility assessment.
#' 
#' @name Peaks.Ste12.Replicate3
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names (yeast chromosomes: "chrI", "chrII", etc.)}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Peak identifiers as character vector}
#' }
#' 
#' @details
#' Ste12 is a transcription factor involved in mating and filamentous growth in
#' yeast. These datasets contain binding sites identified from three biological
#' replicates, which can be used to:
#' \itemize{
#'   \item Demonstrate overlap analysis between replicates
#'   \item Assess reproducibility using \code{\link{findOverlapsOfPeaks}}
#'   \item Generate Venn diagrams with \code{\link{makeVennDiagram}}
#'   \item Perform permutation testing with \code{\link{peakPermTest}}
#' }
#' 
#' @references
#' Lefranois P, Euskirchen GM, Auerbach RK, Rozowsky J, Gibson T, Yellman CM,
#' Gerstein M, Snyder M (2009) Efficient yeast ChIP-Seq using multiplex
#' short-read DNA sequencing. \emph{BMC Genomics} 10:37.
#' \doi{10.1186/1471-2164-10-37}
#' 
#' @seealso
#' \code{\link{Peaks.Ste12.Replicate1}}, \code{\link{Peaks.Ste12.Replicate2}},
#' \code{\link{findOverlapsOfPeaks}}, \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load replicate 3
#' data(Peaks.Ste12.Replicate3)
#' 
#' # Inspect the data
#' Peaks.Ste12.Replicate3
#' length(Peaks.Ste12.Replicate3)
"Peaks.Ste12.Replicate3"

#' Example ChIP-seq peak dataset 1
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing example ChIP-seq peaks
#' for demonstration purposes. This is one of three example peak datasets
#' (\code{peaks1}, \code{peaks2}, \code{peaks3}) that can be used to demonstrate
#' overlap analysis, Venn diagram generation, and multi-set peak comparisons.
#' 
#' @name peaks1
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Peak identifiers as character vector}
#' }
#' 
#' @details
#' This dataset contains 12 example peaks that can be used with \code{peaks2} and
#' \code{peaks3} to demonstrate:
#' \itemize{
#'   \item Multi-set peak overlap analysis with \code{\link{findOverlapsOfPeaks}}
#'   \item Venn diagram generation with \code{\link{makeVennDiagram}}
#'   \item Statistical testing of overlaps with \code{\link{peakPermTest}}
#' }
#' 
#' @seealso
#' \code{\link{peaks2}}, \code{\link{peaks3}}, \code{\link{findOverlapsOfPeaks}},
#' \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(peaks1)
#' 
#' # Inspect the structure
#' peaks1
#' length(peaks1)
#' 
#' # View first few peaks
#' head(peaks1, n = 2)
#' 
#' # Compare with other peak sets
#' data(peaks2)
#' data(peaks3)
#' 
#' # Find overlaps
#' overlaps <- findOverlapsOfPeaks(peaks1, peaks2, peaks3)
"peaks1"

#' Example ChIP-seq peak dataset 2
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing example ChIP-seq peaks
#' for demonstration purposes. This is one of three example peak datasets
#' (\code{peaks1}, \code{peaks2}, \code{peaks3}) that can be used to demonstrate
#' overlap analysis, Venn diagram generation, and multi-set peak comparisons.
#' 
#' @name peaks2
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Peak identifiers as character vector}
#' }
#' 
#' @details
#' This dataset contains example peaks that can be used with \code{peaks1} and
#' \code{peaks3} to demonstrate:
#' \itemize{
#'   \item Multi-set peak overlap analysis with \code{\link{findOverlapsOfPeaks}}
#'   \item Venn diagram generation with \code{\link{makeVennDiagram}}
#'   \item Statistical testing of overlaps with \code{\link{peakPermTest}}
#' }
#' 
#' @seealso
#' \code{\link{peaks1}}, \code{\link{peaks3}}, \code{\link{findOverlapsOfPeaks}},
#' \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(peaks2)
#' 
#' # Inspect the structure
#' peaks2
#' length(peaks2)
#' 
#' # View first few peaks
#' head(peaks2, n = 2)
"peaks2"

#' Example ChIP-seq peak dataset 3
#' 
#' A \code{\link[GenomicRanges]{GRanges}} object containing example ChIP-seq peaks
#' for demonstration purposes. This is one of three example peak datasets
#' (\code{peaks1}, \code{peaks2}, \code{peaks3}) that can be used to demonstrate
#' overlap analysis, Venn diagram generation, and multi-set peak comparisons.
#' 
#' @name peaks3
#' @docType data
#' 
#' @format A \code{\link[GenomicRanges]{GRanges}} object with the following structure:
#' \describe{
#'   \item{seqnames}{Chromosome names}
#'   \item{ranges}{\code{\link[IRanges]{IRanges}} object with start and end positions}
#'   \item{strand}{Strand information ("+", "-", or "*")}
#'   \item{names}{Peak identifiers as character vector}
#' }
#' 
#' @details
#' This dataset contains example peaks that can be used with \code{peaks1} and
#' \code{peaks2} to demonstrate:
#' \itemize{
#'   \item Multi-set peak overlap analysis with \code{\link{findOverlapsOfPeaks}}
#'   \item Venn diagram generation with \code{\link{makeVennDiagram}}
#'   \item Statistical testing of overlaps with \code{\link{peakPermTest}}
#' }
#' 
#' @seealso
#' \code{\link{peaks1}}, \code{\link{peaks2}}, \code{\link{findOverlapsOfPeaks}},
#' \code{\link{makeVennDiagram}}
#' 
#' @keywords datasets
#' 
#' @examples
#' # Load the data
#' data(peaks3)
#' 
#' # Inspect the structure
#' peaks3
#' length(peaks3)
#' 
#' # View first few peaks
#' head(peaks3, n = 2)
"peaks3"

