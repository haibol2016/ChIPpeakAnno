#' LOLA: Locus Overlap Analysis Integration
#' 
#' @description 
#' Re-exported functions from the \code{LOLA} package for enrichment analysis
#' of genomic regions. These functions allow you to test whether your regions
#' of interest are enriched for overlap with reference region sets from public
#' databases (ENCODE, Roadmap Epigenomics, etc.).
#' 
#' The LOLA package provides a complementary approach to \code{\link{peakPermTest}}
#' for enrichment analysis:
#' \itemize{
#'   \item \code{peakPermTest}: Permutation-based testing with bias correction
#'   \item \code{runLOLA}: Database-driven enrichment against many reference sets
#' }
#' 
#' @details
#' 
#' \strong{Key LOLA Functions:}
#' \itemize{
#'   \item \code{loadRegionDB}: Load a pre-built or custom region database
#'   \item \code{runLOLA}: Run enrichment analysis against the database
#' }
#' 
#' \strong{When to Use LOLA vs peakPermTest:}
#' \itemize{
#'   \item Use \code{runLOLA} when: Testing enrichment against many reference
#'         sets, using pre-built databases, need fast analysis
#'   \item Use \code{peakPermTest} when: Comparing two specific peak sets,
#'         need bias correction, want permutation-based testing
#' }
#' 
#' For complete LOLA documentation, see \code{\link[LOLA]{LOLA-package}}.
#' 
#' @name LOLA
#' @seealso \code{\link{peakPermTest}} for permutation-based overlap testing,
#'          \code{\link{buildLOLAregionDB}} for building custom databases
#' @keywords internal
NULL

#' Load a LOLA region database
#' 
#' @description 
#' Loads a LOLA region database from disk. This function is re-exported from
#' the \code{LOLA} package to provide convenient access within ChIPpeakAnno.
#' 
#' @param dbLocation Character string specifying the path to the region database
#'        directory. The directory should contain subdirectories for each
#'        collection, with BED files and metadata.
#' @param ... Additional arguments passed to \code{\link[LOLA]{loadRegionDB}}
#' 
#' @return A LOLA region database object that can be used with
#'         \code{\link{runLOLA}}
#' 
#' @details
#' 
#' The region database directory should have the following structure:
#' \preformatted{
#' regionDB/
#' ├── collection1/
#' │   ├── collection.txt
#' │   └── regions/
#' │       ├── region_set1.bed
#' │       └── region_set1_description.txt
#' ├── collection2/
#' │   └── ...
#' }
#' 
#' @seealso \code{\link[LOLA]{loadRegionDB}} for complete documentation in the
#'          LOLA package
#' @seealso \code{\link{runLOLA}} for running enrichment analysis
#' @seealso \code{\link{buildLOLAregionDB}} for building custom databases
#' 
#' @examples
#' \dontrun{
#' # Load a pre-built LOLA Core database
#' regionDB <- loadRegionDB("/path/to/LOLA_Core")
#' 
#' # Or load a custom database
#' regionDB <- loadRegionDB("/path/to/my_custom_regionDB")
#' }
#' 
#' @importFrom LOLA loadRegionDB
#' @export
loadRegionDB <- function(dbLocation, ...) {
    if (!requireNamespace("LOLA", quietly = TRUE)) {
        stop("LOLA package is required. Install with: ",
             "BiocManager::install('LOLA')", call. = FALSE)
    }
    LOLA::loadRegionDB(dbLocation, ...)
}

#' Run LOLA enrichment analysis
#' 
#' @description 
#' Tests whether query regions are enriched for overlap with reference region
#' sets in a LOLA database. This function is re-exported from the \code{LOLA}
#' package to provide convenient access within ChIPpeakAnno.
#' 
#' The function uses Fisher's exact test to compare observed overlaps between
#' query regions and reference sets against expected overlaps based on the
#' universe (background) regions.
#' 
#' @param userSets A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        containing query regions to test for enrichment
#' @param universe A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        containing all regions that could have been included in \code{userSets}.
#'        This defines the background for statistical testing.
#' @param regionDB A LOLA region database object (from \code{\link{loadRegionDB}})
#' @param ... Additional arguments passed to \code{\link[LOLA]{runLOLA}}
#' 
#' @return A data frame with enrichment results containing:
#' \itemize{
#'   \item \code{pValueLog}: Negative log10 p-value
#'   \item \code{oddsRatio}: Odds ratio for enrichment
#'   \item \code{support}: Number of query regions overlapping the reference set
#'   \item \code{rnk}: Rank based on p-value
#'   \item Additional metadata about the reference sets
#' }
#' 
#' @details
#' 
#' \strong{Universe Selection:}
#' The universe should represent all regions that could have been in your query
#' set. For example:
#' \itemize{
#'   \item If analyzing differential peaks: universe = all peaks tested
#'   \item If analyzing called peaks: universe = all regions with coverage
#'   \item If analyzing specific regions: universe = all testable regions
#' }
#' 
#' \strong{Comparison with peakPermTest:}
#' \itemize{
#'   \item \code{runLOLA}: Fast, database-driven, tests many reference sets
#'   \item \code{peakPermTest}: Permutation-based, bias-corrected, tests one
#'         comparison
#' }
#' 
#' @seealso \code{\link[LOLA]{runLOLA}} for complete documentation in the
#'          LOLA package
#' @seealso \code{\link{loadRegionDB}} for loading databases
#' @seealso \code{\link{peakPermTest}} for permutation-based overlap testing
#' 
#' @examples
#' \dontrun{
#' # Load region database
#' regionDB <- loadRegionDB("/path/to/LOLA_Core")
#' 
#' # Define query regions
#' query_regions <- GRanges("chr1", IRanges(1000000, 2000000))
#' 
#' # Define universe (all testable regions)
#' universe <- all_peaks_from_experiment
#' 
#' # Run enrichment analysis
#' results <- runLOLA(query_regions, universe, regionDB)
#' 
#' # View top enriched sets
#' head(results[order(results$pValueLog, decreasing = TRUE), ])
#' }
#' 
#' @importFrom LOLA runLOLA
#' @export
runLOLA <- function(userSets, universe, regionDB, ...) {
    if (!requireNamespace("LOLA", quietly = TRUE)) {
        stop("LOLA package is required. Install with: ",
             "BiocManager::install('LOLA')", call. = FALSE)
    }
    LOLA::runLOLA(userSets, universe, regionDB, ...)
}


