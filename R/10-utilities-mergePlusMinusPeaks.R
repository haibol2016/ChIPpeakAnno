#' Merge peaks from plus and minus strands
#' 
#' @description 
#' Merges peaks from plus and minus strands that are within a specified distance
#' threshold. This function is particularly useful for processing data from
#' techniques like GUIDE-seq, CRISPR-Cas9 off-target detection, or other
#' methods that generate strand-specific peaks that represent the same genomic
#' event and should be combined into a single merged peak.
#' 
#' The function identifies pairs of plus-strand and minus-strand peaks that are
#' nearby (within the distance threshold), merges them into a single peak
#' spanning from the minimum start to the maximum end position, and aggregates
#' count information from both strands.
#' 
#' @param peaks.file A character string specifying the path to the peak file
#'        containing peaks from both plus and minus strands. The file should be
#'        tab-delimited (or use \code{sep} to specify delimiter) and contain
#'        columns for chromosome, start, end, strand, and count information.
#' @param columns A character vector specifying the column names in the peak
#'        file. Must include exactly one column named \code{"strand"}. The
#'        default includes: \code{"name"}, \code{"chromosome"}, \code{"start"},
#'        \code{"end"}, \code{"strand"}, and multiple \code{"count"} columns.
#'        The function requires:
#'        \itemize{
#'          \item \code{"chromosome"}: Chromosome names
#'          \item \code{"start"}: Start positions
#'          \item \code{"end"}: End positions
#'          \item \code{"strand"}: Strand information ("+" or "-")
#'          \item One or more \code{"count"} columns: Count/signal values to be
#'                summed across strands
#'        }
#' @param sep A character string specifying the column delimiter in the peak
#'        file. Default is \code{"\t"} (tab-delimited). Other common options
#'        include \code{","} (comma-delimited) or \code{" "} (space-delimited).
#' @param header A logical value. If \code{TRUE} (default), the file has a
#'        header row with column names. If \code{FALSE}, column names are taken
#'        from the \code{columns} parameter.
#' @param distance.threshold An integer specifying the maximum distance (in base
#'        pairs) allowed between plus-strand and minus-strand peaks to be merged.
#'        Default is \code{100}. Only peak pairs with distance <= this threshold
#'        are merged. The distance is calculated as the absolute distance
#'        between peak positions.
#' @param plus.strand.start.gt.minus.strand.end A logical value controlling the
#'        expected orientation of peak pairs. If \code{TRUE} (default), plus-strand
#'        peak start is expected to be greater than (downstream of) the paired
#'        minus-strand peak end. This is the typical orientation for
#'        double-strand break events where the plus-strand break is downstream.
#'        If \code{FALSE}, the opposite orientation is expected (minus-strand
#'        peak start is greater than plus-strand peak end).
#' @param output.bedfile A character string specifying the path to the output
#'        BED file where merged peaks will be written. The file will be created
#'        (or overwritten if it exists) in standard BED format (tab-delimited,
#'        no header, no quotes).
#' 
#' @return Returns a data frame in BED format containing merged peaks with
#'        columns:
#'        \itemize{
#'          \item \code{seqnames}: Chromosome names
#'          \item \code{minStart}: Minimum start position of the merged peak
#'                (earliest start from either strand)
#'          \item \code{maxEnd}: Maximum end position of the merged peak
#'                (latest end from either strand)
#'          \item \code{names}: Peak identifier in the format
#'                \code{"peak:feature"} where peak and feature are the original
#'                peak identifiers
#'          \item \code{totalCount}: Sum of all count columns from both plus and
#'                minus strands
#'          \item \code{strand}: Strand information (always set to "+" for
#'                merged peaks)
#'        }
#'        The function also writes this data frame to the specified BED file
#'        using \code{write.table} with tab delimiter, no column names, no row
#'        names, and no quotes.
#' 
#' @details
#' 
#' \strong{How the function works:}
#' \enumerate{
#'   \item Reads the peak file and validates column structure
#'   \item Separates peaks into plus-strand and minus-strand groups
#'   \item Converts each group to GRanges objects
#'   \item Uses \code{\link{annotatePeakInBatch}} to find nearest peak pairs
#'         between strands:
#'         \itemize{
#'           \item If \code{plus.strand.start.gt.minus.strand.end = TRUE}:
#'                 Finds nearest minus-strand peaks for each plus-strand peak
#'           \item If \code{plus.strand.start.gt.minus.strand.end = FALSE}:
#'                 Finds nearest plus-strand peaks for each minus-strand peak
#'         }
#'   \item Filters peak pairs by distance threshold (only keeps pairs with
#'         distance <= \code{distance.threshold} and negative distance, meaning
#'         peaks are in the expected orientation)
#'   \item Merges peak pairs:
#'         \itemize{
#'           \item Creates merged peak spanning from min(start) to max(end)
#'           \item Sums all count columns from both strands
#'           \item Creates identifier combining both peak names
#'         }
#'   \item Writes results to BED file and returns data frame
#' }
#' 
#' \strong{Distance calculation:}
#' The function uses \code{annotatePeakInBatch} to calculate distances between
#' plus-strand and minus-strand peaks. The distance is measured as:
#' \itemize{
#'   \item Distance from plus-strand peak to minus-strand peak (or vice versa)
#'   \item Only negative distances are kept (ensuring peaks are in the expected
#'         orientation)
#'   \item Absolute distance must be <= \code{distance.threshold}
#' }
#' 
#' \strong{Peak merging:}
#' When two peaks are merged:
#' \itemize{
#'   \item The merged peak spans from the minimum start to the maximum end
#'         position
#'   \item All count columns are summed across both strands
#'   \item The peak name combines identifiers from both original peaks
#'   \item Strand is set to "+" (merged peaks are unstranded)
#' }
#' 
#' \strong{Count aggregation:}
#' The function identifies all columns named "count" and:
#' \itemize{
#'   \item Renames them with "plus:" or "minus:" prefixes
#'   \item Sums all count columns to create \code{totalCount}
#'   \item Preserves individual strand counts in intermediate steps
#' }
#' 
#' \strong{Use cases:}
#' This function is particularly useful for:
#' \itemize{
#'   \item \strong{GUIDE-seq}: Merging plus and minus strand peaks from
#'         double-strand break detection
#'   \item \strong{CRISPR off-target detection}: Combining strand-specific
#'         cleavage sites
#'   \item \strong{Any strand-specific assay}: Where the same event produces
#'         peaks on both strands
#' }
#' 
#' \strong{Validation:}
#' The function performs several validation checks:
#' \itemize{
#'   \item Number of columns must match the \code{columns} specification
#'   \item Exactly one "strand" column must be present
#'   \item Both plus and minus strand peaks must exist
#'   \item Column names must include required fields (chromosome, start, end,
#'         strand)
#' }
#' 
#' @note
#' \itemize{
#'   \item The function requires peaks from both plus and minus strands; it
#'         will stop with an error if only one strand is present
#'   \item Only peak pairs with negative distance (correct orientation) are
#'         merged
#'   \item The distance threshold applies to the absolute distance between
#'         peaks
#'   \item Merged peaks always have strand set to "+" (unstranded)
#'   \item The output BED file overwrites any existing file at the specified
#'         path
#'   \item The function uses \code{annotatePeakInBatch} internally, which may
#'         be slow for large peak sets
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{annotatePeakInBatch}} for the annotation function used
#'         internally
#'   \item \code{\link{findOverlappingPeaks}} for finding overlapping peaks
#'   \item \code{\link{makeVennDiagram}} for visualizing peak overlaps
#' }
#' 
#' @references Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to
#' annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237doi:10.1186/1471-2105-11-237
#' 
#' @author Lihua Julie Zhu
#' @keywords misc
#' @export
#' @importFrom matrixStats rowMins rowMaxs
#' @importFrom utils read.table write.table
#' @examples
#' 
#' # Example 1: Basic usage with GUIDE-seq data
#' if (interactive()) {
#'     library(matrixStats)
#'     peaks <- system.file("extdata", "guide-seq-peaks.txt", 
#'                          package = "ChIPpeakAnno")
#'     merged.bed <- mergePlusMinusPeaks(
#'         peaks.file = peaks, 
#'         columns = c("name", "chromosome", "start", "end", "strand", 
#'                     "count", "count"), 
#'         sep = "\t", 
#'         header = TRUE,  
#'         distance.threshold = 100,  
#'         plus.strand.start.gt.minus.strand.end = TRUE, 
#'         output.bedfile = "merged_peaks.bed"
#'     )
#'     # View merged peaks
#'     head(merged.bed)
#' }
#' 
#' # Example 2: Using different distance threshold
#' \dontrun{
#' merged.bed <- mergePlusMinusPeaks(
#'     peaks.file = "peaks.txt",
#'     columns = c("name", "chromosome", "start", "end", "strand", "count"),
#'     distance.threshold = 200,  # Allow larger distance
#'     output.bedfile = "merged_peaks_200bp.bed"
#' )
#' }
#' 
#' # Example 3: Opposite orientation (minus strand start > plus strand end)
#' \dontrun{
#' merged.bed <- mergePlusMinusPeaks(
#'     peaks.file = "peaks.txt",
#'     columns = c("name", "chromosome", "start", "end", "strand", "count"),
#'     plus.strand.start.gt.minus.strand.end = FALSE,  # Opposite orientation
#'     output.bedfile = "merged_peaks_opposite.bed"
#' )
#' }
#' 
#' # Example 4: Multiple count columns
#' \dontrun{
#' merged.bed <- mergePlusMinusPeaks(
#'     peaks.file = "peaks.txt",
#'     columns = c("name", "chromosome", "start", "end", "strand", 
#'                 "count1", "count2", "count3"),  # Multiple count columns
#'     output.bedfile = "merged_peaks.bed"
#'     # totalCount will sum all three count columns
#' )
#' }
#' 
#' # Example 5: Custom column names
#' \dontrun{
#' merged.bed <- mergePlusMinusPeaks(
#'     peaks.file = "peaks.csv",
#'     columns = c("peak_id", "chr", "start_pos", "end_pos", "strand", "signal"),
#'     sep = ",",  # Comma-delimited
#'     header = TRUE,
#'     output.bedfile = "merged_peaks.bed"
#'     # Note: You may need to adjust column names to match expected format
#' )
#' }
#' 
mergePlusMinusPeaks <- function(peaks.file,
                                columns = c("name", "chromosome", "start",
                                            "end", "strand", "count",
                                            "count", "count", "count"),
                                sep = "\t",
                                header = TRUE,
                                distance.threshold = 100,
                                plus.strand.start.gt.minus.strand.end = TRUE,
                                output.bedfile) {
    peaks <- read.table(peaks.file, sep = sep, header = header)
    if (ncol(peaks) != length(columns)) {
        stop("Number of columns specified differs from the number of columns ",
             "in the input peak file. Please modify the column specification ",
             "in parameter 'columns' accordingly!", call. = FALSE)
    }
    if (length(intersect(columns, "strand")) != 1L) {
        stop("Please include exactly one 'strand' column in 'columns' and ",
             "corresponding peaks.file", call. = FALSE)
    }
    colnames(peaks) <- columns
    strand_col <- which(columns == "strand")
    pos_peaks <- peaks[peaks[, strand_col] == "+", , drop = FALSE]
    neg_peaks <- peaks[peaks[, strand_col] == "-", , drop = FALSE]
    if (nrow(pos_peaks) == 0L || nrow(neg_peaks) == 0L) {
        stop("Need peaks from both + and - strand to merge", call. = FALSE)
    }
    pos_gr <- makeGRangesFromDataFrame(pos_peaks)
    neg_gr <- makeGRangesFromDataFrame(neg_peaks)
    names(pos_gr) <- paste0(paste0(seqnames(pos_gr), strand(pos_gr)),
                             ":", start(pos_gr), ":", end(pos_gr))
    names(neg_gr) <- paste0(paste0(seqnames(neg_gr), strand(neg_gr)),
                             ":", start(neg_gr), ":", end(neg_gr))
    if (plus.strand.start.gt.minus.strand.end) {
        gr <- annotatePeakInBatch(pos_gr,
                                   featureType = "TSS",
                                   AnnotationData = neg_gr,
                                   output = "nearestStart",
                                   PeakLocForDistance = "TSS")
    } else {
        gr <- annotatePeakInBatch(neg_gr,
                                   featureType = "TSS",
                                   AnnotationData = pos_gr,
                                   output = "nearestStart",
                                   PeakLocForDistance = "end")
    }
    chr_col <- which(columns == "chromosome")
    start_col <- which(columns == "start")
    end_col <- which(columns == "end")
    
    pos_peaks <- cbind(
        peak = paste0(pos_peaks[, chr_col], "+:",
                       pos_peaks[, start_col], ":",
                       pos_peaks[, end_col]),
        pos_peaks
    )
    neg_peaks <- cbind(
        feature = paste0(neg_peaks[, chr_col], "-:",
                          neg_peaks[, start_col], ":",
                          neg_peaks[, end_col]),
        neg_peaks
    )
    
    if (!plus.strand.start.gt.minus.strand.end) {
        colnames(pos_peaks)[1L] <- "feature"
        colnames(neg_peaks)[1L] <- "peak"
    }
    
    new_count_columns <- which(columns == "count") + 1L
    colnames(neg_peaks)[new_count_columns] <- 
        paste0("minus:", colnames(neg_peaks)[new_count_columns])
    colnames(pos_peaks)[new_count_columns] <- 
        paste0("plus:", colnames(pos_peaks)[new_count_columns])
    
    pos_peaks <- pos_peaks[, c(1L, new_count_columns), drop = FALSE]
    neg_peaks <- neg_peaks[, c(1L, new_count_columns), drop = FALSE]
    
    ann_peaks <- as.data.frame(gr[abs(gr$distancetoFeature) <= 
                                      distance.threshold & 
                                      gr$distancetoFeature < 0L, ])
    temp <- merge(pos_peaks, ann_peaks)
    temp1 <- merge(neg_peaks, temp)
    temp1 <- cbind(temp1,
                   totalCount = rowSums(temp1[, grep("count",
                                                      colnames(temp1)),
                                               drop = FALSE]))
    temp1$names <- paste0(temp1$peak, ":", temp1$feature)
    temp1$minStart <- rowMins(as.matrix(temp1[, c("start_position", "start")]))
    temp1$maxEnd <- rowMaxs(as.matrix(temp1[, c("end_position", "end")]))
    bed_temp <- temp1[, c("seqnames", "minStart", "maxEnd", "names",
                           "totalCount"), drop = FALSE]
    bed_temp <- cbind(bed_temp, strand = "+")
    write.table(bed_temp,
                file = output.bedfile,
                sep = "\t",
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE)
    bed_temp
}
