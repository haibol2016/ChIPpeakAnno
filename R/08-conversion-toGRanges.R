#' Convert data frame to GRanges
#' 
#' @description 
#' Internal helper function to convert a data frame to a GRanges object.
#' Handles column name mapping, strand formatting, coordinate conversion (for
#' BED format), special field processing (thickStart/thickEnd, itemRgb, blocks),
#' and metadata extraction. This function is called by \code{toGRanges} methods
#' after data has been read into a data frame.
#' 
#' @param data A data frame with at least chromosome, start, and end columns.
#'        Column names should match standard genomic format conventions or be
#'        provided via \code{colNames}.
#' @param colNames Character vector of column names. If \code{NULL}, uses
#'        \code{colnames(data)}. Must contain at least \code{"space"} (or
#'        synonym), \code{"start"}, and \code{"end"}. Column name mapping is
#'        case-insensitive and accepts synonyms (e.g., "chr", "chromosome" for
#'        "space").
#' @param format Character string specifying the data format (e.g., "BED",
#'        "GFF", "MACS", etc.). Special processing is applied for BED format:
#'        \itemize{
#'          \item Coordinates are converted from 0-based to 1-based
#'          \item \code{thickStart}/\code{thickEnd} are converted to
#'                \code{thick} (IRanges)
#'          \item \code{itemRgb} is converted from "R,G,B" string to RGB codes
#'          \item \code{blockCount}/\code{blockSizes}/\code{blockStarts} are
#'                converted to \code{blocks} (IRangesList)
#'        }
#' @param ... Additional arguments (currently unused).
#' @return A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with:
#'        \itemize{
#'          \item Standard GRanges fields (seqnames, ranges, strand)
#'          \item Range names (auto-generated if missing/duplicated)
#'          \item Metadata columns from input data
#'          \item Format-specific processed fields (for BED format)
#'        }
#' @keywords internal
df2GRanges <- function(data, colNames = NULL, format = "", ...) {
    if (missing(data)) {
        stop("'data' is required!", call. = FALSE)
    }
    if (!is.data.frame(data) || ncol(data) < 3L) {
        stop("No valid data passed in. ",
             "Expected a data frame as BED format file with at least 3 fields ",
             "in the order of: chromosome, start and end. ",
             "Optional fields are name, score and strand etc. ",
             "Please refer to http://genome.ucsc.edu/FAQ/FAQformat#format1 ",
             "for details.", call. = FALSE)
    }
    if (is.null(colNames)) {
        colNames <- colnames(data)
    }
    colNames_space <-
        tolower(colNames) %in% c("space", "seqnames", "chr", "chrom",
                                 "chromosome", "chromosomes")
    if (sum(colNames_space) == 1L) {
        colNames[colNames_space] <- "space"
    }
    colNames <- gsub("^start$", "start", colNames, ignore.case = TRUE)
    colNames <- gsub("^end$", "end", colNames, ignore.case = TRUE)
    if (!all(c("space", "start", "end") %in% colNames)) {
        stop("colNames must contain space/seqnames, start and end.", 
             call. = FALSE)
    }
    if (length(colNames) < ncol(data)) {
        stop("The length of colNames is less than number of columns of data", 
             call. = FALSE)
    }
    colnames(data) <- colNames[seq_len(ncol(data))]

    getCol <- function(pattern, words, default) {
        ss <- grep(pattern, colnames(data), ignore.case = TRUE)
        if (length(ss) > 1L) {
            stop("Input data has multiple columns for ", words, 
                 " information", call. = FALSE)
        }
        if (length(ss) == 1L) {
            re <- data[, ss, drop = TRUE]
            data[, ss] <<- NULL
        } else {
            re <- default
        }
        re
    }
    ##prepare strand
    strand <- getCol("^strand$", "strand", "*")
    strand <- formatStrand(strand)

    ## Prepare name column
    names <- getCol("^name(s)?$", "names", NA)
    if (any(is.na(names)) || any(duplicated(names))) {
        message("Duplicated or NA names found. Renaming all names by numbers.")
        n <- nrow(data)
        names <- sprintf(paste0("X%0", nchar(as.character(n)), "d"), 
                        seq_len(n))
    }
    names <- make.names(names)

    ##prepare score
    #   score <- getCol("^score$", "score", 1L)
    #   if(length(score)==1) score <- rep(1, nrow(data))
    #   if(all(is.na(score))) score <- rep(1, nrow(data))
    #   score <- as.numeric(as.character(score))

    ## Prepare start, end, seqnames
    start <- data$start
    end <- data$end
    seqnames <- data$space
    if (!is.numeric(start[1L])) {
        start <- as.numeric(as.character(start))
    }
    if (!is.numeric(end[1L])) {
        end <- as.numeric(as.character(end))
    }
    if (!is.character(seqnames[1L])) {
        seqnames <- as.character(data$space)
    }
    ## Trim seqnames
    seqnames <- gsub("^\\s+|\\s+$", "", seqnames)

    gr <- GRanges(seqnames = seqnames,
                  ranges = IRanges(start = start,
                                   end = end,
                                   names = names),
                  strand = strand)
    rm(list = c("start", "end", "names", "strand", "seqnames"))
    metadata <- colnames(data)
    metadata <- metadata[!metadata %in% c("seqnames", "space", "ranges",
                                          "strand", "seqlevels",
                                          "seqlengths", "isCircular",
                                          "genome", "start",
                                          "end", "width", "element")]
    for (col in metadata) {
        mcols(gr)[, col] <- data[, col, drop = TRUE]
    }
    rm(data)
    if (format == "BED") {
        ## BED file is 0-based (start, end], convert to 1-based [start, end]
        start(gr) <- start(gr) + 1L
        if (length(gr$thickStart) > 0L && length(gr$thickEnd) > 0L) {
            if (!(all(gr$thickStart == round(gr$thickStart)) &&
                  all(gr$thickEnd == round(gr$thickEnd)))) {
                stop("This is not a standard BED file. ",
                     "Maybe it is narrowPeak or broadPeak", call. = FALSE)
            }
            gr$thick <- IRanges(gr$thickStart + 1L, gr$thickEnd)
            gr$thickStart <- NULL
            gr$thickEnd <- NULL
        }
        if (length(gr$itemRgb) > 0L) {
            ## itemRgb format: "255,0,0"
            itemRgb <- do.call(rbind, strsplit(as.character(gr$itemRgb), ","))
            gr$itemRgb <- NULL
            if (ncol(itemRgb) == 3L) {
                itemRgb <- apply(itemRgb, 1L, function(.ele) {
                    .ele <- as.numeric(.ele)
                    rgb(.ele[1L], .ele[2L], .ele[3L], maxColorValue = 255L)
                })
                gr$itemRgb <- itemRgb
            } else {
                gr$itemRgb <- NA
            }
        }
        if (length(gr$blockCount) > 0L &&
            length(gr$blockSizes) > 0L &&
            length(gr$blockStarts) > 0L) {
            blocksizes <- strsplit(as.character(gr$blockSizes), ",")
            blockstarts <- strsplit(as.character(gr$blockStarts), ",")
            blocks <- mapply(function(num, sizes, starts) {
                sizes <- as.integer(sizes)[seq_len(num)]
                starts <- as.integer(starts)[seq_len(num)]
                ir <- IRanges(starts + 1L, width = sizes)
                return(ir)
            }, gr$blockCount, blocksizes, blockstarts, SIMPLIFY = FALSE)
            gr$blocks <- IRangesList(blocks)
            gr$blockCount <- NULL
            gr$blockSizes <- NULL
            gr$blockStarts <- NULL
        }
    }
    return(gr)
}

#' Switch column names based on format
#' 
#' @description 
#' Internal helper function to map column names based on data format. Returns
#' the expected column names for each supported format, which are then used to
#' map input data columns to GRanges fields.
#' 
#' @param format Character string specifying the data format. Must be one of:
#'        \code{"BED"}, \code{"GFF"}, \code{"MACS"}, \code{"MACS2"},
#'        \code{"MACS2.broad"}, \code{"narrowPeak"}, \code{"broadPeak"},
#'        \code{"CSV"}, or \code{"others"}.
#' @param colNames Optional character vector of column names for
#'        \code{format = "others"}. Required when format is "others", ignored
#'        otherwise.
#' @return Character vector of column names expected for the specified format.
#'        The order and names correspond to the standard column order for that
#'        format. For \code{format = "others"}, returns \code{colNames} as-is.
#' @keywords internal
switchColNames <- function(format = c("BED", "GFF",
                                     "MACS", "MACS2", "MACS2.broad",
                                     "narrowPeak", "broadPeak", "CSV",
                                     "others"), 
                          colNames = NULL) {
    format <- match.arg(format)
    switch(format,
           BED = c("space", "start", "end", "names",
                   "score", "strand", "thickStart",
                   "thickEnd", "itemRgb", "blockCount",
                   "blockSizes", "blockStarts"),
           GFF = c("space", "source", "names", "start",
                   "end", "score", "strand", "frame", "group"),
           MACS = c("space", "start", "end", "length",
                    "summit", "tags", "-10*log10(pvalue)",
                    "fold_enrichment", "FDR"),
           MACS2 = c("space", "start", "end", "length",
                     "abs_summit", "pileup", "-log10(pvalue)",
                     "fold_enrichment", "-log10(qvalue)", "names"),
           MACS2.broad = c("space", "start", "end", "length",
                           "pileup", "-log10(pvalue)",
                           "fold_enrichment", "-log10(qvalue)", "names"),
           narrowPeak = c("space", "start", "end", "names",
                           "score", "strand", "signalValue",
                           "pValue", "qValue", "peak"),
           broadPeak = c("space", "start", "end", "names",
                         "score", "strand", "signalValue",
                         "pValue", "qValue"),
           others = colNames,
           colNames)
}

#' Convert various data formats to GRanges
#'
#' @description 
#' Converts various genomic data formats (BED, GFF, GTF, MACS output, etc.) or
#' annotation objects (TxDb, EnsDb) to GRanges objects. Supports file paths,
#' connections, data frames, and Bioconductor annotation objects. This is a
#' generic function with methods for different input types (character, connection,
#' data.frame, TxDb, EnsDb).
#' 
#' The function automatically handles format-specific column mappings, coordinate
#' system conversions (0-based BED to 1-based GRanges), and special field
#' processing (e.g., BED blocks, itemRgb, thickStart/thickEnd).
#'
#' @rdname toGRanges
#' @aliases toGRanges toGRanges,data.frame-method toGRanges,connection-method
#' toGRanges,character-method toGRanges,TxDb-method toGRanges,EnsDb-method
#' @param data Input data in one of the following formats:
#'        \itemize{
#'          \item \code{character}: File path to a genomic data file
#'          \item \code{connection}: A readable text-mode connection (see
#'                \code{\link[utils]{read.table}})
#'          \item \code{data.frame}: A data frame with genomic coordinates
#'          \item \code{\link[GenomicFeatures:TxDb-class]{TxDb}}: Annotation
#'                database from GenomicFeatures package
#'          \item \code{\link[ensembldb]{EnsDb}}: Annotation database from
#'                ensembldb package
#'        }
#' @param format A character string specifying the data format. Supported
#'        formats:
#'        \itemize{
#'          \item \code{"BED"}: UCSC BED format (0-based coordinates, converted
#'                to 1-based). Supports standard BED fields including thickStart,
#'                thickEnd, itemRgb, and blocks. See
#'                \url{http://genome.ucsc.edu/FAQ/FAQformat#format1}
#'          \item \code{"GFF"}, \code{"GTF"}: GFF/GTF format files. Uses
#'                \code{rtracklayer::import} for parsing. For Ensembl GFF/GTF
#'                files, it's recommended to use TxDb objects instead.
#'          \item \code{"narrowPeak"}: ENCODE narrowPeak format (BED6+4). Columns:
#'                chromosome, start, end, name, score, strand, signalValue,
#'                pValue, qValue, peak
#'          \item \code{"broadPeak"}: ENCODE broadPeak format (BED6+3). Columns:
#'                chromosome, start, end, name, score, strand, signalValue,
#'                pValue, qValue
#'          \item \code{"MACS"}: MACS1 peak calling output (tab-delimited with
#'                header). Columns: chromosome, start, end, length, summit,
#'                tags, -10*log10(pvalue), fold_enrichment, FDR
#'          \item \code{"MACS2"}: MACS2 peak calling output (tab-delimited with
#'                header). Columns: chromosome, start, end, length,
#'                abs_summit, pileup, -log10(pvalue), fold_enrichment,
#'                -log10(qvalue), name
#'          \item \code{"MACS2.broad"}: MACS2 broad peak output. Columns:
#'                chromosome, start, end, length, pileup, -log10(pvalue),
#'                fold_enrichment, -log10(qvalue), name
#'          \item \code{"CSV"}: Comma-separated values file. Must have columns:
#'                seqnames (or space/chr/chrom/chromosome), start, end, strand
#'                (optional). Header is required.
#'          \item \code{"others"}: Custom format. Requires \code{colNames} to be
#'                specified.
#'        }
#' @param feature A character string specifying the annotation type to extract
#'        when \code{data} is a TxDb or EnsDb object. Options depend on the
#'        object type:
#'        \itemize{
#'          \item For \code{TxDb}: \code{"gene"}, \code{"transcript"},
#'                \code{"exon"}, \code{"CDS"}, \code{"fiveUTR"},
#'                \code{"threeUTR"}, \code{"tRNAs"}, or \code{"geneModel"}
#'                (default is \code{"gene"})
#'          \item For \code{EnsDb}: \code{"gene"}, \code{"transcript"},
#'                \code{"exon"}, or \code{"disjointExons"} (default is
#'                \code{"gene"})
#'        }
#'        Only used when \code{data} is a TxDb or EnsDb object.
#' @param header A logical value indicating whether the file contains column
#'        names as its first line. Default behavior:
#'        \itemize{
#'          \item \code{format = "CSV"}: \code{header = TRUE} (required)
#'          \item \code{format = "MACS"}, \code{"MACS2"}, \code{"MACS2.broad"}:
#'                \code{header = TRUE} (required)
#'          \item Other formats: \code{header = FALSE} by default
#'        }
#'        If missing, the value is determined from the file format.
#' @param comment.char A character vector of length one containing a single
#'        character or an empty string. Lines starting with this character are
#'        treated as comments and ignored. Default is \code{"#"} (standard for
#'        most genomic formats). Use \code{""} to turn off comment interpretation.
#'        For MACS formats, comments starting with \code{"#"} are automatically
#'        handled.
#' @param colNames A character vector specifying column names. Required when
#'        \code{format = "others"}. Must contain at least:
#'        \itemize{
#'          \item Chromosome column: named \code{"space"}, \code{"seqnames"},
#'                \code{"chr"}, \code{"chrom"}, \code{"chromosome"}, or
#'                \code{"chromosomes"} (case-insensitive)
#'          \item \code{"start"}: Start position column
#'          \item \code{"end"}: End position column
#'        }
#'        Optional columns: \code{"strand"}, \code{"name"} or \code{"names"},
#'        \code{"score"}, and any other metadata columns.
#' @param \dots Additional parameters passed to \code{\link[utils]{read.table}}
#'        for file reading (e.g., \code{sep}, \code{quote}, \code{skip},
#'        \code{nrows}, etc.). For BED format, columns beyond 12 are ignored
#'        (set to \code{NULL}).
#' @param OrganismDb An object of class \code{\link[OrganismDbi]{OrganismDb}}.
#'        Only used when \code{data} is a TxDb object and \code{feature =
#'        "geneModel"}. Used to extract gene symbols (SYMBOL) from the
#'        OrganismDb database and add them to the annotation metadata.
#' 
#' @return Returns a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        containing the converted genomic ranges. The object includes:
#'        \itemize{
#'          \item Standard GRanges fields: \code{seqnames}, \code{ranges}
#'                (IRanges with start, end, width), \code{strand}
#'          \item Range names: Automatically generated if missing or duplicated
#'                (format: "X0001", "X0002", etc.)
#'          \item Metadata columns: All additional columns from the input data
#'                (except reserved GRanges column names) are preserved in
#'                \code{mcols(gr)}
#'        }
#'        
#'        \strong{Format-specific additions:}
#'        \itemize{
#'          \item \code{format = "BED"}: 
#'                \itemize{
#'                  \item Coordinates converted from 0-based [start, end) to
#'                        1-based [start, end]
#'                  \item \code{thickStart}/\code{thickEnd} converted to
#'                        \code{thick} (IRanges object)
#'                  \item \code{itemRgb} converted from "R,G,B" string to RGB
#'                        color codes
#'                  \item \code{blockCount}/\code{blockSizes}/\code{blockStarts}
#'                        converted to \code{blocks} (IRangesList object)
#'                }
#'          \item \code{format = "narrowPeak"}: Includes \code{signalValue},
#'                \code{pValue}, \code{qValue}, \code{peak} in metadata
#'          \item \code{format = "broadPeak"}: Includes \code{signalValue},
#'                \code{pValue}, \code{qValue} in metadata
#'          \item \code{format = "MACS"}: Includes MACS-specific fields in
#'                metadata
#'          \item \code{format = "MACS2"}: Includes MACS2-specific fields in
#'                metadata
#'        }
#' 
#' @details
#' 
#' \strong{Supported input types:}
#' \describe{
#'   \item{\code{character} (file path)}{Reads the file and converts based on
#'         \code{format}. Automatically detects column classes for BED/GFF
#'         formats.}
#'   \item{\code{connection}}{Reads from a text connection (e.g., from
#'         \code{textConnection()}). Useful for processing data in memory.}
#'   \item{\code{data.frame}}{Converts directly using column names. Column name
#'         mapping is flexible (case-insensitive, accepts synonyms).}
#'   \item{\code{TxDb}}{Extracts annotation features using
#'         \code{\link{TxDb2GR}}. Supports all standard feature types.}
#'   \item{\code{EnsDb}}{Extracts annotation features using
#'         \code{\link{EnsDb2GR}}. Supports gene, transcript, exon, and
#'         disjointExons.}
#' }
#' 
#' \strong{Column name mapping:}
#' The function automatically maps common column name variations:
#' \itemize{
#'   \item Chromosome: \code{space}, \code{seqnames}, \code{chr}, \code{chrom},
#'         \code{chromosome}, \code{chromosomes} (all mapped to \code{space})
#'   \item Start/End: Case-insensitive matching
#'   \item Strand: Automatically formatted to "+", "-", or "*"
#'   \item Names: \code{name} or \code{names} (both accepted)
#' }
#' 
#' \strong{BED format handling:}
#' BED files use 0-based coordinates [start, end), which are converted to 1-based
#' [start, end] for GRanges:
#' \itemize{
#'   \item \code{start} is incremented by 1
#'   \item \code{thickStart}/\code{thickEnd} are converted to \code{thick}
#'         (IRanges) with coordinates incremented by 1
#'   \item \code{blockStarts} are incremented by 1 when creating \code{blocks}
#'   \item \code{itemRgb} "R,G,B" strings are converted to RGB color codes
#' }
#' 
#' \strong{Name handling:}
#' \itemize{
#'   \item If names are missing: Auto-generated as "X0001", "X0002", etc.
#'   \item If names contain NA or duplicates: Auto-renamed with zero-padded
#'         numbers
#'   \item Names are made valid using \code{make.names()}
#' }
#' 
#' \strong{Metadata preservation:}
#' All columns from the input data (except reserved GRanges column names:
#' \code{seqnames}, \code{space}, \code{ranges}, \code{strand}, \code{start},
#' \code{end}, \code{width}, etc.) are preserved as metadata columns in the
#' returned GRanges object.
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{TxDb2GR}} for TxDb conversion,
#'          \code{\link{EnsDb2GR}} for EnsDb conversion,
#'          \code{\link[rtracklayer]{import}} for GFF/GTF import,
#'          \code{\link[GenomicRanges]{GRanges}} for the output class
#' @keywords misc
#' @exportMethod toGRanges
#' @export toGRanges
#' @importFrom utils read.csv read.table
#' @importFrom grDevices rgb
#' @importFrom rtracklayer import
#' @examples
#' \dontrun{
#' ## Example 1: Convert MACS output file
#' macs <- system.file("extdata", "MACS_peaks.xls", package = "ChIPpeakAnno")
#' macsOutput <- toGRanges(macs, format = "MACS")
#' 
#' ## Example 2: Convert from connection
#' macs <- readLines(macs)
#' macs <- textConnection(macs)
#' macsOutput <- toGRanges(macs, format = "MACS")
#' close(macs)
#' 
#' ## Example 3: Convert BED file
#' bed <- system.file("extdata", "MACS_output.bed", package = "ChIPpeakAnno")
#' bedGR <- toGRanges(bed, format = "BED")
#' 
#' ## Example 4: Convert narrowPeak file (ENCODE format)
#' narrowPeak <- system.file("extdata", "peaks.narrowPeak", 
#'                            package = "ChIPpeakAnno")
#' peaks <- toGRanges(narrowPeak, format = "narrowPeak")
#' 
#' ## Example 5: Convert broadPeak file (ENCODE format)
#' broadPeak <- system.file("extdata", "TAF.broadPeak", 
#'                          package = "ChIPpeakAnno")
#' peaks <- toGRanges(broadPeak, format = "broadPeak")
#' 
#' ## Example 6: Convert CSV file
#' csv <- system.file("extdata", "peaks.csv", package = "ChIPpeakAnno")
#' peaks <- toGRanges(csv, format = "CSV")
#' 
#' ## Example 7: Convert MACS2 output
#' macs2 <- system.file("extdata", "MACS2_peaks.xls", 
#'                      package = "ChIPpeakAnno")
#' peaks <- toGRanges(macs2, format = "MACS2")
#' 
#' ## Example 8: Convert GFF file (uses rtracklayer)
#' gff <- system.file("extdata", "GFF_peaks.gff", package = "ChIPpeakAnno")
#' peaks <- toGRanges(gff, format = "GFF")
#' 
#' ## Example 9: Convert from EnsDb object
#' library(EnsDb.Hsapiens.v75)
#' genes <- toGRanges(EnsDb.Hsapiens.v75, feature = "gene")
#' transcripts <- toGRanges(EnsDb.Hsapiens.v75, feature = "transcript")
#' 
#' ## Example 10: Convert from TxDb object
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' genes <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature = "gene")
#' exons <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature = "exon")
#' 
#' ## Example 11: Convert data.frame
#' macs <- system.file("extdata", "MACS_peaks.xls", package = "ChIPpeakAnno")
#' macs_df <- read.delim(macs, comment.char = "#")
#' peaks <- toGRanges(macs_df)
#' 
#' ## Example 12: Custom format with colNames
#' custom_data <- data.frame(
#'     chr = c("chr1", "chr2"),
#'     start_pos = c(1000, 2000),
#'     end_pos = c(1500, 2500),
#'     strand_info = c("+", "-")
#' )
#' colNames <- c("space", "start", "end", "strand")
#' names(custom_data) <- colNames
#' peaks <- toGRanges(custom_data, format = "others", colNames = colNames)
#' }
#'

#'
setGeneric("toGRanges", function(data, ...) standardGeneric("toGRanges"))
setMethod("toGRanges", "data.frame",
          function(data, colNames = NULL, ...) {
              this.call <- match.call(expand.dots = TRUE)
              this.call[[1]] <- df2GRanges
              if (length(this.call$format) == 0L) {
                  this.call$format <- "data.frame"
              }
              gr <- eval(this.call, parent.frame())
              return(gr)
})


#' Message helper for GTF/GFF files
#' 
#' @description 
#' Internal helper function to display a message recommending the use of
#' TxDb objects for GTF/GFF files from Ensembl. This is called when users
#' import GFF/GTF files directly, suggesting they use \code{makeTxDbFromGFF}
#' instead for better annotation handling.
#' 
#' @param con Character string or connection object representing the file path.
#'        Used in the message to show the recommended code pattern.
#' @keywords internal
message4GTF <- function(con) {
    message("If you are importing files downloaded from Ensembl, ",
            "it will be better to import the files into a TxDb object, ",
            "and then convert to GRanges by toGRanges. ",
            "Here is the sample code:\n",
            "library(GenomicFeatures)\n",
            "txdb <- makeTxDbFromGFF('", con, "')\n",
            "anno <- toGRanges(txdb, format='gene')")
}

#' @importFrom rtracklayer import
#' @importFrom utils read.table
#' @rdname toGRanges
#' @aliases toGRanges,connection-method
setMethod("toGRanges", "connection",
          function(data, format = c("BED", "GFF", "GTF",
                                   "MACS", "MACS2", "MACS2.broad",
                                   "narrowPeak", "broadPeak", "CSV",
                                   "others"),
                   header = FALSE, comment.char = "#", colNames = NULL, ...) {
              format <- match.arg(format)
              if (format %in% c("GFF", "GTF")) {
                  message4GTF("path/to/your/GFF")
                  gr <- import(data, format = format)
                  return(gr)
              }
              if (format %in% c("narrowPeak", "broadPeak")) {
                  data <- read.table(data, header = FALSE,
                                     fill = TRUE, stringsAsFactors = FALSE)
                  data <- data[!grepl("track|browser", data[, 1L]), 
                              seq_len(ncol(data)), drop = FALSE]
                  classes <- c("character", "integer", "integer", "character",
                              "integer", "character", "numeric", "numeric",
                              "numeric", "integer")[seq_len(ncol(data))]
                  for (i in seq_len(ncol(data))) {
                      class(data[, i]) <- mode(data[, i]) <- classes[i]
                  }
              } else if (format == "CSV") {
                  data <- read.csv(data, header = TRUE)
                  if (any(colnames(data) == "X") &&
                      all(!colnames(data) %in% c("names", "name"))) {
                      colnames(data)[colnames(data) == "X"] <- "names"
                  }
                  colNames <- colnames(data)
              } else {
                  if (format %in% c("MACS", "MACS2", "MACS2.broad")) {
                      header <- TRUE
                      comment.char <- "#"
                  }
                  data <- read.table(data, header = header,
                                     comment.char = comment.char,
                                     ...)
              }
              colNames <- switchColNames(format, colNames)
              if (is.null(colNames)) {
                  stop("colNames is required for unknown format.", call. = FALSE)
              }
              gr <- df2GRanges(data, colNames, format, ...)
              return(gr)
          })

#' @rdname toGRanges
#' @aliases toGRanges,TxDb-method
setMethod("toGRanges", "TxDb",
          function(data, feature = c("gene", "transcript", "exon",
                                    "CDS", "fiveUTR", "threeUTR",
                                    "tRNAs", "geneModel"),
                   OrganismDb, ...) {
              feature <- match.arg(feature)
              if (!missing(OrganismDb)) {
                  return(TxDb2GR(data, feature, OrganismDb))
              } else {
                  return(TxDb2GR(data, feature))
              }
          })

#' @rdname toGRanges
#' @aliases toGRanges,EnsDb-method
setMethod("toGRanges", "EnsDb",
          function(data,
                   feature = c("gene", "transcript", "exon", "disjointExons"),
                   ...) {
              feature <- match.arg(feature)
              return(EnsDb2GR(data, feature))
          })

#' @rdname toGRanges
#' @aliases toGRanges,character-method
setMethod("toGRanges", "character",
          function(data, format = c("BED", "GFF", "GTF",
                                   "MACS", "MACS2", "MACS2.broad",
                                   "narrowPeak", "broadPeak", "CSV",
                                   "others"),
                   header = FALSE, comment.char = "#", colNames = NULL, ...) {
              format <- match.arg(format)
              if (format %in% c("GFF", "GTF")) {
                  message4GTF(data)
                  gr <- import(data, format = format)
                  return(gr)
              } else if (format %in% c("narrowPeak", "broadPeak")) {
                  data <- read.table(data, header = FALSE,
                                     fill = TRUE, stringsAsFactors = FALSE)
                  data <- data[!grepl("track|browser", data[, 1L]), 
                              seq_len(ncol(data)), drop = FALSE]
                  classes <- c("character", "integer", "integer", "character",
                              "integer", "character", "numeric", "numeric",
                              "numeric", "integer")[seq_len(ncol(data))]
                  for (i in seq_len(ncol(data))) {
                      class(data[, i]) <- mode(data[, i]) <- classes[i]
                  }
              } else if (format == "CSV") {
                  data <- read.csv(data, header = TRUE)
                  if (any(colnames(data) == "X") &&
                      all(!colnames(data) %in% c("names", "name"))) {
                      colnames(data)[colnames(data) == "X"] <- "names"
                  }
                  colNames <- colnames(data)
              } else {
                  if (format %in% c("MACS", "MACS2", "MACS2.broad")) {
                      header <- TRUE
                      comment.char <- "#"
                  }
                  tab5rows <- read.table(data, header = header,
                                         comment.char = comment.char, ...,
                                         nrows = 5L)
                  classes <- vapply(tab5rows, class, character(1))
                  if (format == "BED") {
                      ## Check class of column 2 and 3
                      if (classes[2L] != "integer" || classes[3L] != "integer") {
                          stop("No valid data passed in. ",
                               "Expected a data frame as BED format file with ",
                               "at least 3 fields in the order of: ",
                               "chromosome, start and end. ",
                               "Optional fields are name, score and strand etc. ",
                               "Column 2 and 3 must be integer. ",
                               "Please refer to ",
                               "http://genome.ucsc.edu/FAQ/FAQformat#format1 ",
                               "for details.", call. = FALSE)
                      }
                      if (!is.na(classes[5L])) {
                          classes[5L] <- "numeric"
                      }
                      classes[1L] <- "character"
                  } else {
                      if (format == "GFF") {
                          ## Check class of column 4 and 5
                          if (classes[4L] != "integer" || classes[5L] != "integer") {
                              stop("No valid data passed in. ",
                                   "Expected a data frame as GFF format file ",
                                   "with 9 fields in the order of: ",
                                   "seqname, source, feature, start, end, ",
                                   "score, strand, frame and group. ",
                                   "Column 4 and 5 must be integer. ",
                                   "Please refer to ",
                                   "http://genome.ucsc.edu/FAQ/FAQformat#format1 ",
                                   "for details.", call. = FALSE)
                          }
                          classes[1L] <- "character"
                      } else {
                          if (format %in% c("MACS", "MACS2", "MACS2.broad")) {
                              ## Do nothing
                          } else {
                              if (header && is.null(colNames)) {
                                  ## format == "others"
                                  colNames <- colnames(tab5rows)
                              }
                          }
                      }
                  }
                  if (format == "BED" && length(classes) > 12L) {
                      classes[13L:length(classes)] <- rep("NULL", 
                                                         length(classes) - 12L)
                  }
                  
                  data <- read.table(data, header = header,
                                     comment.char = comment.char, ...,
                                     colClasses = classes)
                  rm(list = c("tab5rows", "classes"))
              }
              colNames <- switchColNames(format, colNames)
              gr <- df2GRanges(data, colNames, format, ...)
              return(gr)
          })
