#' Class \code{annoGR}
#' 
#' An S4 class extending \code{\link[GenomicRanges]{GRanges}} to represent
#' genomic annotation data (genes, transcripts, exons, etc.) with metadata
#' about the annotation source, creation date, and feature type. Objects of
#' this class can be used by \code{\link{annotatePeakInBatch}} and
#' \code{\link{annoPeaks}} for peak annotation.
#' 
#' @details
#' The \code{annoGR} class adds metadata to standard GRanges objects to track:
#' \itemize{
#'   \item \code{source}: Where the annotation comes from (e.g., "EnsDb.Hsapiens.v86",
#'         "TxDb.Hsapiens.UCSC.hg19.knownGene", or a custom character string)
#'   \item \code{date}: When the annotation object was created (defaults to
#'         \code{Sys.Date()} if not provided)
#'   \item \code{feature}: Type of annotation feature (e.g., "gene", "transcript", "exon")
#'   \item \code{mdata}: Additional metadata as a data frame with columns "name" and
#'         "value" (automatically extracted from TxDb/EnsDb metadata tables if not
#'         provided)
#' }
#' 
#' This class is particularly useful when working with \code{\link[GenomicFeatures]{TxDb}}
#' or \code{\link[ensembldb]{EnsDb}} annotation packages as well as custom GRanges-based
#' annotations, as it preserves information about the annotation source and version.
#' It provides a consistent interface for working with annotation data from different
#' sources.
#' 
#' @name annoGR-class
#' @rdname annoGR
#' @aliases annoGR-class annoGR
#' @docType class
#' 
#' @param ranges An object that can be converted to \code{annoGR}. Supported types:
#'        \itemize{
#'          \item \code{\link[GenomicRanges:GRanges-class]{GRanges}}: Custom annotation
#'                data. If names are missing, they will be automatically generated.
#'          \item \code{\link[GenomicFeatures:TxDb-class]{TxDb}}: Annotation database
#'                from GenomicFeatures package
#'          \item \code{\link[ensembldb:EnsDb-class]{EnsDb}}: Annotation database
#'                from ensembldb package
#'        }
#' @param feature Annotation type. The available options depend on the input type:
#'        \itemize{
#'          \item For \code{GRanges}: Default is \code{"group"}. Can be any character
#'                string to describe the annotation type.
#'          \item For \code{TxDb}: One of \code{"gene"}, \code{"transcript"},
#'                \code{"exon"}, \code{"CDS"}, \code{"fiveUTR"}, \code{"threeUTR"},
#'                \code{"tRNAs"}, or \code{"geneModel"}. Default is \code{"gene"}.
#'                \itemize{
#'                  \item \code{"gene"}: Gene-level annotation (one range per gene)
#'                  \item \code{"transcript"}: Transcript-level annotation
#'                  \item \code{"exon"}: Individual exons with transcript and gene IDs
#'                  \item \code{"CDS"}: Coding sequences
#'                  \item \code{"fiveUTR"}: 5' untranslated regions
#'                  \item \code{"threeUTR"}: 3' untranslated regions
#'                  \item \code{"tRNAs"}: Transfer RNA annotations
#'                  \item \code{"geneModel"}: Combined gene model with CDS, UTRs, and
#'                        non-coding regions. Requires \code{OrganismDb} for gene symbols.
#'                }
#'          \item For \code{EnsDb}: One of \code{"gene"}, \code{"transcript"},
#'                \code{"exon"}, or \code{"disjointExons"}. Default is \code{"gene"}.
#'                \itemize{
#'                  \item \code{"gene"}: Gene-level annotation with gene_id and gene_name
#'                  \item \code{"transcript"}: Transcript-level annotation with tx_id,
#'                        gene_id, and gene_name
#'                  \item \code{"exon"}: Individual exons with exon_id, tx_id, gene_id,
#'                        and gene_name
#'                  \item \code{"disjointExons"}: Disjoint exonic parts (non-overlapping
#'                        exonic regions) linked to single genes
#'                }
#'        }
#' @param date A \code{\link{Date}} object specifying when the annotation was created.
#'        Defaults to \code{Sys.Date()} if not provided.
#' @param source A character string indicating where the annotation comes from.
#'        For \code{TxDb} and \code{EnsDb} objects, defaults to the object name
#'        (deparsed). For \code{GRanges}, can be provided via \code{...}.
#' @param mdata A data frame with columns "name" and "value" containing additional
#'        metadata. For \code{TxDb} and \code{EnsDb} objects, this is automatically
#'        extracted from the database metadata table if not provided.
#' @param OrganismDb An object of \code{\link[OrganismDbi]{OrganismDb}}. Only used
#'        when \code{feature = "geneModel"} for \code{TxDb} objects. It is used to
#'        extract gene symbols (SYMBOL) from the OrganismDb database and add them
#'        to the annotation metadata.
#' @param ... Additional parameters that can be passed when creating \code{annoGR}
#'        from \code{GRanges} objects: \code{feature}, \code{date}, \code{source},
#'        \code{mdata}
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form:
#' \itemize{
#'   \item \code{annoGR(gr)} where \code{gr} is a \code{GRanges} object
#'   \item \code{annoGR(txdb, feature = "gene")} where \code{txdb} is a \code{TxDb} object
#'   \item \code{annoGR(ensdb, feature = "gene")} where \code{ensdb} is an \code{EnsDb} object
#'   \item \code{new("annoGR", ...)} for direct instantiation (not recommended)
#' }
#' 
#' @slot seqnames,ranges,strand,elementMetadata,seqinfo Slots inherit from 
#'        \code{\link[GenomicRanges:GRanges-class]{GRanges}}. The ranges must have
#'        unique names. If names are missing when creating from \code{GRanges}, they
#'        will be automatically generated.
#' @slot source A character string indicating the annotation source (e.g., package
#'        name or custom identifier)
#' @slot date A \code{\link{Date}} object recording when the annotation object was
#'        created
#' @slot feature A character string indicating the annotation type. For \code{TxDb}:
#'        "gene", "exon", "transcript", "CDS", "fiveUTR", "threeUTR", "tRNAs", or
#'        "geneModel". For \code{EnsDb}: "gene", "exon", "transcript", or
#'        "disjointExons"
#' @slot mdata A data frame with columns "name" and "value" containing additional
#'        metadata from the annotation database (e.g., genome version, annotation
#'        version, etc.)
#' 
#' @section Validity:
#' The \code{annoGR} class has the following validity checks:
#' \itemize{
#'   \item The object must contain at least one range (cannot be empty)
#'   \item All ranges must have names (non-NULL)
#'   \item All range names must be unique (no duplicates)
#'   \item If \code{mdata} is provided, it must be a data frame with columns
#'         "name" and "value"
#' }
#' 
#' @section Methods:
#' \describe{
#'   \item{\code{annoGR(ranges, ...)}}{Generic function to create \code{annoGR}
#'         objects from various input types}
#'   \item{\code{info(object)}}{Display information about an \code{annoGR} object,
#'         including source, creation date, feature type, and metadata}
#'   \item{\code{as(object, "GRanges")}}{Coerce \code{annoGR} to \code{GRanges}}
#'   \item{\code{as(object, "annoGR")}}{Coerce \code{GRanges} to \code{annoGR}}
#' }
#' 
#' @author Jianhong Ou, Haibo Liu
#' @seealso \code{\link{annotatePeakInBatch}}, \code{\link{annoPeaks}},
#'          \code{\link[GenomicFeatures]{TxDb}}, \code{\link[ensembldb]{EnsDb}},
#'          \code{\link[GenomicRanges]{GRanges}}
#' @keywords classes
#' @exportClass annoGR
#' @import methods
#' @importFrom ensembldb EnsDb
#' @return Returns an \code{annoGR} object extending \code{GRanges} with additional
#'        slots for source, date, feature type, and metadata. The object can be used
#'        directly with \code{\link{annotatePeakInBatch}} and \code{\link{annoPeaks}}
#'        for peak annotation.
#' 
#' @examples
#' \dontrun{
#' ## Example 1: Create annoGR from EnsDb
#' library(EnsDb.Hsapiens.v79)
#' anno <- annoGR(EnsDb.Hsapiens.v79, feature="gene")
#' info(anno)  # Display annotation information
#' 
#' ## Example 2: Create annoGR from TxDb
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' genes <- annoGR(txdb, feature="gene")
#' transcripts <- annoGR(txdb, feature="transcript")
#' 
#' ## Example 3: Create annoGR from custom GRanges
#' custom_gr <- GRanges("chr1", IRanges(1000, 2000), strand="+")
#' names(custom_gr) <- "custom_feature_1"
#' anno <- annoGR(custom_gr, feature="custom", source="user_annotation")
#' 
#' ## Example 4: Gene model with gene symbols (requires OrganismDb)
#' library(Homo.sapiens)
#' geneModel <- annoGR(txdb, feature="geneModel", OrganismDb=Homo.sapiens)
#' 
#' ## Example 5: Disjoint exons from EnsDb (non-overlapping exonic parts)
#' disjoint_exons <- annoGR(EnsDb.Hsapiens.v79, feature="disjointExons")
#' }
#' 
setClass("annoGR", 
         representation(source = "character",
                        date = "Date",
                        feature = "character",
                        mdata = "data.frame"),
         contains = "GRanges",
         validity = function(object) {
             if (length(object@seqnames) < 1L) {
                 return("annotation is empty")
             }
             if (is.null(names(object@ranges))) {
                 return("annotation must have names")
             }
             if (any(duplicated(names(object@ranges)))) {
                 return("annotation names must be unique (duplicates found)")
             }
             if (!is.null(object@mdata)) {
                 if (!all(colnames(object@mdata) == c("name", "value"))) {
                     return("colnames of mdata must be 'name' and 'value'")
                 }
             }
             TRUE
         })

#' @importFrom S4Vectors Rle DataFrame 
#' @importFrom stats setNames
#' @importFrom GenomeInfoDb Seqinfo
newAnnoGR <- function (seqnames = Rle(), 
                       ranges = IRanges(), 
                       strand = Rle("*", length(seqnames)), 
                       mcols = DataFrame(), 
                       seqlengths = NULL, 
                       seqinfo = NULL,
                       ...) 
{
    if (!is(seqnames, "Rle")) 
        seqnames <- Rle(seqnames)
    if (!is.factor(runValue(seqnames))) 
        runValue(seqnames) <- factor(runValue(seqnames), 
                                     levels = unique(runValue(seqnames)))
    if (!is(ranges, "IRanges")) 
        ranges <- as(ranges, "IRanges")
    if (!is(strand, "Rle")) 
        strand <- Rle(strand)
    if (!is.factor(runValue(strand)) || !identical(levels(runValue(strand)), 
                                                   levels(strand()))) 
        runValue(strand) <- strand(runValue(strand))
    if (any(is.na(runValue(strand)))) {
        warning("missing values in strand converted to \"*\"")
        runValue(strand)[is.na(runValue(strand))] <- "*"
    }
    lx <- max(length(seqnames), length(ranges), length(strand))
    if (lx > 1) {
        if (length(seqnames) == 1) 
            seqnames <- rep(seqnames, lx)
        if (length(ranges) == 1) 
            ranges <- rep(ranges, lx)
        if (length(strand) == 1) 
            strand <- rep(strand, lx)
    }
    if (is.null(seqlengths)) {
        seqlengths <- setNames(rep(NA_integer_, length(levels(seqnames))), 
                               levels(seqnames))
    }
    if (is.null(seqinfo)) {
        seqinfo <- Seqinfo(names(seqlengths), seqlengths)
    }
    runValue(seqnames) <- factor(runValue(seqnames), seqnames(seqinfo))
    if (!is(mcols, "DataFrame")) 
        stop("'mcols' must be a DataFrame object")
    if (ncol(mcols) == 0L) {
        mcols <- mcols(ranges)
        if (is.null(mcols)) {
            rn <- names(ranges)
            removenames <- FALSE
            if (length(rn) != length(ranges)) {
                removenames <- TRUE
                rn <- seq_along(ranges)
            }
            mcols <- DataFrame(row.names = rn)
            if (removenames) {
                rownames(mcols) <- NULL
            }
        }
    }
    if (!is.null(mcols(ranges))) 
        mcols(ranges) <- NULL
    if (!is.null(rownames(mcols))) {
        if (is.null(names(ranges))) 
            names(ranges) <- rownames(mcols)
        rownames(mcols) <- NULL
    }
    new("annoGR", seqnames = seqnames, 
        ranges = ranges, strand = strand, 
        seqinfo = seqinfo, elementMetadata = mcols, ...)
}


newAGR <- function(gr, ...) {
    newAnnoGR(seqnames = seqnames(gr), 
              ranges = ranges(gr), 
              strand = strand(gr), 
              mcols = mcols(gr), 
              seqlengths = seqlengths(gr), 
              seqinfo = seqinfo(gr),
              ...)
}

setGeneric("annoGR", function(ranges, ...) standardGeneric("annoGR"))
setGeneric("info", function(object) standardGeneric("info"))

#' @name coerce
#' @import GenomicRanges
#' @rdname annoGR
#' @aliases coerce,GRanges,annoGR-method
#' coerce,annoGR,GRanges-method
#' @exportMethod coerce
#' @details
#' \strong{Coercion methods:}
#' 
#' \code{annoGR} objects can be coerced to and from \code{GRanges} objects:
#' \itemize{
#'   \item \code{as(annoGR_object, "GRanges")}: Converts \code{annoGR} to \code{GRanges},
#'         preserving all metadata columns but losing the \code{annoGR}-specific slots
#'         (source, date, feature, mdata)
#'   \item \code{as(GRanges_object, "annoGR")}: Converts \code{GRanges} to \code{annoGR}
#'         using default parameters (feature="group", date=Sys.Date())
#' }
#' 
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v79)
#' anno <- annoGR(EnsDb.Hsapiens.v79)
#' 
#' # Convert to GRanges
#' gr <- as(anno, "GRanges")
#' 
#' # Convert back to annoGR
#' anno2 <- as(gr, "annoGR")
#' }
setAs(from = "annoGR", to = "GRanges", function(from) {
    do.call(GRanges, args = append(list(seqnames = seqnames(from), 
                                        ranges = ranges(from),
                                        strand = strand(from),
                                        seqlengths = seqlengths(from),
                                        seqinfo = seqinfo(from)),
                                   as.list(mcols(from))))
})

setAs(from = "GRanges", to = "annoGR", function(from) {
    annoGR(from)
})

#' @rdname annoGR
#' @exportMethod info
#' @param object An \code{annoGR} object
#' @return The \code{info} method prints information about the \code{annoGR} object
#'        to the console, including class, source, creation date, feature type,
#'        and metadata, then returns \code{invisible(NULL)}.
#' @aliases info info,annoGR-method
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v79)
#' anno <- annoGR(EnsDb.Hsapiens.v79)
#' info(anno)  # Display annotation information
#' }
setMethod("info", "annoGR", function(object) {
    cat(class(object), "object;\n")
    cat("# source: ", object@source, "\n")
    cat("# create at: ", format(object@date, "%a %b %d %X %Y %Z"), "\n")
    cat("# feature: ", object@feature, "\n")
    mdata <- object@mdata
    for (i in seq_len(nrow(mdata))) {
        cat("# ", mdata[i, "name"], ": ", mdata[i, "value"],
            "\n", sep = "")
    }
})

#' @rdname annoGR
#' @exportMethod annoGR
#' @aliases annoGR,GRanges-method
#' @details
#' \strong{annoGR,GRanges-method:}
#' 
#' Converts a \code{GRanges} object to \code{annoGR}. This method is useful for
#' creating \code{annoGR} objects from custom annotation data or from GRanges
#' objects obtained from other sources.
#' 
#' \itemize{
#'   \item If \code{names(ranges)} is \code{NULL}, names will be automatically
#'         generated using zero-padded numeric indices
#'   \item The \code{feature} parameter defaults to \code{"group"} but can be
#'         set to any character string to describe the annotation type
#'   \item The \code{date} parameter defaults to \code{Sys.Date()} if not provided
#'   \item Additional parameters (\code{source}, \code{mdata}) can be provided
#'         via \code{...}
#' }
#' 
#' @examples
#' \dontrun{
#' # Create annoGR from custom GRanges
#' custom_gr <- GRanges("chr1", IRanges(1000, 2000), strand="+")
#' names(custom_gr) <- "custom_feature_1"
#' anno <- annoGR(custom_gr, feature="custom", source="user_annotation")
#' }
setMethod("annoGR", "GRanges", 
          function(ranges, feature = "group", date, ...) {
              if (missing("date")) {
                  date <- Sys.Date()
              }
              if (is.null(names(ranges))) {
                  names(ranges) <- make.names(
                      formatC(seq_along(ranges),
                              width = nchar(length(ranges)),
                              flag = "0"))
              }
              newAGR(gr = ranges,
                     date = date, 
                     feature = feature,
                     ...)
          })

#' @importFrom DBI dbGetQuery 
#' @importFrom BiocGenerics dbconn
#' @rdname annoGR
#' @exportMethod annoGR
#' @aliases annoGR,TxDb-method
#' @details
#' \strong{annoGR,TxDb-method:}
#' 
#' Converts a \code{TxDb} object to \code{annoGR} by extracting the specified
#' feature type. This method automatically extracts metadata from the TxDb
#' database and preserves information about the annotation source.
#' 
#' \strong{Feature types and their characteristics:}
#' \itemize{
#'   \item \code{"gene"}: Returns gene-level ranges with \code{gene_id} in metadata.
#'         Names are set to gene IDs.
#'   \item \code{"transcript"}: Returns transcript-level ranges with \code{tx_id},
#'         \code{tx_name}, and \code{gene_id} in metadata. Names are set to transcript IDs.
#'   \item \code{"exon"}: Returns individual exons with \code{exon_id}, \code{tx_name},
#'         and \code{gene_id} in metadata. Names are set to exon IDs.
#'   \item \code{"CDS"}: Returns coding sequences with \code{cds_id}, \code{tx_name},
#'         and \code{gene_id} in metadata. Names are set to CDS IDs.
#'   \item \code{"fiveUTR"}: Returns 5' UTR regions with \code{tx_name} in metadata.
#'         Multiple UTRs per transcript are unlisted.
#'   \item \code{"threeUTR"}: Returns 3' UTR regions with \code{tx_name} in metadata.
#'         Multiple UTRs per transcript are unlisted.
#'   \item \code{"tRNAs"}: Returns transfer RNA annotations using \code{tRNAs()}
#'         function from GenomicFeatures.
#'   \item \code{"geneModel"}: Returns a comprehensive gene model combining CDS,
#'         5' UTRs, 3' UTRs, and non-coding regions. Each element has \code{tx_name}
#'         and \code{feature_type} (CDS, 5UTR, 3UTR, or ncRNA) in metadata. If
#'         \code{OrganismDb} is provided, gene symbols are added to metadata.
#' }
#' 
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' 
#' # Gene-level annotation
#' genes <- annoGR(txdb, feature="gene")
#' 
#' # Transcript-level annotation
#' transcripts <- annoGR(txdb, feature="transcript")
#' 
#' # Gene model with gene symbols (requires OrganismDb)
#' library(Homo.sapiens)
#' geneModel <- annoGR(txdb, feature="geneModel", OrganismDb=Homo.sapiens)
#' }
setMethod("annoGR", "TxDb", 
          function(ranges, feature = c("gene", "transcript", "exon",
                                       "CDS", "fiveUTR", "threeUTR",
                                       "tRNAs", "geneModel"),
                   date, source, mdata, OrganismDb) {
              feature <- match.arg(feature)
              if (missing(mdata)) {
                  mdata <- dbGetQuery(dbconn(ranges), "select * from metadata")
              }
              if (missing(source)) {
                  source <- deparse(substitute(ranges, env = parent.frame()))
              }
              if (missing(date)) {
                  date <- Sys.Date()
              }
              gr <- if (!missing(OrganismDb)) {
                  TxDb2GR(ranges, feature, OrganismDb)
              } else {
                  TxDb2GR(ranges, feature)
              }
              newAGR(gr = gr, 
                     source = source,
                     date = date, 
                     feature = feature, 
                     mdata = mdata)
          })

#' @importFrom DBI dbGetQuery 
#' @importFrom BiocGenerics dbconn
#' @rdname annoGR
#' @exportMethod annoGR
#' @aliases annoGR,EnsDb-method
#' @details
#' \strong{annoGR,EnsDb-method:}
#' 
#' Converts an \code{EnsDb} object to \code{annoGR} by extracting the specified
#' feature type. This method automatically extracts metadata from the EnsDb
#' database and preserves information about the annotation source. The resulting
#' GRanges object is formatted to ensure consistent chromosome naming.
#' 
#' \strong{Feature types and their characteristics:}
#' \itemize{
#'   \item \code{"gene"}: Returns gene-level ranges with \code{gene_id} and
#'         \code{gene_name} in metadata. Names are set to gene IDs.
#'   \item \code{"transcript"}: Returns transcript-level ranges with \code{tx_id},
#'         \code{gene_id}, and \code{gene_name} in metadata. Names are set to
#'         transcript IDs.
#'   \item \code{"exon"}: Returns individual exons with \code{exon_id}, \code{tx_id},
#'         \code{gene_id}, and \code{gene_name} in metadata. Names are set to
#'         exon IDs.
#'   \item \code{"disjointExons"}: Returns disjoint exonic parts (non-overlapping
#'         exonic regions) using \code{exonicParts()} with
#'         \code{linked.to.single.gene.only=TRUE}. This creates non-overlapping
#'         exonic regions where each region is linked to a single gene. Useful for
#'         avoiding double-counting when analyzing exonic coverage.
#' }
#' 
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v79)
#' ensdb <- EnsDb.Hsapiens.v79
#' 
#' # Gene-level annotation
#' genes <- annoGR(ensdb, feature="gene")
#' 
#' # Transcript-level annotation
#' transcripts <- annoGR(ensdb, feature="transcript")
#' 
#' # Disjoint exons (non-overlapping exonic parts)
#' disjoint_exons <- annoGR(ensdb, feature="disjointExons")
#' }
setMethod("annoGR", "EnsDb",
          function(ranges, 
                   feature = c("gene", "transcript", "exon", "disjointExons"),
                   date, source, mdata) {
              feature <- match.arg(feature)
              if (missing(mdata)) {
                  mdata <- dbGetQuery(dbconn(ranges), "select * from metadata")
              }
              if (missing(source)) {
                  source <- deparse(substitute(ranges, env = parent.frame()))
              }
              if (missing(date)) {
                  date <- Sys.Date()
              }
              gr <- EnsDb2GR(ranges, feature)
              newAGR(gr = gr, 
                     source = source,
                     date = date, 
                     feature = feature, 
                     mdata = mdata)
          })
