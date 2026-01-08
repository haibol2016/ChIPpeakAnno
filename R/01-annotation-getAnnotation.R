#' Obtain genomic feature annotations from Ensembl via biomaRt
#' 
#' Retrieve transcription start sites (TSS), exons, UTRs, or transcripts for a
#' specified species from Ensembl using the \code{biomaRt} package. This function
#' queries biomaRt chromosome by chromosome and generates annotations that match
#' your genome assembly, which is the recommended approach for production
#' analyses.
#' 
#' @param mart A \code{\link[biomaRt]{Mart}} object from the \code{biomaRt} package.
#'   Must be a valid gene mart (e.g., "hsapiens_gene_ensembl" for human,
#'   "mmusculus_gene_ensembl" for mouse, "rnorvegicus_gene_ensembl" for rat).
#'   The mart must be a gene mart (dataset ending with "_gene_ensembl") as this
#'   function queries gene-related attributes such as \code{ensembl_gene_id},
#'   \code{chromosome_name}, \code{start_position}, and \code{end_position}.
#'   Create using \code{useMart(biomart = "ensembl", dataset = "species_gene_ensembl")}.
#'   See \code{\link[biomaRt]{useMart}} for details. This parameter is required
#'   and the function will stop if missing or not a valid gene mart object.
#' @param featureType A character string specifying the type of annotation to
#'   retrieve. Must be one of: "TSS" (default), "Exon", "5utr", "3utr",
#'   "ExonPlusUtr", or "transcript". Partial matching is supported via
#'   \code{match.arg()}.
#' 
#' @return A \code{\link[GenomicRanges]{GRanges}} object containing the requested
#'   genomic features. The structure and metadata columns depend on
#'   \code{featureType}:
#'   
#'   \strong{Common elements for all feature types:}
#'   \itemize{
#'     \item \code{seqnames}: Chromosome/contig names (from biomaRt
#'           "chromosome_name")
#'     \item \code{ranges}: \code{\link[IRanges]{IRanges}} with start and end
#'           positions
#'     \item \code{strand}: Strand information (1 for "+", -1 for "-")
#'     \item \code{names}: Feature identifiers (Ensembl IDs) stored in the
#'           IRanges names slot
#'     \item \code{description}: Feature descriptions in metadata columns
#'     \item \code{gene_biotype}: Gene biotype classification (e.g., 
#'           "protein_coding", "lncRNA", "pseudogene") from Ensembl
#'     \item \code{transcript_biotype}: Transcript biotype classification 
#'           (e.g., "protein_coding", "processed_transcript", "nonsense_mediated_decay")
#'           from Ensembl
#'   }
#'   
#'   \strong{Metadata columns by feature type:}
#'   \itemize{
#'     \item \code{"TSS"} (default): \code{description}, \code{gene_biotype}, and
#'           \code{transcript_biotype} metadata columns. The GRanges represents gene boundaries (start_position to
#'           end_position from biomaRt), with names set to
#'           \code{ensembl_gene_id}. For TSS annotation, use the start
#'           coordinate for positive strand genes and end coordinate for
#'           negative strand genes.
#'     \item \code{"Exon"}: \code{description}, \code{gene_biotype}, and
#'           \code{transcript_biotype} metadata columns. Names are set to \code{ensembl_exon_id}. Duplicated exon IDs are removed
#'           with a warning (except for "3utr", "5utr", "ExonPlusUtr").
#'     \item \code{"5utr"}: \code{ensembl_transcript_id}, \code{description},
#'           \code{gene_biotype}, and \code{transcript_biotype} metadata columns. Names are set to
#'           \code{ensembl_transcript_id} (made unique if needed).
#'     \item \code{"3utr"}: \code{ensembl_transcript_id}, \code{description},
#'           \code{gene_biotype}, and \code{transcript_biotype} metadata columns. Names are set to
#'           \code{ensembl_transcript_id} (made unique if needed).
#'     \item \code{"ExonPlusUtr"}: Multiple metadata columns:
#'           \code{ensembl_exon_id}, \code{ensembl_gene_id},
#'           \code{utr5start}, \code{utr5end}, \code{utr3start},
#'           \code{utr3end}, \code{description}, \code{gene_biotype}, and
#'           \code{transcript_biotype}. Names are set to
#'           \code{ensembl_exon_id} (made unique if needed).
#'     \item \code{"transcript"}: \code{ensembl_gene_id}, \code{description},
#'           \code{gene_biotype}, and \code{transcript_biotype} metadata
#'           columns. Names are set to \code{ensembl_transcript_id}.
#'   }
#'   
#'   \strong{Data processing:}
#'   \itemize{
#'     \item Rows with missing coordinates (NA in the third column) are removed
#'     \item Duplicate rows are removed using \code{unique()}
#'     \item Data is sorted by start position (third column)
#'     \item For feature types other than "3utr", "5utr", and "ExonPlusUtr",
#'           duplicated IDs (first column) are removed, keeping only the first
#'           occurrence. A warning is issued if duplicates are found.
#'     \item For "3utr", "5utr", and "ExonPlusUtr", duplicated IDs are allowed
#'           (multiple UTRs or exons per transcript/gene)
#'   }
#' 
#' @note 
#' \itemize{
#'   \item For \code{featureType = "TSS"}: The GRanges represents the full gene
#'         boundaries. The transcription start site is at the \code{start}
#'         coordinate for positive strand genes (strand = 1) and at the
#'         \code{end} coordinate for negative strand genes (strand = -1).
#'   \item The version of the annotation database must match the genome used
#'         for mapping because coordinates may differ between genome releases.
#'         For example, if you are using Mus_musculus.v103 for mapping, you
#'         should use the corresponding Ensembl version (e.g., Ensembl release
#'         103) when creating the Mart object.
#' }
#' @author Lihua Julie Zhu, Jianhong Ou, Kai Hu
#' @references Durinck S. et al. (2005) BioMart and Bioconductor: a powerful
#' link between biological biomarts and microarray data analysis.
#' Bioinformatics, 21, 3439-3440.
#' @keywords misc
#' @export
#' @import GenomicRanges
#' @importFrom biomaRt getBM
#' @examples
#' 
#' if (interactive() || Sys.getenv("USER")=="jou" )
#' {
#'   library(biomaRt)
#'   # Must use a gene mart (dataset ending with "_gene_ensembl")
#'   # To see all available marts, use \code{listMarts()}
#'   # To see all available datasets, use \code{listDatasets(useMart("ensembl"))}
#' 
#'   mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#'   Annotation <- getAnnotation(mart, featureType = "TSS")
#' }
getAnnotation <- function(mart, 
                          featureType=c("TSS", "Exon", 
                                        "5utr", "3utr", 
                                        "ExonPlusUtr", 
                                        "transcript")) {
    featureType <- match.arg(featureType)
    if (missing(mart) || !inherits(mart, "Mart"))
    {
        stop("'mart' must be a valid Mart object from biomaRt. ",
             "Create using useMart() with the appropriate dataset.")
    }
    
    seqnames <- getBM(c("chromosome_name"), mart = mart)
    seqnames <- as.character(seqnames[, "chromosome_name"])
    annotation_data <- lapply(seqnames, function(chr_name){
        attr <- switch(featureType, 
                    TSS=c("ensembl_gene_id", "chromosome_name", 
                            "start_position", "end_position", 
                            "strand", "description", "gene_biotype",
                            "transcript_biotype"),
                    Exon=c("ensembl_exon_id","chromosome_name", 
                            "exon_chrom_start", "exon_chrom_end", 
                            "strand", "description", "gene_biotype",
                            "transcript_biotype"),
                    "5utr"=c("ensembl_transcript_id", "chromosome_name", 
                                "5_utr_start", "5_utr_end", "strand", 
                                "description", "gene_biotype",
                                "transcript_biotype"),
                    "3utr"=c("ensembl_transcript_id", "chromosome_name", 
                                "3_utr_start", "3_utr_end", "strand", 
                                "description", "gene_biotype",
                                "transcript_biotype"),
                    ExonPlusUtr=c('ensembl_exon_id','chromosome_name', 
                                    'exon_chrom_start','exon_chrom_end', 
                                    'strand', 'ensembl_gene_id', 
                                    '5_utr_start', '5_utr_end',
                                    '3_utr_start','3_utr_end', 
                                    'description', 'gene_biotype', 
                                    'transcript_biotype'),
                    transcript=c("ensembl_transcript_id", 
                                    "chromosome_name", 
                                    "transcript_start", "transcript_end", 
                                    'strand', 'description', 
                                    'ensembl_gene_id', 'gene_biotype',
                                    'transcript_biotype')
        )
        getBM(attributes = attr, 
            filters = "chromosome_name",
            values = chr_name, 
            mart = mart)
    })
    annotation_data <- do.call(rbind, annotation_data)
    
    # Remove rows with missing coordinates
    annotation_data <- annotation_data[!is.na(annotation_data[, 3L]), ]
    annotation_data <- unique(annotation_data)
    annotation_data <- annotation_data[order(annotation_data[, 3L]), ]
    
    # Handle duplicated IDs
    duplicated_ids <- annotation_data[duplicated(annotation_data[, 1L]), 1L]
    annotation_data <- annotation_data[!duplicated(annotation_data[, 1L]), ]
    
    if (length(duplicated_ids) > 0L && 
        !featureType %in% c("3utr", "5utr", "ExonPlusUtr"))
    {
        warning("Duplicated IDs found. Only one entry per ID will be returned.\n",
                "Duplicated IDs: ", 
                paste(utils::head(duplicated_ids, 10), collapse = ", "),
                if(length(duplicated_ids) > 10) " ..." else "")
    }
    
    # Create GRanges object based on feature type
    if (featureType == "ExonPlusUtr") {
        result <- GRanges(
            seqnames = as.character(annotation_data[, 2L]),
            IRanges(
                start = as.numeric(as.character(annotation_data[, 3L])),
                end = as.numeric(as.character(annotation_data[, 4L])),
                names = make.names(as.character(annotation_data[, 1L]), 
                                   unique = TRUE)
            ),
            ensembl_exon_id = as.character(annotation_data[, 1L]),
            ensembl_gene_id = as.character(annotation_data[, 6L]),
            utr5start = as.numeric(as.character(annotation_data[, 7L])),
            utr5end = as.numeric(as.character(annotation_data[, 8L])),
            utr3start = as.numeric(as.character(annotation_data[, 9L])),
            utr3end = as.numeric(as.character(annotation_data[, 10L])),
            strand = annotation_data[, 5L],
            description = as.character(annotation_data[, 11L]),
            gene_biotype = as.character(annotation_data[, 12L]),
            transcript_biotype = as.character(annotation_data[, 13L])
        )
    } else if (featureType == "transcript") {
        result <- GRanges(
            seqnames = as.character(annotation_data[, 2L]),
            IRanges(
                start = as.numeric(annotation_data[, 3L]), 
                end = as.numeric(annotation_data[, 4L]), 
                names = as.character(annotation_data[, 1L])
            ), 
            strand = annotation_data[, 5L],   
            description = as.character(annotation_data[, 6L]), 
            ensembl_gene_id = as.character(annotation_data[, 7L]),
            gene_biotype = as.character(annotation_data[, 8L]),
            transcript_biotype = as.character(annotation_data[, 9L])
        )
    } else if (featureType %in% c("5utr", "3utr")) {
        result <- GRanges(
            seqnames = as.character(annotation_data[, 2L]),
            IRanges(
                start = as.numeric(annotation_data[, 3L]), 
                end = as.numeric(annotation_data[, 4L]), 
                names = make.names(as.character(annotation_data[, 1L]),
                                   unique = TRUE)
            ), 
            strand = annotation_data[, 5L],   
            description = as.character(annotation_data[, 6L]), 
            ensembl_transcript_id = as.character(annotation_data[, 1L]),
            gene_biotype = as.character(annotation_data[, 7L]),
            transcript_biotype = as.character(annotation_data[, 8L])
        )
    } else if (featureType == "TSS") {
        result <- GRanges(
            seqnames = as.character(annotation_data[, 2L]),
            IRanges(
                start = as.numeric(annotation_data[, 3L]),
                end = as.numeric(annotation_data[, 4L]), 
                names = as.character(annotation_data[, 1L])
            ), 
            strand = annotation_data[, 5L],   
            description = as.character(annotation_data[, 6L]),
            gene_biotype = as.character(annotation_data[, 7L]),
            transcript_biotype = as.character(annotation_data[, 8L])
        )
    } else { # featureType == "Exon"
        result <- GRanges(
            seqnames = as.character(annotation_data[, 2L]),
            IRanges(
                start = as.numeric(annotation_data[, 3L]),
                end = as.numeric(annotation_data[, 4L]), 
                names = as.character(annotation_data[, 1L])
            ), 
            strand = annotation_data[, 5L],   
            description = as.character(annotation_data[, 6L]),
            gene_biotype = as.character(annotation_data[, 7L]),
            transcript_biotype = as.character(annotation_data[, 8L])
        )
    }
    result
}
