#' Extract gene and transcript information from Ensembl GTF file
#'
#' @description
#' Extracts gene identifiers, transcript identifiers, gene names, and biotype
#' information from an Ensembl GTF file. This function parses the GTF file and
#' extracts the key attributes: gene_id, transcript_id, exon_id, gene_name, gene_biotype,
#' and transcript_biotype.
#'
#' The function reads the GTF file using \code{rtracklayer::import()} and extracts
#' attributes from the metadata columns. For features that don't have certain
#' attributes (e.g., gene features don't have transcript_id), those fields will
#' be set to \code{NA}.
#'
#' @param gtf_file Character string specifying the path to the Ensembl GTF file.
#'   The file can be compressed (.gz) or uncompressed.
#' @param feature_type Optional character vector specifying which feature types
#'   to extract. Common Ensembl feature types include: "gene", "transcript",
#'   "exon", "CDS", "start_codon", "stop_codon", "UTR", "three_prime_utr",
#'   "five_prime_utr". If \code{NULL} (default), all features are returned.
#' @param unique_only Logical indicating whether to return only unique
#'   gene-transcript combinations. If \code{TRUE}, returns one row per unique
#'   combination of gene_id and transcript_id. Default is \code{FALSE}.
#'
#' @return A \code{data.frame} with the following columns:
#'   \itemize{
#' \item \code{seqnames}: Chromosome/contig name
#'     \item \code{start}: Start position (1-based)
#'     \item \code{end}: End position (1-based)
#'     \item \code{strand}: Strand ("+", "-", or "*")
#'     \item \code{gene_id}: Ensembl gene ID (e.g., "ENSG00000139618")
#'     \item \code{transcript_id}: Ensembl transcript ID (e.g., "ENST00000380152")
#'     \item \code{exon_id}: Ensembl exon ID (e.g., "ENST00000380152.1")
#'     \item \code{gene_name}: Gene symbol/name (e.g., "BRCA2")
#'     \item \code{gene_biotype}: Gene biotype (e.g., "protein_coding", "lncRNA")
#'     \item \code{transcript_biotype}: Transcript biotype (e.g., "protein_coding")
#'     \item \code{feature_type}: The GTF feature type (e.g., "gene", "transcript", "exon")
#'   }
#'
#'   If \code{unique_only = TRUE}, only the columns gene_id, transcript_id,
#'   gene_name, gene_biotype, and transcript_biotype are returned, with one row
#'   per unique gene-transcript combination.
#'
#' @details
#' \strong{GTF File Format:}
#' Ensembl GTF files contain annotation information with attributes in the 9th
#' column. This function extracts:
#' \itemize{
#'   \item \code{gene_id}: Always present for gene-related features
#'   \item \code{transcript_id}: Present for transcript, exon, CDS, and UTR features
#'   \item \code{gene_name}: Gene symbol (may be missing for some genes)
#'   \item \code{gene_biotype}: Classification of the gene (protein_coding, lncRNA, etc.)
#'   \item \code{transcript_biotype}: Classification of the transcript
#' }
#'
#' \strong{Feature Types:}
#' Different feature types contain different attributes:
#' \itemize{
#'   \item \code{gene}: Contains gene_id, gene_name, gene_biotype (no transcript_id)
#'   \item \code{transcript}: Contains gene_id, transcript_id, gene_name, gene_biotype, transcript_biotype
#'   \item \code{exon}, \code{CDS}, \code{UTR}: Contain all attributes
#' }
#'
#' \strong{Missing Values:}
#' \itemize{
#'   \item Features without transcript_id (e.g., gene features) will have \code{NA} for transcript_id
#'   \item Features without gene_name will have \code{NA} for gene_name
#'   \item If attributes are missing from the GTF file, they will be \code{NA} in the output
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges mcols seqnames start end strand
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract all information from a GTF file
#' gtf_info <- getGeneIdSymbolBiotypeFromGTF("Homo_sapiens.GRCh38.109.gtf")
#' head(gtf_info)
#'
#' # Extract only transcript features
#' transcript_info <- getGeneIdSymbolBiotypeFromGTF(
#'   "Homo_sapiens.GRCh38.109.gtf",
#'   feature_type = "transcript"
#' )
#'
#' # Get unique gene-transcript combinations
#' unique_genes_transcripts <- getGeneIdSymbolBiotypeFromGTF(
#'   "Homo_sapiens.GRCh38.109.gtf",
#'   unique_only = TRUE
#' )
#'
#' # Extract multiple feature types
#' gene_transcript_info <- getGeneIdSymbolBiotypeFromGTF(
#'   "Homo_sapiens.GRCh38.109.gtf",
#'   feature_type = c("gene", "transcript")
#' )
#' }
getGeneIdSymbolBiotypeFromGTF <- function(gtf_file,
                                          feature_type = NULL,
                                          unique_only = FALSE) {
    if (!file.exists(gtf_file)) {
        stop("GTF file '", gtf_file, "' does not exist!", call. = FALSE)
    }
    
    # Import GTF file
    message("Reading GTF file: ", gtf_file)
    gr <- rtracklayer::import(gtf_file, format = "GTF")
    
    # Filter by feature type if specified
    # If feature_type is NULL, all features are processed
    if (!is.null(feature_type)) {
        if (!all(feature_type %in% as.character(gr$type))) {
            available_types <- unique(as.character(gr$type))
            warning("Some feature types not found in GTF file. ",
                   "Available types: ", paste(available_types, collapse = ", "),
                   call. = FALSE)
        }
        gr <- gr[gr$type %in% feature_type]
    }
    # If feature_type is NULL, gr remains unchanged and all features are included
    
    # Extract attributes from metadata columns
    mcols_data <- GenomicRanges::mcols(gr)
    
    # Initialize result data frame
    result <- data.frame(
        seqnames = as.character(seqnames(gr)),
        start = start(gr),
        end = end(gr),
        strand = as.character(strand(gr)),
        gene_id = character(length(gr)),
        transcript_id = character(length(gr)),
        exon_id = character(length(gr)),
        gene_name = character(length(gr)),
        gene_biotype = character(length(gr)),
        transcript_biotype = character(length(gr)),
        feature_type = as.character(gr$type),
        
        stringsAsFactors = FALSE
    )
    
    # Extract gene_id
    if ("gene_id" %in% colnames(mcols_data)) {
        result$gene_id <- as.character(mcols_data$gene_id)
    } else if ("GeneID" %in% colnames(mcols_data)) {
        result$gene_id <- as.character(mcols_data$GeneID)
    } else {
        warning("gene_id not found in GTF file attributes", call. = FALSE)
        result$gene_id <- NA_character_
    }
    
    # Extract transcript_id
    if ("transcript_id" %in% colnames(mcols_data)) {
        result$transcript_id <- as.character(mcols_data$transcript_id)
    } else if ("TranscriptID" %in% colnames(mcols_data)) {
        result$transcript_id <- as.character(mcols_data$TranscriptID)
    } else {
        # For gene features, transcript_id will be NA
        result$transcript_id <- NA_character_
    }
    
    # Extract exon_id
    if ("exon_id" %in% colnames(mcols_data)) {
        result$exon_id <- as.character(mcols_data$exon_id)
    } else if ("ExonID" %in% colnames(mcols_data)) {
        result$exon_id <- as.character(mcols_data$ExonID)
    } else {
        result$exon_id <- NA_character_
    }
    # Extract gene_name
    if ("gene_name" %in% colnames(mcols_data)) {
        result$gene_name <- as.character(mcols_data$gene_name)
    } else if ("gene_symbol" %in% colnames(mcols_data)) {
        result$gene_name <- as.character(mcols_data$gene_symbol)
    } else if ("Name" %in% colnames(mcols_data)) {
        result$gene_name <- as.character(mcols_data$Name)
    } else {
        result$gene_name <- NA_character_
    }
    
    # Extract gene_biotype
    if ("gene_biotype" %in% colnames(mcols_data)) {
        result$gene_biotype <- as.character(mcols_data$gene_biotype)
    } else if ("gene_type" %in% colnames(mcols_data)) {
        result$gene_biotype <- as.character(mcols_data$gene_type)
    } else {
        result$gene_biotype <- NA_character_
    }
    
    # Extract transcript_biotype
    if ("transcript_biotype" %in% colnames(mcols_data)) {
        result$transcript_biotype <- as.character(mcols_data$transcript_biotype)
    } else if ("transcript_type" %in% colnames(mcols_data)) {
        result$transcript_biotype <- as.character(mcols_data$transcript_type)
    } else {
        result$transcript_biotype <- NA_character_
    }
    
    # Replace empty strings with NA
    result$gene_id[result$gene_id == ""] <- NA_character_
    result$transcript_id[result$transcript_id == ""] <- NA_character_
    result$exon_id[result$exon_id == ""] <- NA_character_
    result$gene_name[result$gene_name == ""] <- NA_character_
    result$gene_biotype[result$gene_biotype == ""] <- NA_character_
    result$transcript_biotype[result$transcript_biotype == ""] <- NA_character_
    
    # If unique_only, return only unique gene-transcript combinations
    if (unique_only) {
        # Keep only rows with both gene_id and transcript_id
        result <- result[!is.na(result$gene_id) & !is.na(result$transcript_id), ]
        
        if (nrow(result) > 0) {
            # Get unique combinations
            unique_cols <- c("gene_id", "transcript_id", "gene_name", 
                           "gene_biotype", "transcript_biotype")
            result <- result[!duplicated(result[, unique_cols, drop = FALSE]), 
                           unique_cols, drop = FALSE]
        } else {
            warning("No rows with both gene_id and transcript_id found. ",
                   "Returning empty data frame.", call. = FALSE)
            result <- result[, unique_cols, drop = FALSE]
        }
    }
    
    message("Extracted ", nrow(result), " features")
    if (unique_only && nrow(result) > 0) {
        message("Unique gene-transcript combinations: ", nrow(result))
    }
    
    return(result)
}

