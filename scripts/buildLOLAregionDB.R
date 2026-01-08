#!/usr/bin/env Rscript
# Script to build a LOLA regionDB from various public data sources
# 
# This script downloads and organizes genomic region sets from:
# - ENCODE transcription factor binding sites
# - ENCODE histone modifications
# - Roadmap Epigenomics data
# - DNase hypersensitive sites
# - Other public sources
#
# Usage:
#   Rscript buildLOLAregionDB.R --output-dir /path/to/regionDB --genome hg19

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(ChIPpeakAnno)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (any(grepl("--output-dir", args))) {
    args[which(grepl("--output-dir", args)) + 1]
} else {
    "LOLA_regionDB"
}
genome <- if (any(grepl("--genome", args))) {
    args[which(grepl("--genome", args)) + 1]
} else {
    "hg19"
}

cat("Building LOLA regionDB in:", output_dir, "\n")
cat("Genome assembly:", genome, "\n\n")

# Create output directory structure
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#' Download and process ENCODE transcription factor binding sites
#' 
#' @param genome Genome assembly (hg19, hg38, mm10, etc.)
#' @param output_dir Output directory for regionDB
#' @return List of processed GRanges objects
downloadENCODE_TFBS <- function(genome = "hg19", output_dir) {
    cat("Downloading ENCODE transcription factor binding sites...\n")
    
    collection_dir <- file.path(output_dir, "ENCODE_TFBS")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    # ENCODE TFBS URLs (these may need to be updated)
    base_url <- switch(genome,
                      hg19 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC",
                      hg38 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/encodeDCC",
                      mm10 = "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/encodeDCC",
                      stop("Unsupported genome: ", genome))
    
    # Download clustered TFBS data
    tfbs_url <- file.path(base_url, "wgEncodeRegTfbsClustered",
                         "wgEncodeRegTfbsClusteredV3.bed.gz")
    
    tryCatch({
        temp_file <- tempfile(fileext = ".bed.gz")
        download.file(tfbs_url, temp_file, quiet = TRUE)
        
        # Read and process
        data <- read.delim(gzfile(temp_file, "r"), header = FALSE,
                          stringsAsFactors = FALSE)
        unlink(temp_file)
        
        if (ncol(data) >= 4) {
            colnames(data)[1:4] <- c("seqnames", "start", "end", "TF")
            
            # Split by transcription factor
            tfs <- unique(data$TF)
            tf_granges <- lapply(tfs, function(tf) {
                tf_data <- data[data$TF == tf, ]
                gr <- GRanges(
                    seqnames = as.character(tf_data$seqnames),
                    ranges = IRanges(start = tf_data$start + 1,
                                    end = tf_data$end),
                    TF = tf
                )
                # Save as BED file
                bed_file <- file.path(regions_dir, paste0("ENCODE_TFBS_", 
                                                         gsub("[^A-Za-z0-9]", "_", tf),
                                                         ".bed"))
                export(gr, bed_file, format = "BED")
                return(gr)
            })
            names(tf_granges) <- tfs
            
            cat("  Downloaded", length(tf_granges), "transcription factor sets\n")
            return(tf_granges)
        }
    }, error = function(e) {
        warning("Failed to download ENCODE TFBS data: ", e$message)
        return(list())
    })
    
    return(list())
}

#' Download ENCODE histone modification data
#' 
#' @param genome Genome assembly
#' @param output_dir Output directory
#' @return List of processed GRanges objects
downloadENCODE_Histone <- function(genome = "hg19", output_dir) {
    cat("Downloading ENCODE histone modification data...\n")
    
    collection_dir <- file.path(output_dir, "ENCODE_Histone")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Common histone marks to download
    histone_marks <- c("H3K4me3", "H3K4me1", "H3K27ac", "H3K27me3",
                      "H3K36me3", "H3K9me3", "H3K9ac")
    
    # Note: This is a simplified example. In practice, you would need to:
    # 1. Query ENCODE API or download from specific URLs
    # 2. Process each mark for different cell types
    # 3. Organize by mark and cell type
    
    cat("  Note: ENCODE histone data requires API access or manual download\n")
    cat("  Please download from: https://www.encodeproject.org/\n")
    
    return(list())
}

#' Create collection metadata file
#' 
#' @param collection_dir Directory for the collection
#' @param collection_name Name of the collection
#' @param description Description of the collection
createCollectionMetadata <- function(collection_dir, collection_name, description) {
    metadata_file <- file.path(collection_dir, "collection.txt")
    metadata <- data.frame(
        collection = collection_name,
        description = description,
        genome = genome,
        stringsAsFactors = FALSE
    )
    write.table(metadata, metadata_file, sep = "\t", quote = FALSE,
               row.names = FALSE)
}

#' Create region set metadata
#' 
#' @param regions_dir Directory containing region sets
#' @param region_name Name of the region set
#' @param description Description
#' @param source Source of the data
#' @param cell_type Cell type (optional)
#' @param factor Factor/mark name (optional)
createRegionMetadata <- function(regions_dir, region_name, description,
                                source, cell_type = NA, factor = NA) {
    metadata_file <- file.path(regions_dir, paste0(region_name, "_description.txt"))
    metadata <- data.frame(
        name = region_name,
        description = description,
        source = source,
        cell_type = ifelse(is.na(cell_type), "", cell_type),
        factor = ifelse(is.na(factor), "", factor),
        genome = genome,
        stringsAsFactors = FALSE
    )
    write.table(metadata, metadata_file, sep = "\t", quote = FALSE,
               row.names = FALSE)
}

#' Download and process data from UCSC
#' 
#' @param genome Genome assembly
#' @param output_dir Output directory
downloadUCSC_Data <- function(genome = "hg19", output_dir) {
    cat("Downloading UCSC data...\n")
    
    collection_dir <- file.path(output_dir, "UCSC")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Example: Download CpG islands
    tryCatch({
        session <- browserSession("UCSC")
        genome(session) <- genome
        
        # CpG islands
        cpg_islands <- track(session, "cpgIslandExt")
        if (length(cpg_islands) > 0) {
            bed_file <- file.path(regions_dir, "UCSC_CpG_Islands.bed")
            export(cpg_islands, bed_file, format = "BED")
            createRegionMetadata(regions_dir, "UCSC_CpG_Islands",
                               "CpG islands from UCSC",
                               "UCSC", factor = "CpG_Islands")
            cat("  Downloaded CpG islands\n")
        }
    }, error = function(e) {
        warning("Failed to download UCSC data: ", e$message)
    })
    
    createCollectionMetadata(collection_dir, "UCSC",
                           "Genomic annotations from UCSC Genome Browser")
}

#' Main function to build regionDB
buildLOLAregionDB <- function(genome = "hg19", output_dir = "LOLA_regionDB") {
    cat("=== Building LOLA regionDB ===\n")
    cat("Genome:", genome, "\n")
    cat("Output directory:", output_dir, "\n\n")
    
    # Download ENCODE TFBS
    encode_tfbs <- downloadENCODE_TFBS(genome, output_dir)
    if (length(encode_tfbs) > 0) {
        collection_dir <- file.path(output_dir, "ENCODE_TFBS")
        createCollectionMetadata(collection_dir, "ENCODE_TFBS",
                               "ENCODE transcription factor binding sites")
    }
    
    # Download ENCODE Histone
    encode_histone <- downloadENCODE_Histone(genome, output_dir)
    if (length(encode_histone) > 0) {
        collection_dir <- file.path(output_dir, "ENCODE_Histone")
        createCollectionMetadata(collection_dir, "ENCODE_Histone",
                               "ENCODE histone modification data")
    }
    
    # Download UCSC data
    downloadUCSC_Data(genome, output_dir)
    
    cat("\n=== RegionDB structure created ===\n")
    cat("To use this database:\n")
    cat("  library(LOLA)\n")
    cat("  regionDB <- loadRegionDB('", output_dir, "')\n", sep = "")
    cat("\nNote: Some data sources require manual download or API access.\n")
    cat("Please refer to the documentation for complete data sources.\n")
}

# Run if executed as script
if (!interactive()) {
    buildLOLAregionDB(genome = genome, output_dir = output_dir)
}

