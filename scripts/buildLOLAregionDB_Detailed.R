#!/usr/bin/env Rscript
# Detailed script to build a comprehensive LOLA regionDB
# 
# This script provides functions to download and organize data from:
# 1. ENCODE Project (TFBS, Histone marks, DNase)
# 2. Roadmap Epigenomics Project
# 3. UCSC Genome Browser annotations
# 4. Custom region sets
#
# Usage:
#   source("buildLOLAregionDB_Detailed.R")
#   buildLOLAregionDB(genome = "hg19", output_dir = "myLOLAregionDB")

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(ChIPpeakAnno)
    library(utils)
})

#' Build LOLA regionDB from multiple sources
#' 
#' @param genome Genome assembly (hg19, hg38, mm10, etc.)
#' @param output_dir Output directory for regionDB
#' @param sources Character vector of sources to include:
#'        "ENCODE_TFBS", "ENCODE_Histone", "ENCODE_DNase", 
#'        "Roadmap", "UCSC", "Custom"
#' @param custom_regions List of GRanges objects for custom regions
#' @param custom_metadata List of metadata data.frames for custom regions
buildLOLAregionDB <- function(genome = "hg19",
                              output_dir = "LOLA_regionDB",
                              sources = c("ENCODE_TFBS", "UCSC"),
                              custom_regions = NULL,
                              custom_metadata = NULL) {
    
    cat("=== Building LOLA regionDB ===\n")
    cat("Genome:", genome, "\n")
    cat("Output directory:", output_dir, "\n")
    cat("Sources:", paste(sources, collapse = ", "), "\n\n")
    
    # Create main output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Process each source
    if ("ENCODE_TFBS" %in% sources) {
        .downloadENCODE_TFBS(genome, output_dir)
    }
    
    if ("ENCODE_Histone" %in% sources) {
        .downloadENCODE_Histone(genome, output_dir)
    }
    
    if ("ENCODE_DNase" %in% sources) {
        .downloadENCODE_DNase(genome, output_dir)
    }
    
    if ("Roadmap" %in% sources) {
        .downloadRoadmap(genome, output_dir)
    }
    
    if ("UCSC" %in% sources) {
        .downloadUCSC(genome, output_dir)
    }
    
    if ("Custom" %in% sources && !is.null(custom_regions)) {
        .addCustomRegions(custom_regions, custom_metadata, output_dir)
    }
    
    cat("\n=== RegionDB build complete ===\n")
    cat("Database location:", output_dir, "\n")
    cat("\nTo use:\n")
    cat("  library(LOLA)\n")
    cat("  regionDB <- loadRegionDB('", normalizePath(output_dir), "')\n", sep = "")
}

#' Download ENCODE transcription factor binding sites
.downloadENCODE_TFBS <- function(genome, output_dir) {
    cat("Processing ENCODE transcription factor binding sites...\n")
    
    collection_dir <- file.path(output_dir, "ENCODE_TFBS")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Determine URL based on genome
    base_url <- switch(genome,
                      hg19 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered",
                      hg38 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/encodeDCC/wgEncodeRegTfbsClustered",
                      mm10 = "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/encodeDCC/wgEncodeRegTfbsClustered",
                      stop("Unsupported genome for ENCODE: ", genome))
    
    tfbs_file <- "wgEncodeRegTfbsClusteredV3.bed.gz"
    tfbs_url <- file.path(base_url, tfbs_file)
    
    tryCatch({
        cat("  Downloading from:", tfbs_url, "\n")
        temp_file <- tempfile(fileext = ".bed.gz")
        download.file(tfbs_url, temp_file, quiet = FALSE)
        
        # Read BED file
        data <- read.delim(gzfile(temp_file, "r"), header = FALSE,
                          stringsAsFactors = FALSE, comment.char = "#")
        unlink(temp_file)
        
        if (ncol(data) >= 4) {
            # BED is 0-based, GRanges is 1-based
            colnames(data)[1:4] <- c("seqnames", "start", "end", "TF")
            data$start <- data$start + 1  # Convert to 1-based
            
            # Split by transcription factor
            tfs <- unique(data$TF)
            cat("  Found", length(tfs), "transcription factors\n")
            
            processed <- 0
            for (tf in tfs) {
                tf_data <- data[data$TF == tf, , drop = FALSE]
                if (nrow(tf_data) > 0) {
                    gr <- GRanges(
                        seqnames = as.character(tf_data$seqnames),
                        ranges = IRanges(start = tf_data$start,
                                        end = tf_data$end),
                        TF = tf
                    )
                    
                    # Clean filename
                    safe_tf <- gsub("[^A-Za-z0-9_]", "_", tf)
                    bed_file <- file.path(regions_dir,
                                         paste0("ENCODE_TFBS_", safe_tf, ".bed"))
                    
                    # Export as BED (rtracklayer will convert back to 0-based)
                    export(gr, bed_file, format = "BED")
                    
                    # Create metadata
                    .createRegionMetadata(regions_dir,
                                         paste0("ENCODE_TFBS_", safe_tf),
                                         paste("ENCODE transcription factor",
                                              "binding sites for", tf),
                                         "ENCODE",
                                         factor = tf)
                    processed <- processed + 1
                }
            }
            cat("  Processed", processed, "transcription factor sets\n")
        }
        
        # Create collection metadata
        .createCollectionMetadata(collection_dir, "ENCODE_TFBS",
                                 paste("ENCODE transcription factor binding",
                                      "sites from clustered data (V3)"))
        
    }, error = function(e) {
        warning("Failed to download ENCODE TFBS: ", e$message, "\n")
        warning("You may need to download manually from ENCODE website\n")
    })
}

#' Download ENCODE histone modification data
.downloadENCODE_Histone <- function(genome, output_dir) {
    cat("Processing ENCODE histone modification data...\n")
    cat("  Note: ENCODE histone data requires API access\n")
    cat("  Please download from: https://www.encodeproject.org/\n")
    cat("  Or use the ENCODE API to query specific experiments\n\n")
    
    collection_dir <- file.path(output_dir, "ENCODE_Histone")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Common histone marks
    histone_marks <- c("H3K4me3", "H3K4me1", "H3K27ac", "H3K27me3",
                      "H3K36me3", "H3K9me3", "H3K9ac", "H3K4me2")
    
    cat("  Expected histone marks:", paste(histone_marks, collapse = ", "), "\n")
    cat("  To add data: Place BED files in:", regions_dir, "\n")
    cat("  Format: ENCODE_Histone_[Mark]_[CellType].bed\n\n")
    
    .createCollectionMetadata(collection_dir, "ENCODE_Histone",
                             "ENCODE histone modification ChIP-seq peaks")
}

#' Download ENCODE DNase hypersensitive sites
.downloadENCODE_DNase <- function(genome, output_dir) {
    cat("Processing ENCODE DNase hypersensitive sites...\n")
    
    collection_dir <- file.path(output_dir, "ENCODE_DNase")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    # DNase data URLs
    base_url <- switch(genome,
                      hg19 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered",
                      hg38 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/encodeDCC/wgEncodeRegDnaseClustered",
                      mm10 = "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/encodeDCC/wgEncodeRegDnaseClustered",
                      stop("Unsupported genome for ENCODE: ", genome))
    
    dnase_file <- "wgEncodeRegDnaseClustered.bed.gz"
    dnase_url <- file.path(base_url, dnase_file)
    
    tryCatch({
        cat("  Downloading from:", dnase_url, "\n")
        temp_file <- tempfile(fileext = ".bed.gz")
        download.file(dnase_url, temp_file, quiet = FALSE)
        
        data <- read.delim(gzfile(temp_file, "r"), header = FALSE,
                          stringsAsFactors = FALSE, comment.char = "#")
        unlink(temp_file)
        
        if (ncol(data) >= 3) {
            colnames(data)[1:3] <- c("seqnames", "start", "end")
            data$start <- data$start + 1
            
            gr <- GRanges(
                seqnames = as.character(data$seqnames),
                ranges = IRanges(start = data$start, end = data$end)
            )
            
            bed_file <- file.path(regions_dir, "ENCODE_DNase_Clustered.bed")
            export(gr, bed_file, format = "BED")
            
            .createRegionMetadata(regions_dir, "ENCODE_DNase_Clustered",
                                 "ENCODE clustered DNase hypersensitive sites",
                                 "ENCODE", factor = "DNase")
            
            cat("  Processed DNase hypersensitive sites\n")
        }
        
        .createCollectionMetadata(collection_dir, "ENCODE_DNase",
                                 "ENCODE DNase hypersensitive sites")
        
    }, error = function(e) {
        warning("Failed to download ENCODE DNase: ", e$message, "\n")
    })
}

#' Download Roadmap Epigenomics data
.downloadRoadmap <- function(genome, output_dir) {
    cat("Processing Roadmap Epigenomics data...\n")
    cat("  Note: Roadmap data requires manual download\n")
    cat("  Or use: https://egg2.wustl.edu/roadmap/web_portal/\n\n")
    
    collection_dir <- file.path(output_dir, "Roadmap_Epigenomics")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    cat("  Expected data structure:\n")
    cat("    - Histone marks (15-state model, 25-state model)\n")
    cat("    - DNase hypersensitive sites\n")
    cat("    - Chromatin states\n")
    cat("  Place processed BED files in:", regions_dir, "\n\n")
    
    .createCollectionMetadata(collection_dir, "Roadmap_Epigenomics",
                             "Roadmap Epigenomics Project data")
}

#' Download UCSC Genome Browser annotations
.downloadUCSC <- function(genome, output_dir) {
    cat("Processing UCSC Genome Browser annotations...\n")
    
    collection_dir <- file.path(output_dir, "UCSC")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    tryCatch({
        session <- browserSession("UCSC")
        genome(session) <- genome
        
        # CpG Islands
        if ("cpgIslandExt" %in% trackNames(session)) {
            cat("  Downloading CpG islands...\n")
            cpg <- track(session, "cpgIslandExt")
            bed_file <- file.path(regions_dir, "UCSC_CpG_Islands.bed")
            export(cpg, bed_file, format = "BED")
            .createRegionMetadata(regions_dir, "UCSC_CpG_Islands",
                                 "CpG islands from UCSC Genome Browser",
                                 "UCSC", factor = "CpG_Islands")
        }
        
        # Repeats (optional - can be large)
        if ("rmsk" %in% trackNames(session)) {
            cat("  Note: RepeatMasker data available but not downloaded\n")
            cat("  (too large for default download)\n")
        }
        
        .createCollectionMetadata(collection_dir, "UCSC",
                                 "Genomic annotations from UCSC Genome Browser")
        
    }, error = function(e) {
        warning("Failed to access UCSC: ", e$message, "\n")
        warning("You may need to install rtracklayer and configure UCSC access\n")
    })
}

#' Add custom region sets
.addCustomRegions <- function(custom_regions, custom_metadata, output_dir) {
    cat("Adding custom region sets...\n")
    
    collection_dir <- file.path(output_dir, "Custom")
    dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
    regions_dir <- file.path(collection_dir, "regions")
    dir.create(regions_dir, showWarnings = FALSE, recursive = TRUE)
    
    if (is.null(custom_regions)) {
        cat("  No custom regions provided\n")
        return()
    }
    
    if (is.list(custom_regions) && all(sapply(custom_regions, inherits, "GRanges"))) {
        for (i in seq_along(custom_regions)) {
            name <- names(custom_regions)[i]
            if (is.null(name) || name == "") {
                name <- paste0("Custom_RegionSet_", i)
            }
            
            gr <- custom_regions[[i]]
            bed_file <- file.path(regions_dir, paste0(name, ".bed"))
            export(gr, bed_file, format = "BED")
            
            # Use provided metadata or create default
            if (!is.null(custom_metadata) && i <= length(custom_metadata)) {
                meta <- custom_metadata[[i]]
            } else {
                meta <- list(name = name,
                           description = paste("Custom region set", i),
                           source = "Custom")
            }
            
            .createRegionMetadata(regions_dir, name,
                                 meta$description %||% paste("Custom set", i),
                                 meta$source %||% "Custom",
                                 meta$cell_type,
                                 meta$factor)
            
            cat("  Added:", name, "\n")
        }
    }
    
    .createCollectionMetadata(collection_dir, "Custom",
                             "User-provided custom region sets")
}

#' Helper function to create collection metadata
.createCollectionMetadata <- function(collection_dir, collection_name, description) {
    metadata_file <- file.path(collection_dir, "collection.txt")
    metadata <- data.frame(
        collection = collection_name,
        description = description,
        stringsAsFactors = FALSE
    )
    write.table(metadata, metadata_file, sep = "\t", quote = FALSE,
               row.names = FALSE, col.names = TRUE)
    cat("  Created collection metadata:", metadata_file, "\n")
}

#' Helper function to create region metadata
.createRegionMetadata <- function(regions_dir, region_name, description,
                                  source, cell_type = NA, factor = NA) {
    metadata_file <- file.path(regions_dir, paste0(region_name, "_description.txt"))
    metadata <- data.frame(
        filename = paste0(region_name, ".bed"),
        name = region_name,
        description = description,
        source = source,
        cell_type = ifelse(is.na(cell_type), "", as.character(cell_type)),
        factor = ifelse(is.na(factor), "", as.character(factor)),
        stringsAsFactors = FALSE
    )
    write.table(metadata, metadata_file, sep = "\t", quote = FALSE,
               row.names = FALSE, col.names = TRUE)
}

#' Helper operator for default values
`%||%` <- function(x, y) if (is.null(x)) y else x

# Example usage
if (FALSE) {
    # Example 1: Build basic regionDB
    buildLOLAregionDB(genome = "hg19",
                      output_dir = "LOLA_regionDB_hg19",
                      sources = c("ENCODE_TFBS", "ENCODE_DNase", "UCSC"))
    
    # Example 2: Add custom regions
    custom_peaks <- GRanges("chr1", IRanges(1000000, 2000000))
    names(custom_peaks) <- "MyCustomPeaks"
    
    buildLOLAregionDB(genome = "hg19",
                      output_dir = "LOLA_regionDB_custom",
                      sources = c("ENCODE_TFBS", "Custom"),
                      custom_regions = list(MyCustomPeaks = custom_peaks),
                      custom_metadata = list(
                          list(name = "MyCustomPeaks",
                              description = "My custom peak set",
                              source = "MyLab")
                      ))
}

