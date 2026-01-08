#!/usr/bin/env Rscript
# Example script: Building a LOLA regionDB
#
# This demonstrates how to use the buildLOLAregionDB functions

# Load required libraries
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(ChIPpeakAnno)
})

# Source the build script
source("scripts/buildLOLAregionDB_Detailed.R")

# Example 1: Build a basic regionDB with ENCODE TFBS and UCSC data
cat("=== Example 1: Basic regionDB ===\n")
buildLOLAregionDB(
    genome = "hg19",
    output_dir = "LOLA_regionDB_basic",
    sources = c("ENCODE_TFBS", "ENCODE_DNase", "UCSC")
)

# Example 2: Add custom region sets
cat("\n=== Example 2: With custom regions ===\n")

# Create some example custom peaks
custom_peaks1 <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(1000000, 2000000, 3000000),
                    end = c(1000500, 2000500, 3000500)),
    name = c("Peak1", "Peak2", "Peak3")
)

custom_peaks2 <- GRanges(
    seqnames = "chr2",
    ranges = IRanges(start = c(500000, 1500000),
                    end = c(500500, 1500500)),
    name = c("PeakA", "PeakB")
)

# Build with custom regions
buildLOLAregionDB(
    genome = "hg19",
    output_dir = "LOLA_regionDB_custom",
    sources = c("ENCODE_TFBS", "Custom"),
    custom_regions = list(
        MyExperiment1 = custom_peaks1,
        MyExperiment2 = custom_peaks2
    ),
    custom_metadata = list(
        list(name = "MyExperiment1",
            description = "Custom peaks from experiment 1",
            source = "MyLab",
            cell_type = "HeLa",
            factor = "MyTF1"),
        list(name = "MyExperiment2",
            description = "Custom peaks from experiment 2",
            source = "MyLab",
            cell_type = "K562",
            factor = "MyTF2")
    )
)

cat("\n=== Examples complete ===\n")
cat("To use the regionDB:\n")
cat("  library(LOLA)\n")
cat("  regionDB <- loadRegionDB('LOLA_regionDB_basic')\n")

