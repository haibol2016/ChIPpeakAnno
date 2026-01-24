# Building a LOLA RegionDB

This directory contains scripts to help build a LOLA (Locus Overlap Analysis) regionDB from various public data sources.

## Overview

LOLA requires a structured database of genomic region sets organized into collections. Each collection contains:
- BED files with genomic regions
- Metadata files describing the region sets

## Quick Start

### Basic Usage

```r
# Source the script
source("scripts/buildLOLAregionDB_Detailed.R")

# Build a basic regionDB
buildLOLAregionDB(
    genome = "hg19",
    output_dir = "LOLA_regionDB",
    sources = c("ENCODE_TFBS", "ENCODE_DNase", "UCSC")
)

# Load and use
library(LOLA)
regionDB <- loadRegionDB("LOLA_regionDB")
```

### Command Line Usage

```bash
Rscript scripts/buildLOLAregionDB.R --output-dir /path/to/regionDB --genome hg19
```

## Data Sources

### 1. ENCODE Transcription Factor Binding Sites

**Source**: ENCODE Project clustered TFBS data
**URL**: http://hgdownload.cse.ucsc.edu/goldenPath/[genome]/encodeDCC/wgEncodeRegTfbsClustered/

**What it includes**:
- Clustered transcription factor binding sites
- Organized by transcription factor
- Multiple cell types and conditions

**Usage**:
```r
buildLOLAregionDB(
    genome = "hg19",
    sources = c("ENCODE_TFBS")
)
```

### 2. ENCODE Histone Modifications

**Source**: ENCODE Project ChIP-seq data
**URL**: https://www.encodeproject.org/

**What it includes**:
- H3K4me3, H3K27ac, H3K27me3, H3K36me3, etc.
- Multiple cell types
- Broad and narrow peaks

**Note**: Requires ENCODE API access or manual download

**Manual download**:
1. Visit https://www.encodeproject.org/
2. Search for specific histone marks and cell types
3. Download BED files
4. Place in: `LOLA_regionDB/ENCODE_Histone/regions/`
5. Format: `ENCODE_Histone_[Mark]_[CellType].bed`

### 3. ENCODE DNase Hypersensitive Sites

**Source**: ENCODE Project DNase-seq data
**URL**: http://hgdownload.cse.ucsc.edu/goldenPath/[genome]/encodeDCC/wgEncodeRegDnaseClustered/

**What it includes**:
- Clustered DNase hypersensitive sites
- Cell type-specific accessibility

**Usage**:
```r
buildLOLAregionDB(
    genome = "hg19",
    sources = c("ENCODE_DNase")
)
```

### 4. Roadmap Epigenomics

**Source**: Roadmap Epigenomics Project
**URL**: http://www.roadmapepigenomics.org/
**Alternative**: https://egg2.wustl.edu/roadmap/web_portal/

**What it includes**:
- 15-state and 25-state chromatin state models
- Histone modification data across 127 cell types
- DNase hypersensitive sites
- DNA methylation data

**Manual download required**:
1. Visit the Roadmap portal
2. Download processed data for your genome
3. Convert to BED format
4. Organize by mark/cell type
5. Place in: `LOLA_regionDB/Roadmap_Epigenomics/regions/`

### 5. UCSC Genome Browser

**Source**: UCSC Genome Browser
**Available tracks**: CpG islands, repeats, genes, etc.

**Usage**:
```r
buildLOLAregionDB(
    genome = "hg19",
    sources = c("UCSC")
)
```

### 6. Custom Region Sets

Add your own region sets:

```r
# Prepare your data
my_peaks <- GRanges("chr1", IRanges(1000000, 2000000))
names(my_peaks) <- "MyExperiment"

# Build with custom regions
buildLOLAregionDB(
    genome = "hg19",
    output_dir = "LOLA_regionDB",
    sources = c("Custom"),
    custom_regions = list(MyExperiment = my_peaks),
    custom_metadata = list(
        list(name = "MyExperiment",
            description = "Peaks from my experiment",
            source = "MyLab",
            cell_type = "HeLa",
            factor = "MyTF")
    )
)
```

## Directory Structure

The script creates the following structure:

```
LOLA_regionDB/
├── ENCODE_TFBS/
│   ├── collection.txt
│   └── regions/
│       ├── ENCODE_TFBS_CTCF.bed
│       ├── ENCODE_TFBS_CTCF_description.txt
│       ├── ENCODE_TFBS_POLR2A.bed
│       └── ...
├── ENCODE_Histone/
│   ├── collection.txt
│   └── regions/
│       └── ...
├── ENCODE_DNase/
│   ├── collection.txt
│   └── regions/
│       └── ENCODE_DNase_Clustered.bed
├── Roadmap_Epigenomics/
│   ├── collection.txt
│   └── regions/
│       └── ...
└── UCSC/
    ├── collection.txt
    └── regions/
        └── UCSC_CpG_Islands.bed
```

## Metadata Format

### Collection Metadata (`collection.txt`)

```tsv
collection	description
ENCODE_TFBS	ENCODE transcription factor binding sites from clustered data (V3)
```

### Region Set Metadata (`*_description.txt`)

```tsv
filename	name	description	source	cell_type	factor
ENCODE_TFBS_CTCF.bed	ENCODE_TFBS_CTCF	ENCODE transcription factor binding sites for CTCF	ENCODE		CTCF
```

## Using the RegionDB

Once built, load and use:

```r
library(LOLA)

# Load the database
regionDB <- loadRegionDB("LOLA_regionDB")

# Define your query regions
query_regions <- GRanges("chr1", IRanges(1000000, 2000000))

# Define universe (all testable regions)
universe <- all_peaks_from_experiment

# Run enrichment analysis
results <- runLOLA(query_regions, universe, regionDB)

# View results
head(results)
```

## Notes and Limitations

1. **ENCODE API**: Some ENCODE data requires API access. See: https://www.encodeproject.org/help/api/

2. **Roadmap Data**: Roadmap data requires manual download and processing

3. **File Sizes**: Some datasets are very large. Consider disk space requirements

4. **Genome Assembly**: Ensure all data uses the same genome assembly

5. **Pre-built Database**: LOLA provides a pre-built core database that may be easier to use:
   ```r
   # Download LOLA Core
   # See: https://code.databio.org/LOLA/
   ```

## Troubleshooting

### Download Failures
- Check internet connection
- Verify URLs are still valid
- Some sources may require authentication

### UCSC Access Issues
- Install and configure `rtracklayer`
- May need to set up UCSC session properly

### Memory Issues
- Large datasets may require significant RAM
- Process in batches if needed

## References

- LOLA Package: https://bioconductor.org/packages/LOLA/
- ENCODE Project: https://www.encodeproject.org/
- Roadmap Epigenomics: http://www.roadmapepigenomics.org/
- UCSC Genome Browser: https://genome.ucsc.edu/

