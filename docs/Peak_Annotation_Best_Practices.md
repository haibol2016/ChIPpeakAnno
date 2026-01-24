# Best Practices for Peak (Genomic Region) Annotation

## Document Scope

This document provides comprehensive best practices for peak (genomic region) annotation, incorporating:

1. **Community-Wide Standards**: Best practices established by major consortia (ENCODE, modENCODE), published guidelines, and community consensus
2. **ChIPpeakAnno-Specific Guidance**: Implementation details and recommendations for using the ChIPpeakAnno package
3. **Factor-Specific Binding Patterns**: Recommendations for different factor types (TFs, histone marks, Pol II, RBPs, etc.) based on their distinct biological binding patterns (see Section 11)
4. **General Bioinformatics Principles**: Reproducibility, documentation, and quality control standards

**Note**: While this document uses ChIPpeakAnno for examples, the principles and standards outlined apply to peak annotation using any tool. Section 1 focuses specifically on community-wide consensus standards that are tool-agnostic.

**Important**: Different biological factors exhibit distinct binding patterns that require tailored annotation approaches. Section 11 provides factor-specific recommendations, and the [Peak Annotation Prioritization Strategies](Peak_Annotation_Prioritization_Strategies.md) document provides detailed prioritization strategies for custom annotation workflows.

## Significance of Peak/Region Annotation in Genomics

Peak annotation—the process of associating genomic regions (peaks) with known biological features such as genes, promoters, enhancers, and other regulatory elements—is a fundamental step in functional genomics analysis. Understanding the significance of this process helps researchers appreciate why proper annotation methodology is critical for biological interpretation.

### Why Peak Annotation Matters

**1. Biological Interpretation of Genomic Data**

Peak annotation transforms raw genomic coordinates into biologically meaningful information:
- **Gene association**: Identifies which genes are potentially regulated by transcription factors, histone modifications, or other chromatin features
- **Regulatory context**: Determines whether peaks are in promoters, enhancers, gene bodies, or intergenic regions
- **Functional prediction**: Enables prediction of biological functions through gene ontology and pathway enrichment analysis

**2. Integration with Other Omics Data**

Properly annotated peaks serve as a bridge between different types of genomic data:
- **RNA-seq integration**: Links binding sites to gene expression changes
- **Epigenomic data**: Connects histone modifications to gene regulation
- **3D chromatin structure**: Associates distal regulatory elements with target genes
- **Genetic variation**: Links disease-associated variants to regulatory regions

**3. Hypothesis Generation and Validation**

Annotation enables researchers to:
- **Generate testable hypotheses**: Identify candidate target genes for experimental validation
- **Prioritize experiments**: Focus on genes with strong binding signals or specific regulatory contexts
- **Validate findings**: Compare annotation results with known biological knowledge and published data

**4. Comparative and Meta-Analysis**

Standardized annotation allows:
- **Cross-study comparisons**: Compare results across different experiments and publications
- **Consortium-level analysis**: Aggregate data from large-scale projects (ENCODE, Roadmap Epigenomics)
- **Reproducibility**: Enable other researchers to reproduce and extend analyses

**5. Clinical and Translational Applications**

In biomedical research, annotation is essential for:
- **Disease mechanism discovery**: Identify regulatory disruptions in disease states
- **Biomarker identification**: Find regulatory regions associated with clinical outcomes
- **Therapeutic target discovery**: Identify genes and pathways for drug development
- **Personalized medicine**: Link genetic variants to regulatory changes

### Common Use Cases

**Transcription Factor Binding Analysis (ChIP-seq, CUT&RUN, CUT&TAG)**:
- Identify direct target genes of transcription factors
- Understand gene regulatory networks
- Discover tissue-specific or condition-specific binding patterns

**Histone Modification Profiling**:
- Characterize chromatin states (active promoters, enhancers, repressed regions)
- Study epigenetic changes during development or disease
- Identify regulatory element activity

**Open Chromatin Analysis (ATAC-seq, DNase-seq)**:
- Map accessible chromatin regions
- Identify active regulatory elements
- Study chromatin accessibility changes

**RNA-Binding Protein and microRNA Binding Localization**:
- **CLIP-seq (Crosslinking and Immunoprecipitation)**: Identify RNA-binding protein (RBP) binding sites on RNA transcripts
  - PAR-CLIP, HITS-CLIP, iCLIP: Different crosslinking methods for RBP binding site identification
  - Requires strand-aware annotation as binding is RNA-specific
  - Annotate to transcript-level features to identify target RNAs
- **RIP-seq (RNA Immunoprecipitation)**: Identify RNAs bound by specific proteins
  - Similar to CLIP-seq but with different crosslinking protocols
  - Transcript-level annotation essential for identifying target RNAs
- **microRNA target prediction**: Identify genomic regions where microRNAs bind
  - microRNA binding sites are typically in 3' UTRs of target genes
  - Requires annotation to 3' UTR regions or gene ends
  - Strand-aware annotation important for microRNA binding sites
- **RNA modification mapping**: Identify sites of RNA modifications (m6A, m5C, etc.)
  - MeRIP-seq, miCLIP: Methods for mapping RNA modifications
  - Transcript-level annotation required to identify modified transcripts
- **Key considerations**:
  - Use transcript-level annotation (`transcripts(TxDb)` or `transcripts(EnsDb)`)
  - Use `ignore.strand = FALSE` for strand-aware methods (CLIP-seq, RIP-seq)
  - Consider alternative splicing when annotating RBP binding sites

**Comparative Studies**:
- Compare binding patterns across conditions (e.g., treatment vs. control)
- Identify condition-specific regulatory changes
- Study evolutionary conservation of regulatory elements

### The Annotation Challenge

Despite its importance, peak annotation presents several challenges:

1. **Multiple possible associations**: A single peak may be near multiple genes or regulatory elements
2. **Distance interpretation**: Determining what constitutes a "functional" association
3. **Context dependency**: The same genomic region may have different functions in different cell types or conditions
4. **Annotation completeness**: Not all regulatory elements are well-annotated, especially enhancers and non-coding RNAs

These challenges underscore the importance of following best practices and using appropriate annotation strategies, as detailed in this document.

### Impact on Downstream Analysis

The quality and methodology of peak annotation directly affects:
- **Enrichment analysis results**: Incorrect annotations lead to misleading functional interpretations
- **Target gene identification**: Poor annotation misses true targets or includes false positives
- **Publication quality**: Incomplete or incorrect annotation undermines research conclusions
- **Reproducibility**: Inconsistent annotation methods prevent comparison across studies

**Best Practice**: Always validate annotation results against known biological knowledge and use multiple annotation approaches when possible to increase confidence in results.

---

## Getting Started

This section provides an overview for users beginning peak annotation analysis:

1. **Parameter Selection**: Consult the [Interactive Parameter Guide](Peak_Annotation_Parameter_Guide.md) to determine optimal parameter settings based on experimental data characteristics
2. **Workflow Implementation**: Refer to [Section 12: Complete Workflow Examples](#12-complete-workflow-examples) for comprehensive code examples
3. **Quality Validation**: Review [Section 7: Quality Control and Validation](#7-quality-control-and-validation) to ensure annotation quality and accuracy

**Common Parameter Settings**:
- **Standard annotation**: `output = "nearestLocation"`, `PeakLocForDistance = "middle"`, `FeatureLocForDistance = "TSS"`
- **Promoter analysis**: `output = "overlapping"`, `bindingRegion = c(-2000, 500)`
- **DNA-seq methods**: Always use `ignore.strand = TRUE`
- **RNA-seq methods**: Use `ignore.strand = FALSE` only if strandedness information is preserved by the sequencing protocol and maintained in the peak data

## Table of Contents

### Part I: Community Standards and Fundamentals
1. [Community Standards and Consensus](#1-community-standards-and-consensus)
2. [Pre-Annotation Considerations](#2-pre-annotation-considerations)

### Part II: Annotation Strategy and Methods
3. [Annotation Strategy Selection](#3-annotation-strategy-selection)
4. [Distance Calculation Methods](#4-distance-calculation-methods)
5. [Handling Overlapping Features](#5-handling-overlapping-features)
6. [Annotation Source Selection](#6-annotation-source-selection)

### Part III: Quality Control and Analysis
7. [Quality Control and Validation](#7-quality-control-and-validation)
8. [Statistical Considerations](#8-statistical-considerations)

### Part IV: Best Practices and Workflows
9. [Reproducibility and Documentation](#9-reproducibility-and-documentation)
10. [Common Pitfalls and Solutions](#10-common-pitfalls-and-solutions)
11. [Recommendations by Use Case](#11-recommendations-by-use-case)
12. [Complete Workflow Examples](#12-complete-workflow-examples)

**Related Documents**:
- **[Interactive Parameter Guide](Peak_Annotation_Parameter_Guide.md)**: Questionnaire-based guide to select optimal parameters for `annotatePeakInBatch()` based on binding factor type, sequencing method, peak width distribution, and metagene patterns
- **[Peak Annotation Prioritization Strategies](Peak_Annotation_Prioritization_Strategies.md)**: Comprehensive guide for factor-specific prioritization strategies, including detailed R code examples for custom prioritization functions for different factor types (TFs, histone marks, Pol II, RBPs, etc.)

---

## 1. Community Standards and Consensus

This section outlines best practices established by the broader research community, including major consortia (ENCODE, modENCODE), published guidelines, and community consensus standards.

### 1.1 Reference Annotation Standards

**Community Consensus**:
- **Use latest genome assemblies**: Always use the most recent, high-quality genome assemblies (e.g., GRCh38/hg38, GRCm39/mm39)
- **Prefer comprehensive annotations**: GENCODE, Ensembl, and RefSeq are widely accepted as gold standards
- **Version consistency**: Document annotation version numbers (e.g., GENCODE v44, Ensembl release 110)
- **Regular updates**: Re-annotate when new annotation versions are released to incorporate latest knowledge

**Recommended Sources** (in order of preference for human):
1. **GENCODE**: Most comprehensive, includes all transcript types, regularly updated
2. **Ensembl**: Comprehensive, good for comparative genomics
3. **RefSeq**: Curated, high-confidence annotations
4. **UCSC Known Genes**: Legacy support, widely compatible

**Best Practice**:
```r
# Use GENCODE or Ensembl for most comprehensive annotation
library(EnsDb.Hsapiens.v110)  # Latest Ensembl release
# OR
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # UCSC (more compatible)
```

### 1.2 Standard Distance Cutoffs

**Community Consensus on genomic region Definitions**:

| Region Type | Standard Definition | Common Variations | Source |
|------------|-------------------|-------------------|--------|
| **Proximal Promoter** | -2000 to +500 bp from TSS | -1000 to +100 bp, -2000 to +200 bp | ENCODE, most studies |
| **Core Promoter** | -500 to +100 bp from TSS | -250 to +50 bp | TATA box region |
| **Extended Promoter** | -5000 to +3000 bp from TSS | -3000 to +2000 bp | For broad factors |
| **5' UTR** | TSS to start codon | Variable | Gene-specific |
| **3' UTR** | Stop codon to polyA site | Variable | Gene-specific |
| **Enhancer** | >5 kb from TSS | >2 kb, >10 kb | Context-dependent |
| **Intergenic** | >5 kb from any gene | >2 kb, >10 kb | Context-dependent |

**ENCODE Standards**:
- Promoter region: **-2000 to +500 bp** from TSS (most commonly used)
- Enhancer region: **>5 kb** from TSS (intergenic or intronic)
- Gene body: From TSS to transcription end site (TES)

**Best Practice**:
```r
# Standard ENCODE promoter definition
promoters <- promoters(TxDb, upstream = 2000, downstream = 500)

# For transcription factors, commonly use:
bindingRegion <- c(-2000, 500)  # Standard promoter
# OR
bindingRegion <- c(-5000, 3000)  # Extended for broad factors
```

### 1.3 Quality Control Metrics and Thresholds

**Community Standards for ChIP-seq Quality**:

**FRiP (Fraction of Reads in Peaks)** - Method-Specific Standards:

| Method | Target | Excellent | Good | Acceptable | Poor | Source |
|--------|--------|-----------|------|------------|------|--------|
| **ChIP-seq** | TFs | >0.05 | 0.01-0.05 | 0.005-0.01 | <0.005 | ENCODE |
| **ChIP-seq** | Histones | >0.10 | 0.05-0.10 | 0.01-0.05 | <0.01 | ENCODE |
| **ATAC-seq** | Open Chromatin | >0.30 | 0.20-0.30 | 0.10-0.20 | <0.10 | Community |
| **Cut&RUN** | Any | >0.40 | 0.30-0.40 | 0.20-0.30 | <0.20 | Community |
| **Cut&Tag** | Any | >0.40 | 0.30-0.40 | 0.20-0.30 | <0.20 | Community |

**NSC (Normalized Strand Cross-correlation)**:
- **ChIP-seq**: >1.05 (excellent), 1.01-1.05 (good), 0.8-1.01 (acceptable), <0.8 (poor) | ENCODE
- **ATAC-seq/Cut&RUN/Cut&Tag**: >1.1 (excellent), higher thresholds expected

**RSC (Relative Strand Cross-correlation)**:
- **ChIP-seq**: >1 (excellent), 0.5-1 (good), 0.25-0.5 (acceptable), <0.25 (poor) | ENCODE
- **ATAC-seq/Cut&RUN/Cut&Tag**: >1.2 (excellent), higher thresholds expected

**Peak Number**:
- Context-dependent: 10K-100K (typical), 5K-10K or 100K-500K (acceptable), 1K-5K or >500K (questionable), <1K (poor)

**Peak Width**:
- Factor-dependent: 200-500 bp (typical for TFs), 100-200 or 500-1000 bp (acceptable), 50-100 or 1-5 kb (wide range), <50 or >5 kb (unusual)

**FRiP (Fraction of Reads in Peaks)**:
- **Definition**: Fraction of mapped reads that fall within called peaks
- **Method-specific considerations**:
  - **ChIP-seq**: Lower thresholds due to high background from cross-linking artifacts
    - ENCODE standard: FRiP > 0.01 (1%) for transcription factors
    - Histone marks: 0.01-0.10 typical
  - **ATAC-seq**: Higher thresholds expected (0.10-0.30+) due to low background
  - **Cut&RUN/Cut&Tag**: Very high thresholds expected (0.20-0.50+) due to extremely low background
- **Important**: Do not apply ChIP-seq thresholds to other methods, as they have different background characteristics

**NSC (Normalized Strand Cross-correlation)**:
- **Definition**: Ratio of fragment-length cross-correlation to background
- **ENCODE standard**: NSC > 1.05 indicates good quality

**RSC (Relative Strand Cross-correlation)**:
- **Definition**: Ratio of fragment-length peak to read-length peak
- **ENCODE standard**: RSC > 1 indicates good quality, RSC > 0.8 acceptable

### 1.4 Annotation Strategy Consensus

**Community Best Practices**:

1. **Use Multiple Annotation Approaches**:
   - Combine peak-centric (nearest gene) and feature-centric (promoter regions) methods
   - Validate annotations across multiple tools when possible
   - Cross-reference with expression data when available

2. **Standard Annotation Workflow**:
   ```
   Step 1: Annotate to nearest TSS (peak-centric)
   Step 2: Identify promoter peaks (feature-centric, -2kb to +500bp)
   Step 3: Classify genomic distribution (promoter, gene body, intergenic)
   Step 4: Integrate with expression data (if available)
   Step 5: Perform enrichment analysis
   ```

3. **Distance Calculation Consensus**:
   - **Peak reference**: Use peak center (middle) for robustness
   - **Feature reference**: Use TSS for gene annotation
   - **Rationale**: Less sensitive to peak width variation, biologically meaningful

4. **Handling Overlapping Features**:
   - **Precedence order**: Promoters > UTRs > Exons > Introns > Intergenic
   - **Multiple annotations**: Report all overlapping features when biologically relevant
   - **Best match**: Use shortest distance or feature precedence rules

### 1.5 Reproducibility Standards

**FAIR Data Principles** (Findable, Accessible, Interoperable, Reusable):

1. **Documentation Requirements**:
   - Genome assembly version and source
   - Annotation database version and release date
   - All parameters used in annotation (distance cutoffs, methods)
   - Tool versions and software dependencies
   - Date of annotation

2. **Standard Metadata**:
   ```r
   # Document in the analysis
   annotation_metadata <- list(
       genome = "GRCh38/hg38",
       annotation_source = "GENCODE v44",
       annotation_date = "2024-01-15",
       tool = "ChIPpeakAnno v3.45.2",
       parameters = list(
           output = "nearestLocation",
           PeakLocForDistance = "middle",
           FeatureLocForDistance = "TSS",
           bindingRegion = c(-2000, 500)
       )
   )
   ```

3. **Version Control**:
   - Use version control (Git) for annotation scripts
   - Tag releases with specific annotation versions
   - Document all software versions in sessionInfo()

4. **Sharing Standards**:
   - Provide annotated peak files in standard formats (BED, GFF)
   - Include metadata files with annotation parameters
   - Share annotation scripts and workflows

### 1.6 Multi-Omics Integration

**Community Best Practice**: Integrate annotation with other data types for biological validation.

**Recommended Integrations**:

1. **RNA-seq Expression Data**:
   - **For activators/transcription factors**: Prioritize annotation to expressed genes
   - **For repressors/suppressors**: Consider genes that are down-regulated or repressed
   - **For context-dependent factors** (can act as both activator and repressor):
     - **Do NOT use expression levels to weight annotation confidence** - suppressed genes will have low expression but binding is still functionally relevant
     - Consider all annotated genes regardless of expression level
     - Use differential expression analysis to identify both up- and down-regulated targets
     - Consider the direction of expression change rather than absolute expression level
   - **Differential expression**: Identify peaks associated with differentially expressed genes (both up- and down-regulated)
   - **Important considerations**:
     - **Expression weighting is problematic for**: Repressors, context-dependent factors, factors that maintain gene expression
     - **Expression weighting may be appropriate for**: Known activators where the focus is on active targets
     - Repressor binding sites may be associated with:
       - Lowly expressed or silenced genes
       - Genes that are down-regulated in your condition
       - Genes with repressive chromatin marks (H3K27me3, H3K9me3)
   
   **Example for repressors**:
   ```r
   # For repressor factors, consider focusing on:
   # 1. Lowly expressed genes
   lowly_expressed <- rownames(dds)[rowMeans(counts(dds)) < 10]
   
   # 2. Down-regulated genes
   res <- results(dds)
   downregulated <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
   
   # 3. Annotate and filter
   annotated <- annotatePeakInBatch(peaks, annoData)
   annotated_repressed <- annotated[annotated$feature %in% downregulated]
   
   # 4. Use appropriate background for enrichment
   # For repressors, consider using all genes or up-regulated genes as background
   enriched_go <- getEnrichedGO(
       annotated_repressed,
       orgAnn = "org.Hs.eg.db",
       universeGene = rownames(dds)  # All genes, not just expressed
   )
   ```
   
   **Example for context-dependent factors**:
   ```r
   # For factors that can activate or repress depending on context:
   # DO NOT filter by expression level - binding is relevant regardless
   
   # 1. Annotate all peaks (do not filter by expression)
   annotated <- annotatePeakInBatch(peaks, annoData)
   
   # 2. Identify both up- and down-regulated targets
   res <- results(dds)
   upregulated <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.05]
   downregulated <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
   
   # 3. Analyze both sets separately
   annotated_up <- annotated[annotated$feature %in% upregulated]
   annotated_down <- annotated[annotated$feature %in% downregulated]
   
   # 4. Enrichment analysis for both
   enriched_go_up <- getEnrichedGO(
       annotated_up,
       orgAnn = "org.Hs.eg.db",
       universeGene = rownames(dds)  # All genes as background
   )
   
   enriched_go_down <- getEnrichedGO(
       annotated_down,
       orgAnn = "org.Hs.eg.db",
       universeGene = rownames(dds)  # All genes as background
   )
   
   # 5. Compare functional categories between activated and repressed targets
   ```

2. **ATAC-seq Data**:
   - Overlap ChIP-seq peaks with accessible chromatin regions
   - Validate promoter/enhancer annotations with accessibility data

3. **Histone Modification Data**:
   - Use H3K4me3 to validate promoter annotations
   - Use H3K27ac to validate enhancer annotations
   - Use H3K4me1 to identify poised enhancers

4. **3D Chromatin Structure (Hi-C, ChIA-PET)**:
   - Link distal peaks to target genes via chromatin loops
   - Validate enhancer-gene associations

**Best Practice Examples**:

**For Activator/Transcription Factors**:
```r
# Integrate with expression data
library(DESeq2)
expressed_genes <- rownames(dds)[rowSums(counts(dds)) > 10]

# Annotate peaks
annotated <- annotatePeakInBatch(peaks, annoData)

# Filter to expressed genes
annotated_expressed <- annotated[annotated$feature %in% expressed_genes]

# Enrichment analysis with expressed genes as background
enriched_go <- getEnrichedGO(
    annotated_expressed,
    orgAnn = "org.Hs.eg.db",
    universeGene = expressed_genes  # Use expressed genes as background
)
```

**For Repressor/Suppressor Factors**:
```r
# For repressors, consider different strategies:

# Strategy 1: Focus on down-regulated genes
res <- results(dds)
downregulated <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]

# Strategy 2: Focus on lowly expressed genes
lowly_expressed <- rownames(dds)[rowMeans(counts(dds)) < 10]

# Strategy 3: Consider all genes (repressors may bind regardless of expression)
annotated <- annotatePeakInBatch(peaks, annoData)

# For enrichment, use all genes as background (not just expressed)
# This is important because repressors may target genes that are not expressed
enriched_go <- getEnrichedGO(
    annotated,
    orgAnn = "org.Hs.eg.db",
    universeGene = rownames(dds)  # All genes, not filtered by expression
)

# Additionally, check for association with repressive marks
# (H3K27me3, H3K9me3) if available
```

**Key Differences**:
- **Activators**: Typically associated with active, expressed genes → filter to expressed genes
- **Repressors**: May bind to silence genes → consider all genes or down-regulated genes
- **Context-dependent factors**: Can activate some genes and repress others → **do not filter by expression level**, analyze both up- and down-regulated targets separately
- **Expression weighting**: 
  - **Appropriate for**: Known activators where the focus is on active targets
  - **NOT appropriate for**: Repressors, context-dependent factors, factors maintaining expression
  - **Reason**: Suppressed genes have low expression but binding is still functionally relevant
- **Background for enrichment**: 
  - Activators: Use expressed genes as background
  - Repressors: Use all genes as background (repressors may target inactive genes)
  - Context-dependent: Use all genes as background, analyze activated and repressed targets separately

### 1.7 Community-Recommended Annotation Tools

**Tool Selection Consensus**:

| Use Case | Recommended Tools | Rationale |
|---------|------------------|-----------|
| **Standard R workflow** | ChIPpeakAnno, ChIPseeker | Comprehensive, well-documented |
| **Custom rules** | UROPA | Flexible, rule-based |
| **Large datasets** | UROPA, BEDTools | Efficient, parallel processing |
| **Web-based (no programming)** | GREAT, PAVIS | User-friendly, accessible |
| **Motif discovery** | HOMER | Specialized algorithms |
| **Multi-tool validation** | Use 2+ tools | Cross-validation improves confidence |

**Best Practice**: Use multiple tools for validation, especially for critical annotations.

### 1.8 Publication Standards

**What to Report in Publications**:

1. **Methods Section Should Include**:
   - Genome assembly version
   - Annotation database and version
   - Annotation tool and version
   - Distance cutoffs used (e.g., promoter = -2kb to +500bp)
   - Number of peaks annotated to each category
   - Annotation statistics (e.g., % in promoters, % intergenic)

2. **Supplementary Materials Should Include**:
   - Complete annotated peak files
   - Annotation scripts or workflows
   - Parameter settings used
   - Quality control metrics

3. **Figures Should Show**:
   - Genomic distribution pie charts
   - Distance to TSS distributions
   - Overlap with known regulatory elements
   - Enrichment analysis results

---

## 2. Pre-Annotation Considerations

### 2.1 Genome Version Matching
**Critical**: Ensure peak coordinates match the genome version of the annotation database.

**Best Practices**:
- **Always verify genome version compatibility** before annotation
- Use the same genome build for peak calling and annotation (e.g., hg19, hg38, mm10, mm39)
- If versions differ, use `rtracklayer::liftOver()` to convert coordinates
- Document the genome version used in the analysis

**Example**:
```r
# Check genome version of peaks
genome(myPeakList)

# Check genome version of annotation
genome(TxDb.Hsapiens.UCSC.hg19.knownGene)

# If mismatch, convert using liftOver
library(rtracklayer)
chain <- import.chain("hg18ToHg19.over.chain")
myPeakList_hg19 <- liftOver(myPeakList, chain)
```

### 2.2 Peak Quality Control
**Before annotation, assess peak quality**:

- **Filter low-quality peaks**: Remove peaks with low scores, fold changes, or p-values
- **Remove blacklisted regions**: Exclude peaks in problematic genomic regions (centromeres, telomeres, ENCODE blacklist)
- **Check peak widths**: Very narrow (<50bp) or very wide (>10kb) peaks may need special handling
- **Verify chromosome naming**: Ensure consistent chromosome naming (chr1 vs 1, chrM vs MT)

**Example**:
```r
# Filter peaks by score
peaks_filtered <- peaks[peaks$score > 10]

# Filter by width
peaks_filtered <- peaks[width(peaks) >= 50 & width(peaks) <= 10000]

# Remove blacklisted regions
library(AnnotationHub)
ah <- AnnotationHub()
blacklist <- ah[["AH5122"]]  # ENCODE blacklist
peaks_clean <- peaks[!overlapsAny(peaks, blacklist)]
```

### 2.3 Data Format Standardization
- **Use GRanges objects**: Convert all data to GRanges format for consistency
- **Standardize chromosome names**: Use `seqlevelsStyle()` to ensure consistent naming
- **Handle strand information**: Decide whether to use strand-aware or strand-agnostic annotation
- **Preserve metadata**: Keep important peak metadata (scores, p-values, etc.) during conversion

**Example**:
```r
# Convert to GRanges
peaks_gr <- toGRanges(peaks_bed, format="BED")

# Standardize chromosome names
library(GenomeInfoDb)
seqlevelsStyle(peaks_gr) <- "UCSC"  # or "NCBI", "Ensembl"
```

---

## 3. Annotation Strategy Selection

**Guidance**: For assistance in selecting the appropriate annotation strategy, consult the [Interactive Parameter Guide](Peak_Annotation_Parameter_Guide.md), which provides systematic guidance for selecting optimal approaches based on data characteristics.

### 3.1 Peak-Centric vs Feature-Centric Methods

#### Peak-Centric Method (Recommended for Most Cases)
**When to use**: When identifying the nearest or overlapping features for each peak is required.

**Options**:
- `output = "nearestLocation"` (default): Find the nearest feature based on distance calculation
- `output = "overlapping"`: Find all features that overlap with peaks
- `output = "both"`: Find both nearest and overlapping features
- `output = "shortestDistance"`: Find features with shortest distance from peak edges

**How `nearestLocation` works**:
1. **Reference point calculation**: Both peaks and features are reduced to single reference points:
   - Peak reference point: Based on `PeakLocForDistance` (start, middle, or end of peak)
   - Feature reference point: Based on `FeatureLocForDistance` (TSS, middle, start, end, or geneEnd)

2. **Distance calculation**: The distance is calculated as the absolute difference between these two reference points:
   - For plus strand features: `distance = PeakLoc - FeatureLoc`
   - For minus strand features: `distance = FeatureLoc - PeakLoc`
   - The sign indicates direction: negative = upstream, positive = downstream

3. **Overlapping features**: When a peak overlaps with a feature, the distance is **still calculated** as the difference between the reference points. For example:
   - If a peak (middle at position 1000) overlaps a gene with TSS at position 500
   - Distance = 1000 - 500 = 500 bp (positive = downstream of TSS)
   - The `insideFeature` column will indicate the overlap relationship ("inside", "overlapStart", "overlapEnd", etc.)
   - **Important**: Overlapping features CAN be selected as "nearest" if they have the smallest distance to the peak's reference point

4. **Selection**: The feature with the smallest absolute distance is selected as the nearest.

**Advantages**:
- Intuitive interpretation
- Each peak gets assigned to features
- Good for general annotation purposes
- Can select overlapping features if they are nearest

**Example**:
```r
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = annoData,
    output = "nearestLocation",
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS"
)
```

#### Feature-Centric Method
**When to use**: When identifying peaks in specific genomic contexts (e.g., promoter regions) is required.

**Options**:
- `output = "upstream"`: Peaks upstream of features
- `output = "downstream"`: Peaks downstream of features
- `output = "inside"`: Peaks completely within features
- `output = "upstream&inside"`: Peaks upstream or inside features

**Advantages**:
- Focuses on specific genomic contexts
- Useful for promoter/enhancer analysis
- Can use `bindingRegion` for precise region definition

**Example**:
```r
# Find peaks in promoter regions (2kb upstream, 500bp downstream of TSS)
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = annoData,
    output = "overlapping",
    FeatureLocForDistance = "TSS",
    bindingRegion = c(-2000, 500)
)
```

### 3.2 Choosing the Right Method

| Research Question | Recommended Method | Parameters |
|-----------------|-------------------|-------------|
| "What genes are near my peaks?" | Peak-centric | `output = "nearestLocation"` |
| "Which peaks are in promoters?" | Feature-centric | `output = "overlapping"`, `bindingRegion = c(-2000, 500)` |
| "Find all genes overlapping peaks" | Peak-centric | `output = "overlapping"` |
| "Identify bidirectional promoters" | Feature-centric | `output = "nearestBiDirectionalPromoters"` |
| "Peaks in gene bodies" | Feature-centric | `output = "inside"` or `bindingType = "fullRange"` |

---

## 4. Distance Calculation Methods

**Guidance**: For detailed guidance on selecting `PeakLocForDistance` and `FeatureLocForDistance` parameters, refer to [Step 3 and Step 7](Peak_Annotation_Parameter_Guide.md#step-3-analyze-peak-width-distribution) in the Parameter Guide.

### 4.1 Understanding Distance Calculation: Overlapping Features

**Key Question**: What is the distance when a peak overlaps with a feature?

**Answer**: The distance is **still calculated** as the difference between reference points, even when peaks overlap features.

**How it works**:

1. **Reference Points**: Both peaks and features are reduced to single coordinate points:
   - Peak point: Determined by `PeakLocForDistance` (e.g., peak middle at position 1000)
   - Feature point: Determined by `FeatureLocForDistance` (e.g., TSS at position 500)

2. **Distance Formula**:
   ```r
   # For plus strand features:
   distance = PeakLoc - FeatureLoc
   
   # For minus strand features:
   # distance = FeatureLoc - PeakLoc
   ```

3. **Example Scenarios**:

   **Scenario A: Peak overlaps gene, peak middle downstream of TSS**
   - Peak: chr1:950-1050 (middle = 1000)
   - Gene: chr1:500-2000 (TSS = 500, plus strand)
   - Distance = 1000 - 500 = **500 bp** (positive = downstream)
   - `insideFeature` = "inside" or "overlapStart"
   - **Result**: This overlapping feature can be selected as "nearest" if it has the smallest distance

   **Scenario B: Peak overlaps gene, peak middle upstream of TSS**
   - Peak: chr1:400-500 (middle = 450)
   - Gene: chr1:500-2000 (TSS = 500, plus strand)
   - Distance = 450 - 500 = **-50 bp** (negative = upstream)
   - `insideFeature` = "overlapStart"
   - **Result**: Negative distance indicates upstream position

   **Scenario C: Peak completely inside gene**
   - Peak: chr1:1000-1100 (middle = 1050)
   - Gene: chr1:500-2000 (TSS = 500, plus strand)
   - Distance = 1050 - 500 = **550 bp** (downstream of TSS)
   - `insideFeature` = "inside"
   - **Result**: Distance is still calculated from peak middle to TSS

4. **Important Points**:
   - **Distance is NOT zero** when overlapping (unless reference points coincide)
   - **Overlapping features CAN be nearest** if they have the smallest distance
   - The `insideFeature` column indicates the overlap relationship separately
   - The `shortestDistance` column shows the actual shortest distance from peak edges to feature edges

5. **Two Distance Metrics**:
   - `distancetoFeature`: Distance between reference points (used for "nearest" selection)
   - `shortestDistance`: Shortest distance from any peak edge to any feature edge (0 if overlapping)

**Visual Example**:
```
Gene:           |==========[TSS]==========|
               500       1000           2000
          
Peak A:      |--| (middle=450)
         400   500
         Distance = 450 - 500 = -50 bp (upstream, overlaps)
         
Peak B:                    |--| (middle=1050)
                        1000  1100
               Distance = 1050 - 500 = 550 bp (downstream, inside)
```

### 4.2 Peak Location for Distance (`PeakLocForDistance`)
    
**Options**:
- `"middle"` (recommended): Use the center of the peak
  - **Best for**: Most general cases, symmetric peaks
  - **Advantage**: Less sensitive to peak width variation
  
- `"start"` (default): Use the 5' end (smaller coordinate)
  - **Best for**: Directional peaks, when peak start is biologically meaningful
  
- `"end"`: Use the 3' end (larger coordinate)
  - **Best for**: When peak end is biologically meaningful
  
- `"endMinusStart"`: Use end for plus strand, start for minus strand
  - **Best for**: Stranded peaks where direction matters

**Recommendation**: Use `"middle"` for most cases as it's more robust to peak width variation.

### 4.3 Feature Location for Distance (`FeatureLocForDistance`)

**Options**:
- `"TSS"` (default): Transcription start site (5' end of gene)
  - **Best for**: Most gene annotation cases
  - **Definition**: Start for plus strand, end for minus strand
  
- `"geneEnd"`: Transcription end site (3' end of gene)
  - **Best for**: Analyzing 3' UTR or gene end regions
  
- `"start"`: Feature start (smaller coordinate)
  - **Best for**: Non-stranded features
  
- `"middle"`: Center of feature
  - **Best for**: When feature center is meaningful
  
- `"end"`: Feature end (larger coordinate)
  - **Best for**: Non-stranded features

**Recommendation**: Use `"TSS"` for gene annotation, as it's the most biologically relevant reference point.

### 4.4 Best Practice: Consistent Distance Calculation

**Recommended combination**:
```r
PeakLocForDistance = "middle"      # Center of peak
FeatureLocForDistance = "TSS"      # Transcription start site
```

**Why this combination**:
- Peak center is less affected by peak width
- TSS is the standard reference for gene annotation
- Consistent with most published analyses
- Biologically meaningful for most ChIP-seq experiments

---

## 5. Handling Overlapping Features

**Guidance**: For guidance on selecting the `select` parameter, refer to [Step 6](Peak_Annotation_Parameter_Guide.md#step-6-handle-multiple-overlapping-features) in the Parameter Guide.

### 5.1 Multiple Overlapping Features

**Problem**: A single peak may overlap with multiple features (e.g., overlapping genes, nested transcripts).

**Solutions**:

1. **Return all overlapping features** (`select = "all"`):
   - **Use when**: Complete information is required
   - **Advantage**: No information loss
   - **Disadvantage**: One peak → multiple rows

2. **Return first overlapping feature** (`select = "first"`):
   - **Use when**: One annotation per peak is required
   - **Advantage**: Straightforward one-to-one mapping
   - **Disadvantage**: May miss important features

3. **Return best match** (`select = "bestOne"` in `annoPeaks`):
   - **Use when**: The most relevant feature based on distance and overlap is required
   - **Advantage**: Prioritizes most relevant annotation using objective criteria
   - **Disadvantage**: Criteria for "best" may not match your specific biological needs
   
   **How "best" is determined**:
   The `select = "bestOne"` option uses a two-tier ranking system:
   
   1. **Primary criterion: Shortest distance** (`distanceToSite`)
      - Distance from peak to the binding site (TSS, gene end, etc., depending on `bindingType`)
      - Shorter distance = higher priority
   
   2. **Secondary criterion: Highest overlap** (Jaccard Index)
      - Calculated as: `width(intersection) / width(union)` of peak and annotation ranges
      - Higher overlap = better match
      - Range: 0 (no overlap) to 1 (complete overlap)
   
   **Selection process**:
   ```r
   # Internally, peaks are ordered by:
   order(peak_name, distanceToSite, -annoScore)
   # Then the first (best) match for each peak is kept
   ```
   
   **Example scenario**:
   - Peak overlaps two genes: Gene A (distance = 100 bp, Jaccard = 0.3) and Gene B (distance = 500 bp, Jaccard = 0.8)
   - Gene A will be selected because distance (100) < distance (500), even though Gene B has higher overlap
   - If distances are equal, the gene with higher Jaccard Index (overlap) is selected
   
   **When to use**:
   - When you need one annotation per peak and want the closest, most overlapping feature
   - When distance to binding site is your primary concern
   - When you want an objective, reproducible selection method
   
   **When NOT to use**:
   - When you need all overlapping features for downstream analysis
   - When biological relevance doesn't correlate with distance/overlap
   - When you want to apply custom selection criteria

**Example**:
```r
# Get all overlapping features
anno_all <- annotatePeakInBatch(peaks, annoData, 
                                output = "overlapping",
                                select = "all")

# Get only first overlapping feature
anno_first <- annotatePeakInBatch(peaks, annoData,
                                  output = "overlapping",
                                  select = "first")

# Get best match (shortest distance + highest overlap)
# Note: select = "bestOne" is available only in annoPeaks, not in annotatePeakInBatch
library(ensembldb)
library(EnsDb.Hsapiens.v75)
annoData <- annoGR(EnsDb.Hsapiens.v75)
anno_best <- annoPeaks(peaks, annoData,
                       bindingType = "startSite",
                       bindingRegion = c(-2000, 500),
                       select = "bestOne")  # Selects best based on distance + overlap
```

**Understanding the Jaccard Index (overlap score)**:
```r
# The Jaccard Index measures overlap between peak and annotation:
# Jaccard = width(intersection) / width(union)
# 
# Example:
# Peak:     chr1:1000-2000 (width = 1000)
# Gene A:   chr1:1500-3000 (width = 1500)
# Intersection: chr1:1500-2000 (width = 500)
# Union:        chr1:1000-3000 (width = 2000)
# Jaccard = 500 / 2000 = 0.25
#
# Higher Jaccard Index = better overlap = more likely to be selected
# (when distances are equal)
```

### 5.2 Feature Precedence Rules

**When using `assignChromosomeRegion()` or similar functions**, you can specify precedence to avoid double-counting:

**Common precedence order**:
```r
precedence <- c("Promoters", "fiveUTRs", "threeUTRs", 
                "Exons", "Introns", "immediateDownstream")
```

**Best Practice**: Define precedence based on biological importance:
1. Promoters (most important for TF binding)
2. UTRs (regulatory regions)
3. Exons (coding regions)
4. Introns (Note: first introns often contain enhancers and regulatory elements)
5. Downstream regions

**Important Note on Introns**:
- **First introns** are particularly important as they frequently contain:
  - Enhancers and regulatory elements
  - Transcription factor binding sites
  - Conserved non-coding sequences
- **Later introns** are less likely to contain functional regulatory elements
- **Consideration**: If your analysis focuses on enhancers or regulatory elements, you may want to:
  - Treat first introns separately from other introns
  - Give first introns higher precedence than later introns
  - Or use a more nuanced precedence: `c("Promoters", "fiveUTRs", "FirstIntrons", "Exons", "OtherIntrons", "threeUTRs")`

**Example**:
```r
distribution <- assignChromosomeRegion(
    peaks,
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
    precedence = c("Promoters", "fiveUTRs", "threeUTRs", 
                   "Exons", "Introns")
)
```

### 5.3 Handling Bidirectional Promoters

**Special case**: Peaks between two divergently transcribed genes.

**Best Practice**: Use dedicated function for bidirectional promoters:
```r
bdp_peaks <- peaksNearBDP(
    peaks,
    annoData,
    bindingRegion = c(-5000, 3000)
)
```

**Why this matters**: Bidirectional promoters are functionally distinct and should be identified separately.

---

## 6. Annotation Source Selection

**Guidance**: For guidance on choosing gene-level versus transcript-level annotation, refer to [Step 5](Peak_Annotation_Parameter_Guide.md#step-5-determine-annotation-level) in the Parameter Guide.

### 6.1 Community Standards for Annotation Sources

**Community Consensus** (from Section 1.1):
- **GENCODE**: Most comprehensive, includes all transcript types (recommended for human)
- **Ensembl**: Comprehensive, good for comparative genomics
- **RefSeq**: Curated, high-confidence annotations
- **UCSC Known Genes**: Legacy support, widely compatible

**Best Practice**: Use the latest version of comprehensive annotation sources (GENCODE or Ensembl) for most analyses.

### 6.2 TxDb vs EnsDb

**TxDb (UCSC-based)**:
- **Advantages**: 
  - Well-established, widely used
  - Good for standard analyses
  - Compatible with many tools
- **Disadvantages**:
  - May have fewer transcripts
  - Less frequently updated
- **Best for**: Standard gene annotation, compatibility

**EnsDb (Ensembl-based)**:
- **Advantages**:
  - More comprehensive transcript annotation
  - Regularly updated
  - Better for alternative splicing analysis
- **Disadvantages**:
  - May have more complexity
  - Some tools may not support
- **Best for**: Comprehensive annotation, alternative transcripts

**Recommendation**: Use EnsDb for comprehensive annotation, TxDb for compatibility.

### 6.3 Annotation Level: Gene vs Transcript

**Gene-level annotation**:
- **Use when**: 
  - You want one annotation per gene for simplicity
  - You're doing gene-level enrichment analysis
  - The gene has a single major transcript/promoter
- **Advantage**: Reduced complexity, minimal redundancy, one-to-one mapping
- **Disadvantage**: May miss transcript-specific regulatory information
- **Function**: `genes(TxDb)` or `genes(EnsDb)`

**Transcript-level annotation**:
- **Use when**: 
  - **Transcript isoforms have different promoters** (very common in mammals)
  - You need to identify which specific transcript/promoter is targeted
  - You're analyzing alternative promoter usage
  - The gene has multiple transcripts with distinct regulatory regions
- **Advantage**: 
  - Captures alternative transcripts and their unique promoters
  - More precise annotation for genes with multiple promoters
  - Can identify transcript-specific regulatory events
- **Disadvantage**: More complex, one peak may map to multiple transcripts
- **Function**: `transcripts(TxDb)` or `transcripts(EnsDb)`

**Important Biological Consideration**:
- **Many genes have multiple transcript isoforms with different promoters**
- In humans, ~50-70% of genes have alternative promoters
- Different promoters may be:
  - Used in different cell types
  - Regulated by different transcription factors
  - Located at different distances from the gene body
- **Implication**: Gene-level annotation may miss important regulatory specificity

**Best Practice Decision Tree**:
```r
# Decision: Gene-level vs Transcript-level

# Use GENE-LEVEL when:
# - Standard gene-level enrichment analysis
# - Gene has single major transcript
# - Initial exploratory analysis
# - Downstream analysis is gene-centric

# Use TRANSCRIPT-LEVEL when:
# - Analyzing promoter-specific regulation
# - Gene has multiple transcripts with different promoters
# - Need to identify alternative promoter usage
# - Cell-type or condition-specific transcript analysis
# - Precision is more important than simplicity
```

**Recommended Approach**:
1. **Start with transcript-level annotation** if you're analyzing promoter regions or transcription start sites
2. **Aggregate to gene-level** for downstream enrichment analysis if needed
3. **Use gene-level** only if you're certain the gene has a single major promoter or if simplicity is the priority

**Example: Transcript-level annotation for promoter analysis**:
```r
# For promoter-focused analysis, use transcript-level
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Get all transcripts (each may have different TSS/promoter)
transcripts_anno <- transcripts(TxDb)

# Annotate peaks to transcripts
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = transcripts_anno,
    output = "nearestLocation",
    FeatureLocForDistance = "TSS"
)

# If needed, aggregate to gene-level later
# (but you'll have lost transcript-specific information)
```

### 6.4 Custom Annotation Sources

**When to use custom annotations**:
- Specialized features (enhancers, TF binding sites, etc.)
- Literature-curated regions
- Experimental data (e.g., ATAC-seq peaks as annotation)

**Best Practice**:
```r
# Create custom annotation from GRanges
custom_anno <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(start = c(1000, 2000), end = c(1500, 2500)),
    strand = c("+", "-"),
    feature_id = c("enhancer1", "enhancer2")
)

# Use as annotation
annotated <- annotatePeakInBatch(peaks, AnnotationData = custom_anno)
```

---

## 7. Quality Control and Validation

**Guidance**: After annotation, utilize the [Validation Checklist](Peak_Annotation_Parameter_Guide.md#validation-checklist) in the Parameter Guide to verify annotation results.

### 7.1 Community Quality Standards

**Reference**: See Section 1.3 for ENCODE quality metrics and thresholds.

**Key Metrics**:

**FRiP Score (Fraction of Reads in Peaks)** - Method-Specific Thresholds:

| Method | Target Type | Excellent | Good | Acceptable | Poor | Notes |
|--------|------------|-----------|------|------------|------|-------|
| **ChIP-seq** | Transcription Factors | >0.05 | 0.01-0.05 | 0.005-0.01 | <0.005 | High background noise expected |
| **ChIP-seq** | Histone Marks | >0.10 | 0.05-0.10 | 0.01-0.05 | <0.01 | Broader peaks, higher signal |
| **ATAC-seq** | Open Chromatin | >0.30 | 0.20-0.30 | 0.10-0.20 | <0.10 | Low background, high signal |
| **Cut&RUN** | Any target | >0.40 | 0.30-0.40 | 0.20-0.30 | <0.20 | Very low background |
| **Cut&Tag** | Any target | >0.40 | 0.30-0.40 | 0.20-0.30 | <0.20 | Very low background |
| **CUT&RUN** | Histone marks | >0.50 | 0.40-0.50 | 0.30-0.40 | <0.30 | Even higher for histones |

**Important Notes**:
- **ChIP-seq**: Lower thresholds due to high background noise from cross-linking artifacts
- **ATAC-seq**: Higher thresholds expected due to low background and high signal-to-noise ratio
- **Cut&RUN/Cut&Tag**: Very high thresholds expected due to extremely low background (no cross-linking)
- **Method-specific considerations**:
  - ChIP-seq FRiP <0.01 may still be usable for TFs if other QC metrics are good
  - ATAC-seq FRiP <0.10 suggests poor library quality or sequencing issues
  - Cut&RUN/Cut&Tag FRiP <0.20 suggests technical problems

**NSC (Normalized Strand Cross-correlation)**:
- **ChIP-seq**: >1.05 indicates good quality
- **ATAC-seq**: >1.1 (higher expected due to better signal)
- **Cut&RUN/Cut&Tag**: >1.1 (higher expected)

**RSC (Relative Strand Cross-correlation)**:
- **ChIP-seq**: >1 indicates good quality, >0.8 acceptable
- **ATAC-seq**: >1.2 (higher expected)
- **Cut&RUN/Cut&Tag**: >1.2 (higher expected)

**Best Practice**: Always use method-specific thresholds. Do not apply ChIP-seq thresholds to ATAC-seq or Cut&RUN/Cut&Tag data.

### 7.2 Annotation Coverage

**Check what fraction of peaks were annotated**:
```r
# Calculate annotation rate
total_peaks <- length(peaks)
annotated_peaks <- sum(!is.na(annotated$feature))
annotation_rate <- annotated_peaks / total_peaks

# Expected rates (typical ChIP-seq):
# - Promoter-binding TFs: 40-60% in promoters
# - Enhancer-binding TFs: 20-40% in intergenic regions
# - Gene body-binding factors: 60-80% in gene bodies
```

**Red flags**:
- <10% annotation rate: Check genome version, annotation source
- >90% annotation rate: May indicate too permissive parameters

### 7.3 Distance Distribution

**Examine distance to features**:
```r
# Plot distance distribution
hist(annotated$distancetoFeature, breaks = 100,
     main = "Distance to Nearest TSS",
     xlab = "Distance (bp)")

# Check for expected patterns:
# - Promoter peaks: should cluster near 0
# - Enhancer peaks: may be further away
```

**Expected patterns**:
- Promoter-binding factors: Most peaks within ±5kb of TSS
- Enhancer-binding factors: More peaks >5kb from TSS
- Gene body factors: Peaks distributed across gene bodies

### 7.4 Genomic Element Distribution

**Validate distribution across genomic elements**:
```r
# Get distribution
distribution <- assignChromosomeRegion(
    peaks,
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene
)

# Visualize
pie1(distribution$percentage)
```
### 7.4.1 Peak Distribution Across Metagene

**Metagene plots** visualize the distribution of peaks relative to genomic features (typically TSS to TES) by averaging across many genes. This provides a global view of where peaks are enriched relative to gene structure.

**What metagene plots show**:
- **X-axis**: Distance from reference point (TSS, gene end, or gene middle)
- **Y-axis**: Peak density or count
- **Patterns reveal**:
  - Promoter-binding factors: Strong enrichment near TSS
  - Gene body-binding factors: Enrichment across gene body
  - Enhancer-binding factors: Enrichment in upstream/downstream regions
  - 3' end-binding factors: Enrichment near transcription end site (TES)

**Using `metagenePlot()` in ChIPpeakAnno**:

```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Basic metagene plot (TSS-centered)
metagenePlot(
    peaks,
    AnnotationData = TxDb,
    PeakLocForDistance = "middle",      # Use peak center
    FeatureLocForDistance = "TSS",      # Center on TSS
    upstream = 5000,                     # 5 kb upstream
    downstream = 5000                    # 5 kb downstream
)

# For full gene body analysis (TSS to TES + flanking regions)
metagenePlot(
    peaks,
    AnnotationData = TxDb,
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS",
    upstream = 5000,                     # 5 kb upstream of TSS
    downstream = 100000                  # Extend well past typical gene end
)

# Compare multiple peak sets
peaks_list <- GRangesList(
    "Condition1" = peaks_cond1,
    "Condition2" = peaks_cond2,
    "Control" = peaks_control
)
metagenePlot(peaks_list, AnnotationData = TxDb)
```

**Alternative: Gene end-centered plot**:
```r
# Center on transcription end site (TES)
metagenePlot(
    peaks,
    AnnotationData = TxDb,
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "geneEnd",   # Center on gene end
    upstream = 5000,                      # 5 kb upstream of TES
    downstream = 5000                     # 5 kb downstream of TES
)
```

**Interpreting metagene plots**:

**Expected patterns** (varies by factor type):

1. **Promoter-binding transcription factors**:
   - Strong peak near TSS (0 position)
   - Sharp enrichment within ±2 kb of TSS
   - Rapid drop-off in gene body
   - Example: Most sequence-specific TFs (e.g., STAT, NF-κB)

2. **Gene body-binding factors**:
   - Moderate enrichment at TSS
   - Sustained signal across gene body
   - Example: Elongation factors, some histone modifiers

3. **Enhancer-binding factors**:
   - Enrichment in upstream regions (>5 kb from TSS)
   - May also show enrichment downstream of gene end
   - Lower signal at TSS
   - Example: Distal enhancer-binding TFs

4. **3' end-binding factors**:
   - Enrichment near transcription end site (TES)
   - Example: Polyadenylation factors, 3' processing factors

5. **Broad histone marks**:
   - H3K4me3: Strong at TSS, extends into gene body
   - H3K36me3: Enriched in gene body, low at TSS
   - H3K27me3: Broad enrichment, may span large regions

**Best Practices for Metagene Plots**:

1. **Region selection**:
   - **TSS-centered**: Use for promoter analysis (upstream = 5 kb, downstream = 5-10 kb)
   - **Full gene**: Use for gene body analysis (upstream = 5 kb, downstream = 5-10 kb)
   - **Gene end-centered**: Use for 3' end analysis

2. **Normalization**:
   - Metagene plots show raw counts/density
   - For comparison across conditions, consider normalizing by:
     - Total number of peaks
     - Total number of reads
     - Number of genes analyzed

3. **Gene filtering**:
   - Consider filtering genes by:
     - Expression level (active vs. inactive genes)
     - Gene length (very short or very long genes may skew results)
     - Gene type (protein-coding vs. non-coding)

4. **Multiple samples**:
   - Use `GRangesList` to compare multiple conditions
   - Overlay plots or use faceting for comparison

**Example: Comprehensive metagene analysis**:
```r
# 1. TSS-centered plot (promoter analysis)
p1 <- metagenePlot(
    peaks,
    AnnotationData = TxDb,
    FeatureLocForDistance = "TSS",
    upstream = 5000,
    downstream = 5000
)

# 2. Full gene body plot
p2 <- metagenePlot(
    peaks,
    AnnotationData = TxDb,
    FeatureLocForDistance = "TSS",
    upstream = 5000,
    downstream = 100000  # Extend past typical gene end
)

# 3. Compare conditions
peaks_comparison <- GRangesList(
    "Treatment" = peaks_treatment,
    "Control" = peaks_control
)
p3 <- metagenePlot(
    peaks_comparison,
    AnnotationData = TxDb,
    upstream = 5000,
    downstream = 5000
)
```

**Common Issues and Solutions**:

1. **Flat/no signal**: 
   - Check if peaks are properly annotated
   - Verify genome version matches annotation
   - Ensure sufficient number of peaks

2. **Unexpected patterns**:
   - Verify `FeatureLocForDistance` matches your biological question
   - Check if gene filtering is appropriate
   - Consider if factor type matches expected pattern

3. **Comparison across conditions**:
   - Normalize by total peak count or read depth
   - Use same annotation and parameters for all conditions
   - Consider statistical testing for differences

**Integration with other analyses**:
- Combine with `genomicElementDistribution()` for categorical distribution
- Use with `binOverGene()` for detailed gene body analysis
- Compare with expression data to link binding to gene activity

### 7.5 Strand Consistency

**Check strand information**:
```r
# Verify strand handling
table(strand(peaks))
table(strand(annotated))

# For strand-aware analysis
annotated_strand <- annotatePeakInBatch(
    peaks,
    annoData,
    ignore.strand = FALSE  # Use strand information
)
```

**Best Practice**: 
- Use `ignore.strand = TRUE` (default) for most ChIP-seq data
- Use `ignore.strand = FALSE` only if peaks are truly stranded

**When peaks are truly stranded** (use `ignore.strand = FALSE`):
- **Important prerequisite**: Strandedness information must be preserved by the sequencing protocol and maintained in the peak data
- **Transcription machinery components**:
  - RNA polymerase II (Pol II) - binds to transcribed strand
  - Transcription elongation factors
  - Transcription initiation factors
- **RNA-binding proteins (RBPs)**:
  - Proteins that bind to nascent RNA transcripts
  - May show strand-specific binding patterns
- **Strand-specific sequencing methods** (RNA-based, not DNA-based):
  - GRO-seq (Global Run-On sequencing) - measures nascent transcription
  - PRO-seq (Precision Run-On sequencing)
  - NET-seq (Native Elongating Transcript sequencing)
  - These methods directly measure RNA transcription, so strand is meaningful
  - **Note**: Only use `ignore.strand = FALSE` if the sequencing protocol preserves strandedness and the peak calling/processing maintains this information

**When peaks are NOT stranded** (use `ignore.strand = TRUE`, default):
- **All DNA-binding ChIP-seq** (including strand-specific library protocols):
  - **Important**: Even with strand-specific library preparation, ChIP-seq measures DNA binding, not RNA
  - DNA-binding proteins bind to double-stranded DNA, so strand information is not biologically meaningful
  - Transcription factors (bind DNA, not strand-specific)
  - Histone modifications (on DNA, not strand-specific)
  - Chromatin remodelers (act on DNA, not strand-specific)
  - Most DNA-binding proteins
  - **Note**: Strand-specific library protocols in ChIP-seq provide sequencing direction information, but this doesn't indicate which strand the protein bound to
- **ATAC-seq**: Open chromatin regions (not strand-specific)
- **Standard ChIP-seq protocols**: Strand information is not meaningful for annotation

**Key Principle**: 
- **DNA-seq methods** (ChIP-seq, ATAC-seq): Strand information is NOT meaningful for standard annotation
  - Proteins bind to double-stranded DNA
  - Even strand-specific library prep doesn't tell you which strand was bound
  - **Standard annotation**: Always use `ignore.strand = TRUE` for DNA-based methods
  - **Exception**: Specialized ML-based binding pattern analysis may use strand as a feature, but this is for downstream pattern recognition, not annotation
- **RNA-seq methods** (GRO-seq, PRO-seq, NET-seq): Strand information IS meaningful **only if preserved**
  - These measure RNA transcription, which is inherently stranded
  - Use `ignore.strand = FALSE` for RNA-based methods **only if strandedness information is preserved by the sequencing protocol and maintained in the peak data**
  - Verify that your peak data contains valid strand information before using `ignore.strand = FALSE`

**Example: Strand-aware annotation for GRO-seq/PRO-seq**:
```r
# For GRO-seq or PRO-seq (RNA-based methods), use strand information
# These methods measure nascent RNA transcription, which is stranded
# IMPORTANT: Only use ignore.strand = FALSE if strandedness is preserved

# First, verify that strand information is present and meaningful
table(strand(groseq_peaks))  # Should show both + and - strands

# Then annotate with strand awareness
annotated_groseq <- annotatePeakInBatch(
    groseq_peaks,
    AnnotationData = annoData,
    output = "nearestLocation",
    ignore.strand = FALSE  # Only if strandedness is preserved by protocol
)

# This ensures peaks are only annotated to genes on the same strand
# GRO-seq signal on + strand only annotates to + strand genes
```

**Example: ChIP-seq (DNA-based, not stranded)**:
```r
# For ALL ChIP-seq data (including Pol II ChIP-seq), ignore strand
# Even Pol II ChIP-seq measures DNA binding, not RNA transcription
annotated_pol2_chip <- annotatePeakInBatch(
    pol2_chip_peaks,
    AnnotationData = annoData,
    output = "nearestLocation",
    ignore.strand = TRUE  # ChIP-seq is DNA-based, strand not meaningful
)

# For standard transcription factors
annotated_tf <- annotatePeakInBatch(
    tf_peaks,
    AnnotationData = annoData,
    output = "nearestLocation",
    ignore.strand = TRUE  # Default - all DNA-binding is not strand-specific
)
```

**Important Distinction**:
- **Pol II ChIP-seq**: Measures DNA binding → use `ignore.strand = TRUE` (default)
- **Pol II GRO-seq/PRO-seq**: Measures RNA transcription → use `ignore.strand = FALSE` **only if strandedness is preserved**
- The method matters, not just the target protein
- Always verify that strand information is meaningful in your data before using `ignore.strand = FALSE`

**Exception: Specialized Binding Pattern Analysis**:
- **Machine learning-based binding pattern analysis**: 
  - For specialized analyses using ML methods to detect binding patterns or motifs
  - Strand information might be used as a feature in pattern recognition
  - However, this is for downstream analysis, not for standard annotation
  - Standard annotation should still use `ignore.strand = TRUE` for DNA-seq
- **Note**: Even in ML approaches, strand information from ChIP-seq reflects sequencing direction, not actual binding strand
- **Best practice**: Use `ignore.strand = TRUE` for annotation, then use strand as a feature in ML models if needed for pattern analysis

---

## 8. Statistical Considerations

### 8.1 Multiple Testing Correction

**When performing enrichment analysis**, always correct for multiple testing:

```r
# GO enrichment with multiple testing correction
enriched_go <- getEnrichedGO(
    annotated,
    orgAnn = "org.Hs.eg.db",
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH"  # Benjamini-Hochberg correction
)
```

**Correction methods**:
- `"BH"` or `"fdr"`: Benjamini-Hochberg (recommended for most cases)
- `"bonferroni"`: Very conservative
- `"BY"`: Benjamini-Yekutieli (for dependent tests)

### 8.2 Permutation Testing for Overlaps

**When comparing peak sets**, use permutation testing:

```r
# Permutation test for overlap significance
perm_test <- peakPermTest(
    peaks1,
    peaks2,
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
    ntimes = 1000  # More permutations = more accurate
)
```

**Best Practice**: Use at least 1000 permutations for reliable p-values.

### 8.3 Background Selection

**For enrichment analysis**, choose appropriate background:

- **All genes in genome**: Standard, unbiased
- **Expressed genes**: More relevant for expression-related factors
- **Peak-associated genes**: For comparing different peak sets

**Example**:
```r
# Use all genes as background (default)
enriched_go <- getEnrichedGO(annotated, orgAnn = "org.Hs.eg.db")

# Use expressed genes as background
expressed_genes <- rownames(expression_data)[expression_data$expressed]
enriched_go <- getEnrichedGO(annotated, orgAnn = "org.Hs.eg.db",
                             universeGene = expressed_genes)
```

---

## 9. Reproducibility and Documentation

### 9.1 FAIR Data Principles

**Community Standard**: Follow FAIR (Findable, Accessible, Interoperable, Reusable) principles for all annotation data.

**Best Practices**:

1. **Findable**:
   - Use persistent identifiers (DOIs) for datasets
   - Include comprehensive metadata
   - Use standard file naming conventions

2. **Accessible**:
   - Store data in public repositories (GEO, ENCODE, Zenodo)
   - Provide clear access instructions
   - Use standard data formats (BED, GFF, BigWig)

3. **Interoperable**:
   - Use standard annotation formats
   - Follow community naming conventions
   - Provide format conversion tools if needed

4. **Reusable**:
   - Include detailed documentation
   - Provide analysis scripts
   - License data appropriately

### 9.2 Documentation Requirements

**Essential Information to Document**:

```r
# Create annotation metadata object
annotation_metadata <- list(
    # Data information
    peak_file = "peaks.bed",
    peak_calling_tool = "MACS2",
    peak_calling_version = "2.2.7.1",
    peak_calling_date = "2024-01-10",
    
    # Genome information
    genome_assembly = "GRCh38/hg38",
    genome_source = "UCSC",
    
    # Annotation information
    annotation_source = "GENCODE",
    annotation_version = "v44",
    annotation_release_date = "2023-10-11",
    annotation_date = "2024-01-15",
    
    # Tool information
    annotation_tool = "ChIPpeakAnno",
    tool_version = "3.45.2",
    R_version = R.version.string,
    
    # Parameters
    parameters = list(
        output = "nearestLocation",
        PeakLocForDistance = "middle",
        FeatureLocForDistance = "TSS",
        bindingRegion = c(-2000, 500),
        maxgap = 0,
        ignore.strand = TRUE
    ),
    
    # Results summary
    total_peaks = length(peaks),
    annotated_peaks = sum(!is.na(annotated$feature)),
    annotation_rate = sum(!is.na(annotated$feature)) / length(peaks)
)

# Save metadata
saveRDS(annotation_metadata, "annotation_metadata.rds")
```

### 9.3 Session Information

**Always capture session information**:

```r
# At the end of your analysis script
session_info <- sessionInfo()
capture.output(session_info, file = "sessionInfo.txt")

# Or use devtools
library(devtools)
session_info()
```

### 9.4 Workflow Documentation

**Create reproducible workflows**:

1. **Script Organization**:
   ```
   project/
   ├── scripts/
   │   ├── 01_peak_calling.R
   │   ├── 02_peak_annotation.R
   │   ├── 03_enrichment_analysis.R
   │   └── 04_visualization.R
   ├── data/
   │   ├── peaks/
   │   └── annotations/
   ├── results/
   │   └── annotated_peaks/
   └── README.md
   ```

2. **README Template**:
   ```markdown
   # Peak Annotation Workflow
   
   ## Overview
   Brief description of the analysis
   
   ## Requirements
   - R version: 4.3.0
   - Bioconductor: 3.18
   - Required packages: [list]
   
   ## Data
   - Peak files: [location]
   - Annotation: [source and version]
   - Genome: [assembly and version]
   
   ## Usage
   1. Run peak calling: `Rscript scripts/01_peak_calling.R`
   2. Run annotation: `Rscript scripts/02_peak_annotation.R`
   ...
   
   ## Parameters
   - Promoter definition: -2000 to +500 bp from TSS
   - Distance calculation: peak middle to TSS
   ...
   ```

### 9.5 Version Control

**Best Practices**:
- Use Git for all analysis scripts
- Tag releases with annotation versions
- Document all parameter changes in commit messages
- Use meaningful branch names (e.g., `annotation-v44-update`)

### 9.6 Sharing Annotated Data

**Standard Formats for Sharing**:

1. **BED Format** (with annotation columns):
   ```r
   # Export annotated peaks to BED
   annotated_bed <- as.data.frame(annotated)
   write.table(annotated_bed, "annotated_peaks.bed", 
               sep = "\t", quote = FALSE, row.names = FALSE)
   ```

2. **GFF Format**:
   ```r
   # Export to GFF
   rtracklayer::export(annotated, "annotated_peaks.gff", format = "GFF")
   ```

3. **Include Metadata File**:
   - JSON or YAML format
   - Include all annotation parameters
   - Include quality metrics

### 9.7 Reproducibility Checklist

**Before Publishing or Sharing**:

- [ ] All software versions documented
- [ ] All parameters documented
- [ ] Genome and annotation versions specified
- [ ] Analysis scripts provided
- [ ] Session information captured
- [ ] Metadata files included
- [ ] Results in standard formats
- [ ] README with clear instructions
- [ ] Code commented and organized
- [ ] Version control used

---

## 10. Common Pitfalls and Solutions

**Guidance**: If issues are encountered, consult the [Troubleshooting section](Peak_Annotation_Parameter_Guide.md#troubleshooting) in the Parameter Guide.

### 10.1 Genome Version Mismatch

**Problem**: Peaks and annotation use different genome versions.

**Symptoms**:
- Very low annotation rate
- Many peaks with `NA` features
- Warnings about sequence levels

**Solution**:
```r
# Check versions
genome(peaks)
genome(annoData)

# Convert if needed
library(rtracklayer)
chain <- import.chain("hg18ToHg19.over.chain")
peaks_converted <- liftOver(peaks, chain)
```

### 10.2 Chromosome Naming Inconsistency

**Problem**: Different chromosome naming conventions (chr1 vs 1, chrM vs MT).

**Symptoms**:
- Warnings about sequence levels
- Peaks not annotated

**Solution**:
```r
# Standardize chromosome names
library(GenomeInfoDb)
seqlevelsStyle(peaks) <- "UCSC"  # or "NCBI", "Ensembl"
seqlevelsStyle(annoData) <- "UCSC"
```

### 10.3 Too Many Unannotated Peaks

**Problem**: High fraction of peaks without annotation.

**Possible causes**:
1. Genome version mismatch
2. Annotation source too restrictive
3. Peaks in non-annotated regions (enhancers, etc.)

**Solutions**:
- Verify genome version
- Try different annotation sources (EnsDb vs TxDb)
- Use larger `maxgap` or `bindingRegion`
- Consider that some peaks may be truly intergenic

### 10.4 Double-Counting in Distribution Analysis

**Problem**: Peaks counted in multiple categories.

**Solution**: Use precedence rules:
```r
distribution <- assignChromosomeRegion(
    peaks,
    TxDb = TxDb,
    precedence = c("Promoters", "fiveUTRs", "threeUTRs", 
                   "Exons", "Introns")
)
```

### 10.5 Memory Issues with Large Datasets

**Problem**: Running out of memory with >100K peaks.

**Solutions**:
- ChIPpeakAnno automatically splits large datasets (>10K peaks)
- Process in batches manually if needed
- Use more efficient annotation sources
- Consider using command-line tools (UROPA) for very large datasets

---

## 11. Recommendations by Use Case

**Guidance**: For complete parameter sets for specific scenarios, refer to the [Parameter Combination Examples](Peak_Annotation_Parameter_Guide.md#parameter-combination-examples) in the Parameter Guide. For detailed factor-specific prioritization strategies, see [Peak Annotation Prioritization Strategies](Peak_Annotation_Prioritization_Strategies.md).

**Important Note**: Different biological factors exhibit distinct binding patterns that require tailored annotation approaches. This section provides factor-specific recommendations based on biological characteristics and binding patterns.

### 11.1 Pointed Binding Regulatory Factors

#### 11.1.1 Transcription Factors (TFs)

**Biological Characteristics:**
- **Peak width**: Method-dependent (ChIP-seq: 200-500 bp; Cut&RUN/Cut&Tag: 100-300 bp; ATAC-seq: 100-200 bp, but can be broad 500-2000 bp around highly transcribed genes)
- **Binding pattern**: Strong TSS/promoter enrichment, distance to regulatory site is critical
- **Relationship**: Often one-to-one relationship with target genes

**Recommended Annotation Strategy:**
- **Primary criterion**: Promoter proximity (within 2kb of TSS)
- **Method**: Use reference point-based internal logic (Method 2) - **do NOT provide `bindingRegion`**
- **Parameters**: `output = "overlapping"`, `maxgap = 2000`, `PeakLocForDistance = "middle"`, `FeatureLocForDistance = "TSS"`, `select = "all"` (for custom prioritization) or `select = "first"` (for simple annotation)

**Recommended approach**:
```r
# Step 1: Annotate with promoter-focused parameters
# Note: NOT providing bindingRegion triggers internal logic (Method 2)
anno_tf <- annotatePeakInBatch(
    peaks, 
    annoData,
    output = "overlapping",           # Find overlapping features
    maxgap = 2000,                     # Search within 2kb of features
    PeakLocForDistance = "middle",     # Use peak center for distance
    FeatureLocForDistance = "TSS",     # Use TSS as reference point
    select = "all"                     # Get all matches for prioritization
)

# Step 2: For simple annotation, use select = "first"
# For custom prioritization, see Peak_Annotation_Prioritization_Strategies.md

# Step 3: Distribution analysis
distribution <- assignChromosomeRegion(peaks, TxDb = TxDb)
# Expected: 40-60% in promoters for typical TFs
```

**Peak Width Considerations:**
- Check peak width distribution: `summary(width(peaks))`
- Narrow peaks (<200 bp): Use `PeakLocForDistance = "start"` or `"middle"`
- Medium peaks (200-500 bp): Use `PeakLocForDistance = "middle"`
- Broad peaks (>500 bp): Use `PeakLocForDistance = "middle"` (required)

#### 11.1.2 Promoter Histone Marks (H3K4me3, H3K4me2)

**Biological Characteristics:**
- **Peak width**: Method-dependent (ChIP-seq: 200-500 bp; Cut&RUN/Cut&Tag: 150-400 bp)
- **Binding pattern**: TSS-focused, but with some gene body signal
- **Note**: H3K4me3 is generally narrower than H3K4me2

**Recommended Annotation Strategy:**
- Similar to TFs but with broader distribution
- Prioritize promoter-proximal regions, but also consider gene body overlaps
- Use `maxgap = 2000` for promoter-focused analysis

**Recommended approach**:
```r
# Similar to TF annotation but may need larger maxgap
anno_h3k4me3 <- annotatePeakInBatch(
    peaks,
    annoData,
    output = "overlapping",
    maxgap = 2000,
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS",
    select = "all"
)

# Distribution analysis
distribution <- assignChromosomeRegion(peaks, TxDb = TxDb)
# Expected: High promoter enrichment, some gene body signal
```

#### 11.1.3 Enhancer Histone Marks (H3K4me1, H3K27ac)

**Biological Characteristics:**
- **Peak width**: Method-dependent (ChIP-seq: 200-1000 bp; Cut&RUN/Cut&Tag: 150-500 bp)
- **Binding pattern**: Intergenic and intronic enhancers, often far from gene bodies
- **Note**: H3K27ac marks active enhancers; H3K4me1 marks poised enhancers

**Recommended Annotation Strategy:**
- **Primary criterion**: Distance to nearest gene (tolerant of large distances)
- **Important**: Many enhancers are intergenic and far from gene bodies, so overlap is often 0
- Use larger `maxgap` (5-10kb) to capture distal enhancers
- Consider intergenic peaks as valid enhancer annotations

**Recommended approach**:
```r
# Use larger maxgap for enhancer analysis
anno_enhancer <- annotatePeakInBatch(
    peaks,
    annoData,
    output = "overlapping",
    maxgap = 10000,                    # Larger gap for distal enhancers
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS",
    select = "all"
)

# Many peaks will be intergenic (upstream/downstream)
intergenic <- anno_enhancer[anno_enhancer$insideFeature %in% 
                            c("upstream", "downstream")]
# This is expected and valid for enhancers
```

### 11.2 Gene Body Binding Factors

#### 11.2.1 Gene Body Histone Marks (H3K36me3, H3K27me3)

**Biological Characteristics:**
- **Peak width**: Method-dependent and mark-dependent
  - H3K36me3 (ChIP-seq): 500-2000 bp (medium to broad)
  - H3K27me3 (ChIP-seq): 2-10 kb (very broad domains)
  - Cut&RUN/Cut&Tag: Generally narrower (300-1000 bp for H3K36me3)
- **Binding pattern**: Gene body enrichment, may span multiple genes
- **Note**: Overlap quality matters more than distance for gene body marks

**Recommended Annotation Strategy:**
- **Primary criterion**: Gene body relationship (`inside` > `overlapStart` > `upstream`)
- **Secondary criterion**: Overlap quality (Jaccard index) - very important for gene body marks
- Use `output = "overlapping"` with `maxgap` to capture gene body overlaps
- Consider center-to-center distance for broad domains

**Recommended approach**:
```r
# Gene body-focused annotation
anno_gene_body <- annotatePeakInBatch(
    peaks,
    annoData,
    output = "overlapping",
    maxgap = 5000,                     # Larger gap for broad domains
    PeakLocForDistance = "middle",      # Required for broad peaks
    FeatureLocForDistance = "middle",  # Use gene center for distance
    select = "all"
)

# Distribution analysis
distribution <- assignChromosomeRegion(peaks, TxDb = TxDb)
# Expected: High gene body enrichment, especially for H3K36me3
# H3K27me3 may show broad domains spanning multiple genes
```

#### 11.2.2 RNA Polymerase II (Pol II)

**Biological Characteristics:**
- **Peak width**: Method-dependent (ChIP-seq: 200-1000 bp; GRO-seq/PRO-seq: strand-specific)
- **Binding pattern**: Primarily in gene body (elongation), but also at promoters (initiation)
- **Critical**: Strand matching is essential for RNA-based methods (GRO-seq, PRO-seq)

**Recommended Annotation Strategy:**
- **Primary criterion**: Gene body relationship (`inside` > `overlapStart`)
- **Critical for RNA-based methods**: Use `ignore.strand = FALSE` for strand-aware annotation
- Prioritize gene body over promoter for elongation marks
- For initiation marks, prioritize promoter regions

**Recommended approach**:
```r
# For ChIP-seq (not strand-specific)
anno_pol2 <- annotatePeakInBatch(
    peaks,
    annoData,
    output = "overlapping",
    maxgap = 5000,
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS",
    ignore.strand = TRUE,              # ChIP-seq is not strand-specific
    select = "all"
)

# For GRO-seq/PRO-seq (strand-specific)
anno_pol2_stranded <- annotatePeakInBatch(
    peaks,
    annoData,
    output = "overlapping",
    maxgap = 5000,
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS",
    ignore.strand = FALSE,             # Strand-aware for RNA-based methods
    select = "all"
)
```

#### 11.2.3 RNA-Binding Proteins (RBPs)

**Biological Characteristics:**
- **Peak width**: Method-dependent (CLIP-seq: 50-200 bp; RIP-seq: 100-500 bp)
- **Binding pattern**: Gene body enrichment, transcript-specific
- **Critical**: Requires transcript-level annotation and strand-aware analysis

**Recommended Annotation Strategy:**
- **Use transcript-level annotation**: `transcripts(TxDb)` or `transcripts(EnsDb)`
- **Strand-aware**: Use `ignore.strand = FALSE` for CLIP-seq, RIP-seq
- Prioritize gene body overlaps over promoter regions
- Consider alternative splicing when annotating RBP binding sites

**Recommended approach**:
```r
# Use transcript-level annotation for RBPs
tx_anno <- transcripts(TxDb)
anno_rbp <- annotatePeakInBatch(
    peaks,
    AnnotationData = tx_anno,
    output = "overlapping",
    maxgap = 2000,
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS",
    ignore.strand = FALSE,             # Strand-aware for RNA-based methods
    select = "all"
)

# Consider exon-level annotation for splicing factors
exon_anno <- exons(TxDb)
anno_exon <- annotatePeakInBatch(
    peaks,
    AnnotationData = exon_anno,
    output = "overlapping",
    ignore.strand = FALSE,
    select = "all"
)
```

### 11.3 Special Cases

#### 11.3.1 3' End Factors (PolyA factors, termination factors)

**Biological Characteristics:**
- **Binding pattern**: Gene end (3' UTR, polyA site) enrichment
- **Distance reference**: Use `FeatureLocForDistance = "geneEnd"`

**Recommended approach**:
```r
anno_3end <- annotatePeakInBatch(
    peaks,
    annoData,
    output = "overlapping",
    maxgap = 2000,
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "geneEnd",  # Use gene end as reference
    select = "all"
)
```

#### 11.3.2 Architectural Proteins (CTCF, Cohesin)

**Biological Characteristics:**
- **Peak width**: Method-dependent (ChIP-seq: 200-1000 bp)
- **Binding pattern**: Intergenic boundaries, domain boundaries, insulators
- **Note**: Often at topological domain boundaries, may be far from genes

**Recommended approach**:
```r
# Use large maxgap for architectural proteins
anno_ctcf <- annotatePeakInBatch(
    peaks,
    annoData,
    output = "overlapping",
    maxgap = 10000,                    # Large gap for boundary elements
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS",
    select = "all"
)

# Many peaks will be intergenic - this is expected
intergenic <- anno_ctcf[anno_ctcf$insideFeature %in% 
                        c("upstream", "downstream")]
```

#### 11.3.3 ATAC-seq Peak Annotation

**Biological Characteristics:**
- **Peak width**: Typically narrow (100-200 bp) but can be broad (500-2000 bp) around highly transcribed genes
- **Binding pattern**: Open chromatin regions (promoters, enhancers, gene bodies)
- **Note**: Peak width is method-dependent and can vary significantly

**Recommended approach**:
```r
# Check peak width first
summary(width(peaks))

# Similar to TF ChIP-seq but check for broad peaks
anno_atac <- annotatePeakInBatch(
    peaks,
    annoData,
    output = "overlapping",
    maxgap = 2000,
    PeakLocForDistance = "middle",     # Adjust based on peak width
    FeatureLocForDistance = "TSS",
    select = "all"
)

# Distribution analysis
distribution <- assignChromosomeRegion(peaks, TxDb = TxDb)
# Expected: High promoter and intergenic enrichment
```

### 11.4 Comparative Analysis (Multiple Conditions)

**Recommended approach**:
```r
# 1. Annotate each condition separately
annotated_cond1 <- annotatePeakInBatch(peaks_cond1, annoData)
annotated_cond2 <- annotatePeakInBatch(peaks_cond2, annoData)

# 2. Compare distributions
dist1 <- assignChromosomeRegion(peaks_cond1, TxDb = TxDb)
dist2 <- assignChromosomeRegion(peaks_cond2, TxDb = TxDb)

# 3. Statistical comparison
overlap_test <- findOverlappingPeaks(peaks_cond1, peaks_cond2)
perm_test <- peakPermTest(peaks_cond1, peaks_cond2, TxDb = TxDb)
```

### 11.5 Factor-Specific Parameter Summary

**Quick Reference Table**:

| Factor Type | Primary Binding | Recommended `maxgap` | `FeatureLocForDistance` | Priority Criterion |
|-------------|----------------|---------------------|----------------------|-------------------|
| **TFs** | Promoters | 2000 | TSS | Promoter proximity |
| **H3K4me3/H3K4me2** | Promoters | 2000 | TSS | Promoter proximity |
| **H3K4me1/H3K27ac** | Enhancers | 5000-10000 | TSS | Distance (tolerant) |
| **H3K36me3** | Gene body | 5000 | middle | Gene body + overlap |
| **H3K27me3** | Broad domains | 5000-10000 | middle | Gene body + overlap |
| **Pol II** | Gene body | 5000 | TSS | Gene body relationship |
| **RBPs** | Gene body | 2000 | TSS | Gene body (transcript-level) |
| **3' End factors** | Gene end | 2000 | geneEnd | Gene end proximity |
| **CTCF/Cohesin** | Boundaries | 10000 | TSS | Intergenic/boundary |
| **ATAC-seq** | Open chromatin | 2000 | TSS | Promoter/enhancer |

**Note**: All recommendations use reference point-based internal logic (Method 2) by **not providing `bindingRegion`**. For region-based annotation using `annoPeaks()`, see [Peak Annotation Prioritization Strategies](Peak_Annotation_Prioritization_Strategies.md) for detailed guidance.

**For detailed prioritization strategies and custom prioritization functions**, refer to [Peak Annotation Prioritization Strategies](Peak_Annotation_Prioritization_Strategies.md).

---

## 12. Complete Workflow Examples

**Guidance**: Before implementing these workflows, consult the [Interactive Parameter Guide](Peak_Annotation_Parameter_Guide.md) to determine optimal parameters for your specific data.

### 12.1 Standard ChIP-seq Annotation Workflow

```r
# Step 1: Load and prepare data
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

peaks <- toGRanges("peaks.bed", format = "BED")
TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Step 2: Prepare annotation
annoData <- genes(TxDb)

# Step 3: Annotate peaks
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = annoData,
    output = "nearestLocation",
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "TSS"
)

# Step 4: Add gene symbols
library(org.Hs.eg.db)
annotated <- addGeneIDs(
    annotated,
    feature_id_type = "entrez_id",
    IDs2Add = "symbol",
    orgAnn = org.Hs.eg.db
)

# Step 5: Distribution analysis
distribution <- assignChromosomeRegion(peaks, TxDb = TxDb)
pie1(distribution$percentage)

# Step 6: Enrichment analysis
enriched_go <- getEnrichedGO(
    annotated,
    orgAnn = "org.Hs.eg.db",
    pvalueCutoff = 0.01
)

# Step 7: Export results
write.csv(as.data.frame(annotated), "annotated_peaks.csv")
```

### 12.2 Promoter-Focused Analysis

```r
# Step 1: Prepare promoter annotation
promoters <- promoters(TxDb, upstream = 2000, downstream = 500)

# Step 2: Annotate to promoters
promoter_peaks <- annotatePeakInBatch(
    peaks,
    AnnotationData = promoters,
    output = "overlapping"
)

# Step 3: For peaks not in promoters, find nearest gene
non_promoter <- peaks[!names(peaks) %in% promoter_peaks$peak]
nearest_genes <- annotatePeakInBatch(
    non_promoter,
    AnnotationData = genes(TxDb),
    output = "nearestLocation"
)

# Step 4: Combine results
all_annotated <- c(promoter_peaks, nearest_genes)
```

### 12.3 Comprehensive Multi-Level Annotation

```r
# Step 1: Gene-level annotation
gene_anno <- annotatePeakInBatch(
    peaks,
    AnnotationData = genes(TxDb),
    output = "nearestLocation"
)

# Step 2: Transcript-level annotation
tx_anno <- annotatePeakInBatch(
    peaks,
    AnnotationData = transcripts(TxDb),
    output = "overlapping",
    select = "first"
)

# Step 3: Exon-level annotation
exon_anno <- annotatePeakInBatch(
    peaks,
    AnnotationData = exons(TxDb),
    output = "overlapping",
    select = "first"
)

# Step 4: Combine information
final_anno <- gene_anno
mcols(final_anno)$transcript_id <- tx_anno$feature[match(names(final_anno), tx_anno$peak)]
mcols(final_anno)$exon_overlap <- !is.na(exon_anno$feature[match(names(final_anno), exon_anno$peak)])
```

---

## Summary: Key Best Practices Checklist

This checklist incorporates both ChIPpeakAnno-specific guidance and community-wide best practices (see Section 1 for community standards).

### Before Annotation
- [ ] Verify genome version compatibility (use latest assemblies: GRCh38/hg38, GRCm39/mm39)
- [ ] Use comprehensive annotation sources (GENCODE, Ensembl, RefSeq)
- [ ] Standardize chromosome naming
- [ ] Filter low-quality peaks (check FRiP, NSC, RSC scores - see Section 1.3)
- [ ] Remove blacklisted regions (ENCODE blacklist)
- [ ] Convert to GRanges format
- [ ] Document annotation source and version

### During Annotation
- [ ] **Consult the [Parameter Guide](Peak_Annotation_Parameter_Guide.md)** to select optimal parameters based on data characteristics
- [ ] Choose appropriate method (peak-centric vs feature-centric) - see [Section 3](#3-annotation-strategy-selection)
- [ ] Use community-standard distance calculation (`middle` + `TSS` - see [Section 4.4](#44-best-practice-consistent-distance-calculation))
- [ ] Use standard distance cutoffs (promoter: -2000 to +500 bp - see [Section 1.2](#12-standard-distance-cutoffs))
- [ ] Set appropriate `maxgap` or `bindingRegion` - see [Parameter Guide](Peak_Annotation_Parameter_Guide.md)
- [ ] Handle overlapping features appropriately (`select` parameter) - see [Section 5](#5-handling-overlapping-features)
- [ ] Use precedence rules for distribution analysis ([Section 5.2](#52-feature-precedence-rules))
- [ ] Consider multi-omics integration if data available ([Section 1.6](#16-multi-omics-integration))

### After Annotation
- [ ] Check annotation coverage rate (expected: 40-60% promoters for TFs) - see [Section 7.2](#72-annotation-coverage)
- [ ] Validate distance distributions - see [Section 7.3](#73-distance-distribution)
- [ ] Examine genomic element distribution - see [Section 7.4](#74-genomic-element-distribution)
- [ ] Create metagene plot to validate patterns - see [Section 7.4.1](#741-peak-distribution-across-metagene)
- [ ] Perform statistical testing with multiple testing correction - see [Section 8.1](#81-multiple-testing-correction)
- [ ] Document all parameters used - see [Section 9.2](#92-documentation-requirements)
- [ ] Capture session information - see [Section 9.3](#93-session-information)

### Quality Assurance
- [ ] Verify expected annotation patterns - see [Section 7.4](#74-genomic-element-distribution)
- [ ] Check for genome version mismatches - see [Section 10.1](#101-genome-version-mismatch)
- [ ] Validate against known biological expectations
- [ ] Compare with published similar analyses
- [ ] Cross-validate with multiple tools if possible - see [Section 1.7](#17-community-recommended-annotation-tools)

### Reproducibility
- [ ] Document all software versions - see [Section 9.3](#93-session-information)
- [ ] Save annotation metadata - see [Section 9.2](#92-documentation-requirements)
- [ ] Provide analysis scripts - see [Section 9.4](#94-workflow-documentation)
- [ ] Use version control - see [Section 9.5](#95-version-control)
- [ ] Follow FAIR data principles - see [Section 9.1](#91-fair-data-principles)

---

## References and Further Reading

### Community Standards and Guidelines

1. **ENCODE Consortium**:
   - ENCODE Project Consortium. (2012) An integrated encyclopedia of DNA elements in the human genome. *Nature* 489:57-74
   - ENCODE Data Standards: https://www.encodeproject.org/data-standards/
   - ENCODE ChIP-seq Guidelines: https://www.encodeproject.org/chip-seq/

2. **ChIP-seq Best Practices Papers**:
   - Landt, S.G. et al. (2012) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. *Genome Research* 22:1813-1831
   - Meyer, C.A. & Liu, X.S. (2014) Identifying and mitigating bias in next-generation sequencing methods for chromatin biology. *Nature Reviews Genetics* 15:709-721

3. **Genomic Annotation Standards**:
   - Harrow, J. et al. (2012) GENCODE: The reference human genome annotation for The ENCODE Project. *Genome Research* 22:1760-1774
   - Zerbino, D.R. et al. (2018) Ensembl 2018. *Nucleic Acids Research* 46:D754-D761

4. **Reproducibility and FAIR Data**:
   - Wilkinson, M.D. et al. (2016) The FAIR Guiding Principles for scientific data management and stewardship. *Scientific Data* 3:160018
   - Sandve, G.K. et al. (2013) Ten simple rules for reproducible computational research. *PLoS Computational Biology* 9:e1003285

### ChIPpeakAnno Package

1. **Primary Reference**:
   - Zhu, L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. *BMC Bioinformatics* 11:237

2. **Package Documentation**:
   - Bioconductor: https://bioconductor.org/packages/ChIPpeakAnno/
   - Vignette: `browseVignettes("ChIPpeakAnno")`

### Additional Resources

1. **Annotation Tools**:
   - UROPA: Open-source peak annotation tool with flexible rule-based system
   - ChIPseeker: R package for ChIP peak annotation and visualization
   - HOMER: Suite of tools for motif discovery and ChIP-seq analysis
   - GREAT: Web-based tool for functional annotation of genomic regions

2. **Quality Control**:
   - ENCODE Quality Metrics: https://www.encodeproject.org/data-standards/terms/
   - Phanstiel, D.H. et al. (2014) MAnorm: a robust model for quantitative comparison of ChIP-seq data sets. *Genome Biology* 15:472

3. **Multi-Omics Integration**:
   - Integrative analysis combining ChIP-seq with RNA-seq, ATAC-seq, and other data types
   - 3D chromatin structure (Hi-C, ChIA-PET) for enhancer-gene linking

---

## Version Information

- **Document Version**: 1.0
- **Last Updated**: 2024
- **Based on**: ChIPpeakAnno v3.45.2
- **Compatible with**: R >= 3.5, Bioconductor >= 3.8

