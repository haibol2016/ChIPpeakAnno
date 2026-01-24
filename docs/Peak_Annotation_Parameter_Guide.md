# Peak Annotation Parameter Selection Guide

## Interactive Questionnaire for `annotatePeakInBatch()` Parameters

This guide helps you determine the optimal parameters for `annotatePeakInBatch()` based on your data characteristics. Answer the questions below step-by-step, and recommendations will be provided for each answer. At the end, you'll have a complete parameter set ready to use.

---

## How to Use This Guide

1. **Work through Steps 1-8** in order, answering each question
2. **Record your answers** - Each step builds on the previous ones
3. **Review the integrated recommendations** - Parameters are suggested based on your answers
4. **Check the scenario examples** - Find a similar scenario to your data
5. **Use the decision tree** (at the end) - For quick reference after you understand the concepts

---

## Step 1: Identify Your Binding Factor Type

**Question 1.1: What type of factor are you studying?**

- [ ] **A. Transcription Factor (TF)** - Sequence-specific DNA-binding protein
- [ ] **B. Histone Modification** - Histone mark (H3K4me3, H3K27ac, etc.)
- [ ] **C. RNA-Binding Protein (RBP)** - Binds to RNA transcripts
- [ ] **D. General Chromatin Binding Protein** - Chromatin remodelers, architectural proteins
- [ ] **E. Transcription Machinery** - RNA polymerase, elongation factors
- [ ] **F. Other/Unknown**

### Recommendations Based on Factor Type

| Answer | Primary Strategy | Initial Parameters | Next Steps |
|--------|-----------------|-------------------|------------|
| **A. TF** | Promoter-focused or nearest gene | `output = "overlapping"` or `"nearestLocation"`<br>`bindingRegion = c(-2000, 500)` | Check metagene plot (Step 4) |
| **B. Histone** | Full gene body or region-specific | `bindingType = "fullRange"` for broad marks<br>`bindingRegion = c(-2000, 500)` for narrow | See Step 1B for specific marks |
| **C. RBP** | Strand-aware, gene body | `output = "overlapping"`<br>`bindingType = "fullRange"`<br>Use transcript-level | Check sequencing method (Step 2) |
| **D. Chromatin** | Full gene body or intergenic | `output = "overlapping"`<br>`bindingType = "fullRange"` | Check metagene plot (Step 4) |
| **E. Transcription Machinery** | Gene body, strand-aware for RNA-seq | `bindingType = "fullRange"`<br>Check if RNA-seq method | Critical: Check Step 2 for strand setting |
| **F. Other** | Exploratory approach | `output = "nearestLocation"`<br>`PeakLocForDistance = "middle"` | Create metagene plot first |

### Step 1B: Histone Mark Specific Recommendations (If you selected B)

**Question 1.1B: Which specific histone mark?**

- [ ] **B1. H3K4me3** - Active promoter mark
- [ ] **B2. H3K27ac** - Active enhancer/promoter mark
- [ ] **B3. H3K36me3** - Gene body mark (elongation)
- [ ] **B4. H3K27me3** - Repressive mark (broad domains)
- [ ] **B5. H3K9me3** - Heterochromatin mark (broad domains)
- [ ] **B6. H3K4me1** - Poised/active enhancer mark
- [ ] **B7. Other histone mark**

**Histone mark-specific parameter recommendations:**

| Mark | Expected Pattern | Recommended Parameters |
|------|-----------------|----------------------|
| **H3K4me3** | Strong at TSS, narrow peaks | `output = "overlapping"`<br>`bindingRegion = c(-2000, 500)`<br>`PeakLocForDistance = "middle"` |
| **H3K27ac** | TSS and enhancers, medium peaks | `output = "overlapping"`<br>`bindingRegion = c(-5000, 3000)`<br>May need enhancer annotation |
| **H3K36me3** | Gene body, low at TSS, broad peaks | `output = "overlapping"`<br>`bindingType = "fullRange"`<br>`PeakLocForDistance = "middle"` (required) |
| **H3K27me3** | Broad domains, very broad peaks | `output = "overlapping"`<br>`bindingType = "fullRange"`<br>`PeakLocForDistance = "middle"` (required) |
| **H3K9me3** | Broad heterochromatin, very broad peaks | `output = "overlapping"`<br>`bindingType = "fullRange"`<br>May be intergenic |
| **H3K4me1** | Enhancers, gene body, medium peaks | `output = "overlapping"`<br>`bindingRegion = c(-5000, 5000)`<br>Check for enhancer enrichment |

---

## Step 2: Identify Your Sequencing Method

**Question 2.1: What sequencing method was used?**

- [ ] **A. ChIP-seq** (Chromatin Immunoprecipitation)
- [ ] **B. ATAC-seq** (Assay for Transposase-Accessible Chromatin)
- [ ] **C. Cut&RUN** (Cleavage Under Targets and Release Using Nuclease)
- [ ] **D. Cut&Tag** (Cleavage Under Targets and Tagmentation)
- [ ] **E. GRO-seq/PRO-seq/NET-seq** (RNA-based methods)
- [ ] **F. Other**

### Integrated Recommendations

| Answer | `ignore.strand` | Quality Thresholds | Combined with Factor Type |
|--------|----------------|-------------------|--------------------------|
| **A. ChIP-seq** | `TRUE` (default) | FRiP > 0.01 (TFs), > 0.05 (histones) | **All DNA-binding factors**: Use `ignore.strand = TRUE`<br>Works for TFs, histones, chromatin proteins |
| **B. ATAC-seq** | `TRUE` | FRiP > 0.20 | Open chromatin regions<br>Use promoter-focused annotation |
| **C. Cut&RUN** | `TRUE` | FRiP > 0.30 | Very low background<br>Higher quality expected |
| **D. Cut&Tag** | `TRUE` | FRiP > 0.30 | Very low background<br>Higher quality expected |
| **E. GRO-seq/PRO-seq** | `FALSE` ⚠️ | Context-dependent | **RNA-based**: Strand IS meaningful<br>Use transcript-level annotation<br>Use `ignore.strand = FALSE` |
| **F. Other** | `TRUE` (unless stranded RNA-based) | Method-specific | Check if method is RNA-based |

**⚠️ Critical Note**: If you selected **E. GRO-seq/PRO-seq/NET-seq**:
- You MUST use `ignore.strand = FALSE`
- Use transcript-level annotation (see Step 5)
- These are RNA-based methods, not DNA-based

**If you selected C. RBP (Step 1) and A. ChIP-seq (Step 2)**:
- RBP ChIP-seq is still DNA-based → use `ignore.strand = TRUE`
- For RBP CLIP-seq (RNA-based), use `ignore.strand = FALSE`

---

## Step 3: Analyze Peak Width Distribution

**Question 3.1: What is the typical peak width in your data?**

**First, calculate your peak width:**
```r
# Check peak width distribution
summary(width(peaks))
hist(width(peaks), breaks = 50, main = "Peak Width Distribution")

# Calculate median
median_width <- median(width(peaks))
cat("Median peak width:", median_width, "bp\n")

# Categorize automatically
if(median_width < 200) {
  cat("Category: Narrow peaks\n")
  category <- "A"
} else if(median_width <= 500) {
  cat("Category: Medium peaks\n")
  category <- "B"
} else if(median_width <= 2000) {
  cat("Category: Broad peaks\n")
  category <- "C"
} else {
  cat("Category: Very broad peaks\n")
  category <- "D"
}
```

- [ ] **A. Narrow peaks** (< 200 bp) - Sharp, focused binding
- [ ] **B. Medium peaks** (200-500 bp) - Typical TF binding
- [ ] **C. Broad peaks** (500-2000 bp) - Histone marks, broad factors
- [ ] **D. Very broad peaks** (> 2000 bp) - Large domains, heterochromatin

### Integrated Recommendations

| Answer | `PeakLocForDistance` | Rationale | Combined with Factor Type |
|--------|---------------------|-----------|--------------------------|
| **A. Narrow** | `"middle"` or `"start"` | Peak is small, both reasonable | Works for all factor types |
| **B. Medium** | `"middle"` ⭐ (recommended) | Most robust to width variation | **Best for TFs, ATAC-seq** |
| **C. Broad** | `"middle"` ⭐⭐ (strongly recommended) | Center less affected by width | **Required for histone marks** |
| **D. Very broad** | `"middle"` ⭐⭐⭐ (required) | Start/end would be misleading | **Required for H3K27me3, H3K9me3** |

**Key Integration**:
- If you selected **B. Histone** (Step 1) and **C or D. Broad/Very broad** (Step 3):
  - You MUST use `PeakLocForDistance = "middle"`
  - Start/end would give incorrect distances for broad peaks

---

## Step 4: Analyze Peak Distribution Over Metagene

**Question 4.1: Where are your peaks distributed relative to genes?**

**First, create a metagene plot:**
```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # Use your TxDb
metagenePlot(peaks, TxDb, upstream = 10000, downstream = 10000)
```

- [ ] **A. Strong enrichment at TSS** (±2 kb) - Promoter-binding
- [ ] **B. Enrichment across gene body** - Gene body-binding
- [ ] **C. Enrichment upstream of TSS** (>5 kb) - Enhancer-binding
- [ ] **D. Enrichment at gene end** - 3' end-binding
- [ ] **E. Broad enrichment** (TSS to TES) - Histone marks, broad factors
- [ ] **F. Mixed pattern** - Multiple binding modes

### Integrated Recommendations

| Answer | `output` | `bindingRegion` / `bindingType` | `FeatureLocForDistance` | Combined Strategy |
|--------|----------|-------------------------------|------------------------|------------------|
| **A. TSS enrichment** | `"overlapping"` or `"nearestLocation"` | `bindingRegion = c(-2000, 500)` | `"TSS"` | **Promoter-focused**<br>Good for TFs, H3K4me3 |
| **B. Gene body** | `"overlapping"` | `bindingType = "fullRange"` | `"TSS"` or `"middle"` | **Gene body binding**<br>Good for H3K36me3, elongation factors |
| **C. Upstream** | `"overlapping"` or `"upstream"` | `bindingRegion = c(-10000, 0)` or `c(-5000, 0)` | `"TSS"` | **Enhancer-binding**<br>May need intergenic annotation |
| **D. Gene end** | `"overlapping"` | `bindingType = "endSite"`<br>`bindingRegion = c(-5000, 3000)` | `"geneEnd"` | **3' end analysis**<br>Polyadenylation factors |
| **E. Broad** | `"overlapping"` | `bindingType = "fullRange"`<br>`bindingRegion = c(-5000, 5000)` | `"TSS"` | **Broad factors**<br>H3K27me3, H3K9me3, chromatin proteins |
| **F. Mixed** | `"nearestLocation"` | Not applicable | `"TSS"` | **Exploratory**<br>Get closest gene, analyze later |

**Integration with Previous Steps**:
- **If Step 1 = A (TF) and Step 4 = A (TSS enrichment)**: Classic promoter-binding TF
  - Use: `output = "overlapping"`, `bindingRegion = c(-2000, 500)`
- **If Step 1 = B (Histone) and Step 4 = E (Broad)**: Broad histone mark
  - Use: `bindingType = "fullRange"`, `PeakLocForDistance = "middle"` (from Step 3)
- **If Step 1 = C (RBP) and Step 4 = B (Gene body)**: RBP binding in gene body
  - Use: `bindingType = "fullRange"`, transcript-level (see Step 5)

---

## Step 5: Determine Annotation Level

**Question 5.1: Do you need transcript-specific information?**

- [ ] **A. Yes** - Gene has multiple transcripts with different promoters
- [ ] **B. No** - Gene has single major transcript, or gene-level is sufficient
- [ ] **C. Unsure** - Start with transcript-level to be safe

### Integrated Recommendations

| Answer | Annotation Level | Function | When to Use |
|--------|-----------------|----------|------------|
| **A. Yes** | Transcript-level | `transcripts(TxDb)` or `transcripts(EnsDb)` | **Required when**:<br>- Multiple transcripts with different promoters<br>- Analyzing alternative promoter usage<br>- RNA-based methods (GRO-seq, PRO-seq) |
| **B. No** | Gene-level | `genes(TxDb)` or `genes(EnsDb)` | **Use when**:<br>- Single major transcript per gene<br>- Gene-level enrichment analysis<br>- Simpler, faster annotation |
| **C. Unsure** | Transcript-level ⭐ (recommended) | `transcripts(TxDb)` | **Best practice**: Start comprehensive<br>Can aggregate to gene-level later |

**Integration with Previous Steps**:
- **If Step 2 = E (GRO-seq/PRO-seq)**: MUST use transcript-level
- **If Step 1 = C (RBP)**: Recommended transcript-level
- **If Step 1 = A (TF) and Step 4 = A (TSS)**: Gene-level usually sufficient
- **If Step 1 = B (Histone)**: Gene-level usually sufficient

**Code:**
```r
# Gene-level (simpler, faster)
annoData <- genes(TxDb)

# Transcript-level (more comprehensive, recommended for RNA-based methods)
annoData <- transcripts(TxDb)
```

---

## Step 6: Handle Multiple Overlapping Features

**Question 6.1: How should overlapping features be handled?**

- [ ] **A. Return all overlapping features** - Complete information, one peak → multiple annotations
- [ ] **B. Return only the nearest feature** - One annotation per peak
- [ ] **C. Return best match** - Closest + most overlapping (only in `annoPeaks`)
- [ ] **D. Return first overlapping** - Simple, one annotation per peak

### Integrated Recommendations

| Answer | `select` Value | Use Case | Integration Notes |
|--------|---------------|----------|------------------|
| **A. All** | `"all"` | Need complete information<br>Downstream analysis can filter | **Use with**: Broad peaks, multiple genes<br>**Example**: Histone marks spanning multiple genes |
| **B. Nearest** | `"first"` (with `output = "nearestLocation"`) | Standard annotation<br>One-to-one mapping | **Use with**: TFs, promoter-focused<br>**Most common choice** |
| **C. Best match** | `"bestOne"` (in `annoPeaks()`, not `annotatePeakInBatch()`) | Want closest + most overlapping | **Requires**: Using `annoPeaks()` instead<br>See "Function Selection" section |
| **D. First** | `"first"` | Simple, quick annotation | **Use with**: Any scenario needing one annotation |

**Integration with Previous Steps**:
- **If Step 1 = B (Histone) and Step 4 = E (Broad)**: Use `select = "all"` (peaks may span multiple genes)
- **If Step 1 = A (TF) and Step 4 = A (TSS)**: Use `select = "first"` (one promoter per peak)
- **If Step 6 = C (Best match)**: You'll need to use `annoPeaks()` instead (see Step 9)

---

## Step 7: Distance Calculation Preferences

**Question 7.1: What reference point do you want to use for features?**

- [ ] **A. Transcription Start Site (TSS)** - Standard for most analyses
- [ ] **B. Gene End (TES)** - For 3' end analysis
- [ ] **C. Gene Middle** - For gene body analysis
- [ ] **D. Feature Start/End** - Non-stranded features

### Integrated Recommendations

| Answer | `FeatureLocForDistance` | Use Case | Integration |
|--------|------------------------|----------|-------------|
| **A. TSS** | `"TSS"` ⭐ (default, recommended) | Most gene annotation<br>Promoter analysis | **Use with**: Steps 1A (TF), 4A (TSS enrichment)<br>**Standard choice** |
| **B. Gene End** | `"geneEnd"` | 3' UTR analysis<br>Polyadenylation sites | **Use with**: Step 4D (Gene end enrichment) |
| **C. Gene Middle** | `"middle"` | Gene body analysis | **Use with**: Step 4B (Gene body), broad factors |
| **D. Start/End** | `"start"` or `"end"` | Non-stranded features<br>Custom annotations | **Use with**: Custom annotation sources |

**Default Recommendation**: Use `"TSS"` unless you have a specific reason (3' end analysis, etc.)

---

## Step 8: Set `maxgap` Parameter

**Question 8.1: How strict should overlap detection be?**

- [ ] **A. Strict overlap only** - Peaks must directly overlap features
- [ ] **B. Allow small gaps** - Peaks within 100-500 bp of features
- [ ] **C. Allow larger gaps** - Peaks within 1-5 kb of features
- [ ] **D. Use bindingRegion instead** - Define specific regions relative to features

### Integrated Recommendations

| Answer | `maxgap` Value | Use Case | Important Note |
|--------|---------------|----------|---------------|
| **A. Strict** | `0` (default) | Only direct overlaps | **Most common** |
| **B. Small gaps** | `100` to `500` | Allow nearby peaks | For slightly offset peaks |
| **C. Large gaps** | `1000` to `5000` | Enhancer analysis<br>Distal binding | For intergenic/enhancer peaks |
| **D. Use bindingRegion** | Not applicable | Precise region definition | **When `bindingRegion` is set, `maxgap` is ignored** |

**Integration with Previous Steps**:
- **If Step 4 = A (TSS enrichment) and you set `bindingRegion = c(-2000, 500)`**: `maxgap` is ignored
- **If Step 4 = C (Upstream/enhancer)**: Consider `maxgap = 1000-5000` if not using `bindingRegion`
- **If Step 4 = A (TSS)**: Usually `maxgap = 0` is fine with `bindingRegion`

---

## Step 9: Function Selection

**Question 9.1: Which function should you use?**

- [ ] **A. Use `annotatePeakInBatch()`** - Standard, comprehensive annotation
- [ ] **B. Use `annoPeaks()`** - Need `select = "bestOne"` or bidirectional promoters

### Decision Guide

| Need | Function | Reason | Integration |
|------|----------|--------|-------------|
| Standard annotation | `annotatePeakInBatch()` ⭐ | More flexible, comprehensive options | **Use for**: Most scenarios |
| Best match selection | `annoPeaks()` with `select = "bestOne"` | Only `annoPeaks()` supports this | **If Step 6 = C (Best match)** |
| Bidirectional promoters | `annoPeaks()` with `bindingType = "nearestBiDirectionalPromoters"` | Specialized function | **If Step 1 = A (TF) and analyzing promoters** |
| Multiple output modes | `annotatePeakInBatch()` | More output options available | **Use for**: Most scenarios |
| Simple nearest gene | Either works | `annotatePeakInBatch()` is more commonly used | **Default choice** |

**Note**: `annotatePeakInBatch()` internally calls `annoPeaks()` for some operations, but provides a more user-friendly interface.

---

## Complete Parameter Sets by Scenario

Based on common combinations of answers, here are ready-to-use parameter sets:

### Scenario 1: Promoter-Binding Transcription Factor

**Answers**: 1.1=A, 2.1=A, 3.1=B, 4.1=A, 5.1=B, 6.1=B, 7.1=A, 8.1=D

```r
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = genes(TxDb),
    output = "overlapping",              # Promoter-focused
    PeakLocForDistance = "middle",       # Medium peaks
    FeatureLocForDistance = "TSS",       # Standard TSS
    bindingRegion = c(-2000, 500),       # Promoter region
    select = "first",                    # One annotation per peak
    ignore.strand = TRUE                 # ChIP-seq is DNA-based
)
```

### Scenario 2: Gene Body Histone Mark (e.g., H3K36me3)

**Answers**: 1.1=B (B3), 2.1=A, 3.1=C, 4.1=B, 5.1=B, 6.1=A, 7.1=A, 8.1=A

```r
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = genes(TxDb),
    output = "overlapping",              # Gene body binding
    PeakLocForDistance = "middle",       # Broad peaks require center
    FeatureLocForDistance = "TSS",       # TSS as reference
    bindingType = "fullRange",           # Full gene range
    select = "all",                      # May have multiple genes
    ignore.strand = TRUE                 # ChIP-seq
)
```

### Scenario 3: Broad Repressive Histone Mark (e.g., H3K27me3)

**Answers**: 1.1=B (B4), 2.1=A, 3.1=D, 4.1=E, 5.1=B, 6.1=A, 7.1=A, 8.1=A

```r
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = genes(TxDb),
    output = "overlapping",              # Broad binding
    PeakLocForDistance = "middle",       # REQUIRED for very broad peaks
    FeatureLocForDistance = "TSS",       # TSS reference
    bindingType = "fullRange",           # Full gene range
    bindingRegion = c(-5000, 5000),      # Extended region
    select = "all",                      # May span multiple genes
    ignore.strand = TRUE                 # ChIP-seq
)
```

### Scenario 4: Enhancer-Binding Factor

**Answers**: 1.1=A, 2.1=A, 3.1=B, 4.1=C, 5.1=B, 6.1=B, 7.1=A, 8.1=C

```r
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = genes(TxDb),
    output = "nearestLocation",          # May be far from genes
    PeakLocForDistance = "middle",       # Medium peaks
    FeatureLocForDistance = "TSS",       # TSS reference
    maxgap = 5000,                       # Allow larger gaps for distal peaks
    select = "first",                    # One annotation
    ignore.strand = TRUE                 # ChIP-seq
)
# Then filter for intergenic or use enhancer databases
intergenic <- annotated[annotated$insideFeature %in% c("upstream", "downstream")]
```

### Scenario 5: RNA-Binding Protein, ChIP-seq

**Answers**: 1.1=C, 2.1=A, 3.1=A, 4.1=B, 5.1=A, 6.1=B, 7.1=A, 8.1=A

```r
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = transcripts(TxDb),  # Transcript-level for RBPs
    output = "overlapping",              # Gene body binding
    PeakLocForDistance = "middle",       # Narrow peaks
    FeatureLocForDistance = "TSS",       # TSS reference
    bindingType = "fullRange",           # Full transcript
    select = "first",                    # One annotation
    ignore.strand = TRUE                 # ChIP-seq (DNA-based, not RNA)
    # Note: For RBP CLIP-seq (RNA-based), use ignore.strand = FALSE
)
```

### Scenario 6: GRO-seq/PRO-seq (RNA-based)

**Answers**: 1.1=E, 2.1=E, 3.1=Any, 4.1=A, 5.1=A, 6.1=B, 7.1=A, 8.1=D

```r
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = transcripts(TxDb),  # Transcript-level
    output = "overlapping",              # TSS regions
    PeakLocForDistance = "middle",       # Standard
    FeatureLocForDistance = "TSS",       # TSS reference
    bindingRegion = c(-2000, 500),       # Promoter
    select = "first",                    # One annotation
    ignore.strand = FALSE                # ⚠️ RNA-based, strand matters!
)
```

### Scenario 7: Bidirectional Promoter Analysis

**Answers**: 1.1=A, 2.1=A, 3.1=B, 4.1=A, 5.1=A, 6.1=A, 7.1=A, 8.1=D

```r
# Must use annoPeaks() for bidirectional promoters
annotated <- annoPeaks(
    peaks,
    annoData = transcripts(TxDb),  # Transcript-level recommended
    bindingType = "nearestBiDirectionalPromoters",
    bindingRegion = c(-5000, 3000),
    select = "all",
    ignore.peak.strand = TRUE
)
```

### Scenario 8: ATAC-seq, Promoter Enrichment

**Answers**: 1.1=F, 2.1=B, 3.1=B, 4.1=A, 5.1=B, 6.1=B, 7.1=A, 8.1=D

```r
annotated <- annotatePeakInBatch(
    peaks,
    AnnotationData = genes(TxDb),
    output = "overlapping",              # Promoter regions
    PeakLocForDistance = "middle",       # Medium peaks
    FeatureLocForDistance = "TSS",       # TSS reference
    bindingRegion = c(-2000, 500),       # Promoter
    select = "first",                    # One annotation
    ignore.strand = TRUE                 # ATAC-seq is DNA-based
)
```

---

## Complete Workflow

**Step-by-step workflow using this guide:**

1. **Prepare your data**:
   ```r
   # Load peaks
   library(ChIPpeakAnno)
   peaks <- toGRanges("peaks.bed", format = "BED")
   
   # Check genome version
   genome(peaks)
   ```

2. **Answer Steps 1-9** of the questionnaire above

3. **Create metagene plot** (for Step 4):
   ```r
   library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # Use your TxDb
   TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
   metagenePlot(peaks, TxDb, upstream = 10000, downstream = 10000)
   ```

4. **Check peak width** (for Step 3):
   ```r
   median_width <- median(width(peaks))
   hist(width(peaks), breaks = 50)
   cat("Median peak width:", median_width, "bp\n")
   ```

5. **Select parameters** based on your answers from Steps 1-9

6. **Prepare annotation data**:
   ```r
   # Based on Step 5 answer
   if(step5_answer == "A" || step5_answer == "C") {
       annoData <- transcripts(TxDb)  # Transcript-level
   } else {
       annoData <- genes(TxDb)  # Gene-level
   }
   ```

7. **Run annotation** with your selected parameters:
   ```r
   annotated <- annotatePeakInBatch(
       peaks,
       AnnotationData = annoData,
       # Add parameters from Steps 1-9
       output = "...",  # From Step 4
       PeakLocForDistance = "...",  # From Step 3
       FeatureLocForDistance = "...",  # From Step 7
       bindingRegion = c(...),  # From Step 4
       select = "...",  # From Step 6
       maxgap = ...,  # From Step 8
       ignore.strand = ...  # From Step 2
   )
   ```

8. **Validate results** (see Validation Checklist below)

9. **Adjust if needed** based on validation results

---

## Validation Checklist

After annotation, validate your parameters:

- [ ] **Annotation rate**: 40-80% of peaks annotated (varies by factor type)
- [ ] **Distance distribution**: Matches expected pattern from metagene plot
- [ ] **Genomic distribution**: Matches biological expectations
- [ ] **Strand handling**: Correct for your method (DNA-seq vs RNA-seq)
- [ ] **Peak width**: Appropriate `PeakLocForDistance` selected

**Code for validation:**
```r
# 1. Check annotation rate
annotation_rate <- sum(!is.na(annotated$feature)) / length(peaks)
cat("Annotation rate:", annotation_rate, "\n")
# Expected: 40-80% for most factors

# 2. Check distance distribution
hist(annotated$distancetoFeature, breaks = 100, 
     main = "Distance to TSS Distribution")
# Should match your metagene plot pattern

# 3. Check genomic distribution
distribution <- assignChromosomeRegion(peaks, TxDb = TxDb)
pie1(distribution$percentage)
# Should match biological expectations for your factor type

# 4. Verify strand handling
table(strand(annotated))  # Should match your ignore.strand setting
# If ignore.strand = TRUE, most should be "*"
# If ignore.strand = FALSE, should have "+" and "-" strands
```

---

## Troubleshooting

**Problem**: Low annotation rate (<20%)
- **Solution**: Check genome version match, verify annotation source, try larger `maxgap` or `bindingRegion`

**Problem**: Unexpected distance distribution
- **Solution**: Verify `PeakLocForDistance` and `FeatureLocForDistance` match your data, check metagene plot

**Problem**: Too many peaks with multiple annotations
- **Solution**: Use `select = "first"` or use `annoPeaks()` with `select = "bestOne"`

**Problem**: Peaks not matching metagene plot pattern
- **Solution**: Adjust `output` and `bindingRegion` based on metagene plot, verify `PeakLocForDistance`

**Problem**: Wrong strand information
- **Solution**: Verify sequencing method - DNA-seq should use `ignore.strand = TRUE`, RNA-seq should use `FALSE`

---

## Parameter Summary Table

| Parameter | Common Values | When to Use | Default |
|-----------|--------------|-------------|---------|
| `output` | `"nearestLocation"` | General annotation, find closest gene | ✅ Default |
| | `"overlapping"` | Promoter/gene body analysis | |
| | `"both"` | Need both nearest and overlapping | |
| `PeakLocForDistance` | `"middle"` ⭐ | Most cases, robust to peak width | `"start"` |
| | `"start"` | Narrow peaks, directional binding | |
| | `"end"` | Narrow peaks, 3' end analysis | |
| `FeatureLocForDistance` | `"TSS"` ⭐ | Standard gene annotation | ✅ Default |
| | `"geneEnd"` | 3' end analysis | |
| | `"middle"` | Gene body analysis | |
| `bindingRegion` | `c(-2000, 500)` | Standard promoter | `c(-5000, 5000)` |
| | `c(-5000, 3000)` | Extended promoter | |
| | `c(-5000, 5000)` | Very broad region | |
| `bindingType` | `"fullRange"` | Gene body binding | `"startSite"` |
| | `"startSite"` | TSS-focused | |
| | `"endSite"` | Gene end-focused | |
| `select` | `"all"` | Need all overlapping features | `"all"` |
| | `"first"` | One annotation per peak | |
| `maxgap` | `0` | Only direct overlaps | ✅ Default |
| | `100-500` | Allow small gaps | |
| | `1000-5000` | Enhancer/distal binding | |
| `ignore.strand` | `TRUE` ⭐ | All DNA-seq methods | ✅ Default |
| | `FALSE` | RNA-seq methods (GRO-seq, PRO-seq) | |

⭐ = Recommended for most cases

---

## Parameter Interaction Guide

**Important parameter interactions:**

1. **`bindingRegion` overrides `maxgap`**:
   - When `bindingRegion` is set, `maxgap` is ignored
   - Use `bindingRegion` for precise region definition

2. **`output` affects `select` behavior**:
   - `output = "nearestLocation"`: `select` determines which nearest feature
   - `output = "overlapping"`: `select` determines which overlapping features

3. **`bindingType` in `annoPeaks()` vs `bindingRegion` in `annotatePeakInBatch()`**:
   - `bindingType` in `annoPeaks()`: "startSite", "endSite", "fullRange", "nearestBiDirectionalPromoters"
   - `bindingRegion` in `annotatePeakInBatch()`: Vector defining region relative to feature
   - Both can be used together in `annotatePeakInBatch()`

4. **Transcript-level vs Gene-level**:
   - Transcript-level: More annotations, may have duplicates, slower
   - Gene-level: Simpler, one annotation per gene, faster

---

## Edge Cases and Special Situations

### Very Large Datasets (>100K peaks)

**Recommendations**:
- ChIPpeakAnno automatically splits large datasets internally (>10K peaks)
- Consider processing in batches if memory is limited
- Use gene-level annotation (faster than transcript-level)

### Mixed Peak Types

**If you have multiple peak types in one dataset**:
- Annotate separately for each type
- Or use `output = "nearestLocation"` to get closest annotation for all

### Non-Standard Genomes

**For non-model organisms**:
- Use custom annotation files (GRanges objects)
- Ensure genome version matches
- May need to create custom TxDb from GFF/GTF

### Low Annotation Rate

**If <20% of peaks are annotated**:
1. Check genome version compatibility
2. Verify chromosome naming (chr1 vs 1)
3. Try larger `maxgap` or `bindingRegion`
4. Check if peaks are in non-annotated regions (enhancers, etc.)

---

## Quick Decision Tree

Use this decision tree for quick parameter selection after understanding the concepts above:

```
START: What factor type? (Step 1)
│
├─ A. Transcription Factor
│   ├─ Sequencing method? (Step 2)
│   │   ├─ ChIP-seq/ATAC-seq/Cut&RUN/Cut&Tag → ignore.strand = TRUE
│   │   └─ GRO-seq/PRO-seq → ignore.strand = FALSE ⚠️
│   ├─ Peak width? (Step 3)
│   │   ├─ < 200 bp → PeakLocForDistance = "middle" or "start"
│   │   ├─ 200-500 bp → PeakLocForDistance = "middle" ⭐
│   │   └─ > 500 bp → PeakLocForDistance = "middle" (required)
│   ├─ Metagene pattern? (Step 4)
│   │   ├─ TSS enrichment → output = "overlapping", bindingRegion = c(-2000, 500)
│   │   ├─ Upstream/enhancer → output = "nearestLocation" or "overlapping", maxgap = 5000
│   │   └─ Mixed → output = "nearestLocation"
│   └─ Annotation level? (Step 5)
│       ├─ Multiple transcripts → transcripts(TxDb)
│       └─ Single transcript → genes(TxDb)
│
├─ B. Histone Modification
│   ├─ Which mark? (Step 1B)
│   │   ├─ H3K4me3 → bindingRegion = c(-2000, 500), narrow peaks
│   │   ├─ H3K27ac → bindingRegion = c(-5000, 3000), medium peaks
│   │   ├─ H3K36me3 → bindingType = "fullRange", broad peaks
│   │   ├─ H3K27me3/H3K9me3 → bindingType = "fullRange", very broad, PeakLocForDistance = "middle" (required)
│   │   └─ H3K4me1 → bindingRegion = c(-5000, 5000), enhancer-focused
│   ├─ Peak width? (Step 3)
│   │   └─ > 500 bp → PeakLocForDistance = "middle" (required) ⭐⭐⭐
│   ├─ Metagene pattern? (Step 4)
│   │   ├─ TSS → bindingRegion = c(-2000, 500)
│   │   ├─ Gene body → bindingType = "fullRange"
│   │   └─ Broad → bindingType = "fullRange", bindingRegion = c(-5000, 5000)
│   └─ ignore.strand = TRUE (all DNA-seq)
│
├─ C. RNA-Binding Protein
│   ├─ Sequencing method? (Step 2)
│   │   ├─ ChIP-seq → ignore.strand = TRUE (DNA-based)
│   │   └─ CLIP-seq → ignore.strand = FALSE (RNA-based) ⚠️
│   ├─ Annotation level → transcripts(TxDb) (recommended)
│   └─ bindingType = "fullRange" (gene body)
│
├─ D. Chromatin Binding Protein
│   ├─ ignore.strand = TRUE (DNA-based)
│   ├─ Peak width > 500 bp? → PeakLocForDistance = "middle"
│   └─ bindingType = "fullRange" (may be intergenic)
│
├─ E. Transcription Machinery
│   ├─ Sequencing method? (Step 2)
│   │   ├─ ChIP-seq → ignore.strand = TRUE
│   │   └─ GRO-seq/PRO-seq → ignore.strand = FALSE ⚠️, transcripts(TxDb)
│   └─ bindingType = "fullRange" (gene body)
│
└─ F. Other/Unknown
    └─ Start with: output = "nearestLocation", PeakLocForDistance = "middle"
       Then check metagene plot and adjust

KEY DECISIONS:
⭐ = Recommended
⚠️ = Critical - don't miss this!
```

---

## References

- See main document: `Peak_Annotation_Best_Practices.md`
- ChIPpeakAnno documentation: `?annotatePeakInBatch`
- Function summary: `ChIPpeakAnno_Function_Summary.md`
