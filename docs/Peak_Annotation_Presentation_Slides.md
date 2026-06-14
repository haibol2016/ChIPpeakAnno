---
marp: true
theme: default
paginate: true
size: 16:9
style: |
  section {
    font-family: 'Arial', sans-serif;
    font-size: 28px;  /* Default Marp font size - adjust as needed */
  }
  h1 {
    color: #2c3e50;
  }
  h2 {
    text-align: center;
    color:rgb(12, 2, 49);
  }
  code {
    background-color: #f4f4f4;
    padding: 2px 4px;
    border-radius: 3px;
  }
  table {
    width: 100%;
    font-size: 0.7em;
    table-layout: fixed;
  }
  table td, table th {
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    padding: 8px 12px;
    text-align: left;
  }
  table th:nth-child(1), table td:nth-child(1) {
    width: 15%;
  }
  table th:nth-child(2), table td:nth-child(2) {
    width: 12%;
  }
  table th:nth-child(3), table td:nth-child(3) {
    width: 25%;
  }
  table th:nth-child(4), table td:nth-child(4) {
    width: 48%;
  }
  img {
    display: block;
    margin: 0 auto;
  }
  .center {
    text-align: center;
  }
---

## Peak Annotation with ChIPpeakAnno  
 
<style scoped>
p, strong {
  text-align: center;
  display: block;
}
</style>

<p><strong>Haibo Liu</strong></p>
<p><strong>Jan 13, 2026</strong></p> 

---

## Why Peak Annotation Matters

**Transforms genomic coordinates into biological meaning:**

- **Gene Association**: Identify regulated genes
- **Regulatory Context**: Promoters, enhancers, gene bodies
- **Functional Prediction**: Enable pathway enrichment
- **Data Integration**: Bridge to RNA-seq, epigenomics, 3D structure

**Challenge**: Different factors require different annotation strategies
![height:6px,width:800px](figures/annotation.png)

---

## Workflow Overview - Part 1

### Core Peak Annotation Pipeline

**Per Replicate (Sample-wise):**
```
1. Import Peaks (toGRanges - per replicate)
2. Quality Assessment (assessPeaks - per replicate)
3. Filter Peaks (filterPeaks - per replicate)
```

**Per Condition:**
```
4. Create Consensus Peaks (findOverlapsOfPeaks/IDRfilter)
5. Preliminary Annotation (getGenomicAnnotation)
6. Hierarchical Annotation (annotateHierarchically)
7. Enhancer Identification (findEnhancers)
```

---

## Workflow Overview - Part 2

### Conditional Steps & Downstream Analysis

**For Multi-Condition Experiments:**
```
8. Differential Peak Analysis (diffBind/DESeq2)
9. Annotate Differential Peaks (ChIPpeakAnno)
```

**Downstream Analysis (on annotated peaks):**
```
10. Peak Overlap Enrichment Analysis 
    (Locus Overlap Analysis (LOLA), ChIPpeakAnno Peak Overlap)
11. Gene Set Enrichment (test_enrichment - MSigDB/Enrich)
12. Motif Analysis (MEME/HOMER)
```
---

## Step 1-3 - Per-Replicate Processing (I)

### Step 1: Import Peaks (Per Replicate)
Use `toGRanges()` to import peak files for each replicate

### Step 2: Quality Assessment (Sample-wise)
**assessPeaks()**: Comprehensive Peak Evaluation **per replicate**

**What it assesses:**
- **Peak Width Distribution**: Identify narrow (TF) vs broad (histone) peaks
- **Score Metrics**: p-value, q-value, signalValue, fold enrichment
- **Per-replicate statistics**: Individual assessment before consensus

---

## Step 1-3 - Per-Replicate Processing (II)

**Key Metrics:**
- Width: 200-500 bp (TFs), 1-5 kb (histones) (ChIPpeakAnno)
- FRiP: >0.05 (TFs), >0.10 (histones) - ENCODE standards (ChIPQC)
- NSC: >1.05, RSC: >1.0 (ChIPQC)

**Output**: Statistical summaries + visualization plots per replicate

Use `assessPeaks()` on each replicate individually or on a list of all replicates for comparison

### Step 3: Filter Peaks (Per Replicate)
Use `filterPeaks()` to filter each replicate based on assessment results (width, scores)

---

## Step 4 - Consensus Peak Creation

**After filtering each replicate, create consensus peaks per condition:**

**Strategies:** 
1. **Intersection**: Most conservative (peaks in ALL replicates)
2. **Union**: Most comprehensive (peaks in ANY replicate)
3. **IDR filtering**: Statistical approach (recommended by ENCODE)
4. **Majority voting**: Peaks in N out of M replicates (e.g., 2 of 3)
5. **Merge overlapping**: Combine nearby peaks across replicates

**Tools**: bedtools, MSPC, DiffBind, pycisTopic, Peak Finder MetaServer
**Recommended**: Use IDR for 2 replicates, majority voting for 3+

Use `IDRfilter()` for 2 replicates, or `findOverlapsOfPeaks()` for 3+ replicates to create consensus peaks or other tools list above

---

## Step 5 - Preliminary Annotation (I)

**Purpose**: Understand peak distribution across genomic features

**Applied to**: Consensus peaks (single condition) OR differential peaks (multi-condition)

**Features annotated:**
- Promoters (with distance bins: ≤1kb, 1-2kb, 2-3kb)
- 5' UTR, 3' UTR
- Exons (first, other)
- Introns (first, other)
- Immediate downstream
- Distal intergenic

---

## Step 5 - Preliminary Annotation (II)

**Priority-based assignment**: When peaks overlap multiple features

**Key insight**: Guides strategy selection for hierarchical annotation

Use `getGenomicAnnotation()` on consensus peaks (single condition) or differential peaks (multi-condition)

---

## Peak Distribution Over Genomic Features

<style scoped>
.figures-container {
  display: flex;
  justify-content: center;
  gap: 20px;
  align-items: center;
}
.figures-container img {
  max-height: 400px;
  max-width: 550px;
  display: block;
  object-fit: contain;
}
</style>

<div class="figures-container">
<img src="figures/metagene.png" alt="Metagene" style="height: 300px; width: 350px;">
<img src="figures/genomicfeaturedistribution.png" alt="Genomic Feature Distribution" style="height: 300px; width: 400px;">
</div>

---

## Histone Mark distribution
<style scoped>
.figures-container {
  display: flex;
  justify-content: center;
  gap: 20px;
  align-items: center;
}
.figures-container img {
  height: 600px;
  width: 800px;
  display: block;
}
</style>

<div class="figures-container">
<img src="figures/TF-histone-distribution.jpg" alt="promoter/gene body" style="height: 400px; width: 300px;">
<img src="figures/The-distribution-of-histone-modifications-across-transcriptional-regulatory-elements-of.png" alt="Genomic Feature Distribution" style="height: 400px; width: 200px;">
</div>

---

## Step 6 - Hierarchical Annotation (I)

### annotateHierarchically(): Factor-Specific Annotation

**Applied to**: Consensus peaks (single condition) OR differential peaks (multi-condition)

**Key Differences**: Multiple strategies + factor-specific prioritization

**How it works:**
1. **Parallel Annotation**: Applies multiple strategies simultaneously
   - Promoter-focused (overlapping, maxgap=3000)
   - Upstream-to-downstream (upstream2downstream)
   - Bidirectional promoters

---

## Step 6 - Hierarchical Annotation (II)

2. **Factor-Specific Prioritization**: 
   - TF: Prioritize promoters
   - H3K27me3: Prioritize gene bodies
   - H3K4me3: Prioritize TSS
   - ...
3. **Special Case Handling**:
   - Bidirectional promoters (preserved)
   - Broad peaks covering features (preserved)

**Result**: Best annotation per peak based on factor biology

Use `annotateHierarchically()` with appropriate `factor_type` on consensus peaks or differential peaks

---

## Step 7 - Enhancer Identification (I)

### findEnhancers(): 3D Chromatin Structure Integration

**Applied to**: Annotated consensus peaks OR annotated differential peaks

**Two approaches:**

1. **With Hi-C data** (preferred):
   - Uses chromatin interaction data
   - Identifies long-range enhancer-gene associations
   - Accounts for 3D proximity

---

## Step 7 - Enhancer Identification (II)

2. **With ENCODE CRE data** (human/mouse):
   - Uses curated cis-regulatory elements
   - Validated enhancer annotations
   - No Hi-C required

**Output**: Peaks associated with enhancers and their target genes

Use `findEnhancers()` with Hi-C data, or use ENCODE CRE annotations directly

---

## Identification of Promoter/Enhancer Peaks using ENCODE/RoadMap Epigenetics database

![width:1200px](figures/browser_tags_peaks.png)

---

## Step 8 - Differential Analysis (Multi-Condition Only)(I)

### Differential Peak Analysis

**When**: Only for multi-condition experiments (e.g., treatment vs control)

**Purpose**: Identify condition-specific peaks

**Conventional Methods:**
- **diffBind**: Designed for ChIP-seq differential analysis
- **DESeq2**: General-purpose, requires count matrices
- **edgeR**: Alternative to DESeq2

---

## Step 8 - Differential Analysis (Multi-Condition Only)(II)

**Alternative Method:**
- Generate consensus peaks using tools of your choice
- Convert consensus peaks (BED) into SAF format
- **featureCounts**: Generate count table
- **DESeq2**: Differential peak analysis

**Output**: Statistically significant differential peaks

Use `diffBind` (recommended for ChIP-seq) or `DESeq2`/`edgeR` for differential analysis, then filter for significant peaks

---

## Steps 10-12 - Downstream Analysis (I)

**Applied to**: Annotated consensus peaks OR annotated differential peaks

**10. Peakset Enrichment Analysis:**
- **LOLA (Local Enrichment Analysis)**: Compare against database peak sets (ENCODE, ROADMAP Epigenomics)
- Tissue/cell-type specific enrichment patterns
- Histone mark co-occurrence patterns (chromatin state)

**11. Gene Set Enrichment Analysis:**
- Extract associated genes from annotations
- MSigDB, Enrich (GO, KEGG, Reactome, ...)
- Hypergeometric testing (Choose the right universe)

---

## Steps 10-12 - Downstream Analysis (I)

**12. Motif Analysis:**
- De novo motif discovery (MEME, HOMER)
- Known motif enrichment (JASPAR, TRANSFAC)
- Motif position analysis (FIMO)

![width:600px](figures/motifanalysis.png)

**Integration**: Link motifs → peaks → genes → pathways

Use `runLOLA()` for enrichment analysis, `getEnrichedGO()`/`getEnrichedPATH()`/`test_enrichment()` for gene set enrichment, and `getAllPeakSequence()` + `write2FASTA()` for motif analysis

---

## MEME Suite for Motif Analysis
![width:800px](figures/meme_suite.png)



## Common Pitfalls & Solutions

### Avoiding Annotation Errors

**Pitfall 1**: Using wrong reference point - ❌ Using "start" for broad peaks - ✅ Use "middle" for peaks >200 bp

**Pitfall 2**: Ignoring factor-specific patterns - ❌ Same strategy for TF and histone marks - ✅ Use `annotateHierarchically()` with `factor_type`

**Pitfall 3**: Not assessing quality first - ❌ Annotating without quality check - ✅ Use `assessPeaks()` before annotation

**Pitfall 4**: Over-filtering - ❌ Too strict filtering loses biological signal - ✅ Use assessment results to guide filtering

**Pitfall 5**: Ignoring reproducibility - ❌ Single replicate analysis - ✅ Compare replicates, use IDR

---

## Choosing the Right Strategy

| Factor Type | Peak Width | Annotation Strategy | Key Parameters |
|------------|------------|---------------------|----------------|
| **TF** | 200-500 bp | Promoter-focused | `factor_type="TF"`, `PeakLocForDistance="middle"` |
| **H3K4me3** | 200-500 bp | TSS-focused | `bindingRegion=c(-2000, 500)` |
| **H3K27ac** | 500-2000 bp | Promoter + Enhancer | `bindingRegion=c(-5000, 3000)` |
| **H3K36me3** | 1-5 kb | Gene body | `bindingType="fullRange"` |
| **H3K27me3** | 5-20 kb | Broad domains | `bindingType="fullRange"`, `PeakLocForDistance="middle"` |
| **Pol II** | Variable | TSS + Gene body | Method-specific (ChIP-seq vs GRO-seq) |

**Key**: Match annotation strategy to biological binding pattern

---

## **New Functions/Features Added**

1. **`getGenomicAnnotation()`**: Quick genomic context assignment (Promoter, UTR, exon, intron, intergenic)
2. **`annotateHierarchically()`**: Intelligent multi-strategy annotation with factor-specific prioritization
3. **Annotation Prioritization System**: 15+ factor-specific functions with Jaccard/distance scoring
4. **Enrichment Integration**: LOLA, MSigDB (20+ species), EnrichR for gene set enrichment
5. **Bidirectional Promoter Analysis**: `annotatePeaksNearBDP()`
6. **GTF Parser**: `getGeneIdSymbolBiotypeFromGTF()` for Ensembl GTF files

---

## Package Development Highlights

**Key Improvements:**
- **Robustness**: Fixed critical bugs in core annotation functions
- **Usability**: Comprehensive documentation and visual guides
- **Functionality**: 15+ new prioritization functions for factor-specific analysis
- **Integration**: Native support for major enrichment databases
- **Performance**: Parallel processing support in hierarchical annotation

**Documentation Efforts:**
- Created detailed parameter explanation diagrams
- Developed best practices guides
- Established workflow documentation
- Enhanced function documentation throughout

---

## Complete Workflow Summary (I)

**Per Replicate:**
```
✅ Quality Assessment (assessPeaks - sample-wise)
   → Filter based on metrics (filterPeaks)
```

**Per Condition:**
```
✅ Create Consensus Peaks (findOverlapsOfPeaks/IDRfilter)
   → Preliminary Annotation (getGenomicAnnotation)
   → Hierarchical Annotation (annotateHierarchically)
   → Enhancer Identification (findEnhancers)
```

**Multi-Condition Only:**
```
✅ Differential Analysis (diffBind/DESeq2)
   → Annotate Differential Peaks
```
---
## Complete Workflow Summary (II)
**Downstream (on annotated peaks):**
```
✅ Enrichment Analysis (LOLA + Peak Overlap)
   → Gene Set Enrichment (MSigDB, Enrich)
   → Motif Analysis (MEME/HOMER)
```
---

## Resources & Conclusion

### Resources

- **ChIPpeakAnno Package**: Bioconductor
- **Documentation**: 
  - `Peak_Annotation_Best_Practices.md`
  - `annotatePeakInBatch_parameter_explanation_Diagram.md`
- **ENCODE Standards**: https://www.encodeproject.org/


### Key Takeaways

1. **Quality first**: Always assess before annotating
2. **Factor-specific**: Match strategy to biology
3. **Hierarchical approach**: Multiple perspectives → best annotation
4. **Integration**: Combine annotation with other omics data


<!--->
### Potential Questions

- **Q**: How to choose between different annotation methods?
  - **A**: Use `getGenomicAnnotation()` first to understand distribution, then choose strategy based on factor type and biological question.

- **Q**: What if I don't have replicates?
  - **A**: Still do quality assessment, but be cautious with filtering. Consider using more lenient thresholds.

- **Q**: How to handle very broad peaks (>10 kb)?
  - **A**: Use `PeakLocForDistance="middle"` and `bindingType="fullRange"`. Consider splitting very large peaks.

- **Q**: Can I use this for ATAC-seq?
  - **A**: Yes! Use `factor_type="ATAC"` and adjust quality thresholds (higher FRiP expected).
<!--->