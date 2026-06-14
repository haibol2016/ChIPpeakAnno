# Peak Annotation Quick Reference Handout
## ChIPpeakAnno Best Practices & Parameter Guide

---

## Table 1: Factor-Specific Annotation Strategies

| Factor Type | Peak Width | `factor_type` | `is_promoter_binding` | `PeakLocForDistance` | `bindingRegion` | `bindingType` |
|------------|------------|---------------|----------------------|---------------------|-----------------|---------------|
| **TF** | 200-500 bp | `"TF"` | `TRUE` | `"middle"` | `c(-2000, 500)` | `"startSite"` |
| **H3K4me3** | 200-500 bp | `"H3K4me3"` | `TRUE` | `"middle"` | `c(-2000, 500)` | `"startSite"` |
| **H3K27ac** | 500-2000 bp | `"H3K27ac"` | `TRUE` | `"middle"` | `c(-5000, 3000)` | `"startSite"` |
| **H3K36me3** | 1-5 kb | `"H3K36me3"` | `FALSE` | `"middle"` ⭐ | `c(-5000, 5000)` | `"fullRange"` |
| **H3K27me3** | 5-20 kb | `"H3K27me3"` | `FALSE` | `"middle"` ⭐ | `c(-5000, 5000)` | `"fullRange"` |
| **H3K9me3** | 5-20 kb | `"H3K9me3"` | `FALSE` | `"middle"` ⭐ | `c(-5000, 5000)` | `"fullRange"` |
| **Pol II (ChIP-seq)** | Variable | `"PolII"` | `TRUE` | `"middle"` | `c(-2000, 500)` | `"startSite"` |
| **Pol II (GRO-seq)** | Variable | `"PolII"` | `TRUE` | `"endMinusStart"` | `c(-2000, 500)` | `"startSite"` |
| **RBP** | 50-200 bp | `"RBP"` | `FALSE` | `"middle"` | N/A | Transcript-level |
| **ATAC-seq** | 100-500 bp | `"ATAC"` | `TRUE` | `"middle"` | `c(-5000, 3000)` | `"startSite"` |

⭐ **Required** for broad peaks

---

## Table 2: Quality Control Thresholds (ENCODE Standards)

| Metric | Transcription Factors | Histone Marks | ATAC-seq | Cut&RUN/Cut&Tag |
|--------|----------------------|---------------|----------|-----------------|
| **FRiP** | >0.05 (excellent)<br/>0.01-0.05 (good) | >0.10 (excellent)<br/>0.05-0.10 (good) | >0.30 (excellent)<br/>0.20-0.30 (good) | >0.40 (excellent)<br/>0.30-0.40 (good) |
| **NSC** | >1.05 | >1.05 | >1.1 | >1.1 |
| **RSC** | >1.0 | >1.0 | >1.2 | >1.2 |
| **Peak Width** | 200-500 bp | 1-5 kb (typical)<br/>5-20 kb (broad) | 100-500 bp | 200-1000 bp |
| **Peak Number** | 10K-100K | 10K-100K | 20K-100K | 10K-50K |

---

## Table 3: `annotatePeakInBatch` Output Modes

| Output Mode | Use Case | Uses `maxgap`? | Uses Reference Points? | Key Parameters |
|------------|----------|---------------|------------------------|----------------|
| `"nearestLocation"` | Find closest gene | ❌ No | ✅ Yes (both) | `PeakLocForDistance`, `FeatureLocForDistance` |
| `"overlapping"` | Promoter/gene body | ✅ Yes | ❌ No | `maxgap` |
| `"both"` | Comprehensive | ✅ Yes | ✅ Yes (nearest part) | All above |
| `"shortestDistance"` | Minimum boundary distance | ❌ No | ❌ No | None |
| `"inside"` | Strict containment | ❌ No | ❌ No | None |
| `"upstream&inside"` | Promoter + gene body | ✅ Yes | ✅ Yes (TSS only) | `maxgap`, `FeatureLocForDistance="TSS"` |
| `"inside&downstream"` | Gene body + 3' region | ✅ Yes | ✅ Yes (geneEnd only) | `maxgap`, `FeatureLocForDistance="geneEnd"` |
| `"upstream"` | Upstream regulatory | ✅ Yes | ✅ Yes (TSS only) | `maxgap`, `FeatureLocForDistance="TSS"` |
| `"downstream"` | 3' regulatory | ✅ Yes | ✅ Yes (geneEnd only) | `maxgap`, `FeatureLocForDistance="geneEnd"` |
| `"upstreamORdownstream"` | Both directions | ✅ Yes | ✅ Yes (TSS & geneEnd) | `maxgap` |
| `"upstream2downstream"` | All regions | ✅ Yes | ❌ No | `maxgap` |
| `"nearestBiDirectionalPromoters"` | Bidirectional promoters | ❌ No | ✅ Yes (bindingType) | `bindingRegion` (required) |

---

## Table 4: Reference Point Selection Guide

| Parameter | Option | When to Use | Example |
|-----------|--------|-------------|---------|
| **`PeakLocForDistance`** | `"start"` | Narrow peaks (<200 bp) with 5' bias | Certain TFs |
| | `"middle"` ⭐ | **Recommended** for most cases, especially peaks ≥200 bp | Most ChIP-seq, ATAC-seq |
| | `"end"` | Narrow peaks with 3' bias | 3' end analysis |
| | `"endMinusStart"` | Strand-specific assays | GRO-seq, PRO-seq |
| **`FeatureLocForDistance`** | `"TSS"` ⭐ | **Default, recommended** for most factors | TFs, H3K4me3, H3K27ac |
| | `"geneEnd"` | 3' regulatory elements | PolyA sites, 3' UTR factors |
| | `"middle"` | Gene body analysis | H3K36me3 |
| | `"start"` | Non-stranded features | Custom annotations |
| | `"end"` | Non-stranded features | Custom annotations |

⭐ **Recommended default**

---

## Table 5: `maxgap` Parameter Values

| Value | Meaning | Use Case |
|-------|---------|----------|
| `-1L` (default) | No gap allowed (strict overlap) | Standard overlap detection |
| `0L` | Allow adjacent intervals | Touching regions |
| `1000L` | Allow 1 kb gap | Promoter-proximal analysis |
| `2000L` | Allow 2 kb gap | Standard promoter analysis |
| `5000L` | Allow 5 kb gap | Extended promoter/enhancer |
| `10000L` | Allow 10 kb gap | Distal regulatory elements |

**Note**: `maxgap` is **ignored** when using `annoPeaks()` (uses `bindingRegion` instead)

---

## Table 6: Common Workflow Functions

| Function | Purpose | Input | Output | Time Estimate | When Used |
|----------|---------|-------|--------|---------------|-----------|
| `toGRanges()` | Import peaks | BED/narrowPeak/MACS files | GRanges | <1 min | Per replicate |
| `assessPeaks()` | Quality assessment | GRanges or list | Assessment object + plots | 2-5 min | **Sample-wise** (per replicate) |
| `filterPeaks()` | Filter by width/scores | GRanges, thresholds | Filtered GRanges | <1 min | **Sample-wise** (per replicate) |
| `findOverlapsOfPeaks()` | Create consensus peaks | List of GRanges (replicates) | OverlappingPeaks object | 2-5 min | **Per condition** (after filtering) |
| `IDRfilter()` | IDR-based consensus | 2 GRanges (replicates) | Filtered GRanges | 2-5 min | **Per condition** (2 replicates) |
| `getGenomicAnnotation()` | Preliminary annotation | GRanges, TxDb/EnsDb | GRanges with feature types | 3-5 min | On **consensus** or **differential** peaks |
| `annotateHierarchically()` | Best annotation | GRanges, EnsDb, factor_type | GRanges with best annotation | 5-10 min | On **consensus** or **differential** peaks |
| `annotatePeakInBatch()` | Single-strategy annotation | GRanges, annotation data | GRanges with annotations | 2-5 min | Alternative to hierarchical |
| `annoPeaks()` | Region-based annotation | GRanges, bindingRegion | GRanges with overlaps | 2-5 min | Alternative to hierarchical |
| `findEnhancers()` | Enhancer identification | GRanges, Hi-C data | GRanges with enhancer links | 2-5 min | On annotated peaks |
| `diffBind()` / `DESeq2` | Differential analysis | Count matrices, conditions | Differential peaks | 5-10 min | **Multi-condition only** |

---

## Table 7: Troubleshooting Common Issues

| Issue | Symptom | Solution |
|-------|---------|----------|
| **No annotations returned** | Empty GRanges | Check seqlevels style match between peaks and annotation |
| **Too many annotations** | Multiple per peak | Use `select="first"` or `keep="best"` |
| **Wrong gene associations** | Unexpected genes | Check `PeakLocForDistance` (use `"middle"` for broad peaks) |
| **Missing promoter peaks** | Expected but not found | Increase `maxgap` or use `bindingRegion` |
| **Slow performance** | Long runtime | Filter peaks first, use parallel processing with BiocParallel |
| **Strand issues** | Wrong strand annotations | Check `ignore.strand` parameter (usually `TRUE` for ChIP-seq) |
| **Memory errors** | Out of memory | Process by chromosome, reduce annotation data size |
| **Factor-specific wrong** | Poor prioritization | Verify `factor_type` matches your factor |

---

## Table 8: Annotation Data Sources

| Source | Package | Function | Use Case |
|--------|---------|----------|----------|
| **EnsDb** | `ensembldb` | `transcripts(EnsDb)` | Recommended, version-matched |
| **TxDb** | `GenomicFeatures` | `transcripts(TxDb)` | Alternative to EnsDb |
| **biomaRt** | `biomaRt` | `getAnnotation()` | On-the-fly download |
| **ENCODE** | Direct download | BED files | Reference datasets |
| **ROADMAP** | Direct download | BED files | Epigenomic reference |

**Best Practice**: Use EnsDb matching your genome assembly version

---

## Table 9: Quick Parameter Selection by Question

| Research Question | Function | Key Parameters |
|-------------------|----------|----------------|
| "Which genes are nearest?" | `annotatePeakInBatch` | `output="nearestLocation"`, `PeakLocForDistance="middle"` |
| "What's the genomic distribution?" | `getGenomicAnnotation` | `tssRegion=c(-3000, 3000)` |
| "Best annotation for my factor?" | `annotateHierarchically` | `factor_type`, `keep="best"` |
| "Peaks in promoters?" | `annotatePeakInBatch` | `output="overlapping"`, `maxgap=2000` |
| "Peaks near TSS?" | `annoPeaks` | `bindingRegion=c(-2000, 500)`, `bindingType="startSite"` |
| "Enhancer associations?" | `findEnhancers` | Hi-C data or ENCODE CRE |
| "All annotations?" | `annotateHierarchically` | `keep="all"` |

---

## Table 10: Code Snippets Quick Reference

### Quality Assessment
```r
assessment <- assessPeaks(peaks_list)
assessment$width_plot
assessment$summary
```

### Filtering
```r
# Using filterPeaks() function
peaks_filtered <- filterPeaks(
  peaks,
  min_width = 100,
  max_width = 5000,
  score_thresholds = list(pValue = c(0.01, NULL))
)

# Alternative: Manual filtering
peaks_filtered <- peaks[
  width(peaks) >= 100 & width(peaks) <= 5000 &
  peaks$pValue < 0.01
]
```

### Preliminary Annotation
```r
genomic_anno <- getGenomicAnnotation(
  peaks, TxDb = EnsDb.Hsapiens.v86,
  tssRegion = c(-3000, 3000)
)
```

### Hierarchical Annotation (TF)
```r
anno <- annotateHierarchically(
  peaks, annoData = transcripts,
  factor_type = "TF",
  is_promoter_binding = TRUE,
  keep = "best"
)
```

### Hierarchical Annotation (Histone)
```r
anno <- annotateHierarchically(
  peaks, annoData = transcripts,
  factor_type = "H3K27me3",
  is_promoter_binding = FALSE,
  keep = "best"
)
```

### Simple Nearest Gene
```r
anno <- annotatePeakInBatch(
  peaks, AnnotationData = transcripts,
  output = "nearestLocation",
  PeakLocForDistance = "middle",
  FeatureLocForDistance = "TSS"
)
```

### Promoter Region Annotation
```r
anno <- annoPeaks(
  peaks, annoData = transcripts,
  bindingType = "startSite",
  bindingRegion = c(-2000, 500)
)
```

---

## Table 11: Workflow Decision Tree

```
Start
  │
  ├─→ Need quick overview?
  │     └─→ getGenomicAnnotation()
  │
  ├─→ Need best annotation?
  │     ├─→ Know factor type?
  │     │     └─→ annotateHierarchically(factor_type=...)
  │     │
  │     └─→ Don't know factor type?
  │           └─→ getGenomicAnnotation() first
  │                 └─→ Then annotateHierarchically()
  │
  ├─→ Need simple nearest gene?
  │     └─→ annotatePeakInBatch(output="nearestLocation")
  │
  └─→ Need region-based (e.g., promoters)?
        └─→ annoPeaks(bindingRegion=...)
```

---

## Table 12: File Format Support

| Format | Extension | Function | Notes |
|--------|-----------|----------|-------|
| BED | `.bed` | `toGRanges(format="BED")` | Standard BED format |
| narrowPeak | `.narrowPeak` | `toGRanges(format="narrowPeak")` | ENCODE standard |
| broadPeak | `.broadPeak` | `toGRanges(format="broadPeak")` | ENCODE standard |
| MACS | `.xls` | `toGRanges(format="MACS")` | MACS output |
| MACS2 | `.xls` | `toGRanges(format="MACS2")` | MACS2 output |
| BEDPE | `.bedpe` | `toGRanges()` | Hi-C interaction format |

---

## Table 13: Common Error Messages & Solutions

| Error Message | Cause | Solution |
|---------------|-------|----------|
| `"seqlevels style not recognized"` | Mismatched chromosome naming | Use `seqlevelsStyle()` to match |
| `"No valid mart object"` | biomaRt connection failed | Check internet, use EnsDb instead |
| `"factor_type must be specified"` | Missing factor type | Provide `factor_type` parameter |
| `"bindingRegion must have length 2"` | Wrong bindingRegion format | Use `c(upstream, downstream)` |
| `"No annotations found"` | No overlaps detected | Check seqlevels, increase maxgap |
| `"PeakLocForDistance must be one of..."` | Typo in parameter | Check spelling: "start", "middle", "end" |

---

## Table 14: Performance Optimization Tips

| Tip | Implementation | Benefit |
|-----|----------------|---------|
| **Filter early** | Filter before annotation | Reduces computation time |
| **Process by chromosome** | Split by `seqnames()` | Reduces memory usage |
| **Use parallel processing** | `BPPARAM` in `annotateHierarchically` | 2-4x speedup |
| **Cache annotation data** | Save EnsDb transcripts | Faster repeated runs |
| **Reduce annotation size** | Use TSS only for initial analysis | Faster for large datasets |
| **Use `select="first"`** | When only need one annotation | Faster than `select="all"` |

---

## Table 15: Integration with Other Tools

| Tool | Purpose | Integration Point |
|------|---------|-------------------|
| **DESeq2/diffBind** | Differential analysis | **Multi-condition only** - before annotation, use peak counts |
| **getEnrichedGO/getEnrichedPATH/test_enrichment** | Gene set enrichment | ChIPpeakAnno native functions for GO, pathway, MSigDB enrichment |
| **MEME/HOMER** | Motif analysis | Extract peak sequences |
| **LOLA** | Enrichment analysis | Use annotated peaks as input |
| **Gviz/ggbio** | Visualization | Use GRanges from annotations |
| **ChIPseeker** | Alternative annotation | Compare results for validation |

---

## Quick Command Reference

### Complete Workflow (Copy-Paste Ready)

#### Single Condition Workflow

```r
library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v86)

# STEP 1-3: Per-Replicate Processing (Sample-wise)
# 1. Import peaks (per replicate)
peaks_rep1 <- toGRanges("rep1.bed")
peaks_rep2 <- toGRanges("rep2.bed")
peaks_rep3 <- toGRanges("rep3.bed")

# 2. Quality assessment (sample-wise)
assessment_rep1 <- assessPeaks(peaks_rep1)
assessment_rep2 <- assessPeaks(peaks_rep2)
assessment_rep3 <- assessPeaks(peaks_rep3)

# 3. Filter peaks (per replicate)
peaks_rep1_filtered <- filterPeaks(
  peaks_rep1,
  min_width = 100,
  max_width = 5000,
  score_thresholds = list(pValue = c(0.01, NULL))
)
peaks_rep2_filtered <- filterPeaks(peaks_rep2, ...)
peaks_rep3_filtered <- filterPeaks(peaks_rep3, ...)

# STEP 4: Create Consensus Peaks (Per Condition)
overlaps <- findOverlapsOfPeaks(
  peaks_rep1_filtered,
  peaks_rep2_filtered,
  peaks_rep3_filtered
)
consensus_peaks <- overlaps$mergedPeaks  # or use IDRfilter for 2 replicates

# STEP 5-7: Annotation (on Consensus Peaks)
# 5. Preliminary annotation
genomic_anno <- getGenomicAnnotation(
  consensus_peaks,
  TxDb = EnsDb.Hsapiens.v86
)

# 6. Hierarchical annotation (TF example)
anno <- annotateHierarchically(
  consensus_peaks,
  annoData = transcripts(EnsDb.Hsapiens.v86),
  factor_type = "TF",
  is_promoter_binding = TRUE,
  keep = "best"
)

# 7. Enhancer identification (optional)
enhancers <- findEnhancers(anno, ...)

# STEP 8-10: Downstream Analysis (on Annotated Peaks)
# 8. Enrichment analysis (LOLA + Peak Overlap)
library(LOLA)
enrichment <- runLOLA(anno, universe, regionDB)

# 9. Gene set and pathway enrichment (ChIPpeakAnno)
# GO terms enrichment
go_enrich <- getEnrichedGO(anno, 
                          orgAnn = "org.Hs.eg.db",
                          feature_id_type = "ensembl_gene_id",
                          maxP = 0.01)

# Pathway enrichment (KEGG/Reactome)
path_enrich <- getEnrichedPATH(anno,
                               orgAnn = "org.Hs.eg.db",
                               pathAnn = "reactome.db",
                               feature_id_type = "ensembl_gene_id")

# MSigDB or Enrich gene set enrichment
library(msigdbr)
collections <- load_terms_from_msigdbr(species = "Homo sapiens",
                                      collection = "H")
universe <- keys(org.Hs.eg.db, keytype = "SYMBOL")
msigdb_enrich <- test_enrichment(anno,
                                 id_column = "gene_symbol",
                                 id_type = "gene_symbol",
                                 collections = collections,
                                 universe = universe)

# 10. Motif analysis (external tools)
# Extract peak sequences for MEME/HOMER
seq <- getAllPeakSequence(consensus_peaks, upstream=20, downstream=20, genome=genome)
write2FASTA(seq, file="peaks.fa")
```

#### Multi-Condition Workflow

```r
# STEP 1-4: Same as single condition (per condition)
# ... create consensus_peaks_condition1 and consensus_peaks_condition2 ...

# STEP 5: Differential Analysis (Multi-Condition Only)
library(diffBind)
dba <- dba(sampleSheet = "samples.csv")
dba <- dba.count(dba)
dba <- dba.analyze(dba)
diff_peaks <- dba.report(dba, th = 1)
diff_peaks_sig <- diff_peaks[abs(diff_peaks$Fold) > 1 & diff_peaks$FDR < 0.05]

# STEP 6-7: Annotation (on Differential Peaks)
# 6. Preliminary annotation
genomic_anno <- getGenomicAnnotation(diff_peaks_sig, TxDb = EnsDb.Hsapiens.v86)

# 7. Hierarchical annotation
anno <- annotateHierarchically(
  diff_peaks_sig,
  annoData = transcripts(EnsDb.Hsapiens.v86),
  factor_type = "TF",
  keep = "best"
)

# STEP 8-10: Downstream Analysis (same as single condition)
```

---

## Contact & Resources

- **ChIPpeakAnno**: Bioconductor package
- **Documentation**: See package vignettes and `Peak_Annotation_Best_Practices.md`
- **ENCODE Standards**: https://www.encodeproject.org/
- **Questions**: Check function help pages with `?functionName`

---

*Last Updated: [Date]*

