# Peak Annotation Workflow - Visual Diagram

## Complete Workflow Overview

```mermaid
graph TD
    A[Import Peaks<br/>toGRanges - Per Replicate] --> B[Quality Assessment<br/>assessPeaks - Sample-wise]
    B --> C{Assessment<br/>Results}
    C -->|Width OK<br/>Scores OK| D[Filter Peaks<br/>filterPeaks - Per Replicate]
    C -->|Issues Found| E[Review & Adjust<br/>Peak Calling]
    E --> A
    D --> F[Create Consensus Peaks<br/>findOverlapsOfPeaks/IDRfilter]
    F --> G{Multi-Condition<br/>Experiment?}
    G -->|No| H[Preliminary Annotation<br/>getGenomicAnnotation]
    G -->|Yes| I[Differential Analysis<br/>diffBind/DESeq2]
    I --> J[Annotate Differential Peaks<br/>getGenomicAnnotation]
    J --> K[Hierarchical Annotation<br/>annotateHierarchically]
    H --> K
    K --> L{Hi-C Data<br/>Available?}
    L -->|Yes| M[Enhancer Identification<br/>findEnhancers with Hi-C]
    L -->|No| N[Enhancer Identification<br/>ENCODE CRE Data]
    M --> O[Enrichment Analysis<br/>LOLA + Peak Overlap]
    N --> O
    O --> P[Gene Set Enrichment<br/>getEnrichedGO/PATH<br/>test_enrichment - MSigDB/Enrich]
    P --> Q[Motif Analysis<br/>MEME/HOMER/JASPAR]
    Q --> R[Integration &<br/>Interpretation]
    
    style A fill:#e1f5ff
    style B fill:#fff4e1
    style D fill:#ffe1e1
    style F fill:#e1ffe1
    style I fill:#f0e1ff
    style K fill:#f0e1ff
    style O fill:#ffe1f0
    style R fill:#e1f5ff
```

## Detailed Step-by-Step Workflow

```mermaid
graph LR
    subgraph "Phase 1: Per-Replicate Processing"
        A1[Import Peaks<br/>toGRanges - Per Replicate] --> A2[assessPeaks<br/>Sample-wise]
        A2 --> A3[Width Distribution]
        A2 --> A4[Score Metrics]
        A3 --> A5[filterPeaks<br/>Per Replicate]
        A4 --> A5
        A5 --> A6[Filtered Peaks<br/>Per Replicate]
    end
    
    subgraph "Phase 2: Consensus Peak Creation"
        A6 --> B1[findOverlapsOfPeaks<br/>or IDRfilter]
        B1 --> B2[Consensus Peaks<br/>Per Condition]
    end
    
    subgraph "Phase 3: Annotation"
        B2 --> C1{Multi-Condition?}
        C1 -->|No| C2[getGenomicAnnotation<br/>on Consensus Peaks]
        C1 -->|Yes| C3[Differential Analysis<br/>diffBind/DESeq2]
        C3 --> C4[Annotate Differential Peaks<br/>getGenomicAnnotation]
        C2 --> C5[annotateHierarchically]
        C4 --> C5
        C5 --> C6[Best Annotation<br/>per Peak]
    end
    
    subgraph "Phase 4: Enhancer Analysis"
        C6 --> D1{Hi-C Available?}
        D1 -->|Yes| D2[findEnhancers<br/>with Hi-C]
        D1 -->|No| D3[ENCODE CRE<br/>Annotations]
        D2 --> D4[Enhancer-Gene<br/>Associations]
        D3 --> D4
    end
    
    subgraph "Phase 5: Downstream Analysis"
        D4 --> E1[Enrichment Analysis<br/>LOLA + Peak Overlap]
        E1 --> E2[Gene Set Enrichment<br/>getEnrichedGO/PATH<br/>test_enrichment]
        E2 --> E3[Motif Analysis]
        E3 --> E4[Final Results]
    end
    
    style A2 fill:#fff4e1
    style B1 fill:#e1ffe1
    style C1 fill:#f0e1ff
    style D2 fill:#ffe1f0
    style E1 fill:#e1f5ff
```

## Decision Tree for Annotation Strategy Selection

```mermaid
graph TD
    Start[Start Annotation] --> Q1{What is the<br/>factor type?}
    
    Q1 -->|Transcription Factor| TF[TF Strategy]
    Q1 -->|Histone Mark| HM{Which histone?}
    Q1 -->|Pol II| POL{Sequencing<br/>method?}
    Q1 -->|RBP| RBP[RNA-binding<br/>Protein]
    Q1 -->|ATAC-seq| ATAC[Open Chromatin]
    
    TF --> TF1[factor_type = 'TF'<br/>is_promoter_binding = TRUE<br/>PeakLocForDistance = 'middle']
    
    HM --> HM1{H3K4me3?}
    HM --> HM2{H3K27ac?}
    HM --> HM3{H3K36me3?}
    HM --> HM4{H3K27me3?}
    
    HM1 --> HM1A[factor_type = 'H3K4me3'<br/>bindingRegion = c-2000, 500<br/>PeakLocForDistance = 'middle']
    HM2 --> HM2A[factor_type = 'H3K27ac'<br/>bindingRegion = c-5000, 3000<br/>PeakLocForDistance = 'middle']
    HM3 --> HM3A[factor_type = 'H3K36me3'<br/>bindingType = 'fullRange'<br/>PeakLocForDistance = 'middle']
    HM4 --> HM4A[factor_type = 'H3K27me3'<br/>bindingType = 'fullRange'<br/>PeakLocForDistance = 'middle']
    
    POL --> POL1{ChIP-seq?}
    POL --> POL2{GRO-seq/PRO-seq?}
    POL1 --> POL1A[Standard TF strategy]
    POL2 --> POL2A[Strand-aware<br/>PeakLocForDistance = 'endMinusStart']
    
    RBP --> RBP1[factor_type = 'RBP'<br/>Transcript-level annotation<br/>ignore.strand = FALSE]
    
    ATAC --> ATAC1[factor_type = 'ATAC'<br/>Promoter + Enhancer focus]
    
    TF1 --> Apply[Apply annotateHierarchically]
    HM1A --> Apply
    HM2A --> Apply
    HM3A --> Apply
    HM4A --> Apply
    POL1A --> Apply
    POL2A --> Apply
    RBP1 --> Apply
    ATAC1 --> Apply
    
    Apply --> Result[Best Annotation<br/>per Peak]
    
    style Start fill:#e1f5ff
    style Apply fill:#f0e1ff
    style Result fill:#e1ffe1
```

## Quality Control Checkpoints

```mermaid
graph TD
    QC1[Peak Import] --> QC2{Width Check}
    QC2 -->|Median < 200bp| QC2A[✓ TF-like peaks]
    QC2 -->|Median 200-2000bp| QC2B[✓ Histone marks]
    QC2 -->|Median > 2000bp| QC2C[⚠ Very broad<br/>Check peak calling]
    
    QC2A --> QC3{Score Check}
    QC2B --> QC3
    QC2C --> QC3
    
    QC3 -->|p-value OK| QC4{FRiP Check}
    QC3 -->|p-value Low| QC3A[⚠ Review filtering]
    
    QC4 -->|FRiP > 0.05 TF<br/>FRiP > 0.10 Histone| QC5{Reproducibility}
    QC4 -->|FRiP Low| QC4A[⚠ Quality issue<br/>Review experiment]
    
    QC5 -->|IDR OK| QC6[✓ Ready for Annotation]
    QC5 -->|IDR Low| QC5A[⚠ Low reproducibility<br/>Review replicates]
    
    QC3A --> QC2
    QC4A --> QC2
    QC5A --> QC2
    
    style QC6 fill:#e1ffe1
    style QC2C fill:#fff4e1
    style QC3A fill:#fff4e1
    style QC4A fill:#ffe1e1
    style QC5A fill:#ffe1e1
```

## Annotation Method Selection Flow

```mermaid
graph LR
    Input[Filtered Peaks] --> Method{Annotation<br/>Method?}
    
    Method -->|Quick Overview| Genomic[getGenomicAnnotation<br/>Feature-centric<br/>Priority-based]
    Method -->|Best Annotation| Hierarchical[annotateHierarchically<br/>Factor-specific<br/>Multiple strategies]
    Method -->|Simple Nearest| Simple[annotatePeakInBatch<br/>output='nearestLocation'<br/>Single strategy]
    Method -->|Region-based| Region[annoPeaks<br/>bindingRegion defined<br/>Region expansion]
    
    Genomic --> Output1[Distribution Summary]
    Hierarchical --> Output2[Best Annotation<br/>+ Priority Scores]
    Simple --> Output3[Nearest Gene]
    Region --> Output4[Region Overlaps]
    
    Output1 --> Next[Downstream Analysis]
    Output2 --> Next
    Output3 --> Next
    Output4 --> Next
    
    style Hierarchical fill:#f0e1ff
    style Next fill:#e1f5ff
```

## ASCII Art Workflow (Alternative Format)

```
┌─────────────────────────────────────────────────────────────────┐
│                    PEAK ANNOTATION WORKFLOW                      │
└─────────────────────────────────────────────────────────────────┘

PHASE 1: PER-REPLICATE PROCESSING (Sample-wise)
═══════════════════════════════════════════════════════════════════
    Import Peaks (toGRanges - Per Replicate)
         │
         ▼
    ┌─────────────────┐
    │  assessPeaks()  │  ← Sample-wise assessment
    │  (per replicate)│     Width, Scores
    └─────────────────┘
         │
         ├─── Width Distribution ────┐
         └─── Score Metrics ─────────┘
         │
         ▼
    ┌─────────────────┐
    │  filterPeaks()  │  ← Per replicate filtering
    │  (per replicate)│
    └─────────────────┘
         │
         ▼
    Filtered Peaks (Per Replicate) ✓

PHASE 2: CONSENSUS PEAK CREATION (Per Condition)
═══════════════════════════════════════════════════════════════════
    Filtered Peaks (All Replicates)
         │
         ▼
    ┌──────────────────────────────┐
    │ findOverlapsOfPeaks()        │
    │ or IDRfilter()               │  ← Create consensus
    └──────────────────────────────┘
         │
         ▼
    Consensus Peaks (Per Condition) ✓

PHASE 3: ANNOTATION
═══════════════════════════════════════════════════════════════════
    Consensus Peaks
         │
         ▼
    ┌─────────────────┐
    │ Multi-Condition?│
    └─────────────────┘
         │
         ├─── NO ──► getGenomicAnnotation()
         │            (on consensus peaks)
         │
         └─── YES ──► Differential Analysis
                      (diffBind/DESeq2)
                      │
                      ▼
                   Annotate Differential Peaks
                      (getGenomicAnnotation)
         │
         ▼
    ┌──────────────────────────────┐
    │ annotateHierarchically()     │
    │  - factor_type specified     │
    │  - Multiple parallel         │
    │    strategies                │
    └──────────────────────────────┘
         │
         ▼
    Best Annotation per Peak ✓

PHASE 4: ENHANCER IDENTIFICATION
═══════════════════════════════════════════════════════════════════
    Annotated Peaks
         │
         ▼
    ┌─────────────────┐
    │  Hi-C Available?│
    └─────────────────┘
         │
         ├─── YES ──► findEnhancers(Hi-C data)
         │
         └─── NO  ──► ENCODE CRE annotations
         │
         ▼
    Enhancer-Gene Associations ✓

PHASE 5: DOWNSTREAM ANALYSIS
═══════════════════════════════════════════════════════════════════
    Annotated Peaks (Consensus or Differential)
         │
         ├───► Enrichment Analysis (LOLA + Peak Overlap)
         │
         ├───► Gene Set Enrichment 
         │    (getEnrichedGO/PATH/test_enrichment - MSigDB/Enrich)
         │
         └───► Motif Analysis
         │
         ▼
    Final Integrated Results ✓
```

## Workflow with Time Estimates

```
┌─────────────────────────────────────────────────────────────────┐
│  STEP                    │  FUNCTION              │  TIME      │
├───────────────────────────┼────────────────────────┼────────────┤
│  PER-REPLICATE:           │                        │            │
│  1. Import & Assess       │  assessPeaks()         │  2-5 min   │
│     (per replicate)       │  (sample-wise)         │            │
│  2. Filter                │  filterPeaks()         │  <1 min    │
│     (per replicate)       │  (per replicate)       │            │
│                           │                        │            │
│  PER-CONDITION:           │                        │            │
│  3. Consensus Peaks       │  findOverlapsOfPeaks() │  2-5 min   │
│     (per condition)       │  or IDRfilter()        │            │
│  4. Preliminary Annotation │  getGenomicAnnotation │  3-5 min   │
│     (on consensus)        │                        │            │
│  5. Hierarchical          │  annotateHierarchically│  5-10 min  │
│     (on consensus)        │                        │            │
│  6. Enhancer ID           │  findEnhancers()       │  2-5 min   │
│                           │                        │            │
│  MULTI-CONDITION ONLY:    │                        │            │
│  7. Differential          │  diffBind/DESeq2      │  5-10 min  │
│  8. Annotate Diff Peaks   │  getGenomicAnnotation │  3-5 min   │
│  9. Hierarchical          │  annotateHierarchically│  5-10 min  │
│     (on differential)     │                        │            │
│                           │                        │            │
│  DOWNSTREAM:              │                        │            │
│  10. Enrichment           │  LOLA/Overlap          │  3-5 min   │
│  11. Gene Set Enrichment  │  getEnrichedGO/PATH/   │  3-5 min   │
│                           │  test_enrichment       │            │
│  12. Motif                │  MEME/HOMER            │  10-30 min │
└───────────────────────────┴────────────────────────┴────────────┘

Total Estimated Time: 
- Single condition: 30-65 minutes
- Multi-condition: 45-90 minutes
(depending on dataset size and number of replicates)
```

