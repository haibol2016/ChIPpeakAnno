# Diagrammatic Explanation of `output` Parameter in `annotatePeakInBatch`

This document provides visual diagrams explaining how each `output` mode works in the `annotatePeakInBatch` function.

## Table of Contents

### Key Concepts
- [Reference Points](#reference-points)
- [Diagram Legend](#diagram-legend)

### Detailed Parameter Explanations
- [1. `PeakLocForDistance` - Reference Point in Peak](#1-peaklocfordistance---reference-point-in-peak)
  - [`PeakLocForDistance = "start"` (Default)](#peaklocfordistance--start-default)
  - [`PeakLocForDistance = "middle"` (Recommended)](#peaklocfordistance--middle-recommended-for-variable-width-peaks)
  - [`PeakLocForDistance = "end"`](#peaklocfordistance--end)
  - [`PeakLocForDistance = "endMinusStart"` (Strand-aware)](#peaklocfordistance--endminusstart-strand-aware)
- [2. `FeatureLocForDistance` - Reference Point in Feature](#2-featurelocfordistance---reference-point-in-feature)
  - [`FeatureLocForDistance = "TSS"` (Default, Recommended)](#featurelocfordistance--tss-default-recommended)
  - [`FeatureLocForDistance = "geneEnd"`](#featurelocfordistance--geneend)
  - [`FeatureLocForDistance = "start"`](#featurelocfordistance--start)
  - [`FeatureLocForDistance = "end"`](#featurelocfordistance--end)
  - [`FeatureLocForDistance = "middle"`](#featurelocfordistance--middle)
- [3. `maxgap` - Maximum Gap Allowed](#3-maxgap---maximum-gap-allowed)
  - [`maxgap = -1L` (Default - No Gap Allowed)](#maxgap--1l-default---no-gap-allowed)
  - [`maxgap = 0L` (Allow Adjacent)](#maxgap--0l-allow-adjacent)
  - [`maxgap = 1000L` (Allow 1kb Gap)](#maxgap--1000l-allow-1kb-gap)
  - [`maxgap = 5000L` (Allow 5kb Gap)](#maxgap--5000l-allow-5kb-gap)
- [4. `bindingRegion` - Region Expansion Offsets](#4-bindingregion---region-expansion-offsets)
- [5. `bindingType` - Region Definition Strategy](#5-bindingtype---region-definition-strategy)
  - [`bindingType = "startSite"` (TSS-based)](#bindingtype--startsite-tss-based)
  - [`bindingType = "endSite"` (Gene End-based)](#bindingtype--endsite-gene-end-based)
  - [`bindingType = "fullRange"` (Full Gene Range)](#bindingtype--fullrange-full-gene-range)
  - [`bindingType = "nearestBiDirectionalPromoters"` (Bidirectional Promoters)](#bindingtype--nearestbidirectionalpromoters-bidirectional-promoters)

### Output Mode Diagrams
- [1. `output = "nearestLocation"` (Default)](#1-output--nearestlocation-default)
- [2. `output = "overlapping"`](#2-output--overlapping)
- [3. `output = "both"`](#3-output--both)
- [4. `output = "shortestDistance"`](#4-output--shortestdistance)
- [5. `output = "inside"`](#5-output--inside)
- [6. `output = "upstream&inside"`](#6-output--upstreaminside)
- [7. `output = "inside&downstream"`](#7-output--insidedownstream)
- [8. `output = "upstream"`](#8-output--upstream)
- [9. `output = "downstream"`](#9-output--downstream)
- [10. `output = "upstreamORdownstream"`](#10-output--upstreamordownstream)
- [11. `output = "upstream2downstream"`](#11-output--upstream2downstream)
- [12. `output = "nearestBiDirectionalPromoters"`](#12-output--nearestbidirectionalpromoters)

### Combined Examples
- [Example 1: `nearestLocation` with Reference Points](#example-1-nearestlocation-with-reference-points)
- [Example 2: `overlapping` with maxgap](#example-2-overlapping-with-maxgap)
- [Example 3: `upstream` with maxgap and FeatureLocForDistance](#example-3-upstream-with-maxgap-and-featurelocfordistance)

### Parameter Interactions
- [When `maxgap` is Used](#when-maxgap-is-used)
- [When Reference Points are Used](#when-reference-points-are-used)

### Quick Reference
- [Quick Reference Table](#quick-reference-table)

### Example Scenarios
- [Scenario 1: Standard Promoter Analysis](#scenario-1-standard-promoter-analysis)
- [Scenario 2: Find Nearest Gene](#scenario-2-find-nearest-gene)
- [Scenario 3: Upstream Enhancers](#scenario-3-upstream-enhancers)
- [Scenario 4: Bidirectional Promoters](#scenario-4-bidirectional-promoters)

### Additional Information
- [Notes](#notes)

---

## Key Concepts

### Reference Points
- **PeakLocForDistance**: Reference point in the peak (start, middle, end, or endMinusStart)
- **FeatureLocForDistance**: Reference point in the feature (TSS, geneEnd, middle, start, or end)
- **maxgap**: Maximum gap (in bp) allowed between peak and feature boundaries

### Diagram Legend
```
[+++++]  = Peak region
[==========]  = Feature/Gene region
  ↑      = PeakLocForDistance (reference point in peak)
  ↑      = FeatureLocForDistance (reference point in feature, e.g., TSS)
  |      = Distance/gap between regions
```

---

## Detailed Parameter Explanations

### 1. `PeakLocForDistance` - Reference Point in Peak

**Purpose**: Defines which point within each peak is used as the reference for distance calculations.

**Options**:

#### `PeakLocForDistance = "start"` (Default)
```
[+++++]
↑
start coordinate (leftmost/smallest coordinate)
```
- Uses the start coordinate of the peak
- For positive strand: leftmost position
- For negative strand: still leftmost position (coordinate-based, not strand-aware)

#### `PeakLocForDistance = "middle"` (Recommended for variable-width peaks)
```
[+++++]
   ↑
  middle = (start + end) / 2
```
- Uses the midpoint of the peak
- Calculated as: `(start + end) / 2`
- **Recommended** for peaks with variable widths (e.g., ChIP-seq peaks)
- More robust than "start" for wide peaks

#### `PeakLocForDistance = "end"`
```
[+++++]
     ↑
     end coordinate (rightmost/largest coordinate)
```
- Uses the end coordinate of the peak
- For positive strand: rightmost position
- For negative strand: still rightmost position (coordinate-based)

#### `PeakLocForDistance = "endMinusStart"` (Strand-aware)
```
When regarding to positive strand feature:
[+++++]
     ↑   [5'=============3']
     Uses peak end

When regarding to negative strand feature:
                     [+++++]
[3'=============5']   ↑
                      Uses peak start
```
- **Strand-aware**: Uses different reference points based on feature strand
- For features on **positive strand**: Uses peak **end**
- For features on **negative strand**: Uses peak **start**
- Useful for strand-specific assays (e.g., GRO-seq, PRO-seq)

**Visual Comparison**:
```
Peak:       [+++++++++++++] (start=100, end=300, middle=200)

"start":    ↑ (position 100)
"middle":          ↑ (position 200)
"end":                   ↑ (position 300)
"endMinusStart":         ↑ (if feature is + strand)
                          or
            ↑ (if feature is - strand)
```

**When is `FeatureLocForDistance` used or NOT used?**
- Used in `output = "nearestLocation"` for distance calculation
- Used in `output = "both"` (for the nearestLocation part)
- Used when calculating `distancetoFeature` in output metadata
- **NOT used** in `output = "overlapping"` (uses boundaries instead)
- **NOT used** in `output = "shortestDistance"` (uses boundaries instead)
- **NOT used** in `output = "inside"` (uses boundary containment instead)

---

### 2. `FeatureLocForDistance` - Reference Point in Feature

**Purpose**: Defines which point within each feature/gene is used as the reference for distance calculations.

**Options**:

#### `FeatureLocForDistance = "TSS"` (Default, Recommended)
```
Positive strand feature:
[5'==========3']
↑
TSS = start (transcription start site)

Negative strand feature:
[3'==========5']
                ↑
                TSS = end (transcription start site)
```
- **Strand-aware**: Uses the transcription start site
- For **positive strand**: TSS = `start` coordinate
- For **negative strand**: TSS = `end` coordinate
- **Standard choice** for most transcription factor and chromatin modifier analyses
- Most regulatory elements are defined relative to TSS

#### `FeatureLocForDistance = "geneEnd"`
```
Positive strand feature:
[5'==========3']
            ↑
            geneEnd = end (transcription end site)

Negative strand feature:
[3'==========5']
↑
geneEnd = start (transcription end site)
```
- **Strand-aware**: Uses the transcription end site (3' end)
- For **positive strand**: geneEnd = `end` coordinate
- For **negative strand**: geneEnd = `start` coordinate
- Useful for 3' regulatory element analysis (polyadenylation sites, 3' UTRs)
- Required for `output = "inside&downstream"`

#### `FeatureLocForDistance = "start"`
```
[5'==========3']
↑
start coordinate (always leftmost, NOT strand-aware)
```
- Always uses the **start** coordinate (leftmost position)
- **NOT strand-aware**: Same for both + and - strand features
- For positive strand: same as TSS
- For negative strand: different from TSS (TSS would be end)

#### `FeatureLocForDistance = "end"`
```
[5'==========3']
            ↑
            end coordinate (always rightmost, NOT strand-aware)
```
- Always uses the **end** coordinate (rightmost position)
- **NOT strand-aware**: Same for both + and - strand features
- For positive strand: same as geneEnd
- For negative strand: different from geneEnd (geneEnd would be start)

#### `FeatureLocForDistance = "middle"`
```
[5'==========3']
        ↑
        middle = (start + end) / 2
```
- Uses the midpoint of the feature
- Calculated as: `(start + end) / 2`
- **NOT strand-aware**: Same calculation for both strands
- Useful for gene body analysis

**Visual Comparison for Positive Strand Gene**:
```
Gene: [5'======================3'] (start=1000, end=5000, middle=3000, TSS=1000, geneEnd=5000)

"TSS":      ↑ (position 1000) - strand-aware
"geneEnd":                    ↑ (position 5000) - strand-aware
"start":    ↑ (position 1000) - NOT strand-aware
"end":                        ↑ (position 5000) - NOT strand-aware
"middle":          ↑ (position 3000) - NOT strand-aware
```

**When is `FeatureLocForDistance used or NOT used`?**
- Used in `output = "nearestLocation"` for distance calculation
- Used in `output = "both"` (for the nearestLocation part)
- Used in `output = "upstream"`, `"downstream"`, `"upstream&inside"`, `"inside&downstream"` to define regions
- Used when calculating `distancetoFeature` in output metadata
- Used to determine `bindingType` when `bindingRegion` triggers `annoPeaks()`
- **NOT used** in `output = "overlapping"` (uses boundaries instead)
- **NOT used** in `output = "shortestDistance"` (uses boundaries instead)
- **NOT used** in `output = "inside"` (uses boundary containment instead)

---

### 3. `maxgap` - Maximum Gap Allowed

**Purpose**: Defines the maximum gap (in base pairs) allowed between peak and feature boundaries for them to be considered as "overlapping".

**Definition of Gap**:
```
Case 1: Overlapping regions (gap = 0)
    [+++++]  Peak
    [==========] Feature
    Gap = 0 (they overlap)

Case 2: Adjacent regions (gap = 0)
    [+++++][==========]
    Peak   Feature
    Gap = 0 (touching, no space between)

Case 3: Disjoint regions (gap > 0)
    [+++++]     [==========]
     Peak       Feature
          |<gap>|
     Gap = distance between nearest boundaries
```

**How Gap is Calculated**:
```
For two disjoint intervals:
- Interval A: [start_A, end_A]
- Interval B: [start_B, end_B]

Gap = max(start_A, start_B) - min(end_A, end_B)

Example:
  Peak: [100, 200]
  Feature: [250, 500]
  Gap = max(100, 250) - min(200, 500) = 250 - 200 = 50 bp
```

**Parameter Values**:

#### `maxgap = -1L` (Default - No Gap Allowed)
```
    [+++++] Peak
    [==========] Feature
    ✓ Returns: Only if they directly overlap (gap = 0)

    [+++++]  [==========]
     Peak     Feature
     ✗ Does not return: Gap > 0
```
- **Strict overlap only**: Peaks must directly overlap features
- No gap allowed between boundaries
- Most restrictive option

#### `maxgap = 0L` (Allow Adjacent)
```
    [+++++] Peak
    [==========] Feature
    ✓ Returns: Overlapping (gap = 0)

    [+++++][==========]
    Peak    Feature
    ✓ Returns: Adjacent (gap = 0)

    [+++++]  [==========]
     Peak     Feature
     ✗ Does not return: Gap = 50 > 0
```
- Allows overlapping OR adjacent intervals
- Gap must be exactly 0

#### `maxgap = 1000L` (Allow 1kb Gap)
```
    [+++++]  [==========]
     Peak     Feature
     |<500bp>|
     ✓ Returns: Gap (500) ≤ maxgap (1000)

    [+++++]        [==========]
     Peak           Feature
     |<--1500bp-->|
     ✗ Does not return: Gap (1500) > maxgap (1000)
```
- Allows peaks within 1kb of features
- Useful for promoter-proximal analysis
- Common values: 500-2000 bp

#### `maxgap = 5000L` (Allow 5kb Gap)
```
    [+++++]              [==========]
     Peak                 Feature
          |<---3000bp--->|
     ✓ Returns: Gap (3000) ≤ maxgap (5000)
```
- Allows peaks within 5kb of features
- Useful for enhancer/distal regulatory element analysis
- Common values: 2000-10000 bp

**Visual Examples with Different maxgap Values**:
```
Peak: [100, 200]
Feature: [400, 600]
Gap = 400 - 200 = 200 bp

maxgap = -1:  ✗ Does not return (gap = 200 > 0)
maxgap = 0:   ✗ Does not return (gap = 200 > 0)
maxgap = 100: ✗ Does not return (gap = 200 > 100)
maxgap = 200: ✓ Returns (gap = 200 ≤ 200)
maxgap = 500: ✓ Returns (gap = 200 ≤ 500)
```

**When is `maxgap` used?** See the [table of when `maxgap` is used](#when-maxgap-is-used) in the Parameter Interactions section for a complete breakdown of which output modes use `maxgap` and how.

**Important Notes**:
1. **Gap vs Distance**: 
   - Gap = distance between **boundaries** (for overlap detection)
   - Distance = distance between **reference points** (for nearestLocation)

2. **When `bindingRegion` is set**: `maxgap` is **ignored** - the function uses `annoPeaks()` which uses `bindingRegion` instead

3. **Common Use Cases**:
   - `maxgap = 0`: Strict overlap only (most common)
   - `maxgap = 2000`: Promoter analysis (2kb upstream/downstream)
   - `maxgap = 5000`: Enhancer analysis (5kb regions)
   - `maxgap = 10000`: Distal regulatory elements (10kb regions)

---

### 4. `bindingRegion` - Region Expansion Offsets

**Purpose**: Defines a region around features by specifying upstream and downstream offsets from a reference point. When provided, triggers the `annoPeaks()` annotation method (region-based) instead of the internal point-based method.

**Format**: A vector with two integer values: `c(upstream_offset, downstream_offset)`

**Requirements**:
- First value: Upstream offset (must be ≤ 0, typically negative)
- Second value: Downstream offset (must be ≥ 1, typically positive)
- Example: `c(-5000, 3000)` means 5kb upstream and 3kb downstream

**How It Works**:
```
Reference Point (TSS or geneEnd)
        ↑
        |
[upstream_offset]  [feature reference]  [downstream_offset]
        |               ↑              |
        |<--negative-->||<--positive-->|
        
Example: bindingRegion = c(-5000, 3000)
         Creates region: [TSS - 5000] to [TSS + 3000]
```

**Visual Example for Positive Strand Gene**:
```
Original Feature:
[5'==========3'] (start=10000, end=15000, TSS=10000)

bindingRegion = c(-2000, 500):
           [5'==========3']
            ↑
           TSS
|<-2000bp->|               |<500bp>|
[Expanded Region: 8000 to 10500]
   [++++++] Peak
[5'-----------==========-----3']          
    Expanded annotation region
    (peaks overlapping this region will be annotated)
```

**Visual Example for Negative Strand Gene**:
```
Original Feature:
[3'==========5'] (start=10000, end=15000, TSS=15000)

bindingRegion = c(-2000, 500):
               [3'==========5']
                            ↑
                           TSS
        |<500bp>|              |<--2000bp-->|
        [Expanded Region: 14500 to 17000]
                                 [+++++++] Peak
        [3'--------==========---------------5'] 
           Expanded annotation region
```

**When `bindingRegion` Triggers `annoPeaks()`**:
```
Conditions:
1. bindingRegion length > 1 (i.e., has 2 values)
2. bindingType can be determined (either provided or auto-determined)

When BOTH conditions are met:
  → Uses annoPeaks() function (Method 1: region-based)
  → maxgap is IGNORED
  → Annotation regions are expanded by bindingRegion

When conditions NOT met:
  → Uses internal logic (Method 2: point-based)
  → bindingRegion is IGNORED
  → maxgap is used instead
```

**Common Use Cases**:
```
1. Promoter analysis:
   bindingRegion = c(-2000, 500)
   bindType = "startSite"
   → 2kb upstream to 500bp downstream of TSS

2. Full gene range with flanking:
   bindingRegion = c(-5000, 5000)
   bindingType = "fullRange"
   → 5kb upstream to 5kb downstream of entire gene
```

**Important Notes**:
1. **Method Selection**: When `bindingRegion` is provided and `bindingType` is determined and valid, annotation switches to `annoPeaks()` which uses region expansion + overlap detection, not point-to-point distances.

2. **Auto-determination of bindingType**: If not explicitly provided, `bindingType` is determined from:
   - `output = "nearestBiDirectionalPromoters"` → `bindingType = "nearestBiDirectionalPromoters"`
   - `output = "overlapping"` and `FeatureLocForDistance = "TSS"` → `bindingType = "startSite"`
   - `output = "overlapping"` and `FeatureLocForDistance = "geneEnd"` → `bindingType = "endSite"`

3. **Interaction with other parameters**:
   - When `bindingRegion` is set and `bindingType` is determined: `maxgap` is ignored
   - `output` must be "overlapping" or "nearestBiDirectionalPromoters" for auto-determination; 
   `bindingType` should be explicitly defined for other output modes to trigger `annoPeaks()` for peak annotation.

---

### 5. `bindingType` - Region Definition Strategy

**Purpose**: Specifies how the binding region is defined relative to features. Determines which reference point is used and how region expansion is constrained. Only used when `bindingRegion` triggers `annoPeaks()`.

**Options**:

#### `bindingType = "startSite"` (TSS-based)

**Description**: Defines binding region relative to the feature start site (TSS for positive strand, gene end for negative strand). Annotation regions are expanded by `bindingRegion`, but downstream expansion is **constrained to not exceed original feature boundaries**.

**Visual Example (Positive Strand)**:
```
Original Feature:
               [5'=======================3'] (start=10000, end=15000, TSS=10000)

bindingRegion = c(-2000, 5000):
               [5'======================3']
               ↑
              TSS
  |<--2000bp-->|<------5000bp---------->|
    [Expanded Region: 8000 to 15000]
[+++++++]
  [5'----------------------------------3']  Expanded annotation region
               ↑ 
               TSS
    (downstream expansion stops at gene end)
```

**Visual Example (Negative Strand)**:
```
Original Feature:
            [3'=======================5'] (start=10000, end=15000, TSS=15000)

bindingRegion = c(-2000, 5000):
            [3'=======================5']
                                     ↑
                                     TSS
            |<------------5000bp---->|<--2000bp--->|
            [5'----------------------------------3']  Expanded annotation region 
                [Expanded Region: 10000 to 17000]
                                              [+++++++] Peak
           (downstream expansion stops at gene start)
```

**Key Characteristics**:
- Uses TSS as reference point (strand-aware)
- Downstream expansion constrained by original feature boundaries
- Best for promoter-proximal binding analysis
- Auto-determined when: `output = "overlapping"` + `FeatureLocForDistance = "TSS"`

---

#### `bindingType = "endSite"` (Gene End-based)

**Description**: Defines binding region relative to the feature end site (gene end for positive strand, TSS for negative strand). Annotation regions are expanded by `bindingRegion`, but upstream expansion is constrained to not exceed original feature boundaries.

**Visual Example (Positive Strand)**:
```
Original Feature:
                       [5'===============3'] (start=10000, end=15000, geneEnd=15000)

bindingRegion = c(-5000, 3000):
                       [5'===============3']
                                        ↑
                                      geneEnd
                         |<-- 5000bp--->|<-3000bp-->|
                       [Expanded Region: 10000 to 18000]
                      [5'----------------------------3']
                                                 [+++++++] Peak
                                        ↑
                                      geneEnd         
           (upstream expansion stops at gene start)
```

**Visual Example (Negative Strand)**:
```

Original Feature:
                       [3'===============5'] (start=10000, end=15000, geneEnd=15000)

bindingRegion = c(-5000, 3000):
                       [3'===============5']
                          ↑
                        geneEnd
           |<-- 3000bp--->|<----5000bp-->|
          [Expanded Region: 10000 to 18000]
          [3'----------------------------5']
        [+++++++] Peak
                          ↑
                        geneEnd         
           (upstream expansion stops at gene start)
```

**Key Characteristics**:
- Uses gene end as reference point (strand-aware)
- Upstream expansion constrained by original feature boundaries
- Best for 3' end regulatory element analysis
- Auto-determined when: `output = "overlapping"` + `FeatureLocForDistance = "geneEnd"`

---

#### `bindingType = "fullRange"` (Full Gene Range)

**Description**: Uses the entire feature range. Annotation regions are expanded by `bindingRegion` in both directions **without constraints**. This allows expansion beyond original feature boundaries.

**Visual Example**:
```
Original Feature:
                [5'==========3'] (start=10000, end=15000)

bindingRegion = c(-5000, 5000):
                [5'==========3']
     |<--5000bp-->|         |<5000bp-->|
     |<-----------===========--------->|
       [Expanded Region: 5000 to 20000]
  [+++++++]                  [+++++++]
```

**Key Characteristics**:
- Uses entire feature range as base
- Expansion in both directions without constraints
- Best for gene body analysis or when you need flanking regions
- Must be explicitly provided (not auto-determined)

**Use Cases**:
- Histone marks that span gene bodies (H3K36me3, H3K27me3)
- Full transcript analysis with flanking regions
- When you need symmetric expansion around entire gene

---

#### `bindingType = "nearestBiDirectionalPromoters"` (Bidirectional Promoters)

**Description**: Identifies peaks near bidirectional promoters using a two-step algorithm.

**Visual Example**:
```
Gene Pair (Divergent):
                  [++++]  [5'==========3']
[3'==========5'] 
   Gene1 (rev)     peak      Gene2 (forward)
        ↑                      ↑
       TSS1                   TSS2
        
Bidirectional Promoter Region:
                        [5'==========3']
[3'==========5']                     
      Gene1                 Gene2
           [<--promoter region-->]
                  [++++] Peak
peak overlaps with bidirectional promoter region
✓ Returns annotations for BOTH Gene1 and Gene2
```

**Key Characteristics**:
- Reports annotations for **both genes** in bidirectional promoter pairs
- `select = "bestOne"` is not supported (always keeps both genes)
- Auto-determined when: `output = "nearestBiDirectionalPromoters"`

**Algorithm Steps**:
1. Create promoter regions using `promoters()` with `bindingRegion` as upstream/downstream distances
2. Find peaks overlapping with the expaned promoter regions
3. Return annotations for both genes if bidirectional promoters exist for a peak; otherwise return the feature whose expanded promoter overlap with peak. 

---

**Summary Table: `bindingType` Options**

| bindingType | Reference Point | Expansion Constraints | Best For | Auto-Determined? |
|------------|----------------|----------------------|----------|------------------|
| `"startSite"` | TSS (strand-aware) | Downstream constrained by gene end | Promoter analysis | Yes (when `output="overlapping"` + `FeatureLocForDistance="TSS"`) |
| `"endSite"` | Gene end (strand-aware) | Upstream constrained by gene start | 3' end analysis | Yes (when `output="overlapping"` + `FeatureLocForDistance="geneEnd"`) |
| `"fullRange"` | Entire feature range | No constraints | Gene body + flanking | No (must be explicit) |
| `"nearestBiDirectionalPromoters"` | TSS (both directions) | Uses `maxTSSDistance` instead | Bidirectional promoters | Yes (when `output="nearestBiDirectionalPromoters"`) |

---

**How `bindingType` Affects Region Expansion**:

```
bindingType = "startSite":
  Annotation regions expanded by bindingRegion
  Downstream expansion constrained by original feature boundaries
  
bindingType = "endSite":
  Annotation regions expanded by bindingRegion
  Upstream expansion constrained by original feature boundaries
  
bindingType = "fullRange":
  Annotation regions expanded by bindingRegion
  NO constraints (can extend beyond feature boundaries)
  
bindingType = "nearestBiDirectionalPromoters":
  Uses promoters() function with bindingRegion
  bindingRegion parameter is ignored for distance filtering
  Uses maxTSSDistance instead
  
```

---

**Interaction with Other Parameters**:

```
When bindingRegion triggers annoPeaks():

1. bindingType determines:
   - Which reference point to use (TSS, geneEnd, or full range)
   - How expansion is constrained
   - Whether to use bidirectional promoter algorithm

2. bindingRegion determines:
   - Upstream/downstream expansion distances

3. FeatureLocForDistance:
   - Used to auto-determine bindingType when not explicit

4. maxgap:
   - IGNORED (bindingRegion is used instead)

5. output:
   - Must be "overlapping" or "nearestBiDirectionalPromoters" for auto-determination
   - Determines which bindingType is auto-selected
```

---

## Output Mode Diagrams

### 1. `output = "nearestLocation"` (Default)

**Description**: Returns the nearest feature based on absolute distance between reference points.

**Formula**: `|PeakLocForDistance - FeatureLocForDistance|`

**Diagram**:
```
Case 1: Peak upstream of feature
    [+++++]  [=====]
       ↑        ↑
       |<--d--->|
     Peak   Feature (TSS)
     ✓ Returns this feature (minimum distance between `peakLocForDistance` and `featureLocforDistance`)

Case 2: Peak overlapping feature
    [+++++] Peak
      [==========]  Feature
       ↑    ↑
       |<d->|
 
     ✓ Returns this feature (distance = 0 or small)

Case 3: Peak downstream of feature
    [+++++]  [=====]
       ↑     ↑
       |<-d->|
             Feature Peak
             ✓ Returns this feature (minimum distance)
```

**Key Points**:
- Uses reference points (PeakLocForDistance and FeatureLocForDistance)
- Returns both overlapping and non-overlapping features
- Always returns the feature with minimum distance

---

### 2. `output = "overlapping"`

**Description**: Returns all features that overlap with peaks, where overlap is defined as gap ≤ maxgap.

**Diagram**:
```
Case 1: Direct overlap (gap = 0)
    [+++++]   Peak 
    [==========] Feature
    ✓ Returns feature (gap = 0 ≤ maxgap)

Case 2: Gap within maxgap
    [+++++]      [=====]
     Peak         Feature
          | <gap>|
     If gap ≤ maxgap: ✓ Returns feature
     If gap > maxgap: ✗ Does not return

Case 3: Peak contained in feature
    [==========] Feature
      [+++++]    Peak
    ✓ Returns feature (gap = 0)

Case 4: Feature contained in peak
    [+++++] Peak
    [==] Feature
    Feature
    ✓ Returns feature (gap = 0)
```

**Key Points**:
- Gap is calculated between peak and feature boundaries (not reference points)
- Uses `maxgap` parameter to allow nearby features
- May return zero, one, or multiple features per peak

---

### 3. `output = "both"`

**Description**: Combines `nearestLocation` and `overlapping` results.

**Diagram**:
```
    [=====][+++++]         
               [=====]
          ↑↑         
     Feature1   Feature2   
     
    Nearest: Feature1 (minimum distance from reference points)
    Overlapping: Feature2 (overlaps with peak, gap ≤ maxgap)
    
    Result: Returns BOTH Feature1 and Feature2
```

**Key Points**:
- Returns nearest feature(s) based on reference points
- Plus any overlapping features (gap ≤ maxgap)
- Most comprehensive annotation option

---

### 4. `output = "shortestDistance"`

**Description**: Returns features with minimum boundary-to-boundary distance.

**Diagram**:
```
Case 1: Peak upstream
   [==========]           [+++++]  [=====]   [=========]
    FeatureC                Peak    FeatureA  FeatureB
                                | d|
     Distance = end(peak) to start(feature)
     ✓ Returns FeatureA if this is the shortest distance

Case 2: Peak downstream
    [=====]  [+++++]     [=========]
    FeatureA   Peak        FeatureB
          |d |
             Distance = end(feature) to start(peak)
             ✓ Returns FeatureA if this is the shortest distance

Case 3: Overlapping
    [+++++]  Peak
    [==========] Feature
    Distance = 0 (overlapping)
    ✓ Returns feature
```

**Key Points**:
- Uses `nearest()` from GenomicRanges
- Considers ALL boundary-to-boundary distances (not reference points)
- Returns both overlapping and non-overlapping features with minimum distance

---

### 5. `output = "inside"`

**Description**: Returns features where peaks are completely contained within (inside) the feature region. Uses `type = "within"` in `findOverlaps()`, which requires the peak to be fully inside the feature boundaries.

**Diagram**:
```
Case 1: Peak completely inside feature
    [==========] Feature
      [+++++]    Peak
    ✓ Returns feature (peak is completely within feature)

Case 2: Peak overlapping but not completely inside
    [+++++]  Peak
    [==========] Feature
    ✗ Does not return (peak extends beyond feature boundaries)

Case 3: Peak outside feature
    [+++++]  [==========]
     Peak     Feature
    ✗ Does not return (peak is outside feature)
```

**Key Points**:
- Uses `findOverlaps()` with `type = "within"`
- Peak must be completely contained within feature boundaries
- Does not use `maxgap` parameter (strict containment only)
- Does not use reference points (uses boundary containment)
- Useful for identifying peaks that are fully within gene bodies or exons

---

### 6. `output = "upstream&inside"`

**Description**: Returns features where peaks are upstream of TSS (within maxgap) OR overlapping with feature.

**Diagram**:
```
    [+++++]     [5'==========3']
     Peak         Feature
          |<maxgap>|
     Region: [upstream by maxgap] + [entire feature]
     
    Case 1: Peak upstream within maxgap
        [+++++]  [5'==========3']
         Peak    Feature
         ✓ Returns feature
    
    Case 2: Peak overlapping feature
        [+++++]  Peak 
        [5'==========3'] Feature
        ✓ Returns feature
    
    Case 3: Peak downstream beyond feature
        [+++++]  [5'==========3']  [+++++]
                Feature      Peak
                ✗ Does not return
```

**Key Points**:
- Defines region: `[TSS - maxgap]` to `[feature end]`
- Returns features if peak overlaps this extended region
- Combines upstream and gene body annotations

---

### 7. `output = "inside&downstream"`

**Description**: Returns features where peaks are downstream of gene end (within maxgap) OR overlapping with feature.

**Requirement**: `FeatureLocForDistance = "geneEnd"`

**Diagram**:
```
    [5'==========3']      [+++++]
    Feature              Peak
               |<maxgap>|
    Region: [feature start] to [geneEnd + maxgap]
    
    Case 1: Peak inside feature
        [5'==========3'] Feature
          [+++++] Peak
        ✓ Returns feature
    
    Case 2: Peak downstream within maxgap
        [5'==========3']  [+++++]
        Feature            Peak
        ✓ Returns feature
    
    Case 3: Peak upstream of feature
        [+++++]  [5'==========3']
         Peak    Feature
         ✗ Does not return
```

**Key Points**:
- Requires `FeatureLocForDistance = "geneEnd"`
- Defines region: `[feature start]` to `[geneEnd + maxgap]`
- Captures 3' regulatory elements and gene body

---

### 8. `output = "upstream"`

**Description**: Returns features where peaks are upstream of TSS (within maxgap).

**Diagram**:
```
    [+++++]        [5'==========3']
     Peak             Feature
          |<maxgap>|
     
    Region: Upstream of TSS, within maxgap distance
    
    Case 1: Peak upstream within maxgap
        [+++++]  [5'==========3']
         Peak    Feature
         ✓ Returns feature
    
    Case 2: Peak overlapping or downstream
        [+++++] Peak
        [5'==========3'] Feature
        ✗ Does not return
    
    Case 3: Peak too far upstream
   [+++++]             [5'==========3']
     Peak                 Feature
         |<--maxgap+-->|
         ✗ Does not return (gap > maxgap)
```

**Key Points**:
- Only returns features where peak is upstream of TSS
- Distance must be ≤ maxgap
- Does not return overlapping or downstream features

---

### 9. `output = "downstream"`

**Description**: Returns features where peaks are downstream of gene end (within maxgap).

**Diagram**:
```
    [5'==========3']  [+++++]
    Feature            Peak
                  |<maxgap>|
    
    Region: Downstream of gene end, within maxgap distance
    
    Case 1: Peak downstream within maxgap
        [5'==========3']  [+++++]
        Feature            Peak
        ✓ Returns feature
    
    Case 2: Peak overlapping or upstream
        [+++++]  Peak
        [==========]  Feature

        ✗ Does not return
    
    Case 3: Peak too far downstream
        [==========]              [+++++]
        Feature                  Peak
                   |<--maxgap+-->|
        ✗ Does not return (gap > maxgap)
```

**Key Points**:
- Only returns features where peak is downstream of gene end
- Distance must be ≤ maxgap
- Does not return overlapping or upstream features

---

### 10. `output = "upstreamORdownstream"`

**Description**: Returns features that are either upstream of TSS OR downstream of gene end (both within maxgap).

**Diagram**:
```
    [+++++]  [==========]  [+++++]
     Peak1   Feature      Peak2
     |<maxgap>|            |<maxgap>|
    
    Upstream region: [TSS - maxgap] to [TSS]
    Downstream region: [geneEnd] to [geneEnd + maxgap]
    
    Case 1: Peak1 upstream
        [+++++]  [==========]
         Peak1   Feature
         ✓ Returns feature
    
    Case 2: Peak2 downstream
        [==========]  [+++++]
        Feature       Peak2
        ✓ Returns feature
    
    Case 3: Peak inside feature
        [==========]
          [+++++]
          Peak
        ✗ Does not return (not upstream or downstream)
```

**Key Points**:
- Returns features in either direction (upstream OR downstream)
- Both directions use maxgap threshold
- Does not return features where peak is inside the gene body

---

### 11. `output = "upstream2downstream"`

**Description**: Returns features where peaks are anywhere between upstream of TSS (within maxgap) and downstream of gene end (within maxgap). This creates a single expanded region that includes the feature plus maxgap upstream and maxgap downstream.

**Diagram**:
```
    [+++++]  [==========]  [+++++]
     Peak1   Feature       Peak2
    |<maxgap>|          |<maxgap>|
    
    Expanded region: [feature start - maxgap] to [feature end + maxgap]
    
    Case 1: Peak1 upstream within maxgap
        [+++++]  [==========]
         Peak1   Feature
         ✓ Returns feature
    
    Case 2: Peak2 downstream within maxgap
        [==========]  [+++++]
        Feature       Peak2
        ✓ Returns feature
    
    Case 3: Peak inside feature
        [==========]
          [+++++]
          Peak
        ✓ Returns feature
    
    Case 4: Peak too far upstream
   [+++++]              [==========]
     Peak                 Feature
            |<--maxga-->|
     ✗ Does not return (gap > maxgap)
```

**Key Points**:
- Creates a single expanded region: `[feature start - maxgap]` to `[feature end + maxgap]`
- Returns features if peak overlaps this expanded region
- Captures peaks upstream, inside, and downstream of features (within maxgap)
- More inclusive than `upstreamORdownstream` (which excludes peaks inside gene body)

---

### 12. `output = "nearestBiDirectionalPromoters"`

**Description**: Uses `annoPeaks()` to find bidirectional promoters. Requires `bindingRegion` 
and determinable `bindingType`. Note: here the definition of bidirectional promoters is
different from the literature's.

**Diagram**:
```
    [3'=========5']  [++++]   [5'==========5']
     Gene1           peak          Gene2
     (rev)          (forward)
             |<--bindingRegion-->|
     
    Region: [TSS1 - upstream_offset] to [TSS2 + downstream_offset]
    
```

**Key Points**:
- **Requires `bindingRegion` parameter** (e.g., `c(-5000, 3000)`)
- Uses `annoPeaks()` function internally
- Reports bidirectional promoters if found in both directions
- Otherwise reports closest promoter in one direction
- Different from literature definition of bidirectional promoters

---

## Combined Examples: How Parameters Work Together

### Example 1: `nearestLocation` with Reference Points
```
Peak: [1000, 1500]
Feature: [5'==========3'] (start=2000, end=5000, TSS=2000)

PeakLocForDistance = "middle" → Peak reference = 1250
FeatureLocForDistance = "TSS" → Feature reference = 2000

Distance = |1250 - 2000| = 750 bp

Result: ✓ Returns this feature (if it's the nearest)
```

### Example 2: `overlapping` with maxgap
```
Peak: [1000, 1500]
             Feature: [5'==========3'] (start=2000, end=5000)

Gap = max(1000, 2000) - min(1500, 5000) = 2000 - 1500 = 500 bp

maxgap = 1000: ✓ Returns (gap 500 ≤ maxgap 1000)
maxgap = 200:  ✗ Does not return (gap 500 > maxgap 200)
```

### Example 3: `upstream` with maxgap and FeatureLocForDistance
```
Peak: [500, 800]
                  Feature: [5'==========3'] (start=2000, end=5000, TSS=2000)

FeatureLocForDistance = "TSS" → TSS = 2000
Peak end = 800
Distance from peak end to TSS = 2000 - 800 = 1200 bp

maxgap = 2000: ✓ Returns (peak is upstream, distance 1200 ≤ maxgap 2000)
maxgap = 1000: ✗ Does not return (distance 1200 > maxgap 1000)
```

---

## Parameter Interactions

### When `maxgap` is Used

| Output Mode | Uses maxgap? | How it's used |
|------------|--------------|---------------|
| `nearestLocation` | ❌ No | Uses reference point distances only |
| `overlapping` | ✅ Yes | Gap between boundaries must be ≤ maxgap |
| `both` | ✅ Yes | For overlapping part only |
| `shortestDistance` | ❌ No | Uses `nearest()` function |
| `inside` | ❌ No | Uses `type = "within"` (strict containment) |
| `upstream&inside` | ✅ Yes | Defines upstream region extent |
| `inside&downstream` | ✅ Yes | Defines downstream region extent |
| `upstream` | ✅ Yes | Maximum distance upstream of TSS |
| `downstream` | ✅ Yes | Maximum distance downstream of gene end |
| `upstreamORdownstream` | ✅ Yes | Both directions use maxgap |
| `upstream2downstream` | ✅ Yes | Both directions use maxgap |
| `nearestBiDirectionalPromoters` | ❌ No | Uses `bindingRegion` instead |

### When Reference Points are Used for Overlap Detection

**Important Distinction:**

This table indicates when reference points are used for **overlap detection** (finding which features match). This is different from **distance calculation**:

- **Overlap Detection**: Determines which features are selected/returned. Some modes use reference points here (e.g., `nearestLocation` uses both, `upstream` uses only `FeatureLocForDistance` for TSS).
- **Distance Calculation**: After overlaps are found, `distancetoFeature` is **always calculated using BOTH** `PeakLocForDistance` and `FeatureLocForDistance` for all output modes (except those using `annoPeaks()`). This happens in STEP2 of the annotation process (lines 866-892 in the code), regardless of what was used for overlap detection.

**Example**: For `output = "upstream"`, only `FeatureLocForDistance` (TSS) is used to define the upstream region for overlap detection, but once overlaps are found, `distancetoFeature` is calculated using both the peak's reference point (`PeakLocForDistance`) and the feature's TSS (`FeatureLocForDistance`).

| Output Mode | Uses PeakLocForDistance? | Uses FeatureLocForDistance? |
|------------|-------------------------|----------------------------|
| `nearestLocation` | ✅ Yes | ✅ Yes (uses prepared reference points to find nearest) |
| `overlapping` | ❌ No | ❌ No (uses boundaries for overlap detection) |
| `both` | ✅ Yes (for nearest part) | ✅ Yes (for nearest part, uses prepared reference points) |
| `shortestDistance` | ❌ No | ❌ No (uses boundaries for nearest detection) |
| `inside` | ❌ No | ❌ No (uses boundary containment) |
| `upstream&inside` | ❌ No | ✅ Yes (uses TSS concept in strand-aware calculations) |
| `inside&downstream` | ❌ No | ✅ Yes (uses geneEnd concept in strand-aware calculations) |
| `upstream` | ❌ No | ✅ Yes (uses TSS concept in strand-aware calculations) |
| `downstream` | ❌ No | ✅ Yes (uses geneEnd concept in strand-aware calculations) |
| `upstreamORdownstream` | ❌ No | ✅ Yes (uses TSS and geneEnd concepts in strand-aware calculations) |
| `upstream2downstream` | ❌ No | ❌ No (uses boundaries, expands by maxgap) |
| `nearestBiDirectionalPromoters` | ❌ No | ✅ Yes (for bindingType determination when using annoPeaks) |

---

## Quick Reference Table

| Output Mode | Primary Use Case | Key Parameters |
|------------|------------------|----------------|
| `nearestLocation` | Find closest gene | PeakLocForDistance, FeatureLocForDistance |
| `overlapping` | Promoter/gene body analysis | maxgap |
| `both` | Comprehensive annotation | PeakLocForDistance, FeatureLocForDistance, maxgap |
| `shortestDistance` | Minimum boundary distance | (none - uses boundaries) |
| `inside` | Peaks completely within features | (none - uses boundary containment) |
| `upstream&inside` | Promoter + gene body | maxgap, FeatureLocForDistance |
| `inside&downstream` | Gene body + 3' region | maxgap, FeatureLocForDistance="geneEnd" |
| `upstream` | Upstream regulatory elements | maxgap, FeatureLocForDistance |
| `downstream` | 3' regulatory elements | maxgap, FeatureLocForDistance |
| `upstreamORdownstream` | Both directions from gene | maxgap, FeatureLocForDistance |
| `upstream2downstream` | Both directions from gene | maxgap, FeatureLocForDistance |
| `nearestBiDirectionalPromoters` | Bidirectional promoter detection | bindingRegion (required) |

---

## Example Scenarios

### Scenario 1: Standard Promoter Analysis
```r
annotatePeakInBatch(peaks, 
                   AnnotationData = genes,
                   output = "overlapping",
                   maxgap = 2000,
                   FeatureLocForDistance = "TSS")
```
**What it does**: Finds genes where peaks overlap or are within 2kb of the gene region (from 2kb upstream of TSS to gene end).

### Scenario 2: Find Nearest Gene
```r
annotatePeakInBatch(peaks,
                   AnnotationData = genes,
                   output = "nearestLocation",
                   PeakLocForDistance = "middle",
                   FeatureLocForDistance = "TSS")
```
**What it does**: Finds the gene with the TSS closest to the peak center, regardless of overlap.

### Scenario 3: Upstream Enhancers
```r
annotatePeakInBatch(peaks,
                   AnnotationData = genes,
                   output = "upstream",
                   maxgap = 5000,
                   FeatureLocForDistance = "TSS")
```
**What it does**: Finds genes where peaks are upstream of TSS within 5kb (does not include overlapping or downstream peaks).

### Scenario 4: Bidirectional Promoters
```r
annotatePeakInBatch(peaks,
                   AnnotationData = genes,
                   output = "nearestBiDirectionalPromoters",
                   bindingRegion = c(-5000, 3000))
```
**What it does**: Uses `annoPeaks()` to find bidirectional promoters within 5kb upstream and 3kb downstream of TSS.

---

## Notes

1. **Gap calculation**: For `overlapping` mode, gap is calculated between the nearest boundaries of peak and feature, not between reference points.

2. **Distance calculation**: For `nearestLocation` mode, distance is calculated between reference points (PeakLocForDistance and FeatureLocForDistance).

3. **Strand awareness**: TSS and geneEnd are strand-aware (TSS = start for + strand, end for - strand).

4. **Multiple features**: Most modes can return multiple features per peak. Use `select` parameter to control which one(s) to return.

5. **Method selection**: When `bindingRegion` is provided and valid, the function uses `annoPeaks()` (Method 1), otherwise uses internal logic (Method 2).

