# Diagrammatic Explanation of `output` Parameter in `annotatePeakInBatch`

This document provides visual diagrams explaining how each `output` mode works in the `annotatePeakInBatch` function.

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

**When is it used?**
- Used in `output = "nearestLocation"` for distance calculation
- Used in `output = "both"` (for the nearestLocation part)
- Used when calculating `distancetoFeature` in output metadata
- **NOT used** in `output = "overlapping"` (uses boundaries instead)
- **NOT used** in `output = "shortestDistance"` (uses boundaries instead)

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

**When is it used?**
- Used in `output = "nearestLocation"` for distance calculation
- Used in `output = "both"` (for the nearestLocation part)
- Used in `output = "upstream"`, `"downstream"`, `"upstream&inside"`, `"inside&downstream"` to define regions
- Used when calculating `distancetoFeature` in output metadata
- Used to determine `bindingType` when `bindingRegion` triggers `annoPeaks()`
- **NOT used** in `output = "overlapping"` (uses boundaries instead)
- **NOT used** in `output = "shortestDistance"` (uses boundaries instead)

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

**When is `maxgap` used?**
- ✅ Used in `output = "overlapping"` - gap between boundaries must be ≤ maxgap
- ✅ Used in `output = "both"` - for the overlapping part
- ✅ Used in `output = "upstream"` - maximum distance upstream of TSS
- ✅ Used in `output = "downstream"` - maximum distance downstream of gene end
- ✅ Used in `output = "upstream&inside"` - defines upstream region extent
- ✅ Used in `output = "inside&downstream"` - defines downstream region extent
- ✅ Used in `output = "upstreamORdownstream"` - both directions
- ✅ Used in `output = "upstream2downstream"` - both directions
- ❌ **NOT used** in `output = "nearestLocation"` - uses reference point distances only
- ❌ **NOT used** in `output = "shortestDistance"` - uses `nearest()` function
- ❌ **NOT used** in `output = "nearestBiDirectionalPromoters"` - uses `bindingRegion` instead

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

## Combined Example: How Parameters Work Together

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

### 5. `output = "upstream&inside"`

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

### 6. `output = "inside&downstream"`

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

### 7. `output = "upstream"`

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

### 8. `output = "downstream"`

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

### 9. `output = "upstreamORdownstream"`

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

### 10. `output = "upstream2downstream"`

**Description**: Returns features upstream of TSS, feature region 
and downstream of gene end (within maxgap).

**Diagram**:
```
    [+++++]  [==========]  [+++++]
     Peak1   Feature       Peak2
    |<maxgap>|          |<maxgap>|
    
```

---

### 11. `output = "nearestBiDirectionalPromoters"`

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

## Parameter Interactions

### When `maxgap` is Used

| Output Mode | Uses maxgap? | How it's used |
|------------|--------------|---------------|
| `nearestLocation` | ❌ No | Uses reference point distances only |
| `overlapping` | ✅ Yes | Gap between boundaries must be ≤ maxgap |
| `both` | ✅ Yes | For overlapping part only |
| `shortestDistance` | ❌ No | Uses `nearest()` function |
| `upstream&inside` | ✅ Yes | Defines upstream region extent |
| `inside&downstream` | ✅ Yes | Defines downstream region extent |
| `upstream` | ✅ Yes | Maximum distance upstream of TSS |
| `downstream` | ✅ Yes | Maximum distance downstream of gene end |
| `upstreamORdownstream` | ✅ Yes | Both directions use maxgap |
| `upstream2downstream` | ✅ Yes | Both directions use maxgap |
| `nearestBiDirectionalPromoters` | ❌ No | Uses `bindingRegion` instead |

### When Reference Points are Used

| Output Mode | Uses PeakLocForDistance? | Uses FeatureLocForDistance? |
|------------|-------------------------|----------------------------|
| `nearestLocation` | ✅ Yes | ✅ Yes |
| `overlapping` | ❌ No | ❌ No (uses boundaries) |
| `both` | ✅ Yes (for nearest) | ✅ Yes (for nearest) |
| `shortestDistance` | ❌ No | ❌ No (uses boundaries) |
| `upstream&inside` | ❌ No | ✅ Yes (for TSS) |
| `inside&downstream` | ❌ No | ✅ Yes (for geneEnd) |
| `upstream` | ❌ No | ✅ Yes (for TSS) |
| `downstream` | ❌ No | ✅ Yes (for geneEnd) |
| `upstreamORdownstream` | ❌ No | ✅ Yes (for TSS and geneEnd) |
| `upstream2downstream` | ❌ No | ✅ Yes (for TSS and geneEnd) |
| `nearestBiDirectionalPromoters` | ❌ No | ✅ Yes (for bindingType determination) |

---

## Quick Reference Table

| Output Mode | Primary Use Case | Key Parameters |
|------------|------------------|----------------|
| `nearestLocation` | Find closest gene | PeakLocForDistance, FeatureLocForDistance |
| `overlapping` | Promoter/gene body analysis | maxgap |
| `both` | Comprehensive annotation | PeakLocForDistance, FeatureLocForDistance, maxgap |
| `shortestDistance` | Minimum boundary distance | (none - uses boundaries) |
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

