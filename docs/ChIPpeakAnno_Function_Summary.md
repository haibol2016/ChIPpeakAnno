# ChIPpeakAnno Package Function Summary

## Package Overview

**ChIPpeakAnno** (version 3.45.2) is a Bioconductor package for batch annotation and analysis of peaks from ChIP-seq, ATAC-seq, and NAD-seq experiments. The package provides comprehensive tools for identifying the closest gene, exon, miRNA, or custom features, retrieving sequences around peaks, and obtaining enriched Gene Ontology (GO) or Pathway terms.

**Key Features:**
- Flexible peak annotation with multiple binding types and distance calculations
- Comprehensive statistical testing (hypergeometric tests, permutation tests)
- Sequence analysis and pattern matching capabilities
- Signal analysis and visualization tools
- Integration with Bioconductor ecosystem (TxDb, EnsDb, BSgenome, biomaRt)

---

## Function Categories and Descriptions

### 1. Peak Annotation Functions

#### `annoPeaks()`
**Purpose:** Annotates peaks with genomic features using flexible binding types.

**Key Parameters:**
- `peaks`: GRanges object containing peaks
- `annoData`: Annotation data (GRanges or annoGR object)
- `bindingType`: "nearestBiDirectionalPromoters", "startSite", "endSite", "fullRange"
- `bindingRegion`: Vector of two integers (upstream, downstream), default c(-5000, 5000)
- `ignore.peak.strand`: Logical, whether to ignore peak strand
- `select`: "all" or "bestOne" - return all annotations or best match

**Return Value:** GRanges object with annotated peaks including feature information, distances, and relationships

**Use Cases:**
- Annotate peaks to TSS with specific upstream/downstream regions
- Find bidirectional promoters
- Annotate to gene ends or full gene ranges

---

#### `annotatePeakInBatch()`
**Purpose:** Obtains distance to nearest TSS, miRNA, exon, or custom features for a list of peaks.

**Key Parameters:**
- `myPeakList`: GRanges object
- `AnnotationData`: GRanges or annoGR object (or mart object for biomaRt queries)
- `output`: "nearestLocation", "overlapping", "both", "shortestDistance", "upstream", "downstream", "nearestBiDirectionalPromoters"
- `PeakLocForDistance`: "start", "middle", "end", "endMinusStart"
- `FeatureLocForDistance`: "TSS", "middle", "start", "end", "geneEnd"
- `maxgap`: Maximum gap allowed between ranges
- `select`: "all", "first", "last", "arbitrary"

**Return Value:** GRanges object with annotation metadata including:
- `feature`: Feature ID
- `insideFeature`: Relationship (upstream, downstream, inside, overlapStart, overlapEnd, includeFeature)
- `distancetoFeature`: Distance to feature
- `shortestDistance`: Shortest distance from either end of peak to either end of feature

**Use Cases:**
- Standard peak annotation workflow
- Finding nearest genes to peaks
- Identifying overlapping features

---

#### `getAnnotation()`
**Purpose:** Obtains TSS, exon, miRNA, transcript, or UTR annotation for specified species using biomaRt.

**Key Parameters:**
- `mart`: biomaRt mart object
- `featureType`: "TSS", "miRNA", "Exon", "5utr", "3utr", "ExonPlusUtr", "transcript"

**Return Value:** GRanges object with annotation data including strand and description

**Use Cases:**
- Fetching up-to-date annotation from Ensembl
- Getting TSS or exon annotations for any supported organism

---

#### `annoGR()`
**Purpose:** Creates annotation GRanges objects from TxDb or EnsDb objects.

**Key Parameters:**
- `ranges`: TxDb, EnsDb, or GRanges object
- `feature`: Annotation type ("gene", "exon", "transcript", "CDS", "fiveUTR", "threeUTR", etc.)
- `source`: Character, annotation source
- `mdata`: Data frame with metadata

**Return Value:** annoGR object (extends GRanges) with annotation metadata

**Use Cases:**
- Converting TxDb/EnsDb to annotation format compatible with annotation functions
- Creating custom annotation objects

---

#### `peaksNearBDP()`
**Purpose:** Identifies peaks near bi-directional promoters with summary statistics.

**Key Parameters:**
- `myPeakList`: GRanges object
- `AnnotationData`: GRanges or annoGR object
- `MaxDistance`: Maximum gap allowed between peak and nearest gene

**Return Value:** List containing:
- `peaksWithBDP`: GRangesList of annotated peaks with bidirectional promoters
- `percentPeaksWithBDP`: Percentage of input peaks containing BDP
- `n.peaks`: Total number of input peaks
- `n.peaksWithBDP`: Number of peaks with BDP

**Use Cases:**
- Identifying enhancer regions
- Finding bidirectional promoter-associated peaks

---

#### `findEnhancers()`
**Purpose:** Finds possible enhancers using DNA interaction data (3C, 5C, HiC).

**Key Parameters:**
- `peaks`: GRanges object
- `annoData`: Annotation data
- `DNAinteractiveData`: DNA interaction data (GRanges with blocks, GInteractions, or BEDPE file)
- `bindingType`: "nearestBiDirectionalPromoters", "startSite", "endSite"
- `bindingRegion`: Vector of two integers

**Return Value:** GRanges object of annotated peaks with enhancer information

**Use Cases:**
- Identifying enhancers using chromosome conformation capture data
- Linking distal peaks to genes through DNA interactions

---

### 2. Overlap Analysis Functions

#### `findOverlapsOfPeaks()`
**Purpose:** Finds overlapping peaks among 2-5 sets of peak ranges (recommended function).

**Key Parameters:**
- `...`: 2-5 GRanges objects
- `maxgap`, `minoverlap`: Overlap detection parameters
- `ignore.strand`: Logical, ignore strand information
- `connectedPeaks`: "merge", "min", or "keepAll" - how to handle connected peaks

**Return Value:** overlappingPeaks object containing:
- `venn_cnt`: VennCounts object for Venn diagram
- `peaklist`: List of overlapping or unique peaks
- `uniquePeaks`: GRanges of all unique peaks
- `mergedPeaks`: GRanges of merged overlapping peaks
- `overlappingPeaks`: List of data frames with overlap annotations

**Use Cases:**
- Comparing multiple ChIP-seq experiments
- Finding common peaks across replicates
- Analyzing peak overlaps between different transcription factors

---

#### `findOverlappingPeaks()`
**Purpose:** Finds overlapping peaks for two peak ranges (deprecated, use `findOverlapsOfPeaks`).

**Key Parameters:**
- `Peaks1`, `Peaks2`: GRanges objects
- `maxgap`, `minoverlap`: Overlap parameters
- `select`: "all", "first", "last", "arbitrary"
- `annotate`: Include overlapFeature and shortestDistance (0 or 1)
- `connectedPeaks`: "min" or "merge"

**Return Value:** List with `OverlappingPeaks` data frame and `MergedPeaks` GRanges

**Use Cases:**
- Legacy code compatibility
- Simple two-way peak comparison

---

#### `makeVennDiagram()`
**Purpose:** Creates Venn diagrams from peak lists and calculates statistical significance of overlaps.

**Key Parameters:**
- `Peaks`: List of GRanges objects (2-5 sets)
- `NameOfPeaks`: Character vector of peak names
- `maxgap`, `minoverlap`: Overlap parameters
- `totalTest`: Total number of tests for hypergeometric test
- `by`: "region", "feature", or "base" - calculation method
- `method`: "hyperG" or "permutation" for p-value calculation
- `TxDb`: TxDb object (for permutation test)
- `plot`: Logical, whether to plot

**Return Value:** List with:
- `p.value`: Data frame with p-values for overlaps
- `vennCounts`: VennCounts object

**Use Cases:**
- Visualizing overlap between experimental replicates
- Statistical testing of peak overlap significance
- Comparing multiple ChIP-seq datasets

---

#### `getVennCounts()`
**Purpose:** Calculates Venn counts for Venn diagram generation (internal function for `makeVennDiagram`).

**Key Parameters:**
- `...`: GRanges objects
- `by`: "region", "feature", or "base"
- `connectedPeaks`: "min", "merge", or "keepAll"

**Return Value:** VennCounts object

**Use Cases:**
- Custom Venn diagram generation
- Overlap counting without visualization

---

### 3. Enrichment Analysis Functions

#### `getEnrichedGO()`
**Purpose:** Performs Gene Ontology enrichment analysis using hypergeometric tests.

**Key Parameters:**
- `annotatedPeak`: GRanges object or vector of feature IDs
- `orgAnn`: Organism annotation package (e.g., "org.Hs.eg.db")
- `feature_id_type`: "ensembl_gene_id", "refseq_id", "gene_symbol", or "entrez_id"
- `maxP`: Maximum p-value threshold
- `minGOterm`: Minimum count in genome for GO term inclusion
- `multiAdjMethod`: Multiple testing correction method (Bonferroni, BH, etc.)
- `condense`: Logical, condense results
- `removeAncestorByPval`: Remove ancestor terms by p-value
- `keepByLevel`: Filter by GO hierarchy level
- `subGroupComparison`: Logical vector for subgroup comparison

**Return Value:** List with three data frames:
- `bp`: Biological process enrichment
- `mf`: Molecular function enrichment
- `cc`: Cellular component enrichment

Each data frame contains: go.id, go.term, Definition, Ontology, count.InDataset, count.InGenome, pvalue, totaltermInDataset, totaltermInGenome

**Use Cases:**
- Functional enrichment analysis of ChIP-seq peaks
- Identifying biological processes associated with binding sites
- Comparing enrichment between experimental groups

---

#### `getEnrichedPATH()`
**Purpose:** Performs pathway enrichment analysis (KEGG/Reactome) using hypergeometric tests.

**Key Parameters:**
- `annotatedPeak`: GRanges object or vector of feature IDs
- `orgAnn`: Organism annotation package
- `pathAnn`: Pathway annotation ("reactome.db" or "KEGGREST")
- `feature_id_type`: Type of feature ID
- `maxP`: Maximum p-value threshold
- `minPATHterm`: Minimum count threshold
- `multiAdjMethod`: Multiple testing correction
- `subGroupComparison`: Logical vector for subgroup comparison

**Return Value:** Data frame with enriched pathways including:
- `path.id`: Pathway ID
- `path.term`: Pathway name
- `count.InDataset`: Count in dataset
- `count.InGenome`: Count in genome
- `pvalue`: P-value from hypergeometric test
- `EntrezID`: Associated Entrez IDs

**Use Cases:**
- Pathway enrichment analysis
- Identifying signaling pathways associated with peaks
- KEGG or Reactome pathway analysis

---

#### `hyperGtest()`
**Purpose:** Internal function for hypergeometric test used by enrichment functions.

**Key Parameters:**
- `alltermcount`: List with GOterm and GOcount for whole genome
- `thistermcount`: List with GOterm and GOcount for peak list
- `totaltermInGenome`: Total GO terms in genome
- `totaltermInPeakList`: Total GO terms in peak list

**Return Value:** List with test results including p-values

**Use Cases:**
- Internal use by enrichment functions
- Custom enrichment testing

---

#### `addAncestors()`
**Purpose:** Adds GO ancestor terms to enrichment results.

**Key Parameters:**
- `go.ids`: Matrix with GO IDs and Entrez IDs
- `ontology`: "bp", "cc", or "mf"

**Return Value:** Matrix with GO IDs including ancestors

**Use Cases:**
- Expanding GO enrichment results with ancestor terms
- Comprehensive GO term coverage

---

#### `enrichmentPlot()`
**Purpose:** Visualizes GO/pathway enrichment results.

**Key Parameters:**
- `res`: Output from `getEnrichedGO()` or `getEnrichedPATH()`
- `n`: Number of terms to plot (default 20)
- `style`: "v" (vertical) or "h" (horizontal)
- `strlength`: Maximum character length for term descriptions
- `orderBy`: "pvalue", "termId", or "none"
- `label_wrap`: Character wrap length for labels

**Return Value:** ggplot object

**Use Cases:**
- Publication-ready enrichment plots
- Visualizing top enriched GO terms or pathways

---

### 4. Sequence Analysis Functions

#### `getAllPeakSequence()`
**Purpose:** Extracts genomic sequences around peaks using BSgenome or biomaRt.

**Key Parameters:**
- `myPeakList`: GRanges object
- `upstream`: Upstream offset from peak start (default 200)
- `downstream`: Downstream offset from peak end (default 200)
- `genome`: BSgenome object or mart object
- `AnnotationData`: Optional GRanges object

**Return Value:** GRanges object with sequences in metadata, including:
- `upstream`: Upstream offset used
- `downstream`: Downstream offset used
- `sequence`: Extracted sequence

**Use Cases:**
- Extracting sequences for motif discovery
- Preparing sequences for PCR or cloning
- Sequence-based analysis

---

#### `summarizePatternInPeaks()`
**Purpose:** Summarizes motif/pattern occurrences and enrichment in peaks.

**Key Parameters:**
- `patternFilePath`: Path to FASTA file containing patterns
- `format`: "fasta" or "fastq"
- `BSgenomeName`: BSgenome object
- `peaks`: GRanges object
- `revcomp`: Logical, search reverse complement (default TRUE)
- `method`: "binom.test" or "permutation.test"
- `expectFrequencyMethod`: "Markov" or "Naive"
- `MarkovOrder`: Order of Markov chain (default 3)
- `bgdForPerm`: "chromosome" or "shuffle" for background
- `nperm`: Number of permutations (default 1000)
- `alpha`: Significance level (default 0.05)

**Return Value:** List with two data frames:
- `motif_enrichment`: Enrichment statistics (patternNum, totalNumPatternWithSameLen, expectedRate, patternRate, pValueBinomTest, cutOffPermutationTest)
- `motif_occurrence`: Detailed occurrence information (motifChr, motifStartInChr, motifEndInChr, motifName, motifPattern, motifStartInPeak, motifEndInPeak, motifFound, motifFoundStrand, peakChr, peakStart, peakEnd, peakWidth, peakStrand)

**Use Cases:**
- Known motif scanning in peaks
- Motif enrichment analysis
- Pattern occurrence statistics

---

#### `countPatternInSeqs()`
**Purpose:** Counts total number of pattern occurrences in sequences.

**Key Parameters:**
- `pattern`: DNAStringSet object
- `sequences`: Vector of sequences

**Return Value:** Numeric, total count of pattern occurrences

**Use Cases:**
- Quick pattern counting
- Simple motif presence/absence testing

---

#### `findMotifsInPromoterSeqs()`
**Purpose:** Finds motif occurrences in promoter regions of input gene list.

**Key Parameters:**
- `patternFilePath1`: File path with known motifs (required)
- `patternFilePath2`: File path with second motif (for paired motifs)
- `BSgenomeName`: BSgenome object
- `txdb`: TxDb object
- `geneIDs`: Entrez gene IDs
- `upstream`, `downstream`: Promoter region (default 5000 each)
- `findPairedMotif`: Logical, find paired motifs (default FALSE)
- `max.distance`, `min.distance`: Distance thresholds for paired motifs
- `motif.orientation`: "both", "motif1UpstreamOfMotif2", "motif2UpstreamOfMoif1"
- `outfile`: Output file path

**Return Value:** GRanges object with motif occurrence information

**Use Cases:**
- Scanning promoter regions for known motifs
- Finding paired motif configurations
- Motif co-occurrence analysis

---

#### `translatePattern()`
**Purpose:** Translates IUPAC nucleotide codes to regex patterns.

**Key Parameters:**
- Pattern string with IUPAC codes

**Return Value:** Translated regex pattern

**Use Cases:**
- Converting IUPAC codes for pattern matching
- Supporting degenerate nucleotide searches

---

#### `write2FASTA()`
**Purpose:** Writes sequences to FASTA format file.

**Key Parameters:**
- `mySeq`: GRanges object with sequence metadata
- `file`: Output file path
- `width`: Maximum letters per line (default 80)

**Return Value:** Writes FASTA file

**Use Cases:**
- Exporting sequences for external tools
- Preparing sequences for motif discovery software

---

### 5. Distribution Analysis Functions

#### `assignChromosomeRegion()`
**Purpose:** Assigns peaks to genomic regions (exon, intron, UTR, promoter, etc.).

**Key Parameters:**
- `peaks.RD`: GRanges object
- `TxDb`: TxDb object (recommended) or individual annotation GRanges
- `proximal.promoter.cutoff`: Upstream/downstream for promoter (default c(2000, 100))
- `immediate.downstream.cutoff`: Upstream/downstream for downstream region (default c(0, 1000))
- `nucleotideLevel`: Logical, nucleotide-centric vs peak-centric view
- `precedence`: Order of precedence for overlapping regions

**Return Value:** List with:
- `percentage`: Named vector with percentage distribution
- `jaccard`: Named vector with Jaccard indices

Regions include: Exons, Introns, fiveUTRs, threeUTRs, Promoter, ImmediateDownstream, Intergenic.Region

**Use Cases:**
- Genomic distribution analysis
- Determining where peaks are located relative to genes
- Publication-ready distribution summaries

---

#### `genomicElementDistribution()`
**Purpose:** Calculates and plots genomic element distribution with customizable categories.

**Key Parameters:**
- `peaks`: GRanges or GRangesList
- `TxDb`: TxDb object
- `nucleotideLevel`: Logical, nucleotide-centric view
- `promoterRegion`: Upstream/downstream for promoter
- `geneDownstream`: Upstream/downstream for downstream region
- `labels`: List of labels for genomic elements
- `labelColors`: Named vector of colors
- `promoterLevel`: List with breaks, labels, colors for promoter sub-regions
- `plot`: Logical, whether to plot

**Return Value:** Invisible list with peaks (annotated) and plot (ggplot object)

**Use Cases:**
- Detailed genomic distribution analysis
- Customizable region definitions
- Publication-quality pie charts and bar plots

---

#### `genomicElementUpSetR()`
**Purpose:** Prepares data for UpSet plots of genomic element distribution.

**Key Parameters:**
- `peaks`: GRanges or GRangesList
- `TxDb`: TxDb object
- `breaks`: List defining genomic element sets (functions or numeric vectors)

**Return Value:** List with:
- `peaks`: Annotated peaks
- `plotData`: Data frame for UpSetR plotting

**Use Cases:**
- Creating UpSet plots for genomic element overlaps
- Complex multi-set visualization

---

#### `binOverFeature()`
**Purpose:** Aggregates peaks over bins from feature sites (TSS, gene end, etc.).

**Key Parameters:**
- `...`: GRanges objects to analyze
- `annotationData`: GRanges or annoGR for annotation
- `select`: "all" or "nearest"
- `radius`: Longest distance to feature site
- `nbins`: Number of bins
- `aroundGene`: Logical, count peaks around features or at site
- `mbins`: Number of bins intra-feature (if aroundGene=TRUE)
- `featureSite`: "FeatureStart", "FeatureEnd", "bothEnd"
- `PeakLocForDistance`: "all", "end", "start", "middle"
- `FUN`: Function for score calculation (sum, mean, median, etc.)
- `errFun`: Function for error bar calculation

**Return Value:** Data frame with bin values

**Use Cases:**
- Binding profile analysis around TSS
- Metagene plots
- Distance-dependent binding analysis

---

#### `binOverGene()`
**Purpose:** Calculates coverage of gene body per gene per bin.

**Key Parameters:**
- `cvglists`: List of RleList objects (coverage data)
- `TxDb`: TxDb object
- `upstream.cutoff`, `downstream.cutoff`: Cutoff lengths
- `nbinsGene`, `nbinsUpstream`, `nbinsDownstream`: Number of bins
- `includeIntron`: Logical, include introns
- `minGeneLen`, `maxGeneLen`: Gene length filters

**Return Value:** List of matrices with coverage data

**Use Cases:**
- Gene body coverage analysis
- Metagene analysis from BigWig files
- Average gene profiles

---

#### `binOverRegions()`
**Purpose:** Calculates coverage of 5'UTR, CDS, and 3'UTR per transcript per bin.

**Key Parameters:**
- `cvglists`: List of RleList objects
- `TxDb`: TxDb object
- `nbinsCDS`, `nbinsUTR`, `nbinsUpstream`, `nbinsDownstream`: Number of bins
- `includeIntron`: Logical, include introns
- `minCDSLen`, `minUTRLen`, `maxCDSLen`, `maxUTRLen`: Length filters

**Return Value:** List of matrices with coverage by region

**Use Cases:**
- Transcript-level coverage analysis
- UTR and CDS-specific analysis
- Ribosome profiling data analysis

---

#### `plotBinOverRegions()`
**Purpose:** Plots coverage of regions from `binOverRegions()` or `binOverGene()` output.

**Key Parameters:**
- `dat`: List of matrices from binOver functions
- `...`: Additional parameters for matplot

**Return Value:** Plots coverage profiles

**Use Cases:**
- Visualizing metagene profiles
- Coverage plots across transcript regions

---

### 6. Signal Analysis Functions

#### `featureAlignedSignal()`
**Purpose:** Extracts signals in given feature ranges.

**Key Parameters:**
- `cvglists`: List of SimpleRleList or RleList objects
- `feature.gr`: GRanges object with identical width
- `upstream`, `downstream`: Upstream/downstream from feature center
- `n.tile`: Number of tiles per feature (default 100)

**Return Value:** List of matrices, each row represents a feature, columns represent tiles

**Use Cases:**
- Extracting ChIP-seq signals around features
- Preparing data for heatmaps
- Signal aggregation around TSS or other features

---

#### `featureAlignedDistribution()`
**Purpose:** Plots distribution of signals in feature ranges.

**Key Parameters:**
- `cvglists`: Output of `featureAlignedSignal()` or list of RleList
- `feature.gr`: GRanges object
- `upstream`, `downstream`: Upstream/downstream regions
- `zeroAt`: Zero point position
- `n.tile`: Number of tiles
- `...`: Additional parameters for matplot

**Return Value:** Invisible matrix of plot data

**Use Cases:**
- Average signal profiles
- Metagene plots
- Signal distribution visualization

---

#### `featureAlignedHeatmap()`
**Purpose:** Creates heatmaps representing signals in given ranges.

**Key Parameters:**
- `cvglists`: Output of `featureAlignedSignal()` or list of RleList
- `feature.gr`: GRanges object
- `upstream`, `downstream`: Upstream/downstream regions
- `zeroAt`: Zero point position
- `n.tile`: Number of tiles
- `annoMcols`: Metadata columns for annotation
- `sortBy`: Sort features by columns or signals
- `color`: Color palette
- `lower.extreme`, `upper.extreme`: Color scale boundaries

**Return Value:** Invisible gList object (grid graphics)

**Use Cases:**
- Publication-quality heatmaps
- Signal visualization around features
- Multi-sample comparison heatmaps

---

#### `featureAlignedExtendSignal()`
**Purpose:** Extends signals around features (similar to featureAlignedSignal with extension).

**Use Cases:**
- Extended signal extraction
- Long-range signal analysis

---

#### `metagenePlot()`
**Purpose:** Creates bar plot for distance to features (metagene plot).

**Key Parameters:**
- `peaks`: GRanges or GRangesList
- `AnnotationData`: TxDb or GRanges
- `PeakLocForDistance`: "middle", "start", "end"
- `FeatureLocForDistance`: "TSS", "middle", "geneEnd"
- `upstream`, `downstream`: Region to plot (default 100000 each)

**Return Value:** ggplot object

**Use Cases:**
- Metagene visualization
- Distance distribution analysis
- Binding profile around features

---

### 7. Statistical Testing Functions

#### `peakPermTest()`
**Purpose:** Performs permutation tests to assess association between two peak lists.

**Key Parameters:**
- `peaks1`, `peaks2`: GRanges objects
- `ntimes`: Number of permutations (default 100)
- `seed`: Random seed
- `mc.cores`: Number of cores for parallel processing
- `maxgap`: Maximum gap for overlap
- `pool`: permPool object (or will be created)
- `TxDb`: TxDb object
- `bindingDistribution`: bindist object
- `bindingType`: "TSS" or "geneEnd"
- `featureType`: "transcript" or "exon"

**Return Value:** permTestResults object from regioneR package

**Use Cases:**
- Statistical significance of peak overlaps
- Comparing experimental conditions
- Validating peak associations

---

#### `preparePool()`
**Purpose:** Prepares permutation pool for statistical testing.

**Key Parameters:**
- `TxDb`: TxDb object
- `template`: GRanges object (template peaks)
- `bindingDistribution`: bindist object
- `bindingType`: "TSS" or "geneEnd"
- `featureType`: "transcript" or "exon"
- `seqn`: Sequence names to filter

**Return Value:** List with:
- `grs`: GRangesList for sampling
- `N`: Numbers to draw from each GRanges

**Use Cases:**
- Preparing background for permutation tests
- Custom statistical testing

---

#### `estFragmentLength()`
**Purpose:** Estimates fragment length from BAM files.

**Key Parameters:**
- `bamfiles`: BAM file paths
- `index`: Index file paths
- `plot`: Logical, plot ACF (default TRUE)
- `lag.max`: Maximum lag for ACF (default 1000)
- `minFragmentSize`: Minimal fragment size to avoid phantom peak (default 100)

**Return Value:** Numeric vector of estimated fragment lengths

**Use Cases:**
- ChIP-seq quality control
- Fragment size estimation for paired-end data
- Cross-correlation analysis for single-end data

---

#### `estLibSize()`
**Purpose:** Estimates library size from BAM files.

**Key Parameters:**
- `bamfiles`: BAM file paths
- `index`: Index file paths

**Return Value:** Numeric vector of library sizes

**Use Cases:**
- Library size estimation
- Normalization calculations
- Quality control

---

### 8. Data Conversion & Utility Functions

#### `toGRanges()`
**Purpose:** Converts various formats (BED, GFF, data.frame, etc.) to GRanges.

**Key Parameters:**
- `data`: Data frame, file path, or other supported format
- `format`: "BED", "GFF", "MACS", "MACS2", "narrowPeak", "broadPeak", "RangedData", etc.
- `colNames`: Column names (if data frame)
- `header`: Logical, file has header
- `skip`: Number of lines to skip

**Return Value:** GRanges object

**Use Cases:**
- Converting peak files to GRanges
- Format conversion for downstream analysis
- Importing external peak files

---

#### `addGeneIDs()`
**Purpose:** Adds gene identifiers (symbol, Entrez ID, etc.) to annotated peaks.

**Key Parameters:**
- `annotatedPeak`: GRanges object or vector of feature IDs
- `orgAnn`: Organism annotation package (e.g., "org.Hs.eg.db")
- `IDs2Add`: Vector of IDs to add ("symbol", "entrez_id", "ensembl", "refseq", etc.)
- `feature_id_type`: Type of input ID
- `silence`: Logical, suppress unmapped ID messages
- `mart`: biomaRt mart object (alternative to orgAnn)

**Return Value:** GRanges (if input is GRanges) or data frame (if input is vector)

**Use Cases:**
- Adding gene symbols to annotated peaks
- Converting between ID types
- Enriching annotation with additional identifiers

---

#### `convert2EntrezID()`
**Purpose:** Converts gene IDs to Entrez IDs.

**Key Parameters:**
- `IDs`: Vector of IDs
- `orgAnn`: Organism annotation package
- `ID_type`: "ensembl_gene_id", "gene_symbol", or "refseq_id"

**Return Value:** Vector of Entrez IDs

**Use Cases:**
- ID conversion for enrichment analysis
- Preparing data for GO analysis

---

#### `reCenterPeaks()`
**Purpose:** Re-centers peaks based on peak centers.

**Key Parameters:**
- `peaks`: GRanges or annoGR object
- `width`: Width of new peaks (default 2000)

**Return Value:** GRanges object with re-centered peaks

**Use Cases:**
- Standardizing peak widths
- Creating fixed-width peaks for analysis

---

#### `mergePlusMinusPeaks()`
**Purpose:** Merges peaks from plus and minus strands within distance threshold.

**Key Parameters:**
- `peaks.file`: Peak file path
- `columns`: Column names
- `distance.threshold`: Maximum gap between plus and minus peaks
- `plus.strand.start.gt.minus.strand.end`: Logical, strand orientation
- `output.bedfile`: Output BED file path

**Return Value:** Data frame in BED format

**Use Cases:**
- Processing GUIDE-seq or similar data
- Merging strand-specific peaks

---

#### `tileGRanges()`
**Purpose:** Tiles GRanges into bins.

**Key Parameters:**
- GRanges object and tiling parameters

**Return Value:** GRangesList of tiled regions

**Use Cases:**
- Creating bins for coverage analysis
- Tiling genome regions

---

#### `tileCount()`
**Purpose:** Counts overlaps in tiled regions.

**Key Parameters:**
- Tiled regions and query GRanges

**Return Value:** Count matrix

**Use Cases:**
- Coverage counting in bins
- Tiled analysis

---

### 9. Visualization Functions

#### `pie1()`
**Purpose:** Draws pie charts with percentage labels.

**Key Parameters:**
- `x`: Vector of numerical quantities
- `labels`: Labels for slices
- `percentage`: Logical, add percentage (default TRUE)
- `rawNumber`: Logical, add raw numbers instead
- `digits`: Significant digits for percentage
- `cutoff`: Minimum percentage to show
- `legend`: Logical, use legend instead of labels

**Return Value:** Pie chart plot

**Use Cases:**
- Genomic distribution visualization
- Category proportion plots

---

#### `cumulativePercentage()`
**Purpose:** Plots cumulative percentage tag allocation in samples.

**Key Parameters:**
- `bamfiles`: BAM file paths
- `gr`: GRanges object (target regions)
- `input`: Which file is input (default 1)
- `binWidth`: Width of each bin (default 1000)

**Return Value:** List of data frames with cumulative percentages

**Use Cases:**
- Normalization assessment
- Bias detection
- Quality control plots

---

### 10. Additional Utility Functions

#### `condenseMatrixByColnames()`
**Purpose:** Condenses matrices by column names.

**Use Cases:**
- Data aggregation
- Matrix manipulation

---

#### `getUniqueGOidCount()`
**Purpose:** Gets unique GO ID counts.

**Use Cases:**
- GO term counting
- Enrichment analysis support

---

#### `getGO()`
**Purpose:** Obtains GO terms for given genes.

**Key Parameters:**
- `all.genes`: Character vector of feature IDs
- `orgAnn`: Organism annotation package
- `ID_type`: Feature ID type
- `writeTo`: Optional output file path

**Return Value:** Invisible table with genes and GO terms

**Use Cases:**
- Retrieving GO annotations for genes
- GO term lookup

---

#### `getGeneSeq()`
**Purpose:** Gets gene sequences.

**Use Cases:**
- Sequence retrieval
- Gene sequence extraction

---

#### `downstreams()`
**Purpose:** Gets downstream regions from features.

**Use Cases:**
- Downstream region definition
- Feature extension

---

#### `egOrgMap()`
**Purpose:** Maps organism names to Entrez Gene database names.

**Use Cases:**
- Organism name conversion
- Database mapping

---

#### `xget()`
**Purpose:** Extended get function for annotation databases.

**Use Cases:**
- Annotation data retrieval
- Database access

---

#### `IDRfilter()`
**Purpose:** Filters peaks using Irreproducible Discovery Rate (IDR).

**Key Parameters:**
- `peaksA`, `peaksB`: GRanges objects
- `bamfileA`, `bamfileB`: BAM file paths
- `maxgap`, `minoverlap`: Overlap parameters
- `singleEnd`: Logical, single-end reads (default TRUE)
- `IDRcutoff`: IDR threshold (default 0.01)

**Return Value:** GRanges object with filtered peaks

**Use Cases:**
- Replicate consistency assessment
- High-confidence peak identification
- IDR-based peak filtering

---

#### `addMetadata()`
**Purpose:** Adds metadata to objects.

**Use Cases:**
- Metadata management
- Object annotation

---

#### `translatePattern()`
**Purpose:** Translates IUPAC nucleotide codes to regex patterns.

**Use Cases:**
- Pattern matching support
- Degenerate sequence searches

---

#### `oligoFrequency()`
**Purpose:** Calculates oligonucleotide frequency.

**Use Cases:**
- Sequence composition analysis
- K-mer counting

---

#### `oligoSummary()`
**Purpose:** Summarizes oligonucleotide patterns.

**Use Cases:**
- Pattern summary statistics
- Sequence analysis

---

#### `summarizeOverlapsByBins()`
**Purpose:** Summarizes overlaps by bins.

**Use Cases:**
- Binned overlap analysis
- Coverage summarization

---

## Common Workflows

### Workflow 1: Basic Peak Annotation
```
1. Load peaks: peaks <- toGRanges("peaks.bed", format="BED")
2. Get annotation: anno <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene)
3. Annotate: annotated <- annotatePeakInBatch(peaks, AnnotationData=anno)
4. Add gene symbols: annotated <- addGeneIDs(annotated, orgAnn="org.Hs.eg.db", IDs2Add="symbol")
5. Distribution: dist <- assignChromosomeRegion(peaks, TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
```

### Workflow 2: Enrichment Analysis
```
1. Annotate peaks: annotated <- annotatePeakInBatch(peaks, AnnotationData=anno)
2. GO enrichment: enrichedGO <- getEnrichedGO(annotated, orgAnn="org.Hs.eg.db", maxP=0.01)
3. Pathway enrichment: enrichedPATH <- getEnrichedPATH(annotated, orgAnn="org.Hs.eg.db", pathAnn="reactome.db")
4. Visualize: enrichmentPlot(enrichedGO)
```

### Workflow 3: Overlap Analysis
```
1. Load multiple peak sets: peaks1 <- toGRanges("set1.bed"); peaks2 <- toGRanges("set2.bed")
2. Find overlaps: overlaps <- findOverlapsOfPeaks(peaks1, peaks2)
3. Venn diagram: makeVennDiagram(overlaps)
4. Statistical test: permTest <- peakPermTest(peaks1, peaks2, TxDb=txdb)
```

### Workflow 4: Signal Analysis
```
1. Load coverage: cvglists <- list(sample1=RleList(...), sample2=RleList(...))
2. Get features: features <- promoters(genes(TxDb), upstream=2000, downstream=2000)
3. Extract signals: signals <- featureAlignedSignal(cvglists, features)
4. Plot heatmap: featureAlignedHeatmap(signals, features)
5. Plot distribution: featureAlignedDistribution(signals, features)
```

### Workflow 5: Motif Analysis
```
1. Extract sequences: seqs <- getAllPeakSequence(peaks, genome=BSgenome, upstream=200, downstream=200)
2. Summarize patterns: motifs <- summarizePatternInPeaks("motifs.fa", peaks=peaks, BSgenomeName=BSgenome)
3. Find in promoters: promoterMotifs <- findMotifsInPromoterSeqs("motif.fa", BSgenomeName=BSgenome, txdb=TxDb, geneIDs=geneIDs)
```

---

## Key Use Cases

1. **Peak Annotation**: Identify nearest genes, genomic regions, and distances
2. **Functional Enrichment**: GO and pathway enrichment analysis
3. **Overlap Analysis**: Compare multiple ChIP-seq experiments
4. **Distribution Analysis**: Determine where peaks are located relative to genes
5. **Signal Analysis**: Extract and visualize ChIP-seq signals around features
6. **Motif Analysis**: Scan for known motifs in peak sequences
7. **Statistical Testing**: Assess significance of overlaps and enrichments
8. **Quality Control**: Estimate fragment lengths, library sizes, and assess reproducibility

---

## Integration with Other Packages

- **GenomicRanges/IRanges**: Core data structures
- **TxDb/EnsDb**: Annotation databases
- **BSgenome**: Genome sequences
- **biomaRt**: Online annotation queries
- **GO.db/org.*.eg.db**: GO and gene annotations
- **regioneR**: Permutation testing
- **VennDiagram**: Venn diagram plotting
- **ggplot2**: Advanced plotting
- **clusterProfiler**: Functional enrichment (via ChIPseeker integration)

---

## Summary Statistics

- **Total Exported Functions**: 58
- **Main Categories**: 10
- **Peak Annotation Functions**: 6
- **Overlap Analysis Functions**: 4
- **Enrichment Analysis Functions**: 5
- **Sequence Analysis Functions**: 6
- **Distribution Analysis Functions**: 7
- **Signal Analysis Functions**: 5
- **Statistical Testing Functions**: 4
- **Data Conversion Functions**: 8
- **Visualization Functions**: 3
- **Utility Functions**: 10

---

## References

Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 11:237. doi:10.1186/1471-2105-11-237

---

*Document generated from ChIPpeakAnno package version 3.45.2*

