# Parallel Computing Framework Comparison: BiocParallel vs future.apply

## Executive Summary

**For ChIPpeakAnno (a Bioconductor package): BiocParallel is the recommended choice.**

## Detailed Comparison

### 1. **BiocParallel** (Recommended for Bioconductor)

#### Advantages:
- ✅ **Bioconductor Standard**: Official parallel computing framework for Bioconductor packages
- ✅ **Ecosystem Integration**: Seamlessly integrates with other Bioconductor packages
- ✅ **Automatic Backend Detection**: Automatically detects and uses appropriate backend
- ✅ **Progress Reporting**: Built-in `progressbar` support via `progressr` package
- ✅ **Error Handling**: Robust error handling and logging capabilities
- ✅ **Resource Management**: Better memory and resource management for large genomic datasets
- ✅ **Consistent API**: `bplapply()`, `bpiterate()`, `bpaggregate()` - consistent interface
- ✅ **Cluster Support**: Excellent support for HPC clusters (SGE, SLURM, LSF, etc.)
- ✅ **BiocCheck Compliance**: Aligns with Bioconductor best practices

#### Disadvantages:
- ⚠️ **Bioconductor Dependency**: Requires Bioconductor installation
- ⚠️ **Slightly More Verbose**: Requires `BiocParallelParam` object setup

#### Code Example:
```r
library(BiocParallel)

# Automatic backend detection
bp <- bpparam()  # Uses default (usually MulticoreParam or SnowParam)

# Or explicit setup
bp <- MulticoreParam(workers = 4, progressbar = TRUE)

# Use in function
results <- bplapply(chr_names, process_chromosome, BPPARAM = bp)
```

---

### 2. **future.apply** (General-Purpose)

#### Advantages:
- ✅ **Flexible Backends**: Supports many backends (multicore, multisession, cluster, remote)
- ✅ **Easy Migration**: Drop-in replacement for base R `lapply()` functions
- ✅ **Future Framework**: Part of the broader `future` ecosystem
- ✅ **Cross-Platform**: Works well on Windows (multisession)
- ✅ **No Bioconductor Dependency**: Can be used in CRAN packages

#### Disadvantages:
- ⚠️ **Manual Setup Required**: User must set `plan()` before use
- ⚠️ **Less Bioconductor Integration**: Not the standard for Bioconductor packages
- ⚠️ **Progress Reporting**: Requires additional setup (e.g., `progressr`)
- ⚠️ **Error Handling**: Less sophisticated than BiocParallel

#### Code Example:
```r
library(future.apply)

# User must set plan
plan(multicore, workers = 4)  # or plan(multisession, workers = 4)

# Use in function
results <- future_lapply(chr_names, process_chromosome)
```

---

### 3. **Hybrid Approach: BiocParallel.FutureParam**

You can combine both frameworks:

```r
library(BiocParallel)
library(future)
library(BiocParallel.FutureParam)

# Use future backend with BiocParallel API
plan(multicore, workers = 4)
bp <- FutureParam()

results <- bplapply(chr_names, process_chromosome, BPPARAM = bp)
```

**Benefits**: BiocParallel API + future backend flexibility

---

## Recommendation for ChIPpeakAnno

### **Use BiocParallel** for the following reasons:

1. **Bioconductor Standard**: ChIPpeakAnno is a Bioconductor package, so using BiocParallel aligns with ecosystem standards

2. **Better User Experience**: 
   - Automatic backend detection (`bpparam()`)
   - No need for users to manually set plans
   - Consistent with other Bioconductor packages

3. **Better Error Handling**: Important for genomic data processing where errors can be costly

4. **Progress Reporting**: Built-in support for progress bars (via `progressr`)

5. **HPC Support**: Better support for cluster environments (common in bioinformatics)

6. **BiocCheck Compliance**: Helps pass Bioconductor package checks

### Implementation Strategy:

```r
# In DESCRIPTION, add to Imports or Suggests:
Suggests: BiocParallel

# In function:
metageneSignalProfile <- function(..., BPPARAM = NULL) {
    # Use BiocParallel if available, fallback to sequential
    if (is.null(BPPARAM)) {
        if (requireNamespace("BiocParallel", quietly = TRUE)) {
            BPPARAM <- BiocParallel::bpparam()
        } else {
            BPPARAM <- BiocParallel::SerialParam()
        }
    }
    
    # Process chromosomes
    chr_results <- BiocParallel::bplapply(
        chr_names, 
        process_chromosome, 
        BPPARAM = BPPARAM
    )
    
    # ... rest of function
}
```

### Migration Path:

1. **Phase 1**: Add BiocParallel support alongside existing `mc.cores` parameter
2. **Phase 2**: Deprecate `mc.cores` in favor of `BPPARAM`
3. **Phase 3**: Remove `mc.cores` support

---

## Performance Comparison

Both frameworks have similar performance characteristics:
- **Multicore**: Similar speedup (3-4x on 4 cores)
- **Memory**: BiocParallel has slightly better memory management
- **Overhead**: Minimal difference (< 5%)

---

## Current State in ChIPpeakAnno

- ✅ `peakPermTest`: Uses `parallel::mclapply` with `mc.cores`
- ✅ `multiStrategyAnnotation`: Uses `future.apply::future_lapply`
- ✅ `metageneSignalProfile`: Uses `parallel::mclapply` with `mc.cores`

**Recommendation**: Standardize on BiocParallel for consistency and better Bioconductor integration.

---

## References

- [BiocParallel Documentation](https://bioconductor.org/packages/BiocParallel/)
- [future.apply Documentation](https://future.apply.futureverse.org/)
- [BiocParallel.FutureParam](https://biocparallel.futureparam.futureverse.org/)
- [Bioconductor Parallel Computing Guide](https://bioconductor.org/help/course-materials/2017/BioC2017/Day2/Workshops/Parallel/Parallel.html)

