#' Generate metagene signal profile with scaled gene bodies
#' 
#' @description 
#' Generates signal profiles along metagenes (TSS to TES with extensions) from
#' either BigWig files (signal intensity) or called peaks (peak coverage
#' density) and annotation databases. All gene bodies are scaled to a fixed
#' length (default 5kb) to allow direct comparison across genes of different
#' lengths, while upstream and downstream regions remain at their original bp
#' scale.
#' 
#' This function is useful for visualizing:
#' \itemize{
#'   \item Signal intensity patterns along genes (from BigWig files)
#'   \item Peak coverage/density patterns along genes (from called peaks)
#'   \item Binding patterns across genes of different lengths (via scaling)
#' }
#' 
#' BigWig files are assumed to contain pre-normalized signal (e.g., RPKM, CPM).
#' 
#' @param bigwigFiles Character vector of BigWig file paths. Optional if
#'        \code{peaks} is provided. BigWig files should contain normalized
#'        signal (e.g., RPKM, CPM). Files will be read using
#'        \code{rtracklayer::import()}.
#' @param peaks A \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
#'        \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object
#'        containing called peaks. Optional if \code{bigwigFiles} is provided.
#'        If a single \code{GRanges} is provided, it will be converted to a
#'        \code{GRangesList} with one element. Peak coverage/density will be
#'        calculated along metagene regions.
#' @param AnnotationData An object of class
#'        \code{\link[GenomicFeatures:TxDb-class]{TxDb}} or
#'        \code{\link[ensembldb:EnsDb-class]{EnsDb}} containing gene
#'        annotations. Required.
#' @param upstream An integer specifying the number of base pairs upstream
#'        from TSS to include in the metagene profile. Default is 2000L.
#' @param downstream An integer specifying the number of base pairs downstream
#'        from TES to include in the metagene profile. Default is 2000L.
#' @param geneBodyLength An integer specifying the fixed length for scaled
#'        gene bodies in base pairs. Default is 5000L (5kb). All gene bodies
#'        will be scaled to this length to allow direct comparison across genes
#'        of different lengths.
#' @param upstreamBins An integer specifying the number of bins for the
#'        upstream region. Default is 20L.
#' @param geneBodyBins An integer specifying the number of bins for the
#'        scaled gene body. Default is 50L.
#' @param downstreamBins An integer specifying the number of bins for the
#'        downstream region. Default is 20L.
#' @param minGeneLength An integer specifying the minimum gene length (in bp)
#'        to include. Genes shorter than this will be filtered out. Default is
#'        \code{NULL} (no filtering).
#' @param peakWeight Optional. For peak-based input, specifies weights for
#'        peaks. Can be:
#'        \itemize{
#'          \item A character string: name of a metadata column in \code{peaks}
#'                containing weights (e.g., "score", "signalValue")
#'          \item A numeric vector: weights for each peak (must match length of
#'                \code{peaks})
#'          \item \code{NULL} (default): uniform weights (count = 1 per peak)
#'        }
#' @param saveData Logical. If \code{TRUE} (default), saves signal/coverage
#'        data matrices to RDS files instead of returning them in the result.
#'        This reduces memory usage for large datasets. Files are saved as
#'        \code{metageneSignalProfile_data_[timestamp].rds} in the current
#'        working directory, or to the path specified in \code{dataFile}.
#' @param dataFile Character string. Path to save the RDS file containing data
#'        matrices. If \code{NULL} (default), a timestamped filename is
#'        generated automatically. Only used if \code{saveData = TRUE}.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'        object specifying the parallel backend to use. If \code{NULL} (default),
#'        uses \code{\link[BiocParallel]{bpparam}()} which automatically detects
#'        the best available backend. For sequential processing, use
#'        \code{BiocParallel::SerialParam()}. For parallel processing, use
#'        \code{BiocParallel::MulticoreParam(workers = n)} or
#'        \code{BiocParallel::SnowParam(workers = n)}. Parallel processing
#'        processes chromosomes in parallel, which can significantly speed up
#'        analysis for large genomes with many chromosomes.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages
#'        during processing. Set to \code{FALSE} to suppress messages.
#' @param ... Additional arguments (currently not used)
#' 
#' @return Returns a list with the following components:
#' \itemize{
#'   \item \code{plot}: A ggplot object showing the metagene profile
#'   \item \code{average}: A data.frame with average signal/coverage per bin
#'         across all genes (one column per input file/peak set)
#'   \item \code{regions}: A \code{GRangesList} with metagene regions used
#'         (scaled gene bodies)
#'   \item \code{binInfo}: A data.frame with bin positions and labels
#'   \item \code{scalingFactors}: A data.frame with original and scaled gene
#'         lengths for reference
#'   \item \code{metadata}: A list with input type, file/peak set names, data
#'         file path (if saved), and other metadata
#' }
#' 
#' \strong{Data Storage:}
#' By default (\code{saveData = TRUE}), signal/coverage data matrices are saved
#' to an RDS file instead of being returned in the result. This reduces memory
#' usage. The file path is stored in \code{metadata$dataFile}. To load the data:
#' \code{data <- readRDS(metadata$dataFile)}
#' 
#' @details
#' 
#' \strong{Input Requirements:}
#' \itemize{
#'   \item Exactly one of \code{bigwigFiles} or \code{peaks} must be provided
#'   \item BigWig files must exist and be readable
#'   \item Peaks must be a valid GRanges or GRangesList object
#'   \item AnnotationData must be a TxDb or EnsDb object
#' }
#' 
#' \strong{Gene Body Scaling:}
#' All gene bodies are scaled to a fixed length (\code{geneBodyLength}) to
#' allow direct comparison across genes of different lengths:
#' \itemize{
#'   \item Original gene body: TSS to TES (variable length)
#'   \item Scaled gene body: TSS to TSS + \code{geneBodyLength} - 1 (fixed
#'         length)
#'   \item Signal/coverage from original gene body is mapped to scaled
#'         coordinates using linear interpolation
#'   \item Upstream and downstream regions remain at original bp scale
#' }
#' 
#' \strong{Signal Extraction (BigWig):}
#' \itemize{
#'   \item BigWig files are read using \code{rtracklayer::import()}
#'   \item Signal is extracted at original genomic coordinates
#'   \item Gene body signal is interpolated to scaled coordinates
#'   \item Signal values are used directly (assumed normalized in BigWig)
#' }
#' 
#' \strong{Peak Coverage (Peaks):}
#' \itemize{
#'   \item Peak coverage is calculated by counting overlaps per bin
#'   \item If \code{peakWeight} is provided, weighted coverage is calculated
#'   \item Coverage is normalized by bin size to get density
#'   \item Gene body coverage is interpolated to scaled coordinates
#' }
#' 
#' \strong{Binning:}
#' The metagene profile is divided into three regions:
#' \itemize{
#'   \item Upstream: \code{upstreamBins} bins covering upstream region
#'   \item Gene body: \code{geneBodyBins} bins covering scaled gene body
#'   \item Downstream: \code{downstreamBins} bins covering downstream region
#' }
#' 
#' @author Jianhong Ou (concept), Implementation by ChIPpeakAnno contributors
#' @seealso \code{\link{metagenePlot}} for peak distance visualization,
#'          \code{\link{featureAlignedSignal}} for signal extraction from
#'          RleList objects, \code{\link{featureAlignedDistribution}} for
#'          average signal plots, \code{\link[BiocParallel]{BiocParallelParam}}
#'          for parallel backend options
#' @importFrom BiocParallel bpparam bplapply SerialParam MulticoreParam bpnworkers
#' @keywords misc
#' @export
#' @import GenomicRanges
#' @importFrom BiocGenerics start end width strand
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle seqlengths
#' @importFrom GenomicFeatures genes promoters
#' @importFrom S4Vectors Rle runLength runValue
#' @importFrom rtracklayer import
#' @importFrom ggplot2 ggplot aes_string geom_line geom_vline theme_bw labs
#' @importFrom reshape2 melt
#' @importFrom stats approx
#' @importFrom reshape2 melt
#' @examples
#' \dontrun{
#' ## Example 1: Using BigWig files with TxDb
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' bigwigFiles <- c("sample1.bw", "sample2.bw")
#' result <- metageneSignalProfile(
#'     bigwigFiles = bigwigFiles,
#'     AnnotationData = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'     upstream = 2000L,
#'     downstream = 2000L
#' )
#' result$plot
#' 
#' ## Example 2: Using peaks with EnsDb
#' library(EnsDb.Hsapiens.v86)
#' peaks <- GRanges("chr1", IRanges(1000000, 2000000))
#' result <- metageneSignalProfile(
#'     peaks = peaks,
#'     AnnotationData = EnsDb.Hsapiens.v86,
#'     upstream = 2000L,
#'     downstream = 2000L
#' )
#' result$plot
#' 
#' ## Example 3: Multiple peak sets with weights
#' peaks_list <- GRangesList(
#'     "Condition1" = peaks1,
#'     "Condition2" = peaks2
#' )
#' result <- metageneSignalProfile(
#'     peaks = peaks_list,
#'     AnnotationData = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'     peakWeight = "score"
#' )
#' result$plot
#' }
metageneSignalProfile <- function(bigwigFiles = NULL,
                                  peaks = NULL,
                                  AnnotationData,
                                  upstream = 2000L,
                                  downstream = 2000L,
                                  geneBodyLength = 5000L,
                                  upstreamBins = 20L,
                                  geneBodyBins = 50L,
                                  downstreamBins = 20L,
                                  minGeneLength = NULL,
                                  peakWeight = NULL,
                                  saveData = TRUE,
                                  dataFile = NULL,
                                  BPPARAM = NULL,
                                  verbose = TRUE,
                                  ...) {
    # Step 1: Validate input
    if (is.null(bigwigFiles) && is.null(peaks)) {
        stop("Either 'bigwigFiles' or 'peaks' must be provided", call. = FALSE)
    }
    if (!is.null(bigwigFiles) && !is.null(peaks)) {
        stop("Only one of 'bigwigFiles' or 'peaks' should be provided",
             call. = FALSE)
    }
    
    # Setup BiocParallel backend
    if (is.null(BPPARAM)) {
        if (requireNamespace("BiocParallel", quietly = TRUE)) {
            BPPARAM <- BiocParallel::bpparam()
        } else {
            if (verbose) {
                message("BiocParallel not available, using sequential processing")
            }
            BPPARAM <- NULL  # Will use sequential processing
        }
    }
    
    input_type <- if (!is.null(bigwigFiles)) "bigwig" else "peaks"
    
    if (input_type == "bigwig") {
        if (!is.character(bigwigFiles)) {
            stop("'bigwigFiles' must be a character vector", call. = FALSE)
        }
        if (length(bigwigFiles) == 0L) {
            stop("'bigwigFiles' must contain at least one file path",
                 call. = FALSE)
        }
        for (bw in bigwigFiles) {
            if (!file.exists(bw)) {
                stop("BigWig file not found: ", bw, call. = FALSE)
            }
        }
        input_names <- basename(bigwigFiles)
    } else {
        if (!inherits(peaks, c("GRanges", "GRangesList"))) {
            stop("'peaks' must be a GRanges or GRangesList object",
                 call. = FALSE)
        }
        if (inherits(peaks, "GRanges")) {
            n <- deparse(substitute(peaks))
            peaks <- GRangesList(peaks)
            names(peaks) <- if (n != "peaks") n else "peaks"
        }
        if (length(peaks) == 0L) {
            stop("'peaks' must contain at least one peak set", call. = FALSE)
        }
        input_names <- names(peaks)
        if (is.null(input_names) || any(input_names == "")) {
            input_names <- paste0("peakSet", seq_along(peaks))
            names(peaks) <- input_names
        }
    }
    
    # Validate AnnotationData
    if (!inherits(AnnotationData, c("TxDb", "EnsDb"))) {
        stop("'AnnotationData' must be a TxDb or EnsDb object", call. = FALSE)
    }
    
    # Validate numeric parameters
    stopifnot(is.numeric(upstream), upstream >= 0L)
    stopifnot(is.numeric(downstream), downstream >= 0L)
    stopifnot(is.numeric(geneBodyLength), geneBodyLength > 0L)
    stopifnot(is.numeric(upstreamBins), upstreamBins > 0L)
    stopifnot(is.numeric(geneBodyBins), geneBodyBins > 0L)
    stopifnot(is.numeric(downstreamBins), downstreamBins > 0L)
    
    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    geneBodyLength <- as.integer(geneBodyLength)
    upstreamBins <- as.integer(upstreamBins)
    geneBodyBins <- as.integer(geneBodyBins)
    downstreamBins <- as.integer(downstreamBins)
    
    # Step 2: Extract genes from annotation
    if (inherits(AnnotationData, "TxDb")) {
        suppressMessages(genes_gr <- genes(AnnotationData))
    } else {
        # EnsDb
        genes_gr <- EnsDb2GR(AnnotationData, feature = "gene")
    }
    
    if (length(genes_gr) == 0L) {
        stop("No genes found in AnnotationData", call. = FALSE)
    }
    
    # Ensure genes have names
    if (is.null(names(genes_gr)) || any(names(genes_gr) == "")) {
        names(genes_gr) <- paste0("gene", seq_along(genes_gr))
    }
    
    # Handle seqlevels compatibility
    if (input_type == "peaks") {
        seql_style <- seqlevelsStyle(genes_gr)
        peaks <- lapply(peaks, function(p) {
            seqlevelsStyle(p) <- seql_style[1L]
            p
        })
        peaks <- GRangesList(peaks)
    }
    
    # Filter by minimum gene length if specified
    if (!is.null(minGeneLength)) {
        stopifnot(is.numeric(minGeneLength), minGeneLength > 0L)
        minGeneLength <- as.integer(minGeneLength)
        gene_lengths <- width(genes_gr)
        genes_gr <- genes_gr[gene_lengths >= minGeneLength]
        if (length(genes_gr) == 0L) {
            stop("No genes remaining after filtering by minGeneLength",
                 call. = FALSE)
        }
    }
    
    # Calculate TSS and TES positions
    genes_gr <- .calculateTSSandTES(genes_gr)
    
    # Store original gene lengths
    original_lengths <- width(genes_gr)
    scaling_factors <- data.frame(
        gene_id = names(genes_gr),
        original_length = original_lengths,
        scaled_length = geneBodyLength,
        scaling_factor = geneBodyLength / original_lengths,
        stringsAsFactors = FALSE
    )
    
    # Step 3: Scale gene bodies to fixed length
    genes_scaled <- .scaleGeneBodies(genes_gr, geneBodyLength)
    
    # Step 4: Create metagene regions
    metagene_regions <- .createMetageneRegions(
        genes_scaled, genes_gr, upstream, downstream,
        upstreamBins, geneBodyBins, downstreamBins
    )
    
    # Step 5: Handle input (BigWig or peaks)
    if (input_type == "bigwig") {
        if (verbose) {
            message("Extracting signal from ", length(bigwigFiles),
                   " BigWig file(s) for ", length(genes_gr), " genes...")
        }
        # For BigWig: Extract signal directly from files per gene (memory efficient)
        signal_matrices <- .extractSignalFromBigWigFiles(
            bigwigFiles, metagene_regions, genes_gr, genes_scaled,
            upstream, downstream, upstreamBins, geneBodyBins, downstreamBins,
            BPPARAM = BPPARAM, verbose = verbose
        )
    } else {
        if (verbose) {
            message("Calculating peak coverage for ", length(genes_gr),
                   " genes from ", length(peaks), " peaks...")
        }
        # Calculate peak coverage
        coverage_list <- .calculatePeakCoverage(
            peaks, metagene_regions, peakWeight
        )
        # Step 6: Extract signal/coverage from scaled regions
        signal_matrices <- .extractSignalFromRegions(
            coverage_list, metagene_regions, genes_gr, genes_scaled,
            upstream, downstream, upstreamBins, geneBodyBins, downstreamBins,
            input_type
        )
    }
    
    # Step 7: Calculate average signal/coverage (before saving data)
    average_signal <- .calculateAverageSignal(signal_matrices, input_names)
    
    # Step 8: Save data to RDS file if requested
    data_file_path <- NULL
    data_saved <- FALSE
    if (saveData) {
        if (verbose) {
            message("Saving data matrices to RDS file...")
        }
        if (is.null(dataFile)) {
            # Generate timestamped filename
            timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
            data_file_path <- paste0("metageneSignalProfile_data_", timestamp, ".rds")
        } else {
            data_file_path <- dataFile
        }
        
        tryCatch({
            saveRDS(signal_matrices, file = data_file_path, compress = TRUE)
            if (verbose) {
                message("Data matrices saved to: ", data_file_path)
            }
            data_saved <- TRUE
            # Clear signal_matrices from memory after saving if successful
            rm(signal_matrices)
            gc(verbose = FALSE)
        }, error = function(e) {
            warning("Failed to save data to RDS file: ", e$message,
                   "\nData will be returned in result instead.", call. = FALSE)
            data_file_path <- NULL
            data_saved <- FALSE
        })
    }
    
    # Step 9: Generate plot
    plot_obj <- .generateMetagenePlot(
        average_signal, upstream, downstream, geneBodyLength,
        upstreamBins, geneBodyBins, downstreamBins, input_type, input_names
    )
    
    # Step 10: Create bin info
    bin_info <- .createBinInfo(
        upstream, downstream, geneBodyLength,
        upstreamBins, geneBodyBins, downstreamBins
    )
    
    # Step 11: Return results
    result <- list(
        plot = plot_obj,
        average = average_signal,
        regions = metagene_regions$scaled_genes,
        binInfo = bin_info,
        scalingFactors = scaling_factors,
        metadata = list(
            inputType = input_type,
            inputNames = input_names,
            upstream = upstream,
            downstream = downstream,
            geneBodyLength = geneBodyLength,
            upstreamBins = upstreamBins,
            geneBodyBins = geneBodyBins,
            downstreamBins = downstreamBins,
            nGenes = length(genes_gr),
            normalization = if (input_type == "bigwig") "assumed normalized" else "peak coverage",
            dataFile = data_file_path,
            dataSaved = data_saved
        )
    )
    
    # Only include data in result if not saved
    if (!data_saved) {
        result$data <- signal_matrices
    }
    
    if (verbose) {
        message("Metagene signal profile analysis complete!")
    }
    
    return(result)
}

# Helper function: Calculate TSS and TES positions
.calculateTSSandTES <- function(genes_gr) {
    # TSS is start for + strand, end for - strand
    # TES is end for + strand, start for - strand
    strand_info <- as.character(strand(genes_gr))
    
    tss_pos <- ifelse(strand_info == "+", start(genes_gr), end(genes_gr))
    tes_pos <- ifelse(strand_info == "+", end(genes_gr), start(genes_gr))
    
    # Store in metadata
    mcols(genes_gr)$TSS <- tss_pos
    mcols(genes_gr)$TES <- tes_pos
    
    return(genes_gr)
}

# Helper function: Scale gene bodies to fixed length
.scaleGeneBodies <- function(genes_gr, geneBodyLength) {
    genes_scaled <- genes_gr
    
    strand_info <- as.character(strand(genes_gr))
    tss_pos <- mcols(genes_gr)$TSS
    tes_pos <- mcols(genes_gr)$TES
    
    # For positive strand: keep TSS, extend TES
    # For negative strand: keep TES, extend TSS (backwards)
    for (i in seq_along(genes_scaled)) {
        if (strand_info[i] == "+") {
            start(genes_scaled[i]) <- tss_pos[i]
            end(genes_scaled[i]) <- tss_pos[i] + geneBodyLength - 1L
        } else {
            end(genes_scaled[i]) <- tes_pos[i]
            start(genes_scaled[i]) <- tes_pos[i] - geneBodyLength + 1L
        }
    }
    
    # Update TSS and TES in metadata
    mcols(genes_scaled)$TSS <- tss_pos
    mcols(genes_scaled)$TES <- ifelse(strand_info == "+",
                                     tss_pos + geneBodyLength - 1L,
                                     tes_pos - geneBodyLength + 1L)
    
    return(genes_scaled)
}

# Helper function: Create metagene regions
.createMetageneRegions <- function(genes_scaled, genes_original,
                                   upstream, downstream,
                                   upstreamBins, geneBodyBins, downstreamBins) {
    # Create upstream regions (from original TSS)
    upstream_regions <- suppressWarnings(
        promoters(genes_original, upstream = upstream, downstream = 0L)
    )
    upstream_regions <- trim(upstream_regions)
    
    # Create downstream regions (from original TES)
    downstream_regions <- suppressWarnings(
        downstreams(genes_original, upstream = 0L, downstream = downstream)
    )
    downstream_regions <- trim(downstream_regions)
    
    # Gene body regions (scaled)
    gene_body_regions <- genes_scaled
    
    # Also need original gene body regions for signal extraction
    gene_body_original <- genes_original
    
    return(list(
        upstream = upstream_regions,
        gene_body_scaled = gene_body_regions,
        gene_body_original = gene_body_original,
        downstream = downstream_regions,
        scaled_genes = genes_scaled
    ))
}

# Helper function: Extract signal directly from BigWig files (memory efficient)
# Processes genes chromosome by chromosome to balance memory and speed
# Supports parallel processing across chromosomes using BiocParallel
.extractSignalFromBigWigFiles <- function(bigwigFiles, metagene_regions,
                                          genes_original, genes_scaled,
                                          upstream, downstream,
                                          upstreamBins, geneBodyBins,
                                          downstreamBins,
                                          BPPARAM = NULL,
                                          verbose = TRUE) {
    n_genes <- length(genes_original)
    total_bins <- upstreamBins + geneBodyBins + downstreamBins
    
    # Initialize signal matrices for each BigWig file
    signal_matrices <- lapply(bigwigFiles, function(bw) {
        matrix(0, nrow = n_genes, ncol = total_bins)
    })
    names(signal_matrices) <- basename(bigwigFiles)
    
    # Group genes by chromosome for efficient processing
    genes_by_chr <- split(seq_len(n_genes), seqnames(genes_original))
    chr_names <- names(genes_by_chr)
    n_chr <- length(chr_names)
    
    if (verbose && n_chr > 1L) {
        if (!is.null(BPPARAM) && requireNamespace("BiocParallel", quietly = TRUE)) {
            if (inherits(BPPARAM, "SerialParam")) {
                message("Processing ", n_chr, " chromosomes (sequential)")
            } else {
                workers <- BiocParallel::bpnworkers(BPPARAM)
                message("Processing ", n_chr, " chromosomes (parallel, ", workers, " workers)")
            }
        } else {
            message("Processing ", n_chr, " chromosomes (sequential)")
        }
    }
    
    # Helper function to process a single chromosome
    process_chromosome <- function(chr) {
        gene_indices <- genes_by_chr[[chr]]
        n_genes_chr <- length(gene_indices)
        
        if (n_genes_chr == 0L) {
            return(list(chr = chr, matrices = NULL))
        }
        
        # Get all regions for genes on this chromosome
        chr_upstream <- metagene_regions$upstream[gene_indices]
        chr_gene_body_orig <- metagene_regions$gene_body_original[gene_indices]
        chr_downstream <- metagene_regions$downstream[gene_indices]
        chr_genes_scaled <- genes_scaled[gene_indices]
        chr_genes_original <- genes_original[gene_indices]
        
        # Combine all regions for this chromosome
        chr_all_regions <- c(chr_upstream, chr_gene_body_orig, chr_downstream)
        chr_all_regions <- chr_all_regions[width(chr_all_regions) > 0]
        
        if (length(chr_all_regions) == 0L) {
            return(list(chr = chr, matrices = NULL))
        }
        
        # Reduce overlapping regions to minimize data read
        chr_regions_reduced <- reduce(chr_all_regions)
        
        # Initialize chromosome-specific matrices
        chr_matrices <- lapply(bigwigFiles, function(bw) {
            matrix(0, nrow = n_genes_chr, ncol = total_bins)
        })
        
        # Read signal for this chromosome from each BigWig file
        for (bw_idx in seq_along(bigwigFiles)) {
            bw <- bigwigFiles[bw_idx]
            tryCatch({
                # Read only this chromosome's regions
                cvg_chr <- import(bw, format = "BigWig",
                                 which = chr_regions_reduced,
                                 as = "RleList")
                
                # Extract signal for each gene on this chromosome
                for (local_idx in seq_len(n_genes_chr)) {
                    # Extract signal for each region of this gene
                    upstream_signal <- .extractRegionSignalFromRleList(
                        cvg_chr, chr_upstream[local_idx],
                        upstreamBins, chr_genes_original[local_idx]
                    )
                    
                    gene_body_signal <- .extractGeneBodySignalFromRleList(
                        cvg_chr, chr_gene_body_orig[local_idx],
                        chr_genes_scaled[local_idx], geneBodyBins
                    )
                    
                    downstream_signal <- .extractRegionSignalFromRleList(
                        cvg_chr, chr_downstream[local_idx],
                        downstreamBins, chr_genes_original[local_idx]
                    )
                    
                    # Store in chromosome-specific matrix
                    chr_matrices[[bw_idx]][local_idx, ] <- c(
                        upstream_signal, gene_body_signal, downstream_signal
                    )
                }
                
                # Clear chromosome coverage from memory
                rm(cvg_chr)
                gc(verbose = FALSE)
                
            }, error = function(e) {
                warning("Failed to read signal for chromosome ", chr,
                       " from ", basename(bw), ": ", e$message,
                       call. = FALSE)
                # Leave as zeros for genes on this chromosome
            })
        }
        
        return(list(chr = chr, gene_indices = gene_indices,
                  matrices = chr_matrices))
    }
    
    # Process chromosomes (parallel or sequential) using BiocParallel
    if (!is.null(BPPARAM) && requireNamespace("BiocParallel", quietly = TRUE) && 
        n_chr > 1L && !inherits(BPPARAM, "SerialParam")) {
        # Parallel processing with BiocParallel
        chr_results <- BiocParallel::bplapply(chr_names, process_chromosome,
                                              BPPARAM = BPPARAM)
    } else {
        # Sequential processing
        chr_results <- lapply(chr_names, process_chromosome)
    }
    
    # Combine results from all chromosomes
    for (result in chr_results) {
        if (is.null(result$matrices)) next
        
        chr <- result$chr
        gene_indices <- result$gene_indices
        
        for (bw_idx in seq_along(bigwigFiles)) {
            signal_matrices[[bw_idx]][gene_indices, ] <- result$matrices[[bw_idx]]
        }
    }
    
    # Clean up
    rm(chr_results)
    gc(verbose = FALSE)
    
    # Set row names
    for (i in seq_along(signal_matrices)) {
        rownames(signal_matrices[[i]]) <- names(genes_original)
    }
    
    return(signal_matrices)
}

# Helper function: Extract signal from RleList for a region
.extractRegionSignalFromRleList <- function(cvg, region, n_bins, gene_original) {
    if (width(region) == 0) {
        return(rep(0, n_bins))
    }
    
    region_seqname <- as.character(seqnames(region))[1]
    if (region_seqname %in% names(cvg)) {
        region_rle <- cvg[[region_seqname]]
        
        # Bin the signal
        bins <- tile(region, n = n_bins)
        bin_signals <- sapply(bins, function(bin) {
            bin_views <- Views(region_rle, ranges(bin))
            mean(viewMeans(bin_views), na.rm = TRUE)
        })
        bin_signals[is.na(bin_signals)] <- 0
        return(bin_signals)
    }
    return(rep(0, n_bins))
}

# Helper function: Extract and interpolate gene body signal from RleList
.extractGeneBodySignalFromRleList <- function(cvg, gene_original, gene_scaled,
                                               n_bins) {
    if (width(gene_original) == 0) {
        return(rep(0, n_bins))
    }
    
    region_seqname <- as.character(seqnames(gene_original))[1]
    if (region_seqname %in% names(cvg)) {
        region_rle <- cvg[[region_seqname]]
        
        # Extract signal values from the original gene body region
        gene_width <- width(gene_original)
        
        if (gene_width > 0) {
            # Create tiles for the gene body to extract signal
            n_tiles <- min(gene_width, 200L)
            gene_tiles <- tile(gene_original, n = n_tiles)
            
            # Extract mean signal per tile
            tile_views <- Views(region_rle, ranges(gene_tiles))
            tile_signals <- viewMeans(tile_views)
            
            # Expand to full gene width by interpolating
            if (n_tiles < gene_width) {
                tile_positions <- seq(1, gene_width, length.out = n_tiles)
                full_positions <- seq_len(gene_width)
                original_signal <- approx(tile_positions, as.numeric(tile_signals),
                                        xout = full_positions,
                                        method = "linear", rule = 2)$y
                original_signal[is.na(original_signal)] <- 0
            } else {
                original_signal <- as.numeric(tile_signals)[seq_len(gene_width)]
            }
        } else {
            original_signal <- numeric(0)
        }
        
        # Interpolate to scaled coordinates
        scaled_width <- width(gene_scaled)
        if (length(original_signal) > 1 && gene_width > 1 && scaled_width > 0) {
            # Map positions: original gene body to scaled gene body
            original_positions <- seq_len(gene_width)
            scaled_positions <- seq(1, scaled_width, length.out = n_bins)
            
            # Map scaled positions back to original positions
            mapped_positions <- (scaled_positions / scaled_width) * gene_width
            mapped_positions <- pmax(1, pmin(gene_width, mapped_positions))
            
            # Interpolate signal
            interpolated <- approx(original_positions, original_signal,
                                 xout = mapped_positions,
                                 method = "linear", rule = 2)$y
            interpolated[is.na(interpolated)] <- 0
            return(interpolated)
        } else if (length(original_signal) > 0) {
            mean_signal <- mean(original_signal, na.rm = TRUE)
            return(rep(if (is.na(mean_signal)) 0 else mean_signal, n_bins))
        }
    }
    return(rep(0, n_bins))
}

# Helper function: Calculate peak coverage
.calculatePeakCoverage <- function(peaks, metagene_regions, peakWeight) {
    # Get weights if provided
    weights_list <- .extractPeakWeights(peaks, peakWeight)
    
    # For each peak set, calculate coverage
    coverage_list <- lapply(seq_along(peaks), function(i) {
        peak_set <- peaks[[i]]
        weights <- weights_list[[i]]
        
        # Calculate coverage as RleList
        # This is a simplified approach - we'll calculate per-bin coverage later
        # For now, return the peaks and weights for later processing
        return(list(peaks = peak_set, weights = weights))
    })
    
    names(coverage_list) <- names(peaks)
    return(coverage_list)
}

# Helper function: Extract peak weights
.extractPeakWeights <- function(peaks, peakWeight) {
    if (is.null(peakWeight)) {
        # Uniform weights
        return(lapply(peaks, function(p) rep(1, length(p))))
    }
    
    weights_list <- lapply(peaks, function(peak_set) {
        if (is.character(peakWeight)) {
            # Extract from metadata column
            if (!peakWeight %in% colnames(mcols(peak_set))) {
                stop("Column '", peakWeight, "' not found in peaks metadata",
                     call. = FALSE)
            }
            return(mcols(peak_set)[[peakWeight]])
        } else if (is.numeric(peakWeight)) {
            # Use provided vector
            if (length(peakWeight) != length(peak_set)) {
                stop("Length of 'peakWeight' must match length of peaks",
                     call. = FALSE)
            }
            return(peakWeight)
        } else {
            stop("'peakWeight' must be a character string or numeric vector",
                 call. = FALSE)
        }
    })
    
    return(weights_list)
}

# Helper function: Extract signal from regions
.extractSignalFromRegions <- function(coverage_list, metagene_regions,
                                     genes_original, genes_scaled,
                                     upstream, downstream,
                                     upstreamBins, geneBodyBins, downstreamBins,
                                     input_type) {
    n_genes <- length(genes_original)
    total_bins <- upstreamBins + geneBodyBins + downstreamBins
    
    signal_matrices <- lapply(seq_along(coverage_list), function(i) {
        cvg <- coverage_list[[i]]
        
        # Initialize matrix: rows = genes, columns = bins
        signal_matrix <- matrix(0, nrow = n_genes, ncol = total_bins)
        
        for (gene_idx in seq_len(n_genes)) {
            # Extract signal for each region
            upstream_signal <- .extractRegionSignal(
                cvg, metagene_regions$upstream[gene_idx],
                upstreamBins, input_type, genes_original[gene_idx]
            )
            
            gene_body_signal <- .extractGeneBodySignal(
                cvg, genes_original[gene_idx], genes_scaled[gene_idx],
                geneBodyBins, input_type
            )
            
            downstream_signal <- .extractRegionSignal(
                cvg, metagene_regions$downstream[gene_idx],
                downstreamBins, input_type, genes_original[gene_idx]
            )
            
            # Combine signals
            signal_matrix[gene_idx, ] <- c(upstream_signal, gene_body_signal,
                                          downstream_signal)
        }
        
        rownames(signal_matrix) <- names(genes_original)
        return(signal_matrix)
    })
    
    names(signal_matrices) <- names(coverage_list)
    return(signal_matrices)
}

# Helper function: Extract signal from a region (upstream/downstream)
.extractRegionSignal <- function(cvg, region, n_bins, input_type, gene_original) {
    if (input_type == "bigwig") {
        # Extract from RleList
        region_seqname <- as.character(seqnames(region))[1]
        if (region_seqname %in% names(cvg)) {
            region_rle <- cvg[[region_seqname]]
            
            # Bin the signal
            if (width(region) > 0) {
                bins <- tile(region, n = n_bins)
                bin_signals <- sapply(bins, function(bin) {
                    bin_views <- Views(region_rle, ranges(bin))
                    mean(viewMeans(bin_views), na.rm = TRUE)
                })
                return(bin_signals)
            }
        }
        return(rep(0, n_bins))
    } else {
        # Calculate peak coverage
        peaks <- cvg$peaks
        weights <- cvg$weights
        
        if (width(region) > 0) {
            bins <- tile(region, n = n_bins)
            bin_counts <- sapply(bins, function(bin) {
                ol <- findOverlaps(bin, peaks)
                if (length(ol) > 0) {
                    sum(weights[subjectHits(ol)])
                } else {
                    0
                }
            })
            # Normalize by bin size (density)
            bin_sizes <- width(bins)
            bin_density <- bin_counts / (bin_sizes / 1000)  # per kb
            return(bin_density)
        }
        return(rep(0, n_bins))
    }
}

# Helper function: Extract and interpolate gene body signal
.extractGeneBodySignal <- function(cvg, gene_original, gene_scaled,
                                   n_bins, input_type) {
    if (input_type == "bigwig") {
        # Extract signal from original gene body
        region_seqname <- as.character(seqnames(gene_original))[1]
        if (region_seqname %in% names(cvg)) {
            region_rle <- cvg[[region_seqname]]
            
            # Extract signal values from the original gene body region
            gene_width <- width(gene_original)
            
            if (gene_width > 0) {
                # Create tiles for the gene body to extract signal
                # Use finer tiles for better interpolation
                n_tiles <- min(gene_width, 200L)
                gene_tiles <- tile(gene_original, n = n_tiles)
                
                # Extract mean signal per tile
                tile_views <- Views(region_rle, ranges(gene_tiles))
                tile_signals <- viewMeans(tile_views)
                
                # Expand to full gene width by interpolating
                if (n_tiles < gene_width) {
                    tile_positions <- seq(1, gene_width, length.out = n_tiles)
                    full_positions <- seq_len(gene_width)
                    original_signal <- approx(tile_positions, as.numeric(tile_signals),
                                            xout = full_positions,
                                            method = "linear", rule = 2)$y
                    original_signal[is.na(original_signal)] <- 0
                } else {
                    original_signal <- as.numeric(tile_signals)[seq_len(gene_width)]
                }
            } else {
                original_signal <- numeric(0)
            }
            
            # Interpolate to scaled coordinates
            scaled_width <- width(gene_scaled)
            if (length(original_signal) > 1 && gene_width > 1 && scaled_width > 0) {
                # Map positions: original gene body to scaled gene body
                original_positions <- seq_len(gene_width)
                scaled_positions <- seq(1, scaled_width, length.out = n_bins)
                
                # Map scaled positions back to original positions
                # scaled_pos / scaled_width corresponds to original_pos / gene_width
                mapped_positions <- (scaled_positions / scaled_width) * gene_width
                mapped_positions <- pmax(1, pmin(gene_width, mapped_positions))
                
                # Interpolate signal
                interpolated <- approx(original_positions, original_signal,
                                     xout = mapped_positions,
                                     method = "linear", rule = 2)$y
                interpolated[is.na(interpolated)] <- 0
            } else if (length(original_signal) > 0) {
                # Single value or very short gene - use mean
                mean_signal <- mean(original_signal, na.rm = TRUE)
                interpolated <- rep(if (is.na(mean_signal)) 0 else mean_signal,
                                   n_bins)
            } else {
                # No signal
                interpolated <- rep(0, n_bins)
            }
            
            return(interpolated)
        }
        return(rep(0, n_bins))
    } else {
        # Calculate peak coverage for gene body
        peaks <- cvg$peaks
        weights <- cvg$weights
        
        # Create bins in original gene body
        original_bins <- tile(gene_original, n = n_bins)
        
        # Calculate coverage per bin
        bin_counts <- sapply(original_bins, function(bin) {
            ol <- findOverlaps(bin, peaks)
            if (length(ol) > 0) {
                sum(weights[subjectHits(ol)])
            } else {
                0
            }
        })
        
        # Normalize by bin size
        bin_sizes <- width(original_bins)
        bin_density <- bin_counts / (bin_sizes / 1000)  # per kb
        
        # Interpolate to scaled coordinates (already binned, so just return)
        return(bin_density)
    }
}

# Helper function: Calculate average signal
.calculateAverageSignal <- function(signal_matrices, input_names) {
    n_bins <- ncol(signal_matrices[[1]])
    
    average_df <- data.frame(
        bin = seq_len(n_bins),
        stringsAsFactors = FALSE
    )
    
    for (i in seq_along(signal_matrices)) {
        avg_signal <- colMeans(signal_matrices[[i]], na.rm = TRUE)
        average_df[[input_names[i]]] <- avg_signal
    }
    
    return(average_df)
}

# Helper function: Generate plot
.generateMetagenePlot <- function(average_signal, upstream, downstream,
                                 geneBodyLength, upstreamBins, geneBodyBins,
                                 downstreamBins, input_type, input_names) {
    # Create x-axis positions
    upstream_positions <- seq(-upstream, 0, length.out = upstreamBins)
    gene_body_positions <- seq(0, geneBodyLength, length.out = geneBodyBins + 1)
    gene_body_positions <- gene_body_positions[-length(gene_body_positions)]
    downstream_start <- geneBodyLength
    downstream_positions <- seq(downstream_start,
                                downstream_start + downstream,
                                length.out = downstreamBins + 1)
    downstream_positions <- downstream_positions[-1]
    
    x_positions <- c(upstream_positions, gene_body_positions,
                    downstream_positions)
    
    # Prepare data for plotting
    plot_data <- data.frame(
        position = x_positions,
        stringsAsFactors = FALSE
    )
    
    for (name in input_names) {
        plot_data[[name]] <- average_signal[[name]]
    }
    
    # Reshape for ggplot
    plot_data_long <- melt(plot_data, id.vars = "position",
                          variable.name = "sample",
                          value.name = "signal")
    
    # Create plot
    p <- ggplot(plot_data_long, aes_string(x = "position", y = "signal",
                                           color = "sample")) +
        geom_line() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
        geom_vline(xintercept = geneBodyLength, linetype = "dashed",
                  color = "gray") +
        theme_bw() +
        labs(
            x = "Position (bp)",
            y = if (input_type == "bigwig") "Normalized Signal" else "Peak Coverage Density",
            title = "Metagene Profile",
            color = "Sample"
        )
    
    return(p)
}

# Helper function: Create bin info
.createBinInfo <- function(upstream, downstream, geneBodyLength,
                          upstreamBins, geneBodyBins, downstreamBins) {
    upstream_positions <- seq(-upstream, 0, length.out = upstreamBins)
    gene_body_positions <- seq(0, geneBodyLength, length.out = geneBodyBins + 1)
    gene_body_positions <- gene_body_positions[-length(gene_body_positions)]
    downstream_start <- geneBodyLength
    downstream_positions <- seq(downstream_start,
                               downstream_start + downstream,
                               length.out = downstreamBins + 1)
    downstream_positions <- downstream_positions[-1]
    
    bin_info <- data.frame(
        bin = seq_len(upstreamBins + geneBodyBins + downstreamBins),
        region = c(rep("upstream", upstreamBins),
                  rep("gene_body", geneBodyBins),
                  rep("downstream", downstreamBins)),
        position = c(upstream_positions, gene_body_positions,
                    downstream_positions),
        stringsAsFactors = FALSE
    )
    
    return(bin_info)
}

