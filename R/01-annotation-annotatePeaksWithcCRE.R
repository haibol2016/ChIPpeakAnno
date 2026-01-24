#' Annotate peaks with ENCODE candidate cis-Regulatory Elements (cCREs)
#' 
#' @description 
#' Annotates genomic peaks by finding overlaps with ENCODE candidate 
#' cis-Regulatory Elements (cCREs). cCREs are genomic regions identified by 
#' ENCODE that show biochemical signatures of regulatory activity. The full list
#' of cCRE signatures includes: CA-CTCF (chromatin accessibility with CTCF binding),
#' CA-H3K4me3 (chromatin accessibility with H3K4me3 histone modification), 
#' CA-only (chromatin accessibility only), 
#' CA-TF (chromatin accessibility with transcription factor binding), 
#' dELS (distal Enhancer-like Sequence), pELS (proximal Enhancer-like Sequence), 
#' PLS (Promoter-like Sequence), and TF (Transcription Factor). In tissue-specific
#' cCRE downloads, a "Low-DNase" signature is included. In the all-in-one cCRE 
#' download, there is a "CA" signature that includes: CA-CTCF, CA-H3K4me3, CA-only,
#' CA-TF. The Silencer sets and MAFF/MAFK dynamic enhancers are not included in 
#' the cCRE download.
#' 
#' \strong{Note:} "Low-DNase" signatures are automatically excluded from annotation 
#' as they are ubiquitous and not informative.
#' 
#' This function can download cCRE data from ENCODE/SCREEN based on species 
#' and tissue type, or accept a pre-downloaded file. It uses \code{findOverlaps} 
#' to identify peaks that overlap with cCREs and adds annotation metadata to 
#' the peaks.
#' 
#' @param peaks A \code{\link[GenomicRanges]{GRanges}} object containing 
#'        peaks to be annotated. Peaks can be pre-annotated by 
#'        \code{\link{annotatePeakInBatch}} or other annotation functions.
#' @param species Character string specifying the species. Must be one of:
#'        \itemize{
#'          \item \code{"Homo sapiens"}: Human cCREs
#'          \item \code{"Mus musculus"}: Mouse cCREs
#'        }
#'        If \code{NULL} (default), uses the global species option set by
#'        \code{\link{setChIPpeakAnnoSpecies}}. If no global option is set,
#'        an error is raised prompting the user to specify species.
#' @param tissue_types Character vector specifying tissue types. If 
#'        \code{"any"} (default), downloads "All Human cCREs" or "All Mouse 
#'        cCREs" depending on species. Otherwise, downloads tissue-specific 
#'        cCREs. Use \code{\link{listAvailablecCREs}} to see available 
#'        tissue types.
#' @param cCRE_file Optional character string. If provided, uses this file path 
#'        or URL instead of downloading from ENCODE. Can be a local file path 
#'        or a URL to a BED file. This takes precedence over \code{species} 
#'        and \code{tissue_types}.
#' @param cache_dir Optional character string specifying a directory for 
#'        caching downloaded cCRE files. If \code{NULL} (default), uses 
#'        \code{tempdir()}. Files are cached to avoid re-downloading.
#' @param ... Additional arguments passed to \code{\link[IRanges]{findOverlaps}}, 
#'        such as \code{maxgap}, \code{minoverlap}, \code{type}, etc.
#' 
#' @return A \code{\link[GenomicRanges]{GRanges}} object containing the input 
#'        peaks with additional metadata columns:
#'        \itemize{
#'          \item \code{cCRE_count}: Integer, number of overlapping cCREs
#'          \item \code{cCRE_coordinates}: Character, comma-separated list of 
#'                overlapping cCRE coordinates in standard format (e.g., 
#'                "chr1:100-1000,chr2:5000-6000")
#'          \item \code{cCRE_types}: Character, comma-separated list of cCRE 
#'                types (e.g., "PLS", "pELS", "dELS", "CA-CTCF", "CA-H3K4me3",
#'                "CA-only", "CA-TF") for overlapping cCREs. Note: "Low-DNase"
#'                signatures are automatically excluded as they are ubiquitous and
#'                not informative.
#'          \item \code{cCRE_accession}: Character, comma-separated list of 
#'                ENCODE accession IDs of overlapping cCREs (if available)
#'        }
#'        Peaks without overlaps will have \code{cCRE_count = 0} and \code{NA} 
#'        for other cCRE-related columns.
#' 
#' @details
#' 
#' \strong{ENCODE cCRE Data:}
#' ENCODE cCREs are identified based on biochemical signatures from various 
#' assays (ChIP-seq, DNase-seq, ATAC-seq, etc.). The complete cCRE catalog 
#' includes the following signature types:
#' \itemize{
#'   \item \strong{CA-CTCF}: Candidate cis-regulatory element marked by chromatin accessibility and CTCF binding
#'   \item \strong{CA-H3K4me3}: Candidate cis-regulatory element marked by chromatin accessibility and H3K4me3 
#'         histone modification
#'   \item \strong{CA-only}: Candidate cis-regulatory element marked by chromatin accessibility only
#'   \item \strong{CA-TF}: Candidate cis-regulatory element marked by chromatin accessibility and transcription factor binding
#'   \item \strong{dELS}: Distal enhancer-like regions far from genes
#'   \item \strong{pELS}: Proximal enhancer-like regions near genes
#'   \item \strong{PLS}: Promoter-like sequence
#'   \item \strong{TF}: Transcription factor binding sites
#'   \item \strong{Low-DNase}: Regions with low DNase signal (automatically excluded from annotation as ubiquitous and not informative)
#' }
#' 
#' \strong{Note:} "Low-DNase" signatures are automatically excluded from 
#' annotation because they are ubiquitous across the genome and provide little 
#' regulatory information. Only informative cCRE signatures are used for 
#' annotation.
#' 
#' \strong{Data Sources:}
#' \itemize{
#'   \item \strong{SCREEN}: https://screen.encodeproject.org/
#'   \item \strong{ENCODE Portal}: https://www.encodeproject.org/
#' }
#' 
#' \strong{Download Behavior:}
#' \itemize{
#'   \item If \code{cCRE_file} is provided, it takes precedence and no 
#'         download occurs
#'   \item If \code{tissue_types == "any"}, downloads the comprehensive 
#'         cCRE catalog for the species
#'   \item If specific tissue types are provided, downloads tissue-specific 
#'         cCRE files
#'   \item Downloaded files are cached in \code{cache_dir} to avoid 
#'         re-downloading
#' }
#' 
#' \strong{Overlap Detection:}
#' Uses \code{findOverlaps(peaks, ccres)} to identify overlaps. By default, 
#' any overlap is considered. Use \code{maxgap} parameter via \code{...} to 
#' allow gaps between peaks and cCREs.
#' 
#' @note
#' \itemize{
#'   \item \strong{Genome Assembly:} cCRE data downloaded from SCREEN is based on
#'         \strong{hg38} for human (\emph{Homo sapiens}) and \strong{mm10} 
#'         for mouse (\emph{Mus musculus}). If your peaks use a different genome
#'         assembly (e.g., mm39 for mouse), the function will automatically perform 
#'         coordinate conversion (liftOver) of cCRE data to match your peaks' genome assembly.
#'   \item \strong{File Format:} cCRE files from SCREEN do not strictly follow
#'         standard BED format. The all-in-one cCRE file format differs from
#'         tissue-specific files. The function uses robust parsing methods to
#'         handle these format variations.
#'   \item Large cCRE files (especially "All cCREs") may require significant 
#'         memory
#'   \item The function preserves all existing metadata columns from the 
#'         input peaks
#'   \item Downloaded files are cached to improve performance in subsequent 
#'         calls
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{listAvailablecCREs}} to discover available tissue types
#'   \item \code{\link{annotatePeakInBatch}} for gene-based annotation
#'   \item \code{\link{findOverlaps}} for overlap detection details
#' }
#' 
#' @author Haibo Liu
#' @keywords misc
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<- genome
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<-
#' @importFrom utils download.file
#' @importFrom rtracklayer import liftOver import.chain
#' @examples
#' 
#' # Example 1: Annotate with all human cCREs
#' \dontrun{
#' data("myPeakList")
#' setChIPpeakAnnoGlobals(genome = "hg38", species = "Homo sapiens")
#' annotated <- annotatePeaksWithcCRE(myPeakList, 
#'                                    species = "Homo sapiens",
#'                                    tissue_types = "any")
#' # Check how many peaks overlap with cCREs
#' table(annotated$cCRE_count > 0)
#' # View cCRE coordinates and types for overlapping peaks
#' overlapping <- annotated$cCRE_count > 0
#' head(annotated$cCRE_coordinates[overlapping])
#' head(annotated$cCRE_types[overlapping])
#' }
#' 
#' # Example 2: List available tissue types
#' \dontrun{
#' available <- listAvailablecCREs("Homo sapiens")
#' head(available)
#' }
annotatePeaksWithcCRE <- function(peaks,
                                  species = NULL,
                                  tissue_types = "any",
                                  cCRE_file = NULL,
                                  cache_dir = NULL,
                                  ...) {
    # Validate inputs
    if (!inherits(peaks, "GRanges")) {
        stop("'peaks' must be a GRanges object", call. = FALSE)
    }
    
    # Use global species option if species not provided
    if (is.null(species)) {
        species <- getChIPpeakAnnoSpecies()
        if (is.null(species)) {
            stop("'species' must be specified. Options: 'Homo sapiens' or 'Mus musculus'. ",
                 "You can set a global default using setChIPpeakAnnoSpecies().",
                 call. = FALSE)
        }
    }
    
    # Validate species
    species <- match.arg(species, c("Homo sapiens", "Mus musculus"))
    
    if (is.null(cache_dir)) {
        cache_dir <- tempdir()
    }
    
    # Detect target genome assembly using hybrid approach (global option + peak detection)
    target_genome <- getTargetGenome(peaks, species)
    
    # cCRE data from SCREEN is hg38 for human, mm10 for mouse
    # Note: mm39 cCRE data is not available from SCREEN, so mm39 peaks require liftOver from mm10
    if (species == "Homo sapiens") {
        ccre_genome <- "hg38"
        need_liftover <- !is.null(target_genome) && target_genome != "hg38"
    } else {
        # Mouse cCRE data is only available in mm10 format
        ccre_genome <- "mm10"
        need_liftover <- !is.null(target_genome) && target_genome != "mm10"
    }
    
    if (need_liftover) {
        message("Detected peaks on ", target_genome, " assembly. ",
               "cCRE data from SCREEN is on ", ccre_genome, ". Will lift over cCRE data to ", 
               target_genome, ".")
    } else if (is.null(target_genome)) {
        # getTargetGenome already checked global option, so if it's NULL here,
        # neither global option nor peak detection worked
        stop("Could not detect genome assembly from peaks and no global genome is set. ",
             "Please either:\n",
             "  1. Set genome metadata on peaks: genome(peaks) <- 'hg38' (or your assembly)\n",
             "  2. Set a global genome option: setChIPpeakAnnoGenome('hg38')\n",
             "  3. Use setChIPpeakAnnoGlobals() to set all global options at once",
             call. = FALSE)
    }
    
    # Load or download cCRE data
    if (!is.null(cCRE_file)) {
        # Use provided file
        message("Make sure the cCRE file is for the same species as the peaks.")
        message("Loading cCRE data from file: ", cCRE_file)
        ccres <- loadcCREsFromFile(cCRE_file)
    } else {
        # Download based on species and tissue type and load the cCRE data into
        # a GRanges object using the `loadcCREsFromFile()` function.
        if (length(tissue_types) == 1 && tissue_types == "any") {
            message("Downloading all cCREs for ", species)
            ccres <- downloadAllcCREs(species, cache_dir)
        } else {
            message("Downloading tissue-specific cCREs for ", 
                   paste(tissue_types, collapse = ", "))
            ccres <- downloadTissueSpecificcCREs(species, tissue_types, cache_dir)
        }
    }
    
    if (length(ccres) == 0) {
        stop("No cCRE data loaded. Returning peaks without cCRE annotation.",
             call. = FALSE)
    }
    
    # Lift over cCRE data if needed (following Roadmap implementation pattern)
    if (need_liftover) {
        message("Lifting over cCRE data from ", ccre_genome, " to ", target_genome, "...")
        chain <- getLiftOverChain(ccre_genome, target_genome, cache_dir)
        if (!is.null(chain)) {
            # liftOver returns a GRangesList (one element per input range)
            # Some ranges may not map, resulting in empty GRanges
            ccres_lifted <- liftOver(ccres, chain)
            # Convert to GRanges, keeping only ranges that successfully mapped
            ccres <- unlist(ccres_lifted)
            if (length(ccres) == 0) {
                stop("No cCRE regions could be lifted over to ", target_genome, 
                     ". This may indicate incompatible genome assemblies.", 
                     call. = FALSE)
            }
            message("Successfully lifted over ", length(ccres), 
                   " cCRE regions to ", target_genome, ".")
        } else {
            warning("Could not obtain liftOver chain file. ",
                   "Proceeding with ", ccre_genome, " cCRE data, which may not match your peaks.",
                   call. = FALSE)
        }
    }
    
    # Match seqlevels styles if needed
    if (length(intersect(seqlevelsStyle(peaks), seqlevelsStyle(ccres))) == 0) {
        message("Matching seqlevels styles between peaks and cCREs")
        seqlevelsStyle(ccres) <- seqlevelsStyle(peaks)[1]
    }
    
    # Find overlaps
    message("Finding overlaps between peaks and cCREs...")
    overlaps <- findOverlaps(peaks, ccres, ...)
    
    # Initialize metadata columns
    mcols(peaks)$cCRE_count <- 0L
    mcols(peaks)$cCRE_coordinates <- NA_character_
    mcols(peaks)$cCRE_types <- NA_character_
    mcols(peaks)$cCRE_accession <- NA_character_
    
    if (length(overlaps) > 0) {
        # Get peak and cCRE indices from overlaps
        peak_hits <- queryHits(overlaps)
        ccre_hits <- subjectHits(overlaps)
        
        # Pre-extract metadata columns and sequence info to avoid repeated 
        # indexing (performance optimization)
        ccre_mcols <- mcols(ccres)
        ccre_seqnames <- as.character(seqnames(ccres))  # Convert to character for faster paste0
        ccre_starts <- start(ccres)
        ccre_ends <- end(ccres)
        
        # Split hits by peak for efficient processing
        # This is much faster than multiple tapply calls
        peak_groups <- split(ccre_hits, peak_hits)
        peak_group_indices <- as.integer(names(peak_groups))
        
        # Count overlaps per peak (vectorized)
        mcols(peaks)$cCRE_count[peak_group_indices] <- 
            as.integer(table(factor(peak_hits, levels = peak_group_indices)))
        
        # Collect coordinates in standard format (chr1:100-1000)
        mcols(peaks)$cCRE_coordinates[peak_group_indices] <- vapply(peak_groups, function(hits) {
            coords <- paste0(ccre_seqnames[hits], ":", ccre_starts[hits], "-", ccre_ends[hits])
            paste(unique(coords), collapse = ",")
        }, character(1), USE.NAMES = FALSE)
        
        # Collect types (if available in metadata) - optimized
        if ("ccre_type" %in% colnames(ccre_mcols)) {
            type_values <- as.character(ccre_mcols$ccre_type)
            types_list <- vapply(peak_groups, function(hits) {
                paste(unique(type_values[hits]), collapse = ",")
            }, character(1), USE.NAMES = FALSE)
            mcols(peaks)$cCRE_types[peak_group_indices] <- types_list
        }
        
        # Collect accessions (if available) - optimized
        # Use id column (name is the same as cCRE_coordinates, so not needed)
        if ("id" %in% colnames(ccre_mcols)) {
            id_values <- as.character(ccre_mcols$id)
            accessions_list <- vapply(peak_groups, function(hits) {
                paste(unique(id_values[hits]), collapse = ",")
            }, character(1), USE.NAMES = FALSE)
            mcols(peaks)$cCRE_accession[peak_group_indices] <- accessions_list
        }
    }
    
    message("Annotation complete. ", sum(mcols(peaks)$cCRE_count > 0), 
           " out of ", length(peaks), " peaks overlap with cCREs.")
    
    return(peaks)
}

#' List available ENCODE cCRE tissue types
#' 
#' @description 
#' Lists available tissue types for ENCODE cCRE data. This helps users 
#' discover valid values for the \code{tissue_types} parameter in 
#' \code{\link{annotatePeaksWithcCRE}}.
#' 
#' @param species Character string specifying the species. Must be one of:
#'        \itemize{
#'          \item \code{"Homo sapiens"}: Human cCREs
#'          \item \code{"Mus musculus"}: Mouse cCREs
#'        }
#'        If \code{NULL} (default), uses the global species option set by
#'        \code{\link{setChIPpeakAnnoSpecies}}. If no global option is set,
#'        an error is raised prompting the user to specify species.
#' 
#' @return A data.frame with columns:
#'        \itemize{
#'          \item \code{tissue_type}: Character, tissue type identifier
#'          \item \code{description}: Character, description of the tissue type
#'          \item \code{biosample_term}: Character, ENCODE biosample term (if available)
#'          \item \code{file_accession}: Character, ENCODE file accession (if available)
#'        }
#' 
#' @details
#' This function queries ENCODE/SCREEN to discover available tissue-specific 
#' cCRE datasets. The returned data.frame can be used to select appropriate 
#' tissue types for \code{annotatePeaksWithcCRE}.
#' 
#' @note
#' \itemize{
#'   \item Requires internet connection to query ENCODE API
#'   \item Results may vary as ENCODE data is updated
#'   \item If API query fails, returns a basic data.frame with common tissue types
#' }
#' 
#' @seealso
#' \code{\link{annotatePeaksWithcCRE}} for using tissue types in annotation
#' 
#' @author Haibo Liu
#' @importFrom jsonlite fromJSON
#' @keywords misc
#' @export
#' @examples
#' 
#' # List available human cCRE tissue types
#' \dontrun{
#' available <- listAvailablecCREs("Homo sapiens")
#' head(available)
#' 
#' # Use a specific tissue type
#' data("myPeakList")
#' # original genome assembly is hg18, for now we set the global genome 
#' # assembly to hg38 to save time.
#' setChIPpeakAnnoGlobals(genome = "hg38", species = "Homo sapiens")
#' annotated <- annotatePeaksWithcCRE(myPeakList,
#'                                     tissue_types = available$tissue_type[1])
#' }
listAvailablecCREs <- function(species = NULL) {
    # Use global species option if species not provided
    if (is.null(species)) {
        species <- getChIPpeakAnnoSpecies()
        if (is.null(species)) {
            stop("'species' must be specified. Options: 'Homo sapiens' or 'Mus musculus'. ",
                 "You can set a global default using setChIPpeakAnnoSpecies().",
                 call. = FALSE)
        }
    }
    
    # Validate species
    species <- match.arg(species, c("Homo sapiens", "Mus musculus"))
    
    # Load SCREEN registry from RData file
    if (species == "Homo sapiens") {
        registry_file <- system.file("data", "human_cCRE_by_cell_tissue.rda", package = "ChIPpeakAnno")
    } else {
        registry_file <- system.file("data", "mouse_cCRE_by_cell_tissue.rda", package = "ChIPpeakAnno")
    }
    
    if (file.exists(registry_file)) {
        # Load the registry data
        env <- new.env()
        load(registry_file, envir = env)
        
        # Get the appropriate registry based on species
        if (species == "Homo sapiens") {
            registry <- env$human_cCRE_by_cell_tissue
        } else {
            registry <- env$mouse_cCRE_by_cell_tissue
        }
        
        if (!is.null(registry) && nrow(registry) > 0) {
            # Extract relevant columns
            # Column names: "Tissue/Biosample", "cCREs (.bed)", "Biosample", "Organ/Tissue", etc.
            tissue_info <- data.frame(
                tissue_type = registry[["Tissue/Biosample"]],
                description = if ("Biosample" %in% colnames(registry)) {
                    # Use Biosample column if available, otherwise use Tissue/Biosample
                    ifelse(!is.na(registry[["Biosample"]]) & registry[["Biosample"]] != "", 
                           registry[["Biosample"]], 
                           registry[["Tissue/Biosample"]])
                } else {
                    registry[["Tissue/Biosample"]]
                },
                download_url = if ("cCREs (.bed)" %in% colnames(registry)) {
                    registry[["cCREs (.bed)"]]
                } else {
                    rep(NA_character_, nrow(registry))
                },
                organ_tissue = if ("Organ/Tissue" %in% colnames(registry)) {
                    registry[["Organ/Tissue"]]
                } else {
                    rep(NA_character_, nrow(registry))
                },
                sample_type = if ("Sample Type" %in% colnames(registry)) {
                    registry[["Sample Type"]]
                } else {
                    rep(NA_character_, nrow(registry))
                },
                stringsAsFactors = FALSE
            )
            
            # Remove rows with empty tissue_type or download_url
            tissue_info <- tissue_info[!is.na(tissue_info$tissue_type) & 
                                       tissue_info$tissue_type != "" &
                                       !is.na(tissue_info$download_url) & 
                                       tissue_info$download_url != "", ]
            
            # Remove duplicates and sort
            tissue_info <- unique(tissue_info)
            tissue_info <- tissue_info[order(tissue_info$tissue_type), ]
            
            return(tissue_info)
        } else {
            stop("No cCRE data found for ", species, call. = FALSE)
        }
    } else {
        stop("SCREEN registry RData file not found: ", registry_file, call. = FALSE)
    }
}

# Internal helper functions

#' Download a cCRE file from URL to cache
#' @param url Character string, URL to download from
#' @param cache_file Character string, path to cache file
#' @return GRanges object with cCRE data, or empty GRanges on failure
#' @keywords internal
downloadcCREFromURL <- function(url, cache_file) {
    # Check cache first
    if (file.exists(cache_file)) {
        message("Using cached cCRE file: ", cache_file)
        return(loadcCREsFromFile(cache_file))
    }
    
    # Download from URL
    message("Downloading cCRE data from SCREEN...")
    message("URL: ", url)
    
    tryCatch({
        # Use binary mode for compressed files, text mode for plain text
        download_mode <- if (grepl("\\.gz$", url)) "wb" else "w"
        download.file(url, cache_file, mode = download_mode, quiet = FALSE)
        if (file.exists(cache_file) && file.info(cache_file)$size > 0) {
            return(loadcCREsFromFile(cache_file))
        } else {
            stop("Downloaded file is empty or missing")
        }
    }, error = function(e) {
        warning("Failed to download cCRE data: ", e$message, 
               ". Please download manually from https://screen.wenglab.org/downloads ",
               "and use the cCRE_file parameter.", call. = FALSE)
        return(GRanges())
    })
}

#' Load cCRE data from a file or URL
#' 
#' @description
#' Loads cCRE data from a file or URL. Handles format variations between
#' all-in-one cCRE files and tissue-specific files, as SCREEN files do not
#' strictly follow standard BED format.
#' 
#' @param file_path Character string, path to cCRE file or URL
#' @return GRanges object with cCRE data
#' @keywords internal
loadcCREsFromFile <- function(file_path) {
    if (!file.exists(file_path) && !grepl("^https?://", file_path)) {
        stop("File not found: ", file_path, call. = FALSE)
    }
    
    # Use custom parser directly for cCRE files
    # rtracklayer::import doesn't handle the non-standard cCRE formats properly
    # (6-column all-in-one format and 11-column tissue-specific format)
    # The custom parser ensures all columns are correctly parsed and preserved
    tryCatch({
            # Read file (handle compressed files)
            if (grepl("\\.gz$", file_path)) {
                con <- gzfile(file_path, "rt")
            } else {
                con <- file(file_path, "rt")
            }
            lines <- readLines(con)
            close(con)
            
            # Filter out header/comment lines and Low-DNase lines in tissue-specific files
            data_lines <- lines[!grepl("^track|^browser|^#|Low-DNase", lines, ignore.case = TRUE)]
            data_lines <- data_lines[data_lines != ""]  # Remove empty lines
            
            if (length(data_lines) == 0) {
                stop("No data lines found in cCRE file", call. = FALSE)
            }
            
            # Parse tab-separated columns and convert to data frame for efficiency
            # This is much faster than multiple vapply calls on the same list
            fields_list <- strsplit(data_lines, "\t", fixed = TRUE)
            
            # Check column count
            n_cols <- unique(lengths(fields_list))
            if (length(n_cols) > 1) {
                stop("Inconsistent column count in cCRE file. Found: ", 
                     paste(n_cols, collapse = ", "), " columns", call. = FALSE)
            }
            
            if (n_cols[1] < 3) {
                stop("cCRE file must have at least 3 columns (chr, start, end). ",
                     "Found: ", n_cols[1], " columns", call. = FALSE)
            }
            
            # Convert to data frame - optimized for large files
            # Pad shorter rows with NA to ensure consistent column count
            max_cols <- max(n_cols)
            # Use vapply to create matrix efficiently (faster than do.call(rbind) for large data)
            fields_matrix <- vapply(fields_list, function(x) {
                length(x) <- max_cols
                x
            }, character(max_cols))
            # Transpose to get rows as expected
            fields_matrix <- t(fields_matrix)
            df <- as.data.frame(fields_matrix, stringsAsFactors = FALSE)
            colnames(df) <- paste0("V", seq_len(ncol(df)))
            
            # Extract required columns (chr, start, end)
            chr <- df$V1
            start_pos <- as.integer(df$V2)
            end_pos <- as.integer(df$V3)
            
            # Determine format and extract strand
            # Format 1: 6 columns (all-in-one): chr, start, end, ID1, ID2, cCRE_type
            # Format 2: 11 columns (tissue-specific): chr, start, end, name, score, strand, 
            #           thickStart, thickEnd, itemRgb, cCRE_type, classification
            if (n_cols[1] == 6) {
                # All-in-one format: no strand column, use "*"
                strand_val <- "*"
            } else if (n_cols[1] >= 6) {
                # Tissue-specific format: strand is column 6
                strand_val <- ifelse(!is.na(df$V6) & df$V6 != ".", df$V6, "*")
            } else {
                strand_val <- "*"
            }
            
            # Create GRanges object
            ccres <- GRanges(
                seqnames = chr,
                ranges = IRanges(start = start_pos + 1L, end = end_pos),  # BED is 0-based, convert to 1-based
                strand = strand_val
            )
            
            # Add metadata columns based on format
            if (n_cols[1] == 6) {
                # All-in-one format: chr, start, end, ID1, ID2, cCRE_type
                mcols(ccres)$name <- paste0(df$V1, ":", df$V2, "_", df$V3)
                mcols(ccres)$id <- df$V5
                mcols(ccres)$ccre_type <- df$V6
            } else if (n_cols[1] == 11) {
                # Tissue-specific format: chr, start, end, name, score, strand, 
                # thickStart, thickEnd, itemRgb, cCRE_type, classification
                mcols(ccres)$name <- paste0(df$V1, ":", df$V2, "_", df$V3)
                mcols(ccres)$id <- df$V4
                mcols(ccres)$ccre_type <- df$V10
            } else {
                stop("Unsupported cCRE file format. Expected 6 or 11 columns. Found: ", 
                     n_cols[1], " columns", call. = FALSE)
            }
            
            message("Loaded ", length(ccres), " cCRE regions from file (using custom parser, ", 
                   n_cols[1], " columns)")
            return(ccres)
    }, error = function(e) {
            stop("Failed to load cCRE file: ", conditionMessage(e), 
                 ". File: ", file_path, call. = FALSE)
    })
}


#' Download all cCREs for a species
#' @param species Character string, species name
#' @param cache_dir Character string, cache directory
#' @return GRanges object with cCRE data
#' @keywords internal
downloadAllcCREs <- function(species, cache_dir) {
    # Load SCREEN registry to get download URL for "any"
    if (species == "Homo sapiens") {
        registry_file <- system.file("data", "human_cCRE_by_cell_tissue.rda", package = "ChIPpeakAnno")
        cache_file <- file.path(cache_dir, "SCREEN_AllHuman_cCREs.bed")
    } else {
        registry_file <- system.file("data", "mouse_cCRE_by_cell_tissue.rda", package = "ChIPpeakAnno")
        cache_file <- file.path(cache_dir, "SCREEN_AllMouse_cCREs.bed")
    }
    
    if (file.exists(registry_file)) {
        env <- new.env()
        load(registry_file, envir = env)
        
        if (species == "Homo sapiens") {
            registry <- env$human_cCRE_by_cell_tissue
        } else {
            registry <- env$mouse_cCRE_by_cell_tissue
        }
        
        # Find the "any" row which contains the URL for all cCREs
        any_row <- registry[registry[["Tissue/Biosample"]] == "any", ]
        if (nrow(any_row) > 0 && "cCREs (.bed)" %in% colnames(any_row)) {
            url <- any_row[["cCREs (.bed)"]][1]
            if (!is.na(url) && url != "") {
                result <- downloadcCREFromURL(url, cache_file)
                if (length(result) > 0) {
                    return(result)
                }
            }
        }
    }
    
    # Fallback: try direct SCREEN URLs if registry not available
    base_url <- "https://downloads.wenglab.org/Registry-V4"
    if (species == "Homo sapiens") {
        url <- paste0(base_url, "/downloads/human/GRCh38-cCREs.bed")
    } else {
        url <- paste0(base_url, "/downloads/mouse/mm10-cCREs.bed")
    }
    
    return(downloadcCREFromURL(url, cache_file))
}

#' Download tissue-specific cCREs
#' @param species Character string, species name
#' @param tissue_types Character vector, tissue type names
#' @param cache_dir Character string, cache directory
#' @return GRanges object with cCRE data
#' @keywords internal
downloadTissueSpecificcCREs <- function(species, tissue_types, cache_dir) {
    # Load SCREEN registry to get download URLs for specific tissues
    if (species == "Homo sapiens") {
        registry_file <- system.file("data", "human_cCRE_by_cell_tissue.rda", package = "ChIPpeakAnno")
    } else {
        registry_file <- system.file("data", "mouse_cCRE_by_cell_tissue.rda", package = "ChIPpeakAnno")
    }
    
    if (!file.exists(registry_file)) {
        warning("SCREEN registry RData file not found. Cannot download tissue-specific cCREs.",
               " Please use tissue_types='any' or provide cCRE_file parameter.",
               call. = FALSE)
        return(GRanges())
    }
    
    env <- new.env()
    load(registry_file, envir = env)
    
    if (species == "Homo sapiens") {
        registry <- env$human_cCRE_by_cell_tissue
    } else {
        registry <- env$mouse_cCRE_by_cell_tissue
    }
    
    if (is.null(registry) || nrow(registry) == 0) {
        warning("Registry is empty. Cannot download tissue-specific cCREs.",
               call. = FALSE)
        return(GRanges())
    }
    
    # Find matching rows for each tissue type
    all_ccres <- GRangesList()
    
    for (tissue in tissue_types) {
        # Try exact match first
        matches <- registry[registry[["Tissue/Biosample"]] == tissue, ]
        
        # If no exact match, try case-insensitive partial match
        if (nrow(matches) == 0) {
            matches <- registry[grepl(tissue, registry[["Tissue/Biosample"]], ignore.case = TRUE), ]
        }
        
        # Also try matching in Biosample column
        if (nrow(matches) == 0 && "Biosample" %in% colnames(registry)) {
            matches <- registry[registry[["Biosample"]] == tissue, ]
            if (nrow(matches) == 0) {
                matches <- registry[grepl(tissue, registry[["Biosample"]], ignore.case = TRUE), ]
            }
        }
        
        if (nrow(matches) > 0 && "cCREs (.bed)" %in% colnames(matches)) {
            # Get the download URL
            urls <- matches[["cCREs (.bed)"]]
            urls <- urls[!is.na(urls) & urls != ""]
            
            if (length(urls) > 0) {
                # Use the first matching URL
                url <- urls[1]
                
                # Create cache file name based on tissue
                safe_tissue_name <- gsub("[^A-Za-z0-9]", "_", tissue)
                cache_file <- file.path(cache_dir, paste0("SCREEN_", safe_tissue_name, "_cCREs.bed"))
                
                # Check cache first
                if (file.exists(cache_file)) {
                    message("Using cached cCRE file for ", tissue, ": ", cache_file)
                    ccres <- loadcCREsFromFile(cache_file)
                } else {
                    # Download from SCREEN
                    message("Downloading cCRE data for ", tissue, " from SCREEN...")
                    message("URL: ", url)
                    
                    tryCatch({
                        # Use binary mode for compressed files, text mode for plain text
                        download_mode <- if (grepl("\\.gz$", url)) "wb" else "w"
                        download.file(url, cache_file, mode = download_mode, quiet = FALSE)
                        if (file.exists(cache_file) && file.info(cache_file)$size > 0) {
                            ccres <- loadcCREsFromFile(cache_file)
                        } else {
                            warning("Downloaded file is empty for tissue: ", tissue, call. = FALSE)
                            next
                        }
                    }, error = function(e) {
                        warning("Failed to download cCRE data for tissue '", tissue, "': ", 
                               e$message, call. = FALSE)
                        next
                    })
                }
                
                if (exists("ccres") && length(ccres) > 0) {
                    all_ccres[[tissue]] <- ccres
                }
            } else {
                warning("No download URL found for tissue: ", tissue, call. = FALSE)
            }
        } else {
            warning("Tissue '", tissue, "' not found in SCREEN registry. ",
                   "Use listAvailablecCREs() to see available tissues.", call. = FALSE)
        }
    }
    
    # Combine all tissue-specific cCREs
    if (length(all_ccres) > 0) {
        combined <- unlist(all_ccres)
        # Remove duplicates if multiple tissues were requested
        combined <- unique(combined)
        message("Loaded ", length(combined), " cCRE regions from ", length(all_ccres), " tissue(s)")
        return(combined)
    } else {
        warning("No cCRE data was successfully downloaded for the requested tissues.",
               call. = FALSE)
        return(GRanges())
    }
}
