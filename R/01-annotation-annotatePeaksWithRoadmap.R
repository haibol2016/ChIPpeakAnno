#' Annotate peaks with Roadmap Epigenomics data
#' 
#' @description 
#' Annotates genomic peaks by finding overlaps with Roadmap Epigenomics data,
#' including DNase hypersensitive sites (promoters and enhancers) and chromatin
#' state segmentations. This function can download tissue/sample-specific data
#' from Roadmap Epigenomics based on epigenome IDs (EIDs) or tissue types.
#' 
#' This function uses data from the NIH Roadmap Epigenomics Mapping Consortium,
#' which provides genome-wide maps of histone modifications, chromatin 
#' accessibility, and 5-mark chromatin states across 127 human epigenomes.
#' 
#' @param peaks A \code{\link[GenomicRanges]{GRanges}} object containing 
#'        peaks to be annotated. Peaks can be pre-annotated by 
#'        \code{\link{annotatePeakInBatch}} or other annotation functions.
#' @param data_type Character string specifying the type of Roadmap data to use.
#'        Must be one of:
#'        \itemize{
#'          \item \code{"promoter"} (default): DNase hypersensitive sites in 
#'                promoter regions
#'          \item \code{"enhancer"}: DNase hypersensitive sites in enhancer regions
#'          \item \code{"chromhmm"}: ChromHMM 15-state chromatin state 
#'                segmentations
#'        }
#' @param EIDs Character vector of epigenome IDs (e.g., "E001", "E017") to use.
#'        If \code{NULL} (default), uses all available epigenomes for the 
#'        selected data type. Use \code{\link{listAvailableRoadmapEpigenomes}} 
#'        to see available EIDs.
#' @param tissue_types Character vector specifying tissue types to filter by.
#'        If \code{NULL} (default), no tissue filtering is applied. Use 
#'        \code{\link{listAvailableRoadmapTissues}} to see available tissue 
#'        types. This parameter filters EIDs by their tissue type before 
#'        downloading data.
#' @param roadmap_file Optional character string. If provided, uses this file 
#'        path or URL instead of downloading from Roadmap. Can be a local file 
#'        path or a URL to a BED file. This takes precedence over \code{EIDs} 
#'        and \code{tissue_types}.
#' @param cache_dir Optional character string specifying a directory for 
#'        caching downloaded Roadmap files. If \code{NULL} (default), uses 
#'        \code{tempdir()}. Files are cached to avoid re-downloading.
#' @param combine_regions Logical. If \code{TRUE} (default), combines regions 
#'        from multiple epigenomes into a single GRanges object before finding 
#'        overlaps. If \code{FALSE}, processes each epigenome separately and 
#'        reports overlaps per epigenome.
#' @param species Character string specifying the species. Currently only 
#'        \code{"Homo sapiens"} is supported (default). Roadmap Epigenomics 
#'        data is only available for human.
#' @param ... Additional arguments passed to \code{\link[IRanges]{findOverlaps}}, 
#'        such as \code{maxgap}, \code{minoverlap}, \code{type}, etc.
#' 
#' @return A \code{\link[GenomicRanges]{GRanges}} object containing the input 
#'        peaks with additional metadata columns:
#'        \itemize{
#'          \item \code{roadmap_count}: Integer, number of overlapping Roadmap 
#'                regions
#'          \item \code{roadmap_coordinates}: Character, comma-separated list of 
#'                Roadmap region coordinates in the format "chr:start-end" 
#'                (e.g., "chr1:10000-20000,chr2:5000-8000")
#'          \item \code{roadmap_EIDs}: Character, comma-separated list of 
#'                epigenome IDs (EIDs) for overlapping regions
#'          \item \code{roadmap_tissues}: Character, comma-separated list of 
#'                tissue types for overlapping regions
#'          \item \code{roadmap_epigenome_names}: Character, comma-separated 
#'                list of standardized epigenome names (e.g., "ES-I3 Cells", 
#'                "H1 Cells", "IMR90 fetal lung fibroblasts Cell Line") for 
#'                overlapping regions. These are descriptive names of the cell 
#'                types or tissues from which the Roadmap data was derived.
#'          \item \code{roadmap_region_type}: Character, indicates the type of 
#'                Roadmap region. For \code{data_type = "promoter"}, this column 
#'                contains "promoter". For \code{data_type = "enhancer"}, this 
#'                column contains "enhancer". For \code{data_type = "chromhmm"}, 
#'                this column is not created (see \code{roadmap_chromhmm_states} 
#'                instead).
#'          \item \code{roadmap_chromhmm_states}: Character, comma-separated 
#'                list of ChromHMM states (e.g., "1_TssA", "7_Enh"). 
#'                \strong{Only present when \code{data_type = "chromhmm"}}. 
#'                ChromHMM states are not available in the promoter or enhancer 
#'                data files because those files contain DNase hypersensitive 
#'                sites classified by accessibility patterns, while ChromHMM 
#'                states are derived from histone modification patterns (H3K4me3, 
#'                H3K27ac, H3K4me1, H3K9me3, H3K36me3) in separate segmentation 
#'                files.
#'        }
#'        Peaks without overlaps will have \code{roadmap_count = 0} and 
#'        \code{NA} for other Roadmap-related columns.
#' 
#' @details
#' 
#' \strong{Roadmap Epigenomics Data:}
#' The Roadmap Epigenomics Project provides:
#' \itemize{
#'   \item \strong{DNase Promoters}: DNase hypersensitive sites in promoter 
#'         regions (111 epigenomes)
#'   \item \strong{DNase Enhancers}: DNase hypersensitive sites in enhancer 
#'         regions (111 epigenomes)
#'   \item \strong{ChromHMM States}: 15-state chromatin state segmentations 
#'         (127 epigenomes, including 16 ENCODE)
#' }
#' 
#' \strong{Species Restriction:}
#' \itemize{
#'   \item This function is \strong{only for human (Homo sapiens) data}
#'   \item The \code{species} parameter must be \code{"Homo sapiens"} (default)
#'   \item Roadmap Epigenomics data is only available for human
#' }
#' 
#' \strong{Genome Assembly:}
#' \itemize{
#'   \item All Roadmap data uses \strong{hg19 (GRCh37)} assembly
#'   \item If your peaks are in a different human assembly (e.g., hg38), the 
#'         function will automatically lift over the Roadmap data to match your 
#'         peaks' genome assembly using \code{rtracklayer::liftOver()}
#'   \item Chain files are automatically downloaded from UCSC if needed
#' }
#' 
#' \strong{Data Sources:}
#' \itemize{
#'   \item URLs are obtained from \code{\link{roadmap_bed_gz_links_metadata}}
#'   \item Original data: \url{https://egg2.wustl.edu/roadmap/web_portal/}
#' }
#' 
#' \strong{Download Behavior:}
#' \itemize{
#'   \item If \code{roadmap_file} is provided, it takes precedence and no 
#'         download occurs
#'   \item If \code{EIDs} is specified, downloads data for those specific 
#'         epigenomes
#'   \item If \code{tissue_types} is specified, filters EIDs by tissue type 
#'         first, then downloads
#'   \item If both are \code{NULL}, uses all available epigenomes for the 
#'         selected data type
#'   \item Downloaded files are cached in \code{cache_dir} to avoid 
#'         re-downloading
#' }
#' 
#' \strong{Overlap Detection:}
#' Uses \code{findOverlaps(peaks, roadmap_regions)} to identify overlaps. By 
#' default, any overlap is considered. Use \code{maxgap} parameter via 
#' \code{...} to allow gaps between peaks and Roadmap regions.
#' 
#' @note
#' \itemize{
#'   \item This function is \strong{only for human data} (Homo sapiens)
#'   \item Roadmap data is on hg19 (GRCh37), but the function automatically 
#'         lifts over to hg38 if your peaks are on hg38
#'   \item Large datasets (especially when using many epigenomes) may require 
#'         significant memory and download time
#'   \item The function preserves all existing metadata columns from the input 
#'         peaks
#'   \item Downloaded files are cached to improve performance in subsequent 
#'         calls
#'   \item Roadmap server may be slow; consider downloading files manually for 
#'         large-scale analyses
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{listAvailableRoadmapEpigenomes}} to discover available 
#'         EIDs and their metadata
#'   \item \code{\link{listAvailableRoadmapTissues}} to see available tissue types
#'   \item \code{\link{roadmap_bed_gz_links_metadata}} for the URL registry
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
#' @importFrom rtracklayer import import.chain liftOver
#' @examples
#' 
#' # Example 1: Annotate with promoter regions from all epigenomes
#' \dontrun{
#' data("myPeakList")
#' # Peaks can be on hg19 or hg38 - function will automatically lift over if needed
#' annotated <- annotatePeaksWithRoadmap(myPeakList, 
#'                                       data_type = "promoter")
#' # Check how many peaks overlap with Roadmap regions
#' table(annotated$roadmap_count > 0)
#' # View coordinates of overlapping Roadmap regions
#' head(annotated$roadmap_coordinates[annotated$roadmap_count > 0])
#' # View EIDs for overlapping peaks
#' head(annotated$roadmap_EIDs[annotated$roadmap_count > 0])
#' }
#' 
#' # Example 2: Annotate with enhancer regions from specific EIDs
#' \dontrun{
#' annotated <- annotatePeaksWithRoadmap(myPeakList,
#'                                       data_type = "enhancer",
#'                                       EIDs = c("E001", "E017", "E003"))
#' }
#' 
#' # Example 3: Annotate with ChromHMM states from brain tissues
#' \dontrun{
#' annotated <- annotatePeaksWithRoadmap(myPeakList,
#'                                       data_type = "chromhmm",
#'                                       tissue_types = "BRAIN")
#' }
#' 
#' # Example 4: Use a pre-downloaded Roadmap file
#' \dontrun{
#' annotated <- annotatePeaksWithRoadmap(myPeakList,
#'                                       roadmap_file = "/path/to/roadmap.bed.gz")
#' }
annotatePeaksWithRoadmap <- function(peaks,
                                    data_type = c("promoter", "enhancer", "chromhmm"),
                                    EIDs = NULL,
                                    tissue_types = NULL,
                                    roadmap_file = NULL,
                                    cache_dir = NULL,
                                    combine_regions = TRUE,
                                    species = "Homo sapiens",
                                    ...) {
    # Validate inputs
    if (!inherits(peaks, "GRanges")) {
        stop("'peaks' must be a GRanges object", call. = FALSE)
    }
    
    data_type <- match.arg(data_type)
    
    # Validate species - only human is supported
    if (species != "Homo sapiens") {
        stop("Roadmap Epigenomics data is only available for human (Homo sapiens). ",
             "Species '", species, "' is not supported.",
             call. = FALSE)
    }
    
    # Detect target genome assembly using hybrid approach (global option + peak detection)
    # Note: Roadmap is human-only, so species is always "Homo sapiens"
    target_genome <- getTargetGenome(peaks, species = "Homo sapiens")
    
    # Roadmap data is always hg19, so we need to lift over if target is different
    need_liftover <- !is.null(target_genome) && target_genome != "hg19"
    
    if (need_liftover) {
        message("Detected peaks on ", target_genome, " assembly. ",
               "Roadmap data is on hg19. Will lift over Roadmap data to ", 
               target_genome, ".")
    } else if (is.null(target_genome)) {
        global_genome <- getChIPpeakAnnoGenome()
        if (!is.null(global_genome)) {
            warning("Could not detect genome assembly from peaks. ",
                   "Using global genome option: ", global_genome, 
                   ". If this is incorrect, set genome(peaks) explicitly or ",
                   "use setChIPpeakAnnoGenome() to change the global option.",
                   call. = FALSE)
            # Use global genome for liftOver decision
            target_genome <- global_genome
            need_liftover <- target_genome != "hg19"
        } else {
            stop("Could not detect genome assembly from peaks and no global genome is set. ",
                 "Please either:\n",
                 "  1. Set genome metadata on peaks: genome(peaks) <- 'hg38' (or your assembly)\n",
                 "  2. Set a global genome option: setChIPpeakAnnoGenome('hg38')\n",
                 "  3. Use setChIPpeakAnnoGlobals() to set all global options at once",
                 call. = FALSE)
        }
    }
    
    if (is.null(cache_dir)) {
        cache_dir <- tempdir()
    }
    
    # Load or download Roadmap data
    if (!is.null(roadmap_file)) {
        # Use provided file
        message("Loading Roadmap data from file: ", roadmap_file)
        roadmap_regions <- loadRoadmapFromFile(roadmap_file)
    } else {
        # Download based on data type, EIDs, and tissue types
        message("Downloading Roadmap ", data_type, " data...")
        roadmap_regions <- downloadRoadmapData(data_type, EIDs, tissue_types, 
                                               cache_dir, combine_regions)
    }
    
    if (length(roadmap_regions) == 0) {
        warning("No Roadmap data loaded. Returning peaks without Roadmap annotation.",
                call. = FALSE)
        # Initialize empty annotation columns and return
        mcols(peaks)$roadmap_count <- 0L
        mcols(peaks)$roadmap_coordinates <- NA_character_
        mcols(peaks)$roadmap_EIDs <- NA_character_
        mcols(peaks)$roadmap_tissues <- NA_character_
        mcols(peaks)$roadmap_epigenome_names <- NA_character_
        if (data_type %in% c("promoter", "enhancer")) {
            mcols(peaks)$roadmap_region_type <- NA_character_
        }
        if (data_type == "chromhmm") {
            mcols(peaks)$roadmap_chromhmm_states <- NA_character_
        }
        return(peaks)
    }
    
    # Handle GRangesList or list if combine_regions = FALSE
    # seqlevelsStyle doesn't work on GRangesList or list, so we need to convert to GRanges
    # This must happen BEFORE any seqlevelsStyle calls
    if (inherits(roadmap_regions, "GRangesList")) {
        # For overlap finding, we need a single GRanges object
        # Even if combine_regions = FALSE, we combine here for seqlevelsStyle and liftOver
        message("Converting GRangesList to GRanges for processing...")
        roadmap_regions <- unlist(roadmap_regions)
    } else if (is.list(roadmap_regions) && !inherits(roadmap_regions, "GRanges")) {
        # Handle plain list (shouldn't happen, but be safe)
        # Try to convert if it's a list of GRanges
        # Filter out NULL entries first
        roadmap_regions <- roadmap_regions[!sapply(roadmap_regions, is.null)]
        if (length(roadmap_regions) > 0 && all(sapply(roadmap_regions, inherits, "GRanges"))) {
            message("Converting list of GRanges to single GRanges object...")
            roadmap_regions <- do.call(c, roadmap_regions)
        } else if (length(roadmap_regions) == 0) {
            stop("roadmap_regions is an empty list. No Roadmap data available.", call. = FALSE)
        } else {
            # Debug: show what we got
            classes <- sapply(roadmap_regions, class)
            stop("roadmap_regions is a list but not a GRangesList or list of GRanges. ",
                 "Expected GRanges, GRangesList, or list of GRanges. ",
                 "Got list with element classes: ", paste(unique(classes), collapse = ", "), 
                 call. = FALSE)
        }
    }
    
    # Ensure roadmap_regions is a GRanges object before proceeding
    # This is critical - all downstream operations require a GRanges object
    if (!inherits(roadmap_regions, "GRanges")) {
        # Try one more time to convert if it's a list
        if (is.list(roadmap_regions)) {
            # Try to extract GRanges from list
            if (length(roadmap_regions) == 1 && inherits(roadmap_regions[[1]], "GRanges")) {
                message("Extracting GRanges from single-element list...")
                roadmap_regions <- roadmap_regions[[1]]
            } else if (all(sapply(roadmap_regions, inherits, "GRanges"))) {
                message("Converting list of GRanges to single GRanges...")
                roadmap_regions <- do.call(c, roadmap_regions)
            }
        }
        
        # Final check
        if (!inherits(roadmap_regions, "GRanges")) {
            stop("roadmap_regions must be a GRanges object after processing. ",
                 "Got: ", paste(class(roadmap_regions), collapse = ", "), 
                 ". Length: ", length(roadmap_regions), 
                 ". This indicates a problem in downloadRoadmapData or loadRoadmapFromFile.",
                 call. = FALSE)
        }
    }
    
    # Lift over Roadmap data if needed
    if (need_liftover) {
        message("Lifting over Roadmap data from hg19 to ", target_genome, "...")
        chain <- getLiftOverChain("hg19", target_genome, cache_dir)
        if (!is.null(chain)) {
            # liftOver returns a GRangesList (one element per input range)
            # Some ranges may not map, resulting in empty GRanges
            roadmap_lifted <- liftOver(roadmap_regions, chain)
            # Convert to GRanges, keeping only ranges that successfully mapped
            roadmap_regions <- unlist(roadmap_lifted)
            if (length(roadmap_regions) == 0) {
                stop("No Roadmap regions could be lifted over to ", target_genome, 
                     ". This may indicate incompatible genome assemblies.", 
                     call. = FALSE)
            }
            message("Successfully lifted over ", length(roadmap_regions), 
                   " Roadmap regions to ", target_genome, ".")
        } else {
            warning("Could not obtain liftOver chain file. ",
                   "Proceeding with hg19 Roadmap data, which may not match your peaks.",
                   call. = FALSE)
        }
    }
    
    # Match seqlevels styles if needed
    # Double-check that roadmap_regions is a GRanges object before calling seqlevelsStyle
    if (!inherits(roadmap_regions, "GRanges")) {
        stop("Internal error: roadmap_regions is not a GRanges object. ",
             "Got: ", paste(class(roadmap_regions), collapse = ", "), 
             ". This should not happen.", call. = FALSE)
    }
    
    peaks_style <- tryCatch(seqlevelsStyle(peaks), error = function(e) NULL)
    roadmap_style <- tryCatch(seqlevelsStyle(roadmap_regions), error = function(e) NULL)
    
    if (is.null(peaks_style) || is.null(roadmap_style)) {
        warning("Could not determine seqlevels styles. Proceeding without matching.", 
                call. = FALSE)
    } else if (length(intersect(peaks_style, roadmap_style)) == 0) {
        message("Matching seqlevels styles between peaks and Roadmap data")
        seqlevelsStyle(roadmap_regions) <- peaks_style[1]
    }
    
    # Find overlaps
    message("Finding overlaps between peaks and Roadmap regions...")
    overlaps <- findOverlaps(peaks, roadmap_regions, ...)
    
    # Initialize metadata columns
    mcols(peaks)$roadmap_count <- 0L
    mcols(peaks)$roadmap_coordinates <- NA_character_
    mcols(peaks)$roadmap_EIDs <- NA_character_
    mcols(peaks)$roadmap_tissues <- NA_character_
    mcols(peaks)$roadmap_epigenome_names <- NA_character_
    if (data_type %in% c("promoter", "enhancer")) {
        mcols(peaks)$roadmap_region_type <- NA_character_
    }
    if (data_type == "chromhmm") {
        mcols(peaks)$roadmap_chromhmm_states <- NA_character_
    }
    
    if (length(overlaps) > 0) {
        # Get peak and Roadmap region indices from overlaps
        peak_hits <- queryHits(overlaps)
        roadmap_hits <- subjectHits(overlaps)
        
        # Mark peaks with overlaps
        peak_indices <- unique(peak_hits)
        
        # Count overlaps per peak (vectorized)
        mcols(peaks)$roadmap_count[peak_indices] <- 
            as.integer(table(factor(peak_hits, levels = peak_indices)))
        
        # Pre-extract metadata columns to avoid repeated indexing (performance optimization)
        roadmap_mcols <- mcols(roadmap_regions)
        roadmap_seqnames <- as.character(seqnames(roadmap_regions))
        roadmap_starts <- start(roadmap_regions)
        roadmap_ends <- end(roadmap_regions)
        
        # Split hits by peak for efficient processing
        # This is much faster than multiple tapply calls
        peak_groups <- split(roadmap_hits, peak_hits)
        peak_group_indices <- as.integer(names(peak_groups))
        
        # Collect Roadmap region coordinates (format: chr:start-end) - vectorized
        coordinates_list <- vapply(peak_groups, function(hits) {
            coords <- paste0(roadmap_seqnames[hits], ":",
                            roadmap_starts[hits], "-",
                            roadmap_ends[hits])
            paste(unique(coords), collapse = ",")
        }, character(1), USE.NAMES = FALSE)
        mcols(peaks)$roadmap_coordinates[peak_group_indices] <- coordinates_list
        
        # Collect EIDs (if available in metadata) - vectorized
        if ("EID" %in% colnames(roadmap_mcols)) {
            eid_values <- as.character(roadmap_mcols$EID)
            eids_list <- vapply(peak_groups, function(hits) {
                paste(unique(eid_values[hits]), collapse = ",")
            }, character(1), USE.NAMES = FALSE)
            mcols(peaks)$roadmap_EIDs[peak_group_indices] <- eids_list
        }
        
        # Collect tissue types (if available) - vectorized
        if ("anatomy" %in% colnames(roadmap_mcols)) {
            anatomy_values <- as.character(roadmap_mcols$anatomy)
            tissues_list <- vapply(peak_groups, function(hits) {
                paste(unique(anatomy_values[hits]), collapse = ",")
            }, character(1), USE.NAMES = FALSE)
            mcols(peaks)$roadmap_tissues[peak_group_indices] <- tissues_list
        }
        
        # Collect epigenome names (if available) - vectorized
        if ("standardized_epigenome_name" %in% colnames(roadmap_mcols)) {
            name_values <- as.character(roadmap_mcols$standardized_epigenome_name)
            names_list <- vapply(peak_groups, function(hits) {
                paste(unique(name_values[hits]), collapse = ",")
            }, character(1), USE.NAMES = FALSE)
            mcols(peaks)$roadmap_epigenome_names[peak_group_indices] <- names_list
        }
        
        # Set region type for promoter/enhancer data
        if (data_type %in% c("promoter", "enhancer")) {
            mcols(peaks)$roadmap_region_type[peak_indices] <- data_type
        }
        
        # Collect ChromHMM states (if available and data_type is chromhmm) - vectorized
        if (data_type == "chromhmm") {
            # ChromHMM states are in the 4th column (name field) of BED files
            # rtracklayer::import loads this into the 'name' metadata column
            if ("name" %in% colnames(roadmap_mcols)) {
                state_values <- as.character(roadmap_mcols$name)
                states_list <- vapply(peak_groups, function(hits) {
                    paste(unique(state_values[hits]), collapse = ",")
                }, character(1), USE.NAMES = FALSE)
                mcols(peaks)$roadmap_chromhmm_states[peak_group_indices] <- states_list
            }
        }
    }
    
    n_overlapping <- sum(mcols(peaks)$roadmap_count > 0)
    message("Annotation complete. ", n_overlapping, 
           " out of ", length(peaks), " peaks overlap with Roadmap regions.")
    
    return(peaks)
}

#' List available Roadmap Epigenomics epigenomes
#' 
#' @description 
#' Lists available epigenomes (EIDs) from Roadmap Epigenomics with their 
#' metadata. This helps users discover valid EIDs for the \code{EIDs} parameter 
#' in \code{\link{annotatePeaksWithRoadmap}}.
#' 
#' @param data_type Character string specifying the data type to filter by.
#'        Must be one of: \code{"promoter"}, \code{"enhancer"}, or 
#'        \code{"chromhmm"}. If \code{NULL} (default), returns all epigenomes.
#' @param tissue_types Character vector specifying tissue types to filter by.
#'        If \code{NULL} (default), no tissue filtering is applied.
#' 
#' @return A data.frame with columns from \code{\link{roadmap_metadata}}:
#'        \itemize{
#'          \item \code{epigenome_id_eid}: Character, epigenome ID (EID)
#'          \item \code{standardized_epigenome_name}: Character, full name
#'          \item \code{anatomy}: Character, tissue type
#'          \item \code{group}: Character, biological group
#'          \item \code{class}: Character, quality class
#'          \item And other metadata columns...
#'        }
#' 
#' @details
#' This function queries the \code{roadmap_bed_gz_links_metadata} dataset to 
#' discover available epigenomes. The returned data.frame can be used to select 
#' appropriate EIDs for \code{annotatePeaksWithRoadmap}.
#' 
#' @seealso
#' \code{\link{annotatePeaksWithRoadmap}} for using EIDs in annotation
#' \code{\link{roadmap_metadata}} for the full metadata
#' 
#' @author Haibo Liu
#' @keywords misc
#' @export
#' @examples
#' 
#' # List all available epigenomes
#' \dontrun{
#' available <- listAvailableRoadmapEpigenomes()
#' head(available)
#' 
#' # List epigenomes available for promoter data
#' promoter_epigenomes <- listAvailableRoadmapEpigenomes(data_type = "promoter")
#' 
#' # List brain epigenomes
#' brain_epigenomes <- listAvailableRoadmapEpigenomes(tissue_types = "BRAIN")
#' }
listAvailableRoadmapEpigenomes <- function(data_type = NULL,
                                          tissue_types = NULL) {
    # Load roadmap metadata
    env <- new.env()
    data("roadmap_metadata", package = "ChIPpeakAnno", envir = env)
    roadmap_metadata <- env$roadmap_metadata
    
    # Filter by data type if specified
    if (!is.null(data_type)) {
        data_type <- match.arg(data_type, c("promoter", "enhancer", "chromhmm"))
        
        # Load roadmap_bed_gz_links_metadata to get available EIDs
        data("roadmap_bed_gz_links_metadata", package = "ChIPpeakAnno", 
             envir = env)
        roadmap_bed_gz_links_metadata <- env$roadmap_bed_gz_links_metadata
        
        # Map data_type to list element name
        element_name <- switch(data_type,
            "promoter" = "DHS_promoter",
            "enhancer" = "DHS_enhancer",
            "chromhmm" = "ChromHMM_15_states_5_coreMarks"
        )
        
        if (element_name %in% names(roadmap_bed_gz_links_metadata)) {
            available_eids <- roadmap_bed_gz_links_metadata[[element_name]]$EID
            roadmap_metadata <- roadmap_metadata[
                roadmap_metadata$epigenome_id_eid %in% available_eids, 
            ]
        } else {
            stop("Data type '", data_type, "' not found in roadmap_bed_gz_links_metadata",
                 call. = FALSE)
        }
    }
    
    # Filter by tissue types if specified
    if (!is.null(tissue_types)) {
        roadmap_metadata <- roadmap_metadata[
            roadmap_metadata$anatomy %in% tissue_types, 
        ]
    }
    
    if (nrow(roadmap_metadata) == 0) {
        warning("No epigenomes found matching the criteria", call. = FALSE)
    }
    
    return(roadmap_metadata)
}

#' List available Roadmap Epigenomics tissue types
#' 
#' @description 
#' Lists available tissue types from Roadmap Epigenomics. This helps users 
#' discover valid values for the \code{tissue_types} parameter in 
#' \code{\link{annotatePeaksWithRoadmap}}.
#' 
#' @return A character vector of unique tissue type names (anatomy categories)
#' 
#' @details
#' This function extracts unique tissue types from \code{roadmap_metadata}.
#' 
#' @seealso
#' \code{\link{annotatePeaksWithRoadmap}} for using tissue types in annotation
#' \code{\link{roadmap_metadata}} for the full metadata
#' 
#' @author Haibo Liu
#' @keywords misc
#' @export
#' @examples
#' 
#' # List all available tissue types
#' \dontrun{
#' tissues <- listAvailableRoadmapTissues()
#' print(tissues)
#' 
#' # Use a specific tissue type
#' annotated <- annotatePeaksWithRoadmap(peaks,
#'                                       data_type = "promoter",
#'                                       tissue_types = "BRAIN")
#' }
listAvailableRoadmapTissues <- function() {
    # Load roadmap metadata
    env <- new.env()
    data("roadmap_metadata", package = "ChIPpeakAnno", envir = env)
    roadmap_metadata <- env$roadmap_metadata
    
    tissues <- unique(roadmap_metadata$anatomy)
    tissues <- tissues[!is.na(tissues)]
    tissues <- sort(tissues)
    
    return(tissues)
}

# Internal helper functions

#' Load Roadmap data from a file or URL
#' @param file_path Character string, path to BED file or URL
#' @return GRanges object with Roadmap data
#' @keywords internal
loadRoadmapFromFile <- function(file_path) {
    if (!file.exists(file_path) && !grepl("^https?://", file_path)) {
        stop("File not found: ", file_path, call. = FALSE)
    }
    
    tryCatch({
        # Use rtracklayer::import directly for BED files
        # This is more robust and handles various BED file variations
        roadmap_regions <- rtracklayer::import(file_path, format = "BED")
        
        # Validate that import returned a GRanges object
        if (!inherits(roadmap_regions, "GRanges")) {
            stop("rtracklayer::import did not return a GRanges object. Got: ", 
                 paste(class(roadmap_regions), collapse = ", "), call. = FALSE)
        }
        
        message("Loaded ", length(roadmap_regions), " Roadmap regions from file")
        return(roadmap_regions)
    }, error = function(e) {
        # Provide diagnostic information if import fails
        error_msg <- conditionMessage(e)
        
        if (file.exists(file_path)) {
            file_size <- file.info(file_path)$size
            if (file_size == 0) {
                stop("Roadmap file is empty: ", file_path, call. = FALSE)
            }
            
            # Try to read first few lines to diagnose the issue
            tryCatch({
                con <- gzfile(file_path, "rt")
                first_lines <- readLines(con, n = 5)
                close(con)
                
                if (length(first_lines) > 0) {
                    # Filter out comment/track lines
                    data_lines <- first_lines[!grepl("^track|^browser|^#", first_lines)]
                    if (length(data_lines) > 0) {
                        n_cols <- length(strsplit(data_lines[1], "\t")[[1]])
                        if (n_cols < 3) {
                            stop("Roadmap file appears to be malformed. ",
                                 "Expected at least 3 columns (chr, start, end) in BED format. ",
                                 "Found ", n_cols, " column(s) in first data line. ",
                                 "File: ", file_path, ". ",
                                 "Error: ", error_msg, call. = FALSE)
                        } else {
                            stop("Failed to load Roadmap file. ",
                                 "File: ", file_path, ". ",
                                 "First data line has ", n_cols, " column(s). ",
                                 "Error: ", error_msg, call. = FALSE)
                        }
                    } else {
                        stop("Roadmap file contains only header/track lines, no data. ",
                             "File: ", file_path, call. = FALSE)
                    }
                }
            }, error = function(e2) {
                # If we can't read the file for diagnosis, just report the error
                stop("Failed to load Roadmap file: ", error_msg, ". ",
                     "File: ", file_path, ". ",
                     "Unable to diagnose file format issue.",
                     call. = FALSE)
            })
        } else {
            stop("Failed to load Roadmap file: ", error_msg, ". ",
                 "File: ", file_path, call. = FALSE)
        }
    })
}

#' Download Roadmap data
#' @param data_type Character string, data type
#' @param EIDs Character vector, epigenome IDs
#' @param tissue_types Character vector, tissue types
#' @param cache_dir Character string, cache directory
#' @param combine_regions Logical, whether to combine regions (currently always TRUE)
#' @return GRanges object with Roadmap data (always combined)
#' @keywords internal
downloadRoadmapData <- function(data_type, EIDs, tissue_types, cache_dir, combine_regions) {
    # Load roadmap_bed_gz_links_metadata
    # Note: roadmap_bed_gz_links_metadata already contains merged metadata from roadmap_metadata
    env <- new.env()
    data("roadmap_bed_gz_links_metadata", package = "ChIPpeakAnno", envir = env)
    roadmap_bed_gz_links_metadata <- env$roadmap_bed_gz_links_metadata
    
    # Map data_type to list element name
    element_name <- switch(data_type,
        "promoter" = "DHS_promoter",
        "enhancer" = "DHS_enhancer",
        "chromhmm" = "ChromHMM_15_states_5_coreMarks"
    )
    
    if (!element_name %in% names(roadmap_bed_gz_links_metadata)) {
        stop("Data type '", data_type, "' not found in roadmap_bed_gz_links_metadata",
             call. = FALSE)
    }
    
    # Get URLs for the selected data type
    url_data <- roadmap_bed_gz_links_metadata[[element_name]]
    
    # Filter by tissue types if specified
    # url_data already contains metadata (including anatomy column) from roadmap_metadata
    if (!is.null(tissue_types)) {
        if (!"anatomy" %in% colnames(url_data)) {
            stop("anatomy column not found in url_data. ",
                 "This indicates a problem with roadmap_bed_gz_links_metadata.", 
                 call. = FALSE)
        }
        url_data <- url_data[url_data$anatomy %in% tissue_types, ]
        if (nrow(url_data) == 0) {
            stop("No epigenomes found for tissue types: ", 
                 paste(tissue_types, collapse = ", "), call. = FALSE)
        }
        message("Filtered to ", nrow(url_data), " epigenomes matching tissue types")
    }
    
    # Filter by EIDs if specified
    if (!is.null(EIDs)) {
        url_data <- url_data[url_data$EID %in% EIDs, ]
        if (nrow(url_data) == 0) {
            stop("No epigenomes found for EIDs: ", 
                 paste(EIDs, collapse = ", "), call. = FALSE)
        }
        message("Filtered to ", nrow(url_data), " epigenomes matching EIDs")
    }
    
    if (nrow(url_data) == 0) {
        stop("No epigenomes selected for download", call. = FALSE)
    }
    
    message("Downloading data for ", nrow(url_data), " epigenome(s)...")
    
    # Download and load each file
    all_regions <- list()
    
    for (i in seq_len(nrow(url_data))) {
        eid <- url_data$EID[i]
        url <- url_data$URL[i]
        
        # Create cache file name
        safe_eid <- gsub("[^A-Za-z0-9]", "_", eid)
        cache_file <- file.path(cache_dir, 
                               paste0("Roadmap_", safe_eid, "_", data_type, ".bed.gz"))
        
        # Initialize regions variable
        regions <- NULL
        regions_loaded <- FALSE
        
        # Check cache first
        if (file.exists(cache_file)) {
            message("Using cached file for ", eid, ": ", cache_file)
            tryCatch({
                regions <- loadRoadmapFromFile(cache_file)
                if (!is.null(regions) && length(regions) > 0) {
                    regions_loaded <- TRUE
                }
            }, error = function(e) {
                warning("Failed to load cached file for EID '", eid, "': ", 
                       e$message, ". Will re-download.", call. = FALSE)
            })
        }
        
        # Download if not cached or cache failed
        if (!regions_loaded) {
            # Download from Roadmap server
            message("Downloading ", eid, " from Roadmap...")
            message("  URL: ", url)
            
            tryCatch({
                download.file(url, cache_file, mode = "wb", quiet = FALSE)
                if (file.exists(cache_file) && file.info(cache_file)$size > 0) {
                    regions <- loadRoadmapFromFile(cache_file)
                    if (!is.null(regions) && length(regions) > 0) {
                        regions_loaded <- TRUE
                    } else {
                        warning("Downloaded file is empty for EID: ", eid, call. = FALSE)
                    }
                } else {
                    warning("Downloaded file is empty or missing for EID: ", eid, call. = FALSE)
                }
            }, error = function(e) {
                warning("Failed to download data for EID '", eid, "': ", 
                       e$message, call. = FALSE)
            })
        }
        
        # Skip to next iteration if regions weren't loaded
        if (!regions_loaded || is.null(regions) || length(regions) == 0) {
            next
        }
        
        # Validate that regions is a GRanges object
        if (!inherits(regions, "GRanges")) {
            warning("Regions for EID '", eid, "' is not a GRanges object. Got: ",
                   paste(class(regions), collapse = ", "), ". Skipping.", call. = FALSE)
            next
        }
        
        # Add EID and metadata to regions
            mcols(regions)$EID <- eid
            if ("anatomy" %in% colnames(url_data)) {
                mcols(regions)$anatomy <- url_data$anatomy[i]
            }
            if ("standardized_epigenome_name" %in% colnames(url_data)) {
                mcols(regions)$standardized_epigenome_name <- url_data$standardized_epigenome_name[i]
            }
            
        # Store regions (combine_regions is handled later)
        all_regions[[eid]] <- regions
    }
    
    if (length(all_regions) == 0) {
        return(GRanges())
    }
    
    # Always combine regions into a single GRanges object
    # The combine_regions parameter is handled at a higher level if needed
    # For overlap finding, we always need a single GRanges object
    if (length(all_regions) > 0) {
        # Filter out any NULL or invalid entries
        all_regions <- all_regions[!sapply(all_regions, is.null)]
        all_regions <- all_regions[sapply(all_regions, inherits, "GRanges")]
        
        if (length(all_regions) > 0) {
            message("Combining regions from ", length(all_regions), " epigenome(s)...")
            
            # Combine using c() - more reliable than do.call(c, ...)
            combined <- all_regions[[1]]
            if (length(all_regions) > 1) {
                for (i in 2:length(all_regions)) {
                    combined <- c(combined, all_regions[[i]])
                }
            }
            
            # Final validation
            if (!inherits(combined, "GRanges")) {
                stop("Failed to combine regions into GRanges. Got: ",
                     paste(class(combined), collapse = ", "), call. = FALSE)
            }
            
            return(combined)
        }
    }
    # Fallback: return empty GRanges if something went wrong
    return(GRanges())
}

#' Get liftOver chain file for genome conversion
#' @param from_genome Character string, source genome (e.g., "hg19")
#' @param to_genome Character string, target genome (e.g., "hg38")
#' @param cache_dir Character string, directory to cache chain files
#' @return Chain object from rtracklayer, or NULL if not available
#' @keywords internal
getLiftOverChain <- function(from_genome, to_genome, cache_dir) {
    # Map genome names to UCSC chain file names
    # Common conversions: hg19 -> hg38, etc.
    chain_name <- paste0(from_genome, "To", 
                        toupper(substring(to_genome, 1, 1)), 
                        substring(to_genome, 2), 
                        ".over.chain")
    
    # UCSC chain file URL base
    base_url <- "http://hgdownload.cse.ucsc.edu/goldenpath/"
    
    # Construct URL for chain file
    chain_url <- paste0(base_url, from_genome, "/liftOver/", chain_name, ".gz")
    
    # Cache file path
    chain_cache <- file.path(cache_dir, paste0(chain_name, ".gz"))
    chain_file <- file.path(cache_dir, chain_name)
    
    # Check if chain file already exists (uncompressed)
    if (file.exists(chain_file)) {
        message("Using cached chain file: ", chain_file)
        tryCatch({
            return(import.chain(chain_file))
        }, error = function(e) {
            message("Error loading cached chain file, will re-download: ", e$message)
        })
    }
    
    # Check if compressed chain file exists
    if (file.exists(chain_cache)) {
        message("Found compressed chain file, decompressing...")
        tryCatch({
            if (requireNamespace("R.utils", quietly = TRUE)) {
                R.utils::gunzip(chain_cache, destname = chain_file, remove = FALSE)
            } else {
                # Use R's built-in gzfile connection
                con_in <- gzfile(chain_cache, "rb")
                con_out <- file(chain_file, "wb")
                writeBin(readBin(con_in, "raw", n = file.info(chain_cache)$size), con_out)
                close(con_in)
                close(con_out)
            }
            if (file.exists(chain_file)) {
                return(import.chain(chain_file))
            }
        }, error = function(e) {
            message("Error decompressing chain file: ", e$message)
        })
    }
    
    # Download chain file
    message("Downloading chain file from UCSC: ", chain_url)
    tryCatch({
        download.file(chain_url, chain_cache, mode = "wb", quiet = FALSE)
        
        # Decompress
        if (requireNamespace("R.utils", quietly = TRUE)) {
            R.utils::gunzip(chain_cache, destname = chain_file, remove = FALSE)
        } else {
            # Use R's built-in gzfile connection
            file_size <- file.info(chain_cache)$size
            if (is.na(file_size) || file_size == 0) {
                stop("Chain file is empty or cannot be read")
            }
            con_in <- gzfile(chain_cache, "rb")
            con_out <- file(chain_file, "wb")
            on.exit(close(con_in), add = TRUE)
            on.exit(close(con_out), add = TRUE)
            writeBin(readBin(con_in, "raw", n = file_size), con_out)
            close(con_in)
            close(con_out)
        }
        
        # Import chain
        if (file.exists(chain_file)) {
            chain_obj <- import.chain(chain_file)
            message("Successfully loaded chain file for ", from_genome, " -> ", to_genome)
            return(chain_obj)
        } else {
            warning("Chain file was not created after download", call. = FALSE)
            return(NULL)
        }
    }, error = function(e) {
        warning("Failed to download chain file from UCSC: ", e$message, 
               ". You may need to manually download the chain file or ",
               "convert your peaks to hg19.", call. = FALSE)
        return(NULL)
    })
}

