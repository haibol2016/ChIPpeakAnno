#' Global Options for ChIPpeakAnno
#' 
#' @description
#' Internal functions to set and retrieve global options for ChIPpeakAnno,
#' including default genome assembly, species, TxDb, and EnsDb objects.
#' These options provide a convenient way to set defaults that are used across
#' multiple functions, while still allowing per-function overrides.
#' 
#' @details
#' Most functions in this file are internal and not exported. Users should use
#' \code{\link{setChIPpeakAnnoGlobals}} to set all global options at once.
#' 
#' @name ChIPpeakAnno-globals
#' @aliases setChIPpeakAnnoGenome getChIPpeakAnnoGenome setChIPpeakAnnoTxDb getChIPpeakAnnoTxDb setChIPpeakAnnoEnsDb getChIPpeakAnnoEnsDb setChIPpeakAnnoSpecies getChIPpeakAnnoSpecies getTargetGenome setChIPpeakAnnoGlobals
#' @keywords internal
NULL

#' Set global genome assembly for ChIPpeakAnno functions
#' 
#' @description
#' Sets a global default genome assembly that will be used by ChIPpeakAnno
#' functions when genome information cannot be detected from the input peaks.
#' This provides a convenient fallback mechanism and ensures consistency
#' across analyses.
#' 
#' @param genome Character string specifying the genome assembly (e.g., "hg38",
#'        "hg19", "mm10", "mm39"). Should be a normalized genome identifier
#'        that matches UCSC naming conventions.
#' 
#' @details
#' The global genome option is used as a fallback when:
#' \itemize{
#'   \item Peaks do not have genome metadata set
#'   \item Genome detection from peaks fails
#'   \item Explicit genome parameter is not provided to a function
#' }
#' 
#' The global option can always be overridden by explicitly setting the genome
#' in individual function calls or by ensuring peaks have proper genome metadata.
#' 
#' @return Invisibly returns the previous value of the global genome option
#' 
#' @keywords internal
#' 
#' @examples
#' # Set global genome to hg38
#' setChIPpeakAnnoGenome("hg38")
#' 
#' # Get the current global genome
#' getChIPpeakAnnoGenome()
#' 
#' # Clear the global genome (set to NULL)
#' setChIPpeakAnnoGenome(NULL)
setChIPpeakAnnoGenome <- function(genome = NULL) {
    if (!is.null(genome) && !is.character(genome)) {
        stop("genome must be NULL or a character string", call. = FALSE)
    }
    
    if (!is.null(genome) && length(genome) != 1) {
        stop("genome must be a single character string", call. = FALSE)
    }
    
    old_genome <- getOption("ChIPpeakAnno.genome")
    
    if (!is.null(genome)) {
        message("Setting global genome assembly to: ", genome)
        options(ChIPpeakAnno.genome = genome)
    } else {
        options(ChIPpeakAnno.genome = NULL)
    }
    
    invisible(old_genome)
}

#' Get global genome assembly for ChIPpeakAnno functions
#' 
#' @description
#' Retrieves the currently set global genome assembly option for peaks annotation.
#' 
#' @return Character string with the genome assembly, or NULL if not set
#' 
#' @keywords internal
#' 
#' @examples
#' # Get current global genome (may be NULL)
#' getChIPpeakAnnoGenome()
#' 
#' # Set and retrieve
#' setChIPpeakAnnoGenome("hg38")
#' getChIPpeakAnnoGenome()
getChIPpeakAnnoGenome <- function() {
    getOption("ChIPpeakAnno.genome")
}

#' Set global TxDb for ChIPpeakAnno functions
#' 
#' @description
#' Sets a global default TxDb (or EnsDb) object that will be used by ChIPpeakAnno
#' functions when a TxDb/EnsDb parameter is not explicitly provided. This provides
#' a convenient way to set defaults that are used across multiple functions, while
#' still allowing per-function overrides.
#' 
#' @param TxDb An object of class \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#'        or \code{\link[ensembldb:EnsDb-class]{EnsDb}} containing genomic
#'        annotation. If NULL, clears the global TxDb option.
#' 
#' @details
#' The global TxDb option is used as a fallback when:
#' \itemize{
#'   \item Functions require a TxDb/EnsDb parameter but it's not provided
#'   \item Functions have \code{TxDb = NULL} or \code{EnsDb = NULL} as default
#' }
#' 
#' The global option can always be overridden by explicitly providing TxDb/EnsDb
#' in individual function calls.
#' 
#' @return Invisibly returns the previous value of the global TxDb option
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' 
#' # Set global TxDb
#' setChIPpeakAnnoTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' 
#' # Get the current global TxDb
#' getChIPpeakAnnoTxDb()
#' 
#' # Clear the global TxDb (set to NULL)
#' setChIPpeakAnnoTxDb(NULL)
#' }
setChIPpeakAnnoTxDb <- function(TxDb = NULL) {
    if (!is.null(TxDb) && !inherits(TxDb, c("TxDb", "EnsDb"))) {
        stop("TxDb must be NULL or an object of class TxDb or EnsDb", call. = FALSE)
    }
    
    old_txdb <- getOption("ChIPpeakAnno.TxDb")
    
    if (!is.null(TxDb)) {
        txdb_name <- if (inherits(TxDb, "TxDb")) {
            "TxDb"
        } else {
            "EnsDb"
        }
        message("Setting global ", txdb_name, " to: ", deparse(substitute(TxDb)))
        options(ChIPpeakAnno.TxDb = TxDb)
    } else {
        options(ChIPpeakAnno.TxDb = NULL)
    }
    
    invisible(old_txdb)
}

#' Get global TxDb for ChIPpeakAnno functions
#' 
#' @description
#' Retrieves the currently set global TxDb (or EnsDb) option for peaks annotation.
#' 
#' @return An object of class \code{TxDb} or \code{EnsDb}, or NULL if not set
#' 
#' @keywords internal
#' 
#' @examples
#' # Get current global TxDb (may be NULL)
#' getChIPpeakAnnoTxDb()
#' 
#' \dontrun{
#' # Set and retrieve
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' setChIPpeakAnnoTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' getChIPpeakAnnoTxDb()
#' }
getChIPpeakAnnoTxDb <- function() {
    getOption("ChIPpeakAnno.TxDb")
}

#' Set global EnsDb for ChIPpeakAnno functions
#' 
#' @description
#' Sets a global default EnsDb object that will be used by ChIPpeakAnno
#' functions when an EnsDb parameter is not explicitly provided. This provides
#' a convenient way to set defaults that are used across multiple functions, while
#' still allowing per-function overrides.
#' 
#' @param EnsDb An object of class \code{\link[ensembldb:EnsDb-class]{EnsDb}}
#'        containing genomic annotation. If NULL, clears the global EnsDb option.
#' 
#' @details
#' The global EnsDb option is used as a fallback when:
#' \itemize{
#'   \item Functions require an EnsDb parameter but it's not provided
#'   \item Functions have \code{EnsDb = NULL} as default
#' }
#' 
#' The global option can always be overridden by explicitly providing EnsDb
#' in individual function calls.
#' 
#' @return Invisibly returns the previous value of the global EnsDb option
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v86)
#' 
#' # Set global EnsDb
#' setChIPpeakAnnoEnsDb(EnsDb.Hsapiens.v86)
#' 
#' # Get the current global EnsDb
#' getChIPpeakAnnoEnsDb()
#' 
#' # Clear the global EnsDb (set to NULL)
#' setChIPpeakAnnoEnsDb(NULL)
#' }
setChIPpeakAnnoEnsDb <- function(EnsDb = NULL) {
    if (!is.null(EnsDb) && !inherits(EnsDb, "EnsDb")) {
        stop("EnsDb must be NULL or an object of class EnsDb", call. = FALSE)
    }
    
    old_ensdb <- getOption("ChIPpeakAnno.EnsDb")
    
    if (!is.null(EnsDb)) {
        message("Setting global EnsDb to: ", deparse(substitute(EnsDb)))
        options(ChIPpeakAnno.EnsDb = EnsDb)
    } else {
        options(ChIPpeakAnno.EnsDb = NULL)
    }
    
    invisible(old_ensdb)
}

#' Get global EnsDb for ChIPpeakAnno functions
#' 
#' @description
#' Retrieves the currently set global EnsDb option for peaks annotation.
#' 
#' @return An object of class \code{EnsDb}, or NULL if not set
#' 
#' @keywords internal
#' 
#' @examples
#' # Get current global EnsDb (may be NULL)
#' getChIPpeakAnnoEnsDb()
#' 
#' \dontrun{
#' # Set and retrieve
#' library(EnsDb.Hsapiens.v86)
#' setChIPpeakAnnoEnsDb(EnsDb.Hsapiens.v86)
#' getChIPpeakAnnoEnsDb()
#' }
getChIPpeakAnnoEnsDb <- function() {
    getOption("ChIPpeakAnno.EnsDb")
}

#' Set global species for ChIPpeakAnno functions
#' 
#' @description
#' Sets a global default species that will be used by ChIPpeakAnno functions
#' when species is not explicitly provided. This provides convenience for
#' users who consistently work with one species.
#' 
#' @param species Character string specifying the species (e.g., "Homo sapiens",
#'        "Mus musculus"). Should match the species names used in ChIPpeakAnno
#'        functions.
#' 
#' @details
#' The global species option is used as a default when:
#' \itemize{
#'   \item Species parameter is not provided to a function
#'   \item Function requires species but user wants to use global default
#' }
#' 
#' The global option can always be overridden by explicitly setting the species
#' in individual function calls.
#' 
#' @return Invisibly returns the previous value of the global species option
#' 
#' @keywords internal
#' 
#' @examples
#' # Set global species to human
#' setChIPpeakAnnoSpecies("Homo sapiens")
#' 
#' # Get the current global species
#' getChIPpeakAnnoSpecies()
#' 
#' # Clear the global species (set to NULL)
#' setChIPpeakAnnoSpecies(NULL)
setChIPpeakAnnoSpecies <- function(species = NULL) {
    if (!is.null(species) && !is.character(species)) {
        stop("species must be NULL or a character string", call. = FALSE)
    }
    
    if (!is.null(species) && length(species) != 1) {
        stop("species must be a single character string", call. = FALSE)
    }
    
    # Validate that species follows Latin binomial format (Genus species)
    if (!is.null(species)) {
        # Latin binomial names typically have: Capitalized genus + lowercase species
        # Pattern: starts with capital letter, space, then lowercase letters
        if (!grepl("^[A-Z][a-z]+ [a-z]+", species)) {
            stop("Species name '", species, "' does not appear to follow Latin binomial ",
                   "format (Genus species, e.g., 'Homo sapiens', 'Mus musculus'). ",
                   "Please use the scientific name format.", call. = FALSE)
        }
    }
    
    old_species <- getOption("ChIPpeakAnno.species")
    
    if (!is.null(species)) {
        message("Setting global species to: ", species)
        options(ChIPpeakAnno.species = species)
    } else {
        options(ChIPpeakAnno.species = NULL)
    }
    
    invisible(old_species)
}

#' Get global species for ChIPpeakAnno functions
#' 
#' @description
#' Retrieves the currently set global species option for peaks annotation.
#' 
#' @return Character string with the species, or NULL if not set
#' 
#' @keywords internal
#' 
#' @examples
#' # Get current global species (may be NULL)
#' getChIPpeakAnnoSpecies()
#' 
#' # Set and retrieve
#' setChIPpeakAnnoSpecies("Homo sapiens")
#' getChIPpeakAnnoSpecies()
getChIPpeakAnnoSpecies <- function() {
    getOption("ChIPpeakAnno.species")
}

#' Detect target genome assembly from peaks or global option
#' 
#' @description
#' Helper function that detects the target genome assembly using a hybrid
#' approach: first checks the global option, then attempts to detect from
#' peaks, and finally falls back to the global option if detection fails.
#' This function normalizes genome names to standard UCSC identifiers.
#' 
#' @param peaks A \code{\link[GenomicRanges]{GRanges}} object containing peaks
#' @param species Optional character string, species name. When provided,
#'        enables species-specific genome normalization (e.g., "Homo sapiens",
#'        "Mus musculus"). When NULL, uses global species option if set, or
#'        works with any genome assembly.
#' 
#' @return Character string with normalized genome assembly (e.g., "hg38", "hg19",
#'         "mm10", "mm39", or any other genome identifier), or NULL if detection
#'         fails and no global option is set
#' 
#' @details
#' Detection priority:
#' \enumerate{
#'   \item Global option (\code{getChIPpeakAnnoGenome()}) - if species is provided,
#'         validates that global genome matches species
#'   \item Detection from peaks metadata (\code{genome(peaks)})
#'   \item Fallback to global option if detection fails
#' }
#' 
#' When \code{species} is provided, the function performs species-specific
#' normalization:
#' \itemize{
#'   \item Human: "hg19", "GRCh37", "Hsapiens.UCSC.hg19" → "hg19"
#'   \item Human: "hg38", "GRCh38", "Hsapiens.UCSC.hg38" → "hg38"
#'   \item Mouse: "mm9", "GRCm37", "Mmusculus.UCSC.mm9" → "mm9"
#'   \item Mouse: "mm10", "GRCm38", "Mmusculus.UCSC.mm10" → "mm10"
#'   \item Mouse: "mm39", "GRCm39", "Mmusculus.UCSC.mm39" → "mm39"
#' }
#' 
#' When \code{species} is NULL, the function attempts to extract any genome
#' identifier from the peaks metadata and returns it as-is, or uses the global
#' option without species validation.
#' 
#' @keywords internal
#' @importFrom GenomeInfoDb genome
getTargetGenome <- function(peaks, species = NULL) {
    # Use global species as default if species not provided
    if (is.null(species)) {
        species <- getChIPpeakAnnoSpecies()
    }
    
    # 1. Check global option first
    global_genome <- getChIPpeakAnnoGenome()
    if (!is.null(global_genome)) {
        # If species is provided, validate global genome matches species
        if (!is.null(species)) {
            if (species == "Homo sapiens" && grepl("^hg", global_genome, ignore.case = TRUE)) {
                return(global_genome)
            } else if (species == "Mus musculus" && grepl("^mm", global_genome, ignore.case = TRUE)) {
                return(global_genome)
            } else {
                warning("species is not human or mouse, using global genome assembly: ", global_genome)
                return(global_genome)
            }
        } else {
            # No species specified, use global genome directly
            return(global_genome)
        }
    }
    
    # 2. Try to detect from peaks
    target_genome <- tryCatch({
        genome_info <- unique(genome(peaks))
        genome_info <- genome_info[!is.na(genome_info)]
        if (length(genome_info) == 0) {
            NULL
        } else {
            genome_str <- paste(genome_info, collapse = " ")
            
            # If species is provided, do species-specific normalization
            if (!is.null(species)) {
                if (species == "Homo sapiens") {
                    if (grepl("hg19|GRCh37|Hsapiens.UCSC.hg19", genome_str, ignore.case = TRUE)) {
                        "hg19"
                    } else if (grepl("hg38|GRCh38|Hsapiens.UCSC.hg38", genome_str, ignore.case = TRUE)) {
                        "hg38"
                    } else {
                        # Try to extract from genome string
                        if (grepl("hg\\d+", genome_str, ignore.case = TRUE)) {
                            regmatches(genome_str, regexpr("hg\\d+", genome_str, ignore.case = TRUE))
                        } else {
                            NULL
                        }
                    }
                } else if (species == "Mus musculus") {
                    if (grepl("mm9|GRCm37|Mmusculus.UCSC.mm9", genome_str, ignore.case = TRUE)) {
                        "mm9"
                    } else if (grepl("mm10|GRCm38|Mmusculus.UCSC.mm10", genome_str, ignore.case = TRUE)) {
                        "mm10"
                    } else if (grepl("mm39|GRCm39|Mmusculus.UCSC.mm39", genome_str, ignore.case = TRUE)) {
                        "mm39"
                    } else {
                        # Try to extract from genome string
                        if (grepl("mm\\d+", genome_str, ignore.case = TRUE)) {
                            regmatches(genome_str, regexpr("mm\\d+", genome_str, ignore.case = TRUE))
                        } else {
                            NULL
                        }
                    }
                } else {
                    # Unknown species, try generic extraction
                    if (grepl("(hg|mm|rn|dm|ce|danRer|galGal|susScr|bosTau)\\d+", genome_str, ignore.case = TRUE)) {
                        regmatches(genome_str, regexpr("(hg|mm|rn|dm|ce|danRer|galGal|susScr|bosTau)\\d+", genome_str, ignore.case = TRUE))
                    } else {
                        # Return as-is if no pattern matches
                        genome_info[1]
                    }
                }
            } else {
                # No species specified, try to extract common genome identifiers
                if (grepl("(hg|mm|rn|dm|ce|danRer|galGal|susScr|bosTau)\\d+", genome_str, ignore.case = TRUE)) {
                    regmatches(genome_str, regexpr("(hg|mm|rn|dm|ce|danRer|galGal|susScr|bosTau)\\d+", genome_str, ignore.case = TRUE))
                } else {
                    # Return first genome info as-is
                    genome_info[1]
                }
            }
        }
    }, error = function(e) NULL)
    
    # Return detected genome (or NULL if detection failed)
    # Note: If global_genome existed, we would have returned in step 1,
    # so no need for fallback here
    return(target_genome)
}

#' Set up global variables for ChIPpeakAnno analysis
#' 
#' @description
#' A convenience function to set all ChIPpeakAnno global options at once.
#' This function provides a single entry point for configuring the default
#' genome assembly, species, TxDb, and EnsDb objects that will be used
#' across ChIPpeakAnno functions.
#' 
#' @param genome Character string specifying the genome assembly (e.g., "hg38",
#'        "hg19", "mm10", "mm39"). Should be a normalized genome identifier
#'        that matches UCSC naming conventions. If \code{NULL}, the global
#'        genome option is cleared.
#' @param species Character string specifying the species (e.g., "Homo sapiens",
#'        "Mus musculus"). Should match the species names used in ChIPpeakAnno
#'        functions. If \code{NULL}, the global species option is cleared.
#' @param TxDb An object of \link[GenomicFeatures:TxDb-class]{TxDb-class} or
#'        \link[ensembldb:EnsDb-class]{EnsDb-class} containing transcript and
#'        gene annotations. If \code{NULL}, the global TxDb option is cleared.
#' @param EnsDb An object of \link[ensembldb:EnsDb-class]{EnsDb-class}
#'        containing genome annotation data (genes, transcripts, exons, etc.).
#'        If \code{NULL}, the global EnsDb option is cleared.
#' 
#' @details
#' This function calls the individual setter functions:
#' \itemize{
#'   \item \code{\link{setChIPpeakAnnoGenome}} for genome assembly
#'   \item \code{\link{setChIPpeakAnnoSpecies}} for species
#'   \item \code{\link{setChIPpeakAnnoTxDb}} for TxDb/EnsDb objects
#'   \item \code{\link{setChIPpeakAnnoEnsDb}} for EnsDb objects
#' }
#' 
#' All parameters are optional. You can set only the options you need, and
#' leave others as \code{NULL} to keep their current values or clear them.
#' 
#' @return Invisibly returns \code{0}
#' 
#' @export
#' 
#' @seealso
#' \code{\link{setChIPpeakAnnoGenome}} for setting genome only,
#' \code{\link{setChIPpeakAnnoSpecies}} for setting species only,
#' \code{\link{setChIPpeakAnnoTxDb}} for setting TxDb only,
#' \code{\link{setChIPpeakAnnoEnsDb}} for setting EnsDb only
#' 
#' @examples
#' # Set up all global options for human hg38 analysis
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(EnsDb.Hsapiens.v86)
#' 
#' setChIPpeakAnnoGlobals(
#'     genome = "hg38",
#'     species = "Homo sapiens",
#'     TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'     EnsDb = EnsDb.Hsapiens.v86
#' )
#' 
#' # Set only genome and species (leave TxDb/EnsDb unchanged)
#' setChIPpeakAnnoGlobals(
#'     genome = "mm10",
#'     species = "Mus musculus"
#' )
#' 
#' # Clear all global options
#' setChIPpeakAnnoGlobals(
#'     genome = NULL,
#'     species = NULL,
#'     TxDb = NULL,
#'     EnsDb = NULL
#' )
setChIPpeakAnnoGlobals <- function(genome = NULL,
                                   species = NULL,
                                   TxDb = NULL,
                                   EnsDb = NULL) {
    setChIPpeakAnnoGenome(genome)
    setChIPpeakAnnoSpecies(species)
    setChIPpeakAnnoTxDb(TxDb)
    setChIPpeakAnnoEnsDb(EnsDb)
    invisible(0)
}

