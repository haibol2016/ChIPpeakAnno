#!/usr/bin/env Rscript
# Web scraping script to extract download links for all .bed.gz files from
# Roadmap Epigenomics DNase promoter BED files directory.
#
# URL: https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_prom/
#
# Usage:
#   Rscript scrape_roadmap_bed_files.R
#   or
#   source("scrape_roadmap_bed_files.R")
#   bed_links <- scrape_roadmap_bed_files()

suppressPackageStartupMessages({
    if (!require("rvest", quietly = TRUE)) {
        stop("Package 'rvest' is required. Install it with: install.packages('rvest')")
    }
    if (!require("httr", quietly = TRUE)) {
        stop("Package 'httr' is required. Install it with: install.packages('httr')")
    }
})

#' Scrape the webpage and extract all .bed.gz file download links
#'
#' @param url The URL of the directory listing page
#' @return Character vector of full URLs to .bed.gz files (excluding plain .bed files)
get_bed_file_links <- function(url) {
    tryCatch({
        # Send GET request to the URL
        cat("Fetching page:", url, "\n")
        
        # Read the HTML page
        page <- read_html(url)
        
        # Find all anchor tags (links)
        links <- page %>%
            html_nodes("a") %>%
            html_attr("href")
        
        # Remove NA values
        links <- links[!is.na(links)]
        
        # Base URL for constructing full URLs
        base_url <- sub("/$", "", url)
        
        bed_urls <- character(0)
        
        for (href in links) {
            if (nchar(href) == 0) next
            
            # Decode URL-encoded characters (e.g., %5F -> _)
            decoded_href <- URLdecode(href)
            
            # Check if it's a .bed.gz file (excluding plain .bed files)
            if (grepl("\\.bed\\.gz$", decoded_href) && !grepl("hg38lift|stateno|expanded|mnemonics", decoded_href)) {
                # Construct full URL
                if (grepl("^https?://", href)) {
                    full_url <- href
                } else {
                    # Handle relative URLs
                    if (startsWith(href, "/")) {
                        # Absolute path from domain root
                        parsed_url <- httr::parse_url(url)
                        full_url <- paste0(parsed_url$scheme, "://", parsed_url$hostname, href)
                    } else {
                        # Relative path
                        full_url <- paste0(base_url, "/", href)
                    }
                }
                
                bed_urls <- c(bed_urls, full_url)
            }
        }
        
        # Remove duplicates while preserving order
        unique_bed_urls <- unique(bed_urls)
        
        # Sort the URLs
        return(sort(unique_bed_urls))
        
    }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n", file = stderr())
        stop("Failed to scrape the webpage")
    })
}

#' Main function to run the scraper
#'
#' @param url The URL to scrape (default: Roadmap Epigenomics DNase promoter, enhancer, and chromhmm directory)
#' @param output_file Optional output file to save links (default: "roadmap_bed_gz_file_links_prom.txt", "roadmap_bed_gz_file_links_enh.txt", "roadmap_bed_gz_file_links_chromhmm.txt")
#' @param verbose Whether to print progress messages (default: TRUE)
#' @return Character vector of .bed.gz file URLs
scrape_roadmap_bed_files <- function(
    url = c("https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_prom/",
    "https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_enh/",
    "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/"),
    verbose = TRUE
) {
    url <- match.arg(url)

    if (verbose) {
        cat(paste(rep("=", 70), collapse = ""), "\n")
        cat(paste("Roadmap Epigenomics ", url, "BED.GZ Files Scraper\n"))
        cat("All files are based on hg19 assembly.\n")
        cat("If you need to use a different assembly, you need to lift over the files to the\n
                desired assembly using `rtracklayer::liftOver()`.\n")
        cat(paste(rep("=", 70), collapse = ""), "\n")
        cat("\n")
    }
    
    # Get all .bed.gz file links
    bed_links <- get_bed_file_links(url)
    
    if (length(bed_links) == 0) {
        cat("No .bed.gz files found on the page.\n")
        return(character(0))
    }
    
    if (verbose) {
        cat("Found", length(bed_links), ".bed.gz file(s):\n\n")
    }
    
    return(bed_links)
}

data("roadmap_metadata")
url <- c("https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_prom/",
"https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_enh/",
"https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/")
verbose <- TRUE
roadmap_bed_gz_links_metadata <- lapply(url, function(.url) {
    urls <- scrape_roadmap_bed_files(.url)
    EID <- gsub(".+(E\\d+).*", "\\1", urls)
    urls <- data.frame(EID = EID, URL = urls, stringsAsFactors = FALSE)
    urls <- merge(urls, roadmap_metadata, by.x = "EID", 
                  by.y = "epigenome_id_eid", all.x = TRUE)
})
names(roadmap_bed_gz_links_metadata) <- c("DHS_promoter", 
"DHS_enhancer", "ChromHMM_15_states_5_coreMarks")
saveRDS(roadmap_bed_gz_links_metadata, "roadmap_bed_gz_links_metadata.rds")

usethis::use_data(roadmap_bed_gz_links_metadata, overwrite = TRUE)

