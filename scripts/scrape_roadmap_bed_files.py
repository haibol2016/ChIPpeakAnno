#!/usr/bin/env python3
"""
Web scraping script to extract download links for all .bed files from
Roadmap Epigenomics DNase promoter BED files directory.

URL: https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_prom/
"""

import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urlparse
import sys
from typing import List, Set


def get_bed_file_links(url: str) -> List[str]:
    """
    Scrape the webpage and extract all .bed file download links.
    
    Args:
        url: The URL of the directory listing page
        
    Returns:
        List of full URLs to .bed files (excluding .bed.gz files)
    """
    try:
        # Send GET request to the URL
        print(f"Fetching page: {url}")
        response = requests.get(url, timeout=30)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Parse the HTML content
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find all anchor tags (links)
        links = soup.find_all('a')
        
        bed_urls = []
        base_url = url.rstrip('/')
        
        for link in links:
            href = link.get('href', '')
            if not href:
                continue
            
            # Decode URL-encoded characters (e.g., %5F -> _)
            from urllib.parse import unquote
            decoded_href = unquote(href)
            
            # Check if it's a .bed file (but not .bed.gz)
            if decoded_href.endswith('.bed') and not decoded_href.endswith('.bed.gz'):
                # Construct full URL
                if href.startswith('http'):
                    full_url = href
                else:
                    # Handle relative URLs
                    full_url = urljoin(base_url + '/', href)
                
                bed_urls.append(full_url)
        
        # Remove duplicates while preserving order
        seen: Set[str] = set()
        unique_bed_urls = []
        for url_item in bed_urls:
            if url_item not in seen:
                seen.add(url_item)
                unique_bed_urls.append(url_item)
        
        return sorted(unique_bed_urls)
    
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the webpage: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing the webpage: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """Main function to run the scraper."""
    url = "https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_prom/"
    
    print("=" * 70)
    print("Roadmap Epigenomics DNase Promoter BED Files Scraper")
    print("=" * 70)
    print()
    
    # Get all .bed file links
    bed_links = get_bed_file_links(url)
    
    if not bed_links:
        print("No .bed files found on the page.")
        sys.exit(1)
    
    print(f"Found {len(bed_links)} .bed file(s):\n")
    
    # Print all links
    for i, link in enumerate(bed_links, 1):
        print(f"{i:3d}. {link}")
    
    print()
    print("=" * 70)
    print(f"Total: {len(bed_links)} .bed file(s)")
    print("=" * 70)
    
    # Optionally save to a file
    output_file = "roadmap_bed_file_links.txt"
    try:
        with open(output_file, 'w') as f:
            for link in bed_links:
                f.write(f"{link}\n")
        print(f"\nLinks saved to: {output_file}")
    except Exception as e:
        print(f"\nWarning: Could not save to file: {e}", file=sys.stderr)
    
    return bed_links


if __name__ == "__main__":
    main()

