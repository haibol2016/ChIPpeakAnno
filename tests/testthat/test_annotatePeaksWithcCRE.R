test_that("annotatePeaksWithcCRE works with custom file", {
    # Create test peaks
    peaks <- GRanges("chr1", IRanges(c(1000, 5000, 10000), width=500))
    
    # Create mock cCRE data (BED format)
    ccres_bed <- data.frame(
        chrom = c("chr1", "chr1", "chr1"),
        start = c(500, 4800, 9800),
        end = c(1500, 5200, 10200),
        name = c("cCRE1", "cCRE2", "cCRE3"),
        score = c(100, 200, 150),
        strand = c("+", "-", "+"),
        stringsAsFactors = FALSE
    )
    
    # Write to temporary file
    temp_bed <- tempfile(fileext = ".bed")
    write.table(ccres_bed, temp_bed, sep = "\t", quote = FALSE, 
               row.names = FALSE, col.names = FALSE)
    
    # Test annotation
    annotated <- annotatePeaksWithcCRE(peaks, cCRE_file = temp_bed)
    
    # Check that metadata columns were added
    expect_true("has_cCRE_overlap" %in% colnames(mcols(annotated)))
    expect_true("cCRE_count" %in% colnames(mcols(annotated)))
    expect_true("cCRE_types" %in% colnames(mcols(annotated)))
    
    # Check that overlaps were found
    expect_true(any(annotated$has_cCRE_overlap))
    expect_true(sum(annotated$has_cCRE_overlap) > 0)
    
    # Clean up
    unlink(temp_bed)
})

test_that("annotatePeaksWithcCRE handles empty overlaps", {
    # Create test peaks that don't overlap with cCREs
    peaks <- GRanges("chr1", IRanges(c(100000, 200000), width=100))
    
    # Create mock cCRE data far away
    ccres_bed <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 2000,
        name = "cCRE1",
        score = 100,
        strand = "+",
        stringsAsFactors = FALSE
    )
    
    temp_bed <- tempfile(fileext = ".bed")
    write.table(ccres_bed, temp_bed, sep = "\t", quote = FALSE, 
               row.names = FALSE, col.names = FALSE)
    
    annotated <- annotatePeaksWithcCRE(peaks, cCRE_file = temp_bed)
    
    # Check that no overlaps were found
    expect_false(any(annotated$has_cCRE_overlap))
    expect_true(all(annotated$cCRE_count == 0))
    
    unlink(temp_bed)
})

test_that("annotatePeaksWithcCRE preserves existing metadata", {
    # Create peaks with existing metadata
    peaks <- GRanges("chr1", IRanges(c(1000, 5000), width=500))
    mcols(peaks)$gene_id <- c("gene1", "gene2")
    mcols(peaks)$distance <- c(100, 200)
    
    # Create mock cCRE data
    ccres_bed <- data.frame(
        chrom = "chr1",
        start = 500,
        end = 1500,
        name = "cCRE1",
        score = 100,
        strand = "+",
        stringsAsFactors = FALSE
    )
    
    temp_bed <- tempfile(fileext = ".bed")
    write.table(ccres_bed, temp_bed, sep = "\t", quote = FALSE, 
               row.names = FALSE, col.names = FALSE)
    
    annotated <- annotatePeaksWithcCRE(peaks, cCRE_file = temp_bed)
    
    # Check that existing metadata is preserved
    expect_true("gene_id" %in% colnames(mcols(annotated)))
    expect_true("distance" %in% colnames(mcols(annotated)))
    expect_equal(annotated$gene_id, c("gene1", "gene2"))
    expect_equal(annotated$distance, c(100, 200))
    
    unlink(temp_bed)
})

test_that("annotatePeaksWithcCRE validates inputs", {
    # Test invalid peaks input
    expect_error(annotatePeaksWithcCRE("not_granges"),
                "'peaks' must be a GRanges object")
    
    # Test invalid species
    peaks <- GRanges("chr1", IRanges(1000, width=500))
    expect_error(annotatePeaksWithcCRE(peaks, species = "invalid"),
                "should be one of")
})

test_that("listAvailablecCREs returns data.frame", {
    # Test human
    result_human <- listAvailablecCREs("Homo sapiens")
    expect_s3_class(result_human, "data.frame")
    expect_true(nrow(result_human) > 0)
    expect_true("tissue_type" %in% colnames(result_human))
    expect_true("description" %in% colnames(result_human))
    
    # Test mouse
    result_mouse <- listAvailablecCREs("Mus musculus")
    expect_s3_class(result_mouse, "data.frame")
    expect_true(nrow(result_mouse) > 0)
    expect_true("tissue_type" %in% colnames(result_mouse))
    
    # Test invalid species
    expect_error(listAvailablecCREs("invalid"),
                "should be one of")
})

test_that("annotatePeaksWithcCRE handles maxgap parameter", {
    # Create peaks
    peaks <- GRanges("chr1", IRanges(c(1000, 5000), width=500))
    
    # Create cCRE data with gap
    ccres_bed <- data.frame(
        chrom = "chr1",
        start = 1600,  # Gap of 100bp from first peak end (1500)
        end = 2000,
        name = "cCRE1",
        score = 100,
        strand = "+",
        stringsAsFactors = FALSE
    )
    
    temp_bed <- tempfile(fileext = ".bed")
    write.table(ccres_bed, temp_bed, sep = "\t", quote = FALSE, 
               row.names = FALSE, col.names = FALSE)
    
    # Without maxgap, no overlap
    annotated_no_gap <- annotatePeaksWithcCRE(peaks, cCRE_file = temp_bed)
    expect_false(annotated_no_gap$has_cCRE_overlap[1])
    
    # With maxgap, should find overlap
    annotated_with_gap <- annotatePeaksWithcCRE(peaks, cCRE_file = temp_bed, 
                                                maxgap = 200)
    expect_true(annotated_with_gap$has_cCRE_overlap[1])
    
    unlink(temp_bed)
})

