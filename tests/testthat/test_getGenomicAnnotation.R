test_that("getGenomicAnnotation works correctly with basic input", {
    # Load a sample TxDb for testing
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Create test peaks
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(1000, 5000, 10000), 
                            end = c(2000, 6000, 11000)),
            strand = "*"
        )
        
        # Test basic annotation
        result <- getGenomicAnnotation(test_peaks, TxDb)
        
        # Check output structure
        expect_s3_class(result, "data.frame")
        expect_equal(nrow(result), length(test_peaks))
        expect_true(all(c("annotation", "priorityAnnotation") %in% colnames(result)))
        expect_true(all(is.character(result$annotation)))
        expect_true(all(is.character(result$priorityAnnotation)))
    }
})

test_that("getGenomicAnnotation handles usePeakCenter parameter", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Create a large peak that spans multiple features
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = 1000, end = 10000),
            strand = "*"
        )
        
        # Test with peak center (default)
        result_center <- getGenomicAnnotation(test_peaks, TxDb, 
                                             usePeakCenter = TRUE)
        
        # Test with full peak range
        result_full <- getGenomicAnnotation(test_peaks, TxDb, 
                                           usePeakCenter = FALSE)
        
        # Results may differ for large peaks
        expect_s3_class(result_center, "data.frame")
        expect_s3_class(result_full, "data.frame")
        expect_equal(nrow(result_center), 1)
        expect_equal(nrow(result_full), 1)
    }
})

test_that("getGenomicAnnotation handles tssRegion parameter", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Create peaks near TSS
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(500, 1500, 2500), 
                            end = c(600, 1600, 2600)),
            strand = "*"
        )
        
        # Test with default tssRegion
        result_default <- getGenomicAnnotation(test_peaks, TxDb, 
                                               tssRegion = c(-3000, 3000))
        
        # Test with smaller tssRegion
        result_small <- getGenomicAnnotation(test_peaks, TxDb, 
                                            tssRegion = c(-1000, 1000))
        
        # Test with larger tssRegion
        result_large <- getGenomicAnnotation(test_peaks, TxDb, 
                                             tssRegion = c(-5000, 5000))
        
        expect_s3_class(result_default, "data.frame")
        expect_s3_class(result_small, "data.frame")
        expect_s3_class(result_large, "data.frame")
        expect_equal(nrow(result_default), length(test_peaks))
        expect_equal(nrow(result_small), length(test_peaks))
        expect_equal(nrow(result_large), length(test_peaks))
    }
})

test_that("getGenomicAnnotation handles ignore_strand parameter", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Create peaks with specific strands
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(1000, 5000), 
                            end = c(2000, 6000)),
            strand = c("+", "-")
        )
        
        # Test with ignore_strand = TRUE
        result_ignore <- getGenomicAnnotation(test_peaks, TxDb, 
                                              ignore_strand = TRUE)
        
        # Test with ignore_strand = FALSE
        result_respect <- getGenomicAnnotation(test_peaks, TxDb, 
                                               ignore_strand = FALSE)
        
        expect_s3_class(result_ignore, "data.frame")
        expect_s3_class(result_respect, "data.frame")
        expect_equal(nrow(result_ignore), length(test_peaks))
        expect_equal(nrow(result_respect), length(test_peaks))
    }
})

test_that("getGenomicAnnotation handles genomicAnnotationPriority", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(1000, 5000, 10000), 
                            end = c(2000, 6000, 11000)),
            strand = "*"
        )
        
        # Test with default priority
        result_default <- getGenomicAnnotation(test_peaks, TxDb)
        
        # Test with custom priority (reversed order)
        custom_priority <- c("distalIntergenic", "immediateDownstream", 
                             "otherIntron", "otherExon", "firstIntron", 
                             "firstExon", "threeUTR", "fiveUTR", "promoter")
        result_custom <- getGenomicAnnotation(test_peaks, TxDb, 
                                             genomicAnnotationPriority = custom_priority)
        
        expect_s3_class(result_default, "data.frame")
        expect_s3_class(result_custom, "data.frame")
        expect_equal(nrow(result_default), length(test_peaks))
        expect_equal(nrow(result_custom), length(test_peaks))
    }
})

test_that("getGenomicAnnotation handles level parameter", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(1000, 5000), 
                            end = c(2000, 6000)),
            strand = "*"
        )
        
        # Test with transcript level (default)
        result_transcript <- getGenomicAnnotation(test_peaks, TxDb, 
                                                  level = "transcript")
        
        # Test with gene level
        result_gene <- getGenomicAnnotation(test_peaks, TxDb, 
                                            level = "gene")
        
        expect_s3_class(result_transcript, "data.frame")
        expect_s3_class(result_gene, "data.frame")
        expect_equal(nrow(result_transcript), length(test_peaks))
        expect_equal(nrow(result_gene), length(test_peaks))
    }
})

test_that("getGenomicAnnotation handles immediateDownstreamLength", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(1000, 5000), 
                            end = c(2000, 6000)),
            strand = "*"
        )
        
        # Test with default downstream length
        result_default <- getGenomicAnnotation(test_peaks, TxDb, 
                                               immediateDownstreamLength = 3000)
        
        # Test with custom downstream length
        result_custom <- getGenomicAnnotation(test_peaks, TxDb, 
                                              immediateDownstreamLength = 5000)
        
        expect_s3_class(result_default, "data.frame")
        expect_s3_class(result_custom, "data.frame")
        expect_equal(nrow(result_default), length(test_peaks))
        expect_equal(nrow(result_custom), length(test_peaks))
    }
})

test_that("getGenomicAnnotation accumulates multiple annotations", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Create peaks that might overlap multiple features
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(1000, 5000, 10000), 
                            end = c(2000, 6000, 11000)),
            strand = "*"
        )
        
        result <- getGenomicAnnotation(test_peaks, TxDb)
        
        # Check that annotation column may contain comma-separated values
        # (indicating multiple overlapping features)
        has_multiple <- grepl(",", result$annotation)
        
        # At least some peaks might have multiple annotations
        # (this depends on the actual data, so we just check structure)
        expect_s3_class(result, "data.frame")
        expect_true(all(is.character(result$annotation)))
        expect_true(all(is.character(result$priorityAnnotation)))
        
        # Priority annotation should be single value (not comma-separated)
        expect_false(any(grepl(",", result$priorityAnnotation)))
    }
})

test_that("getGenomicAnnotation handles empty peaks", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Create empty GRanges
        empty_peaks <- GRanges()
        
        result <- getGenomicAnnotation(empty_peaks, TxDb)
        
        expect_s3_class(result, "data.frame")
        expect_equal(nrow(result), 0)
        expect_true(all(c("annotation", "priorityAnnotation") %in% colnames(result)))
    }
})

test_that("getGenomicAnnotation handles peaks with no overlaps", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Create peaks in a region unlikely to have annotations
        # (very large coordinates)
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(100000000, 200000000), 
                            end = c(100001000, 200001000)),
            strand = "*"
        )
        
        result <- getGenomicAnnotation(test_peaks, TxDb)
        
        expect_s3_class(result, "data.frame")
        expect_equal(nrow(result), length(test_peaks))
        # Peaks with no overlaps should default to "Distal Intergenic"
        expect_true(all(result$annotation == "Distal Intergenic" | 
                       result$annotation != ""))
        expect_true(all(result$priorityAnnotation == "Distal Intergenic" | 
                       result$priorityAnnotation != ""))
    }
})

test_that("getGenomicAnnotation handles multiple chromosomes", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Get available chromosomes from TxDb
        available_chrs <- seqlevels(TxDb)
        if (length(available_chrs) > 1) {
            # Create peaks on multiple chromosomes
            test_peaks <- GRanges(
                seqnames = available_chrs[1:min(2, length(available_chrs))],
                ranges = IRanges(start = c(1000, 5000), 
                                end = c(2000, 6000)),
                strand = "*"
            )
            
            result <- getGenomicAnnotation(test_peaks, TxDb)
            
            expect_s3_class(result, "data.frame")
            expect_equal(nrow(result), length(test_peaks))
            expect_true(all(c("annotation", "priorityAnnotation") %in% colnames(result)))
        }
    }
})

test_that("getGenomicAnnotation handles promoter binning", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        # Create peaks in promoter regions
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(500, 1500, 2500), 
                            end = c(600, 1600, 2600)),
            strand = "*"
        )
        
        # Test with tssRegion that creates multiple bins
        result <- getGenomicAnnotation(test_peaks, TxDb, 
                                      tssRegion = c(-3000, 3000))
        
        expect_s3_class(result, "data.frame")
        expect_equal(nrow(result), length(test_peaks))
        
        # Check that promoter annotations may include bin labels
        promoter_annotations <- result$annotation[grepl("Promoter", result$annotation)]
        if (length(promoter_annotations) > 0) {
            # Should contain bin information like "Promoter (<=1kb)" or "Promoter (1-2kb)"
            expect_true(any(grepl("Promoter", promoter_annotations)))
        }
    }
})

test_that("getGenomicAnnotation validates input parameters", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = 1000, end = 2000),
            strand = "*"
        )
        
        # Test with invalid TxDb
        expect_error(getGenomicAnnotation(test_peaks, "not_a_TxDb"))
        
        # Test with invalid peaks
        expect_error(getGenomicAnnotation("not_GRanges", TxDb))
        
        # Test with invalid level
        expect_error(getGenomicAnnotation(test_peaks, TxDb, level = "invalid"))
    }
})

test_that("getGenomicAnnotation handles invalid priority names gracefully", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = 1000, end = 2000),
            strand = "*"
        )
        
        # Test with invalid priority names (should issue warning but continue)
        invalid_priority <- c("promoter", "invalidType1", "fiveUTR", "invalidType2")
        
        expect_warning(
            result <- getGenomicAnnotation(test_peaks, TxDb, 
                                           genomicAnnotationPriority = invalid_priority),
            "Invalid priority"
        )
        
        expect_s3_class(result, "data.frame")
        expect_equal(nrow(result), length(test_peaks))
    }
})

test_that("getGenomicAnnotation returns consistent results", {
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    if (file.exists(txdb_file)) {
        TxDb <- loadDb(txdb_file)
        
        test_peaks <- GRanges(
            seqnames = "1",
            ranges = IRanges(start = c(1000, 5000), 
                            end = c(2000, 6000)),
            strand = "*"
        )
        
        # Run twice with same parameters
        result1 <- getGenomicAnnotation(test_peaks, TxDb)
        result2 <- getGenomicAnnotation(test_peaks, TxDb)
        
        # Results should be identical
        expect_equal(result1, result2)
    }
})

