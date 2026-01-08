#' Find occurrence of input motifs in promoter regions
#' 
#' @description 
#' Finds occurrences of input motifs in the promoter regions of specified genes.
#' This function identifies motif locations in promoter sequences extracted from
#' a BSgenome object using transcript annotations from a TxDb object. The function
#' can search for single motifs or paired motifs within a specified distance threshold.
#' 
#' The function extracts promoter regions (defined by \code{upstream} and 
#' \code{downstream} parameters relative to TSS) for the specified genes, then
#' searches for motif occurrences using \code{\link{summarizePatternInPeaks}}.
#' For paired motif searches, it uses \code{\link{annotatePeakInBatch}} to find
#' motif pairs that meet distance and orientation constraints.
#' 
#' @param patternFilePath1 Character string specifying the file path containing
#'        a list of known motifs in FASTA format (or other format specified by
#'        \code{format}). Required. Each sequence in the file represents a motif
#'        pattern to search for.
#' @param patternFilePath2 Character string specifying the file path containing
#'        motifs required to be in the flanking regions of the motif(s) in
#'        \code{patternFilePath1}. Required if \code{findPairedMotif = TRUE}.
#'        Used only for paired motif searches.
#' @param BSgenomeName A \code{BSgenome} object containing the reference genome
#'        sequences. For a list of available BSgenome objects, use
#'        \code{BSgenome::available.genomes()}. Examples:
#'        \itemize{
#'          \item \code{BSgenome.Hsapiens.UCSC.hg38} for human hg38
#'          \item \code{BSgenome.Hsapiens.UCSC.hg19} for human hg19
#'          \item \code{BSgenome.Mmusculus.UCSC.mm10} for mouse mm10
#'          \item \code{BSgenome.Celegans.UCSC.ce6} for C. elegans ce6
#'          \item \code{BSgenome.Rnorvegicus.UCSC.rn5} for rat rn5
#'          \item \code{BSgenome.Drerio.UCSC.danRer7} for zebrafish Zv9
#'          \item \code{BSgenome.Dmelanogaster.UCSC.dm3} for Drosophila dm3
#'        }
#'        Required.
#' @param findPairedMotif Logical. If \code{TRUE}, searches for motifs in paired
#'        configuration (requires \code{patternFilePath2}). If \code{FALSE}
#'        (default), searches only for motifs from \code{patternFilePath1}.
#' @param txdb A \code{TxDb} object containing transcript annotations. For
#'        creating and using TxDb objects, see the GenomicFeatures package. For
#'        a list of existing TxDb objects, search for annotation packages starting
#'        with "TxDb" at
#'        \url{http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData}.
#'        Examples:
#'        \itemize{
#'          \item \code{TxDb.Rnorvegicus.UCSC.rn5.refGene} for rat
#'          \item \code{TxDb.Mmusculus.UCSC.mm10.knownGene} for mouse
#'          \item \code{TxDb.Hsapiens.UCSC.hg19.knownGene} for human hg19
#'          \item \code{TxDb.Hsapiens.UCSC.hg38.knownGene} for human hg38
#'          \item \code{TxDb.Dmelanogaster.UCSC.dm3.ensGene} for Drosophila
#'          \item \code{TxDb.Celegans.UCSC.ce6.ensGene} for C. elegans
#'        }
#'        Required.
#' @param geneIDs Integer vector of one or more Entrez gene IDs. For example,
#'        the Entrez ID for EWSR1 is 2130 (see
#'        \url{https://www.genecards.org/cgi-bin/carddisp.pl?gene=EWSR1}).
#'        You can use the \code{\link{addGeneIDs}} function in ChIPpeakAnno to
#'        convert other types of gene IDs to Entrez IDs. Required.
#' @param upstream Integer. Number of base pairs upstream of the TSS to include
#'        in the promoter region. Default is \code{5000L}.
#' @param downstream Integer. Number of base pairs downstream of the TSS to
#'        include in the promoter region. Default is \code{5000L}.
#' @param name.motif1 Character string. Name of the motif from
#'        \code{patternFilePath1} for labeling output columns. Default is
#'        \code{"motif1"}. Used only when \code{findPairedMotif = TRUE}.
#' @param name.motif2 Character string. Name of the motif from
#'        \code{patternFilePath2} for labeling output columns. Default is
#'        \code{"motif2"}. Used only when \code{findPairedMotif = TRUE}.
#' @param max.distance Integer. Maximum allowed distance (in base pairs) between
#'        paired motifs for them to be included in the output. Distance is
#'        calculated based on \code{motif1LocForDistance} and
#'        \code{motif2LocForDistance}. Default is \code{100L}.
#' @param min.distance Integer. Minimum required distance (in base pairs) between
#'        paired motifs for them to be included in the output. Default is
#'        \code{1L}.
#' @param motif.orientation Character string specifying the required relative
#'        orientation between paired motifs. Options:
#'        \itemize{
#'          \item \code{"both"} (default): Any orientation is allowed
#'          \item \code{"motif1UpstreamOfMotif2"}: Motif1 must be located upstream
#'                of motif2 (based on genomic coordinates)
#'          \item \code{"motif2UpstreamOfMoif1"}: Motif2 must be located upstream
#'                of motif1 (based on genomic coordinates)
#'        }
#'        Only used when \code{findPairedMotif = TRUE}.
#' @param ignore.strand Logical. If \code{TRUE}, strand information is ignored
#'        when finding paired motifs (motifs can be on different strands). If
#'        \code{FALSE} (default), paired motifs must be on the same strand.
#'        Only used when \code{findPairedMotif = TRUE}.
#' @param format Character string. Format of the motif files specified in
#'        \code{patternFilePath1} and \code{patternFilePath2}. Default is
#'        \code{"fasta"}. Other formats may be supported depending on
#'        \code{\link{summarizePatternInPeaks}}.
#' @param skip Integer. Number of lines to skip at the beginning of the input
#'        motif files. Default is \code{0L}.
#' @param motif1LocForDistance Character string. Specifies which boundary of
#'        motif1 to use for calculating distance between paired motifs. Options:
#'        \code{"start"} or \code{"end"} (default). Only used when
#'        \code{findPairedMotif = TRUE}. This parameter is passed to
#'        \code{PeakLocForDistance} in \code{\link{annotatePeakInBatch}}.
#' @param motif2LocForDistance Character string. Specifies which boundary of
#'        motif2 to use for calculating distance between paired motifs. Options:
#'        \code{"start"} (default) or \code{"end"}. Only used when
#'        \code{findPairedMotif = TRUE}. This parameter is passed to
#'        \code{FeatureLocForDistance} in \code{\link{annotatePeakInBatch}}.
#' @param outfile Character string. Optional file path to save the search results
#'        as a tab-separated text file. If not provided, results are only returned
#'        as a GRanges object.
#' @param append Logical. If \code{TRUE}, results are appended to the file
#'        specified by \code{outfile}. If \code{FALSE} (default), the file is
#'        overwritten. Only used if \code{outfile} is provided.
#' 
#' @return Returns a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object
#'        containing motif occurrences. The structure depends on
#'        \code{findPairedMotif}:
#'        
#'        \strong{If \code{findPairedMotif = FALSE}:}
#'        Returns a GRanges object with one range per motif occurrence. Metadata
#'        columns include:
#'        \itemize{
#'          \item \code{seqnames}, \code{start}, \code{end}, \code{strand}:
#'                Genomic coordinates of the motif occurrence
#'          \item \code{tx_start}, \code{tx_end}, \code{tx_strand}: Promoter
#'                region boundaries and strand
#'          \item \code{motifName}: Name of the motif (from FASTA header)
#'          \item \code{motifPattern}: The motif sequence pattern
#'          \item \code{motifStartInPeak}, \code{motifEndInPeak}: Motif position
#'                relative to the promoter region start
#'          \item \code{motifFound}: The actual sequence found (may differ from
#'                pattern due to IUPAC ambiguity codes)
#'          \item \code{peakWidth}: Width of the promoter region
#'        }
#'        
#'        \strong{If \code{findPairedMotif = TRUE}:}
#'        Returns a GRanges object with one range per paired motif occurrence.
#'        The range represents the composite motif region (spanning from the
#'        start of the first motif to the end of the second motif). Metadata
#'        columns include:
#'        \itemize{
#'          \item \code{seqnames}, \code{start}, \code{end}, \code{strand}:
#'                Genomic coordinates of the composite motif region
#'          \item \code{compositMotifStart}, \code{compositMotifEnd}: Boundaries
#'                of the composite motif (same as start/end)
#'          \item Columns prefixed with \code{name.motif1}: Information about
#'                motif1 (start, end, width, strand, motifName, motifPattern,
#'                motifFound, etc.)
#'          \item Columns prefixed with \code{name.motif2}: Information about
#'                motif2 (start, end, width, strand, motifName, motifPattern,
#'                motifFound, etc.)
#'          \item \code{tx_start}, \code{tx_end}, \code{tx_strand}: Promoter
#'                region boundaries and strand
#'          \item \code{peakWidth}: Width of the promoter region
#'          \item \code{<name.motif1><name.motif2>relatviePositionTo}: Spatial
#'                relationship between motifs (e.g., "upstream", "downstream",
#'                "inside", "overlapStart", etc.)
#'          \item \code{<name.motif1><name.motif2>distanceto}: Distance between
#'                the motifs (based on \code{motif1LocForDistance} and
#'                \code{motif2LocForDistance})
#'          \item \code{shortestDistance}: Shortest distance between motif
#'                boundaries
#'        }
#' @author Lihua Julie Zhu, Kai Hu
#' @export
#' @importFrom GenomicFeatures transcriptsBy promoters
#' @importFrom S4Vectors mcols
#' @importFrom utils write.table
#' @importFrom matrixStats rowMins rowMaxs
#' @examples
#' 
#' 
#' library("BSgenome.Hsapiens.UCSC.hg38")
#' library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' 
#' patternFilePath1 =system.file("extdata", "motifIRF4.fa", package="ChIPpeakAnno")
#' patternFilePath2 =system.file("extdata", "motifAP1.fa", package="ChIPpeakAnno")
#' pairedMotifs <- findMotifsInPromoterSeqs(patternFilePath1 = patternFilePath1,
#'    patternFilePath2 = patternFilePath2,
#'    findPairedMotif = TRUE,
#'    name.motif1 = "IRF4", name.motif2 = "AP1",
#'    BSgenomeName = BSgenome.Hsapiens.UCSC.hg38,
#'    geneIDs = 7486, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'    outfile = "testPaired.xls")
#' 
#' unPairedMotifs <- findMotifsInPromoterSeqs(patternFilePath1 = patternFilePath1,
#'     BSgenomeName = BSgenome.Hsapiens.UCSC.hg38,
#'    geneIDs = 7486, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'    outfile = "testUnPaired.xls")
#' 
findMotifsInPromoterSeqs <-
  function(patternFilePath1,
           patternFilePath2,
           findPairedMotif = FALSE,
           BSgenomeName,
           txdb,
           geneIDs,
           upstream = 5000L,
           downstream = 5000L,
           name.motif1 = "motif1",
           name.motif2 = "motif2",
           max.distance = 100L, min.distance = 1L,
           motif.orientation = c("both", "motif1UpstreamOfMotif2", 
                                 "motif2UpstreamOfMoif1"),
           ignore.strand = FALSE,
           format = "fasta",
           skip = 0L,
           motif1LocForDistance = "end",
           motif2LocForDistance = "start",
           outfile, append = FALSE) {
    if (missing(patternFilePath1)) {
        stop("Missing required parameter 'patternFilePath1'!", call. = FALSE)
    }
    if (missing(txdb) || !inherits(txdb, "TxDb")) {
        stop("txdb is required as TxDb object!", call. = FALSE)
    }
    if (missing(geneIDs)) {
        stop("'geneIDs' is required as Entrez IDs", call. = FALSE)
    }
    motif.orientation <- match.arg(motif.orientation)
    tx <- transcriptsBy(txdb, by = "gene")
    tx <- tx[names(tx) %in% geneIDs]
    
    peaks <- promoters(tx, upstream = upstream, downstream = downstream)
    
    x1 <- do.call(rbind, 
                  lapply(seq_along(peaks),
                         function(i) {
                             thisPeak <- peaks[[i]]
                             mcols(thisPeak)$gene_id <- names(peaks)[i]
                             x <- summarizePatternInPeaks(patternFilePath = patternFilePath1, 
                                                          format = format,
                                                          skip = skip, 
                                                          BSgenomeName = BSgenomeName, 
                                                          peaks = thisPeak,
                                                          expectFrequencyMethod = "Naive")
                             x$motif_occurrence
                         }
                  )
    )
    
    # original colnames(x1):
      # 1: motifChr
      # 2: motifStartInChr
      # 3: motifEndInChr
      # 4: motifName
      # 5: motifPattern
      # 6: motifStartInPeak
      # 7: motifEndInPeak
      # 8: motifFound
      # 9: motifFoundStrand
      # 10: peakChr
      # 11: peakStart
      # 12: peakEnd
      # 13: peakWidth
      # 14: peakStrand
    colnames(x1)[1] <- "seqnames"
    colnames(x1)[2] <- "start"
    colnames(x1)[3] <- "end"
    colnames(x1)[9] <- "strand"
    colnames(x1)[11] <- "tx_start"
    colnames(x1)[12] <- "tx_end"
    colnames(x1)[14] <- "tx_strand"

    if (!findPairedMotif) {
        if (!missing(outfile))
            write.table(x1, file = outfile, sep = "\t", row.names = FALSE)
        toGRanges(x1)
    }
    else if (missing(patternFilePath2)) {
        stop("Missing required parameter 'patternFilePath2'!", call. = FALSE)
    }
    else {
        x2 <- do.call(rbind, 
                      lapply(seq_along(peaks), 
                             function(i) {
                                 thisPeak <- peaks[[i]]
                                 mcols(thisPeak)$gene_id <- names(peaks)[i]
                                 x <- summarizePatternInPeaks(patternFilePath = patternFilePath2,
                                                              format = format,
                                                              skip = skip, 
                                                              BSgenomeName = BSgenomeName, 
                                                              peaks = thisPeak,
                                                              expectFrequencyMethod = "Naive")
                                 x <- x$motif_occurrence
                             }
                      )
        )
  
        # original colnames(x2):
          # 1: motifChr
          # 2: motifStartInChr
          # 3: motifEndInChr
          # 4: motifName
          # 5: motifPattern
          # 6: motifStartInPeak
          # 7: motifEndInPeak
          # 8: motifFound
          # 9: motifFoundStrand
          # 10: peakChr
          # 11: peakStart
          # 12: peakEnd
          # 13: peakWidth
          # 14: peakStrand
      
        colnames(x2)[1] <- "seqnames"
        colnames(x2)[2] <- "start"
        colnames(x2)[3] <- "end"
        colnames(x2)[9] <- "strand"
        colnames(x2)[11] <- "tx_start"
        colnames(x2)[12] <- "tx_end"
        colnames(x2)[14] <- "tx_strand"

        x1.gr <- toGRanges(x1)
        x2.gr <- toGRanges(x2)
        
        res <- annotatePeakInBatch(x1.gr, AnnotationData = x2.gr,
                                   PeakLocForDistance = motif1LocForDistance,
                                   FeatureLocForDistance = motif2LocForDistance,
                                   ignore.strand = ignore.strand)
        
        temp <- res[res$shortestDistance <= max.distance & 
                    res$shortestDistance >= min.distance, ]
        
        if (motif.orientation == "motif1UpstreamOfMotif2") {
            temp <- temp[temp$insideFeature == "upstream", , drop = FALSE]
        }
        if (motif.orientation == "motif2UpstreamOfMoif1") {
            temp <- temp[temp$insideFeature == "downstream", , drop = FALSE]
        }
        
        temp <- as.data.frame(temp)
        
        m2 <- cbind(names(x2.gr), as.data.frame(x2.gr))
        colnames(m2)[1] <- "feature"
        m2 <- m2[, c(1:2, 7:8, 11)]
        colnames(m2)[3:5] <- paste(name.motif2, colnames(m2)[3:5])
        
        #m2$seqnames <- paste("chr", m2$seqnames, sep = "" )
        
        res <- merge(m2, temp, by = c("feature", "seqnames"))
        colnames(res)[6:14] <- paste(name.motif1, colnames(res)[6:14])
        # colnames(res):
          # 1: feature                  
          # 2: seqnames                 
          # 3: motif2_motifName           
          # 4: motif2_motifPattern         
          # 5: motif2_motifFound           
          # 6: motif1_start                   
          # 7: motif1_end                      
          # 8: motif1_width                    
          # 9: motif1_strand                  
          # 10: motif1_motifName                
          # 11: motif1_motifPattern             
          # 12: motif1_motifStartInPeak        
          # 13: motif1_motifEndInPeak           
          # 14: motif1_motifFound               
          # 15: peakChr                 
          # 16: tx_start                 
          # 17: tx_end                   
          # 18: peakWidth               
          # 19: tx_strand                
          # 20: peak                     
          # 21: start_position          
          # 22: end_position             
          # 23: feature_strand           
          # 24: insideFeature           
          # 25: distancetoFeature        
          # 26: shortestDistance         
          # 27: fromOverlappingOrNearest
        res <- res[, -c(1, 20)]
        
        colnames(res)[grep("insideFeature", colnames(res))] <- 
            paste(name.motif1, name.motif2, sep = "relativePositionTo")
        
        colnames(res)[19:20] <- paste(name.motif2, colnames(res)[19:20])
        colnames(res)[21] <- paste(name.motif2, "strand", sep = "_")
        
        
        tem <- cbind(as.numeric(res[, 5]), as.numeric(res[, 6]), 
                     as.numeric(res[, 19]), as.numeric(res[, 20]))
        res2 <- cbind(seqnames = res[, 1], compositMotifStart = rowMins(tem),
                      compositMotifEnd = rowMaxs(tem), res[, -1])
        
        res2 <- res2[, -grep("InPeak", colnames(res2))]
        
        d.ind <- grep("distancetoFeature", colnames(res2))
        res2 <- cbind(res2[c(1:6, d.ind, 7:ncol(res2)), drop = FALSE])
        res2 <- res2[, -(d.ind + 1)]
        colnames(res2)[grep("distancetoFeature", colnames(res2))] <- 
            paste(name.motif1, name.motif2, sep = "distanceto")
        colnames(res2)[grep("feature_strand", colnames(res2))] <- 
            paste(name.motif2, "strand")
        
        strand.ind <- grep("strand", colnames(res2))[1:3]
        strand <- res2[, strand.ind[1]]
        strand2 <- res2[, strand.ind[3]]
        strand <- ifelse(strand == strand2, as.character(strand), "*")
        res2 <- cbind(res2[, 1:3, drop = FALSE], strand = strand, 
                      res2[, 4:ncol(res2), drop = FALSE])
        if (!missing(outfile)) {
            write.table(res2, file = outfile, sep = "\t", row.names = FALSE)
        }
        colnames(res2)[2:3] <- c("start", "end")
        toGRanges(res2)
    }
  }
