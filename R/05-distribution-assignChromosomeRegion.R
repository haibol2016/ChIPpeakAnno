#import GenomicFeatures
## Jianhong Ou @ Mar.20, 2013


#' Summarize peak distribution over genomic regions
#' 
#' @description 
#' Categorizes peaks into genomic regions including exons, introns, proximal
#' promoters, immediate downstream regions, 5' UTRs, 3' UTRs, and intergenic
#' regions. The function calculates the percentage of peaks (or nucleotides) in
#' each category and computes Jaccard indices (intersection over union) for
#' measuring overlap significance. Supports both peak-centric and
#' nucleotide-centric views.
#' 
#' \strong{Two annotation modes:}
#' \itemize{
#'   \item \strong{Modern mode (recommended):} Uses \code{TxDb} or \code{EnsDb}
#'         objects to extract genomic features automatically. This is the
#'         preferred approach as it ensures consistency and handles all feature
#'         types correctly.
#'   \item \strong{Legacy mode:} Uses individual \code{GRanges} objects for
#'         exon, TSS, utr5, and utr3. This mode is maintained for backward
#'         compatibility but is deprecated. Users are encouraged to use
#'         \code{TxDb} instead.
#' }
#' 
#' @param peaks.RD A \link[GenomicRanges:GRanges-class]{GRanges} object
#'        containing peaks to be categorized. Metadata columns are ignored.
#' @param exon (Legacy mode only) A \code{GRanges} object containing exon
#'        annotations. Must have strand information. This parameter is for
#'        backward compatibility only. \code{TxDb} should be used instead.
#' @param TSS (Legacy mode only) A \code{GRanges} object containing transcription
#'        start site (TSS) annotations. Must have strand information. This
#'        parameter is for backward compatibility only. \code{TxDb} or
#'        \code{EnsDb} should be used instead.
#'        
#'        \strong{Note:} Pre-computed TSS datasets (e.g., \code{TSS.human.GRCh37})
#'        are provided for package examples and testing only. Users should
#'        generate their own annotations using \code{\link{getAnnotation}} or
#'        EnsDb/TxDb packages to match their genome assembly.
#' @param utr5 (Legacy mode only) A \code{GRanges} object containing 5' UTR
#'        annotations. Must have strand information. This parameter is for
#'        backward compatibility only. \code{TxDb} should be used instead.
#' @param utr3 (Legacy mode only) A \code{GRanges} object containing 3' UTR
#'        annotations. Must have strand information. This parameter is for
#'        backward compatibility only. \code{TxDb} should be used instead.
#' @param proximal.promoter.cutoff A named numeric vector with elements
#'        \code{"upstream"} and \code{"downstream"} specifying the boundaries
#'        (in base pairs) for proximal promoter classification. Default is
#'        \code{c(upstream = 2000, downstream = 100)}.
#'        
#'        Peaks that reside within \code{upstream} bases upstream from or
#'        overlap with the transcription start site (within \code{downstream}
#'        bases downstream) are classified as proximal promoters. Peaks that
#'        reside further upstream are classified as intergenic regions
#'        (enhancers).
#' @param immediate.downstream.cutoff A named numeric vector with elements
#'        \code{"upstream"} and \code{"downstream"} specifying the boundaries
#'        (in base pairs) for immediate downstream classification. Default is
#'        \code{c(upstream = 0, downstream = 1000)}.
#'        
#'        Peaks that reside within \code{downstream} bases downstream of gene
#'        end (but not overlapping 3' UTR) are classified as immediate
#'        downstream. Peaks that reside further downstream are classified as
#'        intergenic regions (enhancers).
#' @param nucleotideLevel A logical value. When \code{FALSE} (default),
#'        calculates peak-centric percentages (percentage of peaks in each
#'        category). When \code{TRUE}, calculates nucleotide-centric percentages
#'        (percentage of nucleotides covered by peaks in each category) and uses
#'        disjoined regions to handle overlapping annotations properly.
#' @param precedence A character vector specifying the precedence order for
#'        assigning peaks to categories when overlaps occur. Valid values are:
#'        \code{"Promoters"}, \code{"immediateDownstream"}, \code{"fiveUTRs"},
#'        \code{"threeUTRs"}, \code{"Exons"}, and \code{"Introns"}. Default is
#'        \code{NULL}.
#'        
#'        \strong{Behavior:}
#'        \itemize{
#'          \item If \code{NULL} (default): Double counting is enabled. A peak
#'                overlapping both promoter and 5' UTR will be counted in both
#'                categories.
#'          \item If specified: Peaks are assigned to the first matching category
#'                in the precedence order. For example, if \code{precedence =
#'                c("Promoters", "fiveUTRs", ...)}, a peak overlapping both
#'                promoter and 5' UTR will only be counted in the "Promoters"
#'                category.
#'        }
#'        
#'        \strong{Note:} Precedence only applies when \code{nucleotideLevel =
#'        FALSE}. When \code{nucleotideLevel = TRUE}, overlapping regions are
#'        disjoined and each disjoined segment is assigned to the appropriate
#'        category based on precedence (if specified).
#' @param TxDb An object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}} or
#'        \code{\link[ensembldb:EnsDb-class]{EnsDb}} containing genomic
#'        annotation. When provided, the function uses this to extract all
#'        required features (exons, introns, UTRs, transcripts) automatically.
#'        This is the recommended approach. If \code{NULL}, the function falls
#'        back to legacy mode using individual \code{GRanges} parameters.
#' @return 
#' 
#' \strong{Modern Mode (when \code{TxDb} is provided):}
#' 
#' Returns a list with two named vectors:
#' \itemize{
#'   \item \code{percentage}: A named numeric vector containing the percentage
#'         of peaks (or nucleotides if \code{nucleotideLevel = TRUE}) in each
#'         genomic region category:
#'         \itemize{
#'           \item \code{"Exons"}: Percentage in exon regions
#'           \item \code{"Introns"}: Percentage in intron regions
#'           \item \code{"fiveUTRs"}: Percentage in 5' UTR regions
#'           \item \code{"threeUTRs"}: Percentage in 3' UTR regions
#'           \item \code{"Promoters"}: Percentage in proximal promoter regions
#'                 (within \code{proximal.promoter.cutoff} of TSS)
#'           \item \code{"immediateDownstream"}: Percentage in immediate
#'                 downstream regions (within \code{immediate.downstream.cutoff}
#'                 of gene end)
#'           \item \code{"Intergenic.Region"}: Percentage in intergenic regions
#'                 (enhancers and other non-genic regions)
#'         }
#'   \item \code{jaccard}: A named numeric vector containing Jaccard indices
#'         (intersection over union) for each category. Values range from 0 to 1,
#'         where higher values indicate more significant overlap between peaks
#'         and genomic features.
#' }
#' 
#' \strong{Legacy Mode (when \code{TxDb} is \code{NULL}):}
#' 
#' Returns a named list of percentages only (no Jaccard indices):
#' \itemize{
#'   \item \code{"Exons"}, \code{"Introns"}, \code{"fiveUTRs"},
#'         \code{"threeUTRs"}, \code{"Promoters"}: Same as modern mode
#'   \item \code{"immediate.Downstream"}: Note the dot in the name (different
#'         from modern mode's \code{"immediateDownstream"})
#'   \item \code{"Intergenic.Region"}: Same as modern mode
#' }
#' 
#' \strong{Jaccard Index Calculation (Modern Mode only):}
#' \itemize{
#'   \item \code{nucleotideLevel = FALSE}: Jaccard index = (number of
#'         overlapping peaks) / (total peaks + total features - overlapping
#'         peaks)
#'   \item \code{nucleotideLevel = TRUE}: Jaccard index = (overlapping
#'         nucleotides) / (total nucleotides in union of peaks and features)
#' }
#' @details
#' 
#' \strong{Modern Mode (TxDb/EnsDb):}
#' 
#' When \code{TxDb} or \code{EnsDb} is provided, the function:
#' \itemize{
#'   \item Extracts genomic features automatically: exons, introns, 5' UTRs,
#'         3' UTRs, transcripts, and optionally tRNAs (if available)
#'   \item Creates promoter regions using \code{promoters()} with
#'         \code{proximal.promoter.cutoff}
#'   \item Creates immediate downstream regions using \code{downstreams()} with
#'         \code{immediate.downstream.cutoff}
#'   \item Calculates intergenic regions as gaps between all annotated features
#'   \item Handles sequence level matching: peaks on chromosomes not present in
#'         the annotation are filtered out with a warning
#'   \item Trims ranges to chromosome boundaries using \code{trim()}
#' }
#' 
#' \strong{Legacy Mode:}
#' 
#' When \code{TxDb} is \code{NULL}, the function uses individual \code{GRanges}
#' objects. This mode:
#' \itemize{
#'   \item Uses a sequential annotation approach: first TSS, then UTRs, then exons
#'   \item Calculates percentages based on peak counts only (not nucleotide-level)
#'   \item Does not calculate Jaccard indices
#'   \item Returns a different structure (named list of percentages only)
#' }
#' 
#' \strong{Strand Handling:}
#' 
#' \itemize{
#'   \item If all peaks have strand "*" (unstranded), strand information is
#'         ignored in overlap calculations and intergenic region calculation
#'   \item Otherwise, strand-aware overlaps are performed
#' }
#' 
#' \strong{Precedence and Double Counting:}
#' 
#' When \code{precedence = NULL}:
#' \itemize{
#'   \item Peaks can be assigned to multiple categories if they overlap multiple
#'         feature types
#'   \item Percentages may sum to more than 100%
#' }
#' 
#' When \code{precedence} is specified:
#' \itemize{
#'   \item Peaks are assigned to the first matching category in the precedence
#'         order
#'   \item Percentages sum to 100% (each peak counted exactly once)
#'   \item In nucleotide-level mode, disjoined regions are assigned based on
#'         precedence
#' }
#' 
#' @author Jianhong Ou, Lihua Julie Zhu
#' @seealso \link{genomicElementDistribution}, \link{genomicElementUpSetR},
#' \link{binOverFeature}, \link{binOverGene}, \link{binOverRegions}
#' @references 1. Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package
#' to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237doi:10.1186/1471-2105-11-237
#' 
#' 2. Zhu L.J. (2013) Integrative analysis of ChIP-chip and ChIP-seq dataset.
#' Methods Mol Biol. 2013;1067:105-24. doi: 10.1007/978-1-62703-607-8_8.
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels
#' @importFrom BiocGenerics start end width strand
#' @importFrom GenomicFeatures exons intronsByTranscript fiveUTRsByTranscript 
#' threeUTRsByTranscript transcripts tRNAs
#' @examples
#' 
#' if (interactive() || Sys.getenv("USER")=="jou"){
#'     ##Display the list of genomes available at UCSC:
#'     #library(rtracklayer)
#'     #ucscGenomes()[, "db"]
#'     ## Display the list of Tracks supported by makeTxDbFromUCSC()
#'     #supportedUCSCtables()
#'     ##Retrieving a full transcript dataset for Human from UCSC
#'     ##TranscriptDb <- 
#'     ##     makeTxDbFromUCSC(genome="hg19", tablename="ensGene")
#'     if(require(TxDb.Hsapiens.UCSC.hg19.knownGene)){
#'       TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'       exons <- exons(TxDb, columns=NULL)
#'       fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
#'       Feature.distribution <- 
#'           assignChromosomeRegion(exons, nucleotideLevel=TRUE, TxDb=TxDb)
#'       barplot(Feature.distribution$percentage)
#'       assignChromosomeRegion(fiveUTRs, nucleotideLevel=FALSE, TxDb=TxDb)
#'       data(myPeakList)
#'       assignChromosomeRegion(myPeakList, nucleotideLevel=TRUE, 
#'                              precedence=c("Promoters", "immediateDownstream", 
#'                                           "fiveUTRs", "threeUTRs", 
#'                                           "Exons", "Introns"), 
#'                              TxDb=TxDb)
#'     }
#' }
#' 
assignChromosomeRegion <-
    function(peaks.RD, exon, TSS, utr5, utr3, 
             proximal.promoter.cutoff = c(upstream = 2000, downstream = 100), 
             immediate.downstream.cutoff = c(upstream = 0, downstream = 1000), 
             nucleotideLevel = FALSE, 
             precedence = NULL, TxDb = NULL) {
        ##check inputs
        if (!is.null(TxDb)) {
            if (!inherits(TxDb, c("TxDb", "EnsDb"))) 
                stop("TxDb must be an object of TxDb or similar such as EnsDb, 
                     try\n?TxDb\tto see more info.")
            if (!inherits(peaks.RD, "GRanges")) {
                stop("peaks.RD must be a GRanges object", call. = FALSE)
            }
            if (!all(c("upstream", "downstream") %in%
                    names(proximal.promoter.cutoff))) {
                stop("proximal.promoter.cutoff must contain elements ",
                     "'upstream' and 'downstream'", call. = FALSE)
            }
            if (!all(c("upstream", "downstream") %in%
                    names(immediate.downstream.cutoff))) {
                stop("immediate.downstream.cutoff must contain elements ",
                     "'upstream' and 'downstream'", call. = FALSE)
            }
            if (!is.null(precedence)) {
                if (!all(precedence %in% c("Exons", "Introns", "fiveUTRs", 
                                           "threeUTRs", "Promoters", 
                                           "immediateDownstream"))) 
                    stop("precedence must be a combination of: ",
                         "Exons, Introns, fiveUTRs, threeUTRs, ",
                         "Promoters, immediateDownstream", call. = FALSE)
            }
            ignore.strand <- all(as.character(strand(peaks.RD)) == "*")
            exons <- exons(TxDb, columns=NULL)
            introns <- unique(unlist(intronsByTranscript(TxDb)))
            fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
            threeUTRs <- unique(unlist(threeUTRsByTranscript(TxDb)))
            transcripts <- unique(transcripts(TxDb, columns=NULL))
            options(warn = -1)
            try({
                promoters <- 
                    unique(promoters(TxDb, upstream = proximal.promoter.cutoff["upstream"], 
                                     downstream = proximal.promoter.cutoff["downstream"]))
                immediateDownstream <- 
                    unique(downstreams(transcripts, 
                                 upstream = immediate.downstream.cutoff["upstream"], 
                                 downstream = immediate.downstream.cutoff["downstream"]))
                promoters <- GenomicRanges::trim(promoters)
                immediateDownstream <- GenomicRanges::trim(immediateDownstream)
            })
            # microRNAs <- tryCatch(microRNAs(TxDb), 
            #                       error = function(e) return(NULL))
            tRNAs <- tryCatch(tRNAs(TxDb), error = function(e) return(NULL))
            options(warn = 0)
            annotation <- list(exons, introns, fiveUTRs, threeUTRs, 
                               promoters, immediateDownstream)
            # if (!is.null(microRNAs)) 
            #     annotation <- c(annotation, "microRNAs" = microRNAs)
            if (!is.null(tRNAs)) 
                annotation <- c(annotation, "tRNAs" = tRNAs)
            annotation <- 
                lapply(annotation, function(.anno){mcols(.anno)<-NULL; .anno})
            names(annotation)[1:6] <- 
                c("Exons", "Introns", "fiveUTRs", "threeUTRs", 
                  "Promoters", "immediateDownstream")
            ###clear seqnames, the format should be chr+NUM
            peaks.RD <- formatSeqnames(peaks.RD, exons)
            peaks.RD <- unique(peaks.RD)
            annotation <- GRangesList(annotation)
            newAnno <- c(unlist(annotation))
            if (ignore.strand) {
                newAnno.rd <- newAnno
                strand(newAnno.rd) <- "*"
                newAnno.rd <- reduce(trim(newAnno.rd))
                Intergenic.Region <- gaps(newAnno.rd, end = seqlengths(TxDb))
                Intergenic.Region <- 
                    Intergenic.Region[strand(Intergenic.Region) == "*"]
            } else {
                newAnno.rd <- reduce(trim(newAnno))
                Intergenic.Region <- gaps(newAnno.rd, end = seqlengths(TxDb))
                Intergenic.Region <- 
                    Intergenic.Region[strand(Intergenic.Region) != "*"]
            }
            if (!all(seqlevels(peaks.RD) %in% seqlevels(newAnno))) {
                warning("peaks.RD has sequence levels not in TxDb.",
                        call. = FALSE, immediate. = TRUE)
                sharedlevels <- 
                    intersect(seqlevels(newAnno), seqlevels(peaks.RD))
                peaks.RD <- keepSeqlevels(peaks.RD, sharedlevels, 
                                          pruning.mode = "coarse")
            }
            mcols(peaks.RD) <- NULL
            if (!is.null(precedence)) {
                annotation <- 
                    annotation[unique(c(precedence, names(annotation)))]
            }
            ##    annotation$Intergenic.Region <- peaks.RD
            names(Intergenic.Region) <- NULL
            annotation$Intergenic.Region <- Intergenic.Region
            anno.names <- names(annotation)
            ol.anno <- findOverlaps(peaks.RD, annotation,
                                    ignore.strand = ignore.strand)
            if (nucleotideLevel) {
                ## calculate Jaccard index
                jaccardIndex <- unlist(lapply(annotation, function(.ele) {
                    intersection <- intersect(.ele, peaks.RD, 
                                              ignore.strand = ignore.strand)
                    union <- union(.ele, peaks.RD, ignore.strand = ignore.strand)
                    sum(as.numeric(width(intersection))) /
                        sum(as.numeric(width(union)))
                }))
                jaccardIndex <- jaccardIndex[anno.names]
                names(jaccardIndex) <- anno.names
                jaccardIndex[is.na(jaccardIndex)] <- 0
                
                ## create a new annotations
                newAnno <- unlist(annotation)
                newAnno$source <- rep(names(annotation), lengths(annotation))
                newAnno.disjoin <- disjoin(newAnno, with.revmap = TRUE, 
                                           ignore.strand = ignore.strand)
                if (!is.null(precedence)) {
                    revmap <- cbind(from = unlist(newAnno.disjoin$revmap), 
                                    to = rep(seq_along(newAnno.disjoin), 
                                           lengths(newAnno.disjoin$revmap)))
                    revmap <- revmap[order(revmap[, "to"], revmap[, "from"]), , drop = FALSE]
                    revmap <- revmap[!duplicated(revmap[, "to"]), , drop = FALSE]
                    newAnno.disjoin$source <- newAnno[revmap[, "from"]]$source
                } else {
                    revmap <- unlist(newAnno.disjoin$revmap)
                    newAnno.disjoin <- rep(newAnno.disjoin, lengths(newAnno.disjoin$revmap))
                    newAnno.disjoin$source <- newAnno[revmap]$source
                }
                newAnno.disjoin$revmap <- NULL
                ol.anno <- findOverlaps(peaks.RD, newAnno.disjoin, ignore.strand = ignore.strand)
                queryHits <- peaks.RD[queryHits(ol.anno)]
                subjectHits <- newAnno.disjoin[subjectHits(ol.anno)]
                totalLen <- sum(as.numeric(width(peaks.RD)))
                queryHits.list <- split(queryHits, subjectHits$source)
                lens <- unlist(lapply(queryHits.list, function(.ele) 
                    sum(as.numeric(width(unique(.ele))))))
                percentage <- 100 * lens / totalLen
            } else {
                ##calculate Jaccard index
                ol.anno.splited <- split(queryHits(ol.anno),
                                         anno.names[subjectHits(ol.anno)])
                jaccardIndex <- unlist(lapply(anno.names, function(.name) {
                    union <- length(annotation[[.name]]) + 
                        length(peaks.RD) - 
                        length(unique(subjectHits(findOverlaps(peaks.RD, 
                                                               annotation[[.name]], 
                                                               ignore.strand = ignore.strand))))
                    intersection <- length(ol.anno.splited[[.name]])
                    intersection / union
                }))
                names(jaccardIndex) <- anno.names
                ol.anno <- as.data.frame(ol.anno)
                ####keep the part only annotated in peaks.RD for peaks.RD
                ol.anno.splited <- split(ol.anno, ol.anno[, 2])
                hasAnnoHits <- 
                    do.call(rbind, 
                            ol.anno.splited[names(ol.anno.splited) !=
                                            as.character(length(annotation))])
                hasAnnoHits <- unique(hasAnnoHits[, 1])
                ol.anno <- 
                    ol.anno[!(ol.anno[, 2] == length(annotation) & 
                              (ol.anno[, 1] %in% hasAnnoHits)), ]    
                if (!is.null(precedence)) {
                    ol.anno <- ol.anno[!duplicated(ol.anno[, 1]), ]
                }
                ##calculate percentage
                subjectHits <- anno.names[ol.anno[, 2]]
                counts <- table(subjectHits)
                percentage <- 100 * counts / length(peaks.RD)
            }
            len <- length(anno.names) - length(percentage)
            if (len > 0) {
                tobeadd <- rep(0, len)
                names(tobeadd) <- anno.names[!anno.names %in% 
                                                 names(percentage)]
                percentage <- c(percentage, tobeadd)
            }
            percentage <- percentage[anno.names]
            return(list(percentage = percentage, jaccard = jaccardIndex))
        } else {
            message("Please try to use TxDb next time. Try ",
                    "?TxDb to see more info.")
            annotationList <- list(exon, TSS, utr5, utr3)
            names(annotationList) <- c("Exon", "TSS", "UTR5", "UTR3")
            lapply(annotationList, function(.ele) {
                if (!inherits(.ele, "GRanges")) {
                    stop("Annotation of exon, TSS, utr5, utr3 must ",
                         "be objects of GRanges.", call. = FALSE)
                }
            })
            if (!inherits(peaks.RD, "GRanges")) {
                stop("peaks.RD must be a GRanges object.", call. = FALSE)
            } 
            ann.peaks <- annotatePeakInBatch(peaks.RD, AnnotationData = TSS)
            ann.peaks <- ann.peaks[!is.na(ann.peaks$distancetoFeature)]
            upstream <- 
                ann.peaks[ann.peaks$insideFeature == "upstream" | 
                          (ann.peaks$distancetoFeature < 0 & 
                           ann.peaks$insideFeature == "overlapStart" & 
                           abs(ann.peaks$distancetoFeature) >
                           ann.peaks$shortestDistance) | 
                          ann.peaks$insideFeature == "includeFeature" | 
                          (ann.peaks$distancetoFeature >= 0 & 
                           ann.peaks$insideFeature == "overlapStart" & 
                           ann.peaks$distancetoFeature ==
                           ann.peaks$shortestDistance)]
            
            proximal.promoter.n <- 
                length(upstream[upstream$distancetoFeature >= 
                                -proximal.promoter.cutoff | 
                                upstream$shortestDistance <= 
                                proximal.promoter.cutoff])
            enhancer.n <- length(upstream) - proximal.promoter.n
            
            downstream <- ann.peaks[ann.peaks$insideFeature == "downstream"]
            immediateDownstream.n <- 
                length(downstream[downstream$distancetoFeature <= 
                                  immediate.downstream.cutoff, ])
            enhancer.n <- enhancer.n + 
                nrow(downstream[downstream$distancetoFeature > 
                                immediate.downstream.cutoff, , drop = FALSE])
            
            inside.peaks <- 
                ann.peaks[ann.peaks$insideFeature == "inside" | 
                          ann.peaks$insideFeature ==
                          "overlapEnd" |  
                          (ann.peaks$insideFeature == "overlapStart" & 
                           ann.peaks$distancetoFeature >= 0 & 
                           ann.peaks$distancetoFeature != 
                           ann.peaks$shortestDistance) | 
                          (ann.peaks$insideFeature == "overlapStart" & 
                           ann.peaks$distancetoFeature < 0 & 
                           abs(ann.peaks$distancetoFeature) ==
                           ann.peaks$shortestDistance)]
            
            ann.utr5.peaks <- annotatePeakInBatch(inside.peaks, 
                                                  AnnotationData = utr5)
            
            proximal.promoter.n <- proximal.promoter.n + 
                length(ann.utr5.peaks[ann.utr5.peaks$insideFeature ==
                                     "upstream"])
            
            utr5.n <- length(
                ann.utr5.peaks[ann.utr5.peaks$insideFeature %in% 
                               c("includeFeature", "inside") | 
                               (ann.utr5.peaks$insideFeature == "overlapStart" & 
                                ann.utr5.peaks$distancetoFeature >= 0 & 
                                ann.utr5.peaks$distancetoFeature != 
                                ann.utr5.peaks$shortestDistance) | 
                               (ann.utr5.peaks$insideFeature == "overlapStart" & 
                                ann.utr5.peaks$distancetoFeature < 0 & 
                                abs(ann.utr5.peaks$distancetoFeature) ==
                                ann.utr5.peaks$shortestDistance) | 
                               (ann.utr5.peaks$insideFeature == "overlapEnd" & 
                                ann.utr5.peaks$strand == "+" & 
                                abs(start(ann.utr5.peaks) -
                                    ann.utr5.peaks$end_position) >= 
                                (end(ann.utr5.peaks) -
                                 ann.utr5.peaks$end_position)) | 
                               (ann.utr5.peaks$insideFeature == "overlapEnd" & 
                                ann.utr5.peaks$strand == "-" & 
                                abs(end(ann.utr5.peaks) -
                                    ann.utr5.peaks$start_position) >= 
                                abs(start(ann.utr5.peaks) -
                                    ann.utr5.peaks$start_position))])
            
            proximal.promoter.n <- 
                proximal.promoter.n +  
                length(
                    ann.utr5.peaks[
                        (ann.utr5.peaks$insideFeature == "overlapStart" & 
                         ann.utr5.peaks$distancetoFeature >= 0 & 
                         ann.utr5.peaks$distancetoFeature == 
                         ann.utr5.peaks$shortestDistance) | 
                            (ann.utr5.peaks$insideFeature == "overlapStart" & 
                             ann.utr5.peaks$distancetoFeature < 0 & 
                             abs(ann.utr5.peaks$distancetoFeature) !=
                             ann.utr5.peaks$shortestDistance)])
            
            downstream.utr5 <-
                ann.utr5.peaks[
                    ann.utr5.peaks$insideFeature == "downstream" |
                        (ann.utr5.peaks$insideFeature == "overlapEnd" & 
                         ann.utr5.peaks$strand == "+" & 
                         abs(start(ann.utr5.peaks) -
                             ann.utr5.peaks$end_position) < 
                             (end(ann.utr5.peaks) -
                              ann.utr5.peaks$end_position)) | 
                        (ann.utr5.peaks$insideFeature == "overlapEnd" & 
                         ann.utr5.peaks$strand == "-" & 
                         abs(end(ann.utr5.peaks) -
                             ann.utr5.peaks$start_position) < 
                             abs(start(ann.utr5.peaks) -
                                 ann.utr5.peaks$start_position))] 
            
            ann.utr3.peaks <- annotatePeakInBatch(downstream.utr5, 
                                                  AnnotationData = utr3)
            
            utr3.n <- 
                length(ann.utr3.peaks[ann.utr3.peaks$insideFeature %in% 
                                      c("includeFeature", "overlapStart", 
                                        "overlapEnd", "inside")])
            
            rest.peaks <- ann.utr3.peaks[ann.utr3.peaks$insideFeature %in% 
                                         c("downstream", "upstream")]
            
            ann.rest.peaks <- annotatePeakInBatch(rest.peaks, 
                                                 AnnotationData = exon)
            
            intron.n <- length(ann.rest.peaks[ann.rest.peaks$insideFeature %in%
                                              c("downstream", "upstream")])
            exon.n <- length(ann.rest.peaks) - intron.n
            
            total <- length(peaks.RD) / 100
            
            list("Exons" = exon.n / total, 
                 "Introns" = intron.n / total, 
                 "fiveUTRs" = utr5.n / total, 
                 "threeUTRs" = utr3.n / total, 
                 "Promoters" = proximal.promoter.n / total, 
                 "immediate.Downstream" = immediateDownstream.n / total, 
                 "Intergenic.Region" = enhancer.n / total)
        }
    }
