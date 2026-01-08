#' Aggregate peaks over bins from feature sites
#' 
#' @description 
#' Aggregates peaks into bins relative to feature sites (e.g., TSS, gene ends)
#' and calculates summary statistics for each bin. The function supports multiple
#' peak sets, different binning strategies (distance-based or gene-body-based),
#' and flexible scoring functions. It automatically plots the distribution and
#' returns a data frame with bin values.
#' 
#' The function bins peaks based on their distance to feature sites (TSS or gene
#' ends) and calculates aggregated scores (e.g., counts, sums, means) for each bin.
#' When \code{aroundGene = TRUE}, it includes additional bins within the gene body,
#' allowing visualization of peak distribution across promoter, gene body, and
#' downstream regions.
#' 
#' @param \dots One or more \link[GenomicRanges:GRanges-class]{GRanges} objects
#'        containing peaks to be analyzed. Alternatively, a single
#'        \link[GenomicRanges:GRangesList-class]{GRangesList} can be provided.
#'        Each GRanges object should have a \code{score} metadata column (if missing,
#'        it will be set to 1 for all peaks). Multiple peak sets will be plotted
#'        in separate panels.
#' @param annotationData A \link[GenomicRanges:GRanges-class]{GRanges} or
#'        \link{annoGR} object containing annotation features (e.g., genes, transcripts).
#'        Must be non-empty. If an \code{annoGR} object is provided, it will be
#'        coerced to \code{GRanges}. Duplicate features are automatically removed.
#'        Strand information must be "+", "-", or "*".
#' @param select A character string specifying how to annotate peaks to features.
#'        Options:
#'        \itemize{
#'          \item \code{"all"}: Annotate peaks to all overlapping features (within
#'                \code{radius}). Uses \code{output = "overlapping"} in
#'                \code{annotatePeakInBatch}.
#'          \item \code{"nearest"} (default): Annotate peaks to the nearest feature
#'                only. Uses \code{output = "nearestLocation"} in
#'                \code{annotatePeakInBatch}.
#'        }
#' @param radius An integer specifying the maximum distance (in base pairs) from
#'        the feature site to include peaks. Peaks beyond this distance are excluded.
#'        Default is \code{5000L} (5kb). This parameter defines the upstream and
#'        downstream regions around the feature site.
#' @param nbins An integer specifying the number of bins for the upstream and
#'        downstream regions (each). Default is \code{50L}. Total number of bins
#'        is \code{2 * nbins} when \code{aroundGene = FALSE}, or \code{2 * nbins + mbins}
#'        when \code{aroundGene = TRUE}.
#' @param minGeneLen An integer specifying the minimum gene length (in base pairs)
#'        required for features to be included in the analysis. Features shorter
#'        than this are excluded. Default is \code{1L} (no filtering).
#' @param aroundGene A logical value. When \code{TRUE}, includes additional bins
#'        within the gene body (between feature start and end), creating a three-part
#'        visualization: upstream, gene body, and downstream. When \code{FALSE}
#'        (default), only bins upstream and downstream of the feature site are
#'        created. Gene body bins are normalized by gene length to account for
#'        variable gene sizes.
#' @param mbins An integer specifying the number of bins within the gene body when
#'        \code{aroundGene = TRUE}. Default is \code{nbins}. Gene body values are
#'        normalized by \code{(radius * mbins) / (genelen * nbins)} to account for
#'        variable gene lengths.
#' @param featureSite A character string specifying which feature site to use as
#'        the reference point for distance calculations. Options:
#'        \itemize{
#'          \item \code{"FeatureStart"} (default): Uses the feature start site
#'                (TSS for positive strand, gene end for negative strand)
#'          \item \code{"FeatureEnd"}: Uses the feature end site (gene end for
#'                positive strand, TSS for negative strand)
#'          \item \code{"bothEnd"}: Uses both feature start and end sites. Only
#'                supported when \code{PeakLocForDistance != "all"} and
#'                \code{aroundGene = FALSE}. Creates separate bins for peaks relative
#'                to feature start (upstream) and feature end (downstream).
#'        }
#' @param PeakLocForDistance A character string specifying which part of the peak
#'        to use for distance calculations. Options:
#'        \itemize{
#'          \item \code{"all"} (default): Uses both peak start and end. Peaks that
#'                span multiple bins are distributed across those bins. Not compatible
#'                with \code{featureSite = "bothEnd"}.
#'          \item \code{"start"}: Uses the peak start position
#'          \item \code{"middle"}: Uses the peak center (midpoint between start and end)
#'          \item \code{"end"}: Uses the peak end position
#'        }
#' @param FUN A function or list of functions to be used for aggregating scores
#'        within each bin. Default is \code{sum}. Common options include:
#'        \itemize{
#'          \item \code{sum}: Sum of scores in each bin (default)
#'          \item \code{mean}: Mean score per bin
#'          \item \code{median}: Median score per bin
#'          \item \code{length}: Count of peaks per bin (automatically converted to
#'                \code{sum} with scores set to 1)
#'        }
#'        If a list is provided, each function will be applied to each peak set,
#'        cycling through the list if there are more peak sets than functions.
#' @param errFun A function, list of functions, or numeric value(s) specifying how
#'        to calculate error bars. Default is \code{sd} (standard deviation).
#'        If a function, it will be applied to scores within each bin. If numeric,
#'        that value will be used as the error bar for all bins. If a list is
#'        provided, each element will be applied to each peak set, cycling through
#'        the list if needed.
#' @param xlab A character string or vector specifying x-axis labels. If missing,
#'        defaults to \code{"distance from <featureSite>"} when \code{aroundGene = FALSE},
#'        or \code{"Bins from <featureSite>"} when \code{aroundGene = TRUE}.
#' @param ylab A character string or vector specifying y-axis labels. If missing,
#'        defaults to \code{"Score"}.
#' @param main A character string or vector specifying plot titles. If missing,
#'        defaults to \code{"<peak_set_name> binding over <featureSite>"}.
#' 
#' @return Returns invisibly a data frame (or matrix) with bin values, where:
#'        \itemize{
#'          \item Rows represent bins (ordered from upstream to downstream)
#'          \item Columns represent different peak sets
#'          \item Values are the aggregated scores calculated by \code{FUN}
#'        }
#'        The function also plots the distribution automatically, showing:
#'        \itemize{
#'          \item X-axis: Distance from feature site (or bin numbers when
#'                \code{aroundGene = TRUE})
#'          \item Y-axis: Aggregated scores
#'          \item Error bars: Calculated using \code{errFun} (if not all zeros)
#'          \item Vertical lines: Marking feature boundaries when
#'                \code{aroundGene = TRUE}
#'        }
#'        Multiple peak sets are plotted in separate panels arranged in a grid.
#' @details
#' 
#' \strong{Algorithm Overview:}
#' 
#' The function performs the following steps:
#' \enumerate{
#'   \item Annotates peaks to features using \code{annotatePeakInBatch} with the
#'         specified \code{select} option
#'   \item Filters annotations to keep only those within \code{radius} and with
#'         gene length >= \code{minGeneLen}
#'   \item Calculates distances from peaks to feature sites based on
#'         \code{PeakLocForDistance} and \code{featureSite}
#'   \item Assigns peaks to bins based on their distances:
#'         \itemize{
#'           \item When \code{aroundGene = FALSE}: Creates \code{2 * nbins} bins
#'                 (nbins upstream, nbins downstream)
#'           \item When \code{aroundGene = TRUE}: Creates \code{2 * nbins + mbins}
#'                 bins (nbins upstream, mbins gene body, nbins downstream)
#'         }
#'   \item For gene body bins, normalizes scores by gene length to account for
#'         variable gene sizes
#'   \item Aggregates scores within each bin using \code{FUN}
#'   \item Calculates error bars using \code{errFun}
#'   \item Plots the distribution with error bars
#' }
#' 
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Strand-aware}: Automatically handles strand information, using
#'         TSS for positive strand and gene end for negative strand when
#'         \code{featureSite = "FeatureStart"}
#'   \item \strong{Gene body normalization}: When \code{aroundGene = TRUE}, gene
#'         body bins are normalized by gene length to prevent bias from variable
#'         gene sizes
#'   \item \strong{Peak spanning}: When \code{PeakLocForDistance = "all"}, peaks
#'         that span multiple bins are distributed across those bins
#'   \item \strong{Multiple peak sets}: Supports comparing multiple peak sets in
#'         separate panels
#' }
#' 
#' \strong{Common Use Cases:}
#' \itemize{
#'   \item \strong{Metagene plots}: Visualize peak distribution around TSS or gene
#'         ends with \code{aroundGene = TRUE}
#'   \item \strong{Promoter analysis}: Focus on promoter regions with
#'         \code{aroundGene = FALSE} and \code{featureSite = "FeatureStart"}
#'   \item \strong{3' end analysis}: Analyze peaks near gene ends with
#'         \code{featureSite = "FeatureEnd"}
#'   \item \strong{Peak counting}: Use \code{FUN = length} to count peaks per bin
#'   \item \strong{Signal intensity}: Use \code{FUN = sum} or \code{mean} with
#'         peak scores to analyze signal intensity
#' }
#' 
#' @author Jianhong Ou
#' @seealso \code{\link{annotatePeakInBatch}} for peak annotation,
#'          \code{\link{binOverGene}} for gene-level binning,
#'          \code{\link{featureAlignedDistribution}} for alternative distribution
#'          visualization
#' @keywords misc
#' @export
#' @import GenomicRanges
#' @importFrom BiocGenerics strand start end width score
#' @importFrom graphics segments axis abline
#' @examples
#' \dontrun{
#' ## Example 1: Basic usage - count peaks around TSS
#' bed <- system.file("extdata", "MACS_output.bed", package="ChIPpeakAnno")
#' gr1 <- toGRanges(bed, format="BED", header=FALSE)
#' data(TSS.human.GRCh37)
#' binOverFeature(gr1, annotationData = TSS.human.GRCh37,
#'                radius = 5000, nbins = 10, FUN = length, errFun = 0)
#' 
#' ## Example 2: Metagene plot with gene body
#' binOverFeature(gr1, annotationData = TSS.human.GRCh37,
#'                radius = 5000, nbins = 20, mbins = 20,
#'                aroundGene = TRUE, FUN = sum, errFun = sd)
#' 
#' ## Example 3: Multiple peak sets comparison
#' gr2 <- toGRanges(another_bed_file, format="BED", header=FALSE)
#' binOverFeature(gr1, gr2, annotationData = TSS.human.GRCh37,
#'                radius = 3000, nbins = 15, FUN = mean)
#' 
#' ## Example 4: Using peak scores (e.g., from MACS)
#' ## Assuming gr1 has a 'score' column from MACS
#' binOverFeature(gr1, annotationData = TSS.human.GRCh37,
#'                radius = 5000, nbins = 50, FUN = sum, errFun = sd)
#' 
#' ## Example 5: Gene end analysis
#' binOverFeature(gr1, annotationData = TSS.human.GRCh37,
#'                featureSite = "FeatureEnd", radius = 5000, nbins = 20)
#' }
#' 
binOverFeature <- function(..., annotationData = GRanges(),
                           select = c("all", "nearest"),
                           radius = 5000L, nbins = 50L,
                           minGeneLen = 1L, aroundGene = FALSE, mbins = nbins, 
                           featureSite = c("FeatureStart", "FeatureEnd", 
                                         "bothEnd"),
                           PeakLocForDistance = c("all", "end", 
                                                "start", "middle"), 
                           FUN = sum, errFun = sd, xlab, ylab, main) {
    ###check inputs
    PeaksList <- list(...)
    isGRangesList <- FALSE
    if (length(PeaksList) == 1) {
        if (inherits(PeaksList[[1]], "GRangesList")) {
            PeaksList <- PeaksList[[1]]
            names <- names(PeaksList)
            isGRangesList <- TRUE
        }
    }
    ##save dots arguments names
    if (!isGRangesList) {
        dots <- substitute(list(...))[-1]
        names <- unlist(sapply(dots, deparse))
    }
    
    n <- length(PeaksList)
    if (n == 0L) {
        stop("Missing required argument Peaks!", call. = FALSE)
    } else {
        nr <- ceiling(sqrt(n))
        nc <- ceiling(n / nr)
        op <- par(mfrow = c(nr, nc))
        on.exit(par(op))
    }
    if (missing(annotationData)) {
        stop("Missing required argument 'annotationData'. ",
             "Must be a GRanges or annoGR object.", call. = FALSE)
    }
    if (!inherits(annotationData, c("GRanges", "annoGR")) || 
           length(annotationData) < 1L) {
        stop("'annotationData' must be a non-empty GRanges or annoGR object.",
             call. = FALSE)
    }
    if (inherits(annotationData, "annoGR")) {
        annotationData <- as(annotationData, "GRanges")
    }
    annotationData <- unique(annotationData)
    if (!all(as.character(strand(annotationData)) %in% c("+", "-", "*"))) {
        stop("strands of annotationData must be +, - or *", call. = FALSE)
    }
    select <- match.arg(select)
    featureSite <- match.arg(featureSite)
    PeakLocForDistance <- match.arg(PeakLocForDistance)
    if(!is.list(FUN)) FUN <- list(FUN)
    if(!is.list(errFun)) errFun <- list(errFun)
    lapply(FUN, function(fun) {
        if (mode(fun) != "function") {
            stop("The mode of FUN must be function. ",
                 "The FUN could be any function such as ",
                 "median, mean, sum, length, ...", call. = FALSE)
        }
    })
    lapply(errFun, function(fun) {
        if (mode(fun) != "function" && !is.numeric(fun)) {
            stop("The mode of errFun must be function. ",
                 "The errFun could be any function such as sd", 
                 call. = FALSE)
        }
    })
    
    annotatedPeaksList <- lapply(PeaksList, function(Peaks) {
        if (!inherits(Peaks, "GRanges")) {
            stop("'Peaks' must be a GRanges object", call. = FALSE)
        }
        if (is.null(score(Peaks))) {
            message("score of GRanges object is required for calculation. ",
                    "It will be the input of FUN. Setting score as 1.")
            Peaks$score <- 1L
        }
        # Annotate the peaks
        if (select == "all") {
            annotatedPeaks <- 
                annotatePeakInBatch(Peaks, AnnotationData = annotationData,
                                    output = "overlapping", maxgap = radius, 
                                    select = "all") 
        } else {
            annotatedPeaks <- 
                annotatePeakInBatch(Peaks, AnnotationData = annotationData,
                                    output = "nearestLocation", select = "all")
        }
        # Filter the annotation
        annotatedPeaks <- annotatedPeaks[!is.na(annotatedPeaks$feature_strand)]
        annotatedPeaks <- 
            annotatedPeaks[
                as.numeric(as.character(annotatedPeaks$end_position)) -
                    as.numeric(as.character(annotatedPeaks$start_position)) + 1L >=
                    minGeneLen]
        
        ###if insideFeature==inside or includeFeature, 
        ###featureSite=="bothEnd.intergenic", the annotation should be removed
        #    if(featureSite=="bothEnd.intergenic") 
        #      annotatedPeaks <- 
        #        annotatedPeaks[!annotatedPeaks$insideFeature %in% 
        #                c("inside", "includeFeature"),]
        annotatedPeaks
    })
    
    plotErrBar <- function(x, y, err) {
        yplus <- y + err
        yminus <- y - err
        segments(x, yminus, x, yplus)
        xcoord <- par()$usr[1:2]
        smidge <- 0.015 * (xcoord[2] - xcoord[1]) / 2
        segments(x - smidge, yminus, x + smidge, yminus)
        segments(x - smidge, yplus, x + smidge, yplus)
    }
    
    if (missing(xlab)) {
        xlab <- if (aroundGene) {
            paste("Bins from", featureSite)
        } else {
            paste("distance from", featureSite)
        }
    }
    if (missing(ylab)) ylab <- "Score" 
    if (missing(main)) main <- paste(names, "binding over", featureSite)
    
    binValue <- mapply(function(annotatedPeaks, fun, errfun,
                                xlab.ele, ylab.ele, main.ele) {
        ##genelength
        genelen <- 
            annotatedPeaks$end_position - annotatedPeaks$start_position + 1
        ## change the function if it is length
        if (identical(fun, length)) {
            fun <- sum
            annotatedPeaks$score <- 1
        }
        ###step 1 calculate the distance,
        strand <- annotatedPeaks$feature_strand == "-"
        if (PeakLocForDistance == "all") {
            ##dist, the start distance to feature loc, dist2, 
            ##the end distance to feature loc
            PeakLoc.start <- start(annotatedPeaks)
            PeakLoc.end <- end(annotatedPeaks)
            if (featureSite == "bothEnd") {
                ##bothEnd, TODO
                stop("Cannot handle the combination of ",
                     "PeakLocForDistance=='all' AND featureSite=='bothEnd'",
                     call. = FALSE)
            } else {
                FeatureLoc <-
                    switch(featureSite,
                           FeatureStart = ifelse(strand, 
                                               annotatedPeaks$end_position, 
                                               annotatedPeaks$start_position),
                           FeatureEnd = ifelse(strand, 
                                             annotatedPeaks$start_position, 
                                             annotatedPeaks$end_position),
                           0)
                dist1 <- ifelse(strand, 
                                FeatureLoc - PeakLoc.end, 
                                PeakLoc.start - FeatureLoc)
                dist2 <- ifelse(strand, 
                                FeatureLoc - PeakLoc.start, 
                                PeakLoc.end - FeatureLoc)
                score <- score(annotatedPeaks)
                weight <- (radius * mbins) / (genelen * nbins)
                if (aroundGene) {
                    if (featureSite == "FeatureStart") {
                        ibin1 <- 
                            ifelse(dist1 < 0, 
                                   nbins + floor(dist1 * nbins / radius),
                                   ifelse(dist1 < genelen, 
                                          nbins + floor(dist1 * mbins / genelen), 
                                          nbins + mbins +
                                              floor((dist1 - genelen) *
                                                        nbins / radius)))
                        ibin2 <- 
                            ifelse(dist2 < 0, 
                                   nbins + floor(dist2 * nbins / radius),
                                   ifelse(dist2 < genelen, 
                                          nbins + floor(dist2 * mbins / genelen), 
                                          nbins + mbins +
                                              floor((dist2 - genelen) *
                                                        nbins / radius)))
                        b <- c(-1 * c(nbins:1), 0:(mbins + nbins - 1))
                        distance2 <- -1 * radius
                        distance1 <- genelen + radius
                    } else { ##featureSite == "FeatureEnd"
                        ibin1 <- 
                            ifelse(dist1 >= 0, 
                                   nbins + mbins + floor(dist1 * nbins / radius),
                                   ifelse(dist1 > -1 * genelen, 
                                          nbins + mbins +
                                              floor(dist1 * mbins / genelen), 
                                          nbins +
                                              floor((dist1 + genelen) *
                                                        nbins / radius)))
                        ibin2 <- 
                            ifelse(dist2 >= 0, 
                                   nbins + mbins + floor(dist2 * nbins / radius),
                                   ifelse(dist2 > -1 * genelen, 
                                          nbins + mbins +
                                              floor(dist2 * mbins / genelen), 
                                          nbins +
                                              floor((dist2 + genelen) *
                                                        nbins / radius)))
                        b <- c(-1 * c((mbins + nbins):1), 0:(nbins - 1))
                        distance2 <- -1 * (radius + genelen)
                        distance1 <- radius
                    }
                    numbins <- 2 * nbins + mbins
                    b.type <- c(rep(FALSE, nbins), 
                                rep(TRUE, mbins), 
                                rep(FALSE, nbins))
                } else {
                    ibin1 <- round(nbins + floor(dist1 * nbins / radius))
                    ibin2 <- round(nbins + floor(dist2 * nbins / radius))
                    b <- c(-1 * c(nbins:1), 0:(nbins - 1))
                    numbins <- 2 * nbins
                    distance2 <- -1 * radius
                    distance1 <- radius
                    genelen <- 0
                    minGeneLen <- -1
                    b.type <- rep(FALSE, 2 * nbins)
                }
                ids <- dist2 >= distance2 & 
                    dist1 < distance1 & 
                    genelen > minGeneLen
                score <- score[ids]
                weight <- weight[ids]
                ibin1 <- ibin1[ids]
                ibin2 <- ibin2[ids]
                ibin1[ibin1 < 0] <- 0
                ibin2[ibin2 > numbins] <- numbins
                gps <- lapply(0:(numbins - 1), function(.id) {
                    ##split the scores
                    .idx <- ibin1 <= .id & ibin2 >= .id
                    if (b.type[.id + 1]) {
                        score[.idx] * weight[.idx]
                    } else {
                        score[.idx]
                    }
                })
                names(gps) <- formatC(seq_len(numbins), 
                                      width = nchar(as.character(numbins)), 
                                      flag = "0")
            }
        } else { ##PeakLocForDistance %in% start, middle, end
            PeakLoc <- 
                switch(PeakLocForDistance,
                       middle = round(rowMeans(cbind(start(annotatedPeaks), 
                                                   end(annotatedPeaks)))),
                       start = start(annotatedPeaks),
                       end = end(annotatedPeaks),
                       0)
            if (featureSite == "bothEnd") { ##only consider outside of bothEnd
                FeatureStart <- ifelse(strand, 
                                    annotatedPeaks$end_position, 
                                    annotatedPeaks$start_position)
                FeatureEnd <- ifelse(strand, 
                                  annotatedPeaks$start_position, 
                                  annotatedPeaks$end_position)
                dist1 <- ifelse(strand, 
                                FeatureStart - PeakLoc, 
                                PeakLoc - FeatureStart) ##frome feature start
                dist2 <- ifelse(strand, 
                                FeatureEnd - PeakLoc, 
                                PeakLoc - FeatureEnd) ##from feature end
                score <- score(annotatedPeaks)
                if (aroundGene) {
                    stop("Cannot handle the combination of ",
                         "aroundGene==TRUE AND featureSite=='bothEnd'",
                         call. = FALSE)
                } else {
                    ibin1 <- round(nbins + floor(dist1 * nbins / radius))
                    ibin2 <- round(nbins + floor(dist2 * nbins / radius))
                    score1 <- score[ibin1 < nbins & ibin1 >= 0]
                    ibin1 <- ibin1[ibin1 < nbins & ibin1 >= 0]
                    score2 <- score[ibin2 >= nbins & ibin2 < 2 * nbins]
                    ibin2 <- ibin2[ibin2 >= nbins & ibin2 < 2 * nbins]
                    numbins <- 2 * nbins
                    b <- c(-1 * c(nbins:1), 0:(nbins - 1))
                    b.type <- rep(FALSE, length(b))
                    ##insert the empty bins
                    ibin1 <- formatC(ibin1, width = nchar(numbins), flag = "0")
                    ibin2 <- formatC(ibin2, width = nchar(numbins), flag = "0")
                    gps1 <- split(score1, ibin1)
                    gps2 <- split(score2, ibin2)
                    gps1 <- gps1[formatC(0:numbins, 
                                         width = nchar(numbins), 
                                         flag = "0")]
                    gps2 <- gps2[formatC(0:numbins, 
                                         width = nchar(numbins), 
                                         flag = "0")]
                    names(gps1) <- formatC(0:numbins, 
                                           width = nchar(numbins), 
                                           flag = "0")
                    names(gps2) <- formatC(0:numbins, 
                                           width = nchar(numbins), 
                                           flag = "0")
                    gps <- c(gps1[1:nbins], gps2[(nbins + 1):(2 * nbins)])
                }
            } else { ##featureSite %in% FeatureStart, FeatureEnd
                FeatureLoc <-
                    switch(featureSite,
                           FeatureStart = ifelse(strand, 
                                               annotatedPeaks$end_position, 
                                               annotatedPeaks$start_position),
                           FeatureEnd = ifelse(strand, 
                                             annotatedPeaks$start_position, 
                                             annotatedPeaks$end_position),
                           0)
                dist1 <- ifelse(strand, 
                                FeatureLoc - PeakLoc, 
                                PeakLoc - FeatureLoc)
                weight <- (radius * mbins) / (genelen * nbins)
                if (aroundGene) {
                    if (featureSite == "FeatureStart") {
                        ibin1 <- 
                            ifelse(dist1 < 0, 
                                   nbins + floor(dist1 * nbins / radius),
                                   ifelse(dist1 < genelen, 
                                          nbins + floor(dist1 * mbins / genelen), 
                                          nbins + mbins +
                                              floor((dist1 - genelen) *
                                                        nbins / radius)))
                        b <- c(-1 * c(nbins:1), 0:(mbins + nbins - 1))
                        distance1 <- radius + genelen
                        distance2 <- -1 * radius
                    } else { ##featureSite == "FeatureEnd"
                        ibin1 <- 
                            ifelse(dist1 >= 0, 
                                   nbins + mbins + floor(dist1 * nbins / radius),
                                   ifelse(dist1 > -1 * genelen, 
                                          nbins + mbins +
                                              floor(dist1 * mbins / genelen), 
                                          nbins +
                                              floor((dist1 + genelen) *
                                                        nbins / radius)))
                        b <- c(-1 * c((mbins + nbins):1), 0:(nbins - 1))
                        distance1 <- radius
                        distance2 <- -1 * (radius + genelen)
                    }
                    numbins <- 2 * nbins + mbins
                    b.type <- c(rep(FALSE, nbins), 
                                rep(TRUE, mbins), 
                                rep(FALSE, nbins))
                } else {
                    ibin1 <- round(nbins + floor(dist1 * nbins / radius))
                    numbins <- 2 * nbins
                    b <- c(-1 * c(nbins:1), 0:(nbins - 1))
                    genelen <- 0
                    minGeneLen <- -1
                    distance1 <- radius
                    distance2 <- -1 * radius
                    b.type <- rep(FALSE, 2 * nbins)
                }
                ids <- dist1 >= distance2 & 
                    dist1 < distance1 & 
                    genelen > minGeneLen
                score <- score(annotatedPeaks)[ids]
                weight <- weight[ids]
                ibin1 <- ibin1[ids]
                ibin1[ibin1 < 0] <- 0
                ibin1[ibin1 > numbins] <- numbins
                score[b.type[ibin1 + 1]] <- 
                    score[b.type[ibin1 + 1]] * weight[b.type[ibin1 + 1]]
                ibin1 <- formatC(ibin1, 
                                 width = nchar(numbins), 
                                 flag = "0")
                gps <- split(score, ibin1)
                gps <- gps[formatC(0:numbins, 
                                   width = nchar(numbins),
                                   flag = "0")]
                names(gps) <- formatC(0:numbins, 
                                      width = nchar(numbins), 
                                      flag = "0")
                ##insert the empty bins
                gps <- gps[1:numbins]
            }
        }
        gps <- lapply(gps, function(.ele) {
            if (is.null(.ele[1])) {
                0
            } else {
                .ele
            }
        })
        value <- unlist(lapply(gps, fun))
        std <- if (mode(errfun) == "function") unlist(lapply(gps, errfun)) else 
            rep(errfun, length(gps))
        std[is.na(std)] <- 0
        ##plot the figure
        ylim.min <- min(value[!is.na(value)] - std[!is.na(value)])
        ylim.max <- max(value[!is.na(value)] + std[!is.na(value)])
        ylim.dis <- (ylim.max - ylim.min) / 20
        blabel <- if (aroundGene) {
            c(seq.int(-radius, -radius / nbins, length.out = nbins) + radius / nbins / 2,
              1:mbins,
              seq.int(0, radius - radius / nbins, length.out = nbins) + radius / nbins / 2)
        } else {
            seq.int(-radius, radius - radius / nbins, length.out = 2 * nbins) +
                radius / nbins / 2
        }
        if (aroundGene) {
            plot(b, value, 
                 ylim = c(ylim.min - ylim.dis, ylim.max + ylim.dis),
                 xlab = xlab.ele, 
                 ylab = ylab.ele, 
                 main = main.ele,
                 xaxt = "n")
            b.type.at <- b[which(b.type)]
            b.type.at <- 
                b.type.at[c(1, length(b.type.at))] + c(-.5, .5)
            abline(v = b.type.at, lty = 2)
            b.at <- c(b[1] - .5, b[nbins] + .5, b[nbins + mbins] + .5, b[length(b)] + .5)
            b.label <- c(-1 * radius, "Feature Start", "Feature End", radius)
            axis(1, at = b.at, labels = b.label)
            if (!all(std == 0)) plotErrBar(b, value, std)
        } else {
            plot(blabel, value, 
                 ylim = c(ylim.min - ylim.dis, ylim.max + ylim.dis),
                 xlab = xlab.ele, 
                 ylab = ylab.ele, 
                 main = main.ele)
            if (!all(std == 0)) plotErrBar(blabel, value, std)
        }
        
        names(value) <- blabel
        value
    }, annotatedPeaksList, FUN, errFun, xlab, ylab, main)
    
    colnames(binValue) <- names
    ###output statistics
    return(invisible(binValue))
}
