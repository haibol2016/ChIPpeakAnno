#' Class \code{"bindist"}
#' 
#' @description 
#' An S4 class representing the distribution of binding sites relative to
#' genomic features (e.g., TSS or gene end). Objects of this class store
#' fixed-width bin counts and positions for use in permutation testing to
#' generate background pools that match the observed binding distribution.
#' 
#' 
#' @name bindist-class
#' @rdname bindist
#' @aliases bindist bindist-class bindist-method $,bindist-method
#' $<-,bindist-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("bindist", counts="integer", mids="integer",
#' halfBinSize="integer", bindingType="character", featureType="character")}.
#' @seealso \link{preparePool}, \link{peakPermTest}
#' @exportClass bindist
#' @exportMethod "$" "$<-"
#' @keywords classes

setClass("bindist", 
         representation(
             counts = "integer",
             mids = "integer",
             halfBinSize = "integer",
             bindingType = "character",
             featureType = "character"
         ),
         validity = function(object) {
             re <- TRUE
             if (is.null(object@counts)) {
                 re <- "counts is empty"
             }
             if (is.null(object@mids)) {
                 re <- "mids is empty"
             }
             if (is.null(object@halfBinSize)) {
                 re <- "halfBinSize is empty"
             }
             if (length(object@halfBinSize) > 1L) {
                 re <- "length of halfBinSize must be 1"
             }
             if (length(object@counts) != length(object@mids)) {
                 re <- "the length of mids and counts are not identical"
             }
             if (!object@bindingType %in% c("TSS", "geneEnd")) {
                 re <- "the bindingType must be TSS or geneEnd"
             }
             if (!object@featureType %in% c("transcript", "exon")) {
                 re <- "the featureType must be transcript or exon"
             }
             re
         }
)

setMethod("$", "bindist", function(x, name) slot(x, name))

setReplaceMethod("$", "bindist",
                 function(x, name, value) {
                     slot(x, name, check = TRUE) <- value
                     x
                 }
)
