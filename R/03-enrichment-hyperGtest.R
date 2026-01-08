#' Perform hypergeometric test for enrichment analysis
#' 
#' @description 
#' Performs a hypergeometric test (one-tailed, upper tail) to assess enrichment
#' of terms (e.g., GO terms, pathway terms) in a dataset compared to a background
#' genome. This function is used internally by \code{\link{getEnrichedGO}} and
#' \code{\link{getEnrichedPATH}}.
#' 
#' The test calculates P(X >= q) where X follows a hypergeometric distribution
#' with parameters: m (number of white balls = count of term in genome), 
#' n-m (number of black balls = remaining terms in genome), and
#' k (number of balls drawn = total terms in peak list). The test evaluates
#' whether the observed count of a term in the peak list is significantly
#' higher than expected by chance.
#' 
#' For each term in \code{thistermcount}, the function:
#' \enumerate{
#'   \item Extracts the count of that term in the whole genome from
#'         \code{alltermcount}
#'   \item Uses the count in the peak list from \code{thistermcount}
#'   \item Calculates the p-value using \code{phyper(q-1, m, n-m, k, lower.tail=FALSE)}
#'         where q is the observed count in peak list, m is the count in genome,
#'         n is total terms in genome, and k is total terms in peak list
#' }
#' 
#' @param alltermcount A list with two components:
#'        \itemize{
#'          \item \code{GOterm}: A character vector of GO/term identifiers present
#'                in the whole genome background
#'          \item \code{GOcount}: A numeric vector of counts for each term in
#'                the whole genome (must be same length as \code{GOterm})
#'        }
#'        This represents the background distribution of terms in the genome.
#' @param thistermcount A list with two components:
#'        \itemize{
#'          \item \code{GOterm}: A character vector of GO/term identifiers found
#'                in the peak list
#'          \item \code{GOcount}: A numeric vector of counts for each term in
#'                the peak list (must be same length as \code{GOterm})
#'        }
#'        This represents the observed terms in the dataset of interest.
#'        Typically obtained from \code{\link{getUniqueGOidCount}}.
#' @param totaltermInGenome A single numeric value representing the total number
#'        of GO/term occurrences in the whole genome background. This is the
#'        sum of all term counts across the entire genome, not just the unique
#'        terms in \code{alltermcount}.
#' @param totaltermInPeakList A single numeric value representing the total
#'        number of GO/term occurrences in the peak list. This is the sum of all
#'        term counts in the peak list, not just the unique terms in
#'        \code{thistermcount}.
#' @return Returns a list with 6 components:
#'        \itemize{
#'          \item \code{thisterm}: Character vector of GO/term identifiers from
#'                \code{thistermcount$GOterm} (same order and length)
#'          \item \code{thistermcount}: Numeric vector of counts for each term
#'                in the peak list (from \code{thistermcount$GOcount})
#'          \item \code{thistermtotal}: Numeric vector of counts for each term
#'                in the whole genome (extracted from \code{alltermcount$GOcount}).
#'                If a term from \code{thistermcount} is not found in
#'                \code{alltermcount}, the count will be 0.
#'          \item \code{pvalue}: Numeric vector of p-values from the hypergeometric
#'                test for each term. P-values represent the probability of
#'                observing at least as many occurrences of the term in the peak
#'                list by chance, given the background distribution.
#'          \item \code{totaltermInPeakList}: The input \code{totaltermInPeakList}
#'                value (for reference)
#'          \item \code{totaltermInGenome}: The input \code{totaltermInGenome}
#'                value (for reference)
#'        }
#'        All vectors have the same length (number of terms in
#'        \code{thistermcount}).
#' 
#' @details
#' \strong{Hypergeometric Test Parameters:}
#' 
#' The hypergeometric test models the following scenario:
#' \itemize{
#'   \item An urn contains n balls total
#'   \item m of these balls are "white" (term of interest in genome)
#'   \item n-m are "black" (other terms in genome)
#'   \item We draw k balls without replacement (terms in peak list)
#'   \item We want to know: what's the probability of drawing at least q white
#'         balls (observed count of term in peak list)?
#' }
#' 
#' The p-value is calculated as:
#' \deqn{P(X \geq q) = 1 - P(X \leq q-1)}
#' 
#' where X ~ Hypergeometric(m, n-m, k).
#' 
#' \strong{Interpretation:}
#' \itemize{
#'   \item Low p-values (< 0.05) indicate significant enrichment: the term
#'         appears more frequently in the peak list than expected by chance
#'   \item High p-values indicate no significant enrichment or potential
#'         depletion
#'   \item Multiple testing correction (e.g., FDR, Bonferroni) should be applied
#'         when testing many terms simultaneously
#' }
#' 
#' @note This is an internal function primarily intended for use by
#'       \code{\link{getEnrichedGO}} and \code{\link{getEnrichedPATH}}. Users
#'       typically should not call this function directly, but rather use the
#'       higher-level enrichment functions.
#' 
#' @author Lihua Julie Zhu
#' @seealso \code{\link[stats:phyper]{phyper}} for the underlying statistical
#'          function, \code{\link{getEnrichedGO}} and \code{\link{getEnrichedPATH}}
#'          for high-level enrichment analysis, \code{\link{getUniqueGOidCount}}
#'          for preparing \code{thistermcount} input
#' @references Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) Univariate
#' Discrete Distributions, Second Edition. New York: Wiley
#' @keywords internal
#' @export
#' @importFrom stats phyper
#' @examples
#' \dontrun{
#' ## Example: Test enrichment of GO terms in a peak list
#' 
#' # GO terms found in the peak list (with duplicates)
#' goList <- c("GO:0000075", "GO:0000082", "GO:0000082", "GO:0000122",
#'             "GO:0000122", "GO:0000075", "GO:0000082", "GO:0000082",
#'             "GO:0000122", "GO:0000122", "GO:0000122", "GO:0000122",
#'             "GO:0000075", "GO:0000082", "GO:000012")
#' 
#' # Background: GO terms and their counts in the whole genome
#' alltermcount <- list(GOterm = c("GO:0000075", "GO:0000082", "GO:000012", 
#'                                  "GO:0000122"), 
#'                      GOcount = c(100, 200, 10, 10))
#' 
#' # Get unique GO term counts in the peak list
#' thistermcount <- getUniqueGOidCount(goList)
#' 
#' # Total counts: sum of all GO term occurrences
#' totaltermInPeakList <- 15  # Total GO term occurrences in peak list
#' totaltermInGenome <- 1000  # Total GO term occurrences in genome
#' 
#' # Perform hypergeometric test
#' result <- hyperGtest(alltermcount, thistermcount, 
#'                      totaltermInGenome, totaltermInPeakList)
#' 
#' # View results
#' result$pvalue  # P-values for each term
#' result$thistermcount  # Counts in peak list
#' result$thistermtotal  # Counts in genome
#' }
#' 
hyperGtest <- function(alltermcount, thistermcount, 
                       totaltermInGenome, totaltermInPeakList) {
    n_terms <- length(thistermcount$GOterm)
    pvalue <- numeric(n_terms)
    thistermtotal <- numeric(n_terms)
    
    for (i in seq_len(n_terms)) {
        # m = number of this GO term in the whole genome
        m <- as.numeric(
            alltermcount$GOcount[alltermcount$GOterm == thistermcount$GOterm[i]]
        )
        # q = number of this GO term in the peak list
        q <- as.numeric(thistermcount$GOcount[i])
        # n = total number of GO terms in the whole genome
        n <- as.numeric(totaltermInGenome)
        # k = total number of GO terms in the peak list
        k <- as.numeric(totaltermInPeakList)
        
        # Hypergeometric test: P(X >= q) where X ~ Hypergeometric(m, n-m, k)
        # Using q-1 and lower.tail=FALSE gives P(X >= q)
        pvalue[i] <- phyper(q - 1L, m, n - m, k, 
                           lower.tail = FALSE, log.p = FALSE)
        thistermtotal[i] <- m
    }
    
    list(
        thisterm = thistermcount$GOterm, 
        thistermcount = thistermcount$GOcount,
        thistermtotal = thistermtotal, 
        pvalue = pvalue, 
        totaltermInPeakList = k, 
        totaltermInGenome = n
    )
}

