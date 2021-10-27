#' @title Sorting function based on mut. exclusivity.
#'
#' @param M (dataframe): dataframe of genes vs samples.
#'
#' @return (dataframe) Sorted dataframe of genes vs samples.
#' @examples
#' \donttest{
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @export
memoSort <- function(M) {

    # Input validation --------------------------------------------------------

    checkmate::assertDataFrame(M)


    # memoSort ----------------------------------------------------------------

    geneOrder <- base::sort(rowSums(M), decreasing = TRUE, index.return = TRUE)$ix

    scoreCol <- function(x) {
        score <- 0
        for(i in seq_len(length(x))) {
            if(x[i]) {
                score <- score + 2^(length(x)-i)
            }
        }
        return(score)
    }

    scores <- base::apply(M[geneOrder, ], 2, scoreCol)

    sampleOrder <- base::sort(scores, decreasing = TRUE, index.return = TRUE)$ix

    return(M[geneOrder, sampleOrder])

}
