#' @title Import the driver catalog as determined by HMF.
#'
#' @param pathStats (character): Path to the <sample>.purple.purity.tsv file containing HMF-determined drivers.
#' @return (tibble) Tibble containing PURPLE statistics of the sample.
#' @examples
#' \dontrun{
#'
#' 	importPurityStatsPURPLE(pathStats = '<sample>.purple.purity.tsv')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importPurityStatsPURPLE <- function(pathStats){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathStats, access = 'r')


    # Read driver catalog -----------------------------------------------------

    sprintf('Importing PURPLE purity stats: %s', pathStats) %>% ParallelLogger::logInfo()

    # Clean seqlevels and add chromosome information.
    sample.Stats <- readr::read_tsv(file = pathStats, col_types = readr::cols())

    # Add sample name
    sample.Stats$sample <- base::factor(base::gsub('\\.purple.*', '', base::basename(pathStats)))

    # Return statement --------------------------------------------------------

    sprintf('\tReturning PURPLE purity stats') %>% ParallelLogger::logTrace()

    return(sample.Stats)

}
