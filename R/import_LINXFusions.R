#' @title Import detected fusions by LINX.
#'
#' @param pathFusions (character): Path to the <sample>.linx.fusion.tsv file containing predicted fusions.
#' @return (tibble) Tibble containing the predicted WGS fusions of the sample.
#' @examples
#' \donttest{
#'
#' 	importLINXFusions(pathFusions = '<sample>.linx.fusion.tsv')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importLINXFusions <- function(pathFusions){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathFusions, access = 'r')


    # Read driver catalog -----------------------------------------------------

    sprintf('Importing LINX fusions: %s', pathFusions) %>% ParallelLogger::logInfo()

    # Clean seqlevels and add chromosome information.
    sample.Fusions <- readr::read_tsv(file = pathFusions, col_types = 'ddclcccddlccddddccccccd')

    # Add sample name
    sample.Fusions$sample <- base::factor(base::gsub('\\.linx.*', '', base::basename(pathFusions)))

    # Return statement --------------------------------------------------------

    sprintf('\tReturning LINX Fusions.') %>% ParallelLogger::logTrace()

    return(sample.Fusions)

}
