#' @title Find corresponding center based on CPCT number.
#'
#' @param x (character): CPCT identifiers.
#'
#' @return (character) Character vector of matching centers.
#' @examples
#' \donttest{
#'
#' 	meta$Center <- findCPCTCenter(CPCTId)
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
findCPCTCenter <- function(x){

    # Input validation --------------------------------------------------------

    checkmate::assertCharacter(x)


    # Read sites --------------------------------------------------------------

    siteCodes <- utils::read.delim(system.file("extdata/CPCTSites.txt", package="R2CPCT"), as.is = TRUE, colClasses = 'character')


    # Overlap -----------------------------------------------------------------

    centers <- siteCodes[match(substr(gsub('(CPCT|DRUP|WIDE)..', '', x), 0, 2), siteCodes$Code),]$Site


    # Return matched centers. -------------------------------------------------

    return(centers)

}
