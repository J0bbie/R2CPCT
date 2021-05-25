#' @title Import drivers (LOH / AMP / DEL) detected by LINX.
#'
#' @param pathDrivers (character): Path to the <sample>.linx.drivers.tsv file containing predicted drivers.
#' @return (tibble) Tibble containing the predicted drivers of the sample.
#' @examples
#' \dontrun{
#'
#' 	importLINXDrivers(pathDrivers = '<sample>.linx.drivers.tsv')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importLINXDrivers <- function(pathDrivers){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathDrivers, access = 'r')
    geneInfo <- tibble::as_tibble(S4Vectors::mcols(R2CPCT::GENCODE.v35)) %>% dplyr::select(SYMBOL, ENSEMBL)


    # Read driver catalog -----------------------------------------------------

    sprintf('Importing LINX driver-events: %s', pathDrivers) %>% ParallelLogger::logInfo()

    # Clean seqlevels and add chromosome information.
    sample.Drivers <- readr::read_tsv(file = pathDrivers, col_types = 'ccccccccccccccccccccccccccccccccccccccccccc')

    # Add sample name
    sample.Drivers$sample <- base::factor(base::gsub('\\.linx.*', '', base::basename(pathDrivers)))


    # Clean-up data. ----------------------------------------------------------

    if(nrow(sample.Drivers) > 0){
        # Retrieve the common gene name.
        commonSYMBOL <- limma::alias2SymbolTable(sample.Drivers$gene)
        sample.Drivers$SYMBOL <- base::ifelse(is.na(commonSYMBOL), sample.Drivers$gene, commonSYMBOL)

        # Add ENSEMBL identifier.
        sample.Drivers <- sample.Drivers %>% dplyr::left_join(geneInfo, by = c('SYMBOL' = 'SYMBOL'))

        # Group multiple events per gene.
        sample.Drivers <- sample.Drivers %>%
            dplyr::group_by(gene) %>%
            dplyr::summarise(dplyr::across(dplyr::everything(), list(~ base::paste(base::unique(.), collapse = ', '))))

        colnames(sample.Drivers) <- base::gsub('_1', '', base::colnames(sample.Drivers))

    }else{
        sample.Drivers <- NULL
    }


    # Return statement --------------------------------------------------------

    sprintf('\tReturning LINX driver-events') %>% ParallelLogger::logTrace()

    return(sample.Drivers)

}
