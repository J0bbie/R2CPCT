#' @title Import the driver catalog as determined by HMF.
#'
#' @param pathCatalog (character): Path to the <sample>.driver.catalog.tsv file containing HMF-determined drivers.
#' @return (tibble) Tibble containing the drivers of the sample.
#' @examples
#' \dontrun{
#'
#' 	importdriverCatalogHMF(pathCatalog = '<sample>.driver.catalog.tsv')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importdriverCatalogHMF <- function(pathCatalog){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathCatalog, access = 'r')

    geneInfo <- tibble::as_tibble(S4Vectors::mcols(R2CPCT::GENCODE.v35)) %>% dplyr::select(SYMBOL, ENSEMBL)


    # Read driver catalog -----------------------------------------------------

    sprintf('Importing HMF driver catalog: %s', pathCatalog) %>% ParallelLogger::logInfo()

    # Clean seqlevels and add chromosome information.
    sample.Drivers <- readr::read_tsv(file = pathCatalog, col_types = 'ccccccdddddddldd')

    # Add sample name
    sample.Drivers$sample <- base::factor(base::gsub('\\.driver.*', '', base::basename(pathCatalog)))


    # Clean-up data. ----------------------------------------------------------

    if(nrow(sample.Drivers) > 0){
        # Retrieve the common gene name.
        commonSYMBOL <- limma::alias2SymbolTable(sample.Drivers$gene)
        sample.Drivers$SYMBOL <- base::ifelse(base::is.na(commonSYMBOL), sample.Drivers$gene, commonSYMBOL)

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

    sprintf('\tReturning driver-catalog') %>% ParallelLogger::logTrace()

    return(sample.Drivers)

}
