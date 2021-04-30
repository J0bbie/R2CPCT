#' @title Retrieves the WGS files needed for analysis into the same folder.
#'
#' @description Retrieves the files (grep) based on CPCT identifiers and symlinks them into the same folder for easier manipulation.
#'
#' @param pathHMF (character): Path to main folder containing the WGS files of the data-request (DR).
#' @param pathOutput (character): Path to output folder which will contain the symlinks.
#' @param cpctIds (character): The samples which will be retrieved.
#'
#' @examples
#' \dontrun{
#'
#'  retrieveFilesWGS(
#'  pathHMF = 'DR-071/rawData/WGS/',
#'  pathOutput = 'DR-071/combinedData/',
#'  cpctIds = c('CPCT02124567', 'CPCT02012345123')
#'  )
#'
#' }
#' @return (NULL) Does not return anything, only outputs logging messages.
#' @export
retrieveFilesWGS <- function(pathHMF, pathOutput, cpctIds){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathHMF, access = 'r')
    checkmate::assertAccess(pathOutput, access = 'rw')
    checkmate::assertCharacter(cpctIds)

    # Logging.
    sprintf('Symlinking relevant WGS files into a single folder.') %>% ParallelLogger::logInfo()
    sprintf('\tpathHMF:\t%s', pathHMF) %>% ParallelLogger::logInfo()
    sprintf('\tpathOutput:\t%s', pathOutput) %>% ParallelLogger::logInfo()
    sprintf('\tUnique samples:\t%s', dplyr::n_distinct(cpctIds)) %>% ParallelLogger::logInfo()


    # Retrieve files ----------------------------------------------------------

    files <- base::list()

    # Somatic VCF.
    files$somaticVCF <- base::list.files(pathHMF, pattern = base::paste(base::paste0(cpctIds, '\\.purple.somatic.vcf.gz'), collapse = '|'), full.names = T, recursive = T)

    # Purple / GRIDSS
    files$purple <- list.files(pathHMF, pattern = paste(paste0(cpctIds, '\\.purple\\.[cnv|sv|purity|qc]'), collapse = '|'), full.names = T, recursive = T)
    files$purple <- files$purple[!grepl('range.tsv', files$purple)]

    # Driver files.
    files$driverCatalog <- list.files(pathHMF, pattern = paste(paste0(cpctIds, '\\.driver.catalog.*tsv'), collapse = '|'), full.names = T, recursive = T)

    # LINX.
    files$fusion <- list.files(pathHMF, pattern = paste(paste0(cpctIds, '\\.linx.fusion.tsv'), collapse = '|'), full.names = T, recursive = T)

    # LINX drivers.
    files$linxDriver <- list.files(pathHMF, pattern = paste(paste0(cpctIds, '\\.linx.drivers.tsv'), collapse = '|'), full.names = T, recursive = T)

    # Generate symlinks.
    symInfo <- base::file.symlink(from = base::unique(c(files$somaticVCF, files$purple, files$driverCatalog, files$fusion, files$linxDriver)), to = pathOutput)


    # Return statement --------------------------------------------------------

    sprintf('\tSuccesfully symlinked %s / %s files.', base::sum(symInfo), base::length(symInfo)) %>% ParallelLogger::logInfo()
    return(NULL)

}