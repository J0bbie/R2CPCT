#' @title Retrieves the WGS files needed for analysis into the same folder.
#'
#' @description Retrieves the files (grep) based on CPCT identifiers and symlinks them into the same folder for easier manipulation.
#'
#' @param pathHMF (character): Path to main folder containing the WGS files of the data-request (DR).
#' @param pathOutput (character): Path to output folder which will contain the symlinks.
#' @param cpctIds (character): The samples which will be retrieved.
#'
#' @examples
#' \donttest{
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
  
  # Retrieve all files.  
  sprintf('\tRetrieving an inventory of all relevant files (recursively). this can take some time for folders with MANY files.') %>% ParallelLogger::logInfo()
  data.Files <- tibble::tibble(file = base::list.files(pathHMF, full.names = TRUE, recursive = TRUE, pattern = 'purple|linx|driver'))
  data.Files <- data.Files %>% dplyr::mutate(baseFile = base::basename(file), sample = base::gsub('\\..*', '', baseFile))
  
  # Retrieve the relevant files.
  data.Cohort <- data.Files %>% 
    dplyr::filter(sample %in% cpctIds) %>% 
    dplyr::filter(
      base::grepl(x = baseFile, pattern = 'purple.somatic_annotedWithVEP.vcf.gz|purple.sv.vcf.gz|.purple.cnv.gene.tsv|.purple.cnv.somatic.tsv|.purple.purity.tsv|.purple.qc|driver.catalog|linx.fusion|linx.drivers')
    )
  
  # Generate symlinks.
  symInfo <- base::file.symlink(from = base::unique(data.Cohort$file), to = pathOutput)
  
  
  # Return statement --------------------------------------------------------
  
  sprintf('\tSuccesfully symlinked %s / %s files.', base::sum(symInfo), base::length(symInfo)) %>% ParallelLogger::logInfo()
  return(NULL)
  
}
