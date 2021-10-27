#' @title Run LINX on one or multiple samples to discover fusion and driver-events from PURPLE output.
#'
#' @description For installation instructions, please see the 'howToInstallLINX' vignette.
#' It will output the (relevant) files in the combinedData/ folder.
#'
#' @param pathCombined (character): Path to the combinedData/ folder containing the symlinks to the WGS files of the data-request (DR).
#' @param cpctIds (character): The samples which will be retrieved.
#' @param nThreads (integer): Number of cores to run multiple LINX in parallel.
#' @param dryRun (logical): Only output the tmp. file which stores the commands (TRUE) or run the commands (FALSE).
#' @param pathLINX (character): Path to the folder containing LINX and the respective databases.
#'
#' @examples
#' \donttest{
#'
#'  runLINX('DR-071/combinedData/', c('CPCT02124567', 'CPCT02012345123'), nThreads = 10)
#'
#' }
#' @return (NULL) Does not return anything, only outputs logging messages.
#' @export
performLINX <- function(pathCombined, cpctIds, nThreads, dryRun = TRUE, pathLINX = '/mnt/onco0002/repository/software/LINX_v1.16/'){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathCombined, access = 'rw')
    checkmate::assertCharacter(cpctIds)
    checkmate::assertNumber(nThreads)
    checkmate::assertLogical(dryRun)

    # Logging.
    sprintf('Performing LINX on %s samples', dplyr::n_distinct(cpctIds)) %>% ParallelLogger::logInfo()


    # Generating LINX commands ------------------------------------------------

    sprintf('\tGenerating LINX commands.') %>% ParallelLogger::logInfo()

    commands.LINX <- base::unlist(base::lapply(cpctIds, function(cpctId){

        # Retrieve path to PURPLE/GRIDSS SV file.
        svFile <- base::list.files(pathCombined, full.names = TRUE, pattern = paste0(cpctId, '\\.purple.sv.*vcf.gz$'))

        sprintf('java -jar %s/linx_v1.16.jar
        -sample %s -sv_vcf %s -purple_dir %s -output_dir %s -check_fusions -check_drivers
        -known_fusion_file %s/known_fusion_data.37.csv
        -gene_transcripts_dir %s/ENSEMBLv101/
        -fragile_site_file %s/fragile_sites_hmf.37.csv
        -line_element_file %s/line_elements.37.csv
        -viral_hosts_file %s/viral_host_ref.csv',
                pathLINX, cpctId, svFile, pathCombined, pathCombined, pathLINX, pathLINX, pathLINX, pathLINX, pathLINX) %>% gsub('\n', '', .)
    }))


    # Running LINX ------------------------------------------------------------

    sprintf('\tRunning LINX commands.') %>% ParallelLogger::logInfo()

    tmpFile <- tempfile()
    utils::write.table(commands.LINX, tmpFile, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

    if(dryRun){
        sprintf('Run the bash commands in the following (tmp) file to perform LINX on your samples:\n %s', sprintf('cat %s | parallel -j %s', tmpFile, nThreads)) %>% ParallelLogger::logInfo()
        return(tmpFile)
    }

    command.LINX <- sprintf('cat %s | parallel -j %s', tmpFile, nThreads)
    outputBash <- base::system(paste("/bin/bash -c", base::shQuote(command.LINX)))


    # Cleaning up LINX --------------------------------------------------------

    sprintf('\tCleaning up LINX files.') %>% ParallelLogger::logInfo()
    base::file.remove(base::list.files(pathCombined, pattern = '\\.linx.vis_|\\.linx.viral_inserts.tsv|\\.linx.svs.tsv|\\.linx.links.tsv|linx.version|\\.linx.clusters.tsv|\\.linx.breakend.tsv', full.names = TRUE))


    # Return statement --------------------------------------------------------

    return(NULL)

}
