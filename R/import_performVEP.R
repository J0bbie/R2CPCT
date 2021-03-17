#' @title Generate the commands to run VEP on HMF/CPCT-02 samples.
#'
#' @description Generate the bash commands to run VEP on HMF/CPCT-02 samples.
#' Please see the vignette 'howToInstallVep' on instructions to install VEP.
#'
#' The accompanying bash script is provided in inst/performVEPInParallel.sh.
#' This function will return the command which should be run in bash.
#'
#' By default, it performs annotation based on ensembl-vep (v101) and GENCODE v35 (hg19) gene-structures.
#'
#' @param pathVCF (character): Path to folder containing the CPCT-02 VCF files (combinedData/<sample>.purple.somatic.vcf.gz).
#' @param pathOutput (character): Path to output folder for VEP-annotated VCF files (if left empty, pathVCF will be used).
#' @param cpctIds (character): The samples which will be annotated by VEP.
#' @param fastaFile (character): Path to genome FASTA.
#' @param pathGTF (character): Path to GTF (see vignette on VEP install).
#'
#' @examples
#' \dontrun{
#'
#'  performVEP('DR-071/combinedData/', cpctIds = c('CPCT02124567', 'CPCT02012345123'))
#'
#' }
#' @return (character) Returns the command which should be performed in bash.
#' @export
performVEP <- function(pathVCF, pathOutput = NULL, cpctIds, fastaFile = '/mnt/data/ccbc_environment/general/genomes/hsapiens/hg19_HMF/Homo_sapiens.GRCh37.GATK.illumina.fasta', pathGTF = '/mnt/data/ccbc_environment/general/annotation/hg19/GENCODE/noChrPrefix_VEP_gencode.v35lift37.annotation.gtf.gz'){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathVCF, access = 'r')
    if(is.null(pathOutput)) pathOutput <- pathVCF
    checkmate::assertAccess(pathOutput, access = 'rw')
    checkmate::assertCharacter(cpctIds)
    checkmate::assertAccess(fastaFile, access = 'r')
    checkmate::assertAccess(pathGTF, access = 'r')

    # Logging.
    sprintf('Going to perform VEP on %s samples', dplyr::n_distinct(cpctIds)) %>% ParallelLogger::logInfo()
    sprintf('\tpathVCF:\t%s', pathVCF) %>% ParallelLogger::logTrace()
    sprintf('\tpathOutput:\t%s', pathOutput) %>% ParallelLogger::logTrace()
    sprintf('\tfastaFile:\t%s', fastaFile) %>% ParallelLogger::logTrace()
    sprintf('\tpathGTF:\t%s', pathGTF) %>% ParallelLogger::logTrace()
    sprintf('\tUnique samples:\t%s', dplyr::n_distinct(cpctIds)) %>% ParallelLogger::logTrace()


    # Generate VEP command ----------------------------------------------------

    command.VEP <- base::sprintf('bash %s -v %s -e .purple.somatic.vcf.gz -o %s -f %s -g %s',
                                 system.file('performVEPInParallel.sh', package = 'R2CPCT'),
                                 pathVCF,
                                 pathOutput,
                                 fastaFile,
                                 pathGTF
    )

    sprintf('Paste this command in your terminal and perform using parallel: <cmd> | parallel -j <nThreads>') %>% ParallelLogger::logInfo()
    sprintf(command.VEP) %>% ParallelLogger::logInfo()


    # Return statement --------------------------------------------------------

    return(command.VEP)

}