#' @title Imports all the WGS data of a single CPCT-02 sample.
#'
#' @details Also runs CHORD for HRD predictions and ShatterSeek for chromothripsis samples.
#'
#' @param cpctId (character): Full CPCT identifier of sample, including the T suffix (used to find files).
#' @param inputFolder (character): Path to the folder containing all WGS data of the sample (combinedData).
#' @return (list) List containing the WGS data of the sample.
#' @examples
#' \dontrun{
#'
#' 	importWGSOfSample(cpctId = '<CPCTID>T', 'combinedData/')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importWGSOfSample <- function(cpctId, inputFolder){

    # Input validation --------------------------------------------------------

    checkmate::assertCharacter(cpctId)
    checkmate::assertAccess(inputFolder, access = 'r')

    sprintf('Importing WGS data of: %s', cpctId) %>% ParallelLogger::logInfo()


    # Import data (WGS) -------------------------------------------------------

    # Initialize empty list containing the WGS data of a single sample.
    dataSample <- list()

    # Import VCF files containing the somatic variants.
    dataSample$somaticVariants <- R2CPCT::importSomaticVariantsVEP(pathVCF = base::list.files(inputFolder, full.names = T, pattern = paste0(cpctId, '\\.purple.somatic.vep.vcf.gz$')))

    # Import the VCF files containing the somatic structural variants (SV).
    svFiles <- base::list.files(inputFolder, full.names = T, pattern = paste0(cpctId, '\\.purple.sv.*vcf.gz$'))
    svFiles <- base::ifelse(base::length(svFiles) > 1, svFiles[1], svFiles)

    dataSample$structuralVariants <- R2CPCT::importStructuralVariantsPURPLE(pathSV = svFiles)

    # Import the HMF driver catalog.
    driverPath <- base::list.files(inputFolder, full.names = T, pattern = paste0(cpctId, '\\.driver.catalog.somatic.tsv$'))
    if(S4Vectors::isEmpty(driverPath)) driverPath <- base::list.files(inputFolder, full.names = T, pattern = paste0(cpctId, '\\.driver.catalog.tsv$'))

    dataSample$driverCatalog <- R2CPCT::importdriverCatalogHMF(pathCatalog = driverPath)

    # Import the purity statistics.
    dataSample$purityStats <- R2CPCT::importPurityStatsPURPLE(pathStats = base::list.files(inputFolder, full.names = T, pattern = paste0(cpctId, '\\.purple.purity.tsv$')))

    # Import the copy-number alterations.
    dataSample$copyNumbers <- R2CPCT::importSomaticCopynumberPURPLE(pathCNV = base::list.files(inputFolder, full.names = T, pattern = paste0(cpctId, '\\.purple.cnv.somatic.tsv$')))

    # Add gender to copynumber information.
    dataSample$copyNumbers$gender <- dataSample$purityStats$gender

    # Import LINX - Fusions.
    dataSample$fusionsLINX <- R2CPCT::importLINXFusions(pathFusions = base::list.files(inputFolder, full.names = T, pattern = paste0(cpctId, '\\.linx.fusion.tsv$')))


    # RUN CHORD ---------------------------------------------------------------

    dataSample$CHORD <- R2CPCT::performCHORD(somaticVariants = dataSample$somaticVariants, structuralVariants = dataSample$structuralVariants)


    # Run ShatterSeek ---------------------------------------------------------

    if(!is.null(dataSample$structuralVariants)){
        dataSample$shatterSeek <- R2CPCT::performShatterSeek(copyNumbers = dataSample$copyNumbers, structuralVariants = dataSample$structuralVariants)
    }


    # Return statement --------------------------------------------------------

    return(dataSample)

}

#' @title Imports all the WGS data of multiple samples (cohort).
#'
#' @param cpctIds (character): Full CPCT identifiers of samples, including the T suffixes (used to find files).
#' @param inputFolder (character): Path to the folder containing all WGS data of the samples (combinedData).
#' @param nThreads (integer): Number of threads to use; default is 1.
#' @param performAggregation (logical): Should aggregation be performed? (Is needed for further processing)
#' @return (list) List containing the aggregated WGS data of the underlying samples.
#' @examples
#' \dontrun{
#'
#' 	importWGSOfCohort(cpctId = '<CPCTIDs>', 'combinedData/')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importWGSOfCohort <- function(cpctIds, inputFolder, nThreads = 1, performAggregation = T){

    # Input validation --------------------------------------------------------

    checkmate::assertCharacter(cpctIds)
    checkmate::assertAccess(inputFolder, access = 'r')
    checkmate::assertNumber(nThreads)
    checkmate::assertLogical(performAggregation)

    sprintf('Importing WGS data of %s samples.', dplyr::n_distinct(cpctIds)) %>% ParallelLogger::logInfo()

    # Import data (WGS) -------------------------------------------------------

    data.PerSample <- pbapply::pblapply(cpctIds, function(x){ importWGSOfSample(x, inputFolder) }, cl = nThreads)
    base::names(data.PerSample) <- cpctIds


    # Aggregate the samples ---------------------------------------------------

    if(performAggregation){
        sprintf('Aggregating WGS data of %s samples.', dplyr::n_distinct(cpctIds)) %>% ParallelLogger::logInfo()

        data.AllSamples <- list()

        data.AllSamples$somaticVariants <- base::unlist(GenomicRanges::GRangesList(base::lapply(data.PerSample, function(x) x$somaticVariants)))
        data.AllSamples$structuralVariants <- base::unlist(GenomicRanges::GRangesList(base::lapply(data.PerSample, function(x) x$structuralVariants)))
        data.AllSamples$copyNumbers <- base::unlist(GenomicRanges::GRangesList(base::lapply(data.PerSample, function(x) x$copyNumbers)))
        data.AllSamples$driverCatalog <- dplyr::bind_rows(base::lapply(data.PerSample, function(x) x$driverCatalog))
        data.AllSamples$purityStats <- dplyr::bind_rows(base::lapply(data.PerSample, function(x) x$purityStats))
        data.AllSamples$fusionsLINX <- dplyr::bind_rows(base::lapply(data.PerSample, function(x) x$fusionsLINX))
        data.AllSamples$CHORD <- dplyr::bind_rows(base::lapply(data.PerSample, function(x) x$CHORD$results))
        data.AllSamples$CHORD.contexts <- dplyr::bind_rows(base::lapply(data.PerSample, function(x) tibble::as_tibble(x$CHORD$contexts, rownames = 'sample')))
        data.AllSamples$shatterSeek <- dplyr::bind_rows(base::lapply(data.PerSample, function(x) x$shatterSeek))

    }else{

        # Return non-aggregated samples.
        return(data.PerSample)

    }


    # Return statement --------------------------------------------------------

    return(data.AllSamples)

}
