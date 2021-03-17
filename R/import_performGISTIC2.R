#' @title Detect genes under selection by focal amplification/deletion using the GISTIC2.0 algorithm.
#' @description This can take anywhere between 30 min. up to several hours depending on the number of samples.
#' If you have previously ran GISTIC, use the \link[R2CPCT]{importGISTIC2} function to import and clean the data.
#'
#' @param regions (GRanges): GRanges of CN segments with the following columns: sample, copyNumber, bafCount, sex.
#' @param outputFolder (character): Path to output folder of GISTIC2 results.
#' @param gisticPath (character): Path to gistic2 shell script.
#' @param gisticParameters (character): Optional parameters given to gistic2. See online documentation for options.
#' @param gisticAnno (character): Path to which genome annotation GISTIC should use? Stored in /mnt/data/ccbc_environment/software/general/GISTIC2_2.0.23/
#'
#' @return (list) GISTIC
#' @examples
#' \dontrun{
#'
#' 	gisticOutput <- R2CPCT::performGISTIC2(regions = data.Cohort$copyNumbers)
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family GISTIC
#' @export
performGISTIC2 <- function(regions, outputFolder, gisticPath = '/mnt/data/ccbc_environment/software/general/GISTIC2_2.0.23/gistic2', gisticParameters = '-genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.3 -td 0.3 -savedata 0 -savegene 1 -qvt 0.1 -twoside 0', gisticAnno = '/mnt/data/ccbc_environment/software/general/GISTIC2_2.0.23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat'){

    # Input validation --------------------------------------------------------

    checkmate::checkClass(regions, 'GRanges')

    if(is.null(regions$sample)) stop('Provide the following column containing sample identifiers: sample')
    if(is.null(regions$copyNumber)) stop('Provide the following column containing the absolute copy-numbers of the segment: copyNumber')
    if(is.null(regions$bafCount)) stop('Provide the following column containing the number of heterozygous markers used in segmentation: bafCount')
    if(is.null(regions$gender)) stop('Provide the following column containing the gender of the patient (MALE or FEMALE): gender')

    checkmate::assertAccess(outputFolder, access = 'rw')
    checkmate::assertFileExists(gisticPath, access = 'x')
    checkmate::assertCharacter(gisticParameters, null.ok = T)
    checkmate::assertFileExists(gisticAnno, access = 'r')

    base::sprintf('Generating command for GISTIC 2.0 on %s samples.', dplyr::n_distinct(regions$sample)) %>% ParallelLogger::logInfo()


    # Generate GISTIC2 input --------------------------------------------------

    base::sprintf('\tConverting CN-segments to GISTIC2 format.') %>% ParallelLogger::logInfo()

    # Fix absolute CN below zero.
    regions$copyNumber <- base::ifelse(regions$copyNumber < 0, 0, regions$copyNumber)

    # Add a copy-number to the sex-chromosomes (so a log2(1) -1 equals 0 instead of -1 for the X and Y chromosomes)
    regions[base::grepl('chrX|chrY', GenomeInfoDb::seqnames(regions)) & regions$gender == 'MALE',]$copyNumber <- regions[base::grepl('chrX|chrY', GenomeInfoDb::seqnames(regions)) & regions$gender == 'MALE',]$copyNumber + 1

    # Convert to GISTIC format.
    gistic.regions <- data.frame(regions$sample, GenomeInfoDb::seqnames(regions), IRanges::start(regions), IRanges::end(regions), regions$bafCount, copynumber.log2 = log2(regions$copyNumber) - 1)
    gistic.regions$copynumber.log2 <- base::ifelse(gistic.regions$copynumber.log2 <= -15, -10, gistic.regions$copynumber.log2)

    # Create temporary files for input/output results.
    temp.input <- base::tempfile(pattern = '', fileext = '.txt')

    # Write regions to temp. file.
    utils::write.table(gistic.regions, temp.input, col.names = F, append = F, quote = F, row.names = F, sep = '\t')

    # Generate GISTIC command.
    command.GISTIC <- base::sprintf('%s -b %s -seg %s -refgene %s %s',
                                    base::Sys.which(gisticPath),
                                    outputFolder,
                                    temp.input,
                                    gisticAnno,
                                    if(!is.null(gisticParameters)){ gisticParameters } else{ '' }
    )

    base::sprintf('\tRun this command in your terminal and import the data using R2CPCT::importGISTIC2(\'%s\'):', outputFolder) %>% ParallelLogger::logInfo()
    base::sprintf('\t%s', command.GISTIC) %>% ParallelLogger::logInfo()

}