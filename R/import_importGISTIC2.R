#' @title Imports the GISTIC2 recurrent peak and chromosomal arms output.
#' @description Imports and cleans-up the output from \link[R2CPCT]{performGISTIC2}.
#'
#' @param gisticFolder (character): Path to GISTIC output folder.
#'
#' @examples
#' \dontrun{
#'
#' 	Folder to the GISTIC2 output.
#' 	importGISTIC2('/tmp/etc/etc/')
#'
#' }
#'
#' @return (list) Cleaned-up GISTIC output.
#'
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @export
importGISTIC2 <- function(gisticFolder){

    # Input validation --------------------------------------------------------

    checkmate::checkCharacter(gisticFolder)

    # Required files.
    gisticAllLesionsFile <- base::list.files(gisticFolder, pattern = 'all_lesions.conf.*txt', full.names = T)
    gisticScoresFile <- base::list.files(gisticFolder, pattern = 'scores\\.gistic', full.names = T)
    gisticBroadFile <- base::list.files(gisticFolder, pattern = 'broad_significance_results.txt', full.names = T)
    gisticBroadPerArmFile <- base::list.files(gisticFolder, pattern = 'broad_values_by_arm.txt', full.names = T)

    base::sprintf('Importing GISTIC2 results from: %s', gisticFolder) %>% ParallelLogger::logInfo()

    # Select protein-coding genes only.
    GENCODE.v35 <- R2CPCT::GENCODE.v35[R2CPCT::GENCODE.v35$gene_type == 'protein_coding',]

    # Import recurrently-detected peaks ---------------------------------------

    base::sprintf('\tImporting recurrently-detected peaks.') %>% ParallelLogger::logInfo()

    # Import peaks, clean empty colunms.
    gisticAllLesionsCN <- readr::read_delim(gisticAllLesionsFile, delim = '\t', trim_ws = T)
    gisticAllLesionsCN[grepl('^X', colnames(gisticAllLesionsCN))] <- NULL
    gisticAllLesionsCN <- gisticAllLesionsCN[!grepl('CN values', gisticAllLesionsCN$`Unique Name`),]

    # Annotate the peaks with overlapping genes and drivers -------------------

    # Generate GRanges from discovered wide peaks.
    peakSites <- GenomicRanges::GRanges(base::gsub('\\(.*', '', gisticAllLesionsCN$`Wide Peak Limits`))
    GenomicRanges::mcols(peakSites) <- gisticAllLesionsCN

    # Initialize empty columns (if no overlap).
    peakSites$nGenes.Drivers <- peakSites$overlapGenes.Drivers <- peakSites$nGenes.GENCODE <- peakSites$overlapGenes.GENCODE <- NA

    # Find overlap with genes. Concatenate per overlap.
    overlapGenes.GENCODE <- tibble::as_tibble(IRanges::findOverlaps(peakSites, GENCODE.v35, minoverlap = 10))
    overlapGenes.GENCODE <- overlapGenes.GENCODE %>% dplyr::group_by(queryHits) %>% dplyr::summarise(overlappingGenes = paste(unique(unlist(list(stats::na.omit(GENCODE.v35[subjectHits]$SYMBOL)))), collapse = ', '), nGenes = length(unique(unlist(list(stats::na.omit(GENCODE.v35[subjectHits]$ENSEMBL))))))
    overlapGenes.GENCODE$overlappingGenes <- paste(overlapGenes.GENCODE$overlappingGenes, '(GENCODE v35)')

    overlapGenes.Drivers <- tibble::as_tibble(IRanges::findOverlaps(peakSites, R2CPCT::driverList, minoverlap = 10))
    overlapGenes.Drivers <- overlapGenes.Drivers %>% dplyr::group_by(queryHits) %>% dplyr::summarise(overlappingGenes = paste(unique(unlist(list(stats::na.omit(R2CPCT::driverList[subjectHits]$SYMBOL)))), collapse = ', '), nGenes = length(unique(unlist(list(stats::na.omit(R2CPCT::driverList[subjectHits]$ENSEMBL))))))
    overlapGenes.Drivers$overlappingGenes <- sprintf('%s (Put. Driver; n = %s)', paste(overlapGenes.Drivers$overlappingGenes, ''), overlapGenes.GENCODE[match(overlapGenes.Drivers$queryHits, overlapGenes.GENCODE$queryHits, nomatch = NULL),]$nGenes)

    # Add genes to peaks.
    peaksWithOverlap <- peakSites[overlapGenes.GENCODE$queryHits]

    # GENCODE
    peaksWithOverlap$overlapGenes.GENCODE <- overlapGenes.GENCODE$overlappingGenes
    peaksWithOverlap$nGenes.GENCODE <- overlapGenes.GENCODE$nGenes
    peaksWithOverlap$originalPeak <- overlapGenes.GENCODE$queryHits

    # Drivers
    peaksWithOverlap[peaksWithOverlap$originalPeak %in% overlapGenes.Drivers$queryHits]$overlapGenes.Drivers <- overlapGenes.Drivers$overlappingGenes
    peaksWithOverlap[peaksWithOverlap$originalPeak %in% overlapGenes.Drivers$queryHits]$nGenes.Drivers <- overlapGenes.Drivers$nGenes

    #Find nearest up/downstream gene of peaks without overlapping genes.
    peaksWithoutOverlap <- peakSites[!1:length(peakSites) %in% overlapGenes.GENCODE$queryHits]
    peaksWithoutOverlap$overlapGenes.GENCODE <- GENCODE.v35[IRanges::nearest(peaksWithoutOverlap, GENCODE.v35)]$SYMBOL
    if(length(peaksWithoutOverlap) != 0) peaksWithoutOverlap$overlapGenes.GENCODE <- paste(peaksWithoutOverlap$overlapGenes.GENCODE, '(Nearest)')
    if(length(peaksWithoutOverlap) != 0) peaksWithoutOverlap$nGenes.GENCODE <- 0

    # Combine the annotated peaks with the peaks having no annotations.
    gisticNarrowPeaksWithAnno <- GenomicRanges::sort(c(peaksWithOverlap, peaksWithoutOverlap))

    # Choose the shown annotation of the peak.
    gisticNarrowPeaksWithAnno$overlapGenes.Final <- ifelse(is.na(gisticNarrowPeaksWithAnno$overlapGenes.Drivers), gisticNarrowPeaksWithAnno$overlapGenes.GENCODE, gisticNarrowPeaksWithAnno$overlapGenes.Drivers)
    gisticNarrowPeaksWithAnno$originalPeak <- NULL

    # Reduce the number of genes shown.
    gisticNarrowPeaksWithAnno$overlapGenes.Final <- base::ifelse((stringr::str_count(gisticNarrowPeaksWithAnno$overlapGenes.Final, pattern = ", ") + 1) > 6, base::paste0(stringr::str_count(gisticNarrowPeaksWithAnno$overlapGenes.Final, pattern = ", ") + 1, ' genes ', gsub('.*\\(', '\\(', gisticNarrowPeaksWithAnno$overlapGenes.Final)), gisticNarrowPeaksWithAnno$overlapGenes.Final)


    # Import GISTIC scores ----------------------------------------------------

    base::sprintf('\tImporting all GISTIC2 peaks.') %>% ParallelLogger::logInfo()
    gisticPeakScores <- readr::read_delim(gisticScoresFile, delim = '\t', trim_ws = T)


    # Import broad-level results ----------------------------------------------

    base::sprintf('\tImporting GISTIC2 arm-level results peaks.') %>% ParallelLogger::logInfo()
    if(!S4Vectors::isEmpty(gisticBroadFile)){
        gisticBroadScores <- readr::read_delim(gisticBroadFile, delim = '\t', trim_ws = T)
    }else{
        gisticBroadScores <- NULL
    }

    if(!S4Vectors::isEmpty(gisticBroadPerArmFile)){
        gisticBroadScoresPerArm <- readr::read_delim(gisticBroadPerArmFile, delim = '\t', trim_ws = T) %>%
            reshape2::melt(id.vars = 'Chromosome Arm') %>%
            dplyr::mutate(`Chromosome Arm` = factor(`Chromosome Arm`, levels = gtools::mixedsort(unique(`Chromosome Arm`))))
    }else{
        gisticBroadScoresPerArm <- NULL
    }

    # Return objects ----------------------------------------------------------

    base::sprintf('\tCombining and returning GISTIC2 results.') %>% ParallelLogger::logInfo()

    dataGISTIC <- list()
    dataGISTIC$gisticPeakScores <- gisticPeakScores
    dataGISTIC$gisticBroadScores <- gisticBroadScores
    dataGISTIC$gisticBroadScoresPerArm <- gisticBroadScoresPerArm
    dataGISTIC$gisticNarrowPeaksWithAnno <- gisticNarrowPeaksWithAnno

    return(dataGISTIC)

}
