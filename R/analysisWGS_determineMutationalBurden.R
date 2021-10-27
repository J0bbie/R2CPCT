#' @title Determine the mutational burden of the CPCT-02/HMF samples.
#'
#' @param data.Cohort (list): Results of the \link[R2CPCT]{importWGSOfCohort} function.
#' @param minTAF.Muts (double): Min. TAF for SNV, InDel, MNV in order to count.
#' @param minTAF.SV (double): Min. TAF for structural variants in order to count.
#'
#' @examples
#' \donttest{
#'
#'  data.Cohort <- R2CPCT::importWGSOfCohort(<cpctIds>, <combinedData>)
#'  determineMutationalBurden(data.Cohort)
#'
#' }
#' @return (tibble) Returns a tibble of the dNdS results.
#' @export
determineMutationalBurden <- function(data.Cohort, minTAF.Muts = 0, minTAF.SV = 0){

    # Input validation --------------------------------------------------------

    checkmate::assertClass(data.Cohort, classes = 'list')
    checkmate::assertClass(data.Cohort$somaticVariants, classes = 'VRanges')
    checkmate::assertClass(data.Cohort$structuralVariants, classes = 'GRanges')
    checkmate::assertClass(data.Cohort$copyNumbers, classes = 'GRanges')
    checkmate::assertClass(data.Cohort$purityStats, classes = 'tbl_df')
    checkmate::assertClass(data.Cohort$shatterSeek, classes = 'tbl_df')
    checkmate::assertDouble(minTAF.Muts)
    checkmate::assertDouble(minTAF.SV)

    sprintf('Determining mutational burden for %s samples.', dplyr::n_distinct(data.Cohort$purityStats$sample)) %>% ParallelLogger::logInfo()


    # Calculate TMB and number of muts ----------------------------------------

    # Retrieve the mutations of the cohort.
    mutData <- tibble::as_tibble(S4Vectors::mcols(data.Cohort$somaticVariants))

    # Determine number of SNV, InDel/DelIn and MNV per sample.
    perSample.Muts <- mutData %>% dplyr::group_by(sample, mutType) %>% dplyr::filter(PURPLE_AF >= minTAF.Muts) %>% dplyr::summarise(totalMut = n()) %>% dplyr::ungroup()

    # Determine Genome-wide TMB
    basesInFasta <- 2858674662 # Number of mappable ATCG in reference genome (hg19).
    perSample.TMB <- perSample.Muts %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(
            Genome.TMB = sum(totalMut) / (basesInFasta / 1E6),
            tmbStatus = base::ifelse(Genome.TMB >= 10, 'High TMB (\u226510)', base::ifelse(Genome.TMB >= 5, 'Medium TMB (\u22655-10)',  'Low TMB (0-5)'))
        )

    # Determine number of coding and intragenic mutations per sample.
    perSample.Muts.c <- mutData %>% dplyr::filter(!is.na(ANN.INTRON) | !is.na(ANN.HGVSp), ANN.BIOTYPE == 'protein_coding') %>% dplyr::group_by(sample) %>% dplyr::summarise(totalIntragenicMuts = n()) %>% dplyr::ungroup()
    perSample.Muts.p <- mutData %>% dplyr::filter(!is.na(ANN.HGVSp), ANN.BIOTYPE == 'protein_coding') %>% dplyr::group_by(sample)%>% dplyr::summarise(totalProteinCodingMuts = n()) %>% dplyr::ungroup()

    # Combine the data.
    perSample.Muts <- perSample.Muts %>%
        dplyr::mutate(mutType = paste0('Muts.', mutType)) %>%
        dplyr::mutate(sample = factor(sample, levels = levels(data.Cohort$purityStats$sample))) %>%
        tidyr::complete(sample, mutType) %>%
        tidyr::pivot_wider(names_from = mutType, values_from = totalMut) %>%
        dplyr::left_join(perSample.TMB) %>%
        dplyr::left_join(perSample.Muts.c) %>%
        dplyr::left_join(perSample.Muts.p)


    # Calculate number of SV --------------------------------------------------

    svData <- tibble::as_tibble(S4Vectors::mcols(data.Cohort$structuralVariants))

    perSamples.perSV <- svData %>% dplyr::filter(TAF >= minTAF.SV) %>% dplyr::group_by(sample, SVTYPE) %>% dplyr::filter(!base::duplicated(EVENT)) %>% dplyr::summarise(totalSV = n()) %>% dplyr::ungroup() %>% dplyr::mutate(SVTYPE = paste0('SV.', SVTYPE))
    perSamples.totalSV <- svData %>% dplyr::filter(TAF >= minTAF.SV) %>% dplyr::group_by(sample) %>% dplyr::filter(!base::duplicated(EVENT)) %>% dplyr::summarise(totalSV = n())

    perSamples.perSV <- perSamples.perSV %>%
        dplyr::mutate(sample = factor(sample, levels = levels(data.Cohort$purityStats$sample))) %>%
        tidyr::complete(sample, SVTYPE) %>%
        tidyr::pivot_wider(names_from = SVTYPE, values_from = totalSV) %>%
        dplyr::left_join(perSamples.totalSV)


    # Add mean genome-wide ploidy ---------------------------------------------

    perSample.ploidy <- data.Cohort$purityStats %>%
        dplyr::select(sample, genomePloidy = ploidy, msStatus, wholeGenomeDuplication)


    # Combine data ------------------------------------------------------------

    combinedData <- perSample.Muts %>%
        dplyr::left_join(perSamples.perSV) %>%
        dplyr::left_join(perSample.ploidy) %>%
        # Add CHORD status.
        dplyr::left_join(data.Cohort$CHORD %>% dplyr::select(sample, hr_status, hrd_type)) %>%
        # Add chromothripsis status.
        dplyr::mutate(hasChromothripsis = ifelse(sample %in% (data.Cohort$shatterSeek %>% dplyr::filter(isChromothripsis == 'Yes'))$sample, 'Yes', 'No'))

    # Return statement --------------------------------------------------------

    return(combinedData)

}
