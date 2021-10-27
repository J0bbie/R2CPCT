#' @title Perform a Fisher's Exact test to detect groups enriched with mutants (mutations / deep CNA).
#'
#' @param mutData (tibble): Combinedreport containing the somatic mutations per gene, per sample.
#' @param groupInfo (tibble): Tibble containing grouping information. Needs the following columns: sample, group.
#' @param genesOfInterest (character): List of genes (ENSEMBL/SYMBOL) for which to always look for enrichment (regardless of minRel mutants.)
#' @param minRelMutantsInGroup (numeric): Only check enrichment for genes with >=X rel. mutants in the group.
#'
#' @return (tibble) Returns a tibble of the results from the Fisher's Exact test.
#' @export
detectEnrichmentOfMutants <- function(mutData, groupInfo, genesOfInterest, minRelMutantsInGroup = .1){


    # Input validation --------------------------------------------------------

    checkmate::assertClass(mutData, classes = 'tbl_df')
    checkmate::assertClass(groupInfo, classes = 'tbl_df')
    checkmate::assertTRUE(all(c('sample', 'group') %in% colnames(groupInfo)))

    checkmate::assertCharacter(genesOfInterest)
    checkmate::assertNumber(minRelMutantsInGroup)

    # Logging.
    sprintf('Testing enrichment on %s groups.', dplyr::n_distinct(groupInfo$group)) %>% ParallelLogger::logInfo()


    # Prepate data and groups -------------------------------------------------

    # Determine size of each group.
    groupInfo <- groupInfo %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(totalInGroup = dplyr::n_distinct(sample)) %>%
        dplyr::ungroup()

    ## Min. Muts. (rel.) per group ----
    # They should also not be flagged within FLAGS.
    data.FLAGS <- readr::read_delim(system.file('extdata', 'largeGenes_FLAGS.txt', package = 'R2CPCT'), delim = '\t')

    mutData <- mutData %>%
        # Select relevant samples and genes.
        dplyr::filter(
            isMutant,
            !SYMBOL %in% data.FLAGS$SYMBOL,
            (!is.na(DriverDatabases) | SYMBOL %in% genesOfInterest | ENSEMBL %in% genesOfInterest)
        ) %>%
        # Add grouping info.
        dplyr::inner_join(groupInfo) %>%
        dplyr::group_by(SYMBOL, group) %>%
        dplyr::summarise(
            totalMut = dplyr::n_distinct(sample),
            noMut = totalInGroup - totalMut,
            totalMut.Rel = totalMut / unique(totalInGroup)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct()

    # Select only genes for which a group has >X% mutant samples.
    # Unless users have specified that gene to be looked at.
    discoveryGenes <- mutData %>% dplyr::filter(totalMut.Rel >= minRelMutantsInGroup | SYMBOL %in% genesOfInterest) %>% dplyr::distinct(SYMBOL)

    # Count nr. of mutants per group and gene.
    mutData <- mutData %>%
        dplyr::filter(SYMBOL %in% discoveryGenes$SYMBOL)

    # Complete missing data.
    mutData <- mutData %>%
        tidyr::complete(SYMBOL, group) %>%
        dplyr::group_by(SYMBOL, group) %>%
        dplyr::mutate(
            totalMut = ifelse(is.na(totalMut), 0, totalMut),
            noMut = ifelse(is.na(noMut), unique(groupInfo[groupInfo$group == unique(group),]$totalInGroup), noMut)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::distinct()


    # Perform Fisher's Exact Test per group versus all other groups. ----------

    mutGroups <- pbapply::pblapply(unique(mutData$group), function(k){

        # Perform Fisher Exact per gene.
        base::do.call(base::rbind, base::lapply(unique(mutData$SYMBOL), function(gene){

            geneData <- mutData %>%
                dplyr::filter(SYMBOL == gene) %>%
                dplyr::summarise(
                    k.withMut = totalMut[group == k],
                    k.withoutMut = noMut[group == k],

                    other.withMut = sum(totalMut[group != k]),
                    other.withoutMut = sum(noMut[group != k]),
                )

            test <- data.frame(
                row.names = c(k, 'Other'),
                mut = c(geneData$k.withMut, geneData$other.withMut),
                noMut = c(geneData$k.withoutMut, geneData$other.withoutMut)
            )

            test.df <- data.frame(
                SYMBOL = gene,
                p = stats::fisher.test(test, hybrid = FALSE, alternative = 'greater', simulate.p.value = TRUE)$p.value,
                group = k,
                mutInGroup = test[1,1],
                noMutInGroup = test[1,2],
                mutInOthers = test[2,1],
                noMutInOthers = test[2,2])

            return(test.df)
        }))
    }, cl = 10)

    # Merge results from >2 groups together
    mutGroups <- base::do.call(base::rbind, mutGroups)

    # Check direction and effect size.
    mutGroups <- mutGroups %>% dplyr::mutate(
        effectSize.A = round((mutInGroup / (mutInGroup + noMutInGroup)) * 100, 1),
        effectSize.B = round((mutInOthers / (mutInOthers + noMutInOthers)) * 100, 1),
        effectSize = sprintf('%s%% vs. %s%%', effectSize.A, effectSize.B),
        Direction = ifelse(effectSize.A > effectSize.B, 'Enriched', 'Depleted')
    )

    # Only check for enrichment.
    mutGroups <- mutGroups %>% dplyr::filter(Direction == 'Enriched')

    # Correct for multiple-testing.
    mutGroups$p.adj <- stats::p.adjust(mutGroups$p, method = 'BH')

    # Return statement.
    return(mutGroups)

}