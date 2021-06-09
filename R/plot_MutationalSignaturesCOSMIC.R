#' @title Generates a mutational signature plot (COSMIC signatures).
#'
#' @param mutSigs (list): A list containing the data of the mutational reconstruction (of either the SNV, InDel or DBS) from \code{\link{fitMutSigs}}.
#' @param orderSamples (character): Specific ordering of samples.
#' @param minContrib (numeric): Min. contribution of a signature (in any sample) in order to plot it. Set to 0 to plot all possible signatures.
#' @param minMutants (numeric): Min. number of mutants (SNV, DBS or InDel) on which the fitting should have been performed.
#' Else, it will be flagged as a low-confident reconstruction.
#' @param combineSigs (logical): Combine signatures based on aetiology?
#' @param focusOnSig (character): Should there be a focus on specific signatures? All others will be categorized as 'Other'.
#' @param onlyShowFocus (logical): Should only the focused signatures be shown?
#'
#' @return (ggplot) ggplot object.
#' @export
plotMutationalSignaturesCOSMIC <- function(mutSigs, orderSamples = NULL, minContrib = 5, minMutants = NULL, combineSigs = T, focusOnSig = NULL, onlyShowFocus = F){

    # Input validation --------------------------------------------------------

    checkmate::assertList(mutSigs)
    checkmate::assertCharacter(orderSamples, null.ok = T)
    checkmate::assertNumber(minContrib)
    checkmate::assertNumber(minMutants, null.ok = T)
    checkmate::assertLogical(combineSigs)
    checkmate::assertCharacter(focusOnSig, null.ok = T)
    checkmate::assertLogical(onlyShowFocus)

    sprintf('\tPlotting Mutational Signatures (COSMIC).') %>% ParallelLogger::logInfo()


    # Convert data for plotting -----------------------------------------------

    # Retrieve the rel. contribution of signatures.
    mutSigs.Rel <- mutSigs$relativeContribution

    # Combine COSMIC signatures by proposed aetiology.
    if(combineSigs){

        mutSigs.Rel <- mutSigs.Rel %>%
            dplyr::group_by(sampleId, proposedAetiologyGrouped) %>%
            dplyr::summarise(
                relContribution = base::sum(relContribution),
                proposedAetiology = unique(proposedAetiologyGrouped)
            ) %>%
            dplyr::ungroup()

    }

    if(!is.null(focusOnSig)){

        mutSigs.Rel <- mutSigs.Rel %>%
            dplyr::mutate(proposedAetiology = ifelse(proposedAetiology %in% focusOnSig, proposedAetiology, 'Other signatures')) %>%
            dplyr::group_by(sampleId, proposedAetiology) %>%
            dplyr::summarise(
                relContribution = base::sum(relContribution),
                proposedAetiology = unique(proposedAetiology)
            ) %>%
            dplyr::ungroup()

    }

    # Filter on min. contribution.
    mutSigs.Rel <- mutSigs.Rel %>%
        # Check which signatures have a rel. contribution below the threshold in all samples.
        dplyr::group_by(sampleId, proposedAetiology) %>%
        dplyr::mutate(
            isBelowFilter = all(relContribution < minContrib),
        ) %>%
        dplyr::ungroup() %>%
        # Group filtered signatures together.
        dplyr::group_by(sampleId, proposedAetiology, isBelowFilter) %>%
        dplyr::mutate(
            proposedAetiology = ifelse(isBelowFilter, sprintf('Filtered (<%s%%)', minContrib), proposedAetiology),
        ) %>% dplyr::ungroup() %>%
        dplyr::group_by(sampleId, proposedAetiology) %>%
        dplyr::summarise(
            relContribution = base::sum(relContribution),
        ) %>%
        dplyr::ungroup()

    # Sort samples.
    if(!is.null(orderSamples)) mutSigs.Rel <- mutSigs.Rel %>% dplyr::filter(sampleId %in% orderSamples) %>% dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples))


    # Determine low-mutant samples --------------------------------------------

    if(!is.null(minMutants)){
        mutSigs.Rel <- mutSigs.Rel %>%
            dplyr::mutate(
                lowConfidence = ifelse(sampleId %in% base::colnames(mutSigs$mutMatrix)[base::colSums(mutSigs$mutMatrix) < minMutants], '*', '')
            )
    }

    # Only show the focused signatures.
    if(onlyShowFocus & !is.null(focusOnSig)){

        mutSigs.Rel <- mutSigs.Rel %>%
            dplyr::filter(proposedAetiology != 'Other signatures') %>%
            dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
            tidyr::complete(sampleId, proposedAetiology, fill = list(relContribution = 0))

    }


    # Generate plot -----------------------------------------------------------

    # Color to use.
    if(combineSigs){
        colors <- c(
            'Spontaneous deamination of 5-methylcytosine (clock-like signature) (SBS1)' = '#85660D',
            'Activity of APOBEC family of cytidine deaminases (SBS2, SBS13)' = '#c11d00',
            'Defective homologous recombination DNA damage repair (SBS3, ID6)' = '#0079c2',
            'Tobacco smoking (SBS4, ID3)' = '#16FF32',
            'Tobacco smoking and other mutagens (e.g. acetaldehyde) (DBS2)' = '#139c3c',
            'Unknown (clock-like signature) (SBS5)' = '#EBB584',
            'Defective DNA mismatch repair (>4 signatures)' = '#f9b320',
            'Ultraviolet light exposure (>4 signatures)' = 'yellow',
            'Polimerase eta somatic hypermutation activity (SBS9)' = '#A088C3',
            'Polymerase epsilon exonuclease domain mutations (SBS10a, SBS10b, DBS3)' = '#d1cc71',
            'Temozolomide treatment (SBS11)' = '#91E4A6',
            'Concurrent polymerase epsilon mutation and defective DNA mismatch repair (SBS14)' = '#785EF0',
            'Damage by reactive oxygen species (SBS18)' = '#44AA99',
            'Concurrent POLD1 mutations and defective DNA mismatch repair (SBS20)' = '#CC6677',
            'Aristolochic acid exposure (SBS22)' = '#D869C0',
            'Aflatoxin exposure (SBS24)' = '#68564E',
            'Chemotherapy treatment (SBS25)' = '#9fbfc4',
            'Tobacco chewing (SBS29)' = '#7ED7D1',
            'Defective DNA base excision repair due to NTHL1 mutations (SBS30)' = '#525975',
            'Platinum chemotherapy treatment (SBS31, SBS35, DBS5)' = '#88CCEE',
            'Azathioprine treatment (SBS32)' = '#FFB6C1',
            'Defective DNA base excision repair due to MUTYH mutations (SBS36)' = '#B10DA1',
            'Indirect effect of ultraviolet light (SBS38)' = '#DEA0FD',
            'Haloalkane exposure (SBS42)' = '#90AD1C',
            'Activity of activation-induced cytidine deaminase (AID) (SBS84)' = '#FFD599',
            'Indirect effects of activation-induced cytidine deaminase (AID) (SBS85)' = '#FF6136',
            'Unknown chemotherapy treatment (SBS86)' = 'skyblue',
            'Thiopurine chemotherapy treatment (SBS87)'  = '#A98AAD',
            'Colibactin exposure (E.coli bacteria carrying pks pathogenicity island) (SBS88, ID18)' = '#9BEBEE',
            'Duocarmycin exposure (SBS90)' = '#2B580D',
            'Unknown (possibly related to APOBEC mutagenesis) (DBS11)' = '#a6143b',
            'Slippage during DNA replication of the replicated DNA strand (enriched in cancers with DNA mismatch repair deficiency) (ID1, ID2)'  = '#9c679b',
            'Repair of DNA double strand breaks by non-homologous end-joining mechanisms or mutations in topoisomerase TOP2A (ID8)' = '#76b078',
            'Mutations in topoisomerase TOP2A (ID17)' = '#d49919',
            'Possible Sequencing Artifact (>4 signatures)' = 'grey70',
            'Unknown (>4 signatures)'  = 'grey90',
            'Filtered' = 'grey75',
            'Other signatures' = 'grey75'
        )

        names(colors)[base::length(colors)] <- sprintf('Filtered (<%s%%)', minContrib)

    }else{
        colors <- hues::iwanthue(dplyr::n_distinct(mutSigs.Rel$proposedAetiology))
        names(colors) <- unique(mutSigs.Rel$proposedAetiology)

        colors[names(colors) == 'Other signatures'] <- 'grey75'

    }

    # Sort on signatures (decreasing SBS / DBS numbers).
    mutSigs.Rel <- mutSigs.Rel %>% dplyr::mutate(proposedAetiology = factor(proposedAetiology, levels = names(colors)))

    # Generate the plot.
    plot <- ggplot2::ggplot(mutSigs.Rel, ggplot2::aes(x = sampleId, y = relContribution / 100, fill = proposedAetiology)) +
        ggplot2::geom_bar(stat = 'identity', lwd = .33, color = 'black', width = .8) +
        ggplot2::labs(y = 'Mut. Signatures<br><span style = "font-size:5pt">(Genome-wide)</span>', x = NULL) +
        ggplot2::scale_y_continuous(expand = c(0,0), labels = scales::percent, limits = c(0, ifelse(onlyShowFocus, ceiling(max(mutSigs.Rel$relContribution*2)) / 2 / 100, 1.0001))) +
        ggplot2::scale_fill_manual(values = colors, guide = ggplot2::guide_legend(title = 'Mutational Signatures (COSMIC v3.1)', title.position = 'top', title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5), name = NULL) +
        ggplot2::theme(
            legend.position = 'right',
            legend.direction = 'horizontal',
            text = ggplot2::element_text(size = 7, family = 'Helvetica', face = 'bold'),
            legend.text = ggplot2::element_text(size = 7, family='Helvetica', face = 'plain'),
            axis.text.x = ggplot2::element_blank(),
            axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
            axis.ticks.x = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dashed'),
            panel.grid.minor.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            panel.border = ggplot2::element_rect(fill = NA, colour = NA),
            strip.background = ggplot2::element_rect(colour = 'grey20', fill = 'white')
        )

    if(!is.null(minMutants)) plot <- plot + ggplot2::geom_text(data = mutSigs.Rel %>% dplyr::distinct(sampleId, lowConfidence), ggplot2::aes(y = 0, fill = NULL, label = lowConfidence), size = 5, color = '#EB6453', nudge_y = -.1)


    # Return statement --------------------------------------------------------

    return(plot)

}