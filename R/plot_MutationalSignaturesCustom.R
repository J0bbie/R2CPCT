#' @title Generates a mutational signature plot for custom signatures.
#'
#' @param mutSigs (list): A list containing the data of the mutational reconstruction (of either the SNV, InDel or DBS) from \code{\link{fitMutSigs}}.
#' @param orderSamples (character): Specific ordering of samples.
#' @param minContrib (numeric): Min. contribution of a signature (in any sample) in order to plot it. Set to 0 to plot all possible signatures.
#' @param sortOnMeanContrib (logical): Sort signature order in mean contribution (T) or as leave as supplied (F).
#' @param minMutants (numeric): Min. number of mutants (SNV, DBS or InDel) on which the fitting should have been performed.
#' Else, it will be flagged as a low-confident reconstruction.
#' @param focusOnSig (character): Should there be a focus on specific signatures? All others will be categorized as 'Other'.
#' @param onlyShowFocus (logical): Should only the focused signatures be shown?
#'
#' @return (ggplot) ggplot object.
#' @export
plotMutationalSignaturesCustom <- function(mutSigs, orderSamples = NULL, minContrib = 5, sortOnMeanContrib = TRUE, minMutants = NULL, focusOnSig = NULL, onlyShowFocus = FALSE){

    # Input validation --------------------------------------------------------

    checkmate::assertList(mutSigs)
    checkmate::assertCharacter(orderSamples, null.ok = TRUE)
    checkmate::assertNumber(minContrib)
    checkmate::assertLogical(sortOnMeanContrib)
    checkmate::assertNumber(minMutants, null.ok = TRUE)
    checkmate::assertCharacter(focusOnSig, null.ok = TRUE)
    checkmate::assertLogical(onlyShowFocus)

    sprintf('\tPlotting Mutational Signatures (Custom).') %>% ParallelLogger::logInfo()


    # Convert data for plotting -----------------------------------------------

    # Retrieve the rel. contribution of signatures.
    mutSigs.Rel <- mutSigs$relativeContribution

    if(!is.null(focusOnSig)){

        mutSigs.Rel <- mutSigs.Rel %>%
            dplyr::mutate(Signature = ifelse(Signature %in% focusOnSig, Signature, 'Other signatures')) %>%
            dplyr::group_by(sampleId, Signature) %>%
            dplyr::summarise(
                relContribution = base::sum(relContribution),
                Signature = unique(Signature)
            ) %>%
            dplyr::ungroup()

    }

    # Filter on min. contribution.
    mutSigs.Rel <- mutSigs.Rel %>%
        # Check which signatures have a rel. contribution below the threshold in all samples.
        dplyr::group_by(sampleId, Signature) %>%
        dplyr::mutate(
            isBelowFilter = all(relContribution < minContrib),
        ) %>%
        dplyr::ungroup() %>%
        # Group filtered signatures together.
        dplyr::group_by(sampleId, Signature, isBelowFilter) %>%
        dplyr::mutate(
            Signature = ifelse(isBelowFilter, sprintf('Filtered (<%s%%)', minContrib), Signature),
        ) %>% dplyr::ungroup() %>%
        dplyr::group_by(sampleId, Signature) %>%
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
            dplyr::filter(Signature != 'Other signatures') %>%
            dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
            tidyr::complete(sampleId, Signature, fill = list(relContribution = 0))

    }


    # Generate plot -----------------------------------------------------------

    if(sortOnMeanContrib){
        orderOfSigs <- mutSigs$relativeContribution %>% dplyr::group_by(Signature) %>% dplyr::summarise(meanContrib = mean(relContribution)) %>% dplyr::arrange(-meanContrib) %>% dplyr::pull(Signature)
        orderOfSigs <- c(orderOfSigs, 'Other signatures', sprintf('Filtered (<%s%%)', minContrib))
    }else{
        orderOfSigs <- c(unique(mutSigs$relativeContribution$Signature), 'Other signatures', sprintf('Filtered (<%s%%)', minContrib))
    }

    # Sort order of signatures.
    mutSigs.Rel <- mutSigs.Rel %>% dplyr::mutate(Signature = factor(Signature, levels = orderOfSigs))

    # Colors to use.
    colors <- c('#fe4a49', '#2ab7ca', '#fe9c8f', '#005b96', '#139c3c', '#EBB584', '#f9b320', '#A088C3', '#d1cc71', '#91E4A6', '#785EF0', '#44AA99', '#CC6677', '#D869C0', '#68564E', '#9fbfc4', '#7ED7D1', '#525975', '#88CCEE', '#FFB6C1', '#B10DA1', '#DEA0FD', '#90AD1C', '#FFD599', '#FF6136', '#A98AAD', '#9BEBEE', '#2B580D', '#a6143b', '#9c679b', '#76b078', '#d49919')
    colors <- c(colors[1:(dplyr::n_distinct(orderOfSigs) - 2)], 'grey75', 'grey90')
    names(colors) <- orderOfSigs


    # Generate the plot.
    plot <- ggplot2::ggplot(mutSigs.Rel, ggplot2::aes(x = sampleId, y = relContribution / 100, fill = Signature)) +
        ggplot2::geom_bar(stat = 'identity', lwd = .33, color = 'black', width = .8) +
        ggplot2::labs(y = 'Mut. Signatures<br><span style = "font-size:5pt">(Genome-wide)</span>', x = NULL) +
        ggplot2::scale_y_continuous(expand = c(0,0), labels = scales::percent, limits = c(0, ifelse(onlyShowFocus, ceiling(max(mutSigs.Rel$relContribution*2)) / 2 / 100, 1.0001))) +
        ggplot2::scale_fill_manual(values = colors, guide = ggplot2::guide_legend(title = 'Mutational Signatures (custom)', title.position = 'top', title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5), name = NULL) +
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