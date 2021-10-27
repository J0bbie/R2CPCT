#' @title Generates a mutational signature plot.
#'
#' @param dataOncoplot (tibble): Mutational data from the combinedReport.
#' @param includedSamples (character): Which samples should be included (completes missing samples)?
#' @param sortSamples (character): Specific order of samples ordering (x-axis)?
#' @param sortGenes (character): Specific order of genes (y-axis)?
#' @param showAlleles (logical): Should the mutations be (roughly) separated per paternal/maternal allele (TRUE) or combined (FALSE).
#'
#' @return (ggplot) ggplot object.
#' @export
makeOncoplot.CPCT <- function(dataOncoplot, includedSamples = NULL, sortSamples = NULL, sortGenes = NULL, showAlleles = TRUE){

    # Input validation --------------------------------------------------------

    checkmate::assertClass(dataOncoplot, class = 'tbl_df')
    checkmate::assertCharacter(includedSamples, null.ok = TRUE)
    checkmate::assertCharacter(sortSamples, null.ok = TRUE)
    checkmate::assertCharacter(sortGenes, null.ok = TRUE)
    checkmate::assertLogical(showAlleles, null.ok = TRUE)


    # Generalize consequences -------------------------------------------------

    dataOncoplot <- dataOncoplot %>%
        dplyr::mutate(
            Consequence.Mut.Clean = base::ifelse(base::grepl('Splice', Consequence.Mut), 'Splicing variant', Consequence.Mut),
            Consequence.Mut.Clean = base::ifelse(base::grepl('Stop|Start', Consequence.Mut.Clean), 'Stop/Gain variant', Consequence.Mut.Clean),
            Consequence.Mut.Clean = base::ifelse(base::grepl('Inframe', Consequence.Mut.Clean), 'Inframe insertion/deletion', Consequence.Mut.Clean)
        )

    # Determine structural variants.
    dataOncoplot <- dataOncoplot %>%
        dplyr::group_by(sample, SYMBOL, ENSEMBL) %>%
        dplyr::mutate(
            Consequence.Mut.Clean = base::ifelse(!is.na(Consequence.SV), base::ifelse(!is.na(Consequence.Mut.Clean), 'Multiple coding mutations', 'Structural variant'), Consequence.Mut.Clean)
        ) %>%
        dplyr::ungroup()


    # Complete missing samples ------------------------------------------------

    if(!is.null(includedSamples)){
        dataOncoplot <- dataOncoplot %>%
            dplyr::filter(sample %in% includedSamples) %>%
            dplyr::mutate(sample = factor(sample, levels = includedSamples)) %>%
            tidyr::complete(sample, SYMBOL, fill = list(isMutant = FALSE)) %>%
            dplyr::distinct()
    }


    # memoSort ----------------------------------------------------------------

    memoData <- reshape2::dcast(dataOncoplot, SYMBOL ~ sample, fill = 0, value.var = 'isMutant', fun.aggregate = sum)
    rownames(memoData) <- memoData$SYMBOL; memoData$SYMBOL <- NULL
    memoData[is.na(memoData)] <- 0
    memoData <- R2CPCT::memoSort(memoData)

    dataOncoplot$SYMBOL <- base::factor(dataOncoplot$SYMBOL, levels = base::rev(base::rownames(memoData)))
    dataOncoplot$sample <- base::factor(dataOncoplot$sample, levels = base::colnames(memoData))


    # Override sample / gene order --------------------------------------------

    if(!is.null(sortGenes)) dataOncoplot <- dataOncoplot %>% dplyr::mutate(SYMBOL = base::factor(SYMBOL, levels = sortGenes))
    if(!is.null(sortSamples)) dataOncoplot <- dataOncoplot %>% dplyr::mutate(sample = base::factor(sample, levels = sortSamples))


    # Calculate mutational frequencies per gene -------------------------------

    dataOncoplot.MutFrequencies <- dataOncoplot %>%
        dplyr::filter(isMutant) %>%
        dplyr::group_by(SYMBOL) %>%
        dplyr::summarise(
            'Splicing variants' = base::sum(Consequence.Mut.Clean == 'Splicing variant', na.rm = TRUE),
            'Structural variants' = base::sum(!is.na(Consequence.SV), na.rm = TRUE),
            'Coding variants' = base::sum(!is.na(Consequence.Mut), na.rm = TRUE) - `Splicing variants`,

            'Deep amplifications' = sum(Consequence.CNA == 'Deep Amplification', na.rm = TRUE),
            'Deep deletions' = sum(Consequence.CNA == 'Deep Deletion', na.rm = TRUE),

            totalMuts = `Coding variants` + `Structural variants` + `Splicing variants` + `Deep amplifications` + `Deep deletions`,
            mutSamples = dplyr::n_distinct(sample),
            mutSamples.Rel = mutSamples / dplyr::n_distinct(dataOncoplot$sample)

        ) %>%
        dplyr::ungroup() %>%
        reshape2::melt(id.vars = c('SYMBOL', 'mutSamples', 'mutSamples.Rel', 'totalMuts')) %>%
        dplyr::mutate(
            value.Rel = value / totalMuts,
            SYMBOL = factor(SYMBOL, levels = levels(dataOncoplot$SYMBOL)),
            variable = factor(variable, levels = rev(c('Coding variants', 'Splicing variants', 'Structural variants', 'Deep amplifications', 'Deep deletions')))
        )


    # Determine dN/dS or GISTIC2 evidence -------------------------------------

    dataOncoplot.Evidence <- dataOncoplot %>%
        dplyr::distinct(SYMBOL, dNdS, GISTIC2Peak) %>%
        dplyr::mutate(
            dNdS = ifelse(!is.na(dNdS), 'dN/dS', NA),
            GISTIC2Peak = ifelse(grepl('Amplification', GISTIC2Peak), 'GISTIC2 - Amplification', GISTIC2Peak),
            GISTIC2Peak = ifelse(grepl('Deletion', GISTIC2Peak), 'GISTIC2 - Deletion', GISTIC2Peak),
            SYMBOL = factor(SYMBOL, levels = levels(dataOncoplot$SYMBOL))
        ) %>%
        reshape2::melt(id.vars = 'SYMBOL')


    # Determine allelic consequences ------------------------------------------

    dataOncoplot <- dataOncoplot %>%
        dplyr::mutate(
            Consequence.AlleleA.Mut = 'Neutral',
            Consequence.AlleleB.Mut = 'Neutral',
            Consequence.AlleleA.CNA = 'Neutral',
            Consequence.AlleleB.CNA = 'Neutral',

            # Set Allele A as mutated.
            Consequence.AlleleA.Mut = ifelse(!is.na(Consequence.Mut.Clean), Consequence.Mut.Clean, Consequence.AlleleA.Mut),
            Consequence.AlleleB.Mut = ifelse(Consequence.Mut.Clean == 'Multiple coding mutations', Consequence.Mut.Clean, Consequence.AlleleB.Mut),

            Consequence.AlleleB.CNA = ifelse(Consequence.LOH == 'LOH' | grepl('Dele', Consequence.CNA), 'Heterozygous deletion / LOH', Consequence.AlleleB.CNA),
            Consequence.AlleleA.CNA = ifelse(Consequence.CNA == 'Deep Deletion', 'Homozygous deletion', Consequence.AlleleA.CNA),
            Consequence.AlleleB.CNA = ifelse(Consequence.CNA == 'Deep Deletion', 'Homozygous deletion', Consequence.AlleleB.CNA),


            Consequence.AlleleB.CNA = ifelse(grepl('Amp', Consequence.CNA), 'Amplification', Consequence.AlleleB.CNA),
            Consequence.AlleleA.CNA = ifelse(Consequence.CNA == 'Deep Amplification', 'Deep Amplification', Consequence.AlleleA.CNA),
            Consequence.AlleleB.CNA = ifelse(Consequence.CNA == 'Deep Amplification', 'Deep Amplification', Consequence.AlleleB.CNA)

        )

    # Generate figures --------------------------------------------------------

    # List to keep figures and data.
    figOncoplot <- list()

    # Add transformed oncoplot data.
    figOncoplot$dataOncoplot <- dataOncoplot


    ## Color of the mutations. ----
    colorCNA <- c(
        'Homozygous deletion' = '#223A64',
        'Heterozygous deletion / LOH' = '#94AAD0',
        'Deep Deletion' = '#223A64',
        'Deletion' = '#94AAD0',
        'Neutral' = '#e6e6e6',
        'Amplification' = '#6EC41F',
        'Deep Amplification' = '#177E63'
    )

    colorMuts <- c(
        'Synonymous variant' = 'yellow',
        'Frameshift variant' = 'purple',
        'Inframe insertion/deletion' = '#F341A2', # Pink
        'Stop/Gain variant' = '#4897D8', # Blue
        'Missense variant' = '#000000', # Black
        'Splicing variant' ='#F9A603', # Orange
        'Multiple coding mutations' = 'red', # Red
        'Structural variant' = 'skyblue', # skyblue
        'Initiator codon variant' = 'lightgreen', # Lightgreen
        'Protein altering variant' = '#F0A979'
    )

    ## Figure - Oncoplot ----

    if(showAlleles){
        figOncoplot$Oncoplot <- ggplot2::ggplot(dataOncoplot, aes(x = sample, y = SYMBOL)) +
            ggplot2::geom_tile(color = '#e6e6e6', lwd = .5, fill = 'white', height = .8, width = .7, na.rm = TRUE) +

            ggplot2::geom_tile(aes(color = Consequence.AlleleA.CNA), fill = 'white', lwd = .4, na.rm = TRUE, alpha = 1, height = .4, width = .7, position = ggplot2::position_nudge(y = +.2)) +
            ggplot2::geom_tile(aes(color = Consequence.AlleleB.CNA), fill = 'white', lwd = .4, na.rm = TRUE, alpha = 1, height = .4, width = .7, position = ggplot2::position_nudge(y = -.2)) +

            ggplot2::geom_tile(aes(fill = Consequence.AlleleA.Mut), lwd = 0, na.rm = TRUE, alpha = 1, height = .2, width = .4, position = ggplot2::position_nudge(y = +.2)) +
            ggplot2::geom_tile(aes(fill = Consequence.AlleleB.Mut), lwd = 0, na.rm = TRUE, alpha = 1, height = .2, width = .4, position = ggplot2::position_nudge(y = -.2)) +

            # Colors of mutations.
            ggplot2::scale_fill_manual(values = colorMuts, drop = TRUE) +
            ggplot2::scale_color_manual(values = colorCNA, na.value = 'grey95', drop = TRUE) +
            ggplot2::scale_x_discrete(expand = c(.005,.01)) +
            ggplot2::scale_y_discrete(expand = c(.005,.01)) +
            ggplot2::labs(x = sprintf('Samples <sub>(<i>n</i> = %s)</sub>', dplyr::n_distinct(dataOncoplot$sample)), y = 'Genes') +
            # Legend settings.
            ggplot2::guides(color = guide_legend(title = 'Copynumber Categories', title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +

            ggplot2::guides(fill = guide_legend(title = 'Mutational Categories', title.position = 'top', title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +
            ggplot2::theme(
                legend.position = 'bottom',
                legend.direction = 'horizontal',
                text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
                axis.text.x = ggplot2::element_blank(),
                axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
                axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
                panel.grid = ggplot2::element_blank(),
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.grid.minor.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                panel.border = ggplot2::element_rect(fill = NA, colour = NA),
                strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
                legend.text = ggtext::element_markdown(),
                legend.title = ggtext::element_markdown(),
                legend.key = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank()
            )
    }else{

        figOncoplot$Oncoplot <- ggplot2::ggplot(dataOncoplot, aes(x = sample, y = SYMBOL)) +
            ggplot2::geom_tile(aes(fill = Consequence.CNA), lwd = .2, color = 'grey90', na.rm = TRUE) +
            ggplot2::geom_tile(aes(fill = Consequence.Mut.Clean), lwd = .2, na.rm = TRUE, alpha = 1, height = .5, width = .5) +
            # Colors of mutations.
            ggplot2::scale_fill_manual(values = colorMuts, drop = TRUE) +
            ggplot2::labs(x = sprintf('Samples <sub>(<i>n</i> = %s)</sub>', dplyr::n_distinct(dataOncoplot$sample)), y = 'Genes') +
            # Legend settings.
            ggplot2::guides( fill = guide_legend(title = 'Mutational Categories', title.position = 'top', title.hjust = 0.5, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
            ggplot2::theme(
                legend.position = 'bottom',
                legend.direction = 'horizontal',
                text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
                axis.text.x = ggplot2::element_blank(),
                axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
                axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
                strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
                panel.grid = ggplot2::element_blank(),
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.grid.minor.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
                panel.border = ggplot2::element_rect(fill = NA, colour = NA),
                strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
                legend.text = ggtext::element_markdown(),
                legend.title = ggtext::element_markdown(),
                legend.key = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank()
            )
    }


    ## Figure - Mutational frequencies ----
    figOncoplot$Oncoplot.MutFrequencies <- ggplot2::ggplot(dataOncoplot.MutFrequencies, aes(x = SYMBOL, fill = variable, y = value.Rel)) +
        ggplot2::geom_bar(stat = 'identity', lwd = .33, color = 'black', width = .7) +
        ggplot2::scale_x_discrete(expand = c(0,0)) +
        ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent, expand = c(0,0)) +
        ggplot2::labs(x = NULL, y = 'Mutational category') +
        ggplot2::scale_fill_manual(NULL, values = c('Coding variants' = '#E1B119', 'Splicing variants' = 'pink', 'Structural variants' = 'skyblue', 'Deep amplifications' = '#2EAA59', 'Deep deletions' = '#D4313E'), guide = guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
        ggplot2::coord_flip() +
        ggplot2::theme(
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
            axis.text.x = ggplot2::element_text(size = 6, hjust = .5),
            axis.text.y = ggplot2::element_blank(),
            axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
            axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
            strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            panel.border = ggplot2::element_rect(fill = NA, colour = NA),
            strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
            legend.text = ggtext::element_markdown(),
            legend.title = ggtext::element_markdown(),
            legend.key = ggplot2::element_blank()
        )

    ## Figure - Evidence ----
    figOncoplot$Oncoplot.Evidence <- ggplot2::ggplot(dataOncoplot.Evidence, aes(x = 'Evidence', y = SYMBOL, fill = value)) +
        ggplot2::geom_tile(data = dataOncoplot.Evidence %>% dplyr::filter(variable == 'GISTIC2Peak'), color = 'grey95', size = 0.25, na.rm = TRUE, height = .9) +
        ggplot2::geom_tile(data = dataOncoplot.Evidence %>% dplyr::filter(variable == 'dNdS'), size = 0.25, na.rm = TRUE, alpha = 1, height = .5, width = .5) +
        ggplot2::scale_x_discrete(expand = c(0,0)) +
        ggplot2::scale_y_discrete(expand = c(0,0)) +
        # Colors of mutations.
        ggplot2::scale_fill_manual(NULL, values = c('dN/dS' = 'black', 'GISTIC2 - Amplification' = '#2EAA59', 'GISTIC2 - Deletion' = 'red'), guide = ggplot2::guide_legend(ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
        ggplot2::labs(x = 'Detected by dN/dS and/or GISTIC2', y = NULL) +
        # Change theme.
        ggplot2::theme(
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
            axis.text.x = ggplot2::element_text(size = 6, hjust = .5),
            axis.text.y = ggplot2::element_blank(),
            axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
            axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
            strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            panel.border = ggplot2::element_rect(fill = NA, colour = NA),
            strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
            legend.text = ggtext::element_markdown(),
            legend.title = ggtext::element_markdown(),
            legend.key = ggplot2::element_blank()
        )


    # Return statement --------------------------------------------------------

    return(figOncoplot)

}
