#' @title Performs ShatterSeek analysis.
#'
#' @param copyNumbers (GRanges): VRanges of a single sample from \link[R2CPCT]{importSomaticVariantsVEP}.
#' @param structuralVariants (GRanges): GRanges of a single sample from \link[R2CPCT]{importStructuralVariantsPURPLE}.
#' @return (list) Tibble with the CHORD results.
#' @examples
#' \donttest{
#'
#' 	runCHORD(dataSample$somaticVariants, dataSample$structuralVariants)
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
performShatterSeek <- function(copyNumbers, structuralVariants){

    # Input validation --------------------------------------------------------

    checkmate::assertClass(copyNumbers, 'GRanges')
    checkmate::assertClass(structuralVariants, 'GRanges')

    sprintf('Running ShatterSeek (chromothripsis prediction) for: %s', unique(copyNumbers$sample)) %>% ParallelLogger::logInfo()


    # Convert SV to ShatterSeek input -----------------------------------------

    # Get SV
    structuralVariants.ShatterSeek <- structuralVariants
    structuralVariants.ShatterSeek <- structuralVariants.ShatterSeek[!base::duplicated(structuralVariants.ShatterSeek$EVENT) & structuralVariants.ShatterSeek$SVTYPE != 'SINGLE',]

    # Remove chrY.
    structuralVariants.ShatterSeek <- structuralVariants.ShatterSeek[!base::grepl('Y:', structuralVariants.ShatterSeek$REF) & !base::grepl('Y:', structuralVariants.ShatterSeek$ALT),]

    # Sort the SV.
    structuralVariants.ShatterSeek <- GenomicRanges::sort(structuralVariants.ShatterSeek)

    # Separate INV3 and INV5 SV.
    if(length(structuralVariants.ShatterSeek[structuralVariants.ShatterSeek$SVTYPE == 'INV' & grepl('\\[', structuralVariants.ShatterSeek$ALT)]) > 0) structuralVariants.ShatterSeek[structuralVariants.ShatterSeek$SVTYPE == 'INV' & grepl('\\[', structuralVariants.ShatterSeek$ALT)]$SVTYPE <- 't2tINV'
    if(length(structuralVariants.ShatterSeek[structuralVariants.ShatterSeek$SVTYPE == 'INV' & grepl('\\]', structuralVariants.ShatterSeek$ALT)]) > 0) structuralVariants.ShatterSeek[structuralVariants.ShatterSeek$SVTYPE == 'INV' & grepl('\\]', structuralVariants.ShatterSeek$ALT)]$SVTYPE <- 'h2hINV'

    # Retrieve the strands of each break-end.
    structuralVariants.ShatterSeek$strand1 <- as.character(BiocGenerics::strand(structuralVariants.ShatterSeek))
    structuralVariants.ShatterSeek$strand2 <- as.character(BiocGenerics::strand(structuralVariants[base::match(structuralVariants.ShatterSeek$partner, names(structuralVariants)),]))

    # Convert to ShatterSeek object.
    SV_data <- ShatterSeek::SVs(
        chrom1 = base::gsub('chr', '', base::as.character(GenomeInfoDb::seqnames(structuralVariants.ShatterSeek))),
        pos1 = as.numeric(GenomicRanges::start(structuralVariants.ShatterSeek)),

        chrom2 = base::gsub('.*(\\[|\\])', '', base::gsub(':.*', '', structuralVariants.ShatterSeek$ALT)),
        pos2 = as.integer(base::gsub('(\\[|\\]).*', '', base::gsub('.*:', '', structuralVariants.ShatterSeek$ALT))),

        SVtype = structuralVariants.ShatterSeek$SVTYPE,
        strand1 = structuralVariants.ShatterSeek$strand1,
        strand2 = structuralVariants.ShatterSeek$strand2
    )

    # Convert CN to ShatterSeek object.
    copyNumbers <- GenomeInfoDb::dropSeqlevels(copyNumbers, 'chrY', pruning.mode = 'coarse')

    CN_data <- ShatterSeek::CNVsegs(
        chrom = base::gsub('chr', '', base::as.character(GenomeInfoDb::seqnames(copyNumbers))),
        start = GenomicRanges::start(copyNumbers),
        end = GenomicRanges::end(copyNumbers),
        total_cn = base::round(base::as.numeric(copyNumbers$copyNumber), 0)
    )


    # Perform ShatterSeek -----------------------------------------------------

    output.ShatterSeek <- suppressMessages(suppressWarnings(ShatterSeek::shatterseek(SV.sample = SV_data, seg.sample = CN_data, min.Size = 3)))

    # Convert to tibble and clean-up.
    output.ShatterSeek <- tibble::as_tibble(output.ShatterSeek@chromSummary)
    output.ShatterSeek$sample <- unique(copyNumbers$sample)

    output.ShatterSeek <- output.ShatterSeek %>%
        dplyr::mutate(
            clusterSize_including_TRA = as.numeric(clusterSize_including_TRA),
            number_TRA = as.numeric(number_TRA),
            max_number_oscillating_CN_segments_3_states = as.numeric(max_number_oscillating_CN_segments_3_states),
            pval_exp_cluster = as.numeric(pval_exp_cluster),
            chr_breakpoint_enrichment = as.numeric(chr_breakpoint_enrichment)
        )


    # Determine chromothripsis regions ----------------------------------------

    output.ShatterSeek <- output.ShatterSeek %>% dplyr::mutate(
        isChromothripsis = base::ifelse(
            # Min. 25 SV (excluding translocations)
            (clusterSize_including_TRA - number_TRA) >= 25 &
                # Need CN-oscilating regions.
                (max_number_oscillating_CN_segments_2_states >= 7 | max_number_oscillating_CN_segments_3_states >= 14) &
                # Need to be a large region.
                abs((start - end)) >= 2E7 &
                # Need to be balanced with all types of SV.
                ((pval_exp_cluster <= 0.05 & !is.na(pval_exp_cluster)) | (chr_breakpoint_enrichment <= 0.05 & !is.na(chr_breakpoint_enrichment))), 'Yes', 'No')
    )


    # Return statement --------------------------------------------------------

    return(output.ShatterSeek)
}
