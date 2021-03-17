#' @title Performs CHORD analysis.
#'
#' @param somaticVariants (VRanges): VRanges of a single sample from \link[R2CPCT]{importSomaticVariantsVEP}.
#' @param structuralVariants (GRanges): GRanges of a single sample from \link[R2CPCT]{importStructuralVariantsPURPLE}.
#' @return (list) list with the CHORD contexts and results.
#' @examples
#' \dontrun{
#'
#' 	runCHORD(dataSample$somaticVariants, dataSample$structuralVariants)
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
performCHORD <- function(somaticVariants, structuralVariants){

    # Input validation --------------------------------------------------------

    checkmate::assertClass(somaticVariants, 'VRanges')
    checkmate::assertClass(structuralVariants, 'GRanges')

    sprintf('Running CHORD (HRD prediction) for: %s', unique(somaticVariants$sample)) %>% ParallelLogger::logInfo()


    # Run CHORD ---------------------------------------------------------------

    df.SNV <- somaticVariants[somaticVariants$mutType == 'SNV']
    df.SNV <- data.frame(chrom = GenomeInfoDb::seqnames(df.SNV), pos = GenomicRanges::start(df.SNV), ref = VariantAnnotation::ref(df.SNV), alt = VariantAnnotation::alt(df.SNV)) %>%
        dplyr::mutate(
            chrom = as.character(chrom),
            ref = as.character(ref),
            alt = as.character(alt),
        )

    df.InDel <- somaticVariants[somaticVariants$mutType == 'InDel']
    df.InDel <- data.frame(chrom = GenomeInfoDb::seqnames(df.InDel), pos = GenomicRanges::start(df.InDel), ref = VariantAnnotation::ref(df.InDel), alt = VariantAnnotation::alt(df.InDel)) %>%
        dplyr::mutate(
            chrom = as.character(chrom),
            ref = as.character(ref),
            alt = as.character(alt),
        )

    df.SV <- structuralVariants[!duplicated(structuralVariants$EVENT) & !structuralVariants$SVTYPE %in% c('INS', 'SINGLE'),]
    df.SV <- data.frame(sv_type = df.SV$SVTYPE, sv_len = abs(df.SV$svLen))

    # Retrieve CHORD contexts.
    contexts <- CHORD::extractSigsChord(
        df.snv = df.SNV,
        df.indel = df.InDel,
        df.sv = df.SV,
        sample.name = unique(somaticVariants$sample),
        ref.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    )

    # Retrieve CHORD predictions.
    outputCHORD <- tibble::as_tibble(CHORD::chordPredict(contexts, show.features = T, rf.model = CHORD::CHORD, verbose = T))


    # Return statement --------------------------------------------------------

    return(list(contexts = contexts, results = outputCHORD))
}