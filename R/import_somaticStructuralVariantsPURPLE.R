#' @title Import the somatic structural variants as estimated by GRIDSS / PURPLE.
#'
#' @param pathSV (character): Path to the <sample>.purple.sv.vcf.gz file containing the somatic structural variants as detected by GRIDSS / PURPLE.
#' @param passOnly (logical): Only retain SV which were not found in the PON and/or were inferred by copy-number transition  alone or any of the other filters.
#' @return (GRanges) GRanges object containing the break-ends of the structural variants, so two rows for each SV (except when single break-end).
#' @examples
#' \donttest{
#'
#' 	importStructuralVariantsPURPLE(pathSV = '<sample>.purple.sv.vcf.gz')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importStructuralVariantsPURPLE <- function(pathSV, passOnly = TRUE){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathSV, access = 'r')
    checkmate::assertLogical(passOnly)


    # Read VCF ----------------------------------------------------------------

    sprintf('Importing SV: %s', pathSV) %>% ParallelLogger::logInfo()

    # Clean seqlevels and add chromosome information.
    sample.SV.VCF <- VariantAnnotation::readVcf(file = pathSV, genome = 'hg19') %>%
        R2CPCT::cleanSeqlevels(excludeChr = NULL)

    # Convert to GRanges; retain both single and paired break-ends.
    sample.SV <- c(
        base::suppressWarnings(StructuralVariantAnnotation::breakpointRanges(sample.SV.VCF, unpartneredBreakends = FALSE)),
        base::suppressWarnings(StructuralVariantAnnotation::breakpointRanges(sample.SV.VCF, unpartneredBreakends = TRUE))
    )

    # Only count each event once.
    sample.SV <- sample.SV[!duplicated(sample.SV$sourceId),]

    # Retrieve the information of the SV and remove duplicate fields.
    infoSV <- VariantAnnotation::info(sample.SV.VCF)
    infoSV <- infoSV[!base::toupper(base::colnames(infoSV)) %in% base::toupper(base::colnames(S4Vectors::mcols(sample.SV)))]
    infoSV$sourceId <- base::rownames(sample.SV.VCF)
    infoSV <- infoSV[!duplicated(infoSV$sourceId),]

    # Sort in the same order.
    infoSV <- infoSV[base::order(infoSV$sourceId),]
    sample.SV <- sample.SV[base::order(sample.SV$sourceId),]

    # Merge information.
    S4Vectors::mcols(sample.SV) <- S4Vectors::merge(S4Vectors::mcols(sample.SV), infoSV, by = 'sourceId', all.x = TRUE)

    # Should only SV passing all filters (GRIDSS) be imported?
    if(passOnly){
        totalPriorPass <- base::length(sample.SV)
        sample.SV <- sample.SV[sample.SV$FILTER == 'PASS']
        sprintf('\tRetaining %s / %s PASS-only somatic structural variants.', base::length(sample.SV), totalPriorPass) %>% ParallelLogger::logInfo()
    }

    # If no SV are left post-filtering, return NULL.
    if(base::length(sample.SV) == 0) return(NULL)


    # Determine SV-type -------------------------------------------------------

    # Taken from GRIDSS helper script.
    simpleEventType <- function(gr) {

        gr$svtype <- NULL

        # Remove single partners (INS)
        gr$SVTYPE <- 'Other'
        gr.paired <- gr[!is.na(gr$partner),]

        # For each (paired) SV, retrieve the partner record.
        gr.paired.partner <- StructuralVariantAnnotation::partner(gr[!is.na(gr$partner),])

        # Translocations - Discordant chromosomes.
        gr.paired$SVTYPE <- ifelse(gr.paired$SVTYPE == 'Other' & GenomeInfoDb::seqnames(gr.paired) != GenomeInfoDb::seqnames(gr.paired.partner), 'TRA', gr.paired$SVTYPE)

        # Inversions - Head-tail have the same strands on ref. genome.
        gr.paired$SVTYPE <- ifelse(gr.paired$SVTYPE == 'Other' & BiocGenerics::strand(gr.paired) == BiocGenerics::strand(gr.paired.partner), 'INV', gr.paired$SVTYPE)

        # Insertions - Length of inserted nucleotides should be larger than some magic number.
        gr.paired$SVTYPE <- ifelse(gr.paired$SVTYPE == 'Other' & (gr.paired$insLen >= abs(gr.paired$svLen) * .7) & !is.na(gr.paired$svLen), 'INS', gr.paired$SVTYPE)

        # Deletions - Start of upstream partner in front of start from downstream partner.
        gr.paired$SVTYPE <- ifelse(gr.paired$SVTYPE == 'Other' & xor(BiocGenerics::start(gr.paired) < BiocGenerics::start(gr.paired.partner), BiocGenerics::strand(gr.paired) == "-"), "DEL", gr.paired$SVTYPE)

        # Tandem Duplication - If not any of the other classes.
        gr.paired$SVTYPE <- ifelse(gr.paired$SVTYPE == 'Other', 'DUP', gr.paired$SVTYPE)

        # Re-add the single breakpoints and set as SINGLE.
        gr.svtype <- base::c(gr[is.na(gr$partner),], gr.paired)
        gr.svtype$SVTYPE <- ifelse(gr.svtype$SVTYPE == 'Other', 'SINGLE', gr.svtype$SVTYPE)

        # TODO - Check if SINGLE break-ends aren't insertions as they have a large insLen.

        # Return SV with their type.
        return(gr.svtype)

    }

    # Determine SV-type.
    sample.SV <- simpleEventType(sample.SV)

    # Sort on chromosomal order.
    sample.SV <- GenomicRanges::sort(sample.SV)

    # Add sample identifier.
    sample.SV$sample <- base::factor(base::gsub('\\.purple.sv.*', '', base::basename(pathSV)))


    # Return statement --------------------------------------------------------

    sprintf('\tReturning GRanges') %>% ParallelLogger::logTrace()

    return(sample.SV)

}
