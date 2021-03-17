#' @title Cleans and standardizes chromosomal information.
#'
#' @description Removes unused chromosomes, converts seqlevels to UCSC nomenclature and adds chromosomal information.
#'
#' @param g (GRanges): GRanges object to be cleaned.
#' @param excludeChr (character): Character vector of chromosomes which should be excluded, e.g. 'chrY' if patient is female or chrM. Set NULL to ignore.
#'
#' @examples
#' \dontrun{
#'  gr0 <- GenomicRanges::GRanges(
#'  S4Vectors::Rle(c('chr2', 'chr2', 'chr1', 'chr3'), c(1, 3, 2, 4)), 
#'  IRanges::IRanges(1:10, width=10:1)
#'  )
#' 	gr0.clean <- cleanSeqlevels(gr0)
#'
#' 	GenomeInfoDb::seqinfo(gr0.clean)
#'
#' }
#' @return (GRanges) Returns cleaned-up GRanges.
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family Misc.
#' @export
cleanSeqlevels <- function(g, excludeChr = NULL){

    # Input validation --------------------------------------------------------

    checkmate::assert(
        checkmate::checkClass(g, 'GRanges'),
        checkmate::checkClass(g, 'CollapsedVCF')
    )
    checkmate::assertCharacter(excludeChr, null.ok = T)


    # Start -------------------------------------------------------------------

    sprintf('Cleaning GRanges object using GRCh37 genome') %>% ParallelLogger::logTrace()

    # Remove haplochromosomes and use UCSC nomenclature.
    g <- GenomeInfoDb::keepStandardChromosomes(g, species = 'Homo_sapiens', pruning.mode = 'coarse')
    g <- GenomeInfoDb::renameSeqlevels(g, GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(g), 'UCSC'))
    g <- GenomeInfoDb::sortSeqlevels(g, X.is.sexchrom = T)

    # Use pre-packages chromosomal information or fetch from UCSC.
    futile.logger::flog.trace('Retrieving chromosome lengths')
    suppressMessages(require(BSgenome.Hsapiens.UCSC.hg19)); chrInfo <- GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)

    # Keep only the used chromosomes.
    chrInfo <- suppressWarnings(GenomeInfoDb::keepSeqlevels(chrInfo, GenomeInfoDb::seqlevelsInUse(g), pruning.mode = 'coarse'))
    g <- suppressWarnings(GenomeInfoDb::keepSeqlevels(g, GenomeInfoDb::seqlevels(chrInfo), pruning.mode = 'coarse'))

    # Add to input GRanges object.
    suppressWarnings(GenomeInfoDb::seqinfo(g) <- chrInfo)

    # Trim elements passing the chromosome ends.
    if(!class(g) == 'CollapsedVCF') g <- IRanges::trim(g)
    if(!class(g) == 'CollapsedVCF') g <- g[GenomicRanges::start(g) <= GenomicRanges::end(g)]

    # Remove chromosomes.
    if(!is.null(excludeChr)) g <- GenomeInfoDb::dropSeqlevels(g, excludeChr, pruning.mode = 'coarse')

    return(g)
}