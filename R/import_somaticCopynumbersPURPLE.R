#' @title Import the somatic CNV segments as estimated by PURPLE.
#'
#' @details Will also convert the absolute copynumbers into CN-States for downstream analysis.
#'
#' @param pathCNV (character): Path to the <sample>.purple.cnv.somatic.tsv file containing the somatic CNV regions as detected by PURPLE.
#'
#' @return (GRanges) GRanges object containing the somatic variants and assorted annotations.
#' @examples
#' \dontrun{
#'
#' 	importSomaticCopynumberPURPLE(pathCNV = '<sample>.purple.cnv.somatic.tsv')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importSomaticCopynumberPURPLE <- function(pathCNV){

    # Input validation --------------------------------------------------------

    checkmate::assertAccess(pathCNV, access = 'r')


    # Read CNV ----------------------------------------------------------------

    sprintf('Importing CNV: %s', pathCNV) %>% ParallelLogger::logInfo()

    # Clean seqlevels and add chromosome information.
    sample.CNV <- readr::read_tsv(pathCNV, col_types = 'cddddddcccdddddd') %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
        R2CPCT::cleanSeqlevels()

    # Add sample name
    sample.CNV$sample <- base::factor(base::gsub('\\.purple.*', '', base::basename(pathCNV)))


    # Return statement --------------------------------------------------------

    sprintf('\tReturning GRanges') %>% ParallelLogger::logTrace()

    return(sample.CNV)

}
