#' @title Determine the number of Ti and Tv contexts and Ti/Tv ratio per sample.
#'
#' @param mutMatrix.SNV (matrix): Mutational motif matrix derived from the \link[R2CPCT]{fitMutSigs} function.
#'
#' @examples
#' \donttest{
#'
#'  mutSigs <- R2CPCT::fitMutSigs(data.Cohort$somaticVariants, restrictiveFit = FALSE)
#'  determineTiTv(mutSigs$SNV$mutMatrix)
#'
#' }
#' @return (tibble) Returns a tibble of the Ti/Tv number and ratio per sample.
#' @export
determineTiTv <- function(mutMatrix.SNV){
    
    # Input validation --------------------------------------------------------
    
    checkmate::assertMatrix(mutMatrix.SNV)
    
    sprintf('Calculating Ti/Tv ratios for %s sample(s).', dplyr::n_distinct(base::colnames(mutMatrix.SNV))) %>% ParallelLogger::logInfo()
    
    
    # Determine Ti/Tv mutations and ratio -------------------------------------
    
    data.MutContexts <- mutMatrix.SNV %>%
        reshape2::melt() %>%
        dplyr::select(mutContext = Var1, sample = Var2, value) %>%
        dplyr::mutate(
            mutationalContext = base::gsub('\\[.*]', '', mutContext),
            mutationalType = base::gsub('].*', '', base::gsub('.*\\[', '', mutContext)),
            # Determine CpG mutations.
            mutationalType = base::ifelse(base::grepl('C>T', mutationalType) & base::grepl('G$', mutationalContext), base::paste(mutationalType, '(CpG)'), mutationalType),
            # Determine Ti / Tv mutations.
            mutationalType = base::ifelse(mutationalType %in% c('C>T', 'C>T (CpG)', 'T>C', 'G>A', 'A>T'), base::paste(mutationalType, 'Ti', sep = '\n'), base::paste(mutationalType, 'Tv', sep = '\n'))
        ) %>%
        # Count total mut. contexts per sample.
        dplyr::group_by(sample, mutationalType) %>%
        dplyr::summarise(totalValue = base::sum(value)) %>%
        dplyr::ungroup() %>%
        
        # Determine Ti/Tv ratio per sample.
        dplyr::mutate(TiTvType = ifelse(grepl('\nTv', mutationalType), 'Transversion', 'Transition')) %>%
        dplyr::group_by(sample, TiTvType) %>% dplyr::mutate(totalGroupInSample = sum(totalValue)) %>% dplyr::ungroup() %>%
        dplyr::group_by(sample) %>% dplyr::mutate(TiTvRatio = unique(totalGroupInSample[TiTvType == 'Transition']) / unique(totalGroupInSample[TiTvType == 'Transversion'])) %>% dplyr::ungroup()
    
    
    # Return statement --------------------------------------------------------
    
    return(data.MutContexts)
}
