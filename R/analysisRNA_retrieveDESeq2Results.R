#' @title Retrieve the results from a DESeq2-analysis using IHW and LFC shrinkage.
#'
#' @param DESeq2.Object (DESeqDataSet): DESeqDataSet of a DESeq2 analysis.
#' @param contrast (character): Contrast of the resulst, e.g.: c('testVariable', 'groupA', 'groupb')
#' @param nThreads (integer): Number of threads using during multi-threading
#'
#' @return (tibble) Tibble containing the DESeq2 results.
#' @examples
#' \donttest{
#'
#' DESeq2.allSamples <- DESeq2::DESeqDataSetFromMatrix(countData = countMatrix, colData = samples.meta[match(colnames(countMatrix), samples.meta$sample),], design = ~treatmentType)
#' DESeq2.allSamples <- DESeq2::DESeq(DESeq2.allSamples, test = 'Wald', parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 20))
#' retrieveDESeq2Results(DESeq2.allSamples)
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family RNA-Seq
#' @export
retrieveDESeq2Results <- function(DESeq2.Object, contrast, nThreads = 20){
    
    
    # Input validation --------------------------------------------------------
    
    checkmate::checkClass(DESeq2.Object, 'DESeqDataSet')
    checkmate::checkCharacter(contrast)
    checkmate::assertInt(nThreads)
    
    data('GENCODE.v38', package = 'R2CPCT')
    base::sprintf('Retrieving DESeq2 results using IHW-filtering and LFC-shrinkage (ashr).') %>% ParallelLogger::logInfo()
    
    
    # Get results -------------------------------------------------------------
    
    # Get diff. results.
    results <- DESeq2::results(DESeq2.Object, pAdjustMethod = 'BH', filterFun = IHW::ihw, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = nThreads), tidy = FALSE, contrast = contrast)
    
    # Shrink the LFC.
    results.LFC <- DESeq2::lfcShrink(DESeq2.Object, res = results, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = nThreads), contrast = contrast, type = 'ashr')
    
    # Add ENSEMBL as column.
    results.LFC$ENSEMBL <- BiocGenerics::rownames(results.LFC)
    
    # Re-add the t-statistic
    results.LFC$stat <- results[match(rownames(results), rownames(results.LFC)),]$stat
    
    # Re-add the q weight.
    results.LFC$weight <- results[match(rownames(results), rownames(results.LFC)),]$weight
    
    # Convert to tibble.
    results.LFC <- tibble::as_tibble(results.LFC)
    
    # Add contrast.
    results.LFC$contrast <- paste(contrast, collapse = '_')
    
    # Add the SYMBOL.
    results.LFC <- results.LFC %>% dplyr::left_join(tibble::as_tibble(S4Vectors::mcols(GENCODE.v38)) %>% dplyr::distinct(SYMBOL, ENSEMBL), by = 'ENSEMBL')
    
    
    # Return statement --------------------------------------------------------
    
    return(results.LFC)
    
}
