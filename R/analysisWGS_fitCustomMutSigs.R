#' @title Performs mutational signature fitting against a pre-defined list of (custom) motif matrixes.
#'
#' @details Users need to supply a motif matrix (SNV-only), the names of the columns (signatures) will be used a aetiology.
#'
#' @param dataMuts (VRanges): VRanges containing the mutations which will used as input.
#' @param motifMatrix (tibble): Mutational motif matrix,  with a column termed 'Motif' and then x columns designating the signatures.
#' @param sigType (character): What type of signature should be looked for (motif matrix should match).
#'
#' @examples
#' \donttest{
#'
#'  data.Cohort <- R2CPCT::importWGSOfCohort(<cpctIds>, <combinedData>)
#'  data.MutSigs <- R2CPCT::fitMutSigs(data.Cohort$somaticVariants, motifMatrix)
#'
#' }
#' @return (list) Returns a list of relevant mutational signature output.
#' @export
fitCustomMutSigs <- function(dataMuts, motifMatrix, sigType = 'SBS'){

    # Input validation --------------------------------------------------------

    checkmate::assertClass(dataMuts, classes = 'VRanges')
    checkmate::assertClass(motifMatrix, classes = 'tbl_df')
    checkmate::assertCharacter(sigType, pattern = 'SBS|DBS|ID', len = 1)
    checkmate::assertTRUE('Motif' %in% colnames(motifMatrix))

    sprintf('Performing (custom) mutational signature fitting (%s) on %s unique samples.\nThis can take some minutes.', sigType, dplyr::n_distinct(dataMuts$sample)) %>% ParallelLogger::logInfo()


    # Transforming motif matrix -----------------------------------------------

    motifMatrix.m <- as.matrix(motifMatrix[colnames(motifMatrix) != 'Motif'])


    # Convert mutations to input matrices -------------------------------------

    sprintf('\tConverting input mutations into GRangesLists.') %>% ParallelLogger::logInfo()

    # Convert mutations to correct GRanges for input into MutationalPatterns.
    convertMuts <- function(x, DBS = FALSE){

        S4Vectors::mcols(x) <- S4Vectors::DataFrame(sample = x$sample)

        # Add REF and ALT as column.
        if(DBS){
            x$REF <- Biostrings::DNAStringSet(VariantAnnotation::ref(x))
            x$ALT <- Biostrings::DNAStringSetList(base::lapply(VariantAnnotation::alt(x), Biostrings::DNAStringSet))
        }else{
            x$REF <- VariantAnnotation::ref(x)
            x$ALT <- VariantAnnotation::alt(x)
        }

        # Convert to GRangeslist, split per sample.
        x <- GenomicRanges::GRanges(x)
        x <- GenomicRanges::GRangesList(base::split(x, x$sample))

        # Return.
        return(x)
    }

    # Generate a GRangesList, split per sample, per mutational type.
    if(sigType == 'SBS') inputMuts <- convertMuts(dataMuts[dataMuts$mutType == 'SNV'])
    if(sigType == 'DBS') inputMuts <- convertMuts(dataMuts[dataMuts$mutType == 'MNV' & base::nchar(VariantAnnotation::ref(dataMuts)) == 2 & base::nchar(VariantAnnotation::alt(dataMuts)) == 2], DBS = TRUE)
    if(sigType == 'ID') inputMuts <- convertMuts(dataMuts[dataMuts$mutType == 'InDel'])


    # Retrieve mutational motifs ----------------------------------------------

    sprintf('\tConverting GRangesLists into mutational matrices.') %>% ParallelLogger::logInfo()

    if(sigType == 'SBS'){
        data.mutMatrix <- MutationalPatterns::mut_matrix(inputMuts, ref_genome =  'BSgenome.Hsapiens.UCSC.hg19')
    }

    if(sigType == 'DBS'){
        data.mutMatrix <- MutationalPatterns::get_dbs_context(inputMuts)
        data.mutMatrix <- MutationalPatterns::count_dbs_contexts(data.mutMatrix)
    }

    if(sigType == 'ID'){
        data.mutMatrix <- MutationalPatterns::get_indel_context(inputMuts, ref_genome =  'BSgenome.Hsapiens.UCSC.hg19')
        data.mutMatrix <- MutationalPatterns::count_indel_contexts(data.mutMatrix)
    }

    # Sort on input motifs.
    data.mutMatrix <- data.mutMatrix[base::match(rownames(data.mutMatrix), motifMatrix$Motif),]


    # Perform mutational signature fitting ------------------------------------

    sprintf('\tPerforming mutational signature fitting for %s signatures.', base::ncol(motifMatrix.m)) %>% ParallelLogger::logInfo()

    data.FittedMuts <- MutationalPatterns::fit_to_signatures(data.mutMatrix, motifMatrix.m)


    # Clean-up and add proposed aetologies ------------------------------------

    sprintf('\tCleaning up and calculating the rel. contributions.') %>% ParallelLogger::logInfo()

    cleanSigs <- function(x){

        # Calculate rel. contribution.
        relContribution <- tibble::as_tibble(base::sweep(x, 2, base::colSums(x), "/") * 100, rownames = "Signature") %>%
            # Melt.
            tidyr::pivot_longer(cols = !dplyr::contains('Signature'), names_to = 'sampleId', values_to = 'relContribution')

        return(relContribution)

    }

    data.FittedMuts$relativeContribution <- cleanSigs(data.FittedMuts$contribution)

    # Add the mutational matrices.
    data.FittedMuts$mutMatrix <- data.mutMatrix


    # Return statement --------------------------------------------------------

    return(data.FittedMuts)

}
