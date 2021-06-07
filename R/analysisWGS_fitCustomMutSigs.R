#' @title Performs mutational signature fitting against a pre-defined list of (custom) motif matrixes.
#'
#' @details Users need to supply a motif matrix (SNV-only), the names of the columns (signatures) will be used a aetiology.
#'
#' @param dataMuts (VRanges): VRanges containing the mutations which will used as input.
#' @param motifMatrix (tibble): Mutational motif matrix,  with a column termed 'Motif' and then x columns designating the signatures.
#'
#' @examples
#' \dontrun{
#'
#'  data.Cohort <- R2CPCT::importWGSOfCohort(<cpctIds>, <combinedData>)
#'  data.MutSigs <- R2CPCT::fitMutSigs(data.Cohort$somaticVariants, motifMatrix)
#'
#' }
#' @return (list) Returns a list of relevant mutational signature output.
#' @export
fitCustomMutSigs <- function(dataMuts, motifMatrix){

    # Input validation --------------------------------------------------------

    checkmate::assertClass(dataMuts, classes = 'VRanges')
    checkmate::assertClass(motifMatrix, classes = 'tbl_df')
    checkmate::assertTRUE('Motif' %in% colnames(motifMatrix))

    sprintf('Performing (custom) mutational signature fitting on %s unique samples.\nThis can take some minutes.', dplyr::n_distinct(dataMuts$sample)) %>% ParallelLogger::logInfo()


    # Transforming motif matrix -----------------------------------------------

    motifMatrix.m <- as.matrix(motifMatrix[colnames(motifMatrix) != 'Motif'])


    # Convert mutations to input matrices -------------------------------------

    sprintf('\tConverting input mutations into GRangesLists.') %>% ParallelLogger::logInfo()

    # Convert mutations to correct GRanges for input into MutationalPatterns.
    convertMuts <- function(x){

        S4Vectors::mcols(x) <- S4Vectors::DataFrame(sample = x$sample)

        # Add REF and ALT as column.
        x$REF <- VariantAnnotation::ref(x)
        x$ALT <- VariantAnnotation::alt(x)

        # Convert to GRangeslist, split per sample.
        x <- GenomicRanges::GRanges(x)
        x <- GenomicRanges::GRangesList(base::split(x, x$sample))

        # Return.
        return(x)
    }

    # Generate a GRangesList, split per sample, per mutational type.
    inputMuts <- list()
    inputMuts$SNV <- convertMuts(dataMuts[dataMuts$mutType == 'SNV'])


    # Retrieve mutational motifs ----------------------------------------------

    sprintf('\tConverting GRangesLists into mutational matrices.') %>% ParallelLogger::logInfo()

    data.mutMatrix <- list()

    # SNV (i.e. SBS)
    data.mutMatrix$SNV <- MutationalPatterns::mut_matrix(inputMuts$SNV, ref_genome =  'BSgenome.Hsapiens.UCSC.hg19')

    # Sort on input motifs.
    data.mutMatrix$SNV <- data.mutMatrix$SNV[base::match(rownames(data.mutMatrix$SNV), motifMatrix$Motif),]


    # Perform mutational signature fitting ------------------------------------

    sprintf('\tPerforming mutational signature fitting for %s signatures.', base::ncol(motifMatrix.m)) %>% ParallelLogger::logInfo()

    data.FittedMuts <- list()

    data.FittedMuts$SNV <- MutationalPatterns::fit_to_signatures(data.mutMatrix$SNV, motifMatrix.m)


    # Clean-up and add proposed aetologies ------------------------------------

    sprintf('\tCleaning up and calculating the rel. contributions.') %>% ParallelLogger::logInfo()

    cleanSigs <- function(x){

        # Calculate rel. contribution.
        relContribution <- tibble::as_tibble(base::sweep(x, 2, base::colSums(x), "/") * 100, rownames = "Signature") %>%
            # Melt.
            tidyr::pivot_longer(cols = !dplyr::contains('Signature'), names_to = 'sampleId', values_to = 'relContribution')

        return(relContribution)

    }

    data.FittedMuts$SNV$relativeContribution <- cleanSigs(data.FittedMuts$SNV$contribution)

    # Add the mutational matrices.
    data.FittedMuts$SNV$mutMatrix <- data.mutMatrix$SNV


    # Return statement --------------------------------------------------------

    return(data.FittedMuts)

}