#' @title Performs optimal fitting of the contribution of mutational signatures.
#'
#' @details Performs mutational signature fitting against the COSMIC / Alexandrov signatures (v3.1/June 2020).
#' It does this for SNV, InDels and MNV patterns.
#'
#' @param dataMuts (VRanges): VRanges containing the mutations which will used as input.
#' @param restrictiveFit (logical): Should restrictive (strict) fitting be performed instead of using all signatures, see \link[MutationalPatterns]{fit_to_signatures_bootstrapped}?
#'
#' @examples
#' \dontrun{
#'
#'  data.Cohort <- R2CPCT::importWGSOfCohort(<cpctIds>, <combinedData>)
#'  data.MutSigs <- R2CPCT::fitMutSigs(data.Cohort$somaticVariants)
#'
#' }
#' @return (list) Returns a list of relevant mutational signature output.
#' @export
fitMutSigs <- function(dataMuts, restrictiveFit = F){

  # Input validation --------------------------------------------------------

  checkmate::assertClass(dataMuts, classes = 'VRanges')

  sprintf('Performing mutational signature fitting on %s unique samples.\nThis can take some minutes.', dplyr::n_distinct(dataMuts$sample)) %>% ParallelLogger::logInfo()


  # Retrieve the COSMIC v3.1 signatures -------------------------------------

  sprintf('\tRetrieving COSMIC (v3.1) signature matrices.') %>% ParallelLogger::logInfo()

  mutSigs.COSMIC <- list()
  mutSigs.COSMIC$SNV <- MutationalPatterns::get_known_signatures(muttype = 'snv', source = 'COSMIC', sig_type = 'reference', incl_poss_artifacts = F)
  mutSigs.COSMIC$InDel <- MutationalPatterns::get_known_signatures(muttype = 'indel', source = 'COSMIC', sig_type = 'reference', incl_poss_artifacts = F)
  mutSigs.COSMIC$DBS <- MutationalPatterns::get_known_signatures(muttype = 'dbs', source = 'COSMIC', sig_type = 'reference', incl_poss_artifacts = F)


  # Convert mutations to input matrices -------------------------------------

  sprintf('\tConverting input mutations into GRangesLists.') %>% ParallelLogger::logInfo()

  # Convert mutations to correct GRanges for input into MutationalPatterns.
  convertMuts <- function(x, DBS = F){

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
  inputMuts <- list()
  inputMuts$SNV <- convertMuts(dataMuts[dataMuts$mutType == 'SNV'])
  inputMuts$InDel <- convertMuts(dataMuts[dataMuts$mutType == 'InDel'])
  inputMuts$DBS <- convertMuts(dataMuts[dataMuts$mutType == 'MNV' & base::nchar(VariantAnnotation::ref(dataMuts)) == 2 & base::nchar(VariantAnnotation::alt(dataMuts)) == 2], DBS = T)


  # Retrieve mutational motifs ----------------------------------------------

  sprintf('\tConverting GRangesLists into mutational matrices.') %>% ParallelLogger::logInfo()

  data.mutMatrix <- list()

  # SNV (i.e. SBS)
  data.mutMatrix$SNV <- MutationalPatterns::mut_matrix(inputMuts$SNV, ref_genome =  'BSgenome.Hsapiens.UCSC.hg19')

  # InDel
  data.mutMatrix$InDel <- MutationalPatterns::get_indel_context(inputMuts$InDel, ref_genome =  'BSgenome.Hsapiens.UCSC.hg19')
  data.mutMatrix$InDel <- MutationalPatterns::count_indel_contexts(data.mutMatrix$InDel)

  # DBS
  data.mutMatrix$DBS <- MutationalPatterns::get_dbs_context(inputMuts$DBS)
  data.mutMatrix$DBS <- MutationalPatterns::count_dbs_contexts(data.mutMatrix$DBS)


  # Perform mutational signature fitting ------------------------------------

  sprintf('\tPerforming mutational signature fitting using %s.', ifelse(restrictiveFit, 'restrictive fitting', 'regular fitting method')) %>% ParallelLogger::logInfo()

  data.FittedMuts <- list()

  if(restrictiveFit){
    data.FittedMuts$SNV <- MutationalPatterns::fit_to_signatures_strict(data.mutMatrix$SNV, mutSigs.COSMIC$SNV, max_delta = 0.01)
    data.FittedMuts$InDel <- MutationalPatterns::fit_to_signatures_strict(data.mutMatrix$InDel, mutSigs.COSMIC$InDel, max_delta = 0.01)
    data.FittedMuts$DBS <- MutationalPatterns::fit_to_signatures_strict(data.mutMatrix$DBS, mutSigs.COSMIC$DBS, max_delta = 0.01)

    cleanRestrictiveFit <- function(x){
      # Remove restrictive decay figures.
      x$sim_decay_fig <- NULL

      # Pop the contribution / reconstructed matrices one level down.
      x$contribution <- x$fit_res$contribution
      x$reconstructed <- x$fit_res$reconstructed
      x$fit_res <- NULL

      # Return clean restrictive fit data.
      return(x)

    }

    data.FittedMuts$SNV <- cleanRestrictiveFit(data.FittedMuts$SNV)
    data.FittedMuts$InDel <- cleanRestrictiveFit(data.FittedMuts$InDel)
    data.FittedMuts$DBS <- cleanRestrictiveFit(data.FittedMuts$DBS)

  }else{
    data.FittedMuts$SNV <- MutationalPatterns::fit_to_signatures(data.mutMatrix$SNV, mutSigs.COSMIC$SNV)
    data.FittedMuts$InDel <- MutationalPatterns::fit_to_signatures(data.mutMatrix$InDel, mutSigs.COSMIC$InDel)
    data.FittedMuts$DBS <- MutationalPatterns::fit_to_signatures(data.mutMatrix$DBS, mutSigs.COSMIC$DBS)
  }


  # Clean-up and add proposed aetologies ------------------------------------

  sprintf('\tCleaning up and calculating the rel. contributions.') %>% ParallelLogger::logInfo()

  cleanSigs <- function(x){

    # Calculate rel. contribution.
    relContribution <- tibble::as_tibble(base::sweep(x, 2, base::colSums(x), "/") * 100, rownames = "Signature") %>%
      # Melt.
      tidyr::pivot_longer(cols = !dplyr::contains('Signature'), names_to = 'sampleId', values_to = 'relContribution') %>%
      # Add proposed aetiology.
      dplyr::inner_join(R2CPCT::proposedAetiologyCOSMICv3.1, by = 'Signature')

    return(relContribution)

  }

  data.FittedMuts$SNV$relativeContribution <- cleanSigs(data.FittedMuts$SNV$contribution)
  data.FittedMuts$InDel$relativeContribution <- cleanSigs(data.FittedMuts$InDel$contribution)
  data.FittedMuts$DBS$relativeContribution <- cleanSigs(data.FittedMuts$DBS$contribution)

  # Add the mutational matrices.
  data.FittedMuts$SNV$mutMatrix <- data.mutMatrix$SNV
  data.FittedMuts$InDel$mutMatrix <- data.mutMatrix$InDel
  data.FittedMuts$DBS$mutMatrix <- data.mutMatrix$DBS


  # Return statement --------------------------------------------------------

  return(data.FittedMuts)

}