#' @title Generate a gene-level overview of all overlapping mutations (incl. SV and CNA) and dN/dS + GISTIC2 results.
#'
#' @description The gene-level overview will be performed on a per-sample basis using only protein-coding genes.
#' The combined report will only contain protein-coding genes which have at least a single somatic aberration or were detected in the dN/dS or GISTIC2 analysis.
#'
#' Use the 'selectedGenes' parameter to retrieve your genes-of-interests, even if no aberrations were detected.
#'
#' @param data.Cohort (list): Cohort-wide data obtained from \link[R2CPCT]{importWGSOfCohort}.
#' @param dNdS (list): Output of the dN/dS analysis obtained from \link[R2CPCT]{rundNdS}.
#' @param GISTIC2 (list): Output of the GISTIC2 analysis obtained from \link[R2CPCT]{importGISTIC2}.
#' @param nThreads (integer): Number of cores over which to parallelize (when needed).
#' @param mutantsOnly (logical): Only output records with mutations (CNA, LOH or Muts; TRUE) or return all records (FALSE).
#'
#' @examples
#' \donttest{
#'
#' 	# Generate a combined gene-level report of an imported WGS cohort (which could be a subset).
#' 	generateCombinedReport(data.Cohort, results.Cohort$dNdS, results.Cohort$GISTIC2)
#'
#' }
#'
#' @return (tibble) Returns a per-sample gene-level overview of all reported somatic aberrations on protein-coding genes.
#'
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @export
generateCombinedReport <- function(data.Cohort, dNdS, GISTIC2, nThreads = 40, mutantsOnly = TRUE){
  
  # Input validation --------------------------------------------------------
  
  checkmate::checkList(data.Cohort)
  checkmate::checkList(dNdS)
  checkmate::checkList(GISTIC2)
  checkmate::checkNumber(nThreads)
  checkmate::checkLogical(mutantsOnly)
  
  data('GENCODE.v38', package = 'R2CPCT')
  data('GENCODE.v38.Exons', package = 'R2CPCT')
  data('driverList', package = 'R2CPCT')
  
  base::sprintf('Generating per-sample gene-level overview of various aberrations for %s samples.', dplyr::n_distinct(data.Cohort$somaticVariants$sample)) %>% ParallelLogger::logInfo()
  
  
  # Gather required data. ---------------------------------------------------
  
  somaticData <- list()
  
  # Add gender and ploidy information to determine amplification/deletion status.
  somaticData$somaticVariants <- tibble::as_tibble(S4Vectors::mcols(data.Cohort$somaticVariants))
  somaticData$copyNumbers <- data.Cohort$copyNumbers
  somaticData$dNdS <- dNdS$finalOutput
  somaticData$driverCatalogHMF <- data.Cohort$driverCatalog
  somaticData$driverLINX <- data.Cohort$driverLINX
  somaticData$GISTIC2 <- GISTIC2$gisticNarrowPeaksWithAnno
  somaticData$SV <- data.Cohort$structuralVariants
  purityStats <- data.Cohort$purityStats
  
  
  # Determine overlap of CN with genes --------------------------------------
  
  base::sprintf('Per sample, generating a gene-level overview of purity-corrected exon-overlapping copynumbers (in parallel; can still take a few minutes).') %>% ParallelLogger::logInfo()
  
  combinedReport.Copynumbers <- pbapply::pblapply(base::split(somaticData$copyNumbers, somaticData$copyNumbers$sample), function(copyNumbers.Sample){
    
    # Retrieve overlapping genes (using their exons) per CN-segment (min 10bp overlap).
    overlapCN <- IRanges::findOverlaps(query = GENCODE.v38.Exons, subject = copyNumbers.Sample, minoverlap = 10, select = 'all')
    
    GENCODE.Overlap <- tibble::as_tibble(S4Vectors::mcols(GENCODE.v38.Exons))
    GENCODE.Overlap$chr <- as.character(GenomeInfoDb::seqnames(GENCODE.v38.Exons))
    
    # Per overlapping CN-segment, retrieve the exon(s).
    GENCODE.Overlap <- GENCODE.Overlap[S4Vectors::queryHits(overlapCN),]
    
    # Determine width of overlap of gene with CN-segment.
    overlaps <- IRanges::pintersect(GENCODE.v38.Exons[S4Vectors::queryHits(overlapCN),], copyNumbers.Sample[S4Vectors::subjectHits(overlapCN)])
    GENCODE.Overlap$overlappingWidth <- IRanges::width(overlaps)
    
    # Add the CN-segment copy-number and purity/ploidy-corrected BAF.
    GENCODE.Overlap$copyNumber <- copyNumbers.Sample[S4Vectors::subjectHits(overlapCN)]$copyNumber
    GENCODE.Overlap$BAF <- copyNumbers.Sample[S4Vectors::subjectHits(overlapCN)]$baf
    GENCODE.Overlap$nHetMarkers <- copyNumbers.Sample[S4Vectors::subjectHits(overlapCN)]$bafCount
    
    # Determine mean copy-number and BAF per gene.
    # If overlapping with multiple segments, take the mean of every base-level CN of the overlapping segments.
    # Determine the BAF based on overlapping segments with at least 25 heterozygous markers.
    GENCODE.Overlap.perGene <- GENCODE.Overlap %>%
      dplyr::group_by(SYMBOL, ENSEMBL, chr) %>%
      dplyr::summarise(
        
        # Determine CN.
        copyNumber.Mean = base::ifelse(dplyr::n() == 1, copyNumber, base::mean(base::rep(copyNumber, overlappingWidth))),
        
        # Determine BAF
        BAF.Mean = base::ifelse(any(nHetMarkers >= 25), base::round(base::ifelse(dplyr::n() == 1, BAF, base::mean(base::rep(BAF[nHetMarkers >= 25], nHetMarkers[nHetMarkers >= 25]))), 2), NA),
        BAF.Mean = base::ifelse(BAF.Mean > 1, 1, BAF.Mean),
        BAF.Mean = base::ifelse(BAF.Mean < 0, 0, BAF.Mean),
        
        # Add CN of each exon.
        copyNumbers.Exons = base::list(copyNumber),
        totalExonsInGene = base::length(copyNumber)
        
      ) %>%
      dplyr::ungroup()
    
    # Remove genes present on more than one chromosome (e.g. PAR)
    GENCODE.Overlap.perGene <- GENCODE.Overlap.perGene %>% dplyr::filter(!base::duplicated(base::paste0(ENSEMBL, SYMBOL)))
    
    # Re-add sample identifier.
    GENCODE.Overlap.perGene$sample <- unique(copyNumbers.Sample$sample)
    
    # Add gender and genome-wide ploidy.
    GENCODE.Overlap.perGene <- GENCODE.Overlap.perGene %>% dplyr::left_join(purityStats %>% dplyr::distinct(sample, gender, ploidy), by = c('sample' = 'sample'))
    
    # Determine copy-number status of the overall gene-body (using all overlapping exons).
    GENCODE.Overlap.perGene <- GENCODE.Overlap.perGene %>%
      dplyr::mutate(
        Consequence.CNA = base::ifelse(
          gender == 'MALE' & chr %in% c('chrX', 'chrY'),
          base::as.character(base::cut(copyNumber.Mean, c(-Inf, base::max(c(ifelse((unique(ploidy) - 1) / 2 <= 0.75, ((unique(ploidy) - 1) / 2) / 2, 0.75), (unique(ploidy) - 1) / 3)), (unique(ploidy) - 1) / 1.5, (1.5 * (unique(ploidy) - 1)), min(15, (3 * (unique(ploidy) - 1))), Inf), labels = c('Deep Deletion', 'Deletion', 'Neutral', 'Amplification' ,'Deep Amplification'))),
          base::as.character(base::cut(copyNumber.Mean, c(-Inf, base::max(c(ifelse((unique(ploidy)) / 2 <= 0.75, (unique(ploidy) / 2) / 2, 0.75), (unique(ploidy)) / 3)), (unique(ploidy)) / 1.5, (1.5 * (unique(ploidy))), min(15, (3 * (unique(ploidy)))), Inf), labels = c('Deep Deletion', 'Deletion', 'Neutral', 'Amplification' ,'Deep Amplification')))
        )
      ) %>%
      # Determine LOH status.
      dplyr::mutate(
        Consequence.LOH = base::ifelse(BAF.Mean >= 0.85, 'LOH', NA)
      )
    
    # Determine number of deleted or amplified exons (per gene).
    GENCODE.Overlap.perGene <- GENCODE.Overlap.perGene %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        Consequence.CNA.Exons = base::ifelse(
          gender == 'MALE' & chr %in% c('chrX', 'chrY'),
          base::list(base::as.character(base::cut(base::unlist(copyNumbers.Exons), c(-Inf, base::max(c(ifelse((unique(ploidy) - 1) / 2 <= 0.75, ((unique(ploidy) - 1) / 2) / 2, 0.75), (unique(ploidy) - 1) / 3)), (unique(ploidy) - 1) / 1.5, (1.5 * (unique(ploidy) - 1)), min(15, (3 * (unique(ploidy) - 1))), Inf), labels = c('Deep Deletion', 'Deletion', 'Neutral', 'Amplification', 'Deep Amplification')))),
          base::list(base::as.character(base::cut(base::unlist(copyNumbers.Exons), c(-Inf, base::max(c(ifelse((unique(ploidy)) / 2 <= 0.75, ((unique(ploidy)) / 2) / 2, 0.75), (unique(ploidy)) / 3)), (unique(ploidy)) / 1.5, (1.5 * (unique(ploidy))), min(15, (3 * (unique(ploidy)))), Inf), labels = c('Deep Deletion', 'Deletion', 'Neutral', 'Amplification', 'Deep Amplification'))))
        ),
        
        # Determine number of amplified/deleted exons.
        Consequence.CNA.Exons.Del = base::sum(Consequence.CNA.Exons == 'Deep Deletion'),
        Consequence.CNA.Exons.Amp = base::sum(Consequence.CNA.Exons == 'Deep Amplification'),
        
        # Determine min. and max. CN of all exons.
        copyNumber.ExonMin = base::min(copyNumbers.Exons),
        copyNumber.ExonMax = base::max(copyNumbers.Exons),
        
        copyNumbers.Exons = NULL,
        Consequence.CNA.Exons = NULL,
      ) %>% dplyr::ungroup()
    
    return(GENCODE.Overlap.perGene)
    
  }, cl = nThreads)
  
  # Aggregate the samples.
  combinedReport.Copynumbers <- dplyr::bind_rows(combinedReport.Copynumbers)
  
  combinedReport.Copynumbers <- combinedReport.Copynumbers %>% dplyr::mutate(
    Consequence.CNA.Exons.Amp = base::ifelse(is.na(Consequence.CNA.Exons.Amp), 0, Consequence.CNA.Exons.Amp),
    Consequence.CNA.Exons.Del = base::ifelse(is.na(Consequence.CNA.Exons.Del), 0, Consequence.CNA.Exons.Del)
  )
  
  
  # Generate gene-level overview of somatic mutations -----------------------
  
  base::sprintf('Generating per-sample gene-level overview of protein-coding somatic aberrations.') %>% ParallelLogger::logInfo()
  
  combinedReport.Muts <- somaticData$somaticVariants %>%
    
    # Cleanup HGVS.p
    dplyr::mutate(
      HGVSp = trimws(HGVSp),
      HGVSp = ifelse(is.na(HGVSp) | HGVSp == '', NA, HGVSp)
    ) %>% 
    
    # Only retain the protein-coding somatic mutations.
    dplyr::filter(!is.na(HGVSp), base::grepl('protein_coding', BIOTYPE), Consequence != 'synonymous_variant') %>%
    dplyr::mutate(Existing_variation = as.character(Existing_variation))
  
  
  # Replace all empty string characters with NA.
  combinedReport.Muts[combinedReport.Muts == ''] <- NA
  
  combinedReport.Muts <- combinedReport.Muts %>% 
    # Per sample, summarize the gene-level information.
    dplyr::group_by(Gene, SYMBOL, sample) %>%
    dplyr::summarise(
      # Total number of non-synonymous mutations.
      totalNonSynMutsInGeneInSample = n(),
      
      # Clean-up mutational status.
      Consequence.Mut = base::ifelse(totalNonSynMutsInGeneInSample > 1, 'Multiple coding mutations', Hmisc::capitalize(base::gsub('_', ' ', base::gsub('&.*', '', as.character(Consequence))))),
      Consequence.HGVSp = base::paste0(base::unique(base::gsub('%3D', '=', HGVSp)), collapse = ', '),
      
      # Add additional information.
      CLIN_SIG = base::paste0(base::unique(stats::na.omit(CLIN_SIG)), collapse = ', '),
      IMPACT = base::paste0(base::unique(stats::na.omit(IMPACT)), collapse = ', '),
      mutType = base::paste0(base::unique(stats::na.omit(mutType)), collapse = ', '),
      BIOTYPE = base::paste0(base::unique(stats::na.omit(BIOTYPE)), collapse = ', '),
      Existing_variation = base::paste0(base::unique(stats::na.omit(Existing_variation)), collapse = ', '),
      PURPLE_AF = base::paste0(base::unique(stats::na.omit(PURPLE_AF)), collapse = ', ')
      
    ) %>%
    dplyr::ungroup() %>%
    # Make use of uniform column names.
    dplyr::mutate(ENSEMBL = Gene, Gene = NULL) %>%
    base::droplevels()
  
  # Determine optimal ENSEMBL / SYMBOL.
  geneInfo <- tibble::as_tibble(S4Vectors::mcols(GENCODE.v38)) %>% dplyr::distinct(ENSEMBL.Job = ENSEMBL, SYMBOL.Job = SYMBOL)
  
  combinedReport.Muts <- combinedReport.Muts %>% dplyr::left_join(geneInfo, by = c('ENSEMBL' = 'ENSEMBL.Job'))
  combinedReport.Muts <- combinedReport.Muts %>% dplyr::left_join(geneInfo, by = c('SYMBOL' = 'SYMBOL.Job'))
  
  combinedReport.Muts <- combinedReport.Muts %>% dplyr::mutate(
    ENSEMBL = base::ifelse(!is.na(ENSEMBL.Job), ENSEMBL.Job, ENSEMBL), ENSEMBL.Job = NULL,
    SYMBOL = base::ifelse(!is.na(SYMBOL.Job), SYMBOL.Job, SYMBOL), SYMBOL.Job = NULL
  )
  
  
  # Combine somatic aberrations ---------------------------------------------
  
  base::sprintf('Combining the gene-level/per-sample mutational and CNA data.') %>% ParallelLogger::logInfo()
  
  combinedReport.Final <- combinedReport.Copynumbers %>%
    dplyr::left_join(combinedReport.Muts, by = c('sample', 'ENSEMBL'))
  
  # For non-overlapping records, deduce the correct SYMBOL.
  combinedReport.Final$SYMBOL <- base::ifelse(is.na(combinedReport.Final$SYMBOL.x), combinedReport.Final$SYMBOL.y, combinedReport.Final$SYMBOL.x)
  combinedReport.Final$SYMBOL.x <- NULL; combinedReport.Final$SYMBOL.y <- NULL
  
  combinedReport.Final <- combinedReport.Final %>% base::droplevels()
  
  
  # Add driver status. ------------------------------------------------------
  
  base::sprintf('Adding driver status based on HMF/LINX determination and presence in our custom driver-list.') %>% ParallelLogger::logInfo()
  
  # Add status of HMF driver-determination (determined per sample).
  if(base::nrow(somaticData$driverCatalogHMF) > 0){
    combinedReport.Final <- combinedReport.Final %>% dplyr::left_join(somaticData$driverCatalogHMF %>% dplyr::distinct(ENSEMBL, sample, event.HMF = driver), by = c('ENSEMBL', 'sample'))
  }
  # Add presence within our own combined driver list.
  driverList <- tibble::as_tibble(S4Vectors::mcols(driverList)) %>% dplyr::distinct(ENSEMBL, DriverDatabases = SOURCE)
  combinedReport.Final <- combinedReport.Final %>% dplyr::left_join(driverList, by = c('ENSEMBL'))
  
  
  # Add GISTIC2 status. ------------------------------------------------------
  
  base::sprintf('Assessing which genes were detected within GISTIC2 peaks.') %>% ParallelLogger::logInfo()
  
  gisticGenes.perSample <- tibble::as_tibble(S4Vectors::mcols(somaticData$GISTIC2)) %>%
    dplyr::select(GISTIC2Peak = Unique.Name, overlapGenes.GENCODE, dplyr::contains('CPCT', ignore.case = FALSE), dplyr::contains('DRUP', ignore.case = FALSE), dplyr::contains('WIDE', ignore.case = FALSE)) %>%
    dplyr::mutate(overlappingGenes = base::gsub(' \\(.*', '', overlapGenes.GENCODE)) %>%
    dplyr::mutate(SYMBOL = base::strsplit(overlappingGenes, ', '), overlappingGenes = NULL, overlapGenes.GENCODE = NULL) %>%
    tidyr::unnest(SYMBOL) %>%
    tidyr::pivot_longer(cols = tidyr::starts_with(c('DRUP', 'CPCT', 'WIDE')), names_to = 'sample', values_to = 'Peak.Amplitude') %>%
    dplyr::mutate(
      Peak.Amplitude = ifelse(Peak.Amplitude == 0, 'Low amplitude; t > -0.3', Peak.Amplitude),
      Peak.Amplitude = ifelse(Peak.Amplitude == 1, 'Med. amplitude; -0.3 > t > -1.3', Peak.Amplitude),
      Peak.Amplitude = ifelse(Peak.Amplitude == 2, 'High amplitude; t < -1.3', Peak.Amplitude)
    )
  
  # Add per-sample GISTIC2 data.
  combinedReport.Final <- combinedReport.Final %>% dplyr::left_join(gisticGenes.perSample, by = c('sample', 'SYMBOL'))
  
  
  # Add dN/dS status. -------------------------------------------------------
  
  base::sprintf('Assessing which genes were detected by dN/dS.') %>% ParallelLogger::logInfo()
  
  combinedReport.Final <- combinedReport.Final %>% dplyr::mutate(dNdS = base::ifelse(ENSEMBL %in% somaticData$dNdS$ENSEMBL, 'Significant', NA))
  
  
  # Determine structural variant in coding regions --------------------------
  
  base::sprintf('Determining overlap of structural variants with coding regions.') %>% ParallelLogger::logInfo()
  
  combinedReport.Final <- combinedReport.Final %>%
    dplyr::group_by(sample, SYMBOL, ENSEMBL) %>%
    dplyr::mutate(Consequence.SV = base::ifelse((Consequence.CNA.Exons.Del > 0 | Consequence.CNA.Exons.Amp > 0) & !grepl('Deep', Consequence.CNA), 'Structural Variant', NA)) %>%
    dplyr::ungroup()
  
  
  # Determine mutants. ------------------------------------------------------
  
  base::sprintf('Determining mutants, per gene and sample.') %>% ParallelLogger::logInfo()
  
  combinedReport.Final <- combinedReport.Final %>%
    dplyr::mutate(
      isMutant = base::ifelse(
        base::grepl('Deep', Consequence.CNA)|
          !is.na(Consequence.SV)|
          (!is.na(Consequence.Mut) & Consequence.Mut != 'Splice region variant'), TRUE, FALSE),
      geneId = base::paste(SYMBOL, ENSEMBL, sep = ' - ')
    )
  
  
  # Filter on mutants-only --------------------------------------------------
  
  if(mutantsOnly){
    base::sprintf('Retaining mutants-only.') %>% ParallelLogger::logInfo()
    
    combinedReport.Final <- combinedReport.Final %>%
      dplyr::filter(isMutant | !is.na(event.HMF)) %>%
      base::droplevels()
  }
  
  
  # Return statement --------------------------------------------------------
  
  return(combinedReport.Final)
  
}
