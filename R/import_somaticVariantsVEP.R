#' @title Import VEP-annotated VCF files obtained from HMF (Strelka2/Sage).
#'
#' @description Imports the VCF obtained from HMF and further annotated by VEP.
#' The annotation should have been performed using the workflow available at: https://github.com/J0bbie/VariantAnnotation_VEP
#'
#' @param pathVCF (character): Path to the <CPCT>_post_processed.vcf file containing Strelka2 somatic variants.
#' @param passOnly (logical): Only import variant which passed all Strelka2 and PON filters. This also filters variants based on PON filter (>5 PON occurences).
#' @param gnomADe (double): Maximum allele frequency in the gnomADe database.
#' @param gnomADg (double): Maximum allele frequency in the gnomADg database.
#' @param keepAnnotation (logical): Import and interpret the ANN column?
#'
#' @return (VRanges) VRanges object containing the somatic variants and assorted annotations.
#' @examples
#' \donttest{
#'
#' 	# Path to VEP annotated VCF.
#' 	path <- '<CPCT#R>_<CPCT#T>_post_processed.hg19_multianno.vcf'
#'
#' 	sampleX.somaticVariants <- importCPCT.somaticVariants(path, passOnly = TRUE)
#'
#' 	# View variant annotation.
#' 	do.call(rbind, sampleX.somaticVariants$annotation)
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
importSomaticVariantsVEP <- function(pathVCF, passOnly = TRUE, gnomADe = 0.001, gnomADg = 0.005, keepAnnotation = TRUE){
  
  # Input validation --------------------------------------------------------
  
  checkmate::assertAccess(pathVCF, access = 'r')
  checkmate::assertLogical(passOnly)
  checkmate::assertDouble(gnomADe)
  checkmate::assertDouble(gnomADg)
  checkmate::assertLogical(keepAnnotation)
  
  
  # Read VCF ----------------------------------------------------------------
  
  sprintf('Importing VCF: %s', pathVCF) %>% ParallelLogger::logInfo()
  
  # Clean seqlevels and add chromosome information.
  sample.VCF <- VariantAnnotation::readVcfAsVRanges(x = pathVCF, genome = 'hg19', use.names = TRUE) %>%
    R2CPCT::cleanSeqlevels(excludeChr = NULL)
  
  # Only retain a single record of each variant (of the tumor)
  if(base::length(base::unique(VariantAnnotation::sampleNames(sample.VCF))) > 1){
    sample.VCF <- sample.VCF[!base::grepl('R$', VariantAnnotation::sampleNames(sample.VCF)),]
  }
  
  # Remove non-PASS variants. (PON and variant-caller filtering)
  if(passOnly){
    totalPriorPass <- base::length(sample.VCF)
    sample.VCF <- sample.VCF[matrixStats::rowAlls(VariantAnnotation::softFilterMatrix(sample.VCF))]
    sprintf('\tRetaining %s / %s PASS-only somatic variants.', base::length(sample.VCF), totalPriorPass) %>% ParallelLogger::logInfo()
  }
  
  
  # Convert / clean annotations ---------------------------------------------
  
  if(is.null(sample.VCF$ANN)) stop(sprintf('This file does not contain ANN (annotation) column:\t%s', pathVCF))
  sprintf('\tConverting and cleaning annotations') %>% ParallelLogger::logTrace()
  
  # Convert annotation to tibble.
  varInfo <- tibble::as_tibble(S4Vectors::mcols(sample.VCF)) %>% 
    tidyr::separate(ANN, into = c('Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','TSL','APPRIS','CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','UNIPROT_ISOFORM','SOURCE','GENE_PHENO','DOMAINS','HGVS_OFFSET','CLIN_SIG','SOMATIC','PHENO','PUBMED','HGVSp2','gnomADe','gnomADe_AF','gnomADg','gnomADg_AF','ClinVar','ClinVar_CLNDN','ClinVar_CLNHGVS','ClinVar_CLNSIG'), sep = '\\|') %>% 
    dplyr::mutate(
      gnomADg_AF = as.numeric(gnomADg_AF),
      gnomADe_AF = as.numeric(gnomADe_AF),
      HGVSp2 = NULL
    )
  
  # Remove unused / old annotations (or also present within the ANN fields).
  varInfo$LOF <- NULL
  varInfo$NMD <- NULL
  varInfo$SEC <- NULL
  varInfo$SEW <- NULL
  varInfo$KT <- NULL
  varInfo$DOMAINS <- NULL
  varInfo$SOURCE <- NULL
  
  varInfo <- varInfo[!grepl('^RC_|^RABQ|^RAD|^ClinVar', colnames(varInfo))]
  
  # Add common gene identifier.
  commonSYMBOL <- limma::alias2SymbolTable(varInfo$SYMBOL)
  varInfo$SYMBOL <- base::ifelse(is.na(commonSYMBOL), varInfo$SYMBOL, commonSYMBOL)
  
  # Convert all columns.
  varInfo <- varInfo %>% dplyr::mutate(dplyr::across(.cols = tidyselect::everything(), .fns = utils::type.convert, as.is = TRUE))
  
  # Convert columns with >50 levels into characters.
  varInfo <- varInfo %>% dplyr::mutate(dplyr::across(.cols = base::colnames(varInfo[base::vapply(varInfo, FUN = function(x) base::length(base::levels(x)), FUN.VALUE = 1) >= 50]), .fns = base::as.character))
  
  # Return annotation to the VRanges.
  S4Vectors::mcols(sample.VCF) <- varInfo
  
  # gnomAD filtering --------------------------------------------------------
  
  # Filter on gnoMAD thresholds.
  totalPriorPass <- base::length(sample.VCF)
  sample.VCF <- sample.VCF[(sample.VCF$gnomADe_AF <= gnomADe | is.na(sample.VCF$gnomADe_AF))]
  sample.VCF <- sample.VCF[(sample.VCF$gnomADg_AF <= gnomADg | is.na(sample.VCF$gnomADg_AF))]
  sprintf('\tRetaining %s / %s somatic variants after gnomAD exome and genome filtering.', base::length(sample.VCF), totalPriorPass) %>% ParallelLogger::logInfo()
  
  # Clean dropped levels.
  S4Vectors::mcols(sample.VCF) <- base::droplevels(S4Vectors::mcols(sample.VCF))
  
  
  # Determine mutational type -----------------------------------------------
  
  sample.VCF$mutType <- 'Other'
  sample.VCF$mutType <- ifelse(VariantAnnotation::isSNV(sample.VCF), 'SNV', sample.VCF$mutType)
  sample.VCF$mutType <- ifelse(VariantAnnotation::isIndel(sample.VCF), 'InDel', sample.VCF$mutType)
  sample.VCF$mutType <- ifelse(!VariantAnnotation::isSNV(sample.VCF) & VariantAnnotation::isSubstitution(sample.VCF), 'MNV', sample.VCF$mutType)
  sample.VCF$mutType <- ifelse(VariantAnnotation::isDelins(sample.VCF), 'DelIn', sample.VCF$mutType)
  
  sample.VCF$mutType <- base::factor(sample.VCF$mutType)
  
  
  # Return statement --------------------------------------------------------
  
  # Remove annotation if not needed (reduced size)
  if(!keepAnnotation) S4Vectors::mcols(sample.VCF) <- NULL
  
  # Add sample name
  sample.VCF$sample <- base::factor(base::gsub('\\.purple.*', '', base::basename(pathVCF)))
  
  sprintf('\tReturning VRanges') %>% ParallelLogger::logTrace()
  
  return(sample.VCF)
  
}
