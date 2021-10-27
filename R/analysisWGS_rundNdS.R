#' @title Perform dN/dS analysis on a set of mutations from multiple CPCT-02 samples.
#'
#' @details The dNdS database was made from ENSEMBL v104 (GENCODE v38) using the following code and following the dN/dS buildref vignette:
#'
#'  # Remove non-standard chromosomes, clone-based genes and CDS which cannot be divided by 3.
#'  ENSEMBLv104 <- readr::read_tsv('~/test/mart_export.txt') %>%
#'   dplyr::filter(!grepl('_', `Chromosome/scaffold name`)) %>%
#'   dplyr::filter(`CDS Length` %% 3 == 0) %>%
#'   dplyr::filter(!is.na(`Genomic coding start`)) %>%
#'   dplyr::filter(!grepl('^AC[0-9]', `Gene name`)) %>%
#'   dplyr::filter(!grepl('^AP00|^RP11-', `Gene name`)) %>%
#'   dplyr::filter(!grepl('\\.', `Gene name`)) %>%
#'   dplyr::select(
#'       'gene.id' = `Gene stable ID`,
#'        'gene.name' = `Gene name`,
#'        'cds.id' = `Protein stable ID`,
#'        'chr' = `Chromosome/scaffold name`,
#'        'chr.coding.start' = `Genomic coding start`,
#'        'chr.coding.end' = `Genomic coding end`,
#'        'cds.start' = `CDS start`,
#'        'cds.end' = `CDS end`,
#'        'length' = `CDS Length`,
#'        'strand' = Strand)
#'  
#' write.table(ENSEMBLv104, file = '~/test/mart_export_filtered.txt', sep = '\t', row.names = FALSE, quote = FALSE)
#' 
#' # Generate refCDS database.
#' pathCDS = '~/test/mart_export_filtered.txt'
#' pathFasta = '/mnt/onco0002/repository/general/genomes/hsapiens/hg19_HMF/Homo_sapiens.GRCh37.GATK.illumina.fasta'
#' dndscv::buildref(cdsfile = pathCDS, genomefile = pathFasta, outfile = 'inst/extdata/refCDS_ENSEMBLv104_HMF.rda', excludechrs='MT', useids = TRUE)
#'
#' @param dataMuts (VRanges): VRanges containing the mutations which will be inputted into dN/dS.
#'
#' @examples
#' \donttest{
#'
#'  data.Cohort <- R2CPCT::importWGSOfCohort(<cpctIds>, <combinedData>)
#'  rundNdS(data.Cohort$somaticVariants)
#'
#' }
#' @return (tibble) Returns a tibble of the dNdS results.
#' @export
rundNdS <- function(dataMuts){

    # Input validation --------------------------------------------------------

    checkmate::assertClass(dataMuts, classes = 'VRanges')

    sprintf('Performing dN/dS analysis on %s unique samples.\nThis can take some minutes.', dplyr::n_distinct(dataMuts$sample)) %>% ParallelLogger::logInfo()


    # Perform dN/dS -----------------------------------------------------------

    # Convert mutations to data.frame and remove chr prefix.
    dataMuts.df <- data.frame(sampleID = dataMuts$sample, chr = as.character(GenomeInfoDb::seqnames(dataMuts)), pos = as.numeric(IRanges::start(dataMuts)), ref = VariantAnnotation::ref(dataMuts), mut = VariantAnnotation::alt(dataMuts))
    dataMuts.df <- dataMuts.df %>% dplyr::mutate(chr = base::gsub('chr', '', chr))

    # Perform dN/dS algorithm.
    output.dNdS <- dndscv::dndscv(dataMuts.df, refdb = system.file('extdata/refCDS_ENSEMBLv104_HMF.rda', package = 'R2CPCT'), outp = 3)

    # Remove large (unused) annotation database.
    output.dNdS$annotmuts <- NULL


    # Combine final results ---------------------------------------------------

    output.dNdS$finalOutput <- output.dNdS$sel_cv %>%
        dplyr::left_join(output.dNdS$sel_loc %>% dplyr::select(gene_name, qall_loc, qmis_loc), by = c('gene_name' = 'gene_name')) %>%
        dplyr::filter(qglobal_cv <= 0.1 | qallsubs_cv <= 0.1 | qtrunc_cv <= 0.1 | qmis_cv <= 0.1 | qall_loc <= 0.1 | qmis_loc <= 0.1) %>%
        dplyr::mutate(SYMBOL = gsub('.*:', '', gene_name), ENSEMBL = gsub(':.*', '', gene_name), gene_name = NULL)

    # Retrieve the common gene name.
    commonSYMBOL <- limma::alias2SymbolTable(output.dNdS$finalOutput$SYMBOL)
    output.dNdS$finalOutput$SYMBOL <- base::ifelse(is.na(commonSYMBOL), output.dNdS$finalOutput$SYMBOL, commonSYMBOL)


    # Return statement --------------------------------------------------------

    return(output.dNdS)

}
