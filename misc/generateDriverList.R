library(dplyr)

# Misc. Functions ---------------------------------------------------------

# Internal function.
coalesce_join <- function(x, y, by = NULL, suffix = c(".x", ".y"), join = dplyr::full_join, ...) {
    joined <- join(x, y, by = by, suffix = suffix, ...)
    # names of desired output
    cols <- union(names(x), names(y))

    to_coalesce <- names(joined)[!names(joined) %in% cols]
    suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
    # remove suffixes and deduplicate
    to_coalesce <- unique(substr(
        to_coalesce,
        1,
        nchar(to_coalesce) - nchar(suffix_used)
    ))

    coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
        joined[[paste0(.x, suffix[1])]],
        joined[[paste0(.x, suffix[2])]]
    ))
    names(coalesced) <- to_coalesce

    dplyr::bind_cols(joined, coalesced)[cols]
}


# Import driver databases -------------------------------------------------

# Generate a list of known drivers from various sources.
driverInputs <- list()

# COSMIC - Cancer Gene Census from https://cancer.sanger.ac.uk/cosmic/download.
driverInputs$COSMIC <- readr::read_csv('/mnt/data/ccbc_environment/general/drivers/COSMICv92_CancerGeneCensus.csv')

# Priestley at al. - Pan-cancer whole-genome analyses of metastatic solid tumours. (Nature)
driverInputs$Priestley <- readxl::read_excel('/mnt/data/ccbc_environment/general/drivers/Priestley_Muts.xlsx') %>% dplyr::filter(qallsubs_cv <= 0.1 | qglobal_cv <= 0.1)
driverInputs$Priestley.CNA <- rbind(readxl::read_excel('/mnt/data/ccbc_environment/general/drivers/Priestley_CNA.xlsx', sheet = 1) %>% dplyr::distinct(gene) %>% dplyr::mutate(alteration = 'amplification'), readxl::read_excel('/mnt/data/ccbc_environment/general/drivers/Priestley_CNA.xlsx', sheet = 2) %>% dplyr::distinct(gene) %>% dplyr::mutate(alteration = 'deletion'))

# Dietlein et al. - Identification of cancer driver genes based on nucleotide context. (Nature Genetics)
driverInputs$Dietlein <- readr::read_tsv('/mnt/data/ccbc_environment/general/drivers/dietleinDrivers.txt')

# Martincorena et al. - Universal Patterns of Selection in Cancer and Somatic Tissues. (Cell)
driverInputs$Martincorena <- readr::read_delim('/mnt/data/ccbc_environment/general/drivers/Martincorena_muts.txt', delim = '\t')

# IntOGen - Release date 2020.02.01
driverInputs$IntOGen <- readr::read_delim('/mnt/data/ccbc_environment/general/drivers/IntOGen_Compendium_Cancer_Genes_022020.tsv', delim = '\t')

# Bailey et al. - Comprehensive Characterization of Cancer Driver Genes and Mutations. (Cell)
driverInputs$Bailey <- readr::read_csv('/mnt/data/ccbc_environment/general/drivers/BaileyEtAl2019.csv')

# Custom genes.
driverInputs$Custom <- data.frame(SYMBOL = c('BCOR', 'MN1', 'ADGRB3', 'CRKL', 'DNMT1', 'DNMT3A', 'MET', 'PIK3CG', 'RICTOR', 'SETDB1', 'SLIT2', 'STK11', 'TET1', 'TET2', 'U2AF1', 'AJUBA', 'CUX1', 'STK11', 'ALK', 'CHD9', 'PHOX2B', 'SMARCB1', 'SMARCA4', 'NF1', 'CDKN2B', 'TSC2', 'MEN1', 'SOX11','ZEB1','ZEB2','TWIST1','TWIST2','SNAI1','SNAI2','IKZF2'), Source = 'Custom')


# Combine and add missing identifiers. ------------------------------------

# Import gene annotation
data('GENCODE.v35', package = 'R2CPCT')
geneInfo <- tibble::as_tibble(S4Vectors::mcols(GENCODE.v35)) %>% dplyr::select(ENSEMBL, SYMBOL) %>% dplyr::distinct()

# Combine.
driverInputs$Combined <- dplyr::bind_rows(
    driverInputs$COSMIC %>% dplyr::distinct(SYMBOL = `Gene Symbol`, Source = 'COSMIC'),
    driverInputs$Priestley %>% dplyr::distinct(SYMBOL = gene_name, Source = 'Priestley et al. (Muts)'),
    driverInputs$Priestley.CNA %>% dplyr::distinct(SYMBOL = gene, Source = 'Priestley et al. (CNA)'),
    driverInputs$Dietlein %>% dplyr::distinct(SYMBOL = Gene, Source = 'Dietlein et al.'),
    driverInputs$Martincorena %>% dplyr::distinct(SYMBOL = gene, Source = 'Martincorena et al.'),
    driverInputs$IntOGen %>% dplyr::distinct(SYMBOL, Source = 'IntOGen 02-2020'),
    driverInputs$Bailey %>% dplyr::distinct(SYMBOL = Gene, Source = 'Bailey et al.'),
    unique(driverInputs$Custom)
) %>% dplyr::filter(!grepl('^RP11-|^AC[0-9]|^AP00|^AL[0-9]', SYMBOL))

# Convert aliases to main symbol.
driverInputs$Combined <- driverInputs$Combined %>%
    dplyr::rowwise() %>%
    dplyr::mutate(SYMBOL = limma::alias2Symbol(SYMBOL, species = 'Hs')) %>%
    dplyr::ungroup()

# Add identifiers.
driverInputs$Combined <- driverInputs$Combined %>% dplyr::left_join(geneInfo)

# Search using annotation package.
ad <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, key = trimws(unique((driverInputs$Combined %>% dplyr::filter(is.na(ENSEMBL)))$SYMBOL)), columns = c('SYMBOL', 'ENSEMBL'), keytype="SYMBOL")
ad <- ad %>% dplyr::filter(!duplicated(SYMBOL), !is.na(ENSEMBL))

driverInputs$Combined <- driverInputs$Combined %>% coalesce_join(ad, by = 'SYMBOL')

# Retrieve missing identifiers using BioMart.
missingId <- trimws(unique((driverInputs$Combined %>% dplyr::filter(is.na(ENSEMBL)))$SYMBOL))
bm <- biomaRt::getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    filters = 'hgnc_symbol',
    values = missingId,
    mart = biomaRt::useMart('ensembl', 'hsapiens_gene_ensembl'),
) %>%
    dplyr::distinct(SYMBOL = hgnc_symbol, ENSEMBL = ensembl_gene_id)%>%
    dplyr::filter(!duplicated(SYMBOL), !is.na(ENSEMBL))

driverInputs$Combined <- driverInputs$Combined %>% coalesce_join(bm, by = 'SYMBOL')

# Remove genes without proper ENSEMBL identifier.
driverInputs$Combined <- driverInputs$Combined %>% dplyr::filter(!is.na(ENSEMBL))


# Add chromosomal locations. ----------------------------------------------

x <- GENCODE.v35[GENCODE.v35$ENSEMBL %in% driverInputs$Combined$ENSEMBL,]

locInfo <- tibble::as_tibble(data.frame(
    chrom = GenomeInfoDb::seqnames(x),
    start = IRanges::start(x),
    end = IRanges::end(x),
    strand = BiocGenerics::strand(x),
    ENSEMBL = x$ENSEMBL
)) %>% dplyr::mutate(chrom = as.character(chrom))

# Find missing genes.
bm <- tibble::as_tibble(biomaRt::getBM(
    attributes = c('chromosome_name', 'start_position', 'end_position', 'strand', 'ensembl_gene_id'),
    filters = 'ensembl_gene_id',
    values = (driverInputs$Combined %>% dplyr::filter(!ENSEMBL %in% locInfo$ENSEMBL))$ENSEMBL,
    mart = biomaRt::useMart('ensembl', 'hsapiens_gene_ensembl', host = 'grch37.ensembl.org'),
)) %>%
    dplyr::filter(!duplicated(ensembl_gene_id)) %>%
    dplyr::mutate(chrom = as.character(chromosome_name)) %>%
    dplyr::select(chrom, start = start_position, end = end_position, strand = strand, ENSEMBL = ensembl_gene_id)

bm$strand <- factor(ifelse(bm$strand == 1, '+', '-'))

locInfo <- rbind(locInfo, bm)

driverInputs$Combined <- driverInputs$Combined %>% dplyr::left_join(locInfo, by = c('ENSEMBL' = 'ENSEMBL'))

# Remove genes for which no locus could be found.
driverInputs$Combined <- driverInputs$Combined %>% dplyr::filter(!is.na(chrom))

# Remove genes which are not present in the GENCODE (v35) GTF.
driverInputs$Combined <- driverInputs$Combined %>% dplyr::filter(ENSEMBL %in% GENCODE.v35$ENSEMBL)


# Generate distinct driver list. ------------------------------------------

driverList <- driverInputs$Combined %>%
    dplyr::group_by(chrom, start, end, strand, ENSEMBL, SYMBOL) %>%
    dplyr::summarise(SOURCE = paste(Source, collapse = ', '), nDatabases = dplyr::n_distinct(Source)) %>%
    dplyr::ungroup()

# Convert to GRanges.
driverList <- GenomicRanges::makeGRangesFromDataFrame(driverList, keep.extra.columns = T)

# Retrieve the common gene name.
commonSYMBOL <- limma::alias2SymbolTable(driverList$SYMBOL)
driverList$SYMBOL <- base::ifelse(is.na(commonSYMBOL), driverList$SYMBOL, commonSYMBOL)

# Add genomic location to driverList.
driverList.Gr <- GENCODE.v35[base::match(driverList$ENSEMBL, GENCODE.v35$ENSEMBL),]
S4Vectors::mcols(driverList.Gr) <- driverList
driverList <- driverList.Gr

# Add to R2CPCT
usethis::use_data(driverList, overwrite = T)
