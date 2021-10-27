#' # Import GENCODE v38 annotation.
#' GENCODE.v38.GTF <- rtracklayer::import.gff('/mnt/onco0002/repository/software/ensembl-vep/Plugins/GRCh37/noChrPrefix_gencode.v38lift37.annotation.gtf.bgz')
#' GENCODE.v38.GTF[GENCODE.v38.GTF$transcript_support_level == 'NA' | is.na(GENCODE.v38.GTF$transcript_support_level)]$transcript_support_level <- NA
#' GENCODE.v38.GTF[GENCODE.v38.GTF$gene_name == 'GTF2I',]$gene_id <- 'ENSG00000077809'
#' tmp <- GENCODE.v38.GTF
#' 
#' # Remove pseudogenes and assorted non-interesting types.
#' tmp <- tmp[!grepl('pseudogene|nonsense_mediated_decay|TEC', tmp$gene_type) | tmp$gene_type == 'polymorphic_pseudogene',]
#' 
#' # Remove clone-based genes.
#' tmp <- tmp[!grepl('^AC[0-9]|^AP00|^RP11-', tmp$gene_name)]
#' 
#' # Determine number of unique exons per gene.
#' exonInfo <- tibble::as_tibble(S4Vectors::mcols(tmp[,c('exon_id', 'gene_id', 'transcript_type')]))
#' exonInfo <- exonInfo %>%
#'     dplyr::group_by(gene_id) %>%
#'     dplyr::filter(!grepl('pseudogene|nonsense_mediated_decay|TEC|retained_intron|non_stop_decay', transcript_type) | transcript_type == 'polymorphic_pseudogene') %>%
#'     dplyr::summarize(totalExonsInGene = dplyr::n_distinct(exon_id, na.rm = TRUE))
#' 
#' tmp$totalExonsInGene <- exonInfo[base::match(tmp$gene_id, exonInfo$gene_id),]$totalExonsInGene
#' 
#' # Only retain the correctly-placed genetic elements.
#' tmp <- tmp[tmp$level %in% 1:2,]
#' 
#' # Only retain genes with at least one non-suspect transcript or are single exons.
#' tmp <- tmp[tmp$totalExonsInGene >= 1,]
#' 
#' # Only retain the gene-level elements.
#' tmp <- tmp[tmp$type == 'gene',]
#' 
#' # Remove unused columns.
#' S4Vectors::mcols(tmp) <- S4Vectors::mcols(tmp)[!grepl('remap_|havana_|phase|score|ont|hgnc_id|ccdsid|transcript|exon|protein|gene_status', colnames(S4Vectors::mcols(tmp)))]
#' 
#' # Convert columns.
#' tmp$source <- as.character(tmp$source)
#' tmp$type <- as.character(tmp$type)
#' tmp$ENSEMBL <- as.character(tmp$gene_id); tmp$gene_id <- NULL
#' tmp$SYMBOL <- as.character(tmp$gene_name); tmp$gene_name <- NULL
#' 
#' # Convert to factor.
#' S4Vectors::mcols(tmp) <- apply(S4Vectors::mcols(tmp), 2, function(x) factor(as.character(x)))
#' S4Vectors::mcols(tmp) <- droplevels(S4Vectors::mcols(tmp))
#' GenomeInfoDb::seqlevels(tmp) <- base::paste0('chr', GenomeInfoDb::seqlevels(tmp))
#' tmp$ENSEMBL <- as.character(tmp$ENSEMBL)
#' 
#' # Retrieve the common gene name.
#' tmp$SYMBOL <- as.character(tmp$SYMBOL)
#' commonSYMBOL <- limma::alias2SymbolTable(tmp$SYMBOL)
#' tmp$SYMBOL <- base::ifelse(is.na(commonSYMBOL), tmp$SYMBOL, commonSYMBOL)
#' 
#' # Save to package.
#' GENCODE.v38 <- tmp
#' 
#' usethis::use_data(GENCODE.v38, overwrite = TRUE)
#' @format GRanges object containing the cleaned-up genes from GENCODE v38.
#' @docType data
#'
#' @usage data(GENCODE.v38)
#' @source http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz
#' @keywords datasets
'GENCODE.v38'

#' GENCODE.v38.GTF <- rtracklayer::import.gff('/mnt/onco0002/repository/software/ensembl-vep/Plugins/GRCh37/noChrPrefix_gencode.v38lift37.annotation.gtf.bgz')
#' GENCODE.v38.GTF[GENCODE.v38.GTF$transcript_support_level == 'NA' | is.na(GENCODE.v38.GTF$transcript_support_level)]$transcript_support_level <- NA
#' GENCODE.v38.GTF[GENCODE.v38.GTF$gene_name == 'GTF2I',]$gene_id <- 'ENSG00000077809'
#' tmp <- GENCODE.v38.GTF
#' 
#' # Subset on selected genes.
#' data('GENCODE.v38', package = 'R2CPCT')
#' tmp <- tmp[tmp$gene_id %in% GENCODE.v38$ENSEMBL,]
#' 
#' # Only retain transcripts with at least one non-suspect transcript or are single exons.
#' tmp <- tmp[!grepl('pseudogene|nonsense_mediated_decay|TEC|retained_intron|non_stop_decay', tmp$transcript_type) | tmp$transcript_type == 'polymorphic_pseudogene',]
#' 
#' # Only retain the exon-level elements.
#' tmp <- tmp[tmp$type == 'exon',]
#' 
#' # Remove unused columns.
#' S4Vectors::mcols(tmp) <- S4Vectors::mcols(tmp)[!grepl('remap_|havana_|phase|score|ont|hgnc_id|ccdsid|protein|gene_status|transcript_status', colnames(S4Vectors::mcols(tmp)))]
#' 
#' # Only retain unique exons per gene.
#' tmp.UniqueExons <- GenomicRanges::GRangesList(pbapply::pblapply(base::unique(tmp$gene_id), function(g){
#'   
#'   # Get the original (overlapping) exons of multiple transcripts of the gene.
#'   x <- tmp[tmp$gene_id == g,]
#'   
#'   # Make unique non-overlapping exonic regions.
#'   uniqueExons <- GenomicRanges::reduce(x, with.revmap = TRUE)
#'   
#'   # Retrieve the identifiers of the original exons.
#'   overlappingExons <- base::unlist(
#'     base::lapply(uniqueExons$revmap, function(y){
#'       
#'       overlappingExons <- base::sprintf('%s - %s (%s)', x[y]$exon_id, x[y]$transcript_id, x[y]$transcript_type)
#'       
#'       if(base::length(overlappingExons) > 2){
#'         return(base::paste(overlappingExons[1:2], collapse = ', '))
#'       }else{
#'         return(base::paste(overlappingExons, collapse = ', '))
#'       }
#'     }
#'     )
#'   )
#'   
#'   uniqueExons$overlappingExons <- overlappingExons
#'   uniqueExons$SYMBOL <- base::unique(x$gene_name)
#'   uniqueExons$ENSEMBL <- base::unique(x$gene_id)
#'   uniqueExons$source <- base::paste(base::unique(x$source), collapse = ', ')
#'   uniqueExons$type <- 'exon'
#'   uniqueExons$totalOverlappingExons <- base::unlist(base::lapply(uniqueExons$revmap, dplyr::n_distinct))
#'   uniqueExons$revmap <- NULL
#'   
#'   return(uniqueExons)
#'   
#' }, cl = 60))
#' 
#' tmp.UniqueExons <- base::unlist(tmp.UniqueExons)
#' 
#' # Convert columns to best class.
#' GenomeInfoDb::seqlevels(tmp.UniqueExons) <- base::paste0('chr', GenomeInfoDb::seqlevels(tmp.UniqueExons))
#' 
#' # Retrieve the common gene name.
#' tmp.UniqueExons$SYMBOL <- as.character(tmp.UniqueExons$SYMBOL)
#' commonSYMBOL <- limma::alias2SymbolTable(tmp.UniqueExons$SYMBOL)
#' tmp.UniqueExons$SYMBOL <- base::ifelse(is.na(commonSYMBOL), tmp.UniqueExons$SYMBOL, commonSYMBOL)
#' 
#' # Save to package.
#' GENCODE.v38.Exons <- tmp.UniqueExons
#' 
#' usethis::use_data(GENCODE.v38.Exons, overwrite = TRUE)
#' @format GRanges object containing the cleaned-up non-overlapping exonic regions per gene from GENCODE v38.
#' @docType data
#'
#' @usage data(GENCODE.v38.Exons)
#' @source http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz
#' @keywords datasets
'GENCODE.v38.Exons'

#' See misc/generateDriverList.R
#' @format GRanges object containing the combined drivers.
#' @docType data
#'
#' @usage data(driverList)
#' @source Custom, see misc/generateDriverList.R.
#' @keywords datasets
'driverList'

#' library(rvest)
#' 
#' # SBS Signatures ----------------------------------------------------------
#' 
#' URL <- 'https://cancer.sanger.ac.uk/cosmic/signatures/SBS/'
#' 
#' # Retrieve content.
#' x <- xml2::read_html(URL) %>%
#'   html_nodes(css = '.vignettes_table_plot') %>%
#'   html_children() %>%
#'   html_text(trim = TRUE)
#' 
#' sigs.SBS <- rbind(
#'   do.call(rbind, lapply(x, function(s){
#'     s <- trimws(strsplit(s, '\n')[[1]])
#'     data.frame(Signature = s[[1]], proposedAetiology = gsub('Proposed Aetiology', '', s[[6]]))
#'   }))
#'   , data.frame(Signature = c('SBS27', 'SBS43', 'SBS45','SBS46','SBS47','SBS48','SBS49','SBS50','SBS51','SBS52','SBS53','SBS54','SBS55','SBS56','SBS57','SBS58','SBS59','SBS60'), proposedAetiology = 'Possible Sequencing Artifact')
#' )
#' 
#' 
#' # DBS Signatures ----------------------------------------------------------
#' 
#' URL <- 'https://cancer.sanger.ac.uk/cosmic/signatures/DBS/'
#' 
#' # Retrieve content.
#' x <- xml2::read_html(URL) %>%
#'   html_nodes(css = '.vignettes_table_plot') %>%
#'   html_children() %>%
#'   html_text(trim = TRUE)
#' 
#' sigs.DBS <- do.call(rbind, lapply(x, function(s){
#'   s <- trimws(strsplit(s, '\n')[[1]])
#'   data.frame(Signature = s[[1]], proposedAetiology = gsub('Proposed Aetiology', '', s[[6]]))
#' }))
#' 
#' # InDel Signatures --------------------------------------------------------
#' 
#' URL <- 'https://cancer.sanger.ac.uk/cosmic/signatures/ID'
#' 
#' # Retrieve content.
#' x <- xml2::read_html(URL) %>%
#'   html_nodes(css = '.vignettes_table_plot') %>%
#'   html_children() %>%
#'   html_text(trim = TRUE)
#' 
#' sigs.ID <- do.call(rbind, lapply(x, function(s){
#'   s <- trimws(strsplit(s, '\n')[[1]])
#'   data.frame(Signature = s[[1]], proposedAetiology = gsub('Proposed Aetiology', '', s[[6]]))
#' }))
#' 
#' 
#' # Combine lists -----------------------------------------------------------
#' 
#' proposedAetiologyCOSMICv3.2 <- rbind(sigs.SBS, sigs.DBS, sigs.ID)
#' 
#' # Group signatures based on common proposed aetiology together.
#' proposedAetiologyCOSMICv3.2 <- proposedAetiologyCOSMICv3.2 %>%
#'   dplyr::group_by(proposedAetiology) %>%
#'   dplyr::mutate(proposedAetiologyGrouped = sprintf('%s (%s)', proposedAetiology, ifelse(dplyr::n_distinct(Signature) <= 4, paste(Signature, collapse = ', '), '>4 signatures'))) %>%
#'   dplyr::ungroup()
#' 
#' usethis::use_data(proposedAetiologyCOSMICv3.2)
#' @format Tibble containing the COSMIC (v3.2; March 2021) proposed signatures aetiologies.
#' @docType data
#'
#' @usage data(proposedAetiologyCOSMICv3.2)
#' @source COSMIC
#' @keywords datasets
'proposedAetiologyCOSMICv3.2'


#' # Import protein domains from PROT2HG --------------------------------------
#'
#' proteinDomains <- data.table::fread(input = '~/prot2hg_1938_112019.csv', sep = ';', stringsAsFactors = FALSE)
#'
#' # Select protein-information.
#' proteinDomains <- proteinDomains %>% dplyr::select(protein_ID, gene, prot_start, prot_end, type, feature_name, note, ensembl) %>% dplyr::distinct()
#'
#' # Convert to GRanges.
#' proteinDomains <- GenomicRanges::makeGRangesFromDataFrame(proteinDomains, keep.extra.columns = TRUE, seqnames.field = 'protein_ID', start.field = 'prot_start', end.field = 'prot_end')
#'
#' # Factorize to reduce size.
#' proteinDomains$type <- factor(proteinDomains$type)
#' proteinDomains$gene <- factor(proteinDomains$gene)
#' proteinDomains$ensembl <- factor(proteinDomains$ensembl)
#'
#' usethis::use_data(proteinDomains)
#' @format GRanges containing protein domains and protein-specific sites, derived from PROT2HG (Jan. 2021)
#' @docType data
#'
#' @usage data(proteinDomains)
#' @source proteins
#' @keywords datasets
'proteinDomains'