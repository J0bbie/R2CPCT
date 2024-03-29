% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GENCODE.v38}
\alias{GENCODE.v38}
\title{# Import GENCODE v38 annotation.
GENCODE.v38.GTF <- rtracklayer::import.gff('/mnt/onco0002/repository/software/ensembl-vep/Plugins/GRCh37/noChrPrefix_gencode.v38lift37.annotation.gtf.bgz')
GENCODE.v38.GTF[GENCODE.v38.GTF$transcript_support_level == 'NA' | is.na(GENCODE.v38.GTF$transcript_support_level)]$transcript_support_level <- NA
GENCODE.v38.GTF[GENCODE.v38.GTF$gene_name == 'GTF2I',]$gene_id <- 'ENSG00000077809'
tmp <- GENCODE.v38.GTF}
\format{
GRanges object containing the cleaned-up genes from GENCODE v38.
}
\source{
http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz
}
\usage{
data(GENCODE.v38)
}
\description{
# Remove pseudogenes and assorted non-interesting types.
tmp <- tmp[!grepl('pseudogene|nonsense_mediated_decay|TEC', tmp$gene_type) | tmp$gene_type == 'polymorphic_pseudogene',]
}
\details{
# Remove clone-based genes.
tmp <- tmp[!grepl('^AC[0-9]|^AP00|^RP11-', tmp$gene_name)]

# Determine number of unique exons per gene.
exonInfo <- tibble::as_tibble(S4Vectors::mcols(tmp[,c('exon_id', 'gene_id', 'transcript_type')]))
exonInfo <- exonInfo %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(!grepl('pseudogene|nonsense_mediated_decay|TEC|retained_intron|non_stop_decay', transcript_type) | transcript_type == 'polymorphic_pseudogene') %>%
    dplyr::summarize(totalExonsInGene = dplyr::n_distinct(exon_id, na.rm = TRUE))

tmp$totalExonsInGene <- exonInfo[base::match(tmp$gene_id, exonInfo$gene_id),]$totalExonsInGene

# Only retain the correctly-placed genetic elements.
tmp <- tmp[tmp$level %in% 1:2,]

# Only retain genes with at least one non-suspect transcript or are single exons.
tmp <- tmp[tmp$totalExonsInGene >= 1,]

# Only retain the gene-level elements.
tmp <- tmp[tmp$type == 'gene',]

# Remove unused columns.
S4Vectors::mcols(tmp) <- S4Vectors::mcols(tmp)[!grepl('remap_|havana_|phase|score|ont|hgnc_id|ccdsid|transcript|exon|protein|gene_status', colnames(S4Vectors::mcols(tmp)))]

# Convert columns.
tmp$source <- as.character(tmp$source)
tmp$type <- as.character(tmp$type)
tmp$ENSEMBL <- as.character(tmp$gene_id); tmp$gene_id <- NULL
tmp$SYMBOL <- as.character(tmp$gene_name); tmp$gene_name <- NULL

# Convert to factor.
S4Vectors::mcols(tmp) <- apply(S4Vectors::mcols(tmp), 2, function(x) factor(as.character(x)))
S4Vectors::mcols(tmp) <- droplevels(S4Vectors::mcols(tmp))
GenomeInfoDb::seqlevels(tmp) <- base::paste0('chr', GenomeInfoDb::seqlevels(tmp))
tmp$ENSEMBL <- as.character(tmp$ENSEMBL)

# Retrieve the common gene name.
tmp$SYMBOL <- as.character(tmp$SYMBOL)
commonSYMBOL <- limma::alias2SymbolTable(tmp$SYMBOL)
tmp$SYMBOL <- base::ifelse(is.na(commonSYMBOL), tmp$SYMBOL, commonSYMBOL)

# Save to package.
GENCODE.v38 <- tmp

usethis::use_data(GENCODE.v38, overwrite = TRUE)
}
\keyword{datasets}
