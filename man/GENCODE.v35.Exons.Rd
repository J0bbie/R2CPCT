% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GENCODE.v35.Exons}
\alias{GENCODE.v35.Exons}
\title{GENCODE.v35.GTF <- rtracklayer::import.gff('/mnt/data/ccbc_environment/general/annotation/hg19/GENCODE/noChrPrefix_VEP_gencode.v35lift37.annotation.gtf')
GENCODE.v35.GTF[GENCODE.v35.GTF$transcript_support_level == 'NA' | is.na(GENCODE.v35.GTF$transcript_support_level)]$transcript_support_level <- NA
GENCODE.v35.GTF[GENCODE.v35.GTF$gene_name == 'GTF2I',]$gene_id <- 'ENSG00000077809'
tmp <- GENCODE.v35.GTF}
\format{
GRanges object containing the cleaned-up non-overlapping exonic regions per gene from GENCODE v35.
}
\source{
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
}
\usage{
data(GENCODE.v35.Exons)
}
\description{
# Subset on selected genes.
data('GENCODE.v35', package = 'R2CPCT')
tmp <- tmp[tmp$gene_id %in% GENCODE.v35$ENSEMBL,]
}
\details{
# Only retain transcripts with at least one non-suspect transcript or are single exons.
tmp <- tmp[!grepl('pseudogene|nonsense_mediated_decay|TEC|retained_intron|non_stop_decay', tmp$transcript_type) | tmp$transcript_type == 'polymorphic_pseudogene',]

# Only retain the exon-level elements.
tmp <- tmp[tmp$type == 'exon',]

# Remove unused columns.
S4Vectors::mcols(tmp) <- S4Vectors::mcols(tmp)[!grepl('remap_|havana_|phase|score|ont|hgnc_id|ccdsid|protein|gene_status|transcript_status', colnames(S4Vectors::mcols(tmp)))]

# Only retain unique exons per gene.
tmp.UniqueExons <- GenomicRanges::GRangesList(pbapply::pblapply(base::unique(tmp$gene_id), function(g){

    # Get the original (overlapping) exons of multiple transcripts of the gene.
    x <- tmp[tmp$gene_id == g,]

    # Make unique non-overlapping exonic regions.
    uniqueExons <- GenomicRanges::reduce(x, with.revmap = T)

    # Retrieve the identifiers of the original exons.
    overlappingExons <- base::unlist(
        base::lapply(uniqueExons$revmap, function(y){

            overlappingExons <- base::sprintf('%s - %s (%s)', x[y]$exon_id, x[y]$transcript_id, x[y]$transcript_type)

            if(base::length(overlappingExons) > 2){
                return(base::paste(overlappingExons[1:2], collapse = ', '))
            }else{
                return(base::paste(overlappingExons, collapse = ', '))
            }
        }
        )
    )

    uniqueExons$overlappingExons <- overlappingExons
    uniqueExons$SYMBOL <- base::unique(x$gene_name)
    uniqueExons$ENSEMBL <- base::unique(x$gene_id)
    uniqueExons$source <- base::paste(base::unique(x$source), collapse = ', ')
    uniqueExons$type <- 'exon'
    uniqueExons$totalOverlappingExons <- base::unlist(base::lapply(uniqueExons$revmap, dplyr::n_distinct))
    uniqueExons$revmap <- NULL

    return(uniqueExons)

}, cl = 60))

tmp.UniqueExons <- base::unlist(tmp.UniqueExons)

# Convert columns to best class.
GenomeInfoDb::seqlevels(tmp.UniqueExons) <- base::paste0('chr', GenomeInfoDb::seqlevels(tmp.UniqueExons))

# Retrieve the common gene name.
tmp.UniqueExons$SYMBOL <- as.character(tmp.UniqueExons$SYMBOL)
commonSYMBOL <- limma::alias2SymbolTable(tmp.UniqueExons$SYMBOL)
tmp.UniqueExons$SYMBOL <- base::ifelse(is.na(commonSYMBOL), tmp.UniqueExons$SYMBOL, commonSYMBOL)

# Save to package.
GENCODE.v35.Exons <- tmp.UniqueExons

usethis::use_data(GENCODE.v35.Exons, overwrite = T)
}
\keyword{datasets}
