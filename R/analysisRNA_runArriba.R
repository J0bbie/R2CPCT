#' @title Run Arriba on RNA-Seq samples using structural-variants information derived from WGS.
#'
#' @param data.SV (GRanges): GRanges of SV derived from WGS.
#' @param minTAF (integer): Min. TAF of SV to take into the analysis, set to zero (0) for no filtering.
#' @param bamFolder (character): Path to folder containing the RNA-Seq BAM files (STAR).
#' @param outputFolder (character): Path to folder to write Arriba output.
#'
#' @return (chr) Command to run Arriba.
#' @examples
#' \dontrun{
#'
#' 	runArriba(
#' 	  data.SV = DR71.CohortWGS$structuralVariants,
#' 	  minTAF = 0.05,
#' 	  bamFolder = '/mnt/data2/hartwig/DR71/Oct2020/dataHMF/RNASeq/BAM/',
#' 	  outputFolder = '/mnt/data2/hartwig/DR71/Oct2020/dataHMF/RNASeq/Arriba/'
#' 	)
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family RNA-Seq
#' @export
runArriba <- function(data.SV, minTAF = 0.05, bamFolder, outputFolder){

    # Input validation --------------------------------------------------------

    checkmate::checkClass(data.SV, 'GRanges')
    checkmate::checkNumeric(minTAF)
    checkmate::checkAccess(bamFolder, access = 'r')
    checkmate::checkAccess(outputFolder, access = 'rw')

    base::sprintf('Performing fusion-gene detection using Arriba on %s samples.', dplyr::n_distinct(data.SV$sample)) %>% ParallelLogger::logInfo()


    # Convert SV to Arriba-friendly format ------------------------------------

    base::sprintf('\t- Converting to Arriba-format.') %>% ParallelLogger::logInfo()

    # Remove single-ends.
    data.SV <- data.SV[data.SV$SVTYPE != 'SINGLE',]

    # Clean remaining factors.
    S4Vectors::mcols(data.SV) <- base::droplevels(S4Vectors::mcols(data.SV))

    # Loop per sample
    arriba.SV <- pbapply::pblapply(base::split(data.SV, data.SV$sample), function(x){

        # Only retain partnered SV.
        x <- x[!is.na(base::match(x$partner, base::gsub('.*\\.', '', base::names(x)))),]

        # Format to Arriba standards.
        x.SV <- data.frame(
            'A' = sprintf('%s:%s', base::gsub('chr', '', base::as.character(GenomeInfoDb::seqnames(x))), as.numeric(GenomicRanges::start(x))),
            'B' = sprintf('%s:%s', base::gsub('.*(\\[|\\])', '', base::gsub(':.*', '', x$ALT)), as.integer(base::gsub('(\\[|\\]).*', '', base::gsub('.*:', '', x$ALT)))),
            'A.strand' = as.character(strand(x)),
            'B.strand' = as.character(BiocGenerics::strand(x[base::match(x$partner, base::gsub('.*\\.', '', base::names(x))),])),
            TAF = x$TAF
        )

        # Filter on min. TAF (if requested)
        x.SV <- x.SV[x.SV$TAF >= minTAF,]

        # Filter out SV which had no proper strand orientation.
        x.SV <- x.SV %>% dplyr::filter(A.strand != '*', B.strand != '*')

        return(x.SV)
    }, cl = 1)


    # Run Arriba --------------------------------------------------------------

    base::sprintf('\t- Generating the command to run Arriba.') %>% ParallelLogger::logInfo()

    # Write Arriba-friendly SVs to temp. dir.
    lapply(names(arriba.SV), function(x){
        z.SV <- arriba.SV[names(arriba.SV) == x][[1]]

        write.table(
            x = z.SV,
            file = file.path(tempdir(), sprintf('%s_SV.txt', x)),
            row.names = F, quote = F, sep = '\t', col.names = F
        )
    })

    # Retrieve the BAM files for which SV are present.
    bamFiles <- data.frame(BAM = list.files(path = bamFolder, pattern = 'sorted.*_markDup.bam$', full.names = T)) %>% dplyr::mutate(sample = gsub('_Aligned.sortedByCoord_markDup.bam', '', basename(as.character(BAM))))
    bamFiles <- bamFiles %>% dplyr::filter(sample %in% names(arriba.SV))

    # Generate the Bash commands and write to tmp. file.
    z <- sprintf('/mnt/data/ccbc_environment/software/general/arriba_v2.1.0/arriba -x %s -g /mnt/data/ccbc_environment/general/annotation/hg19/GENCODE/noChrPrefix_VEP_gencode.v35lift37.annotation.gtf -a /mnt/data/ccbc_environment/general/genomes/hsapiens/hg19_HMF/Homo_sapiens.GRCh37.GATK.illumina.fasta -b /mnt/data/ccbc_environment/software/general/arriba_v2.1.0/database/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz -k /mnt/data/ccbc_environment/general/annotation/ChimerDB_4.0/knownRecurrentFusions.txt -o %s -d %s -s reverse', bamFiles$BAM, paste0(outputFolder, '/', bamFiles$sample, '_Arriba_fusions.tsv'), paste0(tempdir(), '/', bamFiles$sample,'_SV.txt'))
    outputFile <- tempfile()
    write.table(z, outputFile, row.names = F, quote = F, sep = '\t', col.names = F)


    # Return statement --------------------------------------------------------

    base::sprintf('Run the following command in Bash:\ncat %s | parallel -j <threads>', outputFile) %>% ParallelLogger::logInfo()

    return(outputFile)

}


#' @title Import Arriba results.
#'
#' @param inputFolder (character): Path to folder containing Arriba output.
#'
#' @return (tibble) Returns a tibble of Arriba-findings.
#' @examples
#' \dontrun{
#'
#' 	importArriba('/mnt/data2/hartwig/DR71/Oct2020/dataHMF/RNASeq/Arriba/')
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family RNA-Seq
#' @export
importArriba <- function(inputFolder){

    # Input validation --------------------------------------------------------

    checkmate::checkAccess(inputFolder, access = 'r')

    base::sprintf('Importing Arriba results from: %s', inputFolder) %>% ParallelLogger::logInfo()


    # Import results ----------------------------------------------------------

    files.Arriba <- list.files(inputFolder, pattern = '_Arriba_fusions.tsv$', full.names = T)

    results.Arriba <- pbapply::pblapply(files.Arriba, function(x){

        # Import.
        z <- readr::read_tsv(x, col_types = 'cccccccccdddddcccccccccccccccc')

        # Add sample-identifier.
        z$sample <- base::gsub('_.*', '', base::basename(x))

        # Clean columns.
        z <- z %>% dplyr::mutate(
            gene1 = `#gene1`,
            `#gene1` = NULL,
            fusion = sprintf('%s-%s', gene1, gene2),
            read_identifiers = NULL
        )

        # Return.
        return(z)
    })

    results.Arriba <- dplyr::bind_rows(results.Arriba)


    # Return statement --------------------------------------------------------

    return(results.Arriba)

}
