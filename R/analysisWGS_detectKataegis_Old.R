#' @title Detect kataegis-like events by using a PCF algorithm.
#'
#' @description Uses piecewise constant fitting and additional hueristics to detect hypermutated (snvs) regions per chromosome.
#'
#' @param snv (GRanges): GRanges object containing all SNV of a single sample.
#' @param min.snvs (numeric): What is the minimum amount of SNV that a (constant fit) segment should have being calculated on. (How many SNV should have a consectutive fit)
#' @param maxMeanDistance (numeric): The maximum mean interchromosomal distance segmented SNV.
#' @param maxk (numeric): Max. amount of consecutive SNV which may be used to discover segmental breakpoints.
#'
#' @return (GRanges) A GRanges object containing all non-overlapping kataegis sites.
#' @examples
#' \dontrun{
#'
#' # Example case.
#' detectKataegis(mutData, min.snvs = 5, maxMeanDistance = 1E3, maxk = 6)
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @family CPCT
#' @export
detectKataegis <- function (snv, min.snvs = 5, maxMeanDistance = 2E3, maxk = 5){

    # Input validation --------------------------------------------------------

    checkmate::checkClass(snv, 'GRanges')
    checkmate::checkNumeric(min.snvs)
    checkmate::checkNumeric(maxMeanDistance)
    checkmate::checkNumeric(maxk)

    checkmate::assert(
        !is.null(GenomeInfoDb::seqinfo(snv))
    )

    if(length(snv) <= 50){
        sprintf('Sample had fewer than 50 mutations, skipping: %s', unique(snv$sample)) %>% ParallelLogger::logInfo()
        return(GenomicRanges::GRanges())
    }

    # Detect kataegis in parallel ---------------------------------------------

    # Calculate kataegis per chromosome.
    snv.perChr <- base::split(snv, GenomeInfoDb::seqnames(snv))

    sprintf('Detecting kataegis per chromosome for sample: %s', unique(snv$sample)) %>% ParallelLogger::logInfo()

    events.kataegis <- base::unlist(GenomicRanges::GenomicRangesList(pbapply::pblapply(snv.perChr, function(muts){

        # Sort before calculating distance.
        muts <- GenomicRanges::sort(muts)

        # Skip chromosomes which do not even have enough SNV to fill a segment.
        if(base::length(muts) >= min.snvs){

            # Calculate the intergenomic distances to next SNV per segment.
            positions.SNV <- IRanges::start(muts)
            muts$distance <- (abs(positions.SNV - c(positions.SNV[-1], max(positions.SNV))))

            if(maxk > length(muts)) maxk <- length(muts)

            # Perform piecewise constant curve fitting.
            maxseg <- ceiling(length(muts) / maxk)
            if(maxseg > 5000) maxseg <- 5000

            seg <- tilingArray::segment(log10(muts$distance), maxseg = maxseg, maxk = maxk)
            segments <- c(1, seg@breakpoints[[maxseg]], length(muts$distance))

            # Discover enriched segments.
            segments.kataegis <- base::unlist(GenomicRanges::GenomicRangesList(lapply(2:length(segments), function(i) {
                j <- i - 1
                index.i <- segments[i]
                index.j <- segments[j]

                # Do the segments have enough snvs in the segment
                if ((index.i - index.j) >= min.snvs) {

                    # Is the mean distance < 1000 bp or spread over a large region.
                    mean.dist <- mean(muts[index.j:index.i, ]$distance)
                    if (mean.dist <= maxMeanDistance) {
                        return(GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqlevelsInUse(muts), ranges = IRanges::IRanges(IRanges::start(muts[index.j]), IRanges::end(muts[index.i]))))
                    }else{
                        return(GenomicRanges::GRanges())
                    }
                }else{
                    return(GenomicRanges::GRanges())
                }
            })))
        }else{
            return(GenomicRanges::GRanges())
        }

        # Reduce overlapping segments.
        segments.kataegis.reduced <- IRanges::reduce(IRanges::reduce(segments.kataegis), min.gapwidth = 50000)
        if(length(segments.kataegis.reduced) > 0) segments.kataegis.reduced$sample <- unique(snv$sample)

        return(segments.kataegis.reduced)
    })))


    # Calculate information on kataegis regions -------------------------------

    if(base::length(events.kataegis) >= 0){

        # Retrieve SNV from sample and kataegis site.
        snvInKataegis <- snv[S4Vectors::subjectHits(IRanges::findOverlaps(events.kataegis, snv)),]
        snvInKataegis <- cleanSeqlevels(snvInKataegis)
        snvInKataegis$ALT <- Biostrings::DNAStringSetList(as.list(VariantAnnotation::alt(snvInKataegis)))
        snvInKataegis$REF <- VariantAnnotation::ref(snvInKataegis)
        snvInKataegis$sample <- 'test'

        # Determine mutational type.
        snvInKataegis$mutType <- MutationalPatterns::mut_type(snvInKataegis)
        mutContext <- Biostrings::DNAStringSet(MutationalPatterns::mut_context(snvInKataegis, 'BSgenome.Hsapiens.UCSC.hg19'))
        mutContext.TpC <- sum(grepl('TCA|TCT|TGA|AGA', mutContext))

        # Determine (APOBEC) signatures of kataegis events.
        sample.signatures <- R2CCBC::performMutSigMatch.Alexandrov(GRangesList(snvInKataegis), versionCOSMIC = 'v3', clust.method = 'none')

        #sample.signatures <- R2CPCT::fitMutSigs(snvInKataegis, restrictiveFit = F)

        # Combine data.
        allKatInfoSample <- tibble::tibble(
            sample = unique(sample$sample),
            nSNV = length(snvInKataegis),
            nCT = table(snvInKataegis$mutType)['C>T'],
            nCG = table(snvInKataegis$mutType)['C>G'],
            nCA = table(snvInKataegis$mutType)['C>A'],
            nTA = table(snvInKataegis$mutType)['T>A'],
            nTC = table(snvInKataegis$mutType)['T>C'],
            nTG = table(snvInKataegis$mutType)['T>G'],
            nKataegisEvents = base::length(events.kataegis),
            mutContext.TpC = mutContext.TpC,
            widthKataegis = sum(IRanges::width(snvInKataegis)),
            sig2.relative = sample.signatures$fittedAlexandrov$relativeContribution[c(2)],
            sig2.absolute = sample.signatures$fittedAlexandrov$contribution[c(2)],
            sig13.relative = sample.signatures$fittedAlexandrov$relativeContribution[c(13)],
            sig13.absolute = sample.signatures$fittedAlexandrov$contribution[c(13)]
        )

        # Combine data.
        cohortKataegis <- list()
        cohortKataegis$locationFoci <- events.kataegis
        cohortKataegis$characteristicsFoci <- allKatInfoSample

        # Return statement --------------------------------------------------------

        # Return the combined list of kataegis events.
        return(cohortKataegis)

    }else{
        return(list(locationFoci = NULL, characteristicsFoci = NULL))
    }
}