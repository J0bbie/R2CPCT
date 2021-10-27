#' @title Generates an overview of the GISTIC2 copy-number landscape.
#'
#' @param data.GISTIC (list): Output from \code{\link{importGISTIC2}}.
#' @param maxG (numeric): Cut-off of the max. G-Score.
#' @param minG (numeric): Cut-off of the min. G-Score.
#' @param maxQ (numeric): Which peaks are plotted (max. q-value).
#' @param sizeTrackCN (numeric): Size of the CN track.
#'
#' @return (NULL) Plots a Circos figure.
#' @examples
#' \donttest{
#'
#'  svg(file = 'overviewCNA.svg', width = 8, height = 8, onefile = TRUE)
#'  plotGISTIC2(data.GISTIC, 2, ,-2, 0.1)
#'  dev.off()
#'
#' }
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @export
plotOverviewGISTIC2 <- function(data.GISTIC, maxG = 3, minG = -2, maxQ = 0.1, sizeTrackCN = 15){
  
  # Input validation --------------------------------------------------------
  
  checkmate::assertList(data.GISTIC)
  checkmate::assertNumeric(maxG)
  checkmate::assertNumeric(minG)
  checkmate::assertNumeric(maxQ)
  checkmate::assertNumeric(sizeTrackCN)
  
  
  # Plot --------------------------------------------------------------------
  
  futile.logger::flog.info('Plotting GISTIC2 profile for %s samples', dplyr::n_distinct(data.GISTIC$gisticBroadScoresPerArm$variable))
  
  # Transform peaks to GRanges and clean-up.
  cnData <- GenomicRanges::makeGRangesFromDataFrame(data.GISTIC$gisticPeakScores, keep.extra.columns = TRUE)
  cnData <- GenomeInfoDb::renameSeqlevels(cnData, GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(cnData), 'UCSC'))
  GenomeInfoDb::seqlevels(cnData) <- base::gsub('chr23', 'chrX', GenomeInfoDb::seqlevels(cnData))
  GenomeInfoDb::seqlevels(cnData) <- base::gsub('chr24', 'chrY', GenomeInfoDb::seqlevels(cnData))
  cnData <- R2CPCT::cleanSeqlevels(cnData, excludeChr = NULL)
  
  # Determines which rows are plotted.
  cnData$`G-score`[cnData$Type == 'Amp' & cnData$`G-score` > maxG ] <- maxG
  cnData$`G-score`[cnData$Type == 'Del'] <- -1 * cnData$`G-score`[cnData$Type == 'Del']
  cnData$`G-score`[cnData$Type == 'Del' & cnData$`G-score` < minG ] <- minG
  
  # GISTIC segments.
  gisticData <- data.frame(chr = GenomeInfoDb::seqnames(cnData), start = IRanges::start(cnData), end = IRanges::end(cnData), gScore = cnData$`G-score`, Type = cnData$Type, meanQ = cnData$`-log10(q-value)`)
  
  
  # Combine similar segments ------------------------------------------------
  
  futile.logger::flog.info('Combining similar segments (might take some time)')
  
  gisticData.GR <- GenomicRanges::makeGRangesFromDataFrame(gisticData, keep.extra.columns = TRUE)
  gisticData.GR$gScore.Rounded <- round(gisticData.GR$gScore, 1)
  gisticData.GR$gScore.Binned <- as.vector(cut(gisticData.GR$gScore, breaks = c(0.01, seq(-100.21, 100, .2))))
  
  gisticData.GR <- GenomicRanges::sort(base::unlist(GenomicRanges::GRangesList(base::lapply(base::split(gisticData.GR, paste(gisticData.GR$gScore.Binned, gisticData.GR$Type)), function(x){
    
    x.Reduced <- GenomicRanges::reduce(x, drop.empty.ranges = TRUE, min.gapwidth = 1000, with.revmap = TRUE)
    x.Reduced$gScore <- mean(x$gScore.Rounded)
    x.Reduced$Type <- unique(x$Type)
    
    return(x.Reduced)
    
  }))))
  
  # Convert back to data.frame for plotting.
  gisticData.df <- data.frame(chr = GenomeInfoDb::seqnames(gisticData.GR), start = IRanges::start(gisticData.GR), end = IRanges::end(gisticData.GR), gScore = gisticData.GR$gScore, Type = gisticData.GR$Type)
  
  peakData <- data.frame(chr = GenomeInfoDb::seqnames(data.GISTIC$gisticNarrowPeaksWithAnno), start = IRanges::start(data.GISTIC$gisticNarrowPeaksWithAnno), end = IRanges::end(data.GISTIC$gisticNarrowPeaksWithAnno), genes = sprintf('Peak %s; %s', trimws(gsub('.*[a-z] ', '', data.GISTIC$gisticNarrowPeaksWithAnno$`Unique Name`)), base::trimws(gsub('\\(.*', '', base::as.character(data.GISTIC$gisticNarrowPeaksWithAnno$overlapGenes.Final)))), type = ifelse(grepl('Ampl', data.GISTIC$gisticNarrowPeaksWithAnno$`Unique Name`), '#016300FF', '#1673B4FF'), qVal = data.GISTIC$gisticNarrowPeaksWithAnno$`q values`, stringsAsFactors = FALSE)
  
  # Clear old plots and settings.
  circlize::circos.clear()
  
  # Start with a small gap after the last chromosome to allow for y-axis notations.
  circlize::circos.par(start.degree = 90, gap.after = c(rep(1, length(GenomeInfoDb::seqlevels(cnData)) - 1), 10))
  circlize::circos.initializeWithIdeogram(ideogram.height = 0.03, species = 'hg19', chromosome.index = GenomeInfoDb::seqlevels(cnData))
  
  # Draw segments.
  circlize::circos.genomicTrackPlotRegion(gisticData.df, numeric.column = 4, ylim = c(minG, maxG), panel.fun = function(region, value, ...) {
    
    # Add segments as lines.
    circlize::circos.genomicRect(
      region, value[[1]],
      ytop = value[[1]],
      ybottom = 0,
      col = ifelse(value[[1]] > 0, '#8DCB8E', '#8F9FE8'),
      border = ifelse(value[[1]] > 0, '#8DCB8E', '#8F9FE8')
    )
    
    
    # Add axis lines.
    circlize::circos.lines(circlize::CELL_META$cell.xlim, c(minG, minG), lty = 2, col = '#00000040')
    circlize::circos.lines(circlize::CELL_META$cell.xlim, c(minG / 2, minG / 2), lty = 2, col = '#00000040')
    circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0, 0), lty = 2, col = '#00000040')
    circlize::circos.lines(circlize::CELL_META$cell.xlim, c(maxG / 2, maxG / 2), lty = 2, col = '#00000040')
    circlize::circos.lines(circlize::CELL_META$cell.xlim, c(maxG, maxG), lty = 2, col = '#00000040')
    
  }, track.height = circlize::uh(sizeTrackCN))
  
  # Add CN axis.
  circlize::circos.yaxis(side = 'left', at = c(minG, minG / 2, 0, maxG / 2, maxG), labels = c(minG, minG / 2 , 0, maxG / 2, maxG), sector.index = circlize::get.all.sector.index()[1], labels.cex = 0.3)
  
  # Add gene identifiers to q-val regions.
  circlize::circos.genomicPosTransformLines(peakData, posTransform = circlize::posTransform.default, horizontalLine = 'top', track.height = 0.1)
  
  circlize::circos.genomicTrackPlotRegion(peakData, ylim = c(0, 1), panel.fun = function(region, value, ...) {
    circlize::circos.genomicText(region, value, y = 1, adj = c(0, 0.5), labels = value[[1]], facing = 'reverse.clockwise', niceFacing = TRUE, posTransform = circlize::posTransform.default, cex = .35, col = value[[2]])
  }, bg.border = NA)
  
  
  # Return statement --------------------------------------------------------
  
  return(NULL)
  
}
