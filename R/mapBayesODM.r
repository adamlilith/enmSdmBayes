#' Make display maps for a species
#'
#' Make display maps for a species.
#' @param species Name of species
#' @param family Family of species
#' @param focus2SpEa Spatial data frame object
#' @param nam1SpEa Spatial object with states, does not have to be cropped to focus2SpEa
#' @param nam2SpEa Spatial object with counties AND counts for each species, does not have to be cropped to focus2SpEa
#' @param mcmc MCMC list object from NIMBLE
#' @param minUnconverged Proportion of unconverged parameters that is acceptable for a model
#' @param targetInstit Institutional code(s)
#' @param isOk TRUE/FALSE if model is sufficiently converged
#' @param outDir Directory to which to save image
#' @param ... Arguments to pass to text()
#' @return Write an image to disk.
#' @export

mapBayesODM <- function(
	species,
	family,
	focus2SpEa,
	nam1SpEa,
	nam2SpEa,
	mcmc,
	minUnconverged,
	targetInstit,
	isOk,
	outDir,
	...
) {

	### formatting
	occCol <- 'forestgreen' # green--map color for occupancy
	occUncertCol <- 'red' # map color for uncertainty in occupancy
	knownOccHomeInstitCol <- 'steelblue3' # map color for counties with known occurrences
	knownOccAwayInstitCol <- 'steelblue1' # map color for counties with known occurrences
	knownOccCol <- '#984ea3' # map color for counties with known occurrences
	detectCol <- 'darkorange' # map color for probability of false absence
	falseAbsCol <- 'darkblue' # map color for probability of false absence
	effortCol <- '#6a51a3' # blue--map color number of collections
	
	lwd <- 0.4 # line width of political borders

	speciesColName <- tolower(gsub(species, pattern=' ', replacement='_'))
	instColName <- paste0(speciesColName, '_', paste(targetInstit, collapse=''))
	
	### assess convergence
	######################

		# if (TRUE) {
			# ignore <- c('z[[]', 'falseDetect[[]', 'constraintFalseDetect[[]')
			# conv <- rhatStats(mcmc$samples, rhatThresh = 1.1, minConv = minUnconverged, ignore = ignore)
		# } else {
			# if (FALSE) conv <- list(); conv$sufficient <- TRUE
		# }

	### prepare plot elements
	#########################

		# centroids of counties with occurrences
		occs2SpEa <- nam2SpEa[nam2SpEa@data[ , speciesColName] > 0, ]
		centsSpEa <- rgeos::gCentroid(occs2SpEa, byid=TRUE)
		
		# counties with occurrences and from/not from home institution
		targetInstit2SpEa <- nam2SpEa[omnibus::isTRUENA(nam2SpEa@data[ , instColName]), ]
		nonTargetInstit2SpEa <- nam2SpEa[nam2SpEa@data[ , speciesColName] > 0 & !omnibus::isTRUENA(nam2SpEa@data[ , instColName]), ]

		# create plotting extent: focal region plus a buffer
		ext <- raster::extent(focus2SpEa)
		extSpEa <- as(ext, 'SpatialPolygons')
		projection(extSpEa) <- raster::projection(nam2SpEa)
		extCentSpEa <- rgeos::gCentroid(extSpEa)
		maxDist_m <- rgeos::gDistance(extCentSpEa, extSpEa, hausdorff=TRUE)
		extSpEa <- rgeos::gBuffer(extSpEa, width=0.05 * maxDist_m)
		ext <- raster::extent(extSpEa)
		extSpEa <- as(ext, 'SpatialPolygons')
		projection(extSpEa) <- raster::projection(nam2SpEa)

		nam1SpEaCrop <- raster::crop(nam1SpEa, extSpEa)

	### plot
	########

		ok <- if (isOk) { 'ok' } else { 'notOk'}
		filename <- paste0(tolower(family), '_', speciesColName)
		png(paste0('./', outDir, '/', ok, '/', filename, '.png'), width=1800, height=1100, res=300)
		
			# layout
			par(bg='white', oma=c(0.7, 0, 1.3, 0), mar=c(0, 2, 0, 0))
			lay <- matrix(c(1, 1, 2, 1, 1, 3), nrow=2, byrow=TRUE)
			layout(lay)

		### occupancy
		
			col <- occCol
			
			plot(nam1SpEaCrop, col='gray90', border=NA)
			
			x <- focus2SpEa@data$psi
			cols <- scales::alpha(col, x / max(x))
			plot(focus2SpEa, col='white', border=NA, add=TRUE)
			plot(focus2SpEa, col=cols, border='gray80', lwd=lwd, add=TRUE)
			
			plot(targetInstit2SpEa, col=knownOccHomeInstitCol, border='gray80', lwd=lwd, add=TRUE)
			plot(nonTargetInstit2SpEa, col=knownOccAwayInstitCol, border='gray80', lwd=lwd, add=TRUE)
			
			plot(nam1SpEaCrop, col=NA, border='gray60', lwd=lwd, add=TRUE)
			
			# legend
			maxPsiNonDetect <- max(focus2SpEa@data$psi[focus2SpEa@data$detect == 0])
			labels <- seq(0, 1, by=1/5) * maxPsiNonDetect
			labels <- sprintf('%0.2f', labels)
			labels[1] <- '0'
			
			width <- 0.04
			height <- 0.9
			
			legTitle <- 'Occupancy\nprobability'
			legCex <- 0.62

			swatches <- list(
				list(swatchAdjY=c(0.03, 0.05), col=knownOccHomeInstitCol, border='gray', labels=targetInstit[1]),
				list(swatchAdjY=c(0, 0.02), col=knownOccAwayInstitCol, border='gray', labels='Other')
			)
			
			legendary::legendGrad('left', inset=-0.01, title=legTitle, col=c('white', col), labels=labels, width=width, height=height, border='gray', boxBorder=NA, adjX=c(0.6, 1), adjY=c(0.09, 0.90), titleAdj=c(1.39, 0.95), labAdj=0.82, boxBg=NA, cex=legCex, lwd=0.5 * lwd, swatches=swatches, pos=2)
			
		### uncertainty in occupancy
		
			col <- occUncertCol
			
			plot(nam1SpEaCrop, col='gray90', border=NA)
			
			x <- focus2SpEa$psi95CI
			cols <- scales::alpha(col, x / max(x))
			plot(focus2SpEa, col='white', border=NA, add=TRUE)
			plot(focus2SpEa, col=cols, border='gray80', lwd=lwd, add=TRUE)

			plot(targetInstit2SpEa, col=knownOccHomeInstitCol, border='gray80', lwd=lwd, add=TRUE)
			plot(nonTargetInstit2SpEa, col=knownOccAwayInstitCol, border='gray80', lwd=lwd, add=TRUE)
			
			plot(nam1SpEaCrop, col=NA, border='gray60', lwd=lwd, add=TRUE)
			# points(centsSpEa, pch=16, cex=0.25)
			
			# legend
			labels <- seq(0, 1, by=1 / 5) * max(x[focus2SpEa@data$detect == 0])
			labels <- sprintf('%0.2f', labels)
			labels[1] <- '0'
			
			legCex <- 0.58
			width <- 0.12
			height <- 0.3
			
			legTitle <- 'Occupancy\nuncertainty\n(95% CI)'

			swatches <- list(
				list(swatchAdjY=c(0.05, 0.09), col=knownOccHomeInstitCol, border='gray', labels=targetInstit[1]),
				list(swatchAdjY=c(0, 0.04), col=knownOccAwayInstitCol, border='gray', labels='Other')
			)
			
			legendary::legendGrad('left', inset=-0.01, title=legTitle, col=c('white', col), labels=labels, width=width, height=0.925, border='gray', boxBorder=NA, adjX=c(0, 0.29), adjY=c(0.12, 0.81), titleAdj=c(0.55, 0.93), labAdj=0.2, boxBg=NA, cex=legCex, lwd=0.5 * lwd, swatches=swatches, pos=2)
			
		### effort
		
			title <- expression(paste('Total '*italic('Asclepias')))
			legTitle <- 'Specimens'
			legCex <- 0.35
			col <- effortCol
			
			plot(nam1SpEaCrop, col='gray90', border=NA)
			
			x <- focus2SpEa@data$effort
			x <- log10(x + 1)
			cols <- scales::alpha(col, x / max(x))
			plot(focus2SpEa, col='white', border=NA, add=TRUE)
			plot(focus2SpEa, col=cols, border='gray80', lwd=lwd, add=TRUE)
			plot(nam1SpEaCrop, col=NA, border='gray60', lwd=lwd, add=TRUE)
			points(centsSpEa, pch=16, cex=0.25)
			
			# legend
			labels <- seq(0, 1, by=0.2) * max(x)
			labels <- round(10^labels - 1, digits=1)
			labels <- sprintf('%.1f', labels)
			labels[1] <- '0'
			
			width <- 0.12
			height <- 0.3
			
			legCex <- 0.58
			width <- 0.12
			height <- 0.3
			
			legTitle <- paste('Total\n', genus)

			legendary::legendGrad('left', inset=-0.01, title=legTitle, col=c('white', col), labels=labels, width=width, height=0.925, border='gray', boxBorder=NA, adjX=c(0, 0.29), adjY=c(0, 0.84), titleAdj=c(0.55, 0.93), labAdj=0.2, boxBg=NA, cex=legCex, lwd=0.5 * lwd, pos=2)

		### titles
		
			main <- substitute(paste(family, ': ', italic(species)), env=list(species=species, family=family))
			mtext(text=main, at=c(0.01), outer=TRUE, cex=1, line=-0.3, adj=0)
			if (!isOk) mtext(text='Insufficient convergence', at=c(0.01), outer=TRUE, cex=1, line=-1.3, adj=0, col='red')
			mtext(text=date(), side=1, at=0.99, outer=TRUE, cex=0.20, line=-0.3, adj=1)

			
			if ('q' %in% rownames(mcmc$summary$all.chains)) {
			
				q <- round(mcmc$summary$all.chains['q', 'Mean'], 2)
				qLow <- round(mcmc$summary$all.chains['q', '95%CI_low'], 2)
				qHigh <- round(mcmc$summary$all.chains['q', '95%CI_upp'], 2)

				msg <- paste0('Probability of mistaken identification of all specimens in a county: ', sprintf('%0.3f', q), ' (95% CI: ', sprintf('%0.3f', qLow), '-', sprintf('%0.3f', qHigh), ')')
				mtext(msg, side=1, at=0.01, outer=TRUE, cex=0.5, line=-0.3, adj=0)
			
			}
			
		dev.off()
	
}
