#' Calculate statistics for a BayesODM model
#'
#' This function calculates statistics for a BayesODM model.
#' @param mcmc MCMC object from \code{trainBayesODM}.
#' @param focusModel SpatialPolygonsDataFrame with BayesODM model output.
#' @param rhat Maximum value of rhat by which to conclude a model parameter is sufficiently converged.
#' @param ... Arguments to pass to \code{rhatStats}.
#' @return
#' @examples
#' @export
bayesODMStats <- function(mcmc, focusModel, rhat = 1.1, ...) {
	
	meta <- attr(focusModel, 'meta')
	
	# # input stats
	# detects <- focusModel@data[ , meta$detect]
	# efforts <- focusModel@data[ , meta$effort]
	
	# output stats
	totalSamples <- (meta$niter - meta$nburnin) / meta$thin
	totalParams <- nrow(mcmc$summary$all.chains)
	proportConverged <- sum(mcmc$summary$all.chains[ , 'rhat'] < rhat, na.rm=TRUE) / totalParams
	
	stats <- data.frame(
	
		totalCounties = nrow(focusModel),
		# occupiedCounties = sum(detects > 0, na.rm=TRUE),
		# countiesWithNonZeroEffort = sum(efforts > 0, na.rm=TRUE),
		# meanDetectsToEfforts = mean(detects / efforts, na.rm=TRUE),
		# maxDetections = max(detects, na.rm=TRUE),
		# maxEfforts = max(efforts, na.rm=TRUE),
		
		meanPsi = mean(focusModel@data$psi, na.rm=TRUE),
		meanPsi95CI = mean(focusModel@data$psi90CI, na.rm=TRUE),
		meanP = mean(focusModel@data$p, na.rm=TRUE),
		meanP95CI = mean(focusModel@data$p90CI, na.rm=TRUE),
		
		niter=meta$niter,
		nburnin=meta$nburnin,
		nchains=meta$nchains,
		thin=meta$thin,
		totalSamples=totalSamples,

		proportConverged = proportConverged,
		rhat = rhat,
		totalParams = totalParams
	
	)

	stats
	
}
