#' Is a set of chains sufficiently converged?
#'
#' Determines if a set of MCMC chains is sufficiently converged based on the proportion of coefficients for which rhat is < a threshold value.
#' @param mcmc An object of class \code{\link[coda]{mcmc.list}}.
#' @param rhatThresh Numeric, value of rhat below which a chain is considered to be converged (default = 1.1).
#' @param minConv Numeric in the range [0, 1], proportion of rhat values to consider a set of chains "sufficiently" converged. Typically this is 0 (default).
#' @param ignore Character or character vector, names of nodes to ignore when calculating rhat. These can include regex expressions (e.g., \code{x[[]}). Default is \code{NULL} (calculate rhat for all nodes).
#' @return List with three elements:
#' \itemize{
#' 		\item Logical (if proportion of rhats is below the threshold).
#' 		\item Proportion of nodes that did not converge.
#'		\item Vector of node names that did not converge.
#'		\item Vector of node names with \code{NA} for rhat value.
#' @export

rhatStats <- function(mcmc, rhatThresh = 1.1, minConv = 0, ignore = NULL) {

	nchains <- length(mcmc)
	
	if (!is.null(ignore)) {
		for (chain in 1:nchains) {
			for (ig in ignore) {
				mcmc[[chain]] <- mcmc[[chain]][ , !grepl(colnames(mcmc[[chain]]), pattern=ig)]
			}
		}
	}
	
	rhat <- wiqid::simpleRhat(mcmc, nchains)

	propUnconv <- sum(rhat >= rhatThresh, na.rm=TRUE) / sum(!is.na(rhat))
	sufficient <- (propUnconv < minConv)
	unconv <- names(rhat[rhat > rhatThresh])
	nas <- names(rhat)[is.na(rhat)]

	returns <- list(sufficient=sufficient, propUnconv=propUnconv, unconv=unconv, nas=nas, rhat=rhat)
	returns

}
