#' Plot a SpatialPointsDataFrame or SpatialPolygonsDataFrame object
#'
#' This function plots a SpatialPointsDataFrame or SpatialPolygonsDataFrame with formatting based on the output from a \code{mcmc.list} object created by \pkg{rjags} or \pkg{R2jags}. It is useful for displaying the output of a Bayesian model where estimated values vary by point/polygon.
#' @param x SpatialPointsDataFrame or SpatialPolygonsDataFrame.
#' @param mcmc Object of class \code{link[coda]{mcmc.list}}.
#" @param param Parameter kin \code{mcmc} used to encode formatting for points/polygons.
#' @param col List of two or more colors to be used as color scale for plotting points. Values will be stretched from the lowest to highest value if \code{stretch = TRUE} (default).
#' @param stretch Logical, if \code{TRUE} (default), then colors will be stretched from the minimum to maximum value.
#' @param naFields Name of field(s) in \code{x} that might have \code{NA} values. If \code{NA}s occur, then it is assumed that these polygons/points were skipped in the model construction process (JAGS cannot handle \code{NA}s), and so points/polygons that correspond to an \code{NA} will be plotted using formatting given by \code{naCol} and \code{naBorder}.
#' @param naCol Fill color to be used for points/polygons that have \code{NA} values.
#' @param naBorder Border color to be used for points/polygons that have \code{NA} values. (Ignored for points unless \code{pch} is 21 through 25 in which case it is used as the point border color.
#' @param ... Arguments to pass to \code{plot}
#' @return Nothing (displays a plot).
#' @seealso
#' @examples
#'	set.seed(123)
#' @export

plotBayesShape <- function(
	x,
	mcmc,
	param,
	col = c('white', 'black'),
	stretch = TRUE,
	naFields = NULL,
	naCol = 'gray',
	naBorder = NA,
	...
) {

	args <- list(...)
	
	nas <- !complete.cases(x@data[ , naField])

}
