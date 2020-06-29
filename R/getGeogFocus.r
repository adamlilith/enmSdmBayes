#' Select a geographic region around the distribution of a species
#'
#' This function takes as arguments a SpatialPolygonsDataFrame with a column representing whether or not a species was detected in each polygon. It then selects a subset of the shape focused around the area of detections.
#' @param x SpatialPolygonsDataFrame
#' @param detect Character, name of column in \code{x} that contains detection data.  Zero values are assumed to represent no detections.
#' @param expand Positive numeric, determines size of the focal region. For example, 0.5 corresponds to a region that comprises the occupied polygons plus a buffer with a width about 50% as large as the distance from the centroid of occupied polygons to the farthest occupied polygon.
#' @param upper Character, name of column in \code{x} that contains names of "upper" level polygons that contain "lower" level polygons which are at the scale at which detection/non-detection is scored. This argument is used to remove all lower-level polygons that fall in the same upper-level polygon but have no detections among them. For example, if upper-level polygons are states/provinces and lower-level polygons are counties with detections, then this could be used to remove all states/provinces with no counties with detections. Note that \code{minLowerPolysPerUpper} must be specified. A value of \code{NULL} (default) will not remove any polygons regardless of the values of \code{minLowerPolysPerUpper}..
#' @param minLowerPolysPerUpper Character, only used if \code{upper} is not \code{NULL}. Minimum number of lower-level polygons necessary to be present in an upper-level polygon if none of them have any detections.
#' @return SpatialPolygonsDataFrame.
#' @examples
#' \donttest{
#' 
#' library(sp)
#' 
#' # polygons of Mexican counties
#' x <- raster::getData('GADM', country='MEX', level=2)
#' 
#' # generate detections
#' set.seed(123)
#' x@data$detect <- 0
#' n <- sum(x@data$NAME_1 == 'Aguascalientes')
#' x@data$detect[x@data$NAME_1 == 'Aguascalientes'] <- rpois(n, 7)
#' 
#' n <- sum(x@data$NAME_1 == 'Jalisco')
#' x@data$detect[x@data$NAME_1 == 'Jalisco'] <- rpois(n, 4)
#' 
#' n <- sum(x@data$NAME_1 == 'Aguascalientes')
#' x@data$detect[x@data$NAME_1 == 'Aguascalientes'] <- rpois(n, 1)
#' 
#' cols <- x@data$detect / max(x@data$detect)
#' cols <- scales::alpha('darkred', cols)
#' plot(x, col=cols, border='gray90')
#' 
#' # generate focal region around detections
#' focus1 <- getGeogFocus(x, detect='detect', expand=0.1)
#' focus2 <- getGeogFocus(x, detect='detect', expand=0.5)
#' focus3 <- getGeogFocus(x, detect='detect', expand=0.1, upper='NAME_1')
#' focus4 <- getGeogFocus(x, detect='detect', expand=0.5, upper='NAME_1')
#' 
#' # plot
#' par(mfrow=c(2, 2))
#' 
#' cols <- focus1@data$detect / max(focus1@data$detect)
#' cols <- scales::alpha('darkred', cols)
#' plot(focus1, col=cols, border='gray90', main='expand=0.1')
#' 
#' cols <- focus2@data$detect / max(focus2@data$detect)
#' cols <- scales::alpha('darkred', cols)
#' plot(focus2, col=cols, border='gray90', main='expand=0.5')
#' 
#' cols <- focus3@data$detect / max(focus3@data$detect)
#' cols <- scales::alpha('darkred', cols)
#' plot(focus3, col=cols, border='gray90', main='expand=0.1 | sans sparse')
#' 
#' cols <- focus4@data$detect / max(focus4@data$detect)
#' cols <- scales::alpha('darkred', cols)
#' plot(focus4, col=cols, border='gray90', main='expand=0.5 | sans sparse')
#' 
#' 
#' }
#' @export

getGeogFocus <- function(
	x,
	detect,
	expand = 0.3,
	upper = NULL,
	minLowerPolysPerUpper = 10
) {

	mollweide <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

	origProj <- raster::projection(x)
	xSpEa <- sp::spTransform(x, sp::CRS(mollweide))

	# get MCP of polygon representing counties with detections
	speciesIndex <- which(xSpEa@data[ , detect] > 0)
	occs2SpEa <- xSpEa[speciesIndex, ]
	mcpSpEa <- enmSdm::mcpFromPolygons(occs2SpEa)
	
	# study region includes MCP and all counties with records
	occsDissolveSpEa <- rgeos::gUnaryUnion(occs2SpEa)
	occsPlusMcpSpEa <- rgeos::gUnion(occsDissolveSpEa, mcpSpEa)

	centSpEa <- rgeos::gCentroid(occsPlusMcpSpEa)
	distFromCent_m <- rgeos::gDistance(occsPlusMcpSpEa, centSpEa, hausdorff=TRUE)

	# study region includes buffer region around this area
	focalAreaSpEa <- rgeos::gBuffer(occsPlusMcpSpEa, width=expand * distFromCent_m)
	
	# get counties in this area
	xSpEa$id <- 1:nrow(xSpEa)
	
	focalAreaSp <- sp::spTransform(focalAreaSpEa, enmSdm::getCRS('wgs84'))
	focalIndices <- sp::over(focalAreaSpEa, xSpEa, returnList=TRUE)
	
	focalIndices <- sp::over(focalAreaSpEa, xSpEa, returnList=TRUE)
	focusSpEa <- x[focalIndices[[1]]$id, ]
	
	# removing states with with no detections that have too-few counties
	if (!is.null(upper)) {
		uppers <- unique(focusSpEa@data[ , upper])
		for (thisUpper in uppers) {
			focusUpper <- focusSpEa[focusSpEa@data[ , upper] == thisUpper, ]
			if (all(focusUpper@data[ , detect] == 0) & nrow(focusUpper) < minLowerPolysPerUpper) {
				focusSpEa <- focusSpEa[-which(focusSpEa@data[ , upper] == thisUpper), ]
			}
		}
	}

	focusSp <- sp::spTransform(focusSpEa, sp::CRS(origProj))
	focusSp

}
