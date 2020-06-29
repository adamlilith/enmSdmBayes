#' Collate number of detections and effort into a GADM SpatialPolygonsDataFrame
#'
#' This function merges information on detections and effort from a data frame into a SpatialPolygonsDataFrame object from GADM.
#' @param df Data frame
#' @param stateProv Name of field in \code{x} with stateProv/province names.
#' @param county Name of field in \code{x} with county names.
#' @param species Name of field in \code{x} with species names.
#' @param speciesName Name of focal species.
#' @param shape SpatialPolygonsDataFrame in GADM level 2 format (www.gadm.org).
#' @examples

matchDfToGadm <- function(
	df,
	stateProv,
	county,
	species,
	speciesName,
	shape
) {

	if (!(speciesName %in% df[ , species])) stop('The species does not appear in this data frame.')
	
	shape@data$detect <- 0
	shape@data$effort <- 0
	
	unfoundAll <- 0 # number of all records that cannot be matched
	unfoundSpecies <- 0 # number of focal species' records that cannot be matched
	
	# replace characters with diacritics... makes matching easier
	shape@data$NAME_1 <- stringi::stri_trans_general(shape@data$NAME_1, 'latin-ascii')
	shape@data$NAME_2 <- stringi::stri_trans_general(shape@data$NAME_2, 'latin-ascii')
	
	df[ , stateProv] <- stringi::stri_trans_general(df[ , stateProv], 'latin-ascii')
	df[ , county] <- stringi::stri_trans_general(df[ , county], 'latin-ascii')
	
	for (i in 1:nrow(df)) {
	
		# match data frame state/province and county to shape's state/province and county
		shapeIndex <- which(shape@data$NAME_1 == df[i, stateProv] & shape@data$NAME_2 == df[i, county])
		
		# no geographic match
		if (length(shapeIndex) == 0) {
		
			unfoundAll <- unfoundAll + 1
			if (!is.na(df[i, species]) && df[i, species] == speciesName) unfoundSpecies <- unfoundSpecies + 1
	
		# geographic match and is the focal species
		} else  {
		
			shape@data$effort[shapeIndex] <- shape@data$effort[shapeIndex] + 1

			if (!is.na(df[i, species]) && df[i, species] == speciesName) {
				shape@data$detect[shapeIndex] <- shape@data$detect[shapeIndex] + 1
			}
			
		}
		
	}
	
	if (unfoundAll > 0 | unfoundSpecies > 0) warning(unfoundAll, ' records were not matched, of which ', unfoundSpecies, ' are of ', speciesName, '.')

	shape

}
