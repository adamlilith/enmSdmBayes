#' Try to return an ISO3 country code from a potentially "messy" country name
#'
#' This function takes as an argument a potentially "messy" set of country names and attempts to return ISO3 codes for each. For example, "messy" names might be a vector \code{c('USA', 'United States', 'U.S.A.', 'U. S. A', 'N. Zealand', 'New Zealand')}.  It will also try to do fuzzy matching if desired.
#' @param x Vector of country names.
#' @param crosswalk Eitehr \code{NULL} (default) or a two-column data frame or matrix that serves as a lookup table for ISO3 codes. If a data frame or matrix is supplied, the first column must be country names and the second must be ISO3 codes.
#' @return Character vector.
#' @examples
#' x <- c('USA', 'United States', 'U.S.A.', 'U. S. A', 'N. Zealand',
#' 'New Zealand', 'North Korea', 'N. Korea', 'Marshall Islands',
#' 'Marshall Is', 'Saint Helena', 'St. Helena', 'Mars')
#' getISO3(x)
#' @export

getISO3 <- function(x, crosswalk = NULL) {

	if (is.null(crosswalk)) {
		warning('TO DO: Uncomment line below to load ISO3 table. Assuming it is already loaded.')
		# data('countryISO3')
		crosswalk <- countryISO3
	}
	
	crosswalk[ , 1] <- .cleanString(crosswalk[ , 1])
	x <- .cleanString(x)
	
	isos <- crosswalk[match(x, crosswalk[ , 1]), 2]
	isos

}

# clean strings
.cleanString <- function(s) {

	s <- tolower(s)
	s <- gsub(s, pattern='st. ', replacement='saint')
	s <- gsub(s, pattern='ste. ', replacement='saint')
	s <- gsub(s, pattern='[.]', replacement='')
	s <- gsub(s, pattern='&', replacement='')
	s <- gsub(s, pattern='\'', replacement='')
	s <- gsub(s, pattern=',', replacement='')
	s <- gsub(s, pattern=' and ', replacement='')
	s <- gsub(s, pattern=' the ', replacement='')
	s <- gsub(s, pattern='[(]', replacement='')
	s <- gsub(s, pattern='[)]', replacement='')
	s <- gsub(s, pattern='[[]', replacement='')
	s <- gsub(s, pattern='[]]', replacement='')
	s <- gsub(s, pattern='[{]', replacement='')
	s <- gsub(s, pattern='[}]', replacement='')
	s <- gsub(s, pattern='-', replacement='')
	s <- gsub(s, pattern=';', replacement='')
	s <- gsub(s, pattern='/', replacement='')
	s <- gsub(s, pattern=' ', replacement='')
	s

}
