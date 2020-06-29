### create example data for enmSdmBayes
### Adam B. Smith | 2020-06

library(BIEN)
library(enmSdm)
library(raster)
library(rgeos)
library(sp)

	### shapefile
	#############

		# get GADM
		can <- raster::getData('GADM', country='CAN', level=2, path='C:/ecology/!Scratch')
		usa <- raster::getData('GADM', country='USA', level=2, path='C:/ecology/!Scratch')
		mex <- raster::getData('GADM', country='MEX', level=2, path='C:/ecology/!Scratch')

		northAmerica <- rbind(can, usa, mex)
		northAmerica <- northAmerica[!(northAmerica@data$NAME_1 %in% c('Alaska', 'Hawaii')), ]

		northAmericaEa <- sp::spTransform(northAmerica, enmSdm::getCRS('albersNA', TRUE))
		centsEa <- rgeos::gCentroid(northAmericaEa, byid=TRUE)
		cents <- sp::spTransform(centsEa, enmSdm::getCRS('wgs84', TRUE))
		centCoords <- sp::coordinates(cents)
		
		long <- -105
		lat <- 54
		
		inWest <- which(centCoords[ , 1] < long & centCoords[ , 2] < lat)
		west <- northAmerica[inWest, ]
		west@data <- west@data[ , c('NAME_0', 'NAME_1', 'NAME_2')]
		names(west@data) <- c('countr', 'state', 'county')
		projection(west) <- enmSdm::getCRS('wgs84')
		
		plot(west)
		
		# get occurrences
		occs <- BIEN::BIEN_occurrence_family(
			family = 'Apocynaceae',
			cultivated = FALSE,
			only.new.world = TRUE,
			observation.type = FALSE,
			all.taxonomy = FALSE,
			native.status = FALSE,
			natives.only = TRUE,
			political.boundaries = FALSE,
			collection.info = FALSE
		)

		occsWest <- occs[which(occs$longitude < long & occs$latitude < lat), ]
		ext <- terra::extract(west, occsWest[ , c('longitude', 'latitude')])
		occsWest <- cbind(occsWest, ext)
		
		# match effort to shapefile
		west$apocyn <- 0
		for (i in 1:nrow(west)) {
		
			these <- which(west$NAME_1[i] == occsWest$NAME_1 & west$NAME_2[i] == occsWest$NAME_2)
			if (length(these) > 0) {
				west$apocyn[i] <- length(these)
			}
		}

		# match species 1 to shapefile
		west$ascAlb <- 0
		occsSpp <- occsWest[which(occsWest$scrubbed_species_binomial == 'Asclepias albicans'), ]
		
		for (i in 1:nrow(west)) {
		
			these <- which(west$NAME_1[i] == occsSpp$NAME_1 & west$NAME_2[i] == occsSpp$NAME_2)
			if (length(these) > 0) {
				west$ascAlb[i] <- length(these)
			}
		}

		# match species 2 to shapefile
		west$ascCal <- 0
		occsSpp <- occsWest[which(occsWest$scrubbed_species_binomial == 'Asclepias californica'), ]
		
		for (i in 1:nrow(west)) {
		
			these <- which(west$NAME_1[i] == occsSpp$NAME_1 & west$NAME_2[i] == occsSpp$NAME_2)
			if (length(these) > 0) {
				west$ascCal[i] <- length(these)
			}
		}

		# match species 3 to shapefile
		west$apoCan <- 0
		occsSpp <- occsWest[which(occsWest$scrubbed_species_binomial == 'Apocynum cannabinum'), ]
		
		for (i in 1:nrow(west)) {
		
			these <- which(west$NAME_1[i] == occsSpp$NAME_1 & west$NAME_2[i] == occsSpp$NAME_2)
			if (length(these) > 0) {
				west$apoCan[i] <- length(these)
			}
		}

		westEa <- sp::spTransform(west, enmSdm::getCRS('albersNA', TRUE))
		west@data$areaKm2 <- rgeos::gArea(westEa, byid=TRUE) / 1000^2
		
		shapefile(west, 'C:/Ecology/Drive/R/enmSdmBayes/R/data/westShape', overwrite=TRUE)
		
	### CSV
	#######

		occsWest <- occsWest[ , c('scrubbed_family', 'scrubbed_species_binomial', 'date_collected', 'NAME_0', 'NAME_1', 'NAME_2')]
		names(occsWest) <- c('family', 'binomial', 'collectDate', 'country', 'stateProv', 'county')
		write.csv(occsWest, 'C:/Ecology/Drive/R/enmSdmBayes/R/data/westCsv.csv', row.names=FALSE)
		