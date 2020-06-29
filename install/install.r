### install file for enmSdmBayes
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2020
###
### last updated 2020-03-31

	# load/install packages
	packs <- c('BIEN', 'coda', 'dismo', 'grpreg', 'grpregOverlap', 'nimble', 'spdep', 'wiqid', 'stringi', 'sp', 'tidyverse', 'raster', 'rgeos', 'scales', 'shiny', 'shinyjs')
	
	for (pack in packs) {
	
		worked <- do.call(require, args=list(package=pack))
		if (!worked) {
			install.packages(pack, repos='https://cloud.r-project.org')
			do.call(require, args=list(package=pack))
		}
	
	}
	
	pack <- c('omnibus', 'legendary', 'birdsEye', 'enmSdm', 'statisfactory')
	
	for (pack in packs) {
	
		worked <- do.call(require, args=list(package=pack))
		if (!worked) {
			remotes::install_github(paste0('adamlilith/', pack))
			do.call(require, args=list(package=pack))
		}
	
	}
