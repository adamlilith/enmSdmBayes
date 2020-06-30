#******************************************
#
#    enmSdmBayes - R Shiny UI
#    by Helgasoft for www.earthSkySea.org
#
#******************************************

# ------- UI -----

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    tags$style(
      HTML("
      #sidebar {
            background-color: #E0FFFF;
      }

      body { font-size: 11px; line-height: 1; }
      .btnContainer { margin-top: 20px; }
      .action-button { background-color:#B2DFEE; font-size: 12px;  padding: 3px 10px; }
      #saveStats { background-color: #B2DFEE; font-size: 12px; padding: 3px 10px;}
      #openFile { margin-bottom: 6px; }
      .selectize-input { padding: 3px 10px; min-height: 30px; }
      .form-control { padding: 3px 10px; height: 30px;  font-size: 11px; }
      .centxt { text-align: center; font-size: 14px; font-weight: bold; margin-bottom:10px;}
      .fsmall { font-size: 10px; }
      .well { padding: 8px 12px; }
      #feedback { font-size: 13px; color: darkslateblue; }
      .f14 { font-size: 14px; margin: 40px; }
      .f12 { font-size: 12px; font-weight:bold; margin: 8px 0px; }
      .pb { margin-right:6px; margin-top:6px; }
      .centr {
        display: flex; justify-content: center; align-items: center; margin-top:6px;
      }
  "))),
  sidebarLayout(

			sidebarPanel(
				id='sidebar',
				div(class='centxt', span('enmSdmBayes'),
            br(),span('Bias-corrected distribution modeling', class='fsmall')),
        fluidRow(
          column(6,	actionButton('openFile', 'Load CSV or Shapefile')),
          column(6,	actionButton('btnHelp', 'Help', icon=icon('question-circle'),
               onclick="window.open('http://www.earthskysea.org/enmsdmbayesland/', '_blank')"), align='right')),
			  conditionalPanel(
			    condition = 'output.panelCountry',
			    selectInput('dcountry', 'Country field', '', multiple=TRUE)
			  ),
				selectInput('dstate',  'State field', ''),
				selectInput('dcounty', 'County field', ''),
				selectInput('dcovar',  'Covariate field (optional)', ''),
				selectInput('deffort', 'Effort field', ''),
				selectInput('ddetect', 'Detections field', ''),
				fluidRow(
          column(6,	numericInput('texpand', 'Model region size', 0.3, min=0.01, max=1)),
				  column(6,	actionButton('modPrep', 'Prepare'), class='btnContainer')),
				#numericInput('tchains', 'blaa', 2),
        fluidRow( column(12, div('Model settings', class='f12'))),
        fluidRow(
				  column(6,	numericInput('tburn', 'Burn-in', 1000, min=1, max=1000000)),
				  column(6,	numericInput('titer', 'Iterations', 2000, min=1, max=1000000))),
				fluidRow(
				  column(6,	numericInput('tthin', 'Thin by', 1, min=1, max=1000)),
				  column(6, actionButton('model', 'Model'), class='btnContainer')),
				fluidRow(
				  column(6,	downloadButton('saveStats', 'Statistics')),
			    column(6,	actionButton('saveFocMod', 'Shapefile', icon=icon('download'))), class='centr'),
        p(), uiOutput('org', class='centxt')
				, width=4
			),

			mainPanel(
			  fluidRow( column(12,
			  #span( class='centr',
  			  shinyjs::hidden(
  			    actionButton('pb1', 'Specimens', class='pb'),
  			    actionButton('pb2', 'Effort', class='pb'),
  	  	    actionButton('pb3', 'Region', class='pb')),
  			  #HTML('&nbsp; &nbsp; ____ &nbsp; &nbsp;'),
			  shinyjs::hidden(
  			    actionButton('pb4', 'Occupancy', class='pb'),
  			    actionButton('pb5', 'Occ Uncert.', class='pb'),
  			    actionButton('pb6', 'Detectability', class='pb'),
    			  actionButton('pb7', 'Detect. Uncert.', class='pb'),
  			    actionButton('pb8', 'stats', class='pb')
  			  ))),
			  fluidRow( column(12,
  			  # plotOutput('plot1', width='100%'), #, height='570px'),
  			  # shinyjs::hidden(
  				  # plotOutput('plot2', width='100%'), #, height='570px'),
  				  # plotOutput('plot3', width='100%'), #, height='570px'),
  				  # plotOutput('plot4', width='100%'), #, height='570px'),
  				  # plotOutput('plot5', width='100%'), #, height='570px'),
  				  # plotOutput('plot6', width='100%'), #, height='570px'),
  				  # plotOutput('plot7', width='100%'), #, height='570px'),
  				  # htmlOutput('plot8', class='f14') #, width='100%', height='570px')
  				# )
  			  plotOutput('plot1', height='900px'),
  			  shinyjs::hidden(
  				  plotOutput('plot2', height='900px'),
  				  plotOutput('plot3', height='900px'),
  				  plotOutput('plot4', height='900px'),
  				  plotOutput('plot5', height='900px'),
  				  plotOutput('plot6', height='900px'),
  				  plotOutput('plot7', height='900px'),
  				  htmlOutput('plot8', class='f14') #, width='100%', height='570px')
  				)
  				,textOutput('feedback')
			  )), width=8
			)
			,position='left'
			,fluid=FALSE
))

# ------- server ------
server <- function(input, output, session) {

  library(shiny)
  library(shinyjs)
  library(enmSdm)
  library(coda)   # for train & after
  library(nimble)
  require(raster)

# ------- globals/utils -----

  srcPath <- ''
  source(paste0(srcPath,'getGeogFocus.r'), local=TRUE)
  # works fine with full country names, but needs the .rda in subfolder 'data'
  #  source(paste0(srcPath,'getISO3.r'), local=TRUE)
  source(paste0(srcPath,'matchDfToGadm.r'), local=TRUE)
  source(paste0(srcPath,'trainBayesODM_pMaxVarByState_psiCarCovar.r'), local=TRUE)
  source(paste0(srcPath,'bayesODMStats.r'), local=TRUE)
  source(paste0(srcPath,'processBayesODM.r'), local=TRUE)
  data('countryISO3', envir=environment() )

  # init vars
  isCSVfile <- FALSE
  fcountry <- fstate <- fcounty <- fdetect <- feffort <- fcovar <- vdetect <- dcsv <- NULL
  shape <- focus <- mcmcModel <- stats <- NULL
  fname <- NULL

  html.table <- function(listo) {
    out = "<h3>Model output statistics</h3> <table style='border-spacing: 10px'>
            <tr><td><b>Name</b></td><td><b>Value</b></td></tr>"
    for(i in 1:length(listo)) out <- paste(out,'<tr><td>',names(listo[i]),'</td><td>',listo[i],'</td>')
    out <- paste(out,'</table>')
    out
  }

  slam <- function(id=NULL) {   # display one plot by name
    stopifnot(!is.null(id))
    for(i in 1:8) { if (grepl(i, id, fixed=TRUE)) shinyjs::show(id) else hide(paste0('plot',i)) }
  }

  shinyLog <- function(actor=NULL, boxID=NULL) {
    stopifnot(!is.null(boxID))
    out <- NULL
    withCallingHandlers({
      shinyjs::html(boxID, '')
      out <- actor
    },
    message = function(m) {
      shinyjs::html(id=boxID, html = paste(m$message, '<br />'), add = FALSE)
    }
    )
    out
  }   # show function messages in Status area

  disable('modPrep'); disable('model');  hide('saveStats'); hide('saveFocMod')
  disable('tburn'); disable('titer'); disable('tthin');

# ------- openFile  -----
  observeEvent(input$openFile, {
    # reset
    fcountry <<- fstate <<- fcounty <<- fdetect <<- feffort <<- fcovar <<- vdetect <<- dcsv <<- NULL
    shape <<- focus <<- mcmcModel <<- stats <<- NULL
    updateSelectInput(session, 'dcountry', label='Country field')
    updateSelectInput(session, 'dstate', label='State field')
    updateSelectInput(session, 'dcounty', label='County field')
    updateSelectInput(session, 'deffort', label='Effort field')
    updateSelectInput(session, 'ddetect', label='Detections field')
    shinyjs::html(id='feedback', html='')
    output$plot1 <- renderPlot({ plot.new() })  # erase area
    for(i in 1:8) {
      id <- paste0('pb',i); hide(id);
      id <- paste0('plot',i); if (i==1) shinyjs::show(id) else hide(id)
    }
    if (input$openFile == 'None') {   # reset dropdowns if file ddl
      updateSelectInput(session, 'dcountry', choices='',selected='' )
      updateSelectInput(session, 'dstate', choices='',selected='' )
      updateSelectInput(session, 'dcounty', choices='',selected='' )
      updateSelectInput(session, 'ddetect', choices='',selected='' )
      updateSelectInput(session, 'deffort', choices='',selected='' )
      updateSelectInput(session, 'dcovar', choices='',selected='' )
      return()
    }

    fname <<- NULL
    tryCatch({ fname <<- file.choose() }, error=function(e) { }) 	# suppress err msg on file Cancel
    if (is.null(fname)) return()

    shinyjs::html(id='feedback', html=paste('Reading file',basename(fname),'...'))

    if (endsWith(tolower(fname), '.csv')) {
      dcsv <<- read.csv(fname, header=TRUE, stringsAsFactors=FALSE)
      isCSVfile <<- TRUE

    } else {
      shape <<- try(raster::shapefile(fname, verbose=TRUE))  # calls rgDal
      if ( class(shape)=='try-error' ) {
        shinyjs::html(id='feedback', html=paste("<div style='color:red'>",shape,"</div>"))
        return()
      }
      dcsv <<- shape@data
      isCSVfile <<- FALSE
    }

    tmp <- basename(fname)
 		shinyjs::html(id='feedback', html=paste('File',ifelse(isCSVfile, tmp, strtrim(tmp, nchar(tmp)-4) ),'processed.' ))
 		disable('modPrep'); disable('model'); hide('saveStats'); hide('saveFocMod')
 		disable('tburn'); disable('titer'); disable('tthin');

    # reset dropdowns
    fchoices <- c('None',unique(names(dcsv)))

      # preset dropdowns with most widely used values for Shape files, csv will be empty
      updateSelectInput(session, 'dcountry',choices=fchoices, selected='NAME_0')
      updateSelectInput(session, 'dstate',  choices=fchoices, selected='NAME_1')
      updateSelectInput(session, 'dcounty', choices=fchoices, selected='NAME_2')
      updateSelectInput(session, 'dcovar',  choices=fchoices, selected='None')
      updateSelectInput(session, 'deffort', choices=fchoices, selected='None')
      updateSelectInput(session, 'ddetect', choices=fchoices, selected='None')

    output$panelCountry <- reactive({ isCSVfile })		# toggles country panel
    outputOptions(output, 'panelCountry', suspendWhenHidden=FALSE)
  })		# open input file (CSV or Shape)
  output$org <- renderUI({
    url <- a("www.earthSkySea.org", href="http://www.earthSkySea.org", class='fsmall')
    tagList(icon('info-circle'), url);
  })
  observeEvent(input$pb1, { slam('plot1') })
  observeEvent(input$pb2, { slam('plot2') })
  observeEvent(input$pb3, { slam('plot3') })
  observeEvent(input$pb4, { slam('plot4') })
  observeEvent(input$pb5, { slam('plot5') })
  observeEvent(input$pb6, { slam('plot6') })
  observeEvent(input$pb7, { slam('plot7') })
  observeEvent(input$pb8, { slam('plot8') })

# ------- parameters -----

      observeEvent(input$dcountry, {
			  if (input$dcountry != '' && input$dcountry != 'None') {
  			  if (is.null(fcountry)) {    # country field selected
  			    fcountry <<- input$dcountry
  			    tmp <- unique(dcsv[, fcountry])
  			    tmp <- tmp[order(tmp)]    # country names from file
  			    notFound <- setdiff(tmp, countryISO3$country) # ?not in official list
  			    if (length(notFound)==0)
  			      selc <- tmp
  			    else {
    			    found <- setdiff(tmp, notFound)     # correct names
    			    notFound <- paste0(notFound, '?')   # mark countries missing from official list
    			    selc <- c(found, notFound)
  			    }
  			    tmp <- c(selc, countryISO3$country)
  			    # initially select all file countries, user will exclude if needed
  			    updateSelectInput(session, 'dcountry', label='Country', choices=tmp, selected=selc )
  			  } else {    # actual countries selected
  			    if (any(endsWith(input$dcountry, '?'))) {
  			      alert(paste('Invalid country',input$dcountry[which(endsWith(input$dcountry,'?'))]))
  			    }
  			  }
			  }
			})
			observeEvent(input$dstate, {
			  if (input$dstate != '' && input$dstate != 'None') {
			    fstate <<- input$dstate
			  }
			})
			observeEvent(input$dcounty, {
			  if (input$dcounty != '' && input$dcounty != 'None') {
 			    fcounty <<- input$dcounty
			  }
			})
			observeEvent(input$deffort, {
			  if (input$deffort != '' && input$deffort != 'None') {
			    feffort <<- input$deffort    # effort field selected
			  }
			})
			observeEvent(input$dcovar, {
			  if (input$dcovar != '' && input$dcovar != 'None') {
			    fcovar <<- input$dcovar    # covar field selected
			  }
			})
			observeEvent(input$ddetect, {
			  if (input$ddetect != '' && input$ddetect != 'None') {
			    if (is.null(fdetect)) {    # detect field selected
			      fdetect <<- input$ddetect
			      if (isCSVfile) {
			        tmp <- unique(na.omit(dcsv[, fdetect]))    # CSV only
			        updateSelectInput(session, 'ddetect', label='Detections', choices=unique(tmp[order(tmp)]) )
			      }
			      enable('modPrep')
			    } else {    # actual detection selected
			      if (isCSVfile)
			        vdetect <<- input$ddetect
			    }
			  }
			})

# ------- modPrep ----

			observeEvent(input$modPrep, {

			  # validate numeric params
			  tmp <- ''
			  if (!(input$texpand >= 0.01 & input$texpand <=1))
			    tmp <- paste(tmp,'<br/>reg.size=',input$texpand)
			  if (!(input$tburn >=1000 & input$tburn<=1000000))
			      tmp <- paste(tmp,'<br/>burn-in=',input$tburn)
			  if (!(input$titer >=1000 & input$titer<=1000000))
			    tmp <- paste(tmp,'<br/>iterations=',input$titer)
			  if (!(input$tthin >=1 & input$tthin<=1000000))
			    tmp <- paste(tmp,'<br/>thin=',input$tthin)
			  if (tmp != '') {
			      shinyjs::html(id='feedback', html=paste("<div style='color:red'>Input errors: <br />",tmp,"</div>"))
			      return()
        }
			  # hide all plot buttons
			  for(id in paste0('pb',1:8)) hide(id)
			  # clean all plots
			  for(id in paste0('plot',1:8)) output[[id]] <- renderPlot({ plot.new() })

			  shinyjs::html(id='feedback', html=paste0('Collating records and preparing modeling region for ',
			                                           ifelse(isCSVfile, vdetect, fdetect),'...'))

			  if (isCSVfile) {     # get shape for CSV
			    # user can redo Prep for diff state/fdetect, need to recover the values
			    feffort <<- input$deffort
			    fstate <<- input$dstate
			    fcounty <<- input$dcounty

			    # get spatial polygons for countries of interest
			    #iso3 <- getISO3(country[1])
			    shape <<- try( raster::getData('GADM', country=input$dcountry[1], level=2) )
			    if ( class(shape)=='try-error') {
			      shinyjs::html(id='feedback', html=paste("<div style='color:red'>cannot open URL, no Internet ?</div>"))
			      shape <<- NULL
			      return()
			    }

			    if (length(input$dcountry) > 1) {
			      for (i in 2:length(input$dcountry)) {
			        #iso3 <- getISO3(input$dcountry[i])
			        shape <<- rbind(
			          shape, raster::getData('GADM', country=input$dcountry[i], level=2)
			        )
			      }
			    }

			    # # remove Alaska and Hawaii to speed things up
			    # shape <<- shape[!(shape@data$NAME_1 %in% c('Alaska', 'Hawaii')), ]

			    # match shape and CSV (not all will match... to be expected)
			    shape <<- try(
			      matchDfToGadm(
  			      df = dcsv,
  			      stateProv = fstate,
  			      county = fcounty,
  			      species = feffort,
  			      speciesName = vdetect,
  			      shape = shape
  			    )
			    )
			    if ( class(shape)=='try-error' ) {
			      shinyjs::html(id='feedback', html=paste("<div style='color:red'>",shape,"</div>"))
			      shape <<- NULL
			      return()
			    }
			    fdetect <<- 'detect'   # all set by matchDfToGadm
			    feffort <<- 'effort'
			    fstate <<- 'NAME_1'
			    fcounty <<- 'NAME_2'

			  }     # CSV file end
			  cat('Params are s:',fstate,'c:',fcounty,'e:',feffort,'fd:',fdetect,'d:',vdetect)

			  # display shape with detections or effort
			  output$plot1 <- renderPlot({
			    scale <- shape@data[fdetect] / max(shape@data[fdetect])
			    scale = as.numeric(scale[,1])
			    cols <- scales::alpha('darkblue', scale)
			    cols[shape@data[fdetect] == 0] <- 'Gainsboro'
			    plot(shape, col=cols, border='gray', lwd=0.4, main=paste('Number of detections of', ifelse(isCSVfile, vdetect, fdetect)) )
			  })
			  slam('plot1'); shinyjs::show('pb1')

			  output$plot2 <- renderPlot({
			    scale <- shape@data[feffort] / max(shape@data[feffort])
			    scale = as.numeric(scale[,1])
			    cols <- scales::alpha('darkred', scale)
			    cols[shape@data[feffort] == 0] <- 'Gainsboro'
			    plot(shape, col=cols, border='gray', lwd=0.4, main=paste('Number of total collections (effort)'))
			  })
			  slam('plot2'); shinyjs::show('pb2')

			  # user-defined settings for selecting focal region
			  expand <- input$texpand  #0.1 # expansion factor... large values increase size of focal region (but also slow the model down)

			  #output is a subset of the SpatialPolygons
			  focus <<- try( getGeogFocus(
			    x = shape,
			    detect = fdetect,
			    expand = expand,
			    upper = fstate #'NAME_1'
			  ) )
			  if ( class(focus)=='try-error' || is.na(focus) ) {
			    shinyjs::html(id='feedback', html=paste("<div style='color:red'>",focus,"</div>"))
			    focus <<- NULL
			    return()
			  }

			  output$plot3 <- renderPlot({			    # plot focal area
			    # plot(shape, col='Gainsboro', border=NA, main='Modeling region')
			    # plot(focus, add=TRUE)
			    plot(focus, col='Gainsboro', border='cornflowerblue', main='Modeling region')
			    plot(shape, border='black', lwd=0.4, add=TRUE)
			  })
			  slam('plot3'); shinyjs::show('pb3')
			  shinyjs::html(id='feedback', html='Model preparation completed. You may implement the model now.')
			  enable('model')
			})

			observe({       # enable buttons only when input done
			  if (!any(c(
			    input$dstate=='None',
      		input$dcounty=='None',
			    #input$dcovar=='',
			    input$deffort=='None',
			    input$ddetect=='None')) &&
			    !any(c(
			      input$dstate=='',
			      input$dcounty=='',
			      input$deffort=='',
			      input$ddetect==''))
			    ) { enable('tburn'); enable('titer'); enable('tthin'); enable('model') }
			})

# ------- modeling ---------------

	observeEvent(input$model, {

	  ##############################
	  ### use area as covariate? ###
	  ##############################

	  # user chooses:
	  if (is.null(fcovar))
	    covariate <- NULL
	  else if (fcovar == '' || fcovar == 'None') {

	    useAreaAsCovariate <- TRUE

  	  if (useAreaAsCovariate) {

  	    covariate <- 'area_km2' # used later

  	    # equal-area PROJ4 string
  	    mollweide <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
  	    focusEqualArea <- sp::spTransform(focus, sp::CRS(mollweide))
  	    area_km2 <- rgeos::gArea(focusEqualArea, byid=TRUE) / 1000^2
  	    focus@data$area_km2 <- area_km2

  	  } else {

  	    covariate <- NULL

  	  }
	  } else
	    covariate <- fcovar

	  ########################
	  ### prepare to model ###
	  ########################

	  # names of columns in focus shape that represent state/province, county, effort, detections, and covariate values
	  # we know what these are because we downloaded the shapes from GADM
	  # but in the route where the user uploads a shapefile we need to ask
	  stateProv <- fstate  #'NAME_1'
	  county <- fcounty  #'NAME_2'
	  detect <- fdetect #ifelse(isCSVfile, 'detect', fdetect)
	  effort <- feffort #ifelse(isCSVfile, 'effort',  feffort)
	  # covariate <- 'area_km2' # already defined above

	  # model settings (user-defined)
	  # note that larger values of nburnin and niter will yield a better model, but can substantially increase model run time
	  niter <- input$titer #220 # for testing
	  nburnin <- input$tburn #20 # for testing
	  nchains <- 4 # for testing
	  thin <- input$tthin  #1 # for testing

	  # niter <- 1100000 # "reasonable" value... could take model several hours
	  # nburnin <- 10000 # "reasonable" value... could take model several hours
	  # nchains <- 4 # "reasonable" value... could take model several hours
	  # thin <- 10 # "reasonable" value... could take model several hours

	  # catch/warn bad values
	  if (nburnin >= niter) warning('The number of burn-in samples must be less than the number of iterations.')
	  nsamples <- (niter - nburnin) / thin
	  if (nsamples < 1000) warning('Your model will have <1000 samples so output may not be reliable.')
	  if (nburnin < 1000) warning('Your model has a burn-in of <1000 samples so output may not be reliable.')

	  #############
	  ### model ###
	  #############

	  # train model
	  # output is a list
	  # mcmcModel$mcmc: model output
	  # mcmcModel$shape: "focus" with model output appended
	  # mcmcModel$meta: model statistics
	  mcmcModel <<- shinyLog(      # to intercept msgs from train-model
	    trainBayesODM_pMaxVarByState_psiCarCovar(
  	    shape = focus,
  	    effort = effort,
  	    detect = detect,
  	    stateProv = stateProv,
  	    county = county,
  	    covariate = covariate,
  	    niter=niter,
  	    nburnin=nburnin,
  	    nchains=nchains,
  	    thin=thin,
  	    na.rm=TRUE,
  	    verbose = TRUE
  	  )
	  , 'feedback')
	  shinyjs::html(id='feedback', html='Model completed.', add=TRUE)

	  mcmc <- mcmcModel$mcmc
	  focusModel <- mcmcModel$shape

	  ######################
	  ### display output ###
	  ######################

	  ### maps of results

	  # generate centroids of counties with known detections
	  occs <- focusModel[focusModel@data[ , detect] > 0, ]
	  cents <- rgeos::gCentroid(occs, byid=TRUE)

	  # probability of occupancy (psi)
	  output$plot4 <- renderPlot({
	  	x <- focusModel@data$psi
	  	x <- x / max(x)
	    cols <- scales::alpha('darkgreen', x)
	    # plot(shape, col='Gainsboro', border=NA, main='Probability of occupancy')
	    plot(focusModel, col=cols, border='gray', lwd=0.4, main='Probability of occupancy')
	    # plot(focusModel, col=cols, border=NA, add=TRUE)
	 #   plot(focusModel, col=cols, border=NA, main='Probability of occupancy')
	    points(cents, pch=16, cex=0.4)
	  })
	  slam('plot4'); shinyjs::show('pb4')

	  # uncertainty occupancy (psi95CI)
	  output$plot5 <- renderPlot({
	    cols <- scales::alpha('darkred', focusModel@data$psi95CI)
	    plot(focusModel, col=cols, border='gray', lwd=0.4, main='Uncertainty in occupancy')
	    points(cents, pch=16, cex=0.4)
	  })
	  slam('plot5'); shinyjs::show('pb5')

	  # probability of detection assuming species is present (p)
	  output$plot6 <- renderPlot({
	  	x <- focusModel@data$p
	  	x <- x / max(x)
	  	cols <- scales::alpha('darkblue', x)
	    plot(focusModel, col=cols, border='gray', lwd=0.4, main='Probability of collection assuming presence')
	    points(cents, pch=16, cex=0.4)
	  })
	  slam('plot6'); shinyjs::show('pb6')

	  # uncertainty in detection (p95CI)
	  output$plot7 <- renderPlot({
	    cols <- scales::alpha('darkorange', focusModel@data$p95CI)
	    plot(focusModel, col=cols, border='gray', lwd=0.4, main='Uncertainty in probability of collection assuming presence')
	    points(cents, pch=16, cex=0.4)
	  })
	  slam('plot7'); shinyjs::show('pb7')

	  ########################
	  ### model statistics ###
	  ########################

	  stats <<- bayesODMStats(mcmc, focusModel)
	  str(stats)

	  tmp=list();i=1; apply(stats, 2, function(z) { tmp[[names(stats)[i]]] <<- ifelse(is.integer(z), z, round(z,4)); i<<-i+1 })
	  #shinyjs::html(id='feedback', html=paste("<div>",html.table(tmp),"</div>"))
	  output$plot8 <- renderUI({ HTML(html.table(tmp)) })
	  slam('plot8'); shinyjs::show('pb8')

	  shinyjs::show('saveStats'); shinyjs::show('saveFocMod')
	  # prepare fname for shape Save
	  fname <<- strtrim(basename(fname), nchar(basename(fname))-ifelse(grepl('.', basename(fname), fixed=TRUE), 4, 0))
	})

	output$saveStats <- downloadHandler(filename='stats.txt', content=function(fileloc) {
	  if (is.null(stats)) return()
    tmp <- data.frame(name=rownames(t(stats)), value=t(stats)[,1])
    rownames(tmp) <- c()
    #write(tmp, file, ncolumns=2, sep='\t\t')
    write.table(tmp, fileloc, row.names=FALSE, quote=FALSE, sep='\t')
	})
	observeEvent(input$saveFocMod, {
	  if (is.null(mcmcModel) || is.null(fname)) return()
	  focusModel <- mcmcModel$shape
	  fileloc <- paste0(path=paste0(srcPath,'data/'), fname, '_focusModel')
	  raster::shapefile(focusModel, filename=fileloc, overwrite=TRUE)
	  shinyjs::html(id='feedback', html=paste("<div style='color:green'>Saved focusModel to ",fileloc,"</div>"))
	})

}

# ------- exec.command ------
shinyApp(ui, server)
