####################################################################################
# helper.R
#
# Various helper functions
# 01jun2017
# Philip Barrett, Washington DC
#
####################################################################################

rg.read <- function( cty = 'USA', start.date = "1960-01-01" ){
## Read and cleans in the data for the specified country
  rfr <- read.csv('data/riskfreerates.csv')
  gth <- read.csv('data/growthrates.csv')
      # Read the data
  cty.gth <- data.frame( date=as.Date(gth$DATE), gth=gth[[cty]] )
  cty.rfr <- data.frame( date=as.Date(rfr$DATE), rfr = 100 * ( (1+rfr[[cty]]/100) ^ .25 - 1 ) )
      # Create country-specific dataframes
  cty.dta <- merge( subset( cty.gth, date >= start.date ),
                    subset( cty.rfr, date >= start.date ) )
      # The country data after the start date
  cty.dta$rmg <- apply( cty.dta[,-1], 1, diff )
      # Create R minus G
  return( cty.dta )
}

rg10.read <- function( cty = 'USA', start.date = "1960-01-01" ){
## Read and cleans in the ten-year data for the specified country
  rfr <- read.csv('data/tenyrrates.csv')
  gth <- read.csv('data/growthrates.csv')
      # Read the data
  cty.rfr <- data.frame( date=as.Date(rfr$DATE), rfr = 100 * ( (1+rfr[[cty]]/100) ^ .25 - 1 ) )
  cty.gth <- data.frame( date=as.Date(gth$DATE), gth=gth[[cty]] )
  # cty.gth <- data.frame( date=as.Date(gth$DATE)[1:(nrow(gth)-40)],
  #                        gth=(exp(filter( log(1+gth[[cty]]/100), rep(1/40,40), sides=1 )[-(1:40)]) - 1) * 100 )
      # Create country-specific dataframes
  cty.dta <- merge( subset( cty.gth, date > start.date ),
                    subset( cty.rfr, date > start.date ) )
      # The country data after the start date
  cty.dta$rmg <- apply( cty.dta[,-1], 1, diff )
      # Create R minus G
  return( cty.dta )
}

hist.read <- function( cty = 'USA', start.year = 1880, ltr=FALSE ){
# Reads up the historical data and returns relevant series
  cty.ifs <- switch( cty,
                 'USA'=111, 'UK'=112, 'GBR'=112, 'FRA'=132,
                 'DEU'=134, 'CAN'=156, 'JPN'=158, 'ITA'=136 )
      # IFS code dictionary
  mauro <- read.csv("data/mauro.csv")
  names(mauro)[1] <- 'ifs'
  mauro.cty <- subset(mauro, ifs==cty.ifs)
      # The Mauro database
  jst.gov <- read.csv('data/JSTgovernmentR2.csv')
  jst.real <- read.csv('data/JSTrealR2.csv')
  jst.mon <- read.csv('data/JSTmoneyR2.csv')
  jst <- merge( merge( jst.gov, jst.real, by=c('ifs','year','country','iso') ),
                jst.mon, by=c('ifs','year','country','iso') )
      # Merge the JST data
  jst.cty <- subset( jst, ifs==cty.ifs )
      # Country-specific subset
  out <- merge( mauro.cty, jst.cty, by=c('ifs','year'), all=TRUE )
  out$gth <- c( NA, ( out$gdp[-1] / out$gdp[-nrow(out)] - 1 ) * 100 )
  out$rfr <- if(ltr) out$ltrate else out$stir
  out$date <- out$year
    # Standard format
  return(subset(out,year>=start.year))
}





