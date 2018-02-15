#####################################################################
# cdsRating.R
#
# Code to compute implied default probabilities as a funciton of ratings.
# 28jun2017
# Philip Barrett, Washington DC
#####################################################################

# #### Get EU data ####
# library(eurostat)
# id <- search_eurostat("General government gross debt")
# dat <- get_eurostat(id$code)
#
# v.cty <- as.character(unique(dat$geo))
# n.cty <- length(v.cty)
# par(mfrow=rep(ceiling(sqrt(n.cty)),2), mar=c(2,2,1,1))
# for( cty in v.cty){
#   with( subset( dat, indic=='GD' & geo == cty & unit == 'PC_GDP' ),
#       plot(time, values, type='l'), lwd=2)
#   title(main=cty, line=-1 )
# }
# par(mfrow=c(1,1))
#

### Read in CDS spreads ###
library(zoo)
cds <- read.csv('data/cds5yr.csv', na.strings = '#N/A', stringsAsFactors = FALSE )[-1,]
cds$q.date <- as.yearqtr(cds$date, format = "%m/%d/%Y")
cds$date <- as.yearmon(cds$date, format = "%m/%d/%Y")
cds[,-c(1,ncol(cds))] <- apply(cds[,-c(1,ncol(cds))], 2, as.numeric)
cds.q <- aggregate( cds[,-ncol(cds)], list(cds$q.date), mean, na.rm=TRUE )
    # Quarterly data
with( cds, plot(date, United.States, type='l', lwd=2) )
with( cds.q, lines(date, United.States, type='l', lwd=2, col='red') )
