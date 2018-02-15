library(plyr)
library(lubridate)
library(ggplot2)
library(ggthemes)

cty <-  'USA' # 'GBR' #   'FRA' # 'CAN' #'ITA' # 'JPN' # 'DEU' #
G7 <- c( 'USA', 'GBR', 'FRA', 'CAN', 'DEU', 'ITA', 'JPN' ) #, 'ESP', 'NLD', 'PRT', 'SWE' )
plot.save <- TRUE # FALSE # TRUE #
rg.dta <- hist.read(cty,start.year)
jst.real <- read.csv('data/JSTrealR2.csv')
jst.real <- ddply( jst.real, 'country', transform, gth=c( NA, ( exp(diff(log(gdp))) - 1 ) * 100 ))
    # Add growth
jst.mon <- read.csv('data/JSTmoneyR2.csv')
    # For other interest rates
boe.rate <- read.csv('data/BOEPRUKA.csv')
boe.rate$year <- year(boe.rate$DATE)
boe.rate <- boe.rate[,3:2]
names(boe.rate)[2] <- 'rfr'
    # Bank rate
fed.rate <- read.csv('data/FEDFUNDS.csv')
fed.rate$year <- year(fed.rate$DATE)
fed.rate <- fed.rate[,3:2]
names(fed.rate)[2] <- 'rfr'
    # Fed dunds rate
dta <- merge( subset(jst.real, iso==cty)[,c('year','gth')],
              subset(jst.mon, iso==cty)[,c('year','stir','ltrate')], by='year')
if(cty=='GBR') dta <- merge( dta, boe.rate, by='year')
if(cty=='USA') dta <- merge( dta, fed.rate, by='year', all = TRUE)
dta <- subset(dta, year >= 1880)
    # The chart data
dta$rmg.1 <- dta$stir - dta$gth
dta$rmg.2 <- dta$ltrate - dta$gth
if(cty=='GBR' | cty =='USA') dta$rmg.3 <- dta$rfr - dta$gth
    # The interest-growth differentials
par.dft <- par()

if(plot.save) pdf(paste0('~/Dropbox/2017/research/debtLimits/charts/rmg_ts_', cty, '.pdf'), width = 10, height=5.5)
par( mar=c(3,4,2,2))
if( cty !='DEU') with( dta, plot( year, rmg.1, type='l', lwd=2, xlab='', ylab='Interest-growth differential' ) )
if( cty =='DEU') with( dta, plot( year, rmg.1, type='l', lwd=2, xlab='', ylab='Interest-growth differential',
                       ylim=c(-60,30) ) )
abline(h=0)
with( dta, lines( year, rmg.2, type='l', col='blue', lty=2, lwd=2 ) )
if(cty=='GBR' ) with( dta, lines( year, rmg.3, type='l', col='red', lty=3, lwd=2 ) )
legend.text <- c('Short-term rate', 'Long-term rate')
if(cty=='GBR') legend.text <- c( legend.text, 'Bank Rate')
# if(cty=='USA') legend.text <- c( legend.text, 'Fed Funds Rate')
legend('topright', legend.text, col=c('black','blue','red'), lty=1:3, lwd=2, bty='n')
if(plot.save) dev.off()

par(par.dft)

sapply( subset( dta, year > 1950), mean, na.rm=TRUE )
sapply( dta, mean, na.rm=TRUE )
quantile( diff(dta$rmg.1), c(.25, .5, .75, .9), na.rm = TRUE)
quantile( diff(subset(dta,year>1950)$rmg.1), c(.25, .5, .75, .8, .9), na.rm = TRUE)

dta.int <- merge( jst.mon[,c('year', 'country', 'iso', 'stir', 'ltrate')], jst.real[, c('year', 'country','gth')], by=c('country','year') )
dta.int <- subset( dta.int, year >= 1880 )
decade.bks <- c( seq( 1879, 2010, by=5 ), 2017 )
decade.labs <- sapply(2:length(decade.bks), function(i) paste( decade.bks[i-(1:0)]+c(1,0), collapse='-' ) )
dta.int$decade <- cut( dta.int$year, breaks = decade.bks, labels = decade.labs )
mean.decade <- ddply( dta.int, .(country, decade, iso), summarize, rmg=mean(stir-gth) )
mean.decade$rmg[mean.decade$rmg < -20] <- NA
names(mean.decade)[1] <- 'Country'
g.1 <-
  ggplot( subset(mean.decade, iso %in% G7), aes( x=decade, y=rmg, group=Country, colour=Country,
                                                 linetype=Country, shape=Country) ) +
  geom_line() + geom_point() +
  geom_hline(yintercept=0) + theme_minimal() +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1.1) )
print(g.1)
ggsave('~/Dropbox/2017/research/debtLimits/charts/rmg_intl.pdf', g.1 )

dta.int$rmg <- dta.int$stir - dta.int$gth
lr.mean <- with( subset(dta.int, iso %in% G7), mean( rmg, na.rm=TRUE, trim=.05) )
mean.60s <- with( subset(dta.int, iso %in% G7 & year %in% (1960 + 0:9) ), mean( rmg, na.rm=TRUE, trim=.05) )
mean.70s <- with( subset(dta.int, iso %in% G7 & year %in% (1970 + 0:9) ), mean( rmg, na.rm=TRUE, trim=.05) )
mean.80s <- with( subset(dta.int, iso %in% G7 & year %in% (1980 + 0:9) ), mean( rmg, na.rm=TRUE, trim=.05) )
mean.90s <- with( subset(dta.int, iso %in% G7 & year %in% (1990 + 0:9) ), mean( rmg, na.rm=TRUE, trim=.05) )
dta.int.wide <- dcast( dta.int, year ~ iso, value.var = 'rmg' )
write.csv(dta.int.wide, 'data/dta_int.csv')


