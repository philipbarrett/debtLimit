#####################################################################
# RminusG.R
#
# Session code to create charts for R minus G and similar
# 04apr2017
# Philip Barrett, Washington DC
#####################################################################

library(ggplot2)
library(reshape2)
library(lubridate)
library(scales)
library(xtable)

rfr <- read.csv('data/riskfreerates.csv')
gth <- read.csv('data/growthrates.csv')
rfr <- rfr[,!(names(rfr)=='ESP')]
gth <- gth[,!(names(gth)=='ESP')]
rfr$DATE <- as.Date(rfr$DATE)
gth$DATE <- as.Date(gth$DATE)
    # Read and clean the data
gth[,-1] <- 100*( ( 1 + gth[,-1]/100 ) ^ 4 - 1 )
    # Annualize the quarterly growth rate

rfr.melt <- melt(rfr,id="DATE", variable.name='Country', value.name='Overnight rate')
ggplot(rfr.melt,aes(x=DATE,y=`Overnight rate`,colour=Country,group=Country)) + geom_line() +
  scale_x_date(labels=date_format("%Y"))

gth.melt <- melt(gth,id="DATE", variable.name='Country', value.name='Growth rate')
ggplot(gth.melt,aes(x=DATE,y=`Growth rate`,colour=Country,group=Country)) + geom_line() +
  scale_x_date(labels=date_format("%Y"))

gth.4 <- data.frame( DATE=gth$DATE, apply( gth[,-1], 2, function(x) stats::filter(x,c(1,1,1,1), sides=1)/4 ) )
gth.4.melt <- melt(gth.4,id="DATE", variable.name='Country', value.name='Annual growth rate')
ggplot(gth.4.melt,aes(x=DATE,y=`Annual growth rate`,colour=Country,group=Country)) + geom_line() +
  scale_x_date(labels=date_format("%Y"))

# usa.r.g <- data.frame( DATE=gth.4$DATE, 'Annual growth rate'=gth.4$USA, 'Riskfree rate'=rfr$USA )
# usa.r.g.melt <- melt(usa.r.g,id="DATE", variable.name='Series', value.name='Value')
# ggplot(usa.r.g.melt,aes(x=DATE,y=`Value`,colour=Series,group=Series)) + geom_line() +
#   scale_x_date(labels=date_format("%Y"))

rminusg <- rfr - gth
rminusg$DATE <- as.Date(rfr$DATE)
    # R minus G
rminusg.melt <- melt(rminusg,id="DATE", variable.name='Country', value.name='R minus G')
ggplot(rminusg.melt,aes(x=DATE,y=`R minus G`,colour=Country,group=Country)) + geom_line() +
  scale_x_date(labels=date_format("%Y"))

cuts <-as.Date(c( '1950/01/01', '1960/01/01', '1970/01/01', '1980/01/01',
                  '1990/01/01', '2000/01/01', '2009/01/01', '2018/01/01'))
labs <- c('50s', '60s', '70s', '80s', '90s', '2000s pre-09', 'post-09')

decades.mu <- aggregate( rminusg[,-1], list(cut(rminusg$DATE, cuts, labs )), mean, na.rm=TRUE )
names(decades.mu)[1] <- 'Decade'
decades.mu.melt <- melt(decades.mu,id="Decade", variable.name='Country', value.name='R minus G')
g1 <- ggplot(decades.mu.melt,aes(x=Decade,y=`R minus G`,colour=Country,group=Country)) +
          geom_line(size=1) + ggtitle('Decadal means')

rminusg.g7 <- rminusg[,-ncol(rminusg)]
    # Drop Sweden
rminusg.g7[,-1] <- rminusg.g7[,-1] #/ 4
    # De-annualize. NOPE
rminusg.g7 <- subset( rminusg.g7, DATE >= 01-01-1960 )
    # Drop the 1950a
print( mean(unlist(rminusg.g7[,-1]), na.rm=T ) )

decades.mu.g7 <- aggregate( rminusg.g7[,-1], list(cut(rminusg.g7$DATE, cuts, labs )), mean, na.rm=TRUE )
names(decades.mu.g7)[1] <- 'Decade'
decades.mu.g7.melt <- melt(decades.mu.g7,id="Decade", variable.name='Country',
                              value.name='Interest growth differential')
g1.g7 <- ggplot(decades.mu.g7.melt,aes(x=Decade,y=`Interest growth differential`,color=Country,group=Country)) +
  geom_line(size=1) + theme_classic() + theme(axis.title.y=element_blank(),
                                              axis.title.x=element_blank() ) # + ggtitle('Decadal means')
ggsave('~/Dropbox/2017/research/debtLimits/charts/rmg.pdf', g1.g7)

sd.out <- function( x, sd.mult = 3, ... ){
# Computes standard deviation without outliers
  sd( x[ abs( x - mean(x, na.rm = TRUE) ) < sd.mult * sd(x, na.rm = TRUE) ], na.rm = TRUE )
}

decades.sd.g7 <- aggregate( rminusg.g7[,-1], list(cut(rminusg.g7$DATE, cuts, labs )), sd, na.rm=TRUE )
decades.sd.out.g7 <- aggregate( rminusg.g7[,-1], list(cut(rminusg.g7$DATE, cuts, labs )), sd.out, na.rm=TRUE )
names(decades.sd.out.g7)[1] <- ''
decades.sd.out.g7$Average <- apply( decades.sd.out.g7[,-1], 1, mean, na.rm=TRUE )
print( xtable(decades.sd.out.g7, digits = 1, label = 'tab:rmg_sd',
              caption = 'Decadal standard deviation of interest-growth differentials in six advanced economies 1960-2016' ),
       file = '~/Dropbox/2017/research/debtLimits/charts/rmg_sd.tex', include.rownames=FALSE )


decades.sd <- aggregate( rminusg[,-1], list(cut(rminusg$DATE, cuts, labs )), sd, na.rm=TRUE)
names(decades.sd)[1] <- 'Decade'
decades.sd.melt <- melt(decades.sd,id="Decade", variable.name='Country', value.name='R minus G sd')
ggplot(decades.sd.melt,aes(x=Decade,y=`R minus G sd`,colour=Country,group=Country)) + geom_line(size=1)

decades.ad <- aggregate( rminusg[,-1], list(cut(rminusg$DATE, cuts, labs )), function(x) mean(abs(diff(x)), na.rm=TRUE) )
names(decades.ad)[1] <- 'Decade'
decades.ad.melt <- melt(decades.ad,id="Decade", variable.name='Country', value.name='R minus G abs diff')
ggplot(decades.ad.melt,aes(x=Decade,y=`R minus G abs diff`,colour=Country,group=Country)) + geom_line(size=1)

decades.aad <- aggregate( rminusg[,-1], list(cut(rminusg$DATE, cuts, labs )), function(x) mean(abs(x-mean(x,na.rm=TRUE)), na.rm=TRUE) )
names(decades.aad)[1] <- 'Decade'
decades.aad.melt <- melt(decades.aad,id="Decade", variable.name='Country', value.name='R minus G aad')
g2 <- ggplot(decades.aad.melt,aes(x=Decade,y=`R minus G aad`,colour=Country,group=Country)) +
        geom_line(size=1) + ggtitle('Decadal AAD')

multiplot(g1,g2)
pdf('~/Dropbox/2017/research/debtLimits/charts/rminusg.pdf')
multiplot(g1,g2)
dev.off()
