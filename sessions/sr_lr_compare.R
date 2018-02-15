l.thresh.q <- lapply( list( CAN=c(.2225, .58), DEU=c(.27,.33),
                            FRA=c(.0852,.175),GBR=c(.086, .169),
                            USA=c(.0475, .19) ), function(x) x * 4 )
    # The four-lag thresholds from the quarterly model

l.thresh.hist <- list( CAN=c(.35, 2.2), DEU=c(.65, 1.02), FRA=c(.45, .92),
                    GBR=c(.16, .46), USA=c(.05, .45), JPN=c(2.18, 2.7),
                    ITA=c(.98, 1.4) )
    # The one-lag thresholds from the historical data

load( 'data/rmg_est.rdata')
pt.q <- 4 * sapply(l.cty.mle.1$est, function(x) diff( x$mu ) )
load( 'data/rmg_hist__1880_est.rdata')
pt.hist <- sapply(l.cty.mle.1$est, function(x) diff( x$mu ) )

nn <- names(pt.hist)[c(1:4,7,5:6)]
pt.q <- sapply(nn, function(x) if( is.null( pt.q[x] ) )NA else pt.q[x] )
thresh.hist <- sapply(nn, function(x) l.thresh.hist[[x]] )
thresh.q <- sapply(nn, function(x) if( is.null( l.thresh.q[[x]] ) )
                                  c(NA,NA) else l.thresh.q[[x]] )

thresh.hist[,'CAN'] <- pt.hist['CAN'] <- NA
    # Drop Canada because it is crazy.
plt.dta <- rbind( pt.hist[nn], pt.q, thresh.hist, thresh.q )

v.bdy <- 90
bdy <- if( v.bdy == 95) c(4,6) else c(3,5) #  For 90% CI set to c(3,5), for 95% set to c(4,6)
plot( 1:ncol(plt.dta), plt.dta[1,], xaxt = "n", xlab='Country',
      ylim=range(plt.dta[c(1,bdy),], na.rm = TRUE), ylab='Interest-growth differential',
      pch=16, col='blue' )
axis(1, at=1:ncol(plt.dta), labels=nn)
abline(h=0)
points( 1:ncol(plt.dta), plt.dta[2,], pch=16, col='red' )
points( 1:ncol(plt.dta), plt.dta[bdy[1],], pch=3, col='blue' )
# points( 1:ncol(plt.dta), plt.dta[4,], pch='O', col='blue' )
points( 1:ncol(plt.dta), plt.dta[bdy[2],], pch=3, col='red' )
# points( 1:ncol(plt.dta), plt.dta[6,], pch='O', col='red' )
legend( 'topleft', c( 'Post-war quarterly data: point estimate',
                      paste0('Post-war quarterly data: ', v.bdy, '% critical value'),
                      '1880-present annual data: point estimate',
                      paste0('1880-present annual data: ', v.bdy, '% critical value') ),
                      bty='n', col=rep(c('red','blue'),each=2), pch=list(16,3))

