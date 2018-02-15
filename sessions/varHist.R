rm(list=ls())
library(VARext)
library(debtLimits)
library(xtable)

#### TO DO: REWORK WITH BETTER OPTIMAL LAG CHOICE ####
#### ALSO: MAKE SURE THAT THE THRESHOLDS ARE CORRECT. THEY ARE NOT (eg. France) ####

ltr <-  TRUE # FALSE #
st.lr <- if(ltr) 'hist_ltr_' else 'hist_'
start.year <- 1880 # 1880 # 1960 # 1880

if( ltr){
  l.thresh.1 <- list( CAN=c(2.05,4),
                      DEU=c(1.85, 2.4),
                      FRA=c(3,3),
                      GBR=c(1.92, 3.2),
                      UK=c(1.92, 3.2),
                      USA=c(1.4, 3),
                      JPN=c(0.95, 1.9),
                      ITA=c(2.3, 2.9) )
  # The approximate thresholds for the one-lag set-up
  l.thresh.opt <- list( CAN=c(1.2, 2.8),
                        DEU=c(1.85, 2.4),
                        FRA=c(3.3, 5),
                        GBR=c(1.4, 2.3),
                        USA=c(1.4, 3),
                        JPN=c(0.95,1.9),
                        ITA=c(4.35, 5) )
  # The approximate thresholds for the optimal-lag set-up
  l.thresh.4 <- list( CAN=c(1.1, 2.5),
                      DEU=c( 1.9, 3.2),
                      FRA=c( 2.72, 3.9),
                      GBR=c(1.4, 2.3),
                      USA=c(.67, 1.2),
                      JPN=c(3,3),
                      ITA=c(3.2, 4.5) )
  # Likewise in the four-lag case.
}else{
  l.thresh.1 <- list( CAN=c(.35, 2.2),
                      DEU=c(.65, 1.02),
                      FRA=c(.45, .92),
                      GBR=c(.16, .46),
                      USA=c(.05, .45),
                      JPN=c(2.18, 2.7),
                      ITA=c(.98, 1.4) )
  # The approximate thresholds for the one-lag set-up
  l.thresh.opt <- list( CAN=c(.35, 2.2),
                        DEU=c(.65, 1.02),
                        FRA=c(.92, 1.67),
                        GBR=c( .2, 1.1),
                        USA=c(-.068,.36),
                        JPN=c(2.18, 2.7),
                        ITA=c(2.67, 3.87) )
  # The approximate thresholds for the optimal-lag set-up
  l.thresh.4 <- list( CAN=c(2.15, NA),
                      DEU=c( .37, .79),
                      FRA=c( .26, .81),
                      GBR=c(-.03, .58),
                      USA=c(-.065,.38),
                      JPN=c(3.75, NA),
                      ITA=c(2.67, 3.87) )
  # Likewise in the four-lag case.
}

# To compute these numbers, start with guesses at zero.  Then run the code and
# look at the 5% (or whatever) thresholds in the Wald charts.  Put these in
# above and re-run.

l.cty.mle.1 <- list()
l.cty.mle.4 <- list()
l.cty.mle.rest.1 <- list()
l.cty.mle.rest.4 <- list()

# cts <- c('USA', 'FRA', 'GBR', 'DEU', 'JPN', 'ITA', 'CAN')
# cts <- c('USA', 'GBR', 'FRA', 'DEU', 'JPN', 'ITA', 'CAN')

# browser()

n.na <- n.obs <- v.lag <- rep(0, length(cts))
names(n.na) <- names(n.obs) <- names(v.lag) <- cts
l.keep <- list()

for(cty in cts){
  rg.dta <- hist.read(cty,start.year, ltr = ltr)
  n.na[cty] <- sum(is.na(rg.dta$gth-rg.dta$rfr))
  n.obs[cty] <- sum(!is.na(rg.dta$gth-rg.dta$rfr))
  # keep <- apply(rg.dta[,c('rfr','gth')], 1, function(x) !any(c(is.na(x), x['gth'] > 40)) )
  keep <- abs( rg.dta$rfr - rg.dta$gth - mean(rg.dta$rfr - rg.dta$gth, na.rm=T) ) < 6 *sd(rg.dta$rfr - rg.dta$gth, na.rm = TRUE )
  keep[is.na(keep)] <- FALSE
  l.keep[[cty]] <- rg.dta$year[keep]
  if(any(!keep)) message('Warning: ', cty, ' ', nrow(rg.dta) - sum(keep), ' years dropped due to missingness' )
  rg.dta <- rg.dta[ keep, ]
  n.var <- 2
  sim <- t( as.matrix( rg.dta[,c('gth','rfr')] ) )
  vv <- VAR(t(sim), ic='AIC', lag.max=5)
  v.lag[cty] <- vv$p
}

print(v.lag)
print(n.obs)
print(n.na)

keep.ranges <- lapply( l.keep,
                        function(x){
                          start <- c(1, which(diff(x) != 1 & diff(x) != 0) + 1)
                          end <- c(start - 1, length(x))
                          return( cbind(x[start],x[end]) )
                        }  )
df.ranges <- data.frame( Country=cts, Observations=n.obs,
                         Years=sapply( keep.ranges, function(x)
                           paste( apply(x, 1, paste, collapse='-'), collapse=', ' ) ) )
print(df.ranges)
print( xtable( df.ranges, digits = 0, label = 'tab:years_kept',
               caption = 'Sample periods for interest and growth rates' ), include.rownames=FALSE,
       file=paste0('~/Dropbox/2017/research/debtLimits/tables/',st.lr,'years_kept.tex' ) )

cts <- c('USA', 'GBR', 'FRA', 'DEU' )

for( cty in cts ){

  it <- 0
  for( lags in c(1,4,v.lag[cty] ) ) { #) ){

    it <- it + 1 ;
    opt.lag <- if( it == 3 ) TRUE else FALSE
    opt.st <- if( opt.lag ) 'opt_' else lags
        # The optimal lag indicator and string

    rg.dta <- hist.read(cty,start.year, ltr = ltr)
    keep <- apply(rg.dta[,c('rfr','gth')], 1, function(x) !any(c(is.na(x), x['gth'] > 40)) )
    if(any(!keep)) message('Warning: ', cty, ' ', lags, ' lags, ',
                           nrow(rg.dta) - sum(keep), ' years dropped due to missingness' )
    rg.dta <- rg.dta[ keep, ]
    n.var <- 2
    sim <- t( as.matrix( rg.dta[,c('gth','rfr')] ) )
    varnames <- c( 'Growth', 'Int. rate' )

    # Estimate VAR via MLE
    l.var.ols <- var.ols( sim, lags )
    l.var.est.mle <- var.mle.est( sim, lags )
    l.var.est.mle.uc <- var.mle.est( sim, lags, cond = FALSE )
    irf.plot( l.var.est.mle, n.pds = 16, varnames = varnames )
    # To check: LR variance for lag > 1
    g <- function(par, lags, n.var=2, theta=0){
      # Imposes restriction that the mean difference in a 2-variable VAR is greater than theta.
      l.var <- par.from(par, n.var, lags)
      mu <- mu.calc(l.var$a, l.var$A)
      return( - mu[2] + mu[1] + theta )
    }
    thet <- if(lags==1) l.thresh.1[[cty]][1] else l.thresh.4[[cty]][1]
    l.var.rest.mle <- var.mle.rest( sim, g, lags, theta=thet )
    l.var.rest.mle.uc <- var.mle.rest( sim, g, lags, cond=TRUE, theta=thet )

    irf.plot( l.var.rest.mle, n.pds = 16, varnames = varnames )
    # Restricted VAR with mu >= 0
    v.mle <- var.mle.se( sim, l.var.est.mle )
    v.mle.uc <- var.mle.se( sim, l.var.est.mle.uc, cond=FALSE )
    # Not working atm
    # v.mle.uc <- var.mle.se( sim, l.var.est.mle.uc, ineff = TRUE )
    v.mle.rest <- var.mle.se( sim, l.var.rest.mle, ineff=TRUE )
    v.mle.rest.uc <- var.mle.se( sim, l.var.rest.mle.uc, cond=FALSE, ineff=TRUE )
    # v.mle.rest.2 <- var.mle.se( sim, l.var.rest.mle.2, ineff=TRUE )
    v.ols <- var.ols.se( sim, lags )
    v.ols.ext <- lapply( v.ols, function(x) {
      out <- matrix( 0, nrow(v.mle), ncol(v.mle) )
      out[1:(n.var*(1+n.var*lags)),1:(n.var*(1+n.var*lags))] <- x
      return(out)
    } )
    print( cbind( sapply(v.ols.ext, function(x) sqrt(diag(x)) ),
                  v.mle=sqrt(diag(v.mle)), v.mle.uc=sqrt(diag(v.mle.uc)),
                  v.mle.rest=sqrt(diag(v.mle.rest)), v.mle.rest.uc=sqrt(diag(v.mle.rest.uc))) )

    # Plot vs. simulations
    sim.fcast.plot(sim, lags, start=lags, l.var = l.var.est.mle, x.vals = rg.dta$year,
                   varnames = varnames, n.split = 15, n.pds=10 )
    pdf(paste0('~/Dropbox/2017/research/debtLimits/charts/fcast_', cty, '_cond_',
               start.year, st.lr, opt.st, 'lag.pdf'))
      sim.fcast.plot(sim, lags, start=lags, l.var = l.var.est.mle, x.vals = rg.dta$year,
                     varnames = varnames, n.split = 15, n.pds=10 )
    dev.off()
    sim.fcast.plot(sim, lags, start=lags, l.var = l.var.est.mle.uc, x.vals = rg.dta$year,
                   varnames = varnames, n.split = 15, n.pds=10 )
    pdf(paste0('~/Dropbox/2017/research/debtLimits/charts/fcast_', cty, '_uncond_',
               start.year, st.lr, opt.st, 'lag.pdf'))
      sim.fcast.plot(sim, lags, start=lags, l.var = l.var.est.mle.uc, x.vals = rg.dta$year,
                     varnames = varnames, n.split = 15, n.pds=10 )
    dev.off()
    sim.fcast.plot(sim, lags, start=lags, l.var = l.var.rest.mle, x.vals = rg.dta$year,
                   varnames = varnames, n.split = 15, n.pds=10 )
    # sim.fcast.plot(sim, lags, start=lags, l.var = l.var.rest.mle.2, x.vals = rg.dta$date,
    #                varnames = varnames,n.split = 39, n.pds=20 )

    # Create latex output
    l.l.var <- list( l.var.ols, l.var.est.mle, l.var.est.mle.uc, l.var.rest.mle, l.var.rest.mle.uc )
    l.m.var <- list( v.ols$ols, v.mle, v.mle.uc, v.mle.rest, v.mle.rest.uc )
    v.cond <- c(T,T,F,T,F)
    v.lhood <- sapply( 1:5, function(i) var_lhood_N( sim, par.to.l(l.l.var[[i]]), lags,
                                                     cond = v.cond[i] ) )
    var.table( l.l.var, l.m.var, file=paste0('./tables/', cty, '_', start.year, st.lr, opt.st, 'lag.tex'),
               specnames = c( 'OLS', 'MLE', 'UMLE', 'Rest MLE', 'Rest UMLE' ), varnames = varnames,
               caption = paste0( lags, '-lag test VAR for ', cty,
                                 '. In restricted MLEs, mean difference is ', round(thet, 4),
                                 if(opt.lag) paste0(' Using AIC opimal lag length ', lags) else NULL ),
               label=paste0('tab:', cty, opt.st), footer=TRUE, v.lhood=v.lhood )

    # Model tests
    lb <- min( diff(l.var.est.mle$mu), diff(l.var.est.mle.uc$mu) )
    ub <- max( -lb, 3 )
    if( cty == 'JPN' || cty == 'ITA' ) ub <- if(lags==1) 3 else 5
    if( cty == 'GBR' || cty == 'FRA' || cty =='CAN' ) ub <- 4
    if( cty == 'DEU' & lags > 1 ) ub <- 4.5
    n.theta <- 20
    v.theta <- seq( lb, max(.2,ub), length.out=n.theta)
        # Range of theta
    drop <- NULL
    # if( cty == 'USA' & lags == v.lag[cty] ) drop <- n.theta
    # if( cty == 'GBR' & lags == v.lag[cty] ) drop <- c(n.theta-1,n.theta)
    # if( cty == 'GBR' & lags == v.lag[cty] ) drop <- n.theta - 7:8
    # if( cty == 'CAN' ) drop <- n.theta - 2 #:0
        # Country-specific adjustments
    # browser()
    # dev.off()

    pdf(paste0('~/Dropbox/2017/research/debtLimits/charts/W_LR_test_', cty, '_', start.year,
               st.lr, opt.st, 'lag.pdf'))
      lr.wald.plot( sim, lags, g, v.theta, xlab='Mean interest-growth differential', drop=drop )
      # Plot the tests
      # if( lags == 1 ){
      #   abline(v=l.thresh.1[[cty]] )
      # }else if( lags == 4 ){
      #   abline(v=l.thresh.4[[cty]][1] )
      # }else{
      #   abline(v=l.thresh.opt[[cty]] )
      # }
      abline( v = 0, lty=2 )
    dev.off()

    if( lags == 1){
      l.cty.mle.1$est[[cty]] <- l.var.est.mle.uc
      l.cty.mle.1$se[[cty]] <- v.mle.uc
      l.cty.mle.1$lhood[cty] <- v.lhood[3]
      l.cty.mle.rest.1$est[[cty]] <- l.var.rest.mle.uc
      l.cty.mle.rest.1$se[[cty]] <- v.mle.rest.uc
      l.cty.mle.rest.1$lhood[cty] <- v.lhood[5]
    }else if( lags == 4 ){
      l.cty.mle.4$est[[cty]] <- l.var.est.mle.uc
      l.cty.mle.4$se[[cty]] <- v.mle.uc
      l.cty.mle.4$lhood[cty] <- v.lhood[3]
      l.cty.mle.rest.4$est[[cty]] <- l.var.rest.mle.uc
      l.cty.mle.rest.4$se[[cty]] <- v.mle.rest.uc
      l.cty.mle.rest.4$lhood[cty] <- v.lhood[5]
    }
  }
}

var.table( l.cty.mle.1$est, l.cty.mle.1$se,
           file=paste0('~/Dropbox/2017/research/debtLimits/tables/all_', st.lr, '1lag.tex'),
           # file=paste0('./tables/all_', start.year, st.lr, '1lag.tex'),
           specnames = c('USA', 'United Kingdom', 'France', 'Germany' ), #, 'Japan', 'Italy', 'Canada'),
           varnames = varnames, n.obs = n.obs[cts],
           caption =paste0( '1-lag VAR for sample of countries. Annual data ', start.year,
                      '-2015. Robust likelihood-based standard errors in parentheses.'),
           label=paste0('tab:all_hist_1lag'), footer=TRUE, v.lhood=l.cty.mle.1$lhood,
           font.size='small' )

var.table( l.cty.mle.4$est, l.cty.mle.4$se,
           file=paste0('~/Dropbox/2017/research/debtLimits/tables/all_', st.lr, '4lag.tex'),
           # file=paste0('./tables/all_', start.year, st.lr, '4lag.tex'),
           specnames = c('USA', 'United Kingdom', 'France', 'Germany' ), #,  'Japan', 'Italy', 'Canada'),
           varnames = varnames, n.obs = n.obs[cts],
           caption =paste0( 'Multi-lag VAR for sample of countries. Annual data ', start.year,
                            '-2015. Robust likelihood-based standard errors in parentheses.'),
           label=paste0('tab:all_hist_4lag'), footer=TRUE, v.lhood=l.cty.mle.4$lhood,
           font.size='small'  )

var.table( l.cty.mle.rest.4$est, l.cty.mle.rest.4$se,
           file=paste0('~/Dropbox/2017/research/debtLimits/tables/all_', st.lr, 'rest_4lag.tex'),
           # file=paste0('./tables/all_', start.year, st.lr, 'rest_4lag.tex'),
           specnames = c('USA', 'United Kingdom', 'France', 'Germany' ), #, 'Japan', 'Italy','Canada'),
           varnames = varnames, n.obs = n.obs[cts],
           caption =paste0( 'Multi-lag restricted VAR for sample of countries. Lon-run means are
                            restricted to be at the 5\\% critical value under the unconditional
                            LR test. Annual data ', start.year,
                            '-2015. Robust likelihood-based standard errors in parentheses'),
           label=paste0('tab:all_hist_rest_4lag'), footer=TRUE, v.lhood=l.cty.mle.rest.4$lhood,
           font.size='small'  )

rmg.all <- sapply( l.cty.mle.1$est, function(x) diff(x$mu) )

save('cts',  'l.cty.mle.1', 'l.cty.mle.4', 'l.cty.mle.rest.1', 'l.cty.mle.rest.4',
     'l.thresh.1', 'l.thresh.4', 'v.lag',
     file = paste0('data/rmg_', st.lr, '_', start.year, '_est.rdata') )
