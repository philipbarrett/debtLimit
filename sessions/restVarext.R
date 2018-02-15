library(VARext)
library(debtLimits)
rm(list=ls())

lr <- FALSE                           # Long run rates flag
hist <- FALSE                         # Historical data flag
st.lr <- if(lr) 'lr_' else (if(hist) 'hist_' else '')

#### TO DO: ADD ITALY ####

l.thresh.1 <- list( CAN=c(.07275, .271)* 4,
                    DEU=c(.375,.375)* 4,
                    FRA=c(.025,.1925)* 4,
                    GBR=c(.1,.325)* 4,
                    USA=c(.0275, .1775)* 4,
                    JPN=c(.0275, .1775)* 4,
                    ITA=c(.0275, .1775) * 4 )
    # The approximate thresholds for the one-lag set-up
l.thresh.4 <- list( CAN=c(.2525, .58)* 4,
                    DEU=c(.26,.33)* 4,
                    FRA=c(.0,.175)* 4,
                    GBR=c(.06, .169)* 4,
                    USA=c(.015, .19)* 4,
                    JPN=c(.0275, .1775)* 4,
                    ITA=c(.0275, .1775)* 4 )
    # Likewise in the four-lag case

l.cty.mle.1 <- list()
l.cty.mle.4 <- list()
l.cty.mle.rest.1 <- list()
l.cty.mle.rest.4 <- list()

cts <- c( 'USA', 'FRA', 'GBR', 'DEU', 'CAN', 'JPN', 'ITA' )
# cts <- c('USA', 'FRA', 'GBR', 'DEU', 'CAN')
# cts <- c('USA', 'GBR', 'CAN')
start.date <- c( 'USA'='1956-01-01', 'FRA'='1956-01-01', 'GBR'='1956-01-01', 'CAN'='1961-01-01', 'DEU'='1970-04-01',
                 'JPN'='1994-01-01', 'ITA'='1991-01-01' )


n.na <- n.obs <- v.lag <- rep(0, length(cts))
names(n.na) <- names(n.obs) <- names(v.lag) <- cts

for(cty in cts){
  rg.dta <- rg.read(cty,start.date[cty])
  n.na[cty] <- sum(is.na(rg.dta$rmg))
  n.obs[cty] <- sum(!is.na(rg.dta$rmg))
  keep <- apply(rg.dta[,c('rfr','gth')], 1, function(x) !any(c(is.na(x), x['gth'] > 40)) )
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

cts <- c( 'USA', 'GBR', 'FRA', 'CAN' )

for( cty in cts ){
  for( lags in c(1,4) ){

    rg.dta <- if(lr) rg10.read( cty ) else ( if(hist) hist.read(cty,start.year) else rg.read( cty, start.date = start.date[cty]  ) )
    rg.dta <- rg.dta[ apply(rg.dta, 1, function(x) !any(is.na(x))), ]
    n.var <- 2
    sim <- t( as.matrix( rg.dta[,c('gth','rfr')] ) ) * 4
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
    sim.fcast.plot(sim, lags, start=lags, l.var = l.var.est.mle, x.vals = rg.dta$date,
                   varnames = varnames, n.split = 15, n.pds=10 )
    pdf(paste0('~/Dropbox/2017/research/debtLimits/charts/fcast_', cty, '_', st.lr, lags, '_lag.pdf'))
    sim.fcast.plot(sim, lags, start=lags, l.var = l.var.est.mle.uc, x.vals = rg.dta$date,
                   varnames = varnames, n.split = 39, n.pds=20 )
    dev.off()
    sim.fcast.plot(sim, lags, start=lags, l.var = l.var.rest.mle, x.vals = rg.dta$date,
                   varnames = varnames, n.split = 39, n.pds=20 )
    # sim.fcast.plot(sim, lags, start=lags, l.var = l.var.rest.mle.2, x.vals = rg.dta$date,
    #                varnames = varnames,n.split = 39, n.pds=20 )

    # Create latex output
    l.l.var <- list( l.var.ols, l.var.est.mle, l.var.est.mle.uc, l.var.rest.mle, l.var.rest.mle.uc )
    l.m.var <- list( v.ols$ols, v.mle, v.mle.uc, v.mle.rest, v.mle.rest.uc )
    v.cond <- c(T,T,F,T,F)
    v.lhood <- sapply( 1:5, function(i) var_lhood_N( sim, par.to.l(l.l.var[[i]]), lags,
                                                     cond = v.cond[i] ) )
    var.table( l.l.var, l.m.var, file=paste0('./tables/', cty, '_', st.lr, lags, 'lag.tex'),
               specnames = c( 'OLS', 'MLE', 'UMLE', 'Rest MLE', 'Rest UMLE' ), varnames = varnames,
               caption = paste0( lags, '-lag test VAR for ', cty,
                          '. In restricted MLEs, mean difference is ', round(thet, 4)),
               label=paste0('tab:', cty, lags, 'lag'), footer=TRUE, v.lhood=v.lhood )

    # Model tests
    lb <- min( diff(l.var.est.mle$mu), diff(l.var.est.mle.uc$mu) )
    ub <- max( -lb-.02, 1.5 )
    if( cty == 'DEU' || cty == 'CAN' ) ub <- 1.5
    # if( cty == 'DEU' & lags > 1 ) ub <- .4
    n.theta <- 20
    pdf(paste0('~/Dropbox/2017/research/debtLimits/charts/W_LR_test_', cty, '_', st.lr, lags, '_lag.pdf'))
    v.theta <- seq( lb, max(.2,ub), length.out=n.theta)
    drop <- NULL
    # if( cty == 'FRA' & lags == 1 ) drop <- c(n.theta-1,n.theta)
    # if( cty == 'FRA' & lags == 4 ) drop <- n.theta
    # if( cty == 'DEU' & lags == 1) drop <- 12:17
        # Country-specific adjustments
    lr.wald.plot( sim, lags, g, v.theta, xlab='Mean interest-growth differential', drop=drop )
        # Plot the tests
    # if(!lr){
    #   if( lags == 1 ){
    #     abline(v=l.thresh.1[[cty]] )
    #   }else{
    #     abline(v=l.thresh.4[[cty]][1] )
    #   }
    # }
    abline( v = 0, lty=2 )
    dev.off()

    # if(lags==4) stop()

    if( lags == 1){
      l.cty.mle.1$est[[cty]] <- l.var.est.mle.uc
      l.cty.mle.1$se[[cty]] <- v.mle.uc
      l.cty.mle.1$lhood[cty] <- v.lhood[3]
      l.cty.mle.rest.1$est[[cty]] <- l.var.rest.mle.uc
      l.cty.mle.rest.1$se[[cty]] <- v.mle.rest.uc
      l.cty.mle.rest.1$lhood[cty] <- v.lhood[5]
    }else{
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
           n.obs=n.obs[cts],
           # file=paste0('./tables/all_', st.lr, '1lag.tex'),
           specnames = c('USA', 'United Kingdom', 'France', 'Canada'), varnames = varnames,
           caption ='1-lag VAR for sample of countries. Annualized quarterly data 1956Q1-2016Q4
                      (from 1961Q1 for Canada). Robust likelihood-based standard errors in parentheses.',
           label=paste0('tab:all_1lag'), footer=TRUE, v.lhood=l.cty.mle.1$lhood,
           font.size='small' )

var.table( l.cty.mle.4$est, l.cty.mle.4$se,
           file=paste0('~/Dropbox/2017/research/debtLimits/tables/all_', st.lr, '4lag.tex'),
           varnames = varnames,
           n.obs=n.obs[cts],
           # file=paste0('./tables/all_', st.lr, '4lag.tex'),
           specnames = c('USA', 'United Kingdom', 'France', 'Canada'),
           caption ='4-lag VAR for sample of countries. Annualized quarterly data 1956Q1-2016Q4
                      (from 1961Q1 for Canada). Robust likelihood-based standard errors in parentheses.',
           label=paste0('tab:all_4lag'), footer=TRUE, v.lhood=l.cty.mle.4$lhood,
           font.size='small'  )

var.table( l.cty.mle.rest.4$est, l.cty.mle.rest.4$se,
           file=paste0('~/Dropbox/2017/research/debtLimits/tables/all_', st.lr, 'rest_4lag.tex'),
           varnames = varnames,
           n.obs=n.obs[cts],
           # file=paste0('./tables/all_', st.lr, 'rest_4lag.tex'),
           specnames = c('USA', 'United Kingdom', 'France', 'Canada'),
           caption ='4-lag restricted VAR for sample of countries. Lon-run means are
                      restricted to be at the 5\\% critical value under the unconditional
                      LR test. Annualized quarterly data 1956Q1-2016Q4
                      (from 1961Q1 for Canada). Robust likelihood-based standard errors in parentheses.',
           label=paste0('tab:all_rest_1lag'), footer=TRUE, v.lhood=l.cty.mle.rest.4$lhood,
           font.size='small'  )

save('cts',  'l.cty.mle.1', 'l.cty.mle.4', 'l.cty.mle.rest.1', 'l.cty.mle.rest.4',
     'l.thresh.1', 'l.thresh.4', file = paste0('data/rmg', st.lr, '_est.rdata') )
