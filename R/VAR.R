#####################################################################
# VAR.R
#
# Contains the functions for estimating and plotting the VAR for R and G
#
# 24apr2017
# Philip Barrett, Washington DC
#####################################################################

rg.read <- function( cty = 'USA', start.date = "1960-01-01" ){
## Read and cleans in the data for the specified country
  rfr <- read.csv('data/riskfreerates.csv')
  gth <- read.csv('data/growthrates.csv')
      # Read the data
  cty.gth <- data.frame( date=as.Date(gth$DATE), gth=gth[[cty]] )
  cty.rfr <- data.frame( date=as.Date(rfr$DATE), rfr = 100 * ( (1+rfr[[cty]]/100) ^ .25 - 1 ) )
      # Create country-specific dataframes
  cty.dta <- merge( subset( cty.gth, date > start.date ),
                    subset( cty.rfr, date > start.date ) )
      # The country data after the start date
  cty.dta$rmg <- apply( cty.dta[,-1], 1, diff )
      # Create R minus G
  return( cty.dta )
}

rg.plot <- function( dta, yr.filter=10, save.pdf=FALSE, cty.name=NULL ){
# Plot the levels and differences of R and G
  n.filter <- 4 * yr.filter
      # The number of periods for the filter

  par.dft <- par('mar')
  if(save.pdf) pdf( paste0('~/Dropbox/2017/research/debtLimits/charts/rg_', cty.name, '.pdf') )
  par(mfrow=c(2,1), mar=c(3,3,3,3))
  plot(dta$date, dta$gth, type='l', lwd=2, main='Nominal GDP growth', xlab='', ylab='' )
  abline(h=0, lwd=.5)
  plot(dta$date, dta$rfr, type='l', lwd=2, main='Nominal interest rate', xlab='', ylab='' )
  abline(h=0, lwd=.5)
  if(save.pdf) dev.off()
  par(mfrow=c(1,1), mar=par.dft)
      # The levels of the series

  if(save.pdf) pdf( paste0('~/Dropbox/2017/research/debtLimits/charts/rmg_', cty.name, '.pdf') )
  par(mar=c(3,3,3,3))
  plot(dta$date, dta$rmg, type='l', lwd=2, xlab='', ylab='Interest-growth differential' )
  lines(dta$date, filter( dta$rmg, rep( 1/n.filter, n.filter ), sides=1 ), col='blue', lwd=2 )
  abline( h=(-2):2, lwd=.5, lty=2)
  legend( 'topright', paste0( yr.filter, ' year rolling average'), bty='n', col='blue', lwd=2 )
  abline(h=0, lwd=.5)
  abline(h=mean(dta$rmg,na.rm=TRUE), lty=2, col='blue')
  if(save.pdf) dev.off()
  par(mfrow=c(1,1), mar=par.dft)
}

var.rg.est <- function( dta ){
# Estimate a VAR based on R and G
  nn <- 2
  cty.var <- VAR(dta[,c('gth', 'rfr')])
      # The VAR
  coeff <- sapply( c('gth', 'rfr'), function(x) cty.var$varresult[[x]]$coeff)
  Sigma <-var(sapply( c('gth', 'rfr'), function(x) cty.var$varresult[[x]]$residuals))
      # Extract the coefficients
  A <- t( coeff[1:nn,] )
  a <- coeff[nn+1,]
  mu <- solve( diag(nn) - A, a )
      # The parameters
  return( list( A=A, a=a, mu=mu, Sigma=Sigma ) )
}

plot.var.rg.resids <- function(dta, yr.filter=2, save.pdf=FALSE, cty.name=NULL){
# Plots the VAR residuals
  n.filter <- 4 * yr.filter
      # The number of periods for the filter
  cty.var <- VAR(dta[,c('gth', 'rfr')])
      # The VAR

  par.dft <- par('mar')
  if(save.pdf) pdf( paste0('~/Dropbox/2017/research/debtLimits/charts/rg_resids_',
                           cty.name, '.pdf') )
  par(mfrow=c(2,1), mar=c(3,3,3,3))
  plot(dta$date[-1], cty.var$varresult$gth$residuals, type='l', lwd=1,
       main='VAR residuals: Nominal GDP', xlab='', ylab='' )
  lines(dta$date[-1], filter( cty.var$varresult$gth$residuals, rep(1/n.filter,n.filter)),
        col='blue', lwd=2 )
  abline(h=0)
  plot(dta$date[-1], cty.var$varresult$rfr$residuals, type='l', lwd=1,
       main='VAR residuals: Nominal interest rate',xlab='', ylab='' )
  lines(dta$date[-1], filter( cty.var$varresult$rfr$residuals, rep(1/n.filter,n.filter)),
        col='blue', lwd=2 )
  abline(h=0)
  if(save.pdf) dev.off()
  par(mfrow=c(1,1), mar=par.dft)

}

plot.level.fcast <- function(dta, l.var, yr.fcast=10, save.pdf=FALSE, cty.name=NULL ){
#  Plot the level of interest rates and growth with forecasts from the VAR
  n.fcast <- max( 20, 4 * yr.fcast )
      # Number of forecast periods (at least 5 yrs)
  n.pds <- nrow(dta)
      # The number of periods in total
  pd.counter <- 1
  fcast.pd <- 0
      # Counters for the loop
  fcast.low <- fcast.up <- fcast.mean <- matrix( NA, 2, n.pds )
  sd.diff <- rep(0,n.pds)
      # Containers for the forecasts
  sel.mat <- matrix( c( 1, -1, -1, 1 ), nn, nn )
      # Weights on the variances and covariances in the difference
  while( pd.counter < n.pds ){
    fcast.low[,pd.counter] <- fcast.up[,pd.counter] <-
      fcast.mean[,pd.counter] <- t(dta[pd.counter,c('gth','rfr')])
        # First period
    while( fcast.pd < n.fcast - 8  & pd.counter < n.pds ){
      fcast.pd <- fcast.pd + 1
      pd.counter <- pd.counter + 1
          # Increment counters
      fcast.mean[,pd.counter] <- l.var$a + l.var$A %*% fcast.mean[,pd.counter-1]
          # Update the mean
      v.sd <- sqrt( diag( Reduce( '+', lapply(1:fcast.pd,
              function(j) ( l.var$A %^% (j-1) ) %*% l.var$Sigma %*% ( t(l.var$A) %^% (j-1) ) ) ) ) )
      fcast.up[,pd.counter] <- fcast.mean[,pd.counter] + v.sd
      fcast.low[,pd.counter] <- fcast.mean[,pd.counter] - v.sd
      sd.diff[pd.counter] <- sqrt( sum( Reduce( '+', lapply(1:fcast.pd,
            function(j) ( l.var$A %^% (j-1) ) %*% l.var$Sigma %*% ( t(l.var$A) %^% (j-1) ) ) )
              * sel.mat ) )
    }
    pd.counter <- pd.counter + 8
    fcast.pd <- 0
  }

  par.dft <- par('mar')
  if(save.pdf) pdf( paste0('~/Dropbox/2017/research/debtLimits/charts/rg_fcasts_',
                           cty.name, '.pdf') )
  par(mfrow=c(2,1), mar=c(3,3,3,3))
  plot(dta$date, dta$gth, type='l', lwd=2, main='Nominal GDP growth', xlab='', ylab='' )
  lines(dta$date, fcast.mean[1,], lwd=2, col='blue' )
  lines(dta$date, fcast.up[1,], lwd=1, col='red', lty=2 )
  lines(dta$date, fcast.low[1,], lwd=1, col='red', lty=2 )
  abline(h=0, lwd=.5)
  abline(h=l.var$mu[1], lwd=.5, lty=2, col='blue')
  plot(dta$date, dta$rfr, type='l', lwd=2, main='Nominal interest rate', xlab='', ylab='' )
  lines(dta$date, fcast.mean[2,], lwd=2, col='blue' )
  lines(dta$date, fcast.up[2,], lwd=1, col='red', lty=2 )
  lines(dta$date, fcast.low[2,], lwd=1, col='red', lty=2 )
  abline(h=0, lwd=.5)
  abline(h=l.var$mu[2], lwd=.5, lty=2, col='blue')
  if(save.pdf) dev.off()
  par(mfrow=c(1,1), mar=par.dft)

  diff.mean <- fcast.mean[2,] - fcast.mean[1,]
  if(save.pdf) pdf( paste0('~/Dropbox/2017/research/debtLimits/charts/rmg_fcasts_',
                           cty.name, '.pdf') )
  plot(dta$date, dta$rfr - dta$gth, type='l', lwd=2, main='Interest-growth differential',
       xlab='', ylab='' )
  lines(dta$date, diff.mean, lwd=2, col='blue' )
  lines(dta$date, diff.mean + sd.diff, lty=2, col='red' )
  lines(dta$date, diff.mean - sd.diff, lty=2, col='red' )
  abline(h=0, lwd=.5)
  abline(h=diff(l.var$mu), lwd=.5, lty=2, col='blue')
  if(save.pdf) dev.off()
}

var.discretize <- function(l.var, dta=NULL, n.int=1e6, n.pts = 3, n.dirs = 8 ){
# Discretizes a VAR object

  nn <- length(l.var$a)
      # Number of variables
  if(is.null(dta)){
    Sig.u <- Reduce( '+',
          lapply( 0:1000, function(i) (l.var$A %^% i) %*% l.var$Sigma %*% ( t(l.var$A) %^% i) ) )
  }else{
    Sig.u <- var(dta[,c('gth','rfr')])
  }
      # The unconditional variance
  q <- .5 / n.pts
  l.0 <- qnorm( seq( .5 + 4 * q / 5, 1 - q / 5, length.out=n.pts ), 0, 1 )
  theta.0 <- seq( 0, 2*pi, length.out=n.dirs+1 )[-(n.dirs+1)]
      # The vector of angles
  zz.0 <- cbind( sin(theta.0), cos(theta.0) )
  pts.0 <- rbind( c(0,0),
                  do.call( 'rbind', lapply( l.0, function(l) l * zz.0 ) ) )
      # The matrix of points
  n.X <- 1 + n.pts * n.dirs
  X <- t( t( chol( Sig.u ) ) %*% t( pts.0 ) ) + rep(1, n.X ) %*% t(l.var$mu)
      # Rescale the points to the mean and unconditional variance
  X[X[,2]<0,2] <- 0
      # Prevent negative rfr
  int.eps <- rmvnorm( n.int, 0 * l.var$a, l.var$Sigma )
      # The vector of integration shocks
  trans <- matrix(0, n.X, n.X)
      # Initiate the transition probability matrix
  for( i in 1:n.X ){
    X.cont <- rep(1,n.int) %*% t( l.var$a + l.var$A %*% X[i,] ) + int.eps
        # The stochastic continuation values of X
    nearest <- nn2( X, X.cont, 1 )
        # The list of nearest neighbours
    trans[i,] <- table(factor(nearest$nn.idx,levels=1:n.X) ) / n.int
  }
  m <- (trans %^% 100)[1,]
      # The unconditional distribution
  rmg <- matrix( X[,'rfr'], nrow(X), ncol(X) ) - matrix( X[,'gth'], nrow(X), ncol(X), byrow=TRUE )
  m.rmg <- m * trans
  d.rmg <- cbind( c(rmg), c(m.rmg) )
  d.rmg <- d.rmg[ order(d.rmg[,1]), ]
      # The unconditional distribution for r minus g
  return( list( X=X, trans=trans, m=m, d.rmg=d.rmg ) )
}

var.data.interp <- function( l.var, dta, plot.on=TRUE ){
# Interprets data as values from a VAR and plots the corresponding series
  s.idx <- nn2( l.var$X, dta[,c('gth','rfr')], 1)$nn.idx
  dta.disc <- data.frame(l.var$X[s.idx,])
  colnames(dta.disc) <- c('gth','rfr')
      # Create the nearest neighbours
  if(plot.on){
    par.dft <- par('mar')
    par(mfrow=c(2,1), mar=c(3,3,3,3))
    plot(dta$date, dta$gth, type='l', lwd=2, main='Nominal GDP growth', xlab='', ylab='' )
    lines(dta$date, dta.disc$gth, col='blue', lwd=1 )
    abline(h=0, lwd=.5)
    plot(dta$date, dta$rfr, type='l', lwd=2, main='Nominal interest rate', xlab='', ylab='' )
    lines(dta$date, dta.disc$rfr, col='blue', lwd=1 )
    abline(h=0, lwd=.5)
    par(mfrow=c(1,1), mar=par.dft)

    plot(dta$date, dta$rfr - dta$gth, type='l', lwd=2, main='Interest-growth differential', xlab='', ylab='' )
    lines(dta$date, dta.disc$rfr - dta.disc$gth, col='blue', lwd=1 )
    abline(h=0, lwd=.5)
  }
  return(list(s.idx=s.idx,dta.disc=dta.disc))
}
