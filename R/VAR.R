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
  Sigma <-summary(cty.var)$covres
      # Extract the coefficients
  A <- t( coeff[1:nn,] )
  a <- coeff[nn+1,]
  mu <- solve( diag(nn) - A, a )
      # The parameters
  return( list( A=A, a=a, mu=mu, Sigma=Sigma, VAR=cty.var ) )
}

var.rest.rg.est <- function( dta, theta=0 ){
# Estimates the restricted VAR with R>=G+theta
  nn <- 2
      # Number of variables
  g <- function(x){
    a <- x[1:nn]
    B <- x[nn+1:(nn^2)]
    mu <- solve( diag(nn) - B, a )
        # the long run averages
    return(-diff(mu)+theta)
  }

  est <- var.rg.est( dta )
      # The unrestricted VAR
  par.init <- c( est$a, est$A, est$Sigma[lower.tri(est$Sigma, TRUE)] )
      # Initial parameters
  Y <- t(as.matrix(dta[,c('gth','rfr')]))
      # The data
  eval_f <- function( par ){
    var_lhood( Y, par )
  }
  opts <- list( algorithm = "NLOPT_LN_COBYLA", maxeval=1e6 )
  sol <- nloptr( x0 = par.init, eval_f = eval_f, eval_g_ineq = g,
                 ub = rep( Inf, length(par.init) ),
                 lb = rep( -Inf, length(par.init) ),
                 opts = opts )
      # The optimum
  if(!(sol$status %in% c(1,3,4))) message("NLOPTR failed to solve. \n", sol$message )
      # Warning message
  a <- sol$solution[1:nn]
  A <- matrix( sol$solution[nn+1:(nn^2)], nn, nn )
  mu <- solve( diag(nn) - A, a )
  Sigma <- 0 * A
  Sigma[ lower.tri(Sigma,TRUE)] <- Sigma[ upper.tri(Sigma,TRUE)] <- sol$solution[-(1:(nn+nn^2))]
  l <- sol$objective
  l.init <- eval_f(par.init)
      # Format the output
  l.ratio <- 2 * nrow(dta) * ( l - l.init )
  l.t.stat <- sapply( c( .9, .95, .975, .99 ), function(x) qchisq(x, 1) )
      # DoF = 1 because we only impose the restriction that R>=G
  names( l.t.stat ) <- c( .9, .95, .975, .99 )
  rest.rej <- l.ratio > l.t.stat
      # The critical values
  return(list( a=a, A=A, mu=mu, Sigma=Sigma, l=l, l.init=l.init,
               status=sol$message, l.ratio=l.ratio, l.t.stat=l.t.stat,
               rest.rej=rest.rej) )
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
  rmg <- matrix( X[,'rfr'], nrow(X), nrow(X) ) - matrix( X[,'gth'], nrow(X), nrow(X), byrow=TRUE )
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

grad.mu <- function( B, A ){
# Returns the numerical derivative of the mean w.r.t the coefficients of A and B
  nn <- nrow(B)
  mm <- nn ^ 2 + nn
  mu <- solve( diag(nn) -B, A )
      # The mean
  m.mu <- matrix( mu, nrow=mm, ncol=length(mu), byrow=TRUE )
      # The matrix of means
  inc <- 1e-06
  m.mu.inc <- 0 * m.mu
      # The incremented value of mu
  for( i in 1:length(A) ){
    A.cpy <- A
    A.cpy[i] <- A.cpy[i] + inc
    mu.new <- solve( diag(nn) - B, A.cpy )
    m.mu.inc[i,] <- mu.new
  }

  for( i in 1:length(B) ){
    B.cpy <- B
    B.cpy[i] <- B.cpy[i] + inc
    mu.new <- solve( diag(nn) - B.cpy, A )
    m.mu.inc[ length(A)+i, ] <- mu.new
  }

  out <- ( m.mu - m.mu.inc ) / inc
      # The output
  reorder <- c( sapply( 1:nn, function(i) c( i + nn * 1:nn, i ) ) )
      # Order should be: coeffs, then const for each equation
  out <- out[ reorder, ]

  return(out)
}

var.table <- function( l.var, file=NULL, add.mean=FALSE, varnames=NULL,
                       label=NULL, caption=NULL, add.se=TRUE, footer=TRUE,
                       start=NULL, end=NULL ){
# Creates a nice regression table for the VAR

  # Set up
  n.v <- length(l.var$varresult)      # Number of variables
  n.col <- 2 + n.v * 2 + add.mean     # Number of columns
  if( is.null(varnames) ) varnames <-
    gsub( '_', '', names(l.var$varresult[[1]]$coefficients[-(n.v+1)]) )
  coeff <- t( sapply( l.var$varresult, function(x) x$coefficients ) )
      # Coefficients
  mu <- solve( (diag(n.v)-coeff[,1:n.v]), coeff[,'const'] )
      # The mean

  if( add.se ){
    ll <- summary(l.var)                # Need for the covariance estimates
    covres <- ll$covres
        # The covariance of the residuals
    se <- t( sapply( ll$varresult, function(x) x$coefficients[,'Std. Error'] ) )
        # The standard error
    vcv <- kronecker( covres, ll$varresult[[1]]$cov.unscaled )
        # The variance-covariance matrix of the estimates.  The sqrt of the diag
        # of this will be se.  We can do this with only the first regression in ll
        # because the data is the same for all of the (as it is a VAR).
    g.h <- grad.mu( coeff[,-(n.v+1)], coeff[,n.v+1] )
    se.mu <- sqrt( diag( t(g.h) %*% vcv %*% g.h ) )
        # The se of the means
  }else{
    covres <- cov( sapply( sim.VAR$varresult, function(x) x$residuals ) )
        # The covariance of the residuals
  }

  # Create top and tail
  head.str <- paste0( '\\begin{table}[htbp] \n\t\\centering \n\t\\begin{tabular}{',
                      paste(rep('c', n.col), collapse =''),'}' )
  header.coeff <- paste0( '\t\t\\hline\\hline\n\t\t \t\t & \\multicolumn{', n.v + 1,
                          '}{c}{Regression coefficients} &' )
  header.mean <- if(!add.mean) '' else paste0( 'Mean &' )
  header.covar <- paste0( '\\multicolumn{', n.v, '}{c}{Covariances} \\\\ \n' )
  coeff.names <- paste0( '\t\t \t\t & Const \t & ',
                         paste0( varnames, collapse=' (-1) \t & ' ), ' (-1)',
                         (if(add.mean) '\t & \t & \t' else '\t & \t' ),
                         paste0( varnames, collapse=' \t & ' ), '\\\\ \\hline \n')
  header <- paste0( header.coeff, header.mean, header.covar, coeff.names )
      # The header
  tail.str <- paste0( '\t \\hline \n \t\\end{tabular}',
                      if(is.null(caption)) '' else paste0('\t\t\\caption{',caption,'}\n'),
                      if(is.null(label)) '' else paste0('\t\t\\label{',label,'}\n'),
                      '\n\\end{table}' )
      # The closing part

  # Create the body
  line.i <- function(i){
  # Create the line(s) in the table for variable number i
    nm <- paste0( varnames[i], '\t &')
    coeff.est <- paste0( round( coeff[i,c( n.v+1, 1:n.v )], 2),
                         collapse = '\t &' )
    if( add.mean ) coeff.est <- paste0( coeff.est, ' \t & ', round( mu[i], 2) )
    sig.est <- paste0( sapply( 1:n.v, function(j) if( i >= j ) round( covres[i,j], 3) else ' \t ' ),
                          collapse = '\t &' )
    if(add.se){
      coeff.se <- paste0( round( se[i,c( n.v+1, 1:n.v )], 3) , collapse = ') \t & (' )
      sig.se <- paste0( rep( ' ', n.v ), collapse = ' \t & ' )
      if( add.mean ) coeff.se <- paste0( coeff.se, ') \t & (', round(se.mu[i],3) )
      se.line <- paste0( '\t & (', coeff.se, ') \t & ', sig.se, ' \\\\\n' )
    }else{
      se.line <- NULL
    }
    return( paste0( nm, coeff.est, '\t &', sig.est, '\\\\\n', se.line ) )
  }
  bdy <- paste0( sapply( 1:n.v, line.i  ), collapse ='\n' )

  if(footer & add.se){
    if( !is.null(start) & !is.null(end) ){
      st.dates <- paste0( ' ', year(start), 'Q', quarter(start), ':',
                          year(end), 'Q', quarter(end) )
    }else{
      st.dates <- NULL
    }
    p.val <- 1 - pnorm( 0, diff(mu), se.rmg(l.var) )
        # p.value of R-G<0 hypothesis test
    st.foot <- paste0( '\t\t\\hline\n\t\t\\multicolumn{', n.col, '}{l}{', nrow(l.var$datamat),
                       ' obs', st.dates , '\\hfill $H_0: \\mathbb E R = \\mathbb E G$ vs. $H_0: \\mathbb E R < \\mathbb E G$ p=',
                       round( p.val, 3 ), '}\\\\\n' )
  }else{
    st.foot <- NULL
  }

  out <- paste0( head.str, header, bdy, st.foot, tail.str )

  if(!is.null(file)) cat(out, file=file)

  return(out)
}

grad.rmg <- function( B, A ){
# Returns the numerical derivative of the mean w.r.t the coefficients of A and B
  nn <- nrow(B)
  mm <- length(B) + length(A)
  rmg <- diff( solve( diag(nn) -B, A ) )
  # The mean
  inc <- 1e-06
  m.rmg.inc <- rep( 0, mm )
  # The incremented value of mu
  for( i in 1:length(A) ){
    A.cpy <- A
    A.cpy[i] <- A.cpy[i] + inc
    rmg.new <- diff( solve( diag(nn) - B, A.cpy ) )
    m.rmg.inc[i] <- rmg.new
  }

  for( i in 1:length(B) ){
    B.cpy <- B
    B.cpy[i] <- B.cpy[i] + inc
    rmg.new <- diff( solve( diag(nn) - B.cpy, A ) )
    m.rmg.inc[ length(A)+i ] <- rmg.new
  }

  out <- ( rmg - m.rmg.inc ) / inc
      # The output
  reorder <- c( sapply( 1:nn, function(i) c( i + nn * 1:nn, i ) ) )
      # Order should be: coeffs, then const for each equation
  return(out[reorder])
}

se.rmg <- function( l.var ){
# Computes the standard error on the difference in means of R and G

  n.v <- length(l.var$varresult)      # Number of variables
  ll <- summary(l.var)                # Need for the covariance estimates
  coeff <- t( sapply( l.var$varresult, function(x) x$coefficients ) )
      # The VAR coefficients
  covres <- ll$covres
      # The covariance of the residuals
  vcv <- kronecker( covres, ll$varresult[[1]]$cov.unscaled )
      # The variance-covariance matrix of the estimates.
  g.h <- grad.rmg( coeff[,-(n.v+1)], coeff[,n.v+1] )
      # The gradient of R minus G
  return( sqrt( t(g.h) %*% vcv %*% g.h ) )
}

var.rg.summy <- function(cty, start="1960-01-01", end="2017-01-01", caption=TRUE){
# Computes the summary of the VAR for R & G
  rg.dta <- rg.read( cty )
  rg.dta <- subset( rg.dta, date >= start )
  rg.dta <- subset( rg.dta, date <= end )
  rg.dta <- rg.dta[!apply(rg.dta,1,function(x)any(is.na(x))),]
  start <- min(rg.dta$date)
  end <- max(rg.dta$date)
      # Load and trim the data
  est <- var.rg.est(rg.dta)
  var.table(est$VAR,
            file=paste0('~/Dropbox/2017/research/debtLimits/paper/',cty,'var_tab.tex'),
            add.mean = TRUE, varnames = c('Nom. GDP', 'Interest rate'),
            caption=if(caption) paste0('VAR results for ',cty) else NULL,
            label=paste0('tab:',cty,'_var'),
            start=start, end=end)
}

var.rg.rest <- function(cty, start="1960-01-01", end="2017-01-01"){
# Computes the restricted VAR for R & G
  rg.dta <- rg.read( cty )
  rg.dta <- subset( rg.dta, date >= start )
  rg.dta <- subset( rg.dta, date <= end )
  rg.dta <- rg.dta[!apply(rg.dta,1,function(x)any(is.na(x))),]
  start <- min(rg.dta$date)
  end <- max(rg.dta$date)
      # Load and trim the data
  est <- var.rest.rg.est(rg.dta)
  return(est)
}

var.rg.rest.theta <- function(cty, start="1960-01-01", end="2017-01-01",
                              theta.range=c(0,4), n.theta=20 ){
# Computes the p-values of the restricted VAR for R & G for a range of theta.
  rg.dta <- rg.read( cty )
  rg.dta <- subset( rg.dta, date >= start )
  rg.dta <- subset( rg.dta, date <= end )
  rg.dta <- rg.dta[!apply(rg.dta,1,function(x)any(is.na(x))),]
  start <- min(rg.dta$date)
  end <- max(rg.dta$date)
      # Load and trim the data
  v.theta <- seq( theta.range[1], theta.range[2], length.out = n.theta )
      # The vector of thetas
  out <- cbind( v.theta,
                sapply( v.theta,
                        function(theta) 1 -
                          pchisq( var.rest.rg.est(rg.dta, theta)$l.ratio, 1 ) ) )
      # Construct the output
  return(out)
}
