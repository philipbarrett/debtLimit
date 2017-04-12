#####################################################################
# plotting.R
#
# Contains the plotting functions
# 22feb2017
# Philip Barrett, Washington DC
#####################################################################

plot.surp <- function(params, non.stoch=TRUE, x.lim=c(0,200), ... ){
# Plots the surplus function
  d <- seq( from=x.lim[1], to=x.lim[2], by = 1 )
      # X values
  surps <- sapply( params$G, function(G)
    sapply(d, surp, coeff=params$v.s.coeff, G=G, tri=params$tri) )
      # Y-values
  y.max <- max(surps)
  d.max <- d[which.max(if(is.null(dim(surps))) surps else surps[,1])]
  if(d.max==0) d.max <- x.lim[2]
  y.min <- if( min(surps) < 0 ) ( if( params$tri ) min(surps) else min(surps[d<d.max]) ) else 0
  y.lim <- c(y.min,y.max)
      # Y range
  plot( x.lim, y.lim, type='n', xlab='Debt', ylab='Surplus', ...  )
  for( i in 1:ncol(surps) ){
    lines( d, surps[,i], lwd=2, col=i )
    if( non.stoch ) abline( 0, params$R[i] - params$G[i], lwd=.5 )
  }
  legend( 'topright', paste0( 'G = ', round( params$G, 3 ) ), lwd=2, col=1:i, bty='n' )
  abline(h=0)
  if( non.stoch ){
    abline( v=sol.nonstoch(params), lwd=.5, col=1:length(params$R) )
  }
}

plot.z <- function( p, d, params, An=c(0), Bn=c(0), Cn=c(0), def=matrix(0), ... ){
# Plots the zed function for all the values of i
  n <- length(params$R)
  n.x <- ceiling( sqrt(n) )
  n.y <- ceiling( n / n.x )
  par(mfrow=c(n.y,n.x))
  global.apx <- p_init_d( params, p, d, An, Bn, Cn, def )
  for( i in 1:n )
    plot.z.i( p, d, params, i, An, Bn, Cn, def, global.apx[i], ... )
  par(mfrow=c(1,1))
}

plot.z.i <- function( p, d, params, i, An, Bn, Cn, def, global, ... ){
# Plots the z function vs. p[i] in state i
  p.seq <- c( seq(0,1e-3,by=1e-5), seq(1e-3,1e-2,by=1e-4), seq(1e-2,1e-1,by=1e-3), seq(1e-1,1,by=1e-2) )
      # The x-values
  y.vals <- t( sapply( p.seq, function(p.i){
    this.p <- p
    this.p[i] <- p.i
    return( zed_2(this.p, d, params, An, Bn, Cn, def )[i,] )
  } ) )
      # The y values
  # if( !exists('ylim')) ylim <- c(0,1)
  plot( p.seq, y.vals[,1], lwd=2, xlab='p', ylab='z', main=paste0( 'i = ', i ), type='l', ... )
  abline( 0, 1, lty=2 )
  abline( v=p[i], lty=2 )
  abline( h=p[i], lty=2 )
      # The level
  abline( v=global, lty=2, col='blue' )
      # The global minimum
  par(new = TRUE)
  y.range <- c( max( min(c(0, y.vals[,2]), na.rm = TRUE), -10 ),
                min( max(y.vals[,2], na.rm = TRUE ), 10 ) )
      # The y limits for the derivative
  plot( p.seq, y.vals[,2], type = "l", col='red', axes = FALSE, bty = "n",
        xlab = "", ylab = "", ylim=y.range )
      # Plot the derivative
  axis(side=4, at = pretty(y.range))
  # mtext("z.p", side=4, line=3)
  abline( h=1, lty=2, col='red' )
      # Axis labelling
}

plot.z.d <- function( p, d, params, An=c(0), Bn=c(0), Cn=c(0), def=matrix(0),
                    d.range=c(0,200), global=NULL, ... ){
# Plots the zed function for all the values of i
  n <- length(params$R)
  n.x <- ceiling( sqrt(n) )
  n.y <- ceiling( n / n.x )
  par(mfrow=c(n.y,n.x))
  # global.apx <- p_init_d( params, p, d, An, Bn, Cn, def )
  for( i in 1:n )
    plot.z.d.i( p, d, params, i, An, Bn, Cn, def, d.range, global, ... )
  par(mfrow=c(1,1))
}

plot.z.d.i <- function( p, d, params, i, An, Bn, Cn, def, d.range=c(0,200), global, ... ){
# Plots the z function vs. d[i] in state i
  d.seq <- seq(d.range[1], d.range[2],length.out=1000)
      # The x-values
  y.vals <- t( sapply( d.seq, function(d.i){
    this.d <- d
    this.d[i] <- d.i
    return( zed_2(p, this.d, params, An, Bn, Cn, def )[i,] )
  } ) )
  # The y values
  plot( d.seq, y.vals[,1], lwd=2, xlab='d', ylab='z', main=paste0( 'i = ', i ), type='l', ... )
  # abline( 0, 1, lty=2 )
  abline( v=d[i], lty=2 )
  abline( h=p[i], lty=2 )
      # The level
  abline( v=global, lty=2, col='blue' )
      # The global minimum
  par(new = TRUE)
  y.range <- c( max( min(c(0, y.vals[,2]), na.rm = TRUE), -10 ),
                min( max(y.vals[,2], na.rm = TRUE ), 10 ) )
      # The y limits for the derivative
  plot( d.seq, y.vals[,2], type = "l", col='red', axes = FALSE, bty = "n",
        xlab = "", ylab = "", ylim=y.range )
      # Plot the derivative
  axis(side=4, at = pretty(y.range))
      # mtext("z.p", side=4, line=3)
  abline( h=1, lty=2, col='red' )
      # Axis labelling
}

plot.q <- function( p, params, An=NULL, def=NULL ){
# Plots the q function for all the values of i
  n <- length(params$R)
  n.x <- ceiling( sqrt(n) )
  n.y <- ceiling( n / n.x )
  par(mfrow=c(n.y,n.x))
  for( i in 1:n )
    plot.q.i( p, params, i, An, def )
  par(mfrow=c(1,1))
}

plot.q.i <- function( p, params, i, An=NULL, def=NULL ){
# Plots the q function vs. p[i] in state i
  p.seq <- seq(0,1,by=.001)
      # The x-values
  n <- length(params$R)
      # Number of states
  if( is.null(An) ) An <- c(0)
  if( is.null(def) ) def <- matrix(0,1,1)
      # Fill in An and def if required
  y.vals <- sapply( p.seq, function(p.i){
    this.p <- p
    this.p[i] <- p.i
    return( q_fn(params$R, this.p, params$trans, params$lambda, params$phi, n,
                 params$cont.type, params$G, An, def ) )
  } )
      # The y values
  plot( range(p.seq), range(y.vals), xlab='p', ylab='q', main=paste0( 'i = ', i ), type='n', ylim=c(0,1) )
  for( j in 1:n ){
    if( j == i ){
      lines( p.seq, y.vals[j,], lwd=2, col='red' )
    }else{
      lines( p.seq, y.vals[j,], lwd=.5, col='black' )
    }
  }
  # abline( h = ( 1 - params$lambda ) * params$phi / params$R[i], lty=2, lwd=1 )
  abline( h = 1 / params$R[i], lty=2, lwd=1 )
  abline( v = p[i], lty=2, lwd=1 )
  legend( 'topright', paste0( 'R=', params$R[i]-1, '\nG=', params$G[i]-1), bty='n' )
}


plot.sol <- function( sol ){
# Plots the outer solution

  par(mfrow=c(2,2))
      # Multiple plots
  n <- length(params$R)

  # Plot Q
  plot( range(sol$d.grid), range(sol$Q), type='n', xlab='Debt value', ylab='',
        main='Debt price' )
  abline(v=sol$d.bar, lty=2, lwd=2, col=1:n )
  abline(v=sol$d.grid, lty=1, lwd=.25 )
  abline(h=c(0,1))
  for( i in 1:n ) lines( sol$d.grid, sol$Q[i,], col=i, lwd=2 )

  # Plot QE
  plot( range(sol$d.grid), range(sol$QE), type='n', xlab='Debt value', ylab='',
        main='Expected continuation debt price' )
  abline(v=sol$d.bar, lty=2, lwd=2, col=1:n )
  abline(v=sol$d.grid, lty=1, lwd=.25 )
  abline(h=c(0,1))
  for( i in 1:n ) lines( sol$d.grid, sol$QE[i,], col=i, lwd=2 )

  # Plot P
  plot( range(sol$d.grid), range(sol$P), type='n', xlab='Debt value', ylab='',
        main='Default probability' )
  abline(v=sol$d.bar, lty=2, lwd=2, col=1:n )
  abline(v=sol$d.grid, lty=1, lwd=.25 )
  abline(h=c(0,1))
  for( i in 1:n ) lines( sol$d.grid, sol$P[i,], col=i, lwd=2 )

  # Plot D.prime
  plot( range(sol$d.grid), range(sol$D.prime), type='n', xlab='Debt value', ylab='',
        main='Expected continuation debt' )
  abline(v=sol$d.bar, lty=2, lwd=2, col=1:n )
  abline(h=sol$d.bar, lty=2, col=1:n )
  abline(v=sol$d.grid, lty=1, lwd=.25 )
  abline(0,1,lty=2)
  for( i in 1:n ) lines( sol$d.grid, sol$D.prime[i,], col=i, lwd=2 )

  par(mfrow=c(1,1))
}

plot.err <- function(sol.err){
# Plots the errors from a solution
  par(mfrow=c(2,2))
  plot.err.i( sol.err$P, main='Default probability' )
  plot.err.i( sol.err$Q, main='Debt price' )
  plot.err.i( sol.err$D.prime, main='Ave. continuation debt' )
  tb <- sapply( c('P','Q','D.prime'), function(x) c('Ave err'=mean(sol.err[[x]]),
                                        'Ave abs err'=mean(abs(sol.err[[x]])) ) )
      # Table of ave and ave abs errors
  d.bar.err <- sol.err$d.bar.range - sol.err$d.bar
      # The errors
  nn <- length(params$R)
      # Number of states
  plot( c(1,nn), range(d.bar.err), type='n', xlab='State #', ylab='Upper-lower bound range',
        xaxt='n', main='Inner vs. outer bounds check')
  axis(1,at=1:nn,labels=format(1:nn,digits = 1))
  abline(h=0,lwd=.5)
  segments(1:nn,d.bar.err[,1],1:nn,d.bar.err[,2])
  epsilon <- 0.02
  segments(1:nn-epsilon,d.bar.err[,1],1:nn+epsilon,d.bar.err[,1])
  segments(1:nn-epsilon,d.bar.err[,2],1:nn+epsilon,d.bar.err[,2])
      # The consistency check for debt limits
  par(mfrow=c(1,1))
  return(tb)
}

plot.err.i <- function(err, ... ){
# Plots an inidividual error distribution
  plot( density( abs(c(err)), from=0 ), xlab='Absolute error', ylab='Density', lwd=2, ... )
  abline( v = quantile( abs(c(err)), c(.5, .75, .9, .95, .99)), lty=2 )
}
