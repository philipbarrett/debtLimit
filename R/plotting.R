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
  y.min <- if( min(surps) < 0 ) ( if( params$tri ) min(surps) else -.25 * y.max ) else 0
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
    abline( v=sol.nonstoch(params), lwd=.5 )
  }
}

plot.z <- function( p, d, params, An=c(0), Bn=c(0), Cn=c(0), def=matrix(0) ){
# Plots the zed function for all the values of i
  n <- length(params$R)
  n.x <- ceiling( sqrt(n) )
  n.y <- ceiling( n / n.x )
  par(mfrow=c(n.y,n.x))
  global.apx <- p.init.d( params, p, d, An, Bn, Cn, def )
  for( i in 1:n )
    plot.z.i( p, d, params, i, An, Bn, Cn, def, global.apx[i] )
  par(mfrow=c(1,1))
}

plot.z.i <- function( p, d, params, i, An, Bn, Cn, def, global=NULL ){
# Plots the z function vs. p[i] in state i
  p.seq <- seq(0,1,by=.001)
      # The x-values
  y.vals <- t( sapply( p.seq, function(p.i){
    this.p <- p
    this.p[i] <- p.i
    return( zed_2(this.p, d, params, An, Bn, Cn, def )[i,] )
  } ) )
      # The y values
  plot( p.seq, y.vals[,1], lwd=2, xlab='p', ylab='z', main=paste0( 'i = ', i ), type='l', ylim=c(0,1) )
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
