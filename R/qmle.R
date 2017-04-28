#####################################################################
# qmle.R
#
# Contains the functions for quasi maximum likelihood of the
# Markov-switching VAR
# 18apr2017
# Philip Barrett, Washington DC
#####################################################################

init.var <- function( date, gr, n.state=3, seed=42 ){
# Creates an initial guess with constant variance and persistent, but varying
# means based on clustering

  set.seed(seed)
      # Fix randomization
  nn <- ncol(gr)
      # Number of columns
  if( is.null( colnames(gr) ) ) colnames(gr) <- c('g','r')
  # Set column names
  cl <- kmeans( gr, n.state )
  # Compute the clustered means
  dta.var <- gr - cl$centers[cl$cluster,]
  # Remove the appropriate clustered mean
  gr.var <- VAR( dta.var )
  # The VAR with means removed
  A <- t(sapply( colnames(gr), function(x) gr.var$varresult[[x]]$coefficient[-(nn+1)] ))
  # Extract the autoregressive coefficients
  a <- ( diag( nn ) - A ) %*% t(cl$centers)
  # Compute the constant term
  sigma <- var( sapply( colnames(gr), function(x) gr.var$varresult[[x]]$residuals ) )
  # The conditional variance
  par <- rbind( a, matrix( c( A, sigma[lower.tri(sigma,TRUE)] ),
                           nrow=nn^2+.5*nn*(nn+1), ncol=n.state ) )
  # The VAR parameter matrix
  freq <- table( cl$cluster[-(nrow(dta.gr))], cl$cluster[-1] )
  P <- freq / apply( freq, 1, sum )
  # Create a guess of the means and transition probabilities
  p.0 <- rep(0,n.state)
  p.0[cl$cluster[1]] <- 1
  # Initialize the first period
  g.par <- c( p.0[-nn], P[,-nn], a, A, sigma[lower.tri(sigma,TRUE)] )
  # The grand parameter vector
  return( list( mu=t(cl$centers), mu=t(cl$centers), a=a, A=A, sigma=sigma,
                par=par, g.par=g.par, p.0=p.0, P=P, cluster=cl$cluster ) )
}


ql.sw.mu <- function(g.par, dta, n.state, nn, par.out=FALSE ){
# Quasi-likelihood function for Markov mean-swtiching
# Takes and in input a grand vector of parameters in the order:
#   p.0[1:(nn-1)], P[,1:(nn-1)], a, A, sigma (lower tri)

  p.0.idx <- 1:(n.state-1)
  P.idx <- n.state-1 + 1:( n.state * (n.state-1) )
  a.idx <- ( n.state + 1 ) * (n.state-1) + 1:(nn*n.state)
  A.idx <- ( n.state + 1 ) * (n.state-1) + nn * n.state + 1:(nn^2)
  sigma.idx <- ( n.state + 1 ) * (n.state-1) + nn * n.state + nn^2 + 1:(.5*(nn+1)*nn)
      # Inidices of the parameters
  p.0.x <- c( g.par[p.0.idx], 1-sum(g.par[p.0.idx]) )
      # Initial p.0
  p.0 <- pmin( pmax( p.0.x, 0 ), 1 )
      # Bounded p.0
  P.x <- matrix( g.par[P.idx], n.state, n.state - 1 )
  P.x <- cbind(P.x, 1 - apply( P.x, 1, sum ) )
      # Initial P
  P <- pmin( pmax( P.x, 0 ), 1 )
      # Bounded P
  a <- matrix( g.par[a.idx], nn, n.state )
  A <- matrix( g.par[A.idx], nn, nn )
      # Constant and persistence matrices
  sigma <- matrix( 0, nn, nn )
  counter <- 0
  for( i in 1:nn ){
    for( j in i:nn ){
      counter <- counter+1
      sigma[i,j] <- sigma[j,i] <- g.par[sigma.idx][counter]
    }
  }
      # Extract the prameters
  penalty <- sum( 100000 * p.0.x * ( p.0.x - 1 ) * ( p.0.x != p.0 ) ) +
    sum( 100000 * P.x * ( P.x - 1 ) * ( P.x != P ) )
  # Penalty function for violating bounds on p
  m.par <- rbind( a, matrix( c( g.par[c(A.idx, sigma.idx)] ),
                             nrow=nn^2+.5*nn*(nn+1), ncol=n.state ) )
  # Arrange non-Markov parameters into a matrix
  l <- msw_var_lhood( t(as.matrix(dta)), p.0, P, m.par ) + penalty
      # The likelihood
  if( par.out ){
    probs <- msw_var_lhood_p( t(as.matrix(dta.gr)), p.0, P, m.par )
        # The filtering and perdiction probabilities
    modal.state <- apply( probs[1:n.state,], 2, which.max )
    return( list( p.0=p.0, P=P, mu=solve(diag(nn)-A,a), a=a, A=A, sigma=sigma,
                  m.par=m.par, penalty=penalty, l=l, probs=probs, modal.state=modal.state ) )
  }
      # Return the parameters
  return(l)
}


msw.var.lhood <- function( Y, p0, P, par, print.level=0 ){
# Returns the forecast and filtering/ probabilities for the quasi-likelihood
# for the Markov switching model

  ll <- 0
  nn <- length(p0)
  mm <- ncol(Y)
      # Problem dimensions
  P.t <- t(P)
  m.p.filter <- matrix( 0, nn, mm )
  m.p.pred <- matrix( 0, nn, mm )
  v.l.inc <- rep(0,mm)
      # Storage
  p.filter <- p0
  p.pred <- P.t %*% p0
  m.p.filter[,1] <- p.filter
  m.p.pred[,1] <- p.pred
      # Initialization
  for( i in 2:mm ){
    l.vec <- one_step_var_lhood( Y[,c(i-1,i)], par, print.level )
    p.vec <- exp( - l.vec )
        # Negative log likelihood and conditional likelihood vector
    l.inc <- sum( p.vec * p.pred )
    v.l.inc[i] <- l.inc
    ll <- ll - log( l.inc ) / ( mm - 1 )
        # Increment likelihood
    p.filter <- p.pred * p.vec / l.inc
    p.pred <- P.t %*% p.filter
        # update MArkov state probabilities
    m.p.filter[,i] <- p.filter
    m.p.pred[,i] <- p.pred
  }
  return( list( ll=ll, m.p.filter=m.p.filter, m.p.pred=m.p.pred, v.l.inc=v.l.inc ) )
}

ms.est.random.starts <- function( date, gr, n.state=3, seed=42, n.start=100, control=NULL ){
# Estimates the model using a sequence of random starts

  set.seed(seed)

  if(is.null(control))
    control <- list( abstol = 1e-12, reltol = 1e-12, trace = 1, maxit=2000 )
      # List of controls
  nn <- ncol(gr)
      # The number of variables
  var.gr <- var( gr )
  mu.gr <- apply( gr, 2, mean )
      # Summary stats of the data
  l <- Inf
      # The likelihood
  par <- init.var( date, gr, n.state, seed )$g.par
      # The initial grand parameter vector

  for( i in 1:n.start ){

    p.0.init <- runif( n.state-1, 0, 1 )
    P.init <- runif( n.state*(n.state-1), 0, 1 )
    a.init <- rmvnorm( n.state, mu.gr, var.gr )
        # Randomize initial guess of parameters
    par[1:((n.state+1)*(n.state-1)+length(a.init))] <- c( p.0.init, P.init, a.init )
        # Replace parameters

    while( is.nan(ql.sw.mu( par, gr, n.state, nn ) ) ){
      p.0.init <- runif( n.state-1, 0, 1 )
      P.init <- runif( n.state*(n.state-1), 0, 1 )
      a.init <- rmvnorm( n.state, mu.gr, var.gr )
      par[1:((n.state+1)*(n.state-1)+length(a.init))] <- c( p.0.init, P.init, a.init )
          # Redo things if we have a failed start
    }

    opt.cand <- optim( par, ql.sw.mu, dta=gr, n.state=n.state, nn=nn,
                       method='BFGS', control=control )
        # The optimization problem
    if( opt.cand$value < l ){
      opt <- opt.cand
      l <- opt$value
    }
        # Update the candidate solution if an improvement
  }
  return(ql.sw.mu( opt$par, dta.gr, n.state, nn, TRUE ))
}




