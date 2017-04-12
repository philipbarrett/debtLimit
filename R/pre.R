#####################################################################
# pre.R
#
# Contains the functions for presolution work (initial guesses)
# 22feb2017
# Philip Barrett, Washington DC
#####################################################################

sol.nonstoch <- function(params){
# Computes the nonstochastic solution
  if(!params$tri)
    return(sol.nonstoch.poly(params))
  return(sol.nonstoch.tri(params))
}

sol.nonstoch.poly <- function(params){
# Polynomial surlpus function
  n <- length(params$R)
  m <- length(params$v.s.coeff)
  out <- rep(0,n)
      # Initialize the output
  for( i in 1:n ){
    z <- c( params$v.s.coeff[2] + ( params$G[i] - 1 ) * 100 * params$v.s.coeff[m],
            params$v.s.coeff[3:(m-1)] )
        # Coefficients of the polynomial function for surpluses (add the
        # function of G as a constant)
    z[2] <- z[2] - ( params$R[i] - params$G[i] ) * params$v.s.coeff[1]
        # Subtract R-G*scaling from the linear ccoefficient
    rt <- polyroot(z)
        # The roots
    re.rt <- Re(rt[ abs(Im(rt)) < 1e-12 ])
        # The real roots
    z.d <- z[-1] * 1:(m-3)
        # Coefficients of the derivative
    re.rt.d <- 1 / params$v.s.coeff[1] * z.d %*% sapply(re.rt, function(x) x ^ (1:(m-3)) )
        # The derivative at the real roots
    out[i] <- max( re.rt[re.rt.d<0] ) * params$v.s.coeff[1]
  }
  return(out)
}

sol.nonstoch.tri <- function(params){
# Triangular surplus function
  n <- length(params$R)
  m <- length(params$v.s.coeff)
  out <- rep(0,n)
      # Initialize the output
  for( i in 1:n ){
    if( params$v.s.coeff[5] + params$v.s.coeff[6]*( params$G[i]-1 )*100 > 0 &
        params$R[i] < params$G[i] ){
        # Case I: Any arbitrarily large debt level can be maintained
      out[i] <- Inf
    }else{
      if( params$v.s.coeff[5] + params$v.s.coeff[6]*( params$G[i]-1 )*100 < 0 &
          params$R[i] > params$G[i] ){
          # Case II: No positive debt level can be maintained
        out[i] <- 0
      }else{
          # Case III: A finite, positive debt level can possibly be maintained
        fn <- function(d) surp_tri( d, params$v.s.coeff, params$G[i] ) -
                              ( params$R[i] - params$G[i] ) * d
        jac <- function(d) d_surp_tri( d, params$v.s.coeff, params$G[i] ) -
                              ( params$R[i] - params$G[i] )
            # Set up the funciton and jacobian for the root of s(d) = (R-G)*d
        s.upper <- params$G[i] * params$v.s.coeff[6] + params$v.s.coeff[5]
        d.0 <- s.upper / ( params$R[i] - params$G[i] )
            # Initial guess
        sol <- nleqslv( d.0, fn, jac )
            # Find the root
        if( abs(sol$fvec) > 1e-06 ){
          if( params$v.s.coeff[3] * ( params$R[i] - params$G[i] ) >
              params$v.s.coeff[5] + params$v.s.coeff[6]*( params$G[i]-1 )*100 ){
          # Case IIIa: If cannot find a solution, check to see if R-G line is
          # above s(b) at the start of the upper flat part.
            out[i] <- 0
          }else{
            warning("Nonstochastic triangular surplus function failng to solve. \nErr=",
                    round(abs(sol$fvec),5) )
              # This should not happen!!
          }
        }else{
          if( sol$x<0 ){
            out[i] <- if( params$R[i] - params$G[i] <= 0 ) Inf else 0
          }else{
            out[i] <- sol$x
          }
        }
      }
    }
  }
  return(out)
}

q.fn <- function( p, params, An, def ){
# Wrapper for q_fn
  return( q_fn(params$R, p, params$trans, params$lambda, params$phi, length(params$R),
               params$cont.type, params$G, An, def ) )
}

max.d.s.lin <- function( params, A, B ){
# Solves max_d s.t. s(d) = A + Bd for each state
  if(!params$tri)
    return(max.d.s.lin.poly(params,A,B))
  return(max.d.s.lin.tri(params,A,B))
}

max.d.s.lin.poly <- function( params, A, B ){
# Polynomial case
  n <- length(params$R)
  m <- length(params$v.s.coeff)
  out <- rep(0,n)
      # Initialize the output
  if(length(A)==1) A <- rep(A,n)
  if(length(B)==1) B <- rep(B,n)
      # Allow for both scalar and vector inputs
  for( i in 1:n ){
    z <- c( params$v.s.coeff[2] + ( params$G[i] - 1 ) * 100 * params$v.s.coeff[m],
            params$v.s.coeff[3:(m-1)] )
        # Coefficients of the polynomial function for surpluses (add the
        # function of G as a constant)
    z[1] <- z[1] - A[i]
    z[2] <- z[2] - B[i] * params$v.s.coeff[1]
        # Subtract A & B from the appropriate coefficients (including some scaling)
    rt <- polyroot(z)
        # The roots
    re.rt <- Re(rt[ abs(Im(rt)) < 1e-12 ])
        # The real roots
    # z.d <- z[-1] * 1:(m-3)
    # # Coefficients of the derivative
    # re.rt.d <- 1 / params$v.s.coeff[1] * z.d %*% sapply(re.rt, function(x) x ^ (1:(m-3)) )
    # # The derivative at the real roots
    out[i] <- max( re.rt ) * params$v.s.coeff[1]
  }
  return(out)
}

max.d.s.lin.tri <- function( params, A, B ){
  warning("Triangular distribution linear intersect incomplete")
  return(0)
}

d.init.p <- function( params, An=c(0), Cn=c(0), def=matrix(0), p=NULL, x.sd=-3 ){
# Computes the initial value of the vector of highest d such that a shock x.sd
# standard deviations below the mean is required to produce stationary debt,
# using the equation:
#     s + x.sd * s.sd + d = d * ((1-lambda)+lambda*qd) / ((1+g)*q)
  if(is.null(p)) p <- 0 * params$R
      # Set the initial probabilities to zero if missing.
  q <- q.fn( p, params, An, def )
      # The price in the current period
  qd <- if( params$cont.type != 'fix' ) q else Cn
      # The continuation price
  if( params$d.tri ) x.sd <- max( x.sd, - sqrt(6) +1e-03 )
      # If using the triangular shock distribution, cannot be more than -sqrt(6)
      # sds below the mean
  A <- - x.sd * params$surp.sd
      # The intercept
  B <- ( 1 - params$lambda + params$lambda * qd ) / ( params$G * q ) - 1
      # The slope
  return( max.d.s.lin( params, A, B ) )
}

p.d.init <- function( params, qe=c(0), qd=c(0), def=matrix(0), p=NULL, x.sd=-2 ){
# Computes an initial guess for (p,d) satisfying z=1 and z' > 1

  if(is.null(p)) p <- 0 * params$R
      # Set the initial probabilities to zero if missing.
  fn <- function(p){
    # The function for which we want to find the root
    d <- d.init.p( params, qe, qd, def, p, x.sd )
        # Find the corresponding d
    z <- zed( p, d, params, qd, qe, def )
        # The vector of: the implied value of p and the gradient wrt p
    return( z - p )
        # Because we want to find the fixed point
  }
  control <- list( maxit=500, allowSingular =TRUE, ftol=1e-10, xtol=1e-12 )
  sol <- nleqslv( p, fn, control = control )
      # The solution object
  p <- sol$x
  d <- d.init.p( params, qe, qd, def, p, x.sd )
  err <- sol$fvec
      # The solution
  return( list( p=p, d=d, err=err ) )
}

p.init.d <- function( params, p, d, An, Bn, Cn, def ){
# Finds a quick initial quess of p as the minimizer of z-p in each dimension
  p.seq <- seq(0,1,by=.001)
      # The x-values
  n <- length(params$R)
      # The number of states
  z <- sapply( 1:n, function(i)
                sapply( p.seq, function(p.i){
                  this.p <- p
                  this.p[i] <- p.i
                  return( zed(this.p, d, params, An, Cn, def )[i] )
                } ) )
      # Z
  z.p <- sapply( 1:n, function(i)
                sapply( p.seq, function(p.i){
                  this.p <- p
                  this.p[i] <- p.i
                  return( zed_2(this.p, d, params, An, Bn, Cn, def )[i,2] )
                } ) )
      # Derivative
  bo.cand <- apply( z.p, 2, function(x) (abs(x-1)<1e-01) )
      # Matrix of candidates with derivative v. close to unity
  p.out <- sapply( 1:n, function(i)
        p.seq[bo.cand[,i]][which.min((z[,i]-p.seq)[bo.cand[,i]])] )
      # Select the lowest score of z-p from the possible candidates
  z.out <- sapply( 1:n, function(i){
                                this.p <- p
                                this.p[i] <- p.out[i]
                                return( zed(this.p, d, params, An, Cn, def )[i] )
                              } )
      # The corresponding z
  out <- sapply( 1:n, function(i) if( z.out[i] - p.out[i] > z[1,i] ) 0 else p.out[i] )
      # If zero does better, replace
  return( out )
}

d.p.init.wrapper <- function(params, An, Bn, Cn, def ){
# Grand wrapper function to find initial guesses for p and d

  ## Set up ##
  maxit.1 <- 40
  maxit.2 <- 40
  nn <- length(params$R)
  p.fp <- p.guess <- rep( 0, nn )
  bo.step.1 <- TRUE
  x.sd <- x.sd.old <- rep( -3, nn )
  it <- 0
  tol <- 1e-06

  ## Step 1: Choose d to guarantee an interior fixed point ##
  while( it < maxit.1 && any(p.fp %in% c(0,1) ) ){
    d.guess <- d.init.p( params, An, Cn, def, rep(0,nn), x.sd )
        # Non-increasing debt level given a shock x.sd std devs from mean
    if( all(d.guess>0) ){
      p.fp <- p_fp( params, rep(0,nn), d.guess, An, Bn, Cn, def )
          # The fixed point computation
      x.sd[p.fp==1] <- x.sd[p.fp==1] * 1.25
      x.sd[p.fp==0] <- x.sd[p.fp==0] * .9
          # Update x.sd
    }else{
      x.sd[d.guess<0] <- x.sd[d.guess<0] *.5
    }
    it <- it + 1
  }

  if(any(p.fp %in% c(0,1)))
    stop('Cannot find guess of d with interior fixed point for p')

  ## Step 2: Iterate over initial points to find something near the tangency condition ##
  it <- 0
  p.guess.old <- p.guess
  d.guess.old <- d.guess
  diff <- 2 * tol
      # Initialize loop variables
  while( it < maxit.2 && diff > tol ){
    p.guess <- p_init_d( params, p.guess.old, d.guess.old, An, Bn, Cn, def )
    d.guess <- d_init_p( params, p.guess, d.guess.old, An, Bn, Cn, def, 400 )
        # Update guesses of p and d
    diff <- max( abs( c( p.guess - p.guess.old, d.guess - d.guess.old ) ) )
        # The difference
    p.guess.old <- p.guess
    d.guess.old <- d.guess
        # Update guesses
    it <- it + 1
        # Increment counter
  }
  return( cbind( p.guess, d.guess ) )
}


