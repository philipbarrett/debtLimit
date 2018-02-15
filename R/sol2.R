#####################################################################
# sol2.R
#
# Trying a simpler way to solve the model
# 07nov2017
# Philip Barrett, Washington DC
#####################################################################

target.i <- function( x, guess, params, i ){
  d <- guess[,2]; p.i <- x[1] ; d[i] <- x[2]
  if( p.i > 1 ){
    ret <- 100000 * c( (p.i-1)^2, 2*(p.i-1) )
  }
  else if( p.i < 0 ){
    ret <- zed_2_ana_0_i( 0, d, params, i-1 )[1] - 100000 * c( p.i^2, 2*p.i )
  }else{
    ret <- zed_2_ana_0_i( p.i, d, params, i-1 ) - c(p.i,1)
  }
  if(d[i] < 0){
    ret <- ret - 1000 * d[i] ^ 2
  }
  return( ret )
}

sol.p.0 <- function( p.guess, d, params ){
# Objective function for finding d s.t. p=0 everywhere
  ret <- zed_2_ana_0(p.guess, d, params)[,1]
  ret[p.guess<0] <- (zed_2_ana_0( 0*p.guess, d, params )[2,1] - 10000 * p.guess[p.guess<0]^2 )[p.guess<0]
  ret[p.guess>0] <- (10000 * (p.guess[p.guess<0]-1)^2)[p.guess>0]
  return(ret)
}

sol.2 <- function(params, init.guess=NULL, maxit=60, gain=.5, tol=1e-05, d.all.init=130, print.level=2 ){
# A simple solution
  n.states <- length(params$R)
  if( is.null(init.guess)){
    d.sol <- nleqslv( rep(d.all.init,n.states), sol.p.0, p.guess=rep(0,n.states), params=params, control=list(allowSingular=TRUE, maxit=400) )
    d.guess <- d.sol$x
    # d.guess <- sol.nonstoch(params)
    # d.guess[d.guess==Inf] <- max(d.guess[d.guess<Inf])
    guess <- cbind( rep(0,n.states), d.guess)
  }else{
    guess <- init.guess
  }
  n.p.alt <- 125 # 250
  err.i <- err <- new.guess <- guess

  for( it in 1:maxit ){
    if( print.level>=2) message('*** iteration ', it, ' ***')
    for( i in 1:n.states){
      # if(it==12 && i==17){
      #   browser()
      # }
      t.i <- target.i( guess[i,], guess=guess, params=params, i=i )
          # Try to see if the solution is finite
      if( all(t.i < Inf) & !any(is.na(t.i) ) ){
        sol <- nleqslv(guess[i,], target.i, guess=guess, params=params, i=i,
                       control = list(allowSingular=TRUE))
        trial <-if(sol$termcd != 1) TRUE else FALSE
        p.alt <- c(t(cbind(seq(guess[i,1], 1, length.out=n.p.alt),
                           seq(guess[i,1], 0, length.out=n.p.alt))))
        counter <- 1
        min.err <- Inf
        new.guess[i,] <- sol$x
        err.i[i,] <- sol$fvec
      }else{
        trial <- TRUE
      }
      min.err <- Inf
      while(trial){
        if(counter > 2*n.p.alt ){
          if( print.level>=2) message('i = ', i, ':  Exhausted p.alt. Min error = ', min.err )
          # browser()
          trial <- FALSE
        }
        t.i <- target.i( guess[i,], guess=guess, params=params, i=i )
        if( all(t.i < Inf) & !any(is.na(t.i) ) ){
          sol <- nleqslv(c( p.alt[counter], guess[i,2]), target.i, guess=guess, params=params, i=i,
                         control = list(allowSingular=TRUE))
          if( max(abs(sol$fvec)) < min.err ) new.guess[i,] <- sol$x
          min.err <- min( min.err, max(abs(sol$fvec)))
          if(sol$termcd == 1 | min.err < 1e-06 ){
            trial <- FALSE
            err.i[i,] <- sol$fvec
          }
        }
        counter <-  counter + 1
      }
      if( print.level>=2) message( 'i = ', i, ', err.i[i,] = ', signif(max(abs(err.i[i,])), 4) )
    }
    diff <- max(abs(new.guess-guess))
    guess <- gain * new.guess + (1-gain) * guess
    if( print.level>=1) message('it = ', it, ', max diff = ', signif(diff, 4))
    if( print.level>=2) message('mean diff = (', signif(mean(new.guess[,1]-guess[,1]), 4), ', ',
                             signif(mean(new.guess[,2]-guess[,2]), 4), ') \n')
    if(diff<tol) break()
  }
  err <- zed_2_ana_0( guess[,1], guess[,2], params ) - cbind(guess[,1],1)

  return(list( p=guess[,1], d=guess[,2], err.i=err.i, diff=diff, err=err, params=params ))
}
