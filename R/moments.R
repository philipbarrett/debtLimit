#####################################################################
# moments.R
#
# Contains the functions for creating the moments of the data to target
# 09jun2017
# Philip Barrett, Washington DC
#####################################################################

target.moms <- function( cty, breaks=NULL, labs=NULL, bo.dta=TRUE ){
  # Creates the target moments for the country
  if( is.null(breaks) )
    breaks <- as.Date( c( '1960/01/01', '1970/01/01', '1980/01/01',
                          '1990/01/01', '2000/01/01', '2009/01/01', '2018/01/01' ) )
  if( is.null(labs) ) labs <- c('60s', '70s', '80s', '90s', '2000s pre-09', 'post-09')
  # Default breaks and labels
  dta <- read.csv(paste0('./data/',cty, '.csv'))
  dta$date <- as.Date( dta$date, "%m/%d/%Y" )
  # Read in and clean the data
  dta$spd.1yr <- ( dta$Int_1y - dta$rfrA ) / 400
  dta$spd.5yr <- ( dta$Int_5y - dta$rfrA ) / 400
  # The spreads over the current risk free rate
  v.moms <- c('spd.1yr','spd.5yr','cnlb_gdp','pb_gdp')
  dec.ave <- do.call( rbind, by( dta[,v.moms], cut(dta$date,breaks,labs),
                                 function(x) apply( x, 2, mean, na.rm=TRUE ) ) )
  s.var <- c(by( dta[,'pb_gdp'], cut(dta$date,breaks,labs), sd, na.rm=TRUE ) )
  r.cor <- c(by( dta, cut(dta$date,breaks,labs),
                 function(x) cor(x$pb_gdp, x$rfrA, use='complete.obs' ) ) )
  g.cor <- c(by( dta, cut(dta$date,breaks,labs),
                 function(x) cor(x$pb_gdp, x$ngdp_pch, use='complete.obs' ) ) )
  moms <- cbind( dec.ave, s.var, r.cor, g.cor )
  # The moments
  if(bo.dta){
    return(list(moms=moms, dta=dta))
  }
  return( moms )
}

mom.series.plot <- function( cty ){
  # Plots the target moment series
  dta <- read.csv(paste0('./data/',cty, '.csv'))
  dta$date <- as.Date( dta$date, "%m/%d/%Y" )
  # Read in and clean the data
  dta$spd.1yr <- ( dta$Int_1y - dta$rfrA ) / 400
  dta$spd.5yr <- ( dta$Int_5y - dta$rfrA ) / 400
  # The spreads over the current risk free rate
  v.moms <- c('spd.1yr','spd.5yr','cnlb_gdp','pb_gdp')

  par(mfrow=c(2,2))
  for(vbl in v.moms){
    plot( dta$date, dta[[vbl]], type='l', lwd=2, main=vbl, col='blue' )
    abline(h=0)
  }
  par(mfrow=c(1,1))
}

surp.fit <- function( sol.o, interp, params, dta, debt='cnlb_gdp' ){
  # Creates the sequence of surpluses required to fit the specified debt series in
  # the interpreted simulation

  pds <- nrow(dta)                            # Number of periods
  s.init <- dta[[debt]][-pds] * ( params$R / params$G )[interp$s.idx] - dta[[debt]][-1]
  # Initial guess, just from debt levels
  s <- s.init
  # Initialize the solution
  f.sim.i <- function(s.i){
    s[i] <- s.i
    return( sim_core( c(1,interp$s.idx), sol.o$d.bar, sol.o$d.grid, sol.o$P, sol.o$Q,
                      params, dta[[debt]][1], TRUE, c(0,s) )[i+1,'d.prime'] -
              dta[[debt]][1+i] )
  }
  # The error on the debt series
  for( i in 1:length(interp$s.idx) ){
    sol <- nleqslv( s[i], f.sim.i )
    s[i] <- sol$x
  }
  return(s)
}

rf.fit <- function( params, dta, rf.debt=600, rf.level=12, v.surp=NULL,
                    d.init=NULL, sol=NULL, sol.o=NULL, init=NULL, print.level=0, ... ){
  # Fits the parameters for the shifts and the lower pat of the surplus function
  # based on a risk-free model.

  ### 1. Solve the risk-free model   ###
  params$v.s.coeff <- c( 0, rf.debt/2, rf.debt, rf.level, rf.level )
  n.X <- length(params$R)
  params$s.shift <- rep(0,n.X)
  # A surplus function that puts no restriction on debt
  if( is.null(init) ) init <- c( 0, rf.debt/4, 0, rep(0,n.X-1) )
  # plot.z( sol$p, sol$d, params, xlim=c(0,max(sol$p)*1.2), ylim=c(0,max(sol$p)*1.2) )
  # plot.sol( sol.o )
  # Plots
  if(is.null(sol)) sol <- sol.wrapper(params, cbind(0,rep(rf.debt,n.X)), plot.on = FALSE )
      # Solve for the (approximate) debt levels

  ### 2. Create the surplus series  ###
  if( is.null(v.surp) ){
    if(is.null(sol.o)) sol.o <- outer.wrapper( sol, params, print_level = 1 )
    # Outer solution and errors
    v.surp <- surp.fit( sol.o, interp, params, dta, ... )
    # Create the first pass at the surplus series
  }else{
    d.max <- 10000
    sol.o <- list( d.bar=sol$d,
                   d.grid=c(0,max(sol$d)), P=matrix(1,ncol=2,nrow=length(params$R)),
                   Q=matrix( 1/params$R, ncol=2, nrow=length(params$R) ) )
  }

  if( is.null(d.init) ) d.init <- dta$cnlb_gdp_lag[1]

  ### 3. Create the error function   ###
  params.cpy <- params
  s.idx.p <- sapply(1:n.X, function(i) sum( interp$s.idx == i ) ) / nrow(interp$s.idx)
      # Probability disctibution of states
  err.fn <- function(x, s.dev=FALSE, aad=FALSE ){
    # The error function
    params.cpy$v.s.coeff <- c( max(x[1],0), max( x[1], min( x[2], rf.debt ) ),
                               rf.debt, min(x[3],rf.level) , rf.level )
        # Lower parameters for surplus function
    params.cpy$s.shift <- c( -sum(x[-(1:3)] * s.idx.p[-1] ) / s.idx.p[1], x[-(1:3)] )
        # Shift parameters average to zero over the sample
    sim <- sim_core( c(1,interp$s.idx), sol.o$d.bar, sol.o$d.grid, sol.o$P, sol.o$Q,
                     params.cpy, d.init, TRUE, c(0,0,v.surp) )
    # Create the simulation (just need surplus function values but this is
    # fine too)
    if(s.dev) return(sd(sim[-1,'eps']))
    if(aad) return(mean(abs(sim[-1,'eps'])))
    return( sum( sim[-1,'eps'] ^ 2 ) )
  }

  ### 4. Minimize the error function ###
  control <- list(trace=print.level, maxit=100000)
  opt <- optim( init, err.fn, control=control )
      ## NEED TO FIT THE SHIFTS BETTER ##
  # control <- list(print_level=print.level, maxeval=1000, tol.abs=1e-04, algorithm="NLOPT_LN_COBYLA")
  # sol <- nloptr( init, err.fn, opts=control )
      # Minimize the error
  out <- list()
  out$v.s.coeff <- c( max(opt$par[1],0), max( opt$par[1], min( opt$par[2], rf.debt ) ),
                         rf.debt, min( opt$par[3], rf.level), rf.level )
  out$s.shift <- c( -sum(opt$par[-(1:3)] * s.idx.p[-1] ) / s.idx.p[1], opt$par[-(1:3)] )
  out$s.shift[s.idx.p==0] <- 0
  out$surp.sd <- err.fn(opt$par, s.dev=TRUE )
  out$aad <- err.fn(opt$par, aad=TRUE )
  out$p <- sol$p
  out$d <- sol$d
  out$v.surp <- v.surp
  out$err <- sol$err
  # out$v.s.coeff <- c( max(sol$solution[1],0), max( sol$solution[1], min( sol$solution[2], rf.debt ) ),
  #                     rf.debt, min( sol$solution[3], rf.level), rf.level )
  # out$s.shift <- c( -sum(sol$solution[-(1:3)] * s.idx.p[-1] ) / s.idx.p[1], sol$solution[-(1:3)] )
  # out$surp.sd <- sqrt( sol$objective )
      # Format output
  return(out)
}

price.diff.1.5 <- function( sol.o, params, sim, dta ){
# Computes the price difference between the 1 and 5 year bonds and compares it
# to the data analogue
  q.dta <- sapply( c(1,5), function(i) ( 1 + dta[,paste0('Int_',i,'y')] / 100 ) ^ -i )
      # The price of the assets in the data
  l.QQ <- lapply( c(1,5),
               function(i){
                 lambda <- 1 - 1 / (4*i)
                 q_hat_mat( sol.o$P, sol.o$d.bar, sol.o$Q, sol.o$Q, sol.o$d.grid, params$R,
                            params$G, params$s.shift, lambda, params$phi, sol.o$e.grid,
                            params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, 0 )
               } )
  q.sim <-sapply( l.QQ, function(QQ)
            apply( sim, 1, function(x)
              approx( sol.o$d.grid, QQ[x['idx'],], x['d.prime'] )$y ) )
      # Create the simulation of prices
  out <- data.frame( date=dta$date, q.dta, diff.dta=apply(q.dta,1,diff),
                     q.sim, diff.sim=apply(q.sim,1,diff) )
      # Create the output
  return(out)
}

price.diff.1.5.err <- function( sol.o, params, sim, dta ){
# Computes the error on the target term premium
  pd <- price.diff.1.5( sol.o, params, sim, dta )
  return( sum( ( pd$diff.dta - pd$diff.sim ) ^ 2, na.rm = TRUE ) )
}

price.diff.min <- function( params, interp, dta, sol, sol.o, h.0=c(25,-4), maxit=50 ){
# Minimizes the price difference error

  #### 1. Set up ####
  trial.diff <- rep( Inf, 3 ) ; trial.diff.tol <- 1e-04
  it <- 1 ; h <- h.0
  v.s.coeff <- params$v.s.coeff ; j <- 2 ;  # Start by reducing the level first
  i.j <- 1 ; i.j.max <- 5 ; err.old <- Inf ; fail.idx <- 0
  sol.out <- sol ; sol.o.out <- sol.o

  while( it < maxit & all(trial.diff > trial.diff.tol ) ){

    #### 2. Create the trial parameters                     ####
    params$v.s.coeff[c(3,5)][j] <- v.s.coeff[c(3,5)][j] + h[j]
    params$v.s.coeff[2] <- min( params$v.s.coeff[c(2,3)])

    #### 3. Solve the model & check solution                ####
    sol <- tryCatch( sol.wrapper( params, cbind(sol.out$p,sol.out$d) ),
                     error=function(e) list(p='fail', err=1) )
    if( any( sol$p == 'fail' ) | max(abs(sol$err)) > 1e-05 ){
        h[j] <- h[j] / 2
        message("\n********************************************")
        message("  Parameter iteration report:    ")
        message("    Iteration ", it )
        message("    j = ", j, ",  i.j = ", i.j )
        message("    Model solution failed" )
        message("    h                 = ", paste0( round(h,2), sep=', ' ) )
        message("    max(abs(sol$err)) = ", max(abs(sol$err)) )
        message("********************************************\n")
        if( h[j] < .1 * h.0[j] ){                           # Escape
          h[j] <- h.0[j]                                    # Reset step length
          j <- if( j==2 ) 1 else j + 1                      # Increment j
          fail.idx <- fail.idx + 1
          if( fail.idx == 3 ){
            message("\n********************************************")
            message("    Failing out")
            message("********************************************")
            return(list( sol=sol.out, sol.o=sol.o.out ))
          }
          i.j <- i.j + 1
        }
      }else{

        #### 4. Create outer solution and simulate            ####
        sol.o <- outer.wrapper( sol, params, Q.init = sol.o.out$Q,
                                d.grid.init=sol.o.out$d.grid )
        v.surp <- surp.fit(sol.o,interp,params,tgt$dta)
            # Create the fitted surpluses
        sim <- sim_core( c(1,interp$s.idx), sol.o$d.bar, sol.o$d.grid,
                         sol.o$P, sol.o$Q, params, tgt$dta$cnlb_gdp_lag[1],
                         TRUE, c(0,0,v.surp) )
            # The simulation

        #### 5. Create prices, measure error                  ####
        err <- price.diff.1.5.err( sol.o, params, sim, dta )

        #### 6. Accept or reject trial parameters             ####
        bo.improve <- err <= err.old + 1e-8
        # trial.diff[j] <-

        #### 7. Create new step or move on to next parameter  ####
        if( bo.improve  ){
          trial.diff[j] <- abs( v.s.coeff[c(3,5)][j] - params$v.s.coeff[c(3,5)][j] )
              # Measure the difference in the updated coefficients.
          v.s.coeff[c(3,5)][j] <- params$v.s.coeff[c(3,5)][j]
              # Store the improved coeffs
          message("\n********************************************")
          message("  Parameter iteration report:    ")
          message("    Iteration ", it )
          message("    j = ", j, ",  i.j = ", i.j )
          message("    bo.improve    = ", if(bo.improve) 'True' else 'False' )
          message("    v.s.coeff     = ", paste0( round(v.s.coeff,2), sep=', ' ) )
          message("    h             = ", paste0( round(h,2), sep=', ' ) )
          message("    err           = ", err )
          message("    err.old       = ", err.old )
          message("    trial.diff[j] = ", round(trial.diff[j]) )
          message("********************************************\n")

          fitted <- rf.fit( params, dta, v.s.coeff[3], v.s.coeff[5], v.surp,
                            sol, sol.o, init=c( params$v.s.coeff[c(1,2,4)], params$s.shift[-1] ) )
              # Refit the lower parameters
          params$v.s.coeff <- v.s.coeff <- fitted$v.s.coeff
          params$s.shift <- fitted$s.shift
          params$surp.sd <- fitted$surp.sd
              # Paste to parameters
          sol.out <- sol ; sol.o.out <- sol.o ; params.out <- params

          message("\n********************************************")
          message("  Refitting lower parameters:    ")
          message("    v.s.coeff     = ", paste0( round(v.s.coeff,2), sep=', ' ) )
          message("    s.shift       = ", paste0( round(params$s.shift,2), sep=', ' ) )
          message("    sd.surp       = ", round( params$surp.sd, 4 ) )
          message("********************************************\n")

          plot.surp( params, TRUE, c(0,1.2*params$v.s.coeff[3]), TRUE,
                     ylim=range( c( params$v.s.coeff[4:5], sim[,'s'] ) ) )
          points( dta$cnlb_gdp_lag, sim[,'s'], pch=16, cex=.75, col=sim[,'idx'] )

          i.j <- 1                                            # Reset j-trial counter
          h[j] <- h.0[j]                                      # Update h[j]
          j <- if( j==2 ) 1 else j + 1                        # Increment j
          err.old <- price.diff.1.5.err( sol.o, params, sim, dta )
              # Store error with new lower parameters

          fail.idx <- 0
        }else{
          message("\n********************************************")
          message("  Parameter iteration report:    ")
          message("    Iteration ", it )
          message("    j = ", j, ",  i.j = ", i.j )
          message("    bo.improve    = ", if(bo.improve) 'True' else 'False' )
          message("    h             = ", paste0( round(h,2), sep=', ' ) )
          message("    err           = ", err )
          message("    err.old       = ", err.old )
          message("********************************************\n")
          h[j] <- h[j] / 2                                    # Decrease the step
          i.j <- i.j + 1                                      # Increase j-trial counter
          if( h[j] < .1 * h.0[j] ){                           # Escape
            h[j] <- h.0[j]                                    # Reset step length
            j <- if( j==2 ) 1 else j + 1                      # Increment j
          }
        }
      }
    it <-  it + 1
  }
  return( list( sol=sol.out, sol.o=sol.o.out ) )
}

