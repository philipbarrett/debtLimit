library(microbenchmark)


## 1. Set up parameters ##
params <- list()
params$trans <- matrix( c(.6,.4,.4,.6), 2, 2 )
params$R <- c( 1.03, 1.02 )
params$G <- c( 1.04, 1.01 )
params$tri <- FALSE
# params$tri <- TRUE
params$v.s.coeff <- c( 10, -5, .75, .18, -.013, 2 )
# params$v.s.coeff <- c( 10, -5, .75, .18, -.014, .0001, 2 )
# params$v.s.coeff <- c( 50, 100, 150, -6, 1, 1 )
params$surp.sd <- 2
params$lambda <- .2 # 0 # .5
params$phi <- .6
params$cont.type <- 'avg'
params$q.e <- c(0)
params$def <- matrix(0,1,1)
params$diff.method <- "numerical"
params$qd.internal <- TRUE

## 2. Plot surplus function(s) ##
plot.surp(params, ylim=c(-5, 15))

## 3. Solve the model ##
d.init <- sol.nonstoch(params)
    # Use the nonstohastic solution as an initial guess for d

## 4. Plot the price functions ##
params$cont.type <- 'avg'
plot.q( c(.2,.2), params )
params$cont.type <- 'low'
plot.q( c(.2,.2), params )

## 5. Create initial guesses ##
params$cont.type <- 'avg'
init.guess.avg <- d.init.p(params)
params$cont.type <- 'low'
init.guess.low <- d.init.p(params)


## 6. Plot the z functions ##
# init.guess <- cbind( rep(.5,length(params$R)), mean( d.init ) )
params$cont.type <- 'avg'
plot.z( c(0,0), init.guess.avg, params )
microbenchmark(zed( c(0,.3), d.init, params, c(0) ))
    # Eval: ~10ms
microbenchmark(zed_2( c(0,.3), d.init, params ))
    # Eval: ~30ms (using numerical differentiation)



## 7. Try solving the model ##
params$cont.type <- 'avg'
sol <- sol.wrapper( params, init.guess, 'core.nl' )
plot.z( sol$p, sol$d, params, sol$qd, sol$qe, sol$def )
plot.q( sol$p, params )
params$cont.type <- 'low'
sol <- sol.wrapper( params, init.guess, 'core.nl' )
plot.z( sol$p, sol$d, params, sol$qd, sol$qe, sol$def )
plot.q( sol$p, params )
