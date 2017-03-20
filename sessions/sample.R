library(microbenchmark)


## 1. Set up parameters ##
params <- list()
# params$trans <- matrix( c(.8,.2,.2,.8), 2, 2 ) # matrix( c(.6,.4,.4,.6), 2, 2 )
# params$R <- c( 1.03, 1.02 )
# params$G <- c( 1.04, 1.01 )

params$trans <- matrix( c(.6,.2,.1,.1,
                          .2,.6,.1,.1,
                          .4,.1,.4,.1,
                          .1,.4,.1,.4), 4, 4, byrow = TRUE ) # matrix( c(.6,.4,.4,.6), 2, 2 )
params$R <- c( 1.03, 1.02, 1.06, 1.02 )
params$G <- c( 1.04, 1.01, 1.05, .99 )

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
params$diff.method <- "num"
An <- 1 / params$R
Bn <- rep( -1, length(params$R) )
Cn <- An
def <- matrix(0,1,1)

## 2. Plot surplus function(s) ##
plot.surp(params, ylim=c(-5, 15))

## 3. Solve the model ##
d.init <- sol.nonstoch(params)
    # Use the nonstohastic solution as an initial guess for d

## 4. Plot the price functions ##
params$cont.type <- 'avg'
plot.q( rep(.05,length(params$R)), params )
params$cont.type <- 'low'
plot.q( rep(.05,length(params$R)), params )
params$cont.type <- 'fix'
plot.q( rep(.05,length(params$R)), params, An )

## 5. Check the price derivatives ##
q.p.ana.avg <- q_d_p(params$R, rep(.2,4), params$trans, params$lambda, params$phi, 4, "avg",
                    "ana", params$G, An, Bn, def )
q.p.num.avg <- q_d_p(params$R, rep(.2,4), params$trans, params$lambda, params$phi, 4, "avg",
                    "num", params$G, An, Bn, def )
print( q.p.ana.avg / q.p.num.avg )
    # The analytical an numerical derivates for the avg case
q.p.ana.fix <- q_d_p(params$R, rep(.2,4), params$trans, params$lambda, params$phi, 4, "fix",
                     "ana", params$G, An, 0*Bn, def )
q.p.num.fix <- q_d_p(params$R, rep(.2,4), params$trans, params$lambda, params$phi, 4, "fix",
                     "num", params$G, An, Bn, def )
print( q.p.ana.fix / q.p.num.fix )
    # The analytical an numerical derivates for the fix case.  By construction,
    # these are only equal when Bn=0.


## 6. Create initial guesses ##
params$cont.type <- 'avg'
init.guess.avg <- p.d.init(params, x.sd = -2)
init.guess.avg <- list( p=rep(0,length(params$R)), d=d.init.p(params, p=rep(0,length(params$R)), x.sd = -1) )
init.guess.avg <- list( p=p.init.d(params, rep(0,length(params$R)),
                                   rep(min(d.init),length(params$R)), An, Bn, Cn, def ),
                        d=rep(min(d.init),length(params$R)) )
params$cont.type <- 'low'
init.guess.low <- d.init.p(params, x.sd = -1)


## 7. Plot the z functions ##
# init.guess <- cbind( rep(.5,length(params$R)), mean( d.init ) )
params$diff.method <- 'num'
params$cont.type <- 'avg'
plot.z( init.guess.avg$p, init.guess.avg$d, params )
params$diff.method <- 'ana'
plot.z( init.guess.avg$p, init.guess.avg$d, params )

microbenchmark(zed( 0*params$R, d.init, params, An, Cn, def ))
    # Eval: ~12ms

params$diff.method <- 'num'
a <- zed_2( rep(0,4), d.init, params, An, Bn, Cn, def )
print(microbenchmark(zed_2( 0 * params$R, d.init, params, An, Bn, Cn, def ), times=10000))
    # Eval: ~48ms (using numerical differentiation)
params$diff.method <- 'ana'
b <- zed_2( rep(0,4), d.init, params, An, Bn, Cn, def )
print(microbenchmark(zed_2( 0 * params$R, d.init, params, An, Bn, Cn, def ), times=10000))
    # Eval: ~33ms (using numerical differentiation)
print(a/b)



## 7. Try solving the model ##
params$cont.type <- 'avg'

sol.i <- init.guess.avg
plot.z( init.guess.avg$p, init.guess.avg$d, params )
for( j in 1:8 ){
  for( i in 1:length(params$R) ){
    p.guess <- sol.i$p
    d.guess <- sol.i$d
    p.guess[i] <- p.init.d( params, sol.i$p, sol.i$d, An, Bn, Cn, def )[i]
    sol.i <- sol.wrapper( params, cbind( p.guess, d.guess ), 'core.nl.i', i=i )
    plot.z( sol.i$p, sol.i$d, params )
  }
}
system.time( sol <- sol.wrapper( params, cbind( sol.i$p, sol.i$d ),
                                 'core.nl' ) )
plot.z( sol$p, sol$d, params )
# plot.q( sol$p, params )
# params$cont.type <- 'low'
# sol <- sol.wrapper( params, cbind( init.guess.avg$p, init.guess.avg$d ), 'core.nl' )
# plot.z( sol$p, sol$d, params, sol$qd, sol$qe, sol$def )
# plot.q( sol$p, params )
