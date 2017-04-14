library(microbenchmark)
library(debtLimits)

## 1. Set up parameters ##
params <- list()
# params$trans <- matrix( c(.8,.2,.2,.8), 2, 2 ) # matrix( c(.6,.4,.4,.6), 2, 2 )
# params$R <- c( 1.03, 1.02 )
# params$G <- c( 1.04, 1.01 )

# params$trans <- matrix( c(.7,.3,
#                           .3, .7), 2, 2, byrow = TRUE ) # matrix( c(.6,.4,.4,.6), 2, 2 )
# params$R <- 1 + c( .03, .03) / 4 #, 1.06, 1.02 )
# params$G <- 1 + c( .02, .02) / 4 #, 1.05, .99 )

load('data/surpEst.rdta')
params <- params.alt
# params$G[2] <- params$G[2] + .005
    # Increase the low-G state a little

params$tri <- FALSE
# params$tri <- TRUE
# params$v.s.coeff <- c( 40, .3, -.15, .2, -.035, 0 )
par(mfrow=c(1,1), mar=rep(4,4))
plot.surp(params)
# params$v.s.coeff <- c( 10, -5, .75, .18, -.014, .0001, 2 )
# params$v.s.coeff <- c( 50, 100, 150, -6, 1, 1 )
# params$surp.sd <- .01
params$lambda <- 0 # 0 # .5
params$phi <- .6
params$cont.type <- 'avg'
params$q.e <- c(0)
params$def <- matrix(0,1,1)
params$diff.method <- "ana"
params$d.tri <- FALSE # TRUE # FALSE      # Triangular distribution for surplus shocks
params$inner.method <- 'all'
An <- 1 / params$R
Bn <- rep( -1, length(params$R) )
Cn <- An
def <- matrix(0,1,1)

## 2. Plot surplus function(s) ##
plot.surp(params)

## 3. Solve the nonstochastic model ##
d.init <- sol.nonstoch(params)
    # The nonstohastic solution

## 4. Plot the price functions ##
params$cont.type <- 'avg'
plot.q( rep(.05,length(params$R)), params )
params$cont.type <- 'low'
plot.q( rep(.05,length(params$R)), params )
params$cont.type <- 'fix'
plot.q( rep(.05,length(params$R)), params, An )

## 5. Check the price derivatives ##
nn <- length(params$R)
q.p.ana.avg <- q_d_p(params$R, rep(.2,nn), params$trans, params$lambda, params$phi, nn, "avg",
                    "ana", params$G, An, Bn, def )
q.p.num.avg <- q_d_p(params$R, rep(.2,nn), params$trans, params$lambda, params$phi, nn, "avg",
                    "num", params$G, An, Bn, def )
print( q.p.ana.avg / q.p.num.avg )
    # The analytical an numerical derivates for the avg case
q.p.ana.fix <- q_d_p(params$R, rep(.2,nn), params$trans, params$lambda, params$phi, nn, "fix",
                     "ana", params$G, An, 0*Bn, def )
q.p.num.fix <- q_d_p(params$R, rep(.2,nn), params$trans, params$lambda, params$phi, nn, "fix",
                     "num", params$G, An, Bn, def )
print( q.p.ana.fix / q.p.num.fix )
    # The analytical an numerical derivates for the fix case.  By construction,
    # these are only equal when Bn=0.


## 6. Create initial guesses ##
params$cont.type <- 'avg'
zed( rep(0,nn), rep(min(d.init)-20,nn), params, An, Cn, def, 1 )
    # Use case with print_level > 0
surp.sd <- params$surp.sd
params$surp.sd <- .1 * surp.sd
zed( rep(0,nn), rep(min(d.init)-20,nn), params, An, Cn, def, 1 )
plot.z( rep(1e-02,nn), rep(min(d.init)-20,nn), params, An, Bn, Cn, def )
    # Just check to see that default probs are lower when surp SD is lower
params$surp.sd <- surp.sd
plot.z( rep(1e-02,nn), rep(min(d.init),nn), params, An, Bn, Cn, def )
plot.z.d( rep(1e-02,nn), rep(min(d.init),nn), params, An, Bn, Cn, def )
d.init.2 <- d_init_p( params, rep(1e-04,nn), rep(min(d.init),nn), An, Bn, Cn, def, max(d.init), 20, 1,max_outer = 4 )
plot.z( rep(.0,nn), d.init.2, params, An, Bn, Cn, def, xlim=c(0,1e-02) )
plot.z.d( rep(1e-02,nn), d.init.2, params, An, Bn, Cn, def )
p.init <- p_init_d( params, rep(0,nn), d.init.2, An, Bn, Cn, def )
init.guess.avg <- list( p=p.init, d=d.init.2 )

# init.guess.avg <- p.d.init(params, x.sd = -2)
# init.guess.avg <- list( p=rep(0,length(params$R)), d=d.init.p(params, p=rep(0,length(params$R)), x.sd = -1) )
# init.guess.avg <- list( p=p.init.d(params, rep(0,length(params$R)),
#                                    rep(min(d.init),length(params$R)), An, Bn, Cn, def ),
#                         d=rep(min(d.init),length(params$R)) )
params$cont.type <- 'low'
init.guess.low <- d.init.p(params, x.sd = -1)


## 7. Plot the z functions ##
# init.guess <- cbind( rep(.5,length(params$R)), mean( d.init ) )
params$diff.method <- 'num'
params$cont.type <- 'avg'
plot.z( init.guess.avg$p, init.guess.avg$d, params )
params$diff.method <- 'ana'
plot.z( init.guess.avg$p, init.guess.avg$d, params )

params$d.tri <- FALSE
params$diff.method <- 'num'
plot.z( init.guess.avg$p, init.guess.avg$d, params )
params$diff.method <- 'ana'
plot.z( init.guess.avg$p, init.guess.avg$d, params )

microbenchmark(zed( 0*params$R, d.init, params, An, Cn, def,0 ))
    # Eval: ~12ms

params$diff.method <- 'num'
a <- zed_2( rep(0,nn), d.init, params, An, Bn, Cn, def )
print(microbenchmark(zed_2( 0 * params$R, d.init, params, An, Bn, Cn, def ), times=10000))
    # Eval: ~48ms (using numerical differentiation)
params$diff.method <- 'ana'
b <- zed_2( rep(0,nn), d.init, params, An, Bn, Cn, def )
print(microbenchmark(zed_2( 0 * params$R, d.init, params, An, Bn, Cn, def ), times=10000))
    # Eval: ~33ms (using numerical differentiation)
print(a/b)


## 7. Try solving the model ##
trans <- params$trans
R <- params$R
G <- params$G

params$cont.type <- 'avg'
params$lambda <- .0
params$surp.sd <- surp.sd * .3
# params$trans <- trans
params$trans <- matrix( c(.98,.02,.02,.98), 2, 2, byrow=TRUE)
# params$trans <- matrix( c(1,0,0,1), 2, 2)

sol.s <- sol.wrapper( params, plot.on=T ) #, do.call( 'cbind', init.guess.avg )
plot.z( sol.s$p, sol.s$d, sol.s$params, xlim=c(0,3e-03), ylim=c(0,3e-03) )

params$lambda <- .0
params$surp.sd <- surp.sd * .45
params$v.s.coeff[2] <- .3
plot.surp(params)
sol.b <- sol.wrapper( params, plot.on=T ) #, do.call( 'cbind', init.guess.avg )
plot.z( sol.b$p, sol.b$d, sol.b$params, xlim=c(0,5e-03), ylim=c(0,5e-03) )
params$lambda <- .95
sol.t <- sol.wrapper( params ) #, do.call( 'cbind', init.guess.avg )
plot.z( sol.t$p, sol.t$d, sol.t$params, xlim=c(0,5e-02), ylim=c(0,5e-02) )

# params$cont.type <- 'fix'
# params$lambda <- .5
# params$it <- 20
# sol.f <- sol.wrapper( params, init.guess = init.guess, An=An, Bn=Bn, Cn=Cn )
#     # This is solving but reporting failure!!!
# plot.z( sol.f$p, sol.f$d, sol.f$params, An, Bn, Cn )
# plot.z( sol.w$p, sol.w$d, sol.w$params )
#
# plot.z.d( sol.f$p, sol.f$d, sol.f$params, An, Bn, Cn )
# plot.z.d( sol.w$p, sol.w$d, sol.w$params )
# plot.z.d( sol.s$p, sol.s$d, sol.w$params )

save( params, file='data/sample.dta')
