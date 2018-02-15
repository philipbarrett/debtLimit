library(microbenchmark)
library(debtLimits)

## 1. Set up parameters ##
# params <- list()
# params$trans <- matrix( c(.8,.2,.2,.8), 2, 2 ) # matrix( c(.6,.4,.4,.6), 2, 2 )
# params$R <- c( 1.03, 1.02 )
# params$G <- c( 1.04, 1.01 )

# params$trans <- matrix( c(.9,.1,
#                           .1, .9), 2, 2, byrow = TRUE ) # matrix( c(.6,.4,.4,.6), 2, 2 )
# params$R <- c( 1.03, 1.02 ) #, 1.06, 1.02 )
# params$G <- c( 1.04, 1.01 )#, 1.05, .99 )
#
# params$tri <- FALSE
# # params$tri <- TRUE
# params$v.s.coeff <- c( 10, -5, .75, .18, -.013, 2 )
# # params$v.s.coeff <- c( 10, -5, .75, .18, -.014, .0001, 2 )
# # params$v.s.coeff <- c( 50, 100, 150, -6, 1, 1 )
# params$surp.sd <- 2
# params$lambda <- 0 # .2 # 0 # .5
# params$phi <- .6
# params$cont.type <- 'avg'
# params$q.e <- c(0)
# params$def <- matrix(0,1,1)
# params$diff.method <- "ana"
# params$d.tri <- FALSE # TRUE #       # Triangular distribution for surplus shocks
# params$inner.method <- 'all'
load( 'data/sample.dta')
lambda <- params$lambda
params$lambda <- 0

An <- 1 / params$R
Bn <- rep( -1, length(params$R) )
Cn <- An
def <- matrix(0,1,1)
nn <- length(params$R)

## 2. Create an approximate solution
plot.surp(params)
sol.w <- sol.wrapper( params )
plot.z( sol.w$p, sol.w$d, sol.w$params, xlim=c(0,2*max(sol.w$p)), ylim=c(0,2*max(sol.w$p)) )

## 3. Create solution grids ##
e.grid <- e_grid_fn(params$surp.sd, 21, params$d.tri )
d.grid <- d_grid_fn(sol.w$d, params$surp.sd, x_sd_mult = 2, n_pts=15 )
Q <- Q_init( d.grid, sol.w$d, params$R )

## 4. Solve for next-period debt ##
d_prime( 0, sol.w$d[1] - 10, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$s.shift, params$lambda,
         e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, print_level=2 )
microbenchmark( d.prime <- d_prime( 0, sol.w$d[1] - 10, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$s.shift,params$lambda,
                                    e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, print_level=0 ) )
    # ~ 25ms

## 5. Solve for expected continuation debt prices ##
q.e <- q_e( d.grid[13], sol.w$d, 1 / params$R, Q, d.grid, params$G, params$s.shift, params$lambda,
            e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, print_level=1 )
microbenchmark(q_e( sol.w$d[4] + 2, sol.w$d, 1 / params$R, Q, d.grid, params$G, params$s.shift, params$lambda,
                   e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, print_level=0 ))
    # ~ 60ms

## 6. Solve for actual debt price ##

qq <- q_hat_fn( d.grid[13], rep(.01,nn), sol.w$d, 1 / params$R, Q, d.grid, params$R,
                 params$G, params$s.shift, params$lambda, params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, print_level=1 )
microbenchmark( qq <- q_hat_fn( d.grid[2], rep(.02,nn), sol.w$d, 1 / params$R, Q, d.grid, params$R,
                 params$G, params$s.shift, params$lambda, params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, print_level=0 ) )
    # ~94ms
QQ <- q_hat_mat( matrix(.01,dim(Q)[1],dim(Q)[2]), sol.w$d, Q, Q, d.grid, params$R,params$G,
                 params$s.shift, params$lambda, params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, 2 )
microbenchmark(QQ <-
         q_hat_mat( matrix(.01,dim(Q)[1],dim(Q)[2]), sol.w$d, Q, Q, d.grid, params$R,params$G,
                    params$s.shift, params$lambda, params$phi, e.grid, params$v.s.coeff,
                    params$tri, matrix(0), FALSE, params$trans, FALSE ))
    # .002s

## 7. Now just solve and plot ##
params$x.mult <- 20
sol.w <- sol.wrapper(params)
sol.o <- outer.wrapper( sol.w, params )
plot.sol(sol.o)
err.o <- outer.err(sol.o, params )
plot.err(err.o)
    # Great!  The solutions line up for lambda = 0 :) :) :)


## 8. Try with lambda != 0 ##
params$cont.type <- 'avg'
# params$lambda <- lambda
params$lambda <- .7
    # For .7 and below we get exact solution, for .8 and above we fail to get convergence.
sol.in <- sol.wrapper(params)
plot.z( sol.w$p, sol.w$d, sol.w$params, xlim=c(0,2*max(sol.w$p)), ylim=c(0,2*max(sol.w$p)) )
plot.z( sol.in$p, sol.in$d, sol.in$params, xlim=c(0,2*max(sol.in$p)), ylim=c(0,2*max(sol.in$p)) )
    # Compare the two solutions
l.p <- lapply( list( sol.w, sol.in), function(x) x$p )
l.d <- lapply( list( sol.w, sol.in), function(x) x$d )
    # Store the default probabilities and debt levels

params$cont.type <- 'fix'
params$d.pts <- 51
params$x.mult <- 40
params$maxit <- 50
sol.out <- outer.wrapper( sol.in, params )
plot.sol(sol.out)
err.out <- outer.err(sol.out, params )
plot.err(err.out)
    # This looks pretty damn close! :)
sol.out.b <- outer.wrapper( sol.w, params )
plot.sol(sol.out.b)
err.out.b <- outer.err(sol.out.b, params )
plot.err(err.out.b)
    # And this doesn't quite work! :) Which is good!

## Try to improve it next:
# l.ABC <- ABC( sol.out, params )
# sol.in.2 <- sol.wrapper(params, An=l.ABC$An, Bn=l.ABC$Bn, Cn=l.ABC$Cn  )
#    # No bueno



## TODO:
## FIGURE OUT HOW TO MAKE THE INNER LOOP CONVERGE FOR SURE
## LOOP
