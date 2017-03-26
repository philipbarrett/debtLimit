## 1. Set up parameters ##
library(microbenchmark)
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
params$d.tri <- FALSE      # Triangular distribution for surplus shocks
An <- 1 / params$R
Bn <- rep( -1, length(params$R) )
Cn <- An
def <- matrix(0,1,1)

## 2. Create an approximate solution
sol.w <- sol.wrapper( params )
plot.z( sol.w$p, sol.w$d, params )

## 3. Create solution grids ##
e.grid <- e_grid_fn(params$surp.sd, 21, params$d.tri )
d.grid <- d_grid_fn(sol.w$d, params$surp.sd, x_sd_mult = 2, n_pts=15 )
Q <- Q_init( d.grid, sol.w$d, params$R )

## 4. Solve for next-period debt ##
d_prime( 0, sol.w$d[1] - 10, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$lambda,
         e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, print_level=2 )
microbenchmark( d.prime <- d_prime( 0, sol.w$d[1] - 10, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$lambda,
                                    e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, print_level=0 ) )
    # ~ 35ms

## 5. Solve for expected continuation debt prices ##
q.e <- q_e( d.grid[13], sol.w$d, 1 / params$R, Q, d.grid, params$G, params$lambda,
            e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, print_level=1 )
microbenchmark(q_e( sol.w$d[4] + 2, sol.w$d, 1 / params$R, Q, d.grid, params$G, params$lambda,
                   e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, print_level=0 ))

## 6. Solve for actual debt price ##
qq <- q_hat_fn( d.grid[13], rep(.01,4), sol.w$d, 1 / params$R, Q, d.grid, params$R,
                 params$G, params$lambda, params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, print_level=1 )
microbenchmark( qq <- q_hat_fn( d.grid[2], rep(.02,4), sol.w$d, 1 / params$R, Q, d.grid, params$R,
                 params$G, params$lambda, params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, print_level=0 ) )
    # ~290ms
QQ <- q_hat_mat( matrix(.01,dim(Q)[1],dim(Q)[2]), sol.w$d, Q, Q, d.grid, params$R,params$G, params$lambda,
                 params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, 2 )
microbenchmark(QQ <-
         q_hat_mat( matrix(.01,dim(Q)[1],dim(Q)[2]), sol.w$d, Q, Q, d.grid, params$R,params$G, params$lambda,
                                params$phi, e.grid, params$v.s.coeff, params$tri,
                                matrix(0), FALSE, params$trans, FALSE ))
    # .02s

## 7. Now just solve and plot ##
sol.o <- outer.wrapper( d.grid, e.grid, sol.w$d, params, An, Bn, Cn )
plot.sol((sol.o))


## TODO:
## ALSO ERROR CHECKS
## THEN EXTRACT AN, BN, CN
## LOOP
