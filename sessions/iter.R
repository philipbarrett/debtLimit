params$d.tri <- TRUE

e.grid <- e_grid_fn(params$surp.sd, 21, params$d.tri )
d.grid <- d_grid_fn(sol.w$d, params$surp.sd, x_sd_mult = 2, n_pts=15 )
Q <- Q_init( d.grid, sol.w$d, params$R )

d_prime( 0, sol.w$d[1] - 10, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$lambda,
         e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, print_level=2 )

microbenchmark( d.prime <- d_prime( 0, sol.w$d[1] - 10, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$lambda,
                                    e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, print_level=0 ) )

q.e <- q_e( d.grid[13], sol.w$d, 1 / params$R, Q, d.grid, params$G, params$lambda,
            e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, print_level=1 )

microbenchmark(q_e( sol.w$d[4] + 2, sol.w$d, 1 / params$R, Q, d.grid, params$G, params$lambda,
                   e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, print_level=0 ))

qq <- q_hat_fn( d.grid[13], rep(.01,4), sol.w$d, 1 / params$R, Q, d.grid, params$R,
                 params$G, params$lambda, params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, print_level=1 )

microbenchmark( qq <- q_hat_fn( d.grid[2], rep(.02,4), sol.w$d, 1 / params$R, Q, d.grid, params$R,
                 params$G, params$lambda, params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, print_level=0 ) )

QQ <- q_hat_mat( matrix(.01,dim(Q)[1],dim(Q)[2]), sol.w$d, Q, Q, d.grid, params$R,params$G, params$lambda,
                 params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, 2 )

microbenchmark(QQ <-
         q_hat_mat( matrix(.01,dim(Q)[1],dim(Q)[2]), sol.w$d, Q, Q, d.grid, params$R,params$G, params$lambda,
                                params$phi, e.grid, params$v.s.coeff, params$tri,
                                matrix(0), FALSE, params$trans, FALSE ))

plot( range(d.grid), range(QQ), type='n' )
abline(v=sol.w$d, lty=2, col=1:4 )
abline(v=d.grid, lty=1, lwd=.25 )
for( i in 1:4 ) lines( d.grid, QQ[i,], col=i, lwd=2 )

ZZ <- ziter( 0 * Q, sol.w$d, Q, Q, d.grid, e.grid, params$tri, matrix(0,1,1), FALSE,
                params, An, Cn, 1 )

system.time(ZZ <- ziter( 0 * Q, sol.w$d, Q, Q, d.grid, e.grid, params$tri, matrix(0,1,1), FALSE,
                         params, An, Cn, 0 ) )

plot( range(d.grid), range(ZZ), type='n' )
abline(v=sol.w$d, lty=2, col=1:4, lwd=2 )
abline(v=d.grid, lty=1, lwd=.25 )
for( i in 1:4 ) lines( d.grid, ZZ[i,], col=i, lwd=2 )
abline(h=sol.w$p, lty=2 )

QQ <- q_hat_mat( ZZ, sol.w$d, Q, Q, d.grid, params$R,params$G, params$lambda,
                 params$phi, e.grid, params$v.s.coeff, params$tri,
                 matrix(0), FALSE, params$trans, 2 )

plot( range(d.grid), range(QQ), type='n' )
abline(v=sol.w$d, lty=2, col=1:4, lwd=2 )
abline(v=d.grid, lty=1, lwd=.25 )
abline(h=0)
for( i in 1:4 ) lines( d.grid, QQ[i,], col=i, lwd=2 )


DD <- d_prime_mat( sol.w$d, QQ, QQ, d.grid, params$G, params$lambda, e.grid, params$trans,
                   params$v.s.coeff, params$tri, matrix(0,1,1), FALSE )

plot( range(d.grid), range(DD), type='n' )
abline(v=sol.w$d, lty=2, col=1:4, lwd=2 )
abline(h=sol.w$d, lty=2, col=1:4, lwd=2 )
abline(v=d.grid, lty=1, lwd=.25 )
for( i in 1:4 ) lines( d.grid, DD[i,], col=i, lwd=2 )

QE <- qe_mat( sol.w$d, QQ, QQ, d.grid, params$G, params$lambda, e.grid, params$trans,
                   params$v.s.coeff, params$tri, matrix(0,1,1), FALSE )

plot( range(d.grid), range(QE), type='n' )
abline(v=sol.w$d, lty=2, col=1:4, lwd=2 )
abline(v=d.grid, lty=1, lwd=.25 )
for( i in 1:4 ) lines( d.grid, QE[i,], col=i, lwd=2 )


## TODO:
## BUILD A SUMMARY CHART OF PRICES AND PROBABILITIES AND DEFAULT THRESHOLDS
## ALSO ERROR CHECKS
## EXTEND GRID ABOVE TOP DEFAULT RATE
