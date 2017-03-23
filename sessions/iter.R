
e.grid <- e_grid_fn(params$surp.sd, 13)
d.grid <- d_grid_fn(sol.w$d, params$surp.sd )
Q <- Q_init( d.grid, sol.w$d, params$R )

d_prime( 0, sol.w$d[1] - 10, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$lambda,
         e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE )

microbenchmark( d.prime <- d_prime( 0, sol.w$d[1] - 10, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$lambda,
                                    e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, verbose = FALSE ) )

q.e <- q_e( sol.w$d[4] + 2, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$lambda,
            e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, verbose=TRUE )

microbenchmark(q_e( sol.w$d[4] + 2, sol.w$d, 1 / params$R[1], Q, d.grid, params$G, params$lambda,
                   e.grid, params$v.s.coeff, params$tri, matrix(0), FALSE, params$trans, verbose=FALSE ))
