rm(list=ls())

## Code to test whether the markov switching process is recovering the DGP
A <- matrix( c( .5, .2, .6, -.2 ), 2, 2, byrow=TRUE )
Sigma <- matrix( c( 1, .1, .1, .5 ), 2, 2, byrow=TRUE )
a <- matrix( c(1,0, 0, 0, 0, 1), 2, 3 )
mu <- solve( A, a )
trans <- matrix( c( .9, .1, 0,
                    .25, .5, .25,
                    0, .3, .7), 3, 3, byrow=TRUE )
mu.ave <- mu %*% ( trans %^% 100 )[1,]
    # The average of averages
m <- 20000
nn <- 2
n.state <- 3

## Simulate the model and check the results
ms <- markov_sim( m, trans, 0, m )
vs <- ms_var(a, A, Sigma, ms, m)
    # The simulated VAR with a Markov-switching mean
vs.mu <- by( t(vs), ms, function(d) apply( d, 2, mean ) )
    # The conditional averages (will not converge)
vs.lead <- vs[,-1]
vs.lag <- vs[,-m]
l.reg <- lapply( 1:n.state, function(i)
    lapply( 1:nn, function(j) lm( vs.lead[j,ms[-1]==i-1] ~ t( vs.lag[,ms[-1]==i-1] ) ) ) )
a.reg <- sapply( l.reg, function(x) sapply( x, function(y) y$coefficients['(Intercept)'] ) )
    # Looks pretty good
A.reg <- lapply( l.reg, function(x) t( sapply( x, function(y) y$coefficients[-1] ) ) )
    # Looks good too
Sigma.reg <- lapply( l.reg, function(x) var( sapply( x, function(y) y$residuals ) ) )
    # Also good

## Check the filtering probabilities ##
p.0 <- 1 * ( 0:(n.state-1) == ms[1] )
m.par <- rbind( a, matrix( c( A, Sigma[lower.tri(Sigma,TRUE)] ), 2*nn^2 - .5*nn*(nn-1), n.state ) )
g.par <- c( p.0[-n.state], trans[,-n.state],  a, A, Sigma[lower.tri(Sigma,TRUE)] )
    # Initialize grand parameter vector
probs.true <- msw_var_lhood_p( vs, p.0, trans, m.par )
    # The filtering and prediction probailities
n.pds <- 150
plot( c(1,n.pds), c(0,1), type='n' )
for( i in 1:n.state ) lines( 1:n.pds, probs.true[i,1:n.pds], col=i, lwd=2 )
for( i in 1:n.state ) lines( 1:n.pds, ms[1:n.pds] == i-1, col=i, lwd=1 )
    # Looks nice

## Now do MLE for the VAR parameters holding fixed the Markov parameters at
## their correct values ##
par.var <- c( a, A, Sigma[lower.tri(Sigma,TRUE)] )
    # The VAR parameters
var.lhood <- function(v.par){
  m.a <- matrix( v.par[1:(nn*n.state)], nn, n.state )
  m.par <- rbind( m.a, matrix( v.par[-(1:(nn*n.state))], length(v.par)- nn * n.state, n.state ) )
  msw_var_lhood( vs, p.0, trans, m.par )
}
# The likelihood for the VAR part
control <- list( abstol = 1e-12, reltol = 1e-12, trace = 1 )
    # Optimization controls
opt.var <- optim( par.var, var.lhood, method='BFGS', control = control )
set.seed(321)
opt.var.2 <- optim( par.var + c( rnorm( nn*n.state + nn^2, 0, 1),
                                 rnorm( .5*nn*(nn+1), 0, .1) ),
                    var.lhood, method='BFGS', control = control )
print( cbind( par.var, c(a.reg, A.reg[[1]], Sigma.reg[[1]][lower.tri(Sigma.reg[[1]],TRUE)]),
              opt.var$par, opt.var.2$par) )
    # Looks good

## Now do MLE for the Markov parameters holding fixed the VAR parameters at
## their correct values ##
v.par <- rbind( a, matrix( c( A, Sigma[lower.tri(Sigma,TRUE)]),
                           nn ^ 2 + .5 * nn * ( nn + 1 ), n.state ) )
    # Matrix of VAR parameters
ms.par.init <- c( p.0[-n.state], trans[,-n.state] )
    # Initialize the Markov parameters
ms.lhood <- function(ms.par){
# Likelihood for Markov parameters
  p.0.x <- ms.par[1:(n.state-1)]
  p.0.x <- c( p.0.x, 1 - sum(p.0.x) )
  p.0.ms <- pmax( pmin( p.0.x, 1 ), 0 )
      # p.0
  P.x <- matrix( ms.par[-(1:(n.state-1))], n.state, n.state-1 )
  P.x <- cbind( P.x, 1 - apply( P.x, 1, sum ) )
  P.ms <- pmax( pmin( P.x, 1 ), 0 )
      # P
  penalty <- 10000 * sum( p.0.x * ( p.0.x - 1 ) * ( p.0.x != p.0.ms ) ) +
                10000 * sum( P.x * ( P.x - 1 ) * ( P.x != P.ms ) )
      # Penalty for exceeding bounds
  msw_var_lhood( vs, p.0.ms, P.ms, v.par ) + penalty
}
opt.ms <- optim( ms.par.init, ms.lhood, method='BFGS', control = control )
set.seed(123)
opt.ms.2 <- optim( ms.par.init + c( rep(0,n.state-1), rnorm( length(ms.par.init) - n.state + 1, 0, 3e-01) ),
                   ms.lhood, method='BFGS', control = control )
print( cbind( ms.par.init, opt.ms$par, opt.ms.2$par) )
    # Not so sure about this.  Not stable or accurate.

## The same, but excluding p.0 from the estimations ##
ms.lhood.P <- function(P.par){
# Likelihood for Markov transition matrix
  P.x <- matrix( P.par, n.state, n.state-1 )
  P.x <- cbind( P.x, 1 - apply( P.x, 1, sum ) )
  P.ms <- pmax( pmin( P.x, 1 ), 0 )
      # P
  penalty <- 10000 * sum( P.x * ( P.x - 1 ) * ( P.x != P.ms ) )
      # Penalty for exceeding bounds
  msw_var_lhood( vs, p.0, P.ms, v.par ) + penalty
}
opt.ms.P <- optim( c(trans[,-n.state]), ms.lhood.P, method='BFGS', control = control )
set.seed(42)
q <- Inf
P.2 <- NULL
for( i in 1:10 ){
  opt.ms.P.2.cand <- optim( runif( n.state*(n.state-1), 0, 1 ), ms.lhood.P,
                       method='BFGS', control = control )
  if( opt.ms.P.2.cand$value < q ){
    opt.ms.P.2 <- opt.ms.P.2.cand
    q <- opt.ms.P.2$value
  }
  P.2 <- cbind( P.2, opt.ms.P.2.cand$par )
}
print( cbind( c( trans[,-n.state] ) , c(opt.ms.P$par), c( opt.ms.P.2$par ) ) )
    # Perhaps multiple random starts are the solution

## Now try p.0 alone ##
ms.lhood.0 <- function(p.0.x){
# Likelihood for Markov parameters
  p.0.x <- c( p.0.x, 1 - sum(p.0.x) )
  p.0.ms <- pmax( pmin( p.0.x, 1 ), 0 )
      # p.0
  penalty <- 10000 * sum( p.0.x * ( p.0.x - 1 ) * ( p.0.x != p.0.ms ) )
      # Penalty for exceeding bounds
  msw_var_lhood( vs, p.0.ms, trans, v.par ) + penalty
}
opt.ms.0 <- optim( p.0[-n.state], ms.lhood.0, method='BFGS', control = control )
set.seed(321)
opt.ms.0.2 <- optim( runif( n.state-1, 0, 1 ), ms.lhood.0, method='BFGS', control = control )
print( cbind( p.0[-n.state], opt.ms.0$par, opt.ms.0.2$par ) )
    # This is ok

## Last, try P with p.0 given by the L/R mean ##




























