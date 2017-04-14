library(vars)
library(MASS)

A <- matrix( c( .5, .2, .6, -.2 ), 2, 2, byrow=TRUE )
Sigma <- matrix( c( 1, .1, .1, .5 ), 2, 2, byrow=TRUE )
m <- 100000
nn <- 2
# Set up

e <- t( mvrnorm( m, rep(0,2), Sigma ) )
# Errors
Y <- matrix( 0, nn, m )
# Initialize the VAR data
for( i in 2:m ){
  Y[,i] <- A %*% Y[,i-1] + e[,i]
}
a.switch <- rep(1,nn)
A.switch <- Sigma.switch <- matrix(1,nn,nn)

par.init <- c( rep(0,nn), c(A), Sigma[lower.tri(Sigma, diag = TRUE)] )
l <- rest_var_lhood( Y, par.init, a.switch, A.switch, Sigma.switch )
    # c++ code for llhood
Sigma.i <- solve(Sigma)
term <- 0
e.hat <- Y[,-1] - A %*% Y[,-m]
for(i in 1:(m-1)) term <- term + t(e.hat[,i]) %*% Sigma.i %*% e.hat[,i]
ll <- .5 * ( nn * log(2*pi) + log( det( Sigma ) ) + term / (m-1) )

vv <- VAR(t(Y))
    # Estimate the var
coeff <- sapply( c('y1', 'y2'), function(x) vv$varresult[[x]]$coeff)

ww <- optim( par.init, function(par) rest_var_lhood( Y, par, a.switch, A.switch, Sigma.switch ) )
