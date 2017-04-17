library(vars)
library(MASS)
library(microbenchmark)
library(nloptr)
library(lubridate)
library(expm)

## 1. TEST THAT FUNCTION EVALUATION WORKS ##
A <- A.vals <- matrix( c( .5, .2, .6, -.2 ), 2, 2, byrow=TRUE )
Sigma <- Sigma.vals <- matrix( c( 1, .1, .1, .5 ), 2, 2, byrow=TRUE )
a <- a.vals <- c(1,0)
m <- 100000
nn <- 2
# Set up

e <- t( mvrnorm( m, rep(0,2), Sigma ) )
# Errors
Y <- matrix( 0, nn, m )
# Initialize the VAR data
for( i in 2:m ){
  Y[,i] <- a + A %*% Y[,i-1] + e[,i]
}
a.switch <- rep(1,nn)
A.switch <- Sigma.switch <- matrix(1,nn,nn)

par.init <- c( rep(0,nn), c(A), Sigma[lower.tri(Sigma, diag = TRUE)] )
microbenchmark( l <- rest_var_lhood( Y, par.init, a.switch, A.switch, Sigma.switch,
                                      a.vals, A.vals, Sigma.vals ) )
    # c++ code for llhood
Sigma.i <- solve(Sigma)
term <- 0
e.hat <- Y[,-1] - A %*% Y[,-m]
for(i in 1:(m-1)) term <- term + t(e.hat[,i]) %*% Sigma.i %*% e.hat[,i]
ll <- .5 * ( nn * log(2*pi) + log( det( Sigma ) ) + term / (m-1) )


## 2. UNRESTRTICTED ESTIMATION ##
vv <- VAR(t(Y))
    # Estimate the var
coeff <- sapply( c('y1', 'y2'), function(x) vv$varresult[[x]]$coeff)
control <- list(maxit = 2000)
print( system.time( unrest <- optim( par.init,
              function(par) rest_var_lhood( Y, par, a.switch, A.switch, Sigma.switch,
                                           a.vals, A.vals, Sigma.vals ), control = control ) ) )
# opts <- list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8, maxeval=5000 )
# print( system.time(  unrest.nl <- nloptr( par.init,
#               function(par) rest_var_lhood( Y, par, a.switch, A.switch, Sigma.switch,
#                                             a.vals, A.vals, Sigma.vals ), opts=opts ) ) )

a.est <- unrest$par[1:nn]
A.est <- matrix(unrest$par[nn+1:nn^2],nn,nn)
Sigma.est <- matrix(0,nn,nn)
Sigma.est[upper.tri(Sigma.est,diag = TRUE)] <- Sigma.est[lower.tri(Sigma.est,diag = TRUE)] <- unrest$par[-(1:(nn+nn^2))]


## 3. TRY WITH SOME VALUES FIXED ##
A.switch[1,2] <- 0
A.vals[1,2] <- .4
par.init.r <- par.init[-5]
microbenchmark( l.r <- rest_var_lhood( Y, par.init.r, a.switch, A.switch, Sigma.switch,
                                        a.vals, A.vals, Sigma.vals ) )
rest <- optim( par.init.r,
                 function(par) rest_var_lhood( Y, par, a.switch, A.switch, Sigma.switch,
                                               a.vals, A.vals, Sigma.vals ), control = control  )
a.est.r <- rest$par[1:nn]
A.est.r <- A.vals
A.est.r[A.switch==1] <- rest$par[nn+1:(nn^2-1)]
Sigma.est.r <- matrix(0,nn,nn)
Sigma.est.r[upper.tri(Sigma.est.r,diag = TRUE)] <- Sigma.est.r[lower.tri(Sigma.est.r,diag = TRUE)] <- rest$par[-(1:(nn+nn^2-1))]

print( cbind( a, var=coeff['const',], a.est, a.est.r ) )

mu.Y <- apply(Y,1,mean)
mu.var <- solve( (diag(nn)-t(coeff[1:nn,])), coeff['const',] )
mu <- solve(diag(nn)-A, a)
mu.est <- solve(diag(nn)-A.est, a.est)
mu.est.r <- solve(diag(nn)-A.est.r, a.est.r)

print( cbind( mu, mu.Y, mu.var, mu.est, mu.est.r ) )

## 4. Try on the US DATA ##
cty <- 'USA'
    # Country
# cuts <-as.Date(c( '1950/01/01', '1960/01/01', '1970/01/01', '1980/01/01',
#                   '1990/01/01', '2000/01/01', '2009/01/01', '2018/01/01'))
# labs <- c('50s', '60s', '70s', '80s', '90s', '2000s pre-09', 'post-09')
#
# decades.mu <- aggregate( dta$rmg, list(cut(dta$DATE, cuts, labs )), mean, na.rm=TRUE )
dta <- read.csv( paste0( 'data/', cty, '.csv' ) )
dta$Date <- as.Date( dta$Date, "%m/%d/%Y" )
    # Read in and clean the data
dta$rfr4 <- dta$rfr / 4
dta$rmg <- dta$rfr4 - dta$ngdp_pch
    # The quarterly risk free rate
with( dta, plot( ngdp_pch, rfr4 ) )
with( dta, lines( ngdp_pch, rfr4, lwd=.5 ) )
plot(dta$Date, dta$rmg, type='l', lwd=2 )
abline( h=c(-1,1), lwd=.5)
    # Basic plots
us.var <- VAR(dta[,c('ngdp_pch', 'rfr4')])
    # The US VAR
fit <- fitted(us.var)
    # The fitterd US VAR
us.coeff <- sapply( c('ngdp_pch', 'rfr4'), function(x) us.var$varresult[[x]]$coeff)
sigma.us <- m / (m-1) * var(sapply( c('ngdp_pch', 'rfr4'), function(x) us.var$varresult[[x]]$residuals))
par.init.us <- c( us.coeff['const',], t(us.coeff[1:nn,]), sigma.us[lower.tri(sigma.us,TRUE)] )
A.switch <- matrix( 1, nn, nn )
    # Remove constraints
l.us <- rest_var_lhood( t(as.matrix(dta[,c('ngdp_pch','rfr4')])), par.init.us,
                a.switch, A.switch, Sigma.switch,
                a.vals, A.vals, Sigma.vals, 1 )
    # The US likelihood
us.est <- optim( par.init.us,
                 function(par) rest_var_lhood( t(as.matrix(dta[,c('ngdp_pch','rfr4')])), par,
                                               a.switch, A.switch, Sigma.switch,
                                               a.vals, A.vals, Sigma.vals ), control = control )

## 5. WITH SOME SWITCHING ##
n.state <- 3
dta.gr <- dta[,c('ngdp_pch','rfr4')]
set.seed(42)
cl <- kmeans( dta.gr, n.state )
a.0 <- cl$centers
freq <- table( cl$cluster[-(nrow(dta.gr))], cl$cluster[-1] )
P.0 <- freq / apply( freq, 1, sum )
    # Create a guess of the means and transition probabilities
p.0 <- rep(0,n.state)
p.0[cl$cluster[1]] <- 1
    # Initialize the first period
var.state <- VAR( dta.gr - cl$centers[cl$cluster,] )
    # Create a VAR with cluster means removed
msw.coeff.0 <- t(sapply( c('ngdp_pch', 'rfr4'), function(x) var.state$varresult[[x]]$coeff[1:nn]))
msw.sigma.0 <- m / (m-1) * var(sapply( c('ngdp_pch', 'rfr4'), function(x) var.state$varresult[[x]]$residuals) )
msw.a.0 <- t( diag(nn) - msw.coeff.0 ) %*% t( cl$centers )
    # Convert to parameters
par.msw.init <- rbind( msw.a.0,
                       matrix( c( msw.coeff.0, msw.sigma.0[lower.tri(msw.sigma.0,TRUE)] ),
                               nrow=nn^2+.5*nn*(nn+1), ncol=n.state ) )
probs <- msw_rest_var_lhood_p( t(as.matrix(dta.gr)), p.0, P.0, par.msw.init,
                                a.switch, A.switch, Sigma.switch,
                                a.vals, A.vals, Sigma.vals )
    # Compute the filtering and prediction probabilities
rg.msw <- cl$centers[ apply(probs[1:n.state,], 2, which.max), ]
    # Convert to cluster centres

par(mfrow=c(2,1))
plot(dta$Date, dta$ngdp_pch, type='l', lwd=2, main='Nominal GDP' )
lines(dta$Date, cl$centers[cl$cluster,1], col='red' )
lines(dta$Date, rg.msw[,1], col='blue' )
plot(dta$Date, dta$rfr4, type='l', lwd=2, main='Nominal interest rate' )
lines(dta$Date, cl$centers[cl$cluster,2], col='red' )
lines(dta$Date, rg.msw[,2], col='blue' )
par(mfrow=c(1,1))

plot( range(dta$Date), c(0,1), type='n' )
for( i in 1:n.state ) lines( dta$Date, probs[i,], lwd=2, col=i )


## 6. ESTIMATE THE SWITCHING PROCESS ##
ql.sw <- function(par){
# Markov switching quasi-likelihood function
  par.p <- par[1:((n.state+1)*(n.state-1))]
  par.z <- par[-(1:((n.state+1)*(n.state-1)))]
      # Extract the markovian and VAR parts separately
  p.init.0 <- c( par.p[1:(n.state-1)], 1 - sum(par.p[1:(n.state-1)] ) )
  p.init <-  pmin( pmax( p.init.0, 0 ), 1 )
  P.0 <- matrix( par.p[n.state-1 + 1:((n.state-1)*n.state)], n.state, n.state - 1 )
  P.0 <- cbind( P.0, 1-apply(P.0,1,sum) )
  P <- pmin( pmax( P.0, 0 ), 1 )
      # The initial and transition probabilities
  penalty <- sum( 100000 * p.init.0 * ( p.init.0 - 1 ) * ( p.init.0 != p.init ) ) +
                      sum( 100000 * P.0 * ( P.0 - 1 ) * ( P.0 != P ) )
      # Penalty function for violating bounds on p
  a.par <- matrix( par.z[1:(nn*n.state)], nn, n.state )
  A.par <- par.z[nn*n.state + 1:(nn^2)]
  sig.par <- tail( par.z, .5 * nn * (nn+1) )
      # Extract parameters
  m.par <- rbind( msw.a.0, matrix( c( A.par, sig.par ), nrow=nn^2+.5*nn*(nn+1), ncol=n.state ) )
      # Arrange into a matrix
  msw_rest_var_lhood( t(as.matrix(dta.gr)), p.init, P, m.par,
                      a.switch, A.switch, Sigma.switch,
                      a.vals, A.vals, Sigma.vals ) + penalty
}
p.init <- (P.0 %^% 100)[1,]
par.ql <- c(p.init[-n.state], P.0[,-n.state], par.msw.init[1:nn,], par.msw.init[-(1:nn),1] )
ql <- ql.sw(par.ql)
    # Evaluate the quasi likelihood
par.ql.2 <- c(p.init[-n.state] + 1, P.0[,-n.state], par.msw.init[1:nn,], par.msw.init[-(1:nn),1] )
    # Test with violation of bounds
ql.2 <- ql.sw(par.ql.2)
    # Debug to see penalty created
control$maxit <- 20000
    # Number of iterations
# lb <- -Inf * par.ql
# lb[1:((n.state+1)*(n.state-1))] <- 0
#     # Lower bound
# ub <- Inf * par.ql
# ub[1:((n.state+1)*(n.state-1))] <- 1
#     # Upper bound
opt.sw <- optim( par.ql, ql.sw, control = control )
    # Maximize the likelihood
# debug(ql.sw)
ql.sw(opt.sw$par)

p.0.opt <- c( opt.sw$par[1:(n.state-1)], 1-sum(opt.sw$par[1:(n.state-1)]) )
P.opt <- matrix( opt.sw$par[n.state-1 + 1:((n.state-1)*n.state)], n.state, n.state - 1 )
P.opt <- cbind( P.opt, 1-apply(P.opt,1,sum) )
par.z.opt <- opt.sw$par[-(1:((n.state+1)*(n.state-1)))]
a.opt <-  matrix( par.z.opt[1:(nn*n.state)], nn, n.state )
A.opt <-  matrix( par.z.opt[nn*n.state+1:(nn^2)], nn, nn )
Sigma.opt <- matrix(0, nn, nn )
Sigma.opt[lower.tri(Sigma.opt, TRUE)] <- Sigma.opt[upper.tri(Sigma.opt, TRUE)] <- tail( par.z.opt, .5 * nn * (nn+1) )
mu.opt <- solve( diag(nn) - A.opt, a.opt )
    # Convert to matrix form
probs.opt <- msw_rest_var_lhood_p( t(as.matrix(dta.gr)), p.0.opt, P.opt,
                                   rbind( a.opt,
                                          matrix(par.z.opt[-(1:(nn*n.state))],
                                                 nrow=nn^2+.5*(nn+1)*nn, ncol=n.state)),
                                     a.switch, A.switch, Sigma.switch,
                                     a.vals, A.vals, Sigma.vals )
    # The optimal probabilities
rg.mu.opt <- t(mu.opt[ , apply(probs.opt[1:n.state,], 2, which.max) ])
    # Time series of means

plot( range(dta$Date), c(0,1), type='n' )
for( i in 1:n.state ) lines( dta$Date, probs.opt[i,], lwd=2, col=i )
    # Plot the state filtering probabilities

par(mfrow=c(2,1))
plot(dta$Date, dta$ngdp_pch, type='l', lwd=2, main='Nominal GDP' )
lines(dta$Date, rg.mu.opt[,1], col='red' )
plot(dta$Date, dta$rfr4, type='l', lwd=2, main='Nominal interest rate' )
lines(dta$Date, rg.mu.opt[,2], col='red' )
par(mfrow=c(1,1))
