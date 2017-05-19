library(vars)
library(MASS)
library(microbenchmark)
library(nloptr)
library(lubridate)
library(expm)
library(debtLimits)
rm(list=ls())

## 1. TEST THAT FUNCTION EVALUATION WORKS ##
A <- matrix( c( .5, .2, .6, -.2 ), 2, 2, byrow=TRUE )
Sigma <- matrix( c( 1, .1, .1, .5 ), 2, 2, byrow=TRUE )
a <- c(1,0)
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

par.init <- c( rep(0,nn), c(A), Sigma[lower.tri(Sigma, diag = TRUE)] )
microbenchmark( l <- var_lhood( Y, par.init ) )
    # c++ code for llhood
Sigma.i <- solve(Sigma)
term <- 0
e.hat <- Y[,-1] - A %*% Y[,-m]
for(i in 1:(m-1)) term <- term + t(e.hat[,i]) %*% Sigma.i %*% e.hat[,i]
ll <- .5 * ( nn * log(2*pi) + log( det( Sigma ) ) + term / (m-1) )
    ## Check that ll == l

## 2. UNRESTRTICTED ESTIMATION ##
vv <- VAR(t(Y))
    # Estimate the var
coeff <- sapply( c('y1', 'y2'), function(x) vv$varresult[[x]]$coeff)
sig.test <- var(sapply( c('y1', 'y2'), function(x) vv$varresult[[x]]$residuals))
par.init <- c( coeff['const', ], t(coeff[1:nn,]), sig.test[lower.tri(sig.test,TRUE)] )
control <- list(maxit = 2000)
print( system.time( unrest <- optim( par.init, var_lhood, Y=Y, method='BFGS', control = control ) ) )
# opts <- list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8, maxeval=5000 )
# print( system.time(  unrest.nl <- nloptr( par.init,
#               function(par) var_lhood( Y, par, a.switch, A.switch, Sigma.switch,
#                                             a.vals, A.vals, Sigma.vals ), opts=opts ) ) )

a.est <- unrest$par[1:nn]
A.est <- matrix(unrest$par[nn+1:nn^2],nn,nn)
Sigma.est <- matrix(0,nn,nn)
Sigma.est[upper.tri(Sigma.est,diag = TRUE)] <- Sigma.est[lower.tri(Sigma.est,diag = TRUE)] <- unrest$par[-(1:(nn+nn^2))]
print(rbind(t(A.est), a.est) - coeff)
    ## SHould be zero

# ## 3. TRY WITH SOME VALUES FIXED ##
# A.switch[1,2] <- 0
# A.vals[1,2] <- .4
# par.init.r <- par.init[-5]
# microbenchmark( l.r <- var_lhood( Y, par.init.r, a.switch, A.switch, Sigma.switch,
#                                         a.vals, A.vals, Sigma.vals ) )
# rest <- optim( par.init.r,
#                  function(par) var_lhood( Y, par, a.switch, A.switch, Sigma.switch,
#                                                a.vals, A.vals, Sigma.vals ), control = control  )
# a.est.r <- rest$par[1:nn]
# A.est.r <- A.vals
# A.est.r[A.switch==1] <- rest$par[nn+1:(nn^2-1)]
# Sigma.est.r <- matrix(0,nn,nn)
# Sigma.est.r[upper.tri(Sigma.est.r,diag = TRUE)] <- Sigma.est.r[lower.tri(Sigma.est.r,diag = TRUE)] <- rest$par[-(1:(nn+nn^2-1))]

print( cbind( a, var=coeff['const',], a.est ) ) #, a.est.r ) )

mu.Y <- apply(Y,1,mean)
mu.var <- solve( (diag(nn)-t(coeff[1:nn,])), coeff['const',] )
mu <- solve(diag(nn)-A, a)
mu.est <- solve(diag(nn)-A.est, a.est)
# mu.est.r <- solve(diag(nn)-A.est.r, a.est.r)

print( cbind( mu, mu.Y, mu.var, mu.est ) ) #, mu.est.r ) )

## 4. Try on the US DATA ##
cty <- 'USA'
    # Country
# cuts <-as.Date(c( '1950/01/01', '1960/01/01', '1970/01/01', '1980/01/01',
#                   '1990/01/01', '2000/01/01', '2009/01/01', '2018/01/01'))
# labs <- c('50s', '60s', '70s', '80s', '90s', '2000s pre-09', 'post-09')
#
# decades.mu <- aggregate( dta$rmg, list(cut(dta$DATE, cuts, labs )), mean, na.rm=TRUE )
dta <- read.csv( paste0( 'data/', cty, '.csv' ) )
dta$date <- as.Date( dta$date, "%m/%d/%Y" )
    # Read in and clean the data
dta$rfr4 <- dta$rfr / 4
dta$rmg <- dta$rfr4 - dta$ngdp_pch
    # The quarterly risk free rate

par.dft <- par('mar')
pdf('~/Dropbox/2017/research/debtLimits/charts/rg_us.pdf')
par(mfrow=c(2,1), mar=c(3,3,3,3))
plot(dta$date, dta$ngdp_pch, type='l', lwd=2, main='Nominal GDP growth', xlab='', ylab='' )
abline(h=0, lwd=.5)
plot(dta$date, dta$rfr4, type='l', lwd=2, main='Nominal interest rate', xlab='', ylab='' )
dev.off()
par(mfrow=c(1,1), mar=par.dft)

rmg.filter <- 4 * 10

par(mfrow=c(1,1))
with( dta, plot( ngdp_pch, rfr4 ) )
with( dta, lines( ngdp_pch, rfr4, lwd=.5 ) )

pdf('~/Dropbox/2017/research/debtLimits/charts/rmg_us.pdf')
plot(dta$date, dta$rmg, type='l', lwd=2, xlab='', ylab='Interest-growth differential' ) #, main='Interest-growth differential' )
lines(dta$date, filter( dta$rmg, rep(1/rmg.filter,rmg.filter), sides=1), col='blue', lwd=2 )
abline( h=(-2):2, lwd=.5, lty=2)
legend( 'topright', paste0( rmg.filter/4, ' year rolling average'), bty='n', col='blue', lwd=2 )
abline(h=0, lwd=.5)
dev.off()
par(mfrow=c(1,1), mar=par.dft)

pdf('~/Dropbox/2017/research/debtLimits/charts/rmg_us4.pdf')
plot(dta$date, 4 * dta$rmg, type='l', lwd=2, xlab='', ylab='' ) #, main='Interest-growth differential' )
lines(dta$date, filter( 4 * dta$rmg, rep(1/rmg.filter,rmg.filter), sides=1), col='blue', lwd=2 )
abline( h=4*(-2):2, lwd=.5, lty=2)
legend( 'topright', paste0( rmg.filter/4, ' year rolling average'), bty='n', col='blue', lwd=2 )
abline(h=0, lwd=.5)
dev.off()

print( mean(4 * dta$rmg))
print( quantile(4 * abs(diff(dta$rmg)), c(.25, .5, .75, .9)))

us.var <- VAR(dta[,c('ngdp_pch', 'rfr4')])
    # The US VAR

fit <- fitted(us.var)
    # The fitterd US VAR
us.coeff <- sapply( c('ngdp_pch', 'rfr4'), function(x) us.var$varresult[[x]]$coeff)
sigma.us <-var(sapply( c('ngdp_pch', 'rfr4'), function(x) us.var$varresult[[x]]$residuals))
par.init.us <- c( us.coeff['const',], t(us.coeff[1:nn,]), sigma.us[lower.tri(sigma.us,TRUE)] )
    # Remove constraints
l.us <- var_lhood( t(as.matrix(dta[,c('ngdp_pch','rfr4')])), par.init.us, 1 )
    # The US likelihood
control$abstol <- control$reltol <- 1e-08
control$trace <- 1
us.est <- optim( par.init.us, var_lhood, Y=t(as.matrix(dta[,c('ngdp_pch','rfr4')])),
                 method='BFGS', control = control )
a.us.est <- us.est$par[1:nn]
A.us.est <- matrix( us.est$par[nn+1:(nn^2)], nn, nn )
Sigma.us.est <- 0 * A.us.est
Sigma.us.est[ lower.tri(Sigma.us.est, TRUE) ] <-
  Sigma.us.est[ upper.tri(Sigma.us.est, TRUE) ] <- tail(us.est$par, .5*nn*(nn+1))
dta.gr <- dta[,c('ngdp_pch','rfr4')]
us.mu <- apply(dta.gr,2,mean)
us.mu.est <- solve( diag(nn) - A.us.est, a.us.est )
print( cbind( us.mu=us.mu,
              us.mu.var=solve( diag(nn) - t(us.coeff[1:nn,]), us.coeff['const',] ),
              us.mu.est=solve( diag(nn) - A.us.est, a.us.est ) ) )
    ## Means line up.  Nice.
fit <- fitted(us.var)
    # Create fitted values
sd.us.var <- sqrt(diag(Sigma.us.est))
    # The standard deviation
par(mfrow=c(2,1))
plot(dta$date, dta$ngdp_pch, type='l', lwd=2, main='Nominal GDP' )
lines(dta$date[-1], fit[,1], col='blue', lwd=2 )
lines(dta$date[-1], fit[,1] + sd.us.var[1], col='red', lwd=1, lty=1 )
lines(dta$date[-1], fit[,1] - sd.us.var[1], col='red', lwd=1, lty=1 )
plot(dta$date, dta$rfr4, type='l', lwd=2, main='Nominal interest rate' )
# lines(dta$date, cl$centers[cl$cluster,2], col='red' )
lines(dta$date[-1], fit[,2], col='blue', lwd=2 )
lines(dta$date[-1], fit[,2] + sd.us.var[2], col='red', lwd=1, lty=1 )
lines(dta$date[-1], fit[,2] - sd.us.var[2], col='red', lwd=1, lty=1 )
par(mfrow=c(1,1))


n.pred <- 4 * 8
pred.1960 <- t(dta.gr[dta$date == "1960-01-01",])
pred.1970 <- t(dta.gr[dta$date == "1970-01-01",])
pred.1980 <- t(dta.gr[dta$date == "1980-01-01",])
pred.1990 <- t(dta.gr[dta$date == "1990-01-01",])
pred.2000 <- t(dta.gr[dta$date == "2000-01-01",])
pred.2010 <- t(dta.gr[dta$date == "2010-01-01",])
for( i in 2:n.pred ){
  pred.1960 <- cbind(pred.1960, a.us.est + A.us.est %*% pred.1960[,i-1])
  pred.1970 <- cbind(pred.1970, a.us.est + A.us.est %*% pred.1970[,i-1])
  pred.1980 <- cbind(pred.1980, a.us.est + A.us.est %*% pred.1980[,i-1])
  pred.1990 <- cbind(pred.1990, a.us.est + A.us.est %*% pred.1990[,i-1])
  pred.2000 <- cbind(pred.2000, a.us.est + A.us.est %*% pred.2000[,i-1])
  pred.2010 <- cbind(pred.2010, a.us.est + A.us.est %*% pred.2010[,i-1])
}

var.sd <- cbind( c(0,0), sapply( 1:(n.pred-1), function(i)
  sqrt( diag( Reduce( '+', lapply( 0:(i-1),
         function(j) ( A.us.est %^% j ) %*% Sigma.us.est %*% ( t(A.us.est) %^% j ) ) ) ) ) ) )
var.sd.u <- sqrt( diag( Reduce( '+', lapply( 0:200,
                    function(j) ( A.us.est %^% j ) %*% Sigma.us.est %*% ( t(A.us.est) %^% j ) ) ) ) )

# par(mfrow=c(2,1))

pdf('~/Dropbox/2017/research/debtLimits/charts/VAR_US_pred_gdp.pdf', width = 8, height = 3)
par(mar=c(2,2,2,2))
plot(dta$date, dta$ngdp_pch, type='l', lwd=2, xlab='', ylab='' ) # , main='Nominal GDP' )
abline(h=us.mu.est[1], col='blue', lwd=1, lty=2 )
abline(h=0, lwd=.5)
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], pred.1960[1,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], pred.1970[1,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], pred.1980[1,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], pred.1990[1,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], pred.2000[1,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], pred.2010[1,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], pred.1960[1,] + var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], pred.1970[1,] + var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], pred.1980[1,] + var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], pred.1990[1,] + var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], pred.2000[1,] + var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], pred.2010[1,] + var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], pred.1960[1,] - var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], pred.1970[1,] - var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], pred.1980[1,] - var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], pred.1990[1,] - var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], pred.2000[1,] - var.sd[1,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], pred.2010[1,] - var.sd[1,], lwd=1, col='red', lty=2 )
dev.off()

pdf('~/Dropbox/2017/research/debtLimits/charts/VAR_US_pred_r.pdf', width = 8, height = 3)
par(mar=c(2,2,2,2))
plot(dta$date, dta$rfr4, type='l', lwd=2, xlab='', ylab='', ylim=c(-0.5,4.5) ) # main='Nominal interest rate',
abline(h=us.mu.est[2], col='blue', lwd=1, lty=2 )
abline(h=0, lwd=.5)
legend('topright', c('Data', 'VAR forecast', '+/- 1 std dev', 'Unconditional mean'),
       lty=c(1,1,2,2), lwd=c(2,2,1,1), col=c('black', 'blue', 'red', 'blue'), bty='n')
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], pred.1960[2,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], pred.1970[2,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], pred.1980[2,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], pred.1990[2,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], pred.2000[2,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], pred.2010[2,], lwd=2, col='blue' )
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], pred.1960[2,] + var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], pred.1970[2,] + var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], pred.1980[2,] + var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], pred.1990[2,] + var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], pred.2000[2,] + var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], pred.2010[2,] + var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], pred.1960[2,] - var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], pred.1970[2,] - var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], pred.1980[2,] - var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], pred.1990[2,] - var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], pred.2000[2,] - var.sd[2,], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], pred.2010[2,] - var.sd[2,], lwd=1, col='red', lty=2 )
dev.off()

# par(mfrow=c(1,1))

sel.mat <- matrix( c( 1, -1, -1, 1 ), nn, nn )
    # Weights on the variances and covariances in the differenc
rmg.var.sd <- c( 0, sapply( 1:(n.pred), function(i)
                    sqrt( sum( Reduce( '+', lapply( 0:(i-1),
                         function(j) ( A.us.est %^% j ) %*% Sigma.us.est %*% ( t(A.us.est) %^% j ) ) ) * sel.mat ) ) ) )

pdf('~/Dropbox/2017/research/debtLimits/charts/VAR_US_pred_rmg.pdf', width = 8, height = 3)
par(mar=c(2,2,2,2))
plot(dta$date, dta$rfr4 - dta$ngdp_pch, type='l', lwd=2, xlab='', ylab='' ) # , main='R minus G'
abline(h=diff(us.mu.est), col='blue', lwd=1, lty=2 )
abline(h=0, lwd=.5)
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], apply(pred.1960,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], apply(pred.1970,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], apply(pred.1980,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], apply(pred.1990,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], apply(pred.2000,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], apply(pred.2010,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], apply(pred.1960,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], apply(pred.1970,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], apply(pred.1980,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], apply(pred.1990,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], apply(pred.2000,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], apply(pred.2010,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], apply(pred.1960,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], apply(pred.1970,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], apply(pred.1980,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], apply(pred.1990,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], apply(pred.2000,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], apply(pred.2010,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
dev.off()

pdf('~/Dropbox/2017/research/debtLimits/charts/VAR_US_pred_rmg_tall.pdf', width = 8, height = 5)
par(mar=c(2,2,2,2))
plot(dta$date, dta$rfr4 - dta$ngdp_pch, type='l', lwd=2, xlab='', ylab='' ) # , main='R minus G'
abline(h=diff(us.mu.est), col='blue', lwd=1, lty=2 )
abline(h=0, lwd=.5)
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], apply(pred.1960,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], apply(pred.1970,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], apply(pred.1980,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], apply(pred.1990,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], apply(pred.2000,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], apply(pred.2010,2,diff), lwd=3, col='blue' )
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], apply(pred.1960,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], apply(pred.1970,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], apply(pred.1980,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], apply(pred.1990,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], apply(pred.2000,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], apply(pred.2010,2,diff) + rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1960-01-01")+1:n.pred-1], apply(pred.1960,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1970-01-01")+1:n.pred-1], apply(pred.1970,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1980-01-01")+1:n.pred-1], apply(pred.1980,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "1990-01-01")+1:n.pred-1], apply(pred.1990,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2000-01-01")+1:n.pred-1], apply(pred.2000,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
lines( dta$date[which(dta$date == "2010-01-01")+1:n.pred-1], apply(pred.2010,2,diff) - rmg.var.sd[-n.pred], lwd=1, col='red', lty=2 )
dev.off()

res.filter <- 16
# par(mfrow=c(2,1))
pdf('~/Dropbox/2017/research/debtLimits/charts/VAR_US_gdp_res.pdf', width = 8, height = 3)
par(mar=c(2,2,2,2))
plot(dta$date[-1], us.var$varresult$ngdp_pch$residuals, type='l', lwd=1,
     xlab='', ylab='' ) # main='VAR residuals: Nominal GDP',
lines(dta$date[-1], filter( us.var$varresult$ngdp_pch$residuals, rep(1/res.filter,res.filter)), col='blue', lwd=2 )
legend('topright', paste0(res.filter, '-quarter moving average'), lwd=2, col='blue', bty='n' )
abline(h=0)
dev.off()

pdf('~/Dropbox/2017/research/debtLimits/charts/VAR_US_r_res.pdf', width = 8, height = 3)
par(mar=c(2,2,2,2))
plot(dta$date[-1], us.var$varresult$rfr4$residuals, type='l', lwd=1,
     xlab='', ylab='' ) # main='VAR residuals: Nominal interest rate',
lines(dta$date[-1], filter( us.var$varresult$rfr4$residuals, rep(1/res.filter,res.filter)), col='blue', lwd=2 )
abline(h=0)
dev.off()
# par(mfrow=c(1,1))

# acf(us.var$varresult$ngdp_pch$residuals)
# acf(us.var$varresult$rfr4$residuals)

## 5. WITH SOME SWITCHING ##
n.state <- 3
dta.gr <- dta[,c('ngdp_pch','rfr4')]
us.init <- init.var( dta$date, dta.gr, n.state )
    # Initialize the parameters
probs <- msw_var_lhood_p( t(as.matrix(dta.gr)), us.init$p.0, us.init$P, us.init$par )
l <- msw_var_lhood( t(as.matrix(dta.gr)), us.init$p.0, us.init$P, us.init$par )
    # Compute the filtering and prediction probabilities conditional on the states
probs.R <- msw.var.lhood( t(as.matrix(dta.gr)), us.init$p.0, us.init$P, us.init$par )
    # Compare to R version.  This looks great.
rg.msw <- t(us.init$mu[ , apply(probs[1:n.state,], 2, which.max) ])
    # Convert to cluster centres

par(mfrow=c(2,1))
plot(dta$date, dta$ngdp_pch, type='l', lwd=2, main='Nominal GDP' )
lines(dta$date, us.init$mu[1,us.init$cluster], col='red' )
lines(dta$date, rg.msw[,1], col='blue' )
plot(dta$date, dta$rfr4, type='l', lwd=2, main='Nominal interest rate' )
lines(dta$date, us.init$mu[2,us.init$cluster], col='red' )
lines(dta$date, rg.msw[,2], col='blue' )
par(mfrow=c(1,1))

plot( range(dta$date), c(0,1), type='n' )
for( i in 1:n.state ) lines( dta$date, probs[i,], lwd=2, col=i )


## 6. ESTIMATE THE SWITCHING PROCESS ##
ql <- ql.sw.mu( us.init$g.par, dta.gr, n.state, nn )
control$maxit <- 20000
    # Number of iterations
p.0 <- (us.init$P %^% 100)[1,]
g.par <- us.init$g.par
g.par[1:2] <- p.0[1:2]
g.par[(n.state+1)*(n.state-1)+1:(nn*n.state)] <- ( diag(nn) - us.init$A ) %*% us.mu
    # Playing
opt.sw <- optim( g.par, ql.sw.mu, dta=dta.gr, n.state=n.state, nn=nn,
                 method='BFGS', control = control )
    # Maximize the likelihood
# debug(ql.sw.mu)
opt.par <- ql.sw.mu( opt.sw$par, dta.gr, n.state, nn, TRUE )
# undebug(ql.sw.mu)


probs.opt <- msw_var_lhood_p( t(as.matrix(dta.gr)), opt.par$p.0, opt.par$P,
                              opt.par$m.par )
    # The optimal probabilities
rg.mu.opt <- t(opt.par$mu[ , apply(probs.opt[1:n.state,], 2, which.max) ])
    # Time series of means

plot( range(dta$date), c(0,1), type='n' )
for( i in 1:n.state ) lines( dta$date, probs.opt[i,], lwd=2, col=i )
    # Plot the state filtering probabilities

par(mfrow=c(2,1))
plot(dta$date, dta$ngdp_pch, type='l', lwd=2, main='Nominal GDP' )
lines(dta$date, rg.mu.opt[,1], col='red' )
plot(dta$date, dta$rfr4, type='l', lwd=2, main='Nominal interest rate' )
lines(dta$date, rg.mu.opt[,2], col='red' )
par(mfrow=c(1,1))

opt.sw.rs <- ms.est.random.starts( dta$date,  as.matrix(dta.gr), n.state, n.start = 10 )

par(mfrow=c(2,1))
plot(dta$date, dta$ngdp_pch, type='l', lwd=2, main='Nominal GDP' )
lines(dta$date, opt.sw.rs$mu[1,opt.sw.rs$modal.state], col='red' )
plot(dta$date, dta$rfr4, type='l', lwd=2, main='Nominal interest rate' )
lines(dta$date, opt.sw.rs$mu[2,opt.sw.rs$modal.state], col='red' )
par(mfrow=c(1,1))
