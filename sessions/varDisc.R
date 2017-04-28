## Code to check the discretization of a VAR
library(debtLimits)

dta <- rg.read( 'USA' )
# dta <- subset(dta, date <= "2015-01-01")
est <- var.rg.est(dta)
disc <- var.discretize( est, n.pts = 5, n.dirs=12 )
# disc <- var.discretize( est, n.pts = 1, n.dirs=4 )

n.sim <- 1e6
n.reg <- 1e5
n.plot <- 150
m.sim <- markov_sim( n.sim, disc$trans, which.max(disc$m)-1, length(disc$m))
sim <- disc$X[m.sim+1,]
  # t(ms_var( matrix(est$a,ncol=1), est$A, est$Sigma, rep(0,n.sim), n.sim ))
colnames(sim) <- c('gth', 'rfr')
var.disc <- VAR(sim)

plot( 1:n.plot, sim[ n.sim - n.plot - 1 + 1:n.plot, 'gth'], type='l', lwd=2,
      ylab='Percent change', xlab='pds' )
lines( 1:n.plot, sim[ n.sim - n.plot - 1 + 1:n.plot, 'rfr'], col='red', lwd=2 )
legend('bottom', c('Simulated growth', 'Simulated risk-free rate'), lwd=2, bty='n', col=c(1,2))
interp <- var.data.interp( disc, dta )

plot( disc$X, cex=disc$m*100 )
# lines( interp$dta.disc, col='blue' )

### Compare coefficients of estimated VAR:
est.disc <- list( a = sapply( var.disc$varresult, function(x) x$coeff['const'] ),
                  A = t( sapply( var.disc$varresult, function(x) x$coeff[c('gth.l1','rfr.l1')] ) ),
                  Sigma = var(sapply( var.disc$varresult, function(x) x$residuals ) ) )
est.disc$mu <- solve( diag(nrow(est.disc$A)) - est.disc$A, est.disc$a )
mu <- apply( dta[,c('gth','rfr')], 2, mean )
Sigma.u <- var(dta[,c('gth','rfr')])
est$Sigma.u <- Reduce( '+', lapply( 0:500, function(i) (est$A %^%i) %*% est$Sigma %*% (t(est$A)%^%i) ) )
est.disc$Sigma.u <- var(sim)

print( '### Mean ###' )
print( rbind( mu, est$mu, est.disc$mu ) )
print( '### Unconditional Variance ###' )
print( rbind( c( Sigma.u), c(est$Sigma.u), c(est.disc$Sigma.u ) ) )
print( '### Persistence ###' )
print( rbind( c(est$A), c(est.disc$A ) ) )
print( '### Conditional Variance ###' )
print( rbind( c(est$Sigma), c(est.disc$Sigma ) ) )

