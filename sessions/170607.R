#### LEARNING POINT: VARIANCE OF R-G STILL TOO HIGH. NEED TO ALLOW FOR AVERAGE MATURITY TO CHANGE ####

rm(list=ls())
library(VARext)
library(debtLimits)

hist <- TRUE
start.year <- 1880

if(hist){
  load('data/rmg_hist__1880_est.rdata')
}else{
  load('data/rmg_est.rdata')
}
    # Load the VAR estimates

#### Set up the stochastic process ####
n.pts <- 2 ; n.dirs <- 4 # (2,5) works
cty <- 'USA' ; lags <- 1
disc <- var.disc(l.cty.mle.rest.1$est[[cty]], n.pts, n.dirs, lb=c(-Inf,0) )

## Try converting to multi-period VAR ##
# maty <- 4
# A.multi <- l.cty.mle.rest.1$est[[cty]]$A %^% maty
# a.multi <- ( diag(2) - A.multi ) %*% l.cty.mle.rest.1$est[[cty]]$mu
# Simga.multi <- 1 / maty * Reduce( '+', lapply( (1:maty)-1,
#                 function(i) (l.cty.mle.rest.1$est[[cty]]$A) %^% i %*%
#                       l.cty.mle.rest.1$est[[cty]]$Sigma %*% t(l.cty.mle.rest.1$est[[cty]]$A %^%i ) ) )
# disc <- var.disc(list(a=a.multi,A=A.multi,Sigma=Simga.multi,mu=mu.calc(a.multi,A.multi)),
#                  n.pts, n.dirs, lb=c(-Inf,0) )

disc <- var.disc(l.cty.mle.rest.4$est[[cty]], n.pts, n.dirs, lb=c(-Inf,0) )
    # The dscretized VAR
set.seed(42)
m.idx <- markov_sim( 1e6, disc$trans, 0, disc$n.X )[-(1:1e5)]
m.sim <- disc$X[m.idx+1, ]
    # Create the simulation
l.var.est <- var.ols( t(m.sim), lags )
    # Check thsurp( d, x[1:5], x[-(1:5)][s.idx], TRUE))e VAR discretization
var.compare <- sapply( list( l.var.est, l.cty.mle.rest.1$est[[cty]] ),
                       function(x) c( a=x$a, mu=x$mu,A=x$A,
                                      Sigma=x$Sigma[lower.tri(x$Sigma,T)]))
var.compare <- rbind( var.compare, cbind( l.var.est$mu, apply(m.sim,2,mean) ) )
var.compare <- cbind( var.compare, apply( var.compare, 1, function(x) abs(diff(x)) ) )
colnames(var.compare) <- c('data', 'disc', 'diff' )
print(var.compare)
    # The comparison to thte data
# pdf('~/Dropbox/2017/research/debtLimits/charts/rg_us_disc.pdf')
plot(disc$X, cex=disc$p.lr * 20, pch=16, ylab='Risk-free rate', xlab='Nominal growth rate')
abline( 0, 1, lty=2 )
# dev.off()
    # Plot
mar.dft <- par('mar') ; par(mfrow=c(4,4), mar=c(2,2,2,2))
for( i in 1:disc$n.X ) barplot(disc$trans[i,], main=paste0( 'state ', i) )
par(mfrow=c(1,1), mar=mar.dft)


#### Create the surplus process ####
params <- list( R=1+disc$X[,2]/100, G=1+disc$X[,1]/100,
                tri=TRUE, trans=disc$trans )


#### Finish the parameter set-up ####
params$surp.sd <- 1 ; params$lambda <- 0 ; # 0 ;
params$phi <- .8 # .8 # .4
params$cont.type <- 'avg' ; params$diff.method <- "ana" ; params$inner.method <- 'all'
params$d.tri <- TRUE # FALSE # TRUE      # Triangulare distribution for surplus shocks
An <- 1 / params$R ; Bn <- rep( -1, length(params$R) ) ; Cn <- An ; def <- matrix(0,1,1)
params$it <- 15 ; params$tol <- .1
# d.init <- sol.nonstoch( params ) # rep( 600, disc$n.X )

### Compute the surplus function from the risk-free solution ###
tgt <- target.moms('USA')
    ### NEED TO UPDATE THE TARGET MOMENTS FOR EACH COUNTRY ###

rg.dta <- if(hist) hist.read(cty,start.year) else rg.read( cty )
interp <- disc.data.interp( disc, filter(rg.dta[,c('gth','rfr')],rep(1/maty,maty), sides = 1)[-(1:(maty-1)),],
                            lags, v.date = rg.dta$date[-(1:(maty-1))],
                            mains=c('Nominal growth rate','Risk free rate',
                                    'Interest-growth differential'))
    # Create the parameters for the surplus function consistent with no risk of
    # default (so solution is guaraneteed).
# params$v.s.coeff <- rf.surp$v.s.coeff
# params$s.shift <- rf.surp$s.shift
# params$surp.sd <- rf.surp$surp.sd
if(hist){
  # params$v.s.coeff <- c(309/4, 334/4, 500/4, -3.88, 20)
  ave.s <- mean(rg.dta$pb[],na.rm=TRUE )
  params$v.s.coeff <- c(60, 150, 350, 0,0)
}else{
  params$v.s.coeff <- c(309, 334, 500, -3.88, 22) # c(324.23, 327.4, 750, -2.88, 20.5)
}

    # c(324.23, 327.4, 645, -2.88, 20)
params$s.shift <- rep(0,disc$n.X) # c( 0.55, 0.09, -1.8, -1, 1.41, 0.08, -0.34, 0.57, 0.94, -0.63, 0.87 )
params$surp.sd <- .25 # 2.3657
if( hist){
  v.surp <- rg.dta$pb
  pb.alt <- (rg.dta$revenue - rg.dta$expenditure) / rg.dta$gdp *100
      #  Ugly join. Should try to improve...
  v.surp[is.na(v.surp)] <- pb.alt[is.na(v.surp)]
  rf.surp <- rf.fit( params, rg.dta, rf.debt = params$surp( d, x[1:5], x[-(1:5)][s.idx], TRUE))surp( d, x[1:5], x[-(1:5)][s.idx], TRUE))oeff[3],
                     rf.level=params$v.s.coeff[5], print.level=1,
                     init = c( params$v.s.coeff[c(1,2,4)], params$s.shift[-1] ),
                     v.surp=v.surp, d.init=rg.dta$d[1] )
}else{
  rf.surp <- rf.fit( params, tgt$dta, rf.debt = params$v.s.coeff[3],
                     rf.level=params$v.s.coeff[5], print.level=1,
                     init = c( params$v.s.coeff[c(1,2,4)], params$s.shift[-1] ) )
}

params$v.s.coeff <- rf.surp$v.s.coeff
s.idx.p <- sapply(1:disc$n.X, function(i) sum( interp$s.idx == i ) ) / nrow(interp$s.idx)
params$s.shift <- rf.surp$s.shift
params$s.shift[s.idx.p==0] <- 0
params$surp.sd <- rf.surp$aad # rf.surp$surp.sd
x.max <- 500
params$d.pts <- 25 ; params$x.mult <- 3


# sol <- sol.wrapper(params, cbind(sol.s$p,sol.s$d), plot.on = FALSE )
params$maxit <- 50
# sol.s <- sol.search( params )

sol <- sol.wrapper(params, cbind(rf.surp$p,rf.surp$d), plot.on = FALSE )
plot.z( sol$p, sol$d, params, xlim=c(0,max(sol$p)*1.2), ylim=c(0,max(sol$p)*1.2) )

# for( b3 in params$v.s.coeff[3] + 1:20 * 10 ){
#   params$v.s.coeff[3] <- b3
#   message('**************************************')
#   message('* params$v.s.coeff[3] = ', params$v.s.coeff[3] )
#   message('**************************************')
#   sol.new <- sol.wrapper(params, cbind(sol$p,sol$d), plot.on = FALSE )
#   plot.z( sol$p, sol$d, params, xlim=c(0,max(sol$p)*1.2), ylim=c(0,max(sol$p)*1.2) )
#   if(max(abs(sol.new$err)) < 1e-06){
#     sol <- sol.new
#   }else{
#     break()
#   }
# }

# params$v.s.coeff[3] <- 840
# sol <- sol.wrapper(params, cbind(rf.surp$p,rf.surp$d), plot.on = FALSE )
# plot.z( sol$p, sol$d, params, xlim=c(0,max(sol$p)*1.2), ylim=c(0,max(sol$p)*1.2) )
# par(mfrow=c(1,1))

par(mfrow=c(1,1))
plot.surp( params, TRUE, c(0,x.max), TRUE ) #, ylim=c(-15,30) )
# points( tgt$dta$cnlb_gdp_lag[-1], rf.surp$v.surp, pch=16, cex=.75, col=interp$s.idx )
points( rg.dta$d, rf.surp$v.surp, pch=16, cex=.75, col=interp$s.idx )
plot( rg.dta$date, v.surp, type='l', lwd=2, col=2)
lines( rg.dta$date, params$s.shift[interp$s.idx], lwd=2 )
abline(h=0)
# points( tgt$dta$cnlb_gdp_lag, sim[,'s'], pch=16, cex=.75, col=interp$s.idx )
# sol <- sol.wrapper(params, cbind(sol$p,sol$d), plot.on = FALSE )
# sol <- sol.wrapper(params, cbind(rf.surp$p,rf.surp$d), plot.on = FALSE )
# sol <- sol.wrapper(params, plot.on = FALSE )

# d.init <- sol.nonstoch(params)
# d.init[d.init==Inf] <- max(d.init[d.init<Inf])
# sol <- sol.wrapper(params, plot.on = FALSE )
plot.z( sol$p, sol$d, params, xlim=c(0,max(sol$p)*1.2), ylim=c(0,max(sol$p)*1.2) )
# sol <- sol.wrapper(params, cbind(sol.s$p,sol.s$d), plot.on = FALSE )

# stop()

# params.hi <- params
# sol.hi <- sol
# for( inc in seq(.005,.5/4,length.out=10)){
#   message( "inc = ", inc )
#   params.hi$R <- params$R + inc / 100
#   sol.hi <-  sol.wrapper(params.hi, cbind(sol.hi$p,sol.hi$d), plot.on = FALSE )
#     # Equivalent to shift of half a percentage point or so
# }

params.hi <- params
inc <- -.01
params.hi$R <- params$R + inc / 100
sol.hi <-  sol.wrapper(params.hi, cbind(sol$p,sol$d), plot.on = FALSE )


stop()

chg <- -diff(l.thresh.1$USA)
elast.hi <- ( sol$d + ( sol$d - sol.hi$d ) * chg / inc )

pdf('~/Dropbox/2017/research/debtLimits/charts/ts_limit.pdf')
plot( tgt$dta$date[-1], sol$d[interp$s.idx] / 4, type='l',
      lwd=2, ylim=range( c(sol$d, elast.hi) )/4,
      xlab='Date', ylab='Debt limit' )
lines( tgt$dta$date[-1],
       elast.hi[interp$s.idx] / 4,
       lwd=2, col='blue' )
legend('bottomright', c('LR R-G at 5% critical value', 'LR R-G at 2.5% critical value'),
       bty='n', lwd=2, col=c('black','blue') )
dev.off()

pdf('~/Dropbox/2017/research/debtLimits/charts/surp_int.pdf')
plot( tgt$dta$date[-1], sol$params$s.shift[interp$s.idx],
      type='l', lwd=2, xlab='Date', ylab='Surplus intercept' )
abline(h=0,lwd=.5)
dev.off()

plot( (params$R - 1)*100, sol$d / 4, xlab='Nominal interest rate', ylab='Debt limit' )
plot( (params$R - params$G)*100, sol$d / 4, xlab='Interest-growth differential',
      ylab='Debt limit' )
plot( (params$G - 1)*100, params$s.shift, xlab='Nominal growth rate', ylab='Surplus intercept' )
plot( (params$R - 1)*100, params$s.shift, xlab='Nominal interest rate', ylab='Surplus intercept' )
plot( (params$R - params$G)*100, params$s.shift, xlab='Interest-growth differential', ylab='Surplus intercept' )
plot( params$s.shift, sol$d / 4, xlab='Surplus intercept', ylab='Debt limit' )

save( tgt, sol, interp, params, elast.hi, file='sol_jul2017.rdata' )

### Compute a discretized simulation
n.sim <- 1e6
n.plot <- 150
n.reg <- 1e6
m.sim <- markov_sim( n.sim, disc$trans, which.max(disc$p.lr)-1, length(disc$p.lr))
sim <- disc$X[m.sim+1,]
colnames(sim) <- c('gth', 'rfr')
varnames <- c( 'Growth', 'Int. rate' )

mu.sim <- apply( sim, 2, mean )
l.var.ols.disc <- var.ols( t(sim), 1 )
l.cty.mle.rest.1$est[[cty]]
var.table( list(data=l.cty.mle.rest.1$est[[cty]], disc=l.var.ols.disc),
           file=paste0('~/Dropbox/2017/research/debtLimits/tables/', cty, '_var_tab_disc.tex'),
           specnames = c( 'Data', 'Simulation' ), varnames = varnames,
           caption = 'Estimated restricted VAR based on data, and VAR estimated from discretized approximation' ,
           label=paste0('tab:', cty, '_var_tab_disc'), footer=TRUE )

### Create the IRFs ###
n.irf <- 20
wt.fn <- function( v, a, b, c ){
# Solves for constants s,t such that v = t*a + s*b + (1-s-t)*c
# Where a, b, c, and v are vectors
  wts <- solve( cbind( a-c, b-c), v-c )
  return( c( wts[1], wts[2], 1-sum(wts) ) )
}
lrv <- var_lr_variance( l.cty.mle.rest.1$est[[cty]]$A, l.cty.mle.rest.1$est[[cty]]$Sigma )
lr.cor <- lrv[1,2] / sqrt( prod( diag( lrv ) ) )
    # The long-run correlation
v.targ <- list( R=disc$X[1,]-c(0,.1), G=disc$X[1,]-c(.1,0), both=disc$X[1,] - c(lr.cor,.1) )
# v.targ <- list( R=disc$X[9,]-c(0,.1), G=disc$X[9,]-c(.1,0), both=disc$X[9,] - .1 * c(lr.cor,1) )
    # Shocks to R, G and both
l.cand.vert <- list( R=c(1,7,16), G=c(1,10,16), both=c(1,9,2) )
# l.cand.vert <- list( R=c(9,2,3), G=c(9,4,5), both=c(9,4,5) )
    # The candidate vertices to put weight on
# shk <- list( R=c(0,1), G=c(1,0), both=1*c(lr.cor,1) )
m.Sigma <- l.cty.mle.rest.1$est[[cty]]$Sigma
cor.innov <- ( m.Sigma[1,2] / sqrt( prod( diag( m.Sigma ) ) ) ) * sqrt( m.Sigma[1,1] )
shk <- list( R=c(0,sqrt(m.Sigma[2,2])), G=c(sqrt(m.Sigma[1,1]),0),
             both=c( cor.innov, sqrt(m.Sigma[2,2])) )
shk.dist <- lapply( shk, function(e) e / apply(disc$X,2,sum) )
    # Need to figure out how to do a shock to the distribution here
p.init <- lapply( shk, function(x) optim( disc$p.lr[-1],
               function(p) sum( ( c(  1-sum(p), p ) %*% disc$X - x - disc$p.lr %*% disc$X ) ^ 2 ) )$par )
p.init <- lapply( p.init, function(p) c(1-sum(p),p) )
l.irf <- list()
for( j in 1:3 ){
  irf.p <- matrix( 0, nrow=n.irf+1, ncol=nrow(disc$X) )
  # irf.p[1,1] <- 1
  # verts <- l.cand.vert[[j]]
  # ww <- wt.fn( v.targ[[j]], disc$X[verts[1],], disc$X[verts[2],], disc$X[verts[3],] )
  irf.p[ 1, ] <- disc$p.lr
      # Initialize the state distribution
  irf.p[ 2, ] <- p.init[[j]]
      # Initialize the shock
  # print( ww )
  for( i in 2:n.irf){
    irf.p[i+1,] <- irf.p[i,] %*% disc$trans
  }
      # Markov chain evolution
  l.irf[[j]] <- irf.p %*% cbind( disc$X * 4, apply(disc$X,1,diff) * 4, sol$d / 4, params$s.shift ) -
      rep( 1, n.irf+1 ) %*% ( irf.p[1,] %*% cbind( disc$X * 4, apply(disc$X,1,diff) * 4, sol$d / 4, params$s.shift ) )
      # The IRF is an average across states
}
v.col <- c('blue', 'red', 'darkgreen')
v.lty <- c( 1, 2, 3)
v.nm <- c('G','R','RmG','limit','surp')
v.l <- c( 'topright', 'topright', 'bottomright', 'bottomright', 'bottomright' )
for( j in 1:5 ){
  X <- 0:n.irf
  Y <- do.call( cbind, lapply(l.irf, function(x) x[,j]) )
  pdf(paste0('~/Dropbox/2017/research/debtLimits/charts/irf_', cty, '_', v.nm[j], '.pdf' ) )
    plot( range(X), range(Y), type='n', xlab='Quarters', ylab='Response')
    # for( i in 1:3 ) lines( X, Y[,i], lwd=2, col=v.col[i] )
    for( i in 1:2 ) lines( X, Y[,i], lwd=2, col=v.col[i], lty=v.lty[i] )
    abline(h=0,lwd=.5)
    # legend(v.l[j], c('R shock', 'G shock', 'Correlated shock'), lwd=2, col=v.col, bty='n')
    legend(v.l[j], c('R shock', 'G shock'), lwd=2, col=v.col, bty='n')
  dev.off()
}

stop()

### TO DO ###
# - Shut down uncertainty

params.0 <- params
params.0$R <- mean(params$R) # rep( mean(params$R), length(params$R) )
params.0$G <- mean(params$G) # rep( mean(params$G), length(params$G) )
params.0$s.shift <- 0 #
params.0$trans <- matrix(1)
params.0$it <- 150
sol.0 <- sol.wrapper( params.0, cbind( 0, max( sol$d ) ), plot.on = TRUE )
plot.z( sol.0$p, sol.0$d, params.0, xlim=c(0,max(sol.0$p)*1.2), ylim=c(0,max(sol.0$p)*1.2) )

stop()

sol.o <- outer.wrapper( sol, params, print_level = 1 )

#### Create the realized shocks ####
mom.series.plot('USA')
v.surp <- surp.fit(sol.o,interp,params,tgt$dta)
sim <- sim_core( c(1,interp$s.idx), sol.o$d.bar, sol.o$d.grid, sol.o$P, sol.o$Q, params,
                 tgt$dta$cnlb_gdp_lag[1], TRUE, c(0,0,v.surp) )
surp.ee.sim <- mean( sim[-1,'eps'] )
surp.sd.sim <- sd( sim[-1,'eps'] )
    # Mean and sd of surplus shock
    # Create the fitted surpluses

sol.err <- outer.err( sol.o, params )
plot.sol( sol.o, xlim=c(min(sol$d)-15,max(sol$d)+15) )
# plot.err( sol.err )

stop()


est <- price.diff.min( params, interp, tgt$dta, sol, sol.o, h.0=c(5,-.5), maxit = 100 )
# est <- price.diff.min( est$sol$params, interp, tgt$dta, est$sol, est$sol.o,
#                        h.0=c(5,-.5), maxit = 100 )
plot.sol(est$sol.o, xlim=c(min(est$sol$d)-15,max(est$sol$d)+15))
sim.est <- sim_core( c(1,interp$s.idx), est$sol.o$d.bar, est$sol.o$d.grid, est$sol.o$P,
                     est$sol.o$Q, est$sol$params, tgt$dta$cnlb_gdp_lag[1], TRUE, c(0,0,v.surp) )

# est <- price.diff.min( params, interp, tgt$dta, maxit = 5 )
stop()


plot( tgt$dta$pb_gdp, sim[,'s'], pch=16, col=interp$s.idx,
      xlab='Data', ylab='Model', main='Primary balance')
abline(0,1,lty=2)

plot( tgt$dta$date, sim[,'d.prime'], type='l', lwd=2, ylim=c(0,max(sol.o$d.bar)) )
lines( tgt$dta$date, tgt$dta$cnlb_gdp_lag, col='blue', lwd=2 )
lines( tgt$dta$date[-1], sim[-1,'d.bar'], col='red', lwd=2 )

plot( tgt$dta$date, tgt$dta$pb_gdp, type='l', lwd=2, ylim=range(c(tgt$dta$pb_gdp,sim[,'s'])) )
lines( tgt$dta$date, sim[,'s'], col='red', lwd=2 )
abline(h=0)

### Create the long debt prices ###
params.lr.1 <- params ; params.lr.1$lambda <- .75 ; params.lr.1$x.mult <- 5
sol.lr.1 <- outer.wrapper( sol, params.lr.1 )
plot.sol(sol.lr.1, xlim=c(min(sol$d)-40,max(sol$d)+10))
ytm <- sol.o$Q ^ - ( 1 - params$lambda )
Q.1 <- q_hat_mat( sol.o$P, sol.o$d.bar, sol.o$Q, sol.o$Q, sol.o$d.grid, params.lr.1$R, params.lr.1$G,
                   params.lr.1$s.shift, params.lr.1$lambda, params.lr.1$phi, sol.o$e.grid,
                   params.lr.1$v.s.coeff, params.lr.1$tri, matrix(0), FALSE, params.lr.1$trans, 2 )
ytm.1 <- Q.1 ^ - ( 1 - params.lr.1$lambda )
    # One-year debt
params.lr.5 <- params ; params.lr.5$lambda <- .95
# sol.lr.5 <- outer.wrapper( sol, params.lr.5 )
# plot.sol(sol.lr.5)
Q.5 <- q_hat_mat( sol.o$P, sol.o$d.bar, sol.o$Q, sol.o$Q, sol.o$d.grid, params.lr.5$R, params.lr.5$G,
                  params.lr.5$s.shift, params.lr.5$lambda, params.lr.1$phi, sol.o$e.grid,
                  params.lr.5$v.s.coeff, params.lr.5$tri, matrix(0), FALSE, params.lr.5$trans, 2 )
ytm.5 <- Q.5 ^ - ( 1 - params.lr.5$lambda )
    # Five-year debt
for( i in 1:disc$n.X){
  plot( range(sol.o$d.grid), range(c(0,spd[i,][spd[i,]<Inf],
                                     spd.1[i,][spd.1[i,]<Inf],
                                     spd.5[i,][spd.5[i,]<Inf])), type='n',
        lwd=2, xlim=sol$d[i]+c(-8,2), main=paste0("i=",i) )
  lines( sol.o$d.grid, spd[i,], type='l', lwd=2, col='black' )
  lines( sol.o$d.grid, spd.1[i,], type='l', lwd=2, col='red' )
  lines( sol.o$d.grid, spd.5[i,], type='l', lwd=2, col='blue' )
  abline(v=sol.o$d.grid, lwd=.5)
}
tp.1 <- ytm.1 - ytm
tp.5 <- ytm.5 - ytm
for( i in 1:disc$n.X){
  plot( range(sol.o$d.grid), range(c(0, tp.1[i,][tp.1[i,]<Inf],
                                     tp.5[i,][tp.5[i,]<Inf]), na.rm = TRUE ), type='n',
        lwd=2, xlim=sol$d[i]+c(-8,2), main=paste0("i=",i) )
  lines( sol.o$d.grid, tp.1[i,], type='l', lwd=2, col='red' )
  lines( sol.o$d.grid, tp.5[i,], type='l', lwd=2, col='blue' )
  abline(v=sol.o$d.grid, lwd=.5)
}
tp.1.sim <- apply( sim, 1, function(x) approx( sol.o$d.grid, tp.1[x['idx'],],
                                              x['d.prime'] )$y )
tp.5.sim <- apply( sim, 1, function(x) approx( sol.o$d.grid, tp.5[x['idx'],],
                                              x['d.prime'] )$y )
tp <- data.frame( date=tgt$dta$date, yr.1=tgt$dta$Int_1y - tgt$dta$Int_3M,
                  yr.1.sim=100*tp.1.sim, yr.5=tgt$dta$Int_5y - tgt$dta$Int_3M,
                  yr.5.sim=100*tp.5.sim, yr.5.1=tgt$dta$Int_5y - tgt$dta$Int_1y,
                  yr.5.1.sim=100*(tp.5.sim-tp.1.sim))
params.l <- params
params.l$v.s.coeff


### THIS ORDERING LOOKS WRONG ###
