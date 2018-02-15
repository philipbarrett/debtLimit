### Trying to make progress with a calibrated solution for the UK ###

### 1. LOAD LIBRARIES ###
rm(list=ls())
library(VARext)
library(debtLimits)
library(zoo)
library(plyr)
filter <- stats::filter

### 2. LOAD DATA ###
cty <- 'UK' # 'USA' # 'FRA' #
start.year <- 1880
rg.dta <- hist.read(cty,start.year)
rmg.mu <- .16 # .46  # NB: This is for the one-lag case
n.pts <- 3 ; n.dirs <- 10

### 3. DATA WRANGLING ###

## 3.1 Fill in missing data (a little approximate) ##
lag.match <- 10                 # Number of preceding periods to use for splicing average
rg.dta$ob.jt <- with( rg.dta, 100*(revenue-expenditure)/gdp )
jt.maro.pb.diff <- filter( rg.dta$ob.jt - rg.dta$pb,
                           rep( 1/lag.match, lag.match ), sides = 1 )
jt.maro.pb.diff <- c( rep(NA,lag.match-1),na.locf(jt.maro.pb.diff))
rg.dta$pb[is.na(rg.dta$pb)] <- rg.dta$ob.jt[is.na(rg.dta$pb)] - jt.maro.pb.diff[is.na(rg.dta$pb)]
    # Carry latest value across NAs
rg.dta$d[is.na(rg.dta$d)] <- rg.dta$debtgdp[is.na(rg.dta$d)] * 100
    # The two data sources line up almost exactly for this series
rg.dta$d.lag <- c(NA, rg.dta$d[-nrow(rg.dta)])

## 3.2 Compress to interval-averages ##
maty <- 4                       # Period length
rg.dta$pd <- rep(1:ceiling(nrow(rg.dta)/maty), each=maty)[1:nrow(rg.dta)]
rg.dta.maty <- ddply( rg.dta, ~pd, summarise, yr.start=min(year), yr.end=max(year),
                      year=min(year), pb=mean(pb,na.rm=TRUE), rfr=mean(rfr,na.rm=TRUE),
                      gth=mean(gth,na.rm=TRUE), rmg=mean(rfr-gth,na.rm=TRUE),
                      d.end=tail(d,1), d=tail(d,1), d.lag=tail(d.lag,1),
                      ie=mean(ie,na.rm=TRUE), ob.jt=mean(ob.jt,na.rm=TRUE),
                      debtgdp=tail(debtgdp,1) )
    # The essential data for estimation + plotting, at maty-year averages

### 4. CHART TO CHECK EVERYTHING IS OK ###
for( dta in list(rg.dta, rg.dta.maty) ){
  par(mfrow=c(3,1))
  with(dta, plot(year, pb, type='l', lwd=1, main='Fiscal Balance',  xlab='',ylab='% GDP',
                    ylim=range(c(pb,pb-ie,ob.jt),na.rm = TRUE)) )
  with(dta, lines(year, pb-ie, type='l', lwd=1, col='red') )
  with(dta, lines(year, ob.jt, type='l', lwd=1, col='darkgreen') )
  abline(h=0)
  legend( 'bottomleft', c('Primary: Mauro+', 'Overall: Mauro', 'Overall: Jorda-Taylor'), lwd=1, col=c('black','red','darkgreen'), bty='n')
  with(dta, plot(year, d, type='l', lwd=1, main='Government debt',  xlab='',ylab='% GDP') )
  with(dta, lines(year, 100*debtgdp, type='l', lwd=1, col='blue' ) )
  legend( 'topleft', c('Mauro', 'Jorda-Taylor'), lwd=1, col=c('black','blue'), bty='n')
  abline(h=0)
  with(dta, plot(year, gth, type='l', lwd=1, main='Nominal growth and interest rates',
                    ylab='', xlab='', col='black', ylim=range(c(gth,-gth),na.rm = TRUE) ) )
  with(dta, lines(year, rfr, type='l', lwd=1, col='red') )
  with(dta, lines(year, rfr - gth, type='l', lwd=1, col='darkgreen') )
  abline(h=0)
  legend( 'bottomleft', c('G', 'R', 'R-G'), lwd=1, col=c('black','red','darkgreen'), bty='n')
  par(mfrow=c(1,1))
}

### 5. ESTIMATION ###

## 5.1 The VAR ##
# Data
sim <- t(rg.dta[,c('gth','rfr')])
sim.maty <- t(rg.dta.maty[,c('gth','rfr')])

# Unrestricted
l.var.ols <- var.ols( sim )
# l.var.ols.postww2 <- var.ols( sim[,rg.dta$year > 1945] )
l.var.ols.maty <- var.ols( sim.maty )
# Restiction function
g <- function(par, lags, n.var=2, theta=0){
  l.var <- par.from(par, n.var, lags)
  mu <- mu.calc(l.var$a, l.var$A)
  return( - mu[2] + mu[1] + theta )
}
# Chart
v.theta <- seq(diff(l.var.ols$mu),max(1,1.5*rmg.mu),length.out = 20)
lr.wald.plot( sim, 1, g, v.theta, xlab='Mean R-G' )
abline(v=c(0,rmg.mu))
lr.wald.plot( sim.maty, 1, g, v.theta, xlab='Mean R-G' )
abline(v=c(0,rmg.mu))

# Restricted estimates version
l.var.rest <- var.mle.rest( sim, g, 1, theta=rmg.mu )
l.var.rest.2 <- list(  )
l.var.rest.maty <- var.mle.rest( sim.maty, g, 1, theta=rmg.mu )
    # Unrestricted, restricted, maty, non-maty

## 5.2 Discretization ##
# X <- var.disc.pts(l.var.rest.maty, 1, n.pts, n.dirs )
X <- var.disc.pts(l.var.rest, 1, n.pts, n.dirs )
disc <- var.disc.1(l.var.rest, X )
disc$p.lr <- (disc$M %^% 100)[1,]
disc.maty <- var.disc(l.var.rest.maty, n.pts, n.dirs, lb=c(-Inf,0) )
mu.dist <- abs(disc$X - rep(1,nrow(disc$X)) %*% t(l.var.rest$mu))
    # Distance from mean
l.sym <- apply( mu.dist, 1, function(x) which( apply( mu.dist, 1, function(y) all.equal( y, x ) ) == TRUE ) )
sym <- c( 1, unlist(sapply( 1:length(l.sym), function(i) l.sym[[i]][l.sym[[i]] != i] )) )
    # The list of rows wihch should be symmetric (assuming the center is in entry #1)
for( i in 1:nrow(disc$X)){
  if( sum(abs(disc$cm.err[,i])) > sum(abs(disc$cm.err[,sym[i]]))){
    disc$M[i,] <- disc$M[sym[i],][sym]
  }
}
# disc$M[1,] <- .5 * ( disc$M[1,] + disc$M[1,][sym] )
    # Enforce symmetry of the discretized matrix
# X <- var.disc.pts(l.var.rest, 1, n.pts, n.dirs, lb =c(-Inf,0) )
# disc <- var.disc.1(l.var.rest, X, M=disc$M )
    # Take accout of lb on R.  NOPE
disc$p.lr <- (disc$M %^% 100)[1,]

m.sim.idx <- markov_sim( 1e6, disc$M, 0, nrow(disc$M) )
m.sim <- disc$X[m.sim.idx+1,]
l.var.m <- var.ols(t(m.sim))
    # Check the resulting VAR estimates

pdf('~/Dropbox/2017/research/debtLimits/charts/rg_uk_disc.pdf')
plot(disc$X, cex=disc$p.lr * 10 * sqrt( n.pts * n.dirs ), pch=16, ylab='Risk-free rate', xlab='Nominal growth rate')
    # Symmetry enforced
abline( 0, 1, lty=2 )
dev.off()

# pdf('~/Dropbox/2017/research/debtLimits/charts/rg_uk_disc.pdf')
interp <- disc.data.interp( disc, rg.dta[,c('gth','rfr')], 1, v.date = rg.dta$year,
                            mains=c('Nominal growth rate','Risk free rate',
                                    'Interest-growth differential'))
# dev.off()
rg.dta$s.idx <- as.factor(interp$s.idx)
# interp.maty <- disc.data.interp( disc.maty, rg.dta.maty[,c('gth','rfr')], 1, v.date = rg.dta.maty$year,
#                                  plot.on = FALSE )
# rg.dta.maty$s.idx <- as.factor(interp.maty$s.idx)

n.sim <- 1e6
n.plot <- 150
n.reg <- 1e6
m.sim <- markov_sim( n.sim, disc$M, which.max(disc$p.lr)-1, length(disc$p.lr))
sim <- disc$X[m.sim+1,]
colnames(sim) <- c('gth', 'rfr')
varnames <- c( 'Growth', 'Int. rate' )

mu.sim <- apply( sim, 2, mean )
l.var.ols.disc <- var.ols( t(sim), 1 )
# l.cty.mle.rest.1$est[[cty]]
var.table( list(data=l.var.rest, disc=l.var.ols.disc),
           file=paste0('~/Dropbox/2017/research/debtLimits/tables/uk_var_tab_disc.tex'),
           specnames = c( 'Data', 'Simulation' ), varnames = varnames,
           caption = 'Estimated restricted VAR based on data, and VAR estimated from
           discretized approximation simulated 100,000 periods' ,
           label=paste0('tab:uk_var_tab_disc'), footer=TRUE )

## 5.3 The primary balance function: linear correlation correct ##
rg.dta$era <- cut( rg.dta$year, c( -Inf, 1913, 1919, 1938, 1946, Inf ),
                   c('Pre-WW1', 'WW1', 'Inter-war', 'WW2', 'Post-WW2') )
d.poly.lm <- lm( pb ~ poly(d.lag,5,raw = TRUE) + rfr + gth,
                 data=rg.dta )
d.seq <- seq(0,max(rg.dta$d.lag,na.rm=TRUE)+40, by=1)
s.pred.0 <- predict( d.poly.lm, newdata = data.frame( d.lag=d.seq, rfr=0, gth=0) )
s.pred.mu <- predict( d.poly.lm, newdata = data.frame( d.lag=d.seq, rfr=mean(rg.dta$rfr),
                                                      gth=mean(rg.dta$gth)))
    # Predictions for s from a polynomial approximation
c.min.fn <- function( x ){
  coeff <- c( 0, 0, x )
  surp.fn.pred <- sapply( d.seq, surp, coeff=coeff, shift=0, tri=TRUE)
  # surp.fn.pred <- sapply( d.seq[d.seq<max(rg.dta$d.lag,na.rm=TRUE)],
  #                         surp, coeff=coeff, shift=0, tri=TRUE)
  return(sqrt(mean((surp.fn.pred-s.pred.0)^2)))
}
c.opt <- optim( c(150,-20,2), c.min.fn )
    # The best fit to the polynomial

paper.fn <- function( x, data ){
  return(x[1] + x[2] * ( data$d.lag - x[3] ) ^ 2 * (1*(data$d.lag<x[3]))
         + x[4]*data$rfr + x[5]*data$gth)
}
paper.min.fn <- function( x, data ){
  return(sqrt(mean(( paper.fn(x,data)- data$pb )^2, na.rm=TRUE )))
}
paper.coeff.init <- c(c.opt$par[3], (c.opt$par[2]-c.opt$par[3])/c.opt$par[1]^2,
                      c.opt$par[1], d.poly.lm$coefficients['rfr'], d.poly.lm$coefficients['gth'])
paper.min <- optim( paper.coeff.init, paper.min.fn, data=rg.dta )
paper.pred <- paper.fn( paper.min$par, data.frame( d.lag=d.seq, rfr=mean(rg.dta$rfr),
                                                   gth=mean(rg.dta$gth)) )


coeff.init <- c(0,0,c.opt$par)
shift.init <- d.poly.lm$coefficients['rfr'] * disc$X[,2] +
  d.poly.lm$coefficients['gth'] * disc$X[,1]
if(any(is.na(shift.init))) shift.init <- rep(0,nrow(disc$X))
    # No constant coefficient from the polynomial regression, as already in the
    # intercept of the nonlinear surplus function
plot( d.seq, s.pred.mu, type='l' )
lines( d.seq, s.pred.0, lty=2 )
abline(h=0)
lines( d.seq, sapply( d.seq, surp, coeff=coeff.init, shift=0, tri=TRUE), lty=2, lwd=2, col='blue' )
    # The polynomial and approxiamted function
lines(d.seq, paper.pred, col='red')
par(new=T)
hist(rg.dta$d.lag[-1], xlab='', ylab='', axes=F, main='',
     col=rgb(0,0,1,0.25), lty=0)
    # Plot the polynomial surplus function plus the debt density
rg.dta$surp <- sapply( 1:nrow(rg.dta),
                       function(i) surp( rg.dta$d.lag[i], coeff.init,
                                         shift.init[rg.dta$s.idx[i]], TRUE ) )
# surp.sd <- sd(rg.dta$surp - rg.dta$pb, na.rm = TRUE)
surp.sd <- with( subset(rg.dta, era=='Post-WW2'), sd(surp-pb, na.rm = TRUE))
with( rg.dta, plot( year, pb, type='l' ) )
with( rg.dta, lines( year, surp, type='l', col='red' ) )
with( rg.dta, lines( year, shift.init[s.idx], type='l', col='blue' ) )
abline(h=0)

## 5.4 Plotting the primary balance function ##
pars <- paper.min$par[1:3]
s.inc <- 2
coeff.init <- c(0,0,pars[3], pars[1]+pars[2]*pars[3]^2, pars[1] + s.inc)
shift.init <- paper.min$par[4] * disc$X[,2] + paper.min$par[5] * disc$X[,1]

params <- list( R=1+disc$X[,2]/100, G=1+disc$X[,1]/100, tri=TRUE, surp.sd=surp.sd,
                trans=disc$M, v.s.coeff=coeff.init, s.shift=shift.init)
plot.surp( params, x.lim = c(0,max(rg.dta$d.lag,na.rm = TRUE)), non.stoch = TRUE,
           vert=TRUE, leg = FALSE )
with(rg.dta, points( d.lag, pb, col=s.idx, pch=16) )
ns.b <- sol.nonstoch(params)


### 6. MODEL SOLUTION ###
params$surp.sd <- 1
params$lambda <- 0 ; params$phi <- .9 # .8 # .4
params$cont.type <- 'avg' ; params$diff.method <- "ana" ; params$inner.method <- 'all'
params$d.tri <- TRUE ; params$maxit <- 50
# # sol.s <- sol.search( params )
# sol <- sol.wrapper(params, plot.on=TRUE )
# plot.z(sol$p,sol$d,sol$params)

## 6.1 Solve model ##
# sol <- sol.small( params, 30, d.guess=TRUE )
sol <- sol.2(params, gain=.05, maxit=2000, d.all.init = 85 )

# stop()

# sol.b <- sol.2.way( params, 30 )
# params$surp.sd <- 1.1
# sol.sd <- sol.small( params, 30, init.guess = cbind(sol$p,sol$d), d.guess=TRUE )
params$surp.sd <- sol$params$surp.sd ; params$G <- sol$params$G - 5e-03 # A 0.5pp decrease in nominal growth
sol.G <- sol.2(params,init.guess = cbind(sol$p,sol$d), gain=.05, maxit=1500, print.level = 1 )
params$G <- sol$params$G ; params$s.shift <- sol$params$s.shift + .5 # A 1pp increase in surpluses
sol.shift <- sol.2(params,init.guess = cbind(sol$p,sol$d), gain=.05, maxit=1500, print.level = 1 )
params$s.shift <- sol$params$s.shift ; params$R <- sol$params$R + 5e-03 # A 1pp increase in interest rates
sol.R <-  sol.2(params,init.guess = cbind(sol.G$p,sol.G$d), gain=.05, maxit=1500, print.level = 1 )

sd.RG.seq <- seq( 0, .45, by=.05)
sol.RG <- list(sol)
for( i in 2:length(sd.RG.seq) ){
  params <- sol$params
  this.RG.red <- sd.RG.seq[i]
  params$R <- (1-this.RG.red) * sol$params$R + this.RG.red * sol$params$R[1]
  params$G <- (1-this.RG.red) * sol$params$G + this.RG.red * sol$params$G[1]
  gain <- if( i>10) .02 else .05
  sol.RG[[i]] <- sol.2(params,init.guess = cbind(sol.RG[[i-1]]$p,sol.RG[[i-1]]$d), gain=.05,
                       maxit=2000, print.level = 1 )
}

params <- sol$params
p.lr <- (sol$params$trans %^% 100)[1,]

pdf('~/Dropbox/2017/research/debtLimits/charts/ave_limit_RG_var.pdf')
plot( sd.RG.seq*100, sapply(sol.RG, function(x) sum( x$d * (x$params$trans %^% 100)[1,] ) ), type='l',  lwd=2,
      xlab='Percent reduction in interest and growth rate std dev', ylab='Average debt limit')
dev.off()


sol.maty <- list( sol )
maty.seq <- c(1, 2, 5, 6, 8, 9, 10 )
for( i in 2:length(maty.seq) ){
  this.maty <- maty.seq[i]
  maty.gps <- rep( 1:(ceiling(n.sim/this.maty)), each=this.maty )[1:n.sim]
  sim.maty <- ddply( as.data.frame(cbind(sim, maty.gps)), 'maty.gps', function(x) apply( x, 2, mean ) )
  l.var.maty <- var.ols(t(sim.maty[,1:2]))
  X.maty <- var.disc.pts(l.var.maty, 1, n.pts, n.dirs )
  disc.maty <- var.disc.1(l.var.maty, X.maty )
  params <- sol$params ; params$trans <- disc.maty$M
  params$R <- 1 + X.maty[,2] / 100 ;params$G <- 1 + X.maty[,1] / 100 ;
  params$surp.sd <- sol$params$surp.sd / sqrt(this.maty)
  gain <- if(this.maty>5) .01 else .02
  maxit <- if(i>=3) 20000 else 5000
  sol.maty[[i]] <- sol.2(params,init.guess = cbind(sol.maty[[i-1]]$p,sol.maty[[i-1]]$d), gain=gain, maxit=maxit,
                         print.level = 1 )
  # if( i == 4 ){
  #   sol.maty[[i]]$d[2:3] <- 160
  #   sol.maty[[i]] <- sol.2(params,init.guess = cbind(sol.maty[[i]]$p,sol.maty[[i]]$d), gain=gain, maxit=4000,
  #                          print.level = 1 )
  # }
}
params <- sol$params
ns.sol <- sol.nonstoch(params)
# stop()


pdf('~/Dropbox/2017/research/debtLimits/charts/ave_limit_maty.pdf')
plot( maty.seq, sapply(sol.maty, function(x) sum( x$d * (x$params$trans %^% 100)[1,] ) ), type='l',  lwd=2,
      xlab='Period length', ylab='Average debt limit')
dev.off()
plot( maty.seq, sapply(sol.maty, function(x) sqrt(sum( (x$d - sum(x$d * (x$params$trans %^% 100)[1,]))^2  *
                                                         (x$params$trans %^% 100)[1,] ) )), type='l', lwd=2,
      xlab='Calibration horizon', ylab='Std dev of debt limit')
plot( maty.seq, sapply( sol.maty, function(x) lm( x$d ~ x$params$R )$coeff )[2,]/100, type='l', lwd=2,
      xlab='Calibration horizon', ylab='Debt limit interest rate elasticity' )
plot( maty.seq, sapply( sol.maty, function(x) lm( x$d ~ x$params$G )$coeff )[2,]/100, type='l', lwd=2,
      xlab='Calibration horizon', ylab='Debt limit growth rate elasticity' )
pdf('~/Dropbox/2017/research/debtLimits/charts/rmg_elasty_maty.pdf')
plot( maty.seq, sapply( sol.maty, function(x){ rmg <- x$params$R - x$params$G; (lm( x$d ~ rmg )$coeff)[2]/100}), type='l', lwd=2,
      xlab='Calibration horizon', ylab='Debt limit interest-grwoth differential elasticity' )
dev.off()

plot(rg.dta$year, sol$d[rg.dta$s.idx], type='l')
n.smooth <- 5
pdf('~/Dropbox/2017/research/debtLimits/charts/ts_limits_uk.pdf')
plot(rg.dta$year, filter(sol$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), type='l', lwd=2,
     ylim=c(80,110), xlab='Year', ylab='Estimated debt limit' )
lines(rg.dta$year, filter(sol.G$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), col='blue', lwd=2, lty=2 )
points(rg.dta$year, filter(sol.G$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), col='blue', pch=1, cex=.75 )
lines(rg.dta$year, filter(sol.R$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), col='darkgreen', lwd=2, lty=4 )
points(rg.dta$year, filter(sol.R$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), col='darkgreen', pch=3, cex=.75 )
lines(rg.dta$year, filter(sol.RG[[2]]$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), col='red', lwd=2, lty=2 )
points(rg.dta$year, filter(sol.RG[[2]]$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), col='red', pch=2, cex=.75 )
# lines(rg.dta$year, filter(sol.shift$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), col='red', lwd=2 )
# lines(rg.dta$year, filter(sol.sd$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), col='yellow', lwd=2 )
legend('topright', c('Baseline', 'L/R growth 0.5pp lower', 'L/R interest rate 0.5pp higher',
                    'Growth and interest st dev. 5% lower'), #, 'Surpluses 0.5pp higher'),
       lwd=2, col=c('black','blue','darkgreen','red','yellow'), lty=c(1,2,4,1), pch=c(NA, 1, 3, 2 ), bty='n')
dev.off()
# plot(rg.dta$year, filter(sol.sd$d[rg.dta$s.idx], rep(1,n.smooth)/n.smooth, sides = 1 ), type='l' )

decade.bk <- seq(1879,2020, by=10)
rg.dta$decade <- cut(rg.dta$year, decade.bk,
                     labels = sapply(2:length(decade.bk), function(i) paste( decade.bk[i-(1:0)]+c(1,0), collapse='-' ) ) )
rg.dta$debt.limit <- sol$d[rg.dta$s.idx]
# rg.dta$debt.shift <- sol.shift$d[rg.dta$s.idx]
rg.dta$debt.R <- sol.R$d[rg.dta$s.idx]
rg.dta$debt.G <- sol.G$d[rg.dta$s.idx]
dec.dta <- ddply( rg.dta, 'decade', function(x) c( d.limit=mean(x$debt.limit), R.limit=mean(x$debt.R), G.limit=mean(x$debt.G),
                                                   R=mean(x$rfr), G=mean(x$gth), RmG=mean(x$rfr-x$gth ) ) )
with( dec.dta, plot( 1:length(d.limit), d.limit, type='l', xaxt='n', ylab='Debt limit', xlab='', lwd=2, ylim=c(80,100) ) )
with( dec.dta, lines( 1:length(R.limit), R.limit, col= 'blue', lwd=2, lty=2 ) )
with( dec.dta, lines( 1:length(G.limit), G.limit, col= 'darkgreen', lwd=2, lty=4 ) )
axis(1, 1:nrow(dec.dta), labels=dec.dta$decade )

with( dec.dta, plot( 1:length(d.limit), d.limit, type='l', xaxt='n', ylab='Debt limit', xlab='', lwd=2 ) )
with( dec.dta, plot( 1:length(G), G, type='l',lwd=2, xaxt='n', ylab='', xlab='', ylim=range(c(G,-G) ) ) )
axis(1, 1:nrow(dec.dta), labels=dec.dta$decade )
with( dec.dta, lines( 1:length(R), R, type='l', xaxt='n', ylab='', xlab='', col='blue', lwd=2  ) )
with( dec.dta, lines( 1:length(RmG), RmG, type='l', xaxt='n', ylab='', xlab='', col='red', lwd=2  ) )
legend('topleft', c('Growth', 'Interest rate', 'Differential'),
       lwd=2, col=c('black','blue','red'), bty='n')
abline(h=0)


paper.pred <- paper.fn( paper.min$par, data.frame( d.lag=d.seq, rfr=mean(rg.dta$rfr),
                                                   gth=mean(rg.dta$gth)) )
    # No constant coefficient from the polynomial regression, as already in the
    # intercept of the nonlinear surplus function
abc <- c(coeff.init[5], (coeff.init[4]-coeff.init[5]) / coeff.init[3]^2 , coeff.init[3])

pdf('~/Dropbox/2017/research/debtLimits/charts/surp_est.pdf')
plot( d.seq, s.pred.mu, type='l', lwd=2, xlab='Lagged debt/GDP ratio', ylab='Surplus' )
lines( d.seq, s.pred.0, lty=2 )
abline(h=0)
lines( d.seq, sapply( d.seq, surp, coeff=coeff.init, shift=0, tri=TRUE), lty=2, lwd=2, col='blue' )
    # The polynomial and approximated function
par(new=T)
hist(rg.dta$d.lag[-1], xlab='', ylab='', axes=F, main='',
     col=rgb(0,0,1,0.25), lty=0)
legend(x=170, y=22, c('Fifth-order polynomial', 'Bounded quadratic', 'Debt density'), lwd=2, lty=c(1,2, NA),
       fill=c(NA,NA,rgb(0,0,1,0.25)), col=c('black','blue', NA), border=NA, bty='n' )
# Plot the polynomial surplus function plus the debt density
dev.off()

### Create the IRFs ###
n.irf <- 20
wt.fn <- function( v, a, b, c ){
  # Solves for constants s,t such that v = t*a + s*b + (1-s-t)*c
  # Where a, b, c, and v are vectors
  wts <- solve( cbind( a-c, b-c), v-c )
  return( c( wts[1], wts[2], 1-sum(wts) ) )
}
lrv <- var_lr_variance( l.var.rest$A, l.var.rest$Sigma )
lr.cor <- lrv[1,2] / sqrt( prod( diag( lrv ) ) )
    # The long-run correlation
m.Sigma <- l.var.rest$Sigma
cor.innov <- ( m.Sigma[1,2] / sqrt( prod( diag( m.Sigma ) ) ) ) * sqrt( m.Sigma[1,1] )

shk <- list( R=c(0,sqrt(m.Sigma[2,2])), G=c(sqrt(m.Sigma[1,1]),0),
             both=c( cor.innov, sqrt(m.Sigma[2,2])) )
shk.dist <- lapply( shk, function(e) e / apply(disc$X,2,sum) )
    # Need to figure out how to do a shock to the distribution here

v.targ <- lapply( shk, function(x) disc$X[1,] + x )
# v.targ <- list( R=disc$X[9,]-c(0,.1), G=disc$X[9,]-c(.1,0), both=disc$X[9,] - .1 * c(lr.cor,1) )
    # Shocks to R, G and both
n.pts <- c(2,3,3)
    # Beause trying to solve for three points for R gives an error
l.cand.vert <- lapply( 1:3, function(i) c( 1, order(apply((t(disc$X[-1,]) - v.targ[[i]])^2,2,sum))[1:(n.pts[i]-1)]+1 ) )
# l.cand.vert <- list( R=c(9,2,3), G=c(9,4,5), both=c(9,4,5) )
# The candidate vertices to put weight on
# shk <- list( R=c(0,1), G=c(1,0), both=1*c(lr.cor,1) )

p.init <- list()
for( i in 1:3){
  this.p <- 0 * disc$p.lr
  if(n.obs[i]==2){
    p.vals <- solve( disc$X[l.cand.vert[[i]],] %*% t(disc$X[l.cand.vert[[i]],]), disc$X[l.cand.vert[[i]],] %*% v.targ[[i]], )
    this.p[l.cand.vert[[i]]] <- p.vals
  }else{
    opt <- optim( c(.5,.5),
                  function(p.v){
                    p.v <- c( 1-sum(p.v), p.v )
                    t.p <- 0 * disc$p.lr
                    t.p[ l.cand.vert[[i]] ] <- p.v
                    sum( ( t.p %*% disc$X - v.targ[[i]] ) ^ 2 )
                  } )
    this.p[l.cand.vert[[i]]] <- c( 1-sum(opt$par),opt$par )
  }
  p.init[[i]] <- this.p
}
names(p.init) <- names(v.targ)
shk.check <- sapply( p.init, function(p) p %*% disc$X - disc$X[1,] ) - do.call( cbind, shk )

l.irf <- list()
for( j in 1:3 ){
  irf.p <- matrix( 0, nrow=n.irf+1, ncol=nrow(disc$X) )
  irf.p[ 1, ] <- disc$p.lr
      # Initialize the state distribution
  irf.p[ 2, ] <- p.init[[j]]
      # Initialize the shock
  for( i in 2:n.irf){
    irf.p[i+1,] <- irf.p[i,] %*% disc$M
  }
      # Markov chain evolution
  l.irf[[j]] <- irf.p %*% cbind( disc$X, apply(disc$X,1,diff), sol$d, params$s.shift ) -
    rep( 1, n.irf+1 ) %*% ( irf.p[1,] %*% cbind( disc$X, apply(disc$X,1,diff), sol$d, params$s.shift ) )
      # The IRF is an average across states
}
v.col <- c('blue', 'red', 'darkgreen')
v.lty <- c( 1, 2, 3)
v.nm <- c('G','R','RmG','limit','surp')
v.l <- c( 'topright', 'topright', 'bottomright', 'bottomright', 'bottomright' )
for( j in 1:5 ){
  X <- 0:n.irf
  Y <- do.call( cbind, lapply(l.irf, function(x) x[,j]) )
  pdf(paste0('~/Dropbox/2017/research/debtLimits/charts/irf_uk_', v.nm[j], '.pdf' ) )
  plot( range(X), range(Y[,1:2]), type='n', xlab='Years', ylab='Response' )
  # for( i in 1:3 ) lines( X, Y[,i], lwd=2, col=v.col[i] )
  for( i in 1:2 ) lines( X, Y[,i], lwd=2, col=v.col[i], lty=v.lty[i] )
  abline(h=0,lwd=.5)
  # legend(v.l[j], c('R shock', 'G shock', 'Correlated shock'), lwd=2, col=v.col, bty='n')
  legend(v.l[j], c('R shock', 'G shock'), lwd=2, col=v.col, bty='n')
  dev.off()
}

elasty <- do.call( cbind, lapply(l.irf, function(x) x[,4]) )[2,1:2] / c(l.irf[[1]][2,2], l.irf[[2]][2,1])

##### Plot of plot.z.i here ####
pdf('~/Dropbox/2017/research/debtLimits/charts/sol_err_uk.pdf', width = 10, height = 15)
par(mfrow=c(16,4),mar=c(1,1,1,1),  mgp = c(2, .5, 0) )
An=c(0); Bn=c(0); Cn=c(0); def=matrix(0);
# for( i in 1:length(params$R)) plot.z.i( sol$p, sol$d, sol$params, i, An, Bn, Cn, def, 2 )
for( i in 1:length(params$R)) plot.z.i( sol$p, sol$d, sol$params, i, An, Bn, Cn, def, 2,
                         cex.main=.5, cex.lab = .5, cex.axis = .5,
                         xlim=c(0,max(sol$p) * 1.2), ylim=c(0,1)) #, lwd=1 )
dev.off()
par(mfrow=c(1,1))

#### Debt demand curves here? ####
