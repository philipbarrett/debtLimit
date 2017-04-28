## Code to calculate a calibrated solution with a small number of states

library(debtLimits)

## Create the exogenous stocahstic process for R & G ##
cty <- 'USA'
rg.dta <- rg.read( cty )
# dta <- subset(dta, date <= "2015-01-01")
est <- var.rg.est(rg.dta)
disc <- var.discretize( est, n.pts = 1, n.dirs=4 )
gkmoq.init <- nrow(disc$X)
    # An ok small discretization
plot( disc$X, cex=disc$m*100 )
interp <- var.data.interp( disc, rg.dta )
rg.dta$x <- interp$s.idx
rg.dta$x.gth <- interp$dta.disc[,'gth']
rg.dta$x.rfr <- interp$dta.disc[,'rfr']
plot( disc$d.rmg )
n.states <- nrow(disc$X)

## Create the basic surplus function ##
dta.all <- read.csv( paste0( 'data/', cty, '.csv' ) )
dta.all$date <- as.Date( dta.all$date, "%m/%d/%Y" )
dta <- merge( dta.all, rg.dta )
for( nm in c('x','x.gth','x.rfr') ) dta[[nm]] <- as.factor(dta[[nm]])
n.poly <- 2
mod <- lm( pb_gdp ~ poly(cnlb_gdp_lag,n.poly) + gth, data = dta)

## Convert to form used in model ##
x.max <- 800
gkmoq.init <- c( 4, 0, -.2080, .0032, -.0000122)
    # First entry is scaling
v.s.coeff <- gkmoq.init
shift.opt <- optim( rep(0,n.states),
              function(par) sum( ( mapply( surp, dta$cnlb_gdp_lag, par[dta$x],
                                           MoreArgs = list(coeff=v.s.coeff,
                                                           tri=FALSE) ) -
                                     dta$pb_gdp ) ^ 2, na.rm = TRUE ) )
    # Compute the optimal level shifts given the polynomial coefficients
params <- list(v.s.coeff=v.s.coeff, R=1+disc$X[,'rfr']/100, G=1+disc$X[,'gth']/100,
               s.shift=shift.opt$par, tri=FALSE, trans=disc$trans )
plot.surp( params, TRUE, c(0,x.max) )
points( dta[,c('cnlb_gdp_lag','pb_gdp')], col=dta$x )
    # TO ADD: function to plot the fitted surplus function

## Now try solving the model ##
params$surp.sd <- .05 # Starting small
params$lambda <- 0
params$phi <- .4
params$cont.type <- 'avg'
params$diff.method <- "ana"
params$inner.method <- 'all'
params$d.tri <- FALSE # TRUE # FALSE      # Triangular distribution for surplus shocks
An <- 1 / params$R
Bn <- rep( -1, length(params$R) )
Cn <- An
def <- matrix(0,1,1)

# d.init <- rep( max( sol.nonstoch(params) ), length( params$R ) )
d.init <- sol.nonstoch(params)

for( i in 1:length(params$R) ){
  xx <- sol.search.i( params, cbind(1e-05,d.init), i, An, Bn, Cn, def,
                    print.level = 1 )
  plot.z( xx$p, xx$d, params, An, Bn, Cn, def, xlim=c(0,2*xx$p[i]), ylim=c(0,2*xx$p[i]) )
}



plot.z( rep(0,length(params$R)), d.init, params ) #, xlim=c(0,3e-03), ylim=c(0,3e-03) )
plot.z.d( rep(0,length(params$R)), d.init, params ) #, xlim=c(0,3e-03), ylim=c(0,3e-03) )

params$it <- 30
params$tol <- 0.25
sol.s <- sol.search(params, plot.on = TRUE)

params$surp.sd <- .075 # Try increasing the variance
sol.t <- sol.search(params, cbind( sol.s$p, sol.s$d ), plot.on = TRUE)
params$surp.sd <- .1 # This is about the highest this will go
sol.u <- sol.search(params, cbind( sol.t$p, sol.t$d ), plot.on = TRUE)


params$tol <- 1e-6
sol <- sol.wrapper( params, plot.on=T, init.guess = cbind( sol.u$p, sol.u$d ) )
plot.surp(params, x.lim = c(0, 1.5 * max(sol$d) ) )
abline(v=sol$d,lty=2)
points( dta[,c('cnlb_gdp_lag','pb_gdp')], col=dta$x )
plot.z( sol$p, sol$d, sol$params , xlim=c(0,1e-04), ylim=c(0,1e-04) )
plot.z( sol$p, sol$d, sol$params )

