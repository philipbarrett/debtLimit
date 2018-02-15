## Code to calculate a calibrated solution with a small number of states

library(debtLimits)
library(lubridate)
library(xtable)

## Create the exogenous stocahstic process for R & G ##
cty <- 'USA'
rg.dta <- rg.read( cty )
# dta <- subset(dta, date <= "2015-01-01")
est <- var.rg.est(rg.dta)
disc <- var.discretize( est, n.pts = 1, n.dirs=4 )
# disc <- var.discretize( est, n.pts = 1, n.dirs=16 )
gkmoq.init <- nrow(disc$X)
    # An ok small discretization
var.table(est$VAR, file='~/Dropbox/2017/research/debtLimits/paper/usvar_tab.tex',
          add.mean = TRUE, varnames = c('Nom. GDP', 'Interest rate'),
          caption='VAR results for the US', label='tab:us_var')
    # The US VAR table

pdf('~/Dropbox/2017/research/debtLimits/charts/rg_us_disc.pdf')
plot( disc$X, cex=disc$m*15, pch=16, xlab='Quarterly nominal growth rate',
      ylab='Quarterly nominal interest rate' )
abline(0,1,lty=2)
legend('topright', 'R=G', bty='n', lty=2)
dev.off()
interp <- var.data.interp( disc, rg.dta )
rg.dta$x <- interp$s.idx
rg.dta$x.gth <- interp$dta.disc[,'gth']
rg.dta$x.rfr <- interp$dta.disc[,'rfr']
plot( disc$d.rmg )
n.states <- nrow(disc$X)

rmg.mu <- diff(est$mu)
rmg.sd <- se.rmg(est$VAR)
xx <- seq( rmg.mu - 4 * rmg.sd, rmg.mu + 4 * rmg.sd, length.out = 101 )
plot( xx, sapply(xx, function(x) 1 - pnorm( x, mean=rmg.mu, sd=rmg.sd ) ), type='l' )
print( paste0( "p-value of r-g >= 0: ", round( 1 - pnorm( 0, rmg.mu, rmg.sd ), 3 ) ) )


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

m.surp <- data.frame( 'Nom int rate'=(params$R-1)*100, 'Nom gth rate'=100*(params$G-1),
                      "Interest-growth differential"=(params$R-params$G)*100,
                      'Constant'=shift.opt$par - min(shift.opt$par) )
print.xtable(
  xtable( m.surp, label = 'tab:surp_shift', align = rep('c',ncol(m.surp)+1),
          caption =
            'State-dependent constants in surplus function (relative to lowest state)'),
  sanitize.colnames.function=function(x)gsub("\\."," ",x ),
  include.rownames = FALSE, file = '~/Dropbox/2017/research/debtLimits/paper/surp_shift.tex' )

## Now try solving the model ##
params$surp.sd <- .05 # Starting small
params$lambda <- 0
params$phi <- .8 # .4
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

# for( i in 1:length(params$R) ){
#   xx <- sol.search.i( params, cbind(1e-05,d.init), i, An, Bn, Cn, def,
#                     print.level = 1 )
#   plot.z( xx$p, xx$d, params, An, Bn, Cn, def, xlim=c(0,2*xx$p[i]), ylim=c(0,2*xx$p[i]) )
# }

plot.z( rep(0,length(params$R)), d.init, params ) #, xlim=c(0,3e-03), ylim=c(0,3e-03) )
plot.z.d( rep(0,length(params$R)), d.init, params ) #, xlim=c(0,3e-03), ylim=c(0,3e-03) )

params$it <- 30
params$tol <- 0.25
sol.s <- sol.search(params) #, plot.on = TRUE)

params$surp.sd <- .075 # Try increasing the variance
sol.t <- sol.search(params, cbind( sol.s$p, sol.s$d )) #, plot.on = TRUE)
# params$surp.sd <- 1 # No bueno.  Sends d~=0
sol.u <- sol.search(params, cbind( sol.t$p, sol.t$d )) #, plot.on = TRUE)


params$tol <- 1e-6
sol <- sol.wrapper( params, init.guess = cbind( sol.t$p, sol.t$d ) ) #, plot.on=T )
plot.surp(params, x.lim = c(0, 1.5 * max(sol$d) ) )
abline(v=sol$d,lty=2)
points( dta[,c('cnlb_gdp_lag','pb_gdp')], col=dta$x )
plot.z( sol$p, sol$d, sol$params , xlim=c(0,2e-04), ylim=c(0,2e-04) )
plot.z( sol$p, sol$d, sol$params )

params.d <- params
dec <- 1e-4
params.d$R <- params$R - dec
sol.d <- sol.wrapper( params.d, init.guess = cbind( sol$p, sol$d ) )
plot.z( sol.d$p, sol.d$d, sol.d$params , xlim=c(0,2e-04), ylim=c(0,2e-04) )
dd.dR <- (sol.d$d -sol$d) / dec / 100 / 4
    # Divide by 100 to get the response to a 1pp change in the L/R annual interest rate

dl.ts <- data.frame(date=dta$date, limit=sol$d[dta$x] )
annual.ratio.decade <- by( dl.ts$limit, 10 * year(dl.ts$date) %/% 10, mean ) / 4

params.no.RG.var <- params
params.no.RG.var$R <- rep( sum(params$R * disc$m), length(params$R) )
params.no.RG.var$G <- rep( sum(params$G * disc$m), length(params$G) )
sol.no.RG.var.init <- sol.search( params.no.RG.var ) #, plot.on=T )
sol.no.RG.var <- sol.wrapper( params.no.RG.var, # plot.on=T,
                              init.guess = cbind( sol.no.RG.var.init$p, sol.no.RG.var.init$d ) )

plot.z( sol.no.RG.var$p, sol.no.RG.var$d, sol.no.RG.var$params, An, Bn, Cn, def,
        xlim=c( 0, 2*max(sol.no.RG.var$p)), ylim=c( 0, 2*max(sol.no.RG.var$p)) )
plot.z( sol.no.RG.var$p, sol.no.RG.var$d, sol.no.RG.var$params, An, Bn, Cn, def )

pdf('~/Dropbox/2017/research/debtLimits/charts/limits_us.pdf')
plot(dta$date, dta$cnlb_gdp / 4, type='l', lwd=2, ylim=c(0,max(sol.no.RG.var$d/4)),
     xlab='', ylab='Debt/Annual GDP ratio' )
lines(dta$date, sol$d[dta$x] / 4, col='blue', lwd=2 )
lines(dta$date, sol.no.RG.var$d[dta$x] / 4, col='red', lwd=2 )
legend('bottomright',
       c('Data (excludes SSA)', 'Estimated debt limit',
         'Estimated debt limit for fixed R, G' ), bty='n',
       col=c('black','blue','red'), lwd=2)
dev.off()

### Compute a discretized simulation
n.sim <- 1e6
n.plot <- 150
n.reg <- 1e6
m.sim <- markov_sim( n.sim, disc$trans, which.max(disc$m)-1, length(disc$m))
sim <- disc$X[m.sim+1,]
colnames(sim) <- c('gth', 'rfr')

mu.sim <- apply( sim, 2, mean )
sim.VAR <- VAR(sim[1:n.reg,])
var.table(sim.VAR, file='~/Dropbox/2017/research/debtLimits/paper/usvar_tab_disc.tex',
          add.mean = TRUE, varnames = c('Nom. GDP', 'Interest rate'),
          caption=paste0('VAR results for the discretized simulation of the US. Results from ',
                         n.reg, ' periods'),
          label='tab:us_var_disc', add.se = FALSE )

