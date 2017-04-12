
## 1. Load parameters ##
rm(list=ls())
library(debtLimits)
library(expm)
params <- list()
load('data/surpEst.rdta')
params <- params.alt
params$tri <- FALSE # TRUE
params$lambda <- .95 # 0 # .5
    # This actually works better with high lambda because the surplus risk is less problematic
params$phi <- .6
params$cont.type <- 'avg'
params$q.e <- c(0)
params$def <- matrix(0,1,1)
params$diff.method <- "ana"
params$d.tri <- FALSE # TRUE # FALSE      # Triangular distribution for surplus shocks
params$inner.method <- 'all'
An <- 1 / params$R
Bn <- rep( -1, length(params$R) )
Cn <- An
def <- matrix(0,1,1)
nn <- length(params$R)

## 2. Mess with the measured prameters ##
plot.surp(params)
params$surp.sd <- params$surp.sd * .1
    # Reduce the volatility of the shocks
params$trans <- matrix( c(.99,  .005, .005,
                          .005, .99,  .005,
                          .005, .005, .99 ), nrow = 3, ncol=3, byrow = TRUE )
    # To suppress the problems with state 2 in particular. There, the expost
    # risk to R-G is all downside.
# params$G <- rep( mean(params$G), nn )
# params$R <- rep( mean(params$R), nn )
params$d.tri <- FALSE
plot.surp(params)

## 3. Solve the model ##
sol <- sol.wrapper(params )
plot.z( sol$p, sol$d, sol$params, xlim=c(0,min(1,2*max(sol$p))), ylim=c(0,min(1,2*max(sol$p))) )

## 4. Try solving the outer solution ##
