#####################################################################
# surpEst.R
#
# Estimates the surplus function
# 07apr2017
# Philip Barrett, Washington DC
#####################################################################

library(debtLimits)
library(lubridate)
library(zoo)

## 1. Read in the data ##
rm(list=ls())
cty <- 'USA'
    # Country
dta <- read.csv( paste0( 'data/', cty, '.csv' ) )
dta$Date <- as.Date( dta$Date, "%m/%d/%Y" )
    # Read in data
# dta$pub_D.yrend <- mapply( function(x,y) if(month(x)==10) y else NA, dta$Date, dta$Pub_D )
# dta$pub_D.yrstart <- mapply( function(x,y) if(month(x)==1) y else NA, dta$Date, dta$Pub_D )
#     # Add year end/start debt (only one is true)
# dta$OB_GDP4 <- dta$OB_GDP / 4
#     # Compute the overall balance as a share of annual GDP
# dta$Pub_D.imp <- dta$Pub_D.imp.1 <- dta$Pub_D
# dta$Pub_D.imp[-1] <- dta$Pub_D.imp.1[-1] <- NA
# for( i in 2:nrow(dta) ){
#   dta$Pub_D.imp[i] <- ( dta$Pub_D.imp[i-1]  ) / ( 1 + dta$ngdp_pch[i] /100 ) - dta$OB_GDP4[i] # + dta$Int_Exp_2[i]
#   if(month(dta$Date[i])==10){
#     dta$Pub_D.imp.1[i] <- dta$Pub_D[i]
#   }else{
#     dta$Pub_D.imp.1[i] <- ( dta$Pub_D.imp.1[i-1]  ) / ( 1 + dta$ngdp_pch[i] /100 ) - dta$OB_GDP4[i]
#   }
# }
    # Compute implied debt (assumes debt is year-end) face value
dta$ngdp_pch4 <- filter( dta$ngdp_pch, c(1,1,1,1), sides=1 ) / 4

par(mfrow=c(2,1))
with( dta, plot(Date, pb_gdp, type='l' ) )
abline(h=0, col='blue')
with( dta, plot(Date, cnlb_gdp, type='l' ) )
par(mfrow=c(1,1))
with( dta, plot(cnlb_gdp_lag, pb_gdp, col=1+1*( ngdp_pch>mean(ngdp_pch) ), pch=16 ) )
with( dta, lines(cnlb_gdp_lag, pb_gdp, lwd=.5 ) )


## 2. Estimate the unconstrained surplus model ##
scaling <- 80
x.max <- 400
n.poly <- 3
n.states <- 2
set.seed(42)
kk <- kmeans( cbind(1+dta$rfr/400,1+dta$ngdp_pch/100), n.states)
plot( 1+dta$rfr/400, 1+dta$ngdp_pch/100, col=kk$cluster,
      ylab='Growth', xlab='Risk free rate')
points(kk$centers, col=1:n.states, pch=15, cex=2 )
abline(0,1,lty=2)
    # Discretize the states
dta$x <- as.factor(kk$cluster)
mod <- lm( pb_gdp ~ poly(cnlb_gdp_lag, n.poly, raw = TRUE), # + poly(ngdp_pch, 2, raw=TRUE),
           data=subset(dta), na.action = na.exclude )
pred <- predict(mod, dta)
v.s.coeff <- c( scaling, mod$coefficients[1:(n.poly+1)] * scaling ^ (0:(n.poly)) )
shift <- rep(0, n.states) # c(0, tail(mod$coefficients, n.states-1))
params <- list(v.s.coeff=v.s.coeff, R=rep(mean(dta$rfr)/4,n.states),
               G=rep(mean(dta$ngdp_pch),n.states), s.shift=shift, tri=FALSE)
pred.scale <- mapply( surp, dta$cnlb_gdp_lag, params$s.shift[dta$x], MoreArgs = list(coeff=v.s.coeff, tri=FALSE) )
plot.surp(params, x.lim = c(0,x.max), ylim=c(-20,4))
with( dta, points( cnlb_gdp_lag, pb_gdp ) ) #, col=x ) )
# par(new = T)
# hist(dta$Pub_D, axes=F, xlab=NA, ylab=NA, xlim=c(0,x.max), main='', col=rgb(0,0,1,alpha=0.2),
#      breaks = 20  )

## 3. Estimate the alternative surplus model ##
gkmoq.init <- c( 1, 4, -.2080, .0032, -.0000122)
    # The GKMOQ estimates
params.gkmoq <- params
params.gkmoq$v.s.coeff <- c( scaling, .25 * gkmoq.init[-1] * scaling ^ (0:3) )
    # Introduce some scaling and make quarterly
params.gkmoq$s.shift <- rep(params.gkmoq$v.s.coeff[2],n.states)
params.gkmoq$v.s.coeff[2] <- 0
    # Make the intercept zero - put all the difference in the shift
plot.surp(params.gkmoq,x.lim = c(0,250), ylim=c(-3,3) )
    # plot the function
gkmo.opt <- optim( params.gkmoq$s.shift,
             function(par) sum( ( mapply( surp, dta$cnlb_gdp_lag, par[dta$x],
                                          MoreArgs = list(coeff=params.gkmoq$v.s.coeff,
                                                                  tri=params$tri) ) -
                                    dta$OB_GDP4 ) ^ 2, na.rm = TRUE ) )
    # Substitute in the Ghosh-Kim-Mendoza-Ostry-Qureshi params
params.alt <- params.gkmoq
params.alt$s.shift <- gkmo.opt$par
pred.alt <- mapply( surp, dta$cnlb_gdp_lag, params.alt$s.shift[dta$x],
                    MoreArgs = list(coeff=params.gkmoq$v.s.coeff, tri=params$tri) )
    # Try minimizing absolute error
plot.surp(params.alt, ylim=c(-4,4), x.lim=c(0,x.max))
with( dta, points( cnlb_gdp_lag, pb_gdp ), col=x )
# with( dta, points( Pub_D.lag, OB_GDP4, col=x ) )
# with( subset( dta, ngdp_pch > mean(ngdp_pch) ), points( Pub_D.lag, OB_GDP4, col='blue' ) )
# with( subset( dta, ngdp_pch < mean(ngdp_pch) ), points( Pub_D.lag, OB_GDP4, col='red' ) )


## 4. Make some general plots ##
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(dta$Date[!is.na(dta$pub_D.yrend)], na.omit(dta$pub_D.yrend), type='l', xlab='', ylab='', lwd=2,
     main='Debt-GDP ratio', ylim=c(0,max(dta$pub_D.yrend,na.rm = TRUE)))
lines(dta$Date, dta$Pub_D.imp, lwd=2, col='blue' )
lines(dta$Date, dta$Pub_D.imp.1, lwd=2, col='red' )
plot(dta$Date, dta$ngdp_pch4, type='l', xlab='', ylab='',
     main='Four-quarter nominal GDP growth')
plot(dta$Date, dta$OB_GDP4, type='l', xlab='', ylab='',
     main='Quarterly surplus-Annual GDP ratio')
lines( dta$Date, pred.scale, col='red' )
lines( dta$Date, pred, col='blue' )
lines( dta$Date, pred + sd(mod$residuals), col='red', lty=2 )
lines( dta$Date, pred - sd(mod$residuals), col='red', lty=2 )
lines( dta$Date, pred.alt, col='green' )
plot(dta$Date, dta$rfr, type='l', xlab='', ylab='', main='Risk free rate')
lines( dta$Date, dta$rfr + dta$cds_5y / 100, col='blue', lty=2 )
par(mfrow=c(1,1), mar=rep(3,4))

## 5. Create a transition matrix ##
freq <- table( kk$cluster[-(nrow(dta))], kk$cluster[-1] )
params.alt$trans <- freq / apply( freq, 1, sum )
    # Brief transition probabilities
params.alt$surp.sd <- sd( dta$OB_GDP4 - pred.alt, na.rm = TRUE )
params.alt$surp.sd <- sd(mod$residuals)
    # Take the variance from the well-fitted model
save(params.alt, file='data/surpEst.rdta')
