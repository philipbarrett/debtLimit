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
cty <- 'USA'
    # Country
dta <- read.csv( paste0( 'data/', cty, '.csv' ) )
dta$Date <- as.Date( dta$Date, "%m/%d/%Y" )
    # Read in data
dta$pub_D.yrend <- mapply( function(x,y) if(month(x)==10) y else NA, dta$Date, dta$Pub_D )
dta$pub_D.yrstart <- mapply( function(x,y) if(month(x)==1) y else NA, dta$Date, dta$Pub_D )
    # Add year end/start debt (only one is true)
dta$OB_GDP4 <- dta$OB_GDP / 4
    # Compute the overall balance as a share of annual GDP
dta$Pub_D.imp <- dta$Pub_D.imp.1 <- dta$Pub_D
dta$Pub_D.imp[-1] <- dta$Pub_D.imp.1[-1] <- NA
for( i in 2:nrow(dta) ){
  dta$Pub_D.imp[i] <- ( dta$Pub_D.imp[i-1]  ) / ( 1 + dta$ngdp_pch[i] /100 ) - dta$OB_GDP4[i] # + dta$Int_Exp_2[i]
  if(month(dta$Date[i])==10){
    dta$Pub_D.imp.1[i] <- dta$Pub_D[i]
  }else{
    dta$Pub_D.imp.1[i] <- ( dta$Pub_D.imp.1[i-1]  ) / ( 1 + dta$ngdp_pch[i] /100 ) - dta$OB_GDP4[i]
  }
}
    # Compute implied debt (assumes debt is year-end) face value
dta$ngdp_pch4 <- filter( dta$ngdp_pch, c(1,1,1,1), sides=1 ) / 4

## 2. Estimate the unconstrained surplus model ##
scaling <- 60
x.max <- 200
dta$Pub_D.lag <- c( NA, dta$Pub_D.imp.1[-1] )
n.poly <- 3
mod <- lm( OB_GDP4 ~ poly(Pub_D.lag, n.poly, raw = TRUE) + ngdp_pch,
           data=subset(dta, Date > "1940-01-01"), na.action = na.exclude )
pred <- predict(mod, dta)
v.s.coeff <- c( scaling, mod$coefficients[1:(n.poly+1)] * scaling ^ (0:(n.poly)), mod$coefficients['ngdp_pch'] )
kk <- kmeans( cbind(1+dta$rfr/400,1+dta$ngdp_pch/100), 3)
plot( 1+dta$rfr/400, 1+dta$ngdp_pch/100, xlab='rfr', ylab='ngdp_pch')
points(ff$centers, col='red', pch=8 )
params <- list(v.s.coeff=v.s.coeff, R=kk$centers[,1], G=kk$centers[,2], tri=FALSE)
pred.scale <- mapply( surp, dta$Pub_D.lag, 1+dta$ngdp_pch/100, MoreArgs = list(coeff=v.s.coeff, tri=FALSE) )
plot.surp(params, xlim=c(0,x.max))
with( subset( dta, ngdp_pch > mean(ngdp_pch) ), points( Pub_D.lag, OB_GDP4, col='blue' ) )
with( subset( dta, ngdp_pch < mean(ngdp_pch) ), points( Pub_D.lag, OB_GDP4, col='red' ) )
par(new = T)
hist(dta$Pub_D, axes=F, xlab=NA, ylab=NA, xlim=c(0,x.max), main='', col=rgb(0,0,1,alpha=0.2),
     breaks = 20  )

## 3. Estimate the alternative surplus model ##
gkmoq.init <- c( 1, 4, -.2080, .0032, -.0000122, 0)
    # The GKMOQ estimates
params.gkmoq <- params
params.gkmoq$v.s.coeff <- c( scaling, .25 * gkmoq.init[2:5] * scaling ^ (0:3), 0 )
    # Introduce some scaling and make quarterly
plot.surp(params.gkmoq)
    # plot the function
gkmo.opt <- optim( params.gkmoq$v.s.coeff[c(2,6)],
             function(par) sum( ( mapply( surp, dta$Pub_D.lag, 1+dta$ngdp_pch/100,
                                          MoreArgs = list(coeff=c(params.gkmoq$v.s.coeff[1],
                                                                  par[1], params.gkmoq$v.s.coeff[3:5],
                                                                  par[2]),
                                                                  tri=params$tri) ) -
                                    dta$OB_GDP4 ) ^ 2, na.rm = TRUE ) )
    # Substitute in the Ghosht-Kim-Mendoza-Ostry params
params.alt <- params.gkmoq
params.alt$v.s.coeff[c(2,6)] <- gkmo.opt$par
pred.alt <- mapply( surp, dta$Pub_D.lag, 1+dta$ngdp_pch/100,
                      MoreArgs = list(coeff=params.alt$v.s.coeff, tri=params.alt$tri) )
    # Try minimizing absolute error
plot.surp(params.alt)
with( subset( dta, ngdp_pch > mean(ngdp_pch) ), points( Pub_D.lag, OB_GDP4, col='blue' ) )
with( subset( dta, ngdp_pch < mean(ngdp_pch) ), points( Pub_D.lag, OB_GDP4, col='red' ) )


## 4. Make some general plots ##
mar.dft <- par('mar')
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(dta$Date[!is.na(dta$pub_D.yrend)], na.omit(dta$pub_D.yrend), type='l', xlab='', ylab='', lwd=2,
     main='Debt-GDP ratio', ylim=c(0,max(dta$pub_D.yrend,na.rm = TRUE)))
lines(dta$Date, dta$Pub_D.imp, lwd=2, col='blue' )
lines(dta$Date, dta$Pub_D.imp.1, lwd=2, col='red' )
plot(dta$Date, dta$ngdp_pch4, type='l', xlab='', ylab='',
     main='Four-quarter nominal GDP growth')
plot(dta$Date, dta$OB_GDP4, type='l', xlab='', ylab='',
     main='Quarterly surplus-Annual GDP ratio')
lines( dta$Date, pred, col='red' )
lines( dta$Date, pred.scale, col='blue' )
lines( dta$Date, pred + sd(mod$residuals), col='red', lty=2 )
lines( dta$Date, pred - sd(mod$residuals), col='red', lty=2 )
lines( dta$Date, pred.alt, col='green' )
plot(dta$Date, dta$rfr, type='l', xlab='', ylab='', main='Risk free rate')
lines( dta$Date, dta$rfr + dta$cds_5y / 100, col='blue', lty=2 )
par(mfrow=c(1,1), mar=mar.dft)

## 5. Create a transition matrix ##
freq <- table( kk$cluster[-(nrow(dta))], kk$cluster[-1] )
params.alt$trans <- freq / apply( freq, 1, sum )
    # Brief transition probabilities
params.alt$surp.sd <- sd( dta$OB_GDP4 - pred.alt, na.rm = TRUE )
save(params.alt, file='data/surpEst.rdta')
