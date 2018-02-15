# Makes tables for the Mueller-Watson approach
library(xtable)

cts.a <- c('USA', 'UK', 'FRA', 'DEU' )
cts.q <- c('USA', 'UK', 'FRA', 'CAN' )
cts <- union( cts.a, cts.q )
n.cts <- length(cts)
cts.name <- c( 'USA'='United States', 'UK'='United Kingdom', 'FRA'='France',
               'DEU'='Germany', 'JPN'='Japan', 'ITA'='Italy', 'CAN'='Canada')

dir.a <- '/home/philip/Dropbox/2017/research/debtLimits/spectral/export/annual/'
dir.q <- '/home/philip/Dropbox/2017/research/debtLimits/spectral/export/quarterly/'

l.ci.a <- list() ; l.ci.q <- list() ; l.d.a <- list() ; l.d.q <- list()

for( cty in cts ){
  if(cty %in% cts.a){
    dta.in.a <- read.csv(paste0(dir.a, 'ci_lr_rmg_', tolower(cty), '_90.csv'), header = TRUE)
        # Read data in
    l.ci.a[[cty]] <- matrix( dta.in.a[ c(5:6,1:4), ncol(dta.in.a)], ncol=2, nrow=3, byrow = TRUE )
    l.d.a[[cty]] <- read.csv(paste0(dir.a, 'll_rmg_', tolower(cty), '.csv'), header = FALSE)
  }
  if(cty %in% cts.q){
    dta.in.q <- read.csv(paste0(dir.q, 'ci_lr_rmg_', tolower(cty), '_90.csv'), header = TRUE)
    l.ci.q[[cty]] <- matrix( dta.in.q[ c(5:6,1:4), ncol(dta.in.q)], ncol=2, nrow=3, byrow = TRUE )
    l.d.q[[cty]] <- read.csv(paste0(dir.q, 'll_rmg_', tolower(cty), '.csv'), header = FALSE)
  }
}

l.mid.up.a <- lapply( l.ci.a, function(x) cbind( apply( x, 1, mean ), x[,2] ) )
l.mid.up.q <- lapply( l.ci.q, function(x) cbind( apply( x, 1, mean ), x[,2] ) )
    # The midpoints and upper boundaries

m.d.a <- cbind( d=l.d.a[[1]][,1], do.call(cbind, lapply(l.d.a, function(x) x[,2]) ) )
m.d.q <- cbind( d=l.d.q[[1]][,1], do.call(cbind, lapply(l.d.q, function(x) x[,2]) ) )
d.seq <- seq( -1, 1, by=.001 )
d.opt.a <- sapply( cts.a,
                   function(x) d.seq[which.max(predict(lm( paste0( x, ' ~ poly(d,2)'),
                                 as.data.frame(m.d.a) ), data.frame(d=d.seq) ))] )
d.opt.q <- sapply( cts.q,
                   function(x) d.seq[which.max(predict(lm( paste0( x, ' ~ poly(d,2)'),
                               as.data.frame(m.d.q) ), data.frame(d=d.seq) ))] )

df.a <- data.frame( rbind( do.call( cbind, l.mid.up.a ), rep( d.opt.a, each=2 ) ) )
rownames(df.a) <- c('I(0)', 'Bayes', 'Bayes superset', 'ML fractional integration')
print( xtable( round(df.a,1) ),
       file='/home/philip/Dropbox/2017/research/debtLimits/tables/spectral_a.tex' )
df.q <- data.frame( rbind( do.call( cbind, l.mid.up.q ), rep( d.opt.q, each=2 ) ) )
rownames(df.q) <- c('I(0)', 'Bayes', 'Bayes superset', 'ML fractional integration')
print( xtable( round(df.q,1) ),
       file='/home/philip/Dropbox/2017/research/debtLimits/tables/spectral_q.tex' )

i.cts.a <- which( cts %in% cts.a )
i.cts.q <- which( cts %in% cts.q )

pdf('~/Dropbox/2017/research/debtLimits/charts/ci_compare.pdf', width = 10, height=8)

plot( i.cts.a, sapply(l.mid.up.a, function(x)x[1,1]), pch=16, col='black', ylim=c(-4,4), xlim=c(.5,n.cts+.5),
      xaxt = "n", xlab='', ylab='Long-run interest-growth differential', cex=1.5 )
axis(1, at=1:n.cts, labels=cts.name[cts])
abline(h=0)
points( i.cts.a, sapply(l.mid.up.a, function(x)x[1,2]), pch=2, cex=1.5, col='black', lwd=2)
# points( i.cts.a, sapply(l.mid.up.a, function(x)x[2,1]), pch=16, col='darkgreen')
# points( i.cts.a, sapply(l.mid.up.a, function(x)x[2,2]), pch='+', col='darkgreen')
load('data/rmg_est.rdata')
l.thresh.4 <- list( CAN=c(.2525, .58)* 4,
                    DEU=c(.26,.33)* 4,
                    FRA=c(.0,.175)* 4,
                    UK=c(.06, .169)* 4,
                    USA=c(.015, .19)* 4,
                    JPN=c(.0275, .1775)* 4,
                    ITA=c(.0275, .1775)* 4 )
l.cty.mle.4$est$UK <- l.cty.mle.4$est$GBR
points( i.cts.q, sapply( l.cty.mle.4$est[cts.q], function(x) diff(x$mu) ), pch=16, col='blue', cex=1.5)
points( i.cts.q, sapply( l.thresh.4[cts.q], function(x) x[1] ), pch=2, lwd=2, cex=1.5, col='blue')
load('data/rmg_hist__1880_est.rdata')
l.thresh.4.hist <- list( CAN=c(2.15, NA),
                         DEU=c( .37, .79),
                         FRA=c( .26, .81),
                         UK=c(-.03, .58),
                         USA=c(-.065,.38),
                         JPN=c(3.75, NA),
                         ITA=c(2.67, 3.87) )
l.cty.mle.4$est$UK <- l.cty.mle.4$est$GBR
points( i.cts.a, sapply( l.cty.mle.4$est[cts.a], function(x) diff(x$mu) ), pch=16, col='red', cex=1.5)
points( i.cts.a, sapply( l.thresh.4.hist[cts.a], function(x) x[1] ), pch=2, lwd=2, cex=1.5, col='red')
legend('topleft', c('VAR point estimate', 'VAR 95% critical value',
                    'MW I(0) midpoint', 'MW I(0) 95% critical value',
                    # 'MW Bayes midpoint', 'MW Bayes 95% critical value',
                    'VAR point estimate: quarterly data 1956-2016', 'VAR 95% critical value: quarterly data 1956-2016'),
       pch=rep(c(16,2), 4), col=rep(c('red','black','blue'), each=2), bty='n', cex=1.2, lwd=2, lty=NA)
dev.off()
