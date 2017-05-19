## Code to analyze the restricted VARs for the various countries ##

v.cty <- c('USA', 'GBR', 'CAN', 'FRA', 'DEU')
n.cty <- length(v.cty)
theta.range <- c(0,2)
l.p <- lapply( v.cty, function(x) var.rg.rest.theta( x, theta.range = theta.range,
                                                     start="1970-01-01" ) )
v.theta <- l.p[[1]][,'v.theta']
m.p <- sapply( l.p, function(x) x[,2] )
colnames(m.p) <- v.cty
    # Extract the values

pdf('~/Dropbox/2017/research/debtLimits/charts/l_ratio_p_vals.pdf')
  plot( theta.range, range(c(0,m.p)), type='n', xlab='Mean interest-growth differential',
        ylab='p-value' )
  v.col <- c('black', 'red', 'blue', 'darkgreen', 'brown' )
  for(i in 1:n.cty){
    if( v.cty[i] != 'GBR' ){
      lines( v.theta, m.p[,i], lwd=2, col=v.col[i] )
    } # Skip the UK.  Answer looks wrong.
  }
  legend('topright', v.cty, lwd=2, col=v.col, bty='n')
  abline(h=c( .025, .05, .1), lty=2)
dev.off()
