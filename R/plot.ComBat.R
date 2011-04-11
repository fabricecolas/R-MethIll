plot.ComBat <-
function(x, ...){
   gamma.hat <- attr(x, 'gamma.hat')
   gamma.bar <- attr(x, 'gamma.bar')
   delta.hat <- attr(x, 'delta.hat')
   a.prior <- attr(x, 'a.prior')
   b.prior <- attr(x, 'b.prior')
   par.prior <- attr(x, 'par.prior')
   t2 <- attr(x, 't2')

   if (par.prior){
      par(mfrow=c(2,2))
      tmp <- density(gamma.hat[1,])
      plot(tmp,  type='l', main="Density Plot")
      xx <- seq(min(tmp$x), max(tmp$x), length=100)
      lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
      qqnorm(gamma.hat[1,])   
      qqline(gamma.hat[1,], col=2)   
   
      tmp <- density(delta.hat[1,])
      invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
      tmp1 <- density(invgam)
      plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
      lines(tmp1, col=2)
      qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')   
      lines(c(0,max(invgam)),c(0,max(invgam)),col=2)   
      title('Q-Q Plot')
   }
   else
      cat('\nNo plot for non-parametric adjustment.\n')
}

