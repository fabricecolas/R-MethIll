ComBatItSol <-
function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
   n <- apply(!is.na(sdat),1,sum)
   g.old <- g.hat
   d.old <- d.hat
   change <- 1
   while(change>conv){
      g.new <- ComBatPostMean(g.hat,g.bar,n,d.old,t2)
      sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=TRUE)
      d.new <- ComBatPostVar(sum2,n,a,b)
      change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
      g.old <- g.new
      d.old <- d.new
   }
   adjust <- rbind(g.new, d.new)
   rownames(adjust) <- c("g.star","d.star")
   adjust
}

