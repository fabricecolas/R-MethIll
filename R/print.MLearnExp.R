`print.MLearnExp` <- 
function(x, latex=FALSE, rel2H0=FALSE, ...){
   m <- data.frame(attr(x, 'results'))
   cNameMeas <- grep('(Recall)|(Precision)|(F1)',colnames(m))
   out <- m[1:(nrow(m)-attr(x,'NH0')-1),]
   muH0 <- colMeans(m[is.na(m$DataSet),cNameMeas])
   if(!rel2H0){
      out <- rbind(out, muH0=NA)
      out['muH0', cNameMeas] <- muH0
   }
   else
      out[,cNameMeas] <-  100*(as.matrix(out[,cNameMeas]) / rep(1,nrow(out)) %*% t(muH0) -1)
   if(latex){
      id <- grep('(Recall)|(Precision)|(F1)|(Time)',colnames(out))
      out2 <- cbind(data.frame(out[,-id]),apply(out[,id],c(1,2), as.numeric))
      id2 <- grep('(Time)',colnames(out2))
      print(xtable(out2[,-id2], caption=paste2('Machine learning classification accuracy outults for experiment:',
         paste( colnames(out2[,-id]), collapse=', ')), digits=2), floating=FALSE, tabular.environment='longtable')
      return(invisible(out))
   }
   return(out)
}
