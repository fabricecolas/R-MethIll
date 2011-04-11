NAProbes <-
function(x, latex=FALSE){
   ## identify probes presenting NA-values
   naSumByProbe <- apply(is.na(exprs(x)), 1, sum)
   ## table of the number of NA-probes by chromosome
   naByProbe <- cbind(fData(x), hasNA=0)
   naByProbe[row.names(naByProbe) %in% names(which(naSumByProbe>0)), 'hasNA'] <- 1
   hasNA <- sum(naByProbe$hasNA)/nrow(exprs(x))
   naByChr <- table(naByProbe[,c('chr','hasNA')])
   naByChr <- naByChr[order(as.numeric(row.names(naByChr))),]
   if(is.null(dim(naByChr)))
      naByChr <- cbind(naByChr, isNA=0)   
   pval <- chisq.test(naByChr, simulate.p.value=TRUE, B=2000)$p.value
   naRate <- table(is.na(exprs(x)))[["TRUE"]]/(nrow(exprs(x))*ncol(exprs(x)))
   res <- structure(naSumByProbe, byChr=naByChr, pval=pval, hasNA=hasNA, naRate=naRate, class='NAProbes')
   return(res)
}

