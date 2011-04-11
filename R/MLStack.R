MLStack <-
function(Class, DataSet, FeatSel='rand', NFeat=as.character(as.integer(2^(.5+3:9))), Algo=c('svmI'), ...){
   ### Cartesian product of the different parameters
   res <- expand.grid(FeatSel=FeatSel, NFeat=NFeat, Algo=Algo, DataSet=DataSet, Class=Class, Precision=NA,
      Recall=NA, F1=NA, ...)
   ### Get {Precision, Recall, F1} for each Class
   res <- reshape(res, direction='wide', timevar='Class', idvar=colnames(res)[!(colnames(res) %in% c('Precision','Recall','F1'))])
   ### Append additional logging columns
   res <- cbind(res, maF1=NA, Time=NA, DataMD5=NA, Date=NA)
   ### get md5 of each experiment, and take unique instances
   res <- cbind(res, ExpMD5=sapply(apply(res[, c('FeatSel', 'NFeat', 'Algo', 'DataSet')], 1, paste2), digest))
   res <- res[match(unique(res$ExpMD5),res$ExpMD5),]
   ### reorder experiments by increasing number of features, and by algorithm
   res <- res[order(res$NFeat, res$Algo),]
   return(res)
}

