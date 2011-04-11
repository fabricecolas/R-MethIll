MLearnExp <-
function(x, zdat, cl='ccstatus', mc.cores=2, NH0=100){
   if(mc.cores >= 2) require(multicore)
   tmp <- list()
   if(class(zdat) == 'MethIll') zdat <- attr(zdat, 'eSet')
   tmp[['beta']] <- t(exprs(zdat))
   tmp[['meth']] <- t(assayData(zdat)[['methylated']])
   tmp[['unmeth']] <- t(assayData(zdat)[['unmethylated']])
   ### matrix|Class 
   tmp <- lapply(tmp, function(ww) cbind(ww[row.names(pData(zdat)),], as.data.frame(pData(zdat)[,cl])))
   NC <- ncol(tmp[['beta']])
   colnames(tmp[['beta']])[NC] <- colnames(tmp[['unmeth']])[NC] <- colnames(tmp[['meth']])[NC] <- 'class'
   ### MLearn 
   m <- apply(x, 1, MLearnFun, dat=tmp, mc.cores=mc.cores)
   m <- t(m)
   ### H0
   mH0 <- NULL
   for(i in 6013:(6013+NH0)){
      set.seed(i)
      y <- attr(confuMat2(table(pData(datDP)[,'ccstatus'], sample(pData(datDP)[,'ccstatus']))),'measure')
      if(is.null(mH0)) mH0 <- y
      else mH0 <- rbind(mH0, y)
   }
   ### bind mH0 and m
   m <- rbind(m, matrix(NA, nrow(mH0), ncol(m)))
   m[(nrow(m)-nrow(mH0)+1):nrow(m),colnames(mH0)] <- mH0
   m <- as.data.frame(m)
   for(i in grep('(Recall)|(Precision)|(F1)|(Time)',colnames(m)))
      m[,i] <- as.numeric(as.character(m[,i]))
   row.names(m) <- 1:nrow(m)
   out <- structure('MLearnExp', results=m, class=cl, NH0=NH0, betaMD5=digest(tmp[['beta']]), methMD5=digest(tmp[['meth']]), unmeth=digest(tmp[['unmeth']]),
      date=format(Sys.time(), "%Y-%M-%d %H:%M:%S"), mc.cores=mc.cores, class='MLearnExp')
   return(out)
}

