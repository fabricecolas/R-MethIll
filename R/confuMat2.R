confuMat2 <-
function(x, clSet=NULL){
   x1 <- x 
   if(class(x) != 'table'){
      require(MLInterfaces)
      x1 <- confuMat(x)
   }
   if(is.null(clSet))
      clSet <- colnames(x1)
   NClass <- length(clSet)
   cName <- c(apply(expand.grid(c('Precision','Recall','F1'),clSet),1,paste,collapse='.'), 'maF1')
   m <- matrix(NA, 1 , length(cName), dimnames=list(NULL, cName))
   if(dim(x1)[1] < NClass){
      x1 <- rbind(x1, matrix(0, NClass-dim(x1)[1], ncol(x1)))
      row.names(x1)[which(!(row.names(x1) %in% clSet))] <- colnames(x1)[which(!(clSet %in% row.names(x1)))]
   }
   if(dim(x1)[2] < NClass){
      x1 <- cbind(x1, matrix(0, nrow(x1), NClass-dim(x1)[2]))
      colnames(x1)[which(!(colnames(x1) %in% clSet))] <- row.names(x1)[which(!(clSet %in% colnames(x1)))]
   }
   #else{
   #   x1 <- reshape(data.frame(x1), direction='wide', timevar='predicted', idvar='given')
   #   row.names(x1)<- x1[,1]
   #   x1 <- data.matrix(x1)[,-1]
   #   colnames(x1) <- gsub('Freq.', '', colnames(x1))
   #}
   res <- x1[clSet,clSet]
   meas <- list()
   for(i in (1:nrow(x1)-1)){
      tp <-     x1[  1+i ,  1+i]
      tn <- sum(x1[-(1+i),-(1+i)])
      fp <- sum(x1[-(1+i),  1+i])
      fn <- sum(x1[  1+i ,-(1+i)])
      prec <- tp/(tp+fp)
      if(is.nan(prec)) prec <- 0
      rec <- tp/(tp+fn)
      if(is.nan(rec)) rec <- 0
      f1 <-2 * prec * rec / (prec + rec)
      if(is.nan(f1)) f1 <- 0
      m[1,paste2('F1','.',clSet[i+1])] <- f1
      m[1,paste2('Precision','.',clSet[i+1])] <- prec
      m[1,paste2('Recall','.',clSet[i+1])] <- rec
   }
   m[1,'maF1'] <- mean(m[1,colnames(m)[grep('^F1',colnames(m))]])
   return(structure(res, measure=m, class='confuMat2'))
}

