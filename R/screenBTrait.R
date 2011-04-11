screenBTrait <-
function(x, btrait='SampleLabel', nnodes=1, csv=TRUE){
   if(class(x) == 'MethyLumiSet'){
      datE <- exprs(x)
      y <- pData(x)[,btrait]
      y.val <- levels(factor(y))
      y.fact <- ifelse(y == y.val[[1]], 1, ifelse(y == y.val[[2]], 2, NA))
   }
   else if(class(x) == 'MethIll')
      return(screenBTrait(attr(x,'eSet'), btrait=btrait, nnodes=nnodes, csv=csv))
   else ### stops based on conditioning
      stop('x must be either a MethIll or a MethyLumiSet.')
   if (length(y.val) > 2) 
      stop("The sample trait y contains more than 2 levels. Please input a binary variable y")
   if (length(y.val) == 1) 
      stop("The sample trait y is constant. Please input a binary sample trait with some variation.")
   if (length(y.fact) != ncol(datE)) 
      stop("the length of the sample trait y does not equal the number of rows of exprs(x)")
   fun.stats <- function(ww){m <- cor.test(y.fact, datE[ww,], use="p"); rcorr <- Hmisc::rcorr.cens(datE[ww,],
      y.fact,outx=TRUE)[[1]] ; k.p <- kruskal.test(datE[ww,]~factor(y.fact))$p.value ; list('cor.p'=m$estimate,
      'cor.pvalStudent'=m$p.value, 'cor.qvalStudent'=NA,'rcorr.cens'=rcorr, 'kruskal.pval'=k.p, 'kruskal.q'=NA)}
   if(nnodes > 1){
      require(snow)
      cl <- makeCluster(nnodes, type='MPI')
      res <- t(parSapply(cl, row.names(datE), fun.stats))
      stopCluster(cl)
   }
   else
      res <- t(sapply(row.names(datE), fun.stats))
   res[,'kruskal.q'] <- qvalue(as.numeric(res[,'kruskal.pval']))$qvalues
   res[,'cor.qvalStudent'] <- qvalue(as.numeric(res[,'cor.pvalStudent']))$qvalues
   res <- as.data.frame(apply(res,2,as.numeric))
   row.names(res) <- row.names(datE)
   res <- cbind(res, fData(x)[row.names(res),])
   res <- res[,!(colnames(res) %in% c('TargetID','ProbeID_A','ProbeID_B'))]
   res <- res[order(res[,'kruskal.pval']),]
   p.bonf <- .05/nrow(exprs(x))
   if(csv)
      write.table(res, file='standard_screening_binary_trait.csv')
   res <- structure('standard screening binary trait', res=res, btrait=btrait, nnodes=nnodes, p.bonf=p.bonf,
      p.name=c('cor.pvalStudent','kruskal.pval'), class='screenBTrait')
   return(invisible(res))
}

