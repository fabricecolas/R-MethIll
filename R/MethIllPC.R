MethIllPC <-
function(x, G) {
   G.lev <- levels(factor(G))
   # module eigengenes (PCs), 
   pc <- data.frame(matrix(NA, nrow=dim(x)[[1]], ncol=length(G.lev), dimnames=list(list(), paste2("PC", G.lev)))) 
   # percent variance explained by the first 5 PCs of a module
   varE <- data.frame(matrix(NA, nrow=5, ncol=length(G.lev)))
   for(i in 1:length(G.lev)){
      print(i)   
      # in the following, rows are genes and columns are samples
      xG <- t(x[,as.character(G) == G.lev[i]])
      if ( length(xG)==3 && !is.null(names(xG)) && names(xG)[1]=="data" )
         xG <- xG$data
      xG <- t(scale(t(xG)))
      xSvd <- svd(xG)
      varE[,i] <- (xSvd$d[1:5])^2/sum(xSvd$d^2)
      # first principal component
      pc[,i] <- xSvd$v[,1]
      signh1 <- sign(sum(cor(pc[,i],  t(xG))))
      if (signh1 != 0)  
         pc[,i] <- signh1* pc[,i]
   }
   return(list(PrinComps=pc, varexplained=varE))
}

