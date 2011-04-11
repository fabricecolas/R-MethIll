expDesign <-
function(x, gr1=NULL, gr2=NULL, latex=FALSE, exclude=NULL){
   if(class(x) == 'MethIll')
      return(expDesign(attr(x, 'eSet'), gr1=gr1, gr2=gr2, latex=latex))
   else if(class(x) == 'MethyLumiSet')
      return(expDesign(pData(x), gr1=gr1, gr2=gr2, latex=latex))
   else{
      if(is.null(gr1) & is.null(gr2))
         gr1 <- colnames(x)
      if(!is.null(gr1) & is.null(gr2))
         m <- t(combn(gr1, 2))
      if(!is.null(gr1) & !is.null(gr2))
         m <- as.matrix(expand.grid(gr1,gr2))
      res <- list()
      for(i in 1:nrow(m)){
         p <- chisq.test(x[,m[i,1]], x[,m[i,2]], simulate.p.value=TRUE, B=2000)$p.value
         res[[i]] <- table(x[,c(m[i,1], m[i,2])], exclude=exclude)
         if(is.null(exclude)){
            row.names(res[[i]])[is.na(row.names(res[[i]]))] <- 'NA'
            colnames(res[[i]])[is.na(colnames(res[[i]]))] <- 'NA'
         }
         names(res)[i] <- paste2(m[i,1],'-', m[i,2])
         # new
         res[[i]] <- res[[i]][hclust(dist(res[[i]]))$order, hclust(dist(t(res[[i]])))$order]
         if(latex){
            cap <- paste2('Joint distribution of \\textbf{', m[i,1], '} and \\textbf{',m[i,2],'}, $p_{\\chi^2}=',
               signif(p, digits=-1),'$.')
            print(xtable(res[[i]], cap=cap), table.placement="!ht",tabular.environment="longtable",
               latex.environments=c("center", "footnotesize"), floating=FALSE)
         }
         else
            print(res[[i]])
      }
      return(invisible(res))
   }
}

