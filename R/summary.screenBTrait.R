summary.screenBTrait <-
function(object, latex=FALSE, ...){
   p.bonf <- attr(object,'p.bonf')
   res <- attr(object,'res')
   out <- res[res[,'kruskal.pval'] < p.bonf,]
   if(latex & nrow(out) > 0)
      print(xtable(out,  digits=c(0,2,-1,-1,2,-1,-1,0,0,0), caption='Standard screening binary trait.'),  floating=FALSE,
         tabular.environment="longtable", latex.environments=c("center", "footnotesize"))
   else if(latex & nrow(out) == 0)
      cat('\\begin{itemize}\\item no probe with $p$ value of Kruskal and Wallis test below
         $p_{bonf}$\\end{itemize}')
   else
      print(out)
   return(invisible(out))
}

