print.MethIll <-
function(x, latex=FALSE,...){
   rd <- attr(x, 'round_metrics')
   rd_n <- nrow(attr(x, 'metrics'))
   if(latex){
      cap <- paste2('At round ', rd_n, ', ', sum(rd$Z.K.IAC< -3),' sample(s) with ')
      cap <- paste2(cap, 'connectivity score(s) lower than $z_i^{kIAC}<-3$ and ',sum(rd$Z.K.IAC> -3 & rd$Z.K.IAC< -2)) 
      cap <- paste2(cap, ' with scores between $-3<z_i^{kIAC}<-2$, $i\\in \\{1..N\\}$ referring to the samples.') 
      print(xtable(as.matrix(rd), cap=cap), floating=FALSE, tabular.environment="longtable",
         lalatex.environments=c("center", "footnotesize"))
   }
   else
      return(rd)
}

