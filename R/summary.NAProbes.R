summary.NAProbes <-
function(object, latex=FALSE, ...){
   hasNA <- signif(100*attr(object,'hasNA'), 2)
   naRate <- signif(100*attr(object,'naRate'), 2)
   if(latex){
      caption <- paste2('This table reports the number of probes in each chromosome that present none (0)/at
         least one (1) missing value. With $p_{\\chi^2}=',signif(attr(object,'pval'), digits=2),'$, we show that there
         is no apparent chromosome bias in missingness. We also note that ',hasNA,'\\% of all probes show at least
         one missing value, whereas for the assay, the missing value rate is ',naRate,'\\%.')
      print(xtable(attr(object,'byChr'), cap=caption, label='tab:nachr'), floating=FALSE,
         tabular.environment="longtable", latex.environments=c("center", "footnotesize"))
   }
   else
      print(attr(object,'byChr'))
}

