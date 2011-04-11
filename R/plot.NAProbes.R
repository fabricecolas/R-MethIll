plot.NAProbes <-
function(x, file='NAProbesBarplot', ...){
   ## histogram showing the distribution of NA by probs
   f <- paste2(file, '.pdf')
   pdf(f)
   barplot(log10(table(x[1:length(x)])), bty='n', axes=FALSE, border=0)
   axis(2, at=c(0,2,4), label=10^c(0,2,4))
   graphics.off()
   return(invisible(f))
}

