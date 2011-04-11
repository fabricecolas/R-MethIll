plot.screenBTrait <-
function(x, p=NULL, ...){
   pCol <- attr(x,'p.name')
   if(is.null(p))
      p <- attr(x, 'p.bonf')
   df <- cbind(na.omit(print(x)),y=NA)
   df[df$kruskal.pval > p,pCol] <- 1
   df[,pCol] <- -1 * log10(df[,pCol])
   chrList <- unique(df$chr)
   chrList <- c(sort(chrList[grep('^[0-9]$',chrList)]), 
                sort(chrList[grep('^[1-2][0-9]$',chrList)]))
   chrList <- na.omit(c(chrList, unique(df$chr)[!(unique(df$chr) %in% chrList)]))
   chrSeq <- 1:length(chrList)
   df$kruskal.pval <- df$kruskal.pval-min(df$kruskal.pval) 
   df$kruskal.pval <- .90*df$kruskal.pval/max(df$kruskal.pval)
   for(i in chrSeq)
      df[df$chr %in% chrList[i],'y'] <- i+df[df$chr %in% chrList[i],'kruskal.pval']
   df <- df[order(df$loc),]
   cSel <- brewer.pal(9,'RdBu')[7]
   fName <- 'screenBTrait.pdf'
   ## get the location range for each chromosome (serves to plotting) 
   chrLocRange <- t(sapply(chrList, function(zz) return(range(df[df$chr == zz,'loc'], na.rm=TRUE))))
   ## redirect output to pdf
   pdf(file=fName, paper='letter')
   plot(NULL, xlim=c(min(chrLocRange), max(chrLocRange)), ylim=range(chrSeq), bty='n', axes=FALSE, xlab='position',
      ylab='chromosome', main='association w/ SampleLabel', sub=paste2('sites whose p<', signif(p,2)), cex.main=.8)
   axis(2, at=chrSeq, labels=chrList, las=2, cex=.8) ; axis(1, las=2,cex=.8)
   sapply(chrSeq, function(zz){lines(x=df[df$chr==chrList[zz],'loc'], y=df[df$chr==chrList[zz],'y'], type='s')})
   ## close pdf device
   graphics.off()
   return(invisible(fName))
}

