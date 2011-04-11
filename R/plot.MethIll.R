plot.MethIll <-
function(x, covariate=NULL, ...){
   ## expression and sample/pheno data 
   datE <- exprs(attr(x, 'eSet'))
   pdat <- pData(attr(x,'eSet'))
   ## 
   mIAC <- attr(x, 'IAC')
   ADJ <- attr(x, 'ADJ')
   PC <- attr(x, 'PC')
   meanIAC <- attr(x, 'meanIAC')
   cex <- attr(x, 'cex')
   Z.K.IAC <- attr(x, 'Z.K.IAC')
   Z.CC <- attr(x, 'Z.CC')
   c1 <- attr(x, 'c1')
   batch <- attr(x,'batch') 
   if(is.null(covariate))
      covariate <- attr(x, 'covariate')
   if(is.null(covariate))
      covariate <- colnames(pdat)
   strata <- attr(x, 'strata')
   asfactor <- attr(x, 'asfactor')
   whichround <- print(x)[nrow(print(x)),'Round']
   model <- attr(x, 'model')
   fitMeasure <- attr(x, 'fit')

   ## build univariate linear model formula
   m.terms <- tr <- union(covariate, c(batch, strata))
   m.terms[tr %in% asfactor] <- lpaste(tr[tr %in% asfactor], "factor(",")")

   # contribution of each trait wrt PC
   tr.aov <- MethIllAOV(PC, pdat, m.terms)
                        
   # individual contribution of trait-levels wrt PC for factors
   tr.aovL <- list()
   for(f in asfactor){
      lev <- sort(names(which(table(pdat[, f])>1)))
      rt <- ifelse(any(!(tr %in% f)), paste2('+', paste(m.terms[!(tr %in% f)], collapse='+')), '')
      tr.aovL[[f]] <- sapply(lev, function(ww) MethIllAOV(PC, pdat, paste2("factor(",f,"=='",ww,"')", rt)))
      names(tr.aovL[[f]]) <- lpaste(lev, 'factor(',')')
   } 
            
   ## START MULTI-ROW/COLUMN PLOT
   f <- list()
   cSelect <- colorRampPalette(rev(brewer.pal(12,'Paired')))
   if(!file.exists('fig/'))
     dir.create('fig/')
   f <- paste2('fig/',format(Sys.time(),"%Y-%m-%d--%H-%M-%S-"), "connectivity.pdf")

   pdf(file=f, width=7, height=10)

   ## LAYOUT 
   n.traits <- length(tr)
   n.plots <- 4+n.traits
   n.col <- 3
   m.layout <- t(sapply(1:(1+length(tr)), function(ww) rep(ww, n.col)))
   h <- c(5, rep(1/8, n.traits)) 
   i <- 0
   while((n.plots%%(i+.1)) != n.plots){
      m.layout <- rbind(m.layout, (max(m.layout)+1):(max(m.layout)+n.col))
      h <- c(h, 5)
      i <- i+n.col
   }
   layout(m.layout, heights=h, width=rep(1, n.col))

   ## HIERARCHICAL CLUSTERING
   par(mar=c(0, 4.1, 3, 1.1))
   plot(c1, cex=.8*cex, main=paste2("mean IAC = ", signif(meanIAC,3)), xlab="", ylab="1 - IAC", bty='n', las=1)
   m <- apply(as.matrix(pdat[c1$order,tr]), 2, function(ww) as.numeric(as.factor(ww)))   
   ## IMAGE BY TRAIT 
   par(mar=c(0, 5.3, 0, 2.4))
   sapply(tr, function(ww) {image(as.matrix(m[,ww]), col=cSelect(length(unique(pdat[,ww]))), axes=FALSE, new=TRUE);
      text(0, 0, ww, cex=cex*.8, pos=2)})
   par(mar=c(3.9, 4.1, 2.5, 1.4))

   cSelName <- covariate[1]
   cSelVal <- cSelect(length(unique(m[,covariate[1]])))[m[,cSelName]]
   ## CONNECTIVITY DISTRIBUTION
   plot(Z.K.IAC, main="Connectivity distribution", ylab="Z.K.IAC",type="n", xaxt="n", xlab="", bty='n', las=1)
   mtext(paste2('Intra-Array Connectivity score (', cSelName,')'), cex=cex, line=0.2)
   text(Z.K.IAC, labels=row.names(pdat), cex=cex, col=cSelVal)
   abline(h=c(-2,-3))

   ## CLUSTER COEFFICIENT DISTRIBUTION
   plot(Z.CC, main="ClusterCoef distribution", ylab="Z.CC", type="n", xaxt="n", xlab="", bty='n', las=1)
   mtext(paste2('Cluster coefficient score (', cSelName,')'), cex=cex, line=0.2)
   text(Z.CC, labels=row.names(pdat), cex=cex, col=cSelVal)
   abline(h=c(2, 3, -2, -3))

   ## CONNECTIVITY VERSUS CLUSTER COEFFICIENT
   #plot(Z.K.IAC, Z.CC, main="Connectivity vs ClusterCoef", xlab="Z.K.IAC", ylab="Z.CC", col=cSelVal, bty='n', las=1)
   #abline(lm(Z.CC~Z.K.IAC),col="red")
   #mtext(paste2("r = ",signif(cor.test(Z.K.IAC,Z.CC)$estimate,3), " p = ",
   #   signif(cor.test(Z.K.IAC,Z.CC)$p.value,3)), line=0.1, cex=cex)

   ## UNIVARIATE ANOVA BARPLOT
   tr.aov[tr.aov==0]=1e-300
   barplot(-log(tr.aov, 10), names.arg=gsub('^factor','f',names(tr.aov)), las=3, cex.names=cex, 
      ylab="-log10 p-val", main=paste2("ANOVA of PC"), bty='n', border=1, las=2)
   mtext('univariate', cex=cex, line=0.2)
   abline(h=-log(.05,10), col="blue", lwd=2)
   ## BONFFERRONI
   abline(h=-log(.05/length(tr), 10), col="red", lwd=2)

   ## ANOVA P-VALUE FOR EACH FACTOR LEVEL OF THE TRAITS  
   for(v in names(tr.aovL)){
      aovL <- tr.aovL[[v]]
      aovL[aovL==0] <- 1e-300
      plot(-log(aovL,10), type="n", xlab='', ylab="-log10 p-val", main="ANOVA of PC", xaxt="n", bty='n', las=1)
      mtext(paste2('univariate -- ', v), cex=cex, line=0.2)
      cSelVal <- cSelect(length(unique(m[,v])))
      text(-log(aovL,10), labels=gsub("(^factor[(])|([)]$)",'',names(aovL)), col=cSelVal, cex=cex)
      abline(h=-log(.05,10), col="blue", lwd=2)
      abline(h=-log(.05/length(aovL),10), col="red", lwd=2)
   } 
   graphics.off()
   return(f)
}

