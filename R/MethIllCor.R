MethIllCor <-
function(zz, gdata='beta'){
   ## calculates intra array connectivity
   print("Calculating sample correlations...")
   if(gdata == 'beta')
      datE <- exprs(attr(zz, 'eSet'))
   else
      datE <- assayData(attr(zz, 'eSet')[[gdata]])
   pdat <- pData(attr(zz, 'eSet'))

   attr(zz, 'IAC') <- mIAC <- IAC(datE)

   ## fundamental network concepts
   attr(zz, 'FNC') <- FNC <- WGCNA::fundamentalNetworkConcepts(attr(mIAC,'ADJ'))
   ## scaled connectivity 
   attr(zz, 'kIAC2') <- kIAC2 <- FNC$ScaledConnectivity
   attr(zz, 'PC.all') <- PC.all <- MethIllPC(t(datE), rep("gold", nrow(datE)))
   attr(zz, 'PC') <- PC <- PC.all$PrinComps$PCgold
   ## connectivity and cluster coefficient score
   attr(zz,'Z.K.IAC') <- Z.K.IAC <- (kIAC2-mean(kIAC2))/sd(kIAC2)
   attr(zz,'Z.CC') <- Z.CC <- (FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)

   ## metrics log-data-frame
   m <- rbind(attr(zz,'metrics'), NA)
   m[nrow(m),c('Round','Samples')] <- c(nrow(m), ncol(datE))
   m[nrow(m),'Mean_IAC'] <- attr(zz,'meanIAC') <- mean(mIAC[upper.tri(mIAC)], na.rm=TRUE)
   attr(zz,'metrics') <- m

   ## hierarchical clustering 
   attr(zz, 'c1') <- c1 <- hclust(as.dist(1-mIAC), method="average")

   ## metrics for that round: make, store as current and in the log, print 
   tr <- union(attr(zz,'covariate'),c(attr(zz,'batch'), attr(zz,'strata')))
   rd <- data.frame(pdat[,tr], K.IAC=signif(kIAC2, 3), Z.CC=signif(Z.CC, 3), Z.K.IAC=signif(Z.K.IAC, 3))
   attr(zz,'round_log')[[nrow(m)]] <- attr(zz,'round_metrics') <- rd[order(rd[,"Z.K.IAC"], decreasing=TRUE),]
   return(zz)
}

