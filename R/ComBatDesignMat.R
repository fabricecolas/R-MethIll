ComBatDesignMat <-
function(x, batch, covariate){
   batch <- as.factor(x[,which(colnames(x) == batch)])
   cat("Found", nlevels(batch), 'batches\n')
   design <- ComBatBuildDesign(batch, start=1)
   cat("Found", length(covariate), 'covariate(s)\n')
   for (j in covariate)
      design <- ComBatBuildDesign(as.factor(x[,j]), des=design)
   return(design)
}

