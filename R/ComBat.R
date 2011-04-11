ComBat <-
function(x, pdat, batch, covariate, par.prior=TRUE){
   pdat <- as.matrix(pdat)
   ## DESIGN MATRIX
   design <- ComBatDesignMat(pdat, batch, covariate)   

   ## LEVELS FOR EACH BATCH, NUMBER AND LENGTH OF BATCHES, ARRAY
   batches <- NULL
   b.fact <- as.factor(pdat[,batch])
   for (i in 1:nlevels(b.fact))
      batches <- append(batches, list((1:length(b.fact))[b.fact==levels(b.fact)[i]]))

   n.batch <- length(batches)
   n.batches <- sapply(batches, length)
   n.array <- sum(n.batches)
   
   ## Standardize Data across genes
   cat('Standardizing Data across genes\n')
   B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(x))
   grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
   var.pooled <- ((x-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)

   stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
   if(!is.null(design)){
      tmp <- design
      tmp[,c(1:n.batch)] <- 0
      stand.mean <- stand.mean+t(tmp%*%B.hat)
   }   
   s.data <- (x-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

   ## Get regression batch effect parameters
   cat("Fitting L/S model and finding priors\n")
   batch.design <- design[,1:n.batch]
   gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
   delta.hat <- NULL
   for (i in batches)
      delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=TRUE))

   ## Find Priors
   gamma.bar <- apply(gamma.hat, 1, mean)
   t2 <- apply(gamma.hat, 1, var)
   a.prior <- apply(delta.hat, 1, ComBatAPrior)
   b.prior <- apply(delta.hat, 1, ComBatBPrior)

   ## Find EB batch adjustments
   gamma.star <- delta.star <- NULL
   if(par.prior){
      cat("Finding parametric adjustments\n")
      for (i in 1:n.batch){
         temp <- ComBatItSol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
         gamma.star <- rbind(gamma.star,temp[1,])
         delta.star <- rbind(delta.star,temp[2,])
      }
   }
   else{
      cat("Finding nonparametric adjustments\n")
      for (i in 1:n.batch){
         temp <- ComBatIntEprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
         gamma.star <- rbind(gamma.star,temp[1,])
         delta.star <- rbind(delta.star,temp[2,])
      }
   }

   ### Normalize the Data ###
   cat("Adjusting the Data\n")
   bayesdata <- s.data
   j <- 1
   for (i in batches){
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
      j <- j+1
   }

   bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean

   res <- structure(bayesdata, gamma.hat=gamma.hat, delta.hat=delta.hat, gamma.bar=gamma.bar, a.prior=a.prior, b.prior=b.prior,
      par.prior=par.prior, t2=t2, class='ComBat')

   return(res)
}

