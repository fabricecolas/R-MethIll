ComBatInt <-
function(x, pdat, batch, covariate, strata){
   if(class(x) == 'ExpressionSet'){
      covs <- covariate[!(covariate %in% union(batch,strata))]
      out <- lapply(c('methylated','unmethylated'), function(ww) ComBatInt(assayData(x)[[ww]], pdat=pData(x),
         batch=batch, covariate=covs, strata))
      llog <- list(methylated=out[[1]]$llog, unmethylated=out[[2]]$llog)
      assayData(x)[['methylated']] <- out[[1]]$data[1:nrow(out[[1]]$data),1:ncol(out[[1]]$data)]
      assayData(x)[['unmethylated']] <- out[[2]]$data[1:nrow(out[[2]]$data),1:ncol(out[[2]]$data)]
      exprs(x) <- assayData(x)[['methylated']]/(assayData(x)[['methylated']]+assayData(x)[['unmethylated']]+100)
      storageMode(x) <-'lockedEnvironment'
   }
   else{
      llog <- list()
      m <- matrix(NA,1,2,dimnames=list(list(), c('strata','batch')))
      if(!is.null(strata) & !is.null(batch))
         m <- na.omit(expand.grid(strata=unique(pdat[,strata]), batch=batch))
      else if(is.null(strata) & !is.null(batch)){
         m <-  matrix(NA,length(batch),2,dimnames=list(list(), c('strata','batch')))
         m[,2] <- batch
      }
      ## for each combination of 'strata' and 'batch'
      for(k in 1:nrow(m)){
         batchk <- as.character(m[k,'batch'])
         if(is.null(batch))
            batchk <- NULL 
         ## eventually subset to strata
         pdatk <- pdat
         if(!is.na(m[k,1]))
            pdatk <- pdat[pdat[,strata] %in% m[k,1],]
         ## we drop those samples in the given strata, which have 'unique' batch
         btab <- table(pdatk[,as.character(m[k,2])])
         pdatk <- pdatk[pdatk[,as.character(m[k,2])] %in% names(btab[btab > 1]),]
         ## subset the samples
         snk <- row.names(pdatk)
         llog[[k]] <- tmp <- tryCatch(ComBat(x[,snk], pdat=pdatk, batch=batchk, covariate=covariate),
            error=function(e) NULL) 
         if(!is.null(tmp))
            x[,snk]<- tmp[1:nrow(tmp),1:ncol(tmp)] 
      }
   }
   return(list(data=x,llog=llog))
}

