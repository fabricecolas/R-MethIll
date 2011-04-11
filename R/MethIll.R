MethIll <-
function(x, batch=NULL, covariate=NULL, strata=NULL, asfactor=NULL, filter=NULL, latex=FALSE, impute=TRUE){
   if(class(x) != 'MethIll'){
      metrics <- matrix(NA, 0, 4, dimnames=list(list(), c("Group", "Round", "Samples", "Mean_IAC")))
      ## .. MAKE A LOCAL COPY OF THE EXPR SET 
      storageMode(x) <- 'lockedEnvironment'
      y <- x
      storageMode(y) <- 'environment'
      ## .. FILTER DATA ..
      ## select probes having a variance different to zero
      probeVar <- apply(exprs(y), 1, var, na.rm=TRUE)
      probeVar <- ifelse(is.na(probeVar), 0, probeVar)
      y <- y[probeVar>0,]
      print(paste2("Note: ",sum(probeVar == 0)," probe sets had 0 variance"))
      log <- list()
   }
   else if(class(x) == 'MethIll' & !is.null(filter))
      return(MethIllIACFilter(x, filter=filter))
   else{
      y <- attr(x, 'eSet')
      storageMode(y) <- 'environment'
      if(is.null(asfactor))
         asfactor <- attr(x, 'asfactor')
      attr(x,'batch') <- batch 
      attr(x,'covariate') <- covariate
      attr(x,'strata') <- strata 
      log <- attr(x,'log') 
   }
   if(is.null(asfactor))
      asfactor <- union(batch, covariate)
   else
      asfactor <- asfactor[(asfactor %in% union(covariate, c(batch,strata)))]
   ## .. IMPUTATION: IF NA OCCUR AND IMPUTE IS TRUE
   for(i in c('methylated','unmethylated')){
      hasNA <- sum(is.na(assayData(y)[[i]])) > 0
      if (hasNA & !impute) 
         stop(paste2("Missing data in assay data ", i," ...")) 
      else if (hasNA & impute){
         print(paste2("Perform imputation of assay data ", i))
         assayData(y)[[i]] <- impute.knn(assayData(y)[[i]])$data
      }
   }
   ## .. CANCEL BATCH EFFECT AND UPDATE
   if(!is.null(batch)){
      cbat <- ComBatInt(y, batch=batch, covariate=covariate, strata=strata)
      attr(x, 'eSet') <- cbat$data
      attr(x, 'log') <- c(attr(x, 'log'), cbat$log)
      attr(x, 'asfactor') <- asfactor
   }
   ## LOCK THE ENVIRONMENT OF THE EXPRESSION SET
   if(class(x) != 'MethIll'){
      storageMode(y) <-'lockedEnvironment'
      ## .. INIT STRUCTURE
      x <- structure("A sample Network", batch=batch, covariate=covariate, strata=strata, asfactor=asfactor,
         metrics=metrics, impute=impute, IAC=NULL, IAC.history=NULL, eSet=y, eSetOrig=y, Z.K.IAC=NULL, Z.CC=NULL,
         Z.MAR=NULL, log=log, round_log=list(), cex=.6, R.Version=R.Version(), class='MethIll')
   }
   ## .. MAKE A 1ST ROUND OF SAMPLE CORRELATIONS AND RETURN.. 
   return(MethIllCor(x))
}

