MethIllIACFilter <-
function(x, filter){
   sSel <- print(x)$Z.K.IAC>filter
   if(any(sSel)){
      sName <- row.names(print(x)[sSel,])
      ## FILTER OUT THE SAMPLES
      storageMode(attr(x, 'eSet')) <- 'environment'
      attr(x, 'eSet') <- attr(x, 'eSet')[,sName]
      storageMode(attr(x, 'eSet')) <- 'lockedEnvironment'
      ## RECALCULATE CORRELATION BETWEEN THE SAMPLES
      x <- MethIllCor(x)
   }
   return(x)
}

