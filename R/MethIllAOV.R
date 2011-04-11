MethIllAOV <-
function(PC, pdat, r.term){
   ## access directly the variable by their name
   ## we do it in a local env to avoid concurrency with local vars 
   if('ww' %in% colnames(pdat))
      stop('The column name \'ww\' is already in use for the local env ')
   attach(pdat, warn.conflicts = FALSE)
   res <- sapply(r.term, function(ww) anova(lm(as.formula(paste2("PC~", ww))))[1,5])
   detach(pdat)
   res <- ifelse(is.nan(res), 1, res)
   return(res)
}

