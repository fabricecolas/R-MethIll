summary.confuMat2 <-
function(object, ...){ 
   cat('\n\t\t\tmaF1:', round(attr(object, 'maF1'), d=2))
   for(i in names(attr(object,'meas'))){
      out <- attr(object,'meas')[[i]]
      cat('\n. ', i, ':\t', round(out$F1, d=2),'\t( prec: ', round(out$prec, d=2), ',\trec:', round(out$rec, d=2), ')')
   }
   cat('\n')
   return(invisible(attributes(object)))
}

