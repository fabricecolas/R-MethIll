lpaste <-
function(x, pre="pdat[,'", suf="']"){ 
   x <- lapply(x, function(ww) paste(pre, ww, suf, sep=''))
   return(unlist(x))
}

