IAC <-
function(zz){
   mIAC <- cor(zz, method="p",use="p")
   ## mIAC <- stats::cor(zz, method="p",use="p")
   diag(mIAC) <- 0
   ADJ <- ((1+mIAC)/2)^2
   diag(ADJ) <- 0
   res <- structure(mIAC, ADJ=ADJ, class='IAC')
   return(res)
}

