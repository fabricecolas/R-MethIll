ComBatBuildDesign <-
function(x, des=NULL, start=2){
   tmp <- matrix(0,length(x),nlevels(x)-start+1)
   for (i in 1:ncol(tmp))
      tmp[,i] <- x==levels(x)[i+start-1]
   return(cbind(des,tmp))
}

