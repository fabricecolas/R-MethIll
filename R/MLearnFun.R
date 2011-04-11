MLearnFun <-
function(x, dat, mc.cores){
   ### data: true classification or permuted-H0 hypothesis
   dat <- dat[[x['DataSet']]]
   ### feature selection
   fsFun <- function(formula, data) formula 
   if(x['FeatSel'] == 'rand'){
      set.seed(6013) # identical definition of random set of features accross classification tasks
      dat <- dat[,c(sample(colnames(dat)[!(colnames(dat) %in% 'class')])[1:as.numeric(x['NFeat'])],'class')]
      x['DataMD5'] <- digest(dat)
   }
   else if(x['FeatSel'] == 'fs.absT'){ 
      fsFun <- fs.absT(as.numeric(x['NFeat']))
      x['DataMD5'] <- 'To Do'
   }
   set.seed(6013) # identical split in the cross validation
   t1 <- Sys.time()
   xClassifier <- MLearn(as.factor(class)~., dat , get(x['Algo']), xvalSpec("LOO", fsFun=fsFun), mc.cores=mc.cores)
   t2 <- Sys.time()
   out <- attr(confuMat2(xClassifier, clSet=sort(unique(dat[,'class']))), 'measure')
   x[colnames(out)] <- out
   x['Time'] <- as.numeric(difftime(t2, t1, units='secs')) 
   x['Date'] <- format(Sys.time(), "%Y-%M-%d %H:%M:%S")
   return(x)
}

