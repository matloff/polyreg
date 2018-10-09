
# for each column in df, if factor then replace by dummies, all but last
# class, else just copy column

factorsToDummies <- function(df) 
{
   require(dummies)
   outDF <- data.frame(rep(0,nrow((df))))  # filler start
   for (i in 1:ncol(df)) {
      dfi <- df[,i]
      if (!is.factor(dfi)) {
         outDF <- cbind(outDF,dfi) 
         names(outDF)[ncol(outDF)] <- names(df)[i]
      } else {
         dumms <- factorToDummies(dfi,names(df)[i])
         outDF <- cbind(outDF,dumms)
      }
   }
   outDF[,1] <- NULL  # delete filler
   outDF
}

factorToDummies <- function (f,fname) 
{
    n <- length(f)
    fl <- levels(f)
    ndumms <- length(fl) - 1
    dms <- matrix(nrow = n, ncol = ndumms)
    for (i in 1:ndumms) dms[, i] <- as.integer(f == fl[i])
    colnames(dms) <- paste(fname,'.', fl[-(length(fl))], sep = "")
    dms
}


