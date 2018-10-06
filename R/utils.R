
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
         dumms <- dummy(dfi)
         outDF <- cbind(outDF,dumms[,-ncol(dumms),drop=FALSE])
      }
   }
   outDF[,1] <- NULL  # delete filler
   outDF
}

