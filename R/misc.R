
# converts each df column in cols to factor, returns new df

toFactors <- function(df,cols)
{  
   for (i in cols) {
      df[,i] <- as.factor(df[,i])
   }
   df
}

# converts a single column to a matrix of dummy variables; column names
# will be concatenation of the original column names and levels(the
# column); column must be a factor

toDummies <- function(df,col)
{
   require(dummies)
   dumms <- dummies::dummy(df[,col])
   namesDumms <- paste0(names(df)[col],levels(df[,col]))
   tmp <- cbind(df[,-col],dumms)
   nctmp <- ncol(tmp)
   ncdumms <- ncol(dumms)
   names(tmp)[(nctmp-ncdumms+1):nctmp] <- namesDumms
   tmp
}

