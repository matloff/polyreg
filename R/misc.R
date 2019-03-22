
# converts each df column in cols to factor, returns new df
# for getPoly(), and by extension polyFit() and FSR():
# should be used on any categorical variable stored as an integer
# is optional for binary variables
# is optional for categorical variables stored as characters
#' @export
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
# SHOULD NOT BE USED with getPoly() or FSR()
#' @export
toDummies <- function(df,col)
{
   # require(dummies)
   # requireNamespace(dummies)
   dumms <- dummies::dummy(df[,col])
   namesDumms <- paste0(names(df)[col],levels(df[,col]))
   tmp <- cbind(df[,-col],dumms)
   nctmp <- ncol(tmp)
   ncdumms <- ncol(dumms)
   names(tmp)[(nctmp-ncdumms+1):nctmp] <- namesDumms
   tmp
}

