
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

