
# converts each df column in cols to factor, returns new df

toFactors <- function(df,cols)
{  
   for (i in cols) {
      df[,i] <- as.factor(df[,i])
   }
   df
}

