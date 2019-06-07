
# time series routines

##############################  mov(): moving average ##########################

# the standard moving-average routines seem to include the given element
# itself, misleading in a prediction context

mov <- function(x,lag)
{
   lx <- length(x)
   res <- vector(length=lx-lag)
   for (i in 1:(lx-lag)) {                                                                         s <- i
      e <- i+lag-1
      res[i] <- mean(x[s:e])
   }                                                                                            res
}

