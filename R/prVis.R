
# two-dimensional visualization of the X data in classification problems,
# similar in spirit to ordinary PCA and t-sne, color coded by Y values
# (class IDs in classification case, subinternvals of Y in regression
# case)

# t-sne, e.g. in the Rtsne package, applied in dimension k, attempts to
# find a k-dimensional manifold for which most of the data are "near";
# for visualization purposes, typically k = 2, which is assumed heree

# the idea here is to expand the data with polynomial terms, using
# getPoly(), then apply PCA to the result

# typically these methods are applied only to a subsample of the data,
# due both to the lengthy computation time and the "black screen
# problem" (having a lot of points fills the screen, rendering the plot
# useless)

# arguments:
 
#    xy:  data frame, Y in last column; factor means it's a classification
#         problem
#    deg:  degree of polynomial expansion
#    nSubSam:  number of rows to randomly select; 0 means get all
#    nIntervals: in regression case, number of intervals to use for
#                partioning Y range
#    saveOutputs: if TRUE, return list with gpOut = output of getPoly(), 
#                  prout = output of prcomp()

prVis <- function(xy,deg=2,nSubSam=2000,nIntervals=NULL,saveOutputs=FALSE)
{ 
  if (nSubSam > 0) xy <- xy[sample(1:nrow(xy),nSubSam),]
  ydata <- xy[,ncol(xy)]
  if (is.null(nIntervals) && !is.factor(ydata))
     stop('Y must be a factor for classif.; set nIntervals for regress.')
  if (!is.factor(ydata)) {
    rng <- range(ydata)
    increm <- (rng[2] - rng[1]) / nIntervals
    ydata <- round((ydata - rng[1]) / increm)
  }

  xdata <- xy[,-ncol(xy), drop=FALSE]
  xdata <- as.matrix(xdata)
  polyMat <- getPoly(xdata, deg)$xdata
  x.pca <- prcomp(polyMat,center=TRUE)
  plot(x.pca$x[,1:2], col=ydata, pch=15, cex=0.5)
  if (saveOutputs) 
     return(list(gpOut=polyMat,prout=x.pca))
}

