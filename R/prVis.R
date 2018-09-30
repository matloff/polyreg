
# two-dimensional visualization of the X data in classification problems,
# similar in spirit to ordinary PCA and t-sne, color coded by Y values
# (class IDs in classification case, subinternvals of Y in regression
# case)

# t-sne, e.g. in the Rtsne package, applied in dimension k, attempts to
# find a k-dimensional manifold for which most of the data are "near";
# for visualization purposes, typically k = 2, which is assumed here

# the idea here is to expand the data with polynomial terms, using
# getPoly(), then apply PCA to the result

# typically these methods are applied only to a subsample of the data,
# due both to the lengthy computation time and the "black screen
# problem" (having a lot of points fills the screen, rendering the plot
# useless)

# arguments:
 
#    xy:  data frame
#    labels:  if TRUE, last column is Y for a classification problem;
#             must be an R factor, unless nIntervals is non-NULL, in
#             which case Y will be discretized to make labels
#    deg:  degree of polynomial expansion
#    scale:  if TRUE, first call scale() on the X data
#    nSubSam:  number of rows to randomly select; 0 means get all
#    nIntervals: in regression case, number of intervals to use for
#                partioning Y range to create labels
#    saveOutputs: if TRUE, return list with gpOut = output of getPoly(), 
#                 prout = output of prcomp()
#    cex: argument to R plot(), controlling point size

prVis <- function(xy,labels=FALSE,deg=2,scale=FALSE,nSubSam=0,nIntervals=NULL,
   saveOutputs=FALSE,cex=0.5)
{  
  nrxy <- nrow(xy)
  ncxy <- ncol(xy)

  rns <- row.names(xy)
  if (scale) {
     if (labels) {
        xy[,-ncxy] <- scale(xy[,-ncxy])
     } else xy <- scale(xy)
     row.names(xy) <- rns
  }

  if (nSubSam < nrxy && nSubSam > 0)  
     xy <- xy[sample(1:nrxy,nSubSam),]

  if (labels) {
     ydata <- xy[,ncxy]
     if (is.null(nIntervals) && !is.factor(ydata))
        stop('Y must be a factor for classif.; set nIntervals for regress.')
     if (!is.null(nIntervals)) {
       rng <- range(ydata)
       increm <- (rng[2] - rng[1]) / nIntervals
       ydata <- round((ydata - rng[1]) / increm)
       ydata <- as.factor(ydata)
     }
     xdata <- xy[,-ncxy, drop=FALSE]
  } else xdata <- xy

  xdata <- as.matrix(xdata)
  polyMat <- getPoly(xdata, deg)$xdata
  x.pca <- prcomp(polyMat,center=TRUE)
  if (labels)  {
     plot(x.pca$x[,1:2], col=ydata, pch=15, cex=0.5) 
   } else plot(x.pca$x[,1:2], pch=15, cex=0.5)
  
  if (saveOutputs) 
     return(list(gpOut=polyMat,prout=x.pca))
}

# intended to be used when a plot produced by prVis() is on the screen;
# chooses np points at random from the PCA output, writing their row
# numbers on the plot; these are the numbers from the full dataset, even
# if nSubSam > 0; the argument savedPrVisOut is the return value of
# prVis()

addRowNums <- function(np,savedPrVisOut) 
{
   pcax <- savedPrVisOut$prout$x[,1:2]
   if(is.null(row.names(pcax))) 
      stop('no row names')
   npcax <- nrow(pcax)
   tmp <- sample(1:npcax,np,replace=FALSE)
   rowNames <- row.names(pcax[tmp,])
   print('highlighted rows:')
   sorted <- sort(as.numeric(rowNames))
   for (i in 1:length(rowNames)) {
      rn <- rowNames[i]
      print(sorted[i])
      coords <- pcax[rn,]
      text(coords[1],coords[2],rn)
   }
}

