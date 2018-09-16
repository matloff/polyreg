
# image classification of B&W images, based on transitions from 0 to 1
# or vice versa

#####################  transitFit() and helpers  ########################

transitFit <- function() 
{

}

#####################  transits() and helpers  ########################

# arguments:

# img:  image data frame, one row per pixel row, class ID in last col
# nr, nc:  numbers of rows and cols of pixels in this image
# maxFlips: max number of anticipated  0 -> 1 transitions for this,
#           image or col names pattern (see below)
# maxBright: max pixel value 
# newCases: if TRUE, we are predicting new cases
 
# value:
 
# data frame; each row is the transformed version of the corresponding
# row in the original, showing how many rows had i 0 -> 1 transitions,
# for various values of i, and for horiz. and vert. directions; the last
# column has the original class IDs, unless newCases is TRUE

# details on maxFlips: if maxFlips is, say, 10, it may be that no row or
# col in the image has anymore than, say, 4, transitions; this will
# produce columns of all-0s in the output, which are deleted; but we
# must take care to delete the same ones when predicting new cases; this
# is done by saving the column names of the transformed training data,
# and providing them as maxFlips in predicting new data

transits <- function(img,nr,nc,maxFlips=10,maxBright=255,newCases=FALSE) 
{
   nColsImg <- ncol(img)
   if (!( 
          (!newCases && nColsImg == nr*nc+1) || 
           (newCases && nColsImg == nr*nc)  
        )) stop('incorrect number of columns in input')
   # convert to 0,1
   inc <- if (!newCases) img[,-(nr*nc+1)] else img
   inc <- round(inc / maxBright)

   ## works for small number of rows, not large; "sapply() ATTEMPTS to
   ## create a matrix..."
   inct <- as.data.frame(t(inc))
   tmp <- sapply(inct,getNewX,nr,maxFlips,newCases)
   xformedImg <- as.data.frame(t(tmp))
   if (!newCases) {
      xformedImg <- cbind(xformedImg,img[,nr*nc+1])
      names(xformedImg)[ncol(xformedImg)] <- 'y'
   }

   # may have some 0 cols, due to overly conservative setting of
   # maxFlips; in new cases setting, need to use the cols found during
   # the training set phase
   if (!newCases) {
      all0 <- function(col)  all(col == 0)
      xformedImg1 <- xformedImg[,-ncol(xformedImg)]
      zeroCols <- which(apply(xformedImg1,2,all0))
   } else {
      allNames <- names(xformedImg)
      zeroCols <- setdiff(allNames,maxFlips)
   }
      xformedImg[zeroCols] <- NULL
   xformedImg
}

getNewX <- function(oneImg,nr,maxFlips,newCases) 
{
   # get transition counts
   if (newCases) maxFlips <- length(maxFlips)
   inc <- matrix(oneImg,byrow=TRUE,nrow=nr)
   tmp <- getFlips(inc,maxFlips)
   if (length(tmp) > maxFlips) 
      stop('maxFlips set too small')

   # form new X vector for this image
   maxFlips1 <- maxFlips + 1
   h <- rep(0,maxFlips1)
   names(h) <- paste0('h',as.character(0:maxFlips))
   h[names(tmp$h)] <- tmp$h
   v <- rep(0,maxFlips1)
   names(v) <- paste0('v',as.character(0:maxFlips))
   v[names(tmp$v)] <- tmp$v
   c(h,v)
}

# finds the transition counts for this image, 'inc'; returns an R list,
# one element per number of flips found (horiz. or vert.)

getFlips <- function(inc,maxFlips) 
{
   # flips <- function(x) x[-length(x)]!= x[-1]
   flips <- function(x) x[-1] == 1 & x[-length(x)] == 0
   flipCount <- function(x) sum(flips(x))
   horiz <- apply(inc,1,flipCount) 
   th <- table(horiz)
   h <- as.vector(th)
   names(h) <- names(th) 
   names(h) <- paste0('h',names(h))
   vert <- apply(inc,2,flipCount) 
   tv <- table(vert)
   v <- as.vector(tv)
   names(v) <- names(tv) 
   names(v) <- paste0('v',names(v))
   list(h = h, v = v)
}
