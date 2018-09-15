
# image classification of B&W images, based on transitions from 0 to 1
# or vice versa

# arguments:

# img:  image data frame, one row per pixel row, class ID in last col
# nr, nc:  numbers of rows and cols of pixels in one image
# deg:  desired polynomial degree
# maxBright: max pixel value 
 
# value:
 
# S3 object of class 'transits', to be fed into predict.transits()

transits <- function(img,nr,nc,deg,maxFlips=10,maxBright=255) 
{
   # convert to 0,1
   inc <- img[,-nc]
   inc <- round(inc / maxBright)

   ## works for small number of rows, not large; "sapply() ATTEMPTS to
   ## create a matrix..."
   ## inct <- as.data.frame(t(inc))
   ## tmp <- sapply(inct,getNewX,nr,maxFlips)
   inct <- as.data.frame(t(inc))
   tmp <- lapply(inct,getNewX,nr,maxFlips)
   browser()
#       tmp <- NULL
#       for (i in 1:nrow(inc)) {
#          rw <- as.matrix(inc[i,])
#          tmp <- cbind(tmp,getNewX(rw,nr,maxFlips))
#       }
   xformedImg <- cbind(t(tmp),img[,nc])
   names(xformedImg)[nc] <- 'y'

   # now need to do the regression

}

getNewX <- function(oneImg,nr,maxFlips) 
{
   # get transition counts
   inc <- matrix(oneImg,byrow=TRUE,nrow=nr)
   tmp <- getFlips(inc,maxFlips)

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
   flips <- function(x) x[-length(x)]!= x[-1]
   flipCount <- function(x) sum(flips(x))
   horiz <- apply(inc,1,flipCount) / 2
   th <- table(horiz)
   h <- as.vector(th)
   names(h) <- names(th) 
   names(h) <- paste0('h',names(h))
   vert <- apply(inc,2,flipCount) / 2
   tv <- table(vert)
   v <- as.vector(tv)
   names(v) <- names(tv) 
   names(v) <- paste0('v',names(v))
   list(h = h, v = v)
}
