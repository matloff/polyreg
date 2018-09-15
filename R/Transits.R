
# image classification of B&W images, based on transitions from 0 to 1
# or vice versa

# arguments:

# img:  image data frame, one row per pixel row, class ID in last col
# deg:  desired polynomial degree
# maxBright: max pixel value 
 
# value:
 
# S3 object of class 'transits', to be fed into predict.transits()

transits <- function(img,deg,maxBright=255) 
{
   nr <- nrow(img)
   nc <- ncol(img)

   # convert to 0,1
   inc <- img[,-nc]
   inc <- round(inc / maxBright)

   tmp <- apply(inc,1,getNewX)

   xformedImg <- cbind(tmp,img)

}

getNewX <- function(oneImg,nr) 
{
   # get transition counts
   inc <- matrix(matrix(oneImg,byrow=TRUE,nrow=nr)
   tmp <- getFlips(inc)

   # form new X vector for this image 
   tmph <- sapply(tmp$h,function(z) z)
   tmpv <- sapply(tmp$v,function(z) z)
   c(tmph,tmpv)
}

# finds the transition counts for this image, 'inc'; returns an R list,
# one element per number of flips found (horiz. or vert.)

getFlips <- function(inc) 
{
   flips <- function(x) x[-length(x)]!= x[-1]
   flipCount <- function(x) sum(flips(x))
   horiz <- apply(inc,1,flipCount)
   vert <- apply(inc,2,flipCount)
   list(h = as.list(table(horiz)), v = as.list(table(vert)))
}
