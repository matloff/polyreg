Tests on nytaxi
================

Data preprocessing
------------------

``` r
library(readr)
library(chron)
library(caret)
# train.csv is downloaded from https://www.kaggle.com/c/nyc-taxi-trip-duration/data
nytaxi <- read_csv("~/Desktop/project/train.csv")
nytaxi$pickupDate <- as.Date(nytaxi$pickup_datetime)
nytaxi$pickupTime <- as.factor(format(nytaxi$pickup_datetime, "%H"))
nytaxi$pickupWeekday <- ifelse(is.weekend(nytaxi$pickupDate), 0, 1)
nytaxi$dropoffDate <- as.Date(nytaxi$dropoff_datetime)
nytaxi$dropoffTime <- as.factor(format(nytaxi$dropoff_datetime, "%H"))
nytaxi$dropoffWeekday <- ifelse(is.weekend(nytaxi$dropoffDate), 0, 1)
nytaxi <- nytaxi[, -c(3,4)]
nytaxi <- subset(nytaxi, select = -c(pickupDate,dropoffDate))
# divide time to 4 classes (00-05, 06-11, 12-17, 18-23)
# create 3 dummy variables
check1 = nytaxi$pickupTime == "00" | nytaxi$pickupTime == "01" | nytaxi$pickupTime == "02" | nytaxi$pickupTime == "03" | nytaxi$pickupTime == "04" | nytaxi$pickupTime == "05"
check2 = nytaxi$pickupTime == "06" | nytaxi$pickupTime == "07" | nytaxi$pickupTime == "08" | nytaxi$pickupTime == "09" | nytaxi$pickupTime == "10" | nytaxi$pickupTime == "11"
check3 = nytaxi$pickupTime == "12" | nytaxi$pickupTime == "13" | nytaxi$pickupTime == "14" | nytaxi$pickupTime == "15" | nytaxi$pickupTime == "16" | nytaxi$pickupTime == "17"
nytaxi$pickupTime1 <- ifelse(check1, 1, 0)
nytaxi$pickupTime2 <- ifelse(check2, 1, 0)
nytaxi$pickupTime3 <- ifelse(check3, 1, 0)
# same for drop off time
check1 = nytaxi$dropoffTime == "00" | nytaxi$dropoffTime == "01" | nytaxi$dropoffTime == "02" | nytaxi$dropoffTime == "03" | nytaxi$dropoffTime == "04" | nytaxi$dropoffTime == "05"
check2 = nytaxi$dropoffTime == "06" | nytaxi$dropoffTime == "07" | nytaxi$dropoffTime == "08" | nytaxi$dropoffTime == "09" | nytaxi$dropoffTime == "10" | nytaxi$dropoffTime == "11"
check3 = nytaxi$dropoffTime == "12" | nytaxi$dropoffTime == "13" | nytaxi$dropoffTime == "14" | nytaxi$dropoffTime == "15" | nytaxi$dropoffTime == "16" | nytaxi$dropoffTime == "17"
nytaxi$dropoffTime1 <- ifelse(check1, 1, 0)
nytaxi$dropoffTime2 <- ifelse(check2, 1, 0)
nytaxi$dropoffTime3 <- ifelse(check3, 1, 0)
nytaxi <- subset(nytaxi, select = -c(pickupTime,dropoffTime))
y <- nytaxi$trip_duration
nytaxi <- subset(nytaxi, select = -c(trip_duration))
nytaxi <- cbind(nytaxi, y)
nytaxi <- subset(nytaxi, select = -c(id))
nytaxi$store_and_fwd_flag <- ifelse(nytaxi$store_and_fwd_flag=="N", 1, 0)
str(nytaxi)
```

    ## 'data.frame':    1458644 obs. of  16 variables:
    ##  $ vendor_id         : int  2 1 2 2 2 2 1 2 1 2 ...
    ##  $ passenger_count   : int  1 1 1 1 1 6 4 1 1 1 ...
    ##  $ pickup_longitude  : num  -74 -74 -74 -74 -74 ...
    ##  $ pickup_latitude   : num  40.8 40.7 40.8 40.7 40.8 ...
    ##  $ dropoff_longitude : num  -74 -74 -74 -74 -74 ...
    ##  $ dropoff_latitude  : num  40.8 40.7 40.7 40.7 40.8 ...
    ##  $ store_and_fwd_flag: num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ pickupWeekday     : num  1 0 1 1 0 0 1 0 1 1 ...
    ##  $ dropoffWeekday    : num  1 0 1 1 0 0 1 0 1 1 ...
    ##  $ pickupTime1       : num  0 1 0 0 0 0 0 0 0 0 ...
    ##  $ pickupTime2       : num  0 0 1 0 0 0 0 1 0 0 ...
    ##  $ pickupTime3       : num  1 0 0 0 1 0 0 0 0 0 ...
    ##  $ dropoffTime1      : num  0 1 0 0 0 0 0 0 0 0 ...
    ##  $ dropoffTime2      : num  0 0 0 0 0 0 0 1 0 0 ...
    ##  $ dropoffTime3      : num  1 0 1 0 1 0 0 0 0 0 ...
    ##  $ y                 : int  455 663 2124 429 435 443 341 1551 255 1225 ...

Use polyreg
-----------

``` r
library(polyreg)
polyE1 <- xvalPoly(nytaxi,3,2,"lm",0.8)
```

    ## Warning in predict.lm(object$fit, plm.newdata): prediction from a rank-
    ## deficient fit may be misleading

    ## Warning in predict.lm(object$fit, plm.newdata): prediction from a rank-
    ## deficient fit may be misleading

``` r
polyE1
```

    ## [1] 609.7732 607.0927 708.0398

``` r
polyE2 <- xvalPoly(nytaxi,3,2,"lm",0.8,TRUE)
polyE2
```

    ## [1]  636.9757 676.1902 588.2270

Use nnet
--------

``` r
library(nnet)
set.seed(500)
n <- nrow(nytaxi)
ntrain <- round(0.8*n)
trainidxs <- sample(1:n, ntrain, replace = FALSE)
nynn <- nnet(y~., data=nytaxi[trainidxs,],size=20,maxit=10000,decay=.001)
npred <- predict(nynn, nytaxi[-trainidxs, ])
mean(abs(nytaxi[-trainidxs,]$y - npred))
```

    ## # weights:  341
    ## initial  value 29109540357252.925781 
    ## final  value 29108245202694.953125 
    ## converged

    ## [1] 962.2623

``` r
nynn1 <- nnet(y~., data=nytaxi[trainidxs,],size=30,maxit=10000,decay=.005)
npred1 <- predict(nynn1, nytaxi[-trainidxs, ])
mean(abs(nytaxi[-trainidxs,]$y - npred1))
```

    ## # weights:  511
    ## initial  value 29109678474153.464844 
    ## final  value 29108245203792.468750 
    ## converged

    ## [1] 962.2623
