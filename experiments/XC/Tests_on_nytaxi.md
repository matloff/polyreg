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
    
``` r
t <- system.time(e3 <- xvalPoly(nytaxi,4,2,"lm",0.8,TRUE))
e3
```
    ## [1] 636.9757 676.1902 588.2270 543.6791

``` r
t
```
    ## user   system  elapsed
    ## 627.854  268.115 1042.937

NN
--------

``` r
xvalNnet(nytaxi, 10, FALSE)
```

    ## # weights:  171
    ## initial  value 41177886297900.554688 
    ## final  value 41176425758745.000000 
    ## converged
    ## [1] 1024.766

``` r
xvalNnet(nytaxi, 10, TRUE)
```

    ## # weights:  171
    ## initial  value 41179492639648.429688 
    ## final  value 39846815903483.242188 
    ## converged
    ## [1] 697.8903

``` r
xvalNnet(nytaxi, 20, TRUE)
```
    ## # weights:  341
    ## initial  value 41181754445990.875000 
    ## final  value 39846815903465.554688 
    ## converged
    ## [1] 697.8905
``` r
xvalNnet(nytaxi, 10, TRUE, scaleXMat = TRUE)
```
    ## # weights:  171
    ## initial  value 41178685261225.015625 
    ## iter  10 value 39656113416532.195312
    ## iter  20 value 39494341317865.414062
    ## iter  30 value 39447970266546.703125
    ## iter  40 value 39421261147416.734375
    ## iter  50 value 39396487665756.843750
    ## iter  60 value 39387215509856.609375
    ## iter  70 value 39383057160969.226562
    ## iter  80 value 39381409673662.734375
    ## iter  90 value 39373094490773.351562
    ## iter 100 value 39351902000363.765625
    ## final  value 39351902000363.765625 
    ## stopped after 100 iterations
    ## [1] 572.457
    
``` r
xvalDnet(nytaxi[,1:(ncol(nytaxi)-1)], nytaxi[,ncol(nytaxi)], hidden=c(5,5))
```
    ####loss on step 10000 is : 414133.440000
    ####loss on step 20000 is : 329143.740000
    ####loss on step 30000 is : 571835.260000
    ####loss on step 40000 is : 37790387.760000
    ####loss on step 10000 is : 572719.220000
    ####loss on step 20000 is : 660127.845000
    ####loss on step 30000 is : 669182.675000
    ####loss on step 40000 is : 411117.860000
    ## [1] 1024.766
    
``` r
xvalDnet(nytaxi[,1:(ncol(nytaxi)-1)], nytaxi[,ncol(nytaxi)], hidden=c(2,2,2))
```
    ####loss on step 10000 is : 702229.290000
    ####loss on step 20000 is : 583618.890000
    ####loss on step 30000 is : 582522.900000
    ####loss on step 40000 is : 573078.175000
    ####loss on step 10000 is : 390625.595000
    ####loss on step 20000 is : 669002.045000
    ####loss on step 30000 is : 471985.960000
    ####loss on step 40000 is : 523480.440000
    ####loss on step 10000 is : 503311.300000
    ####loss on step 20000 is : 482499.385000
    ####loss on step 30000 is : 36635487.370000
    ####loss on step 40000 is : 467624.255000
    ## [1] 1024.766
