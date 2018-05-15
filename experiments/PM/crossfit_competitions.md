Crossfit Data
================

``` r
library(kerasformula)
library(polyreg)
```

Below find publically-available data from the Crossfit annual open, an amateur athletics competition consisting of five workouts each year. "Rx", as in "perscription", denotes the heaviest weights and most complex movements. In this analysis, we restrict the predictors to height, weight, age, region, and performance in the first of the five rounds. We also restrict analysis who did all five rounds "Rx" and reported that other data. The analysis is repeated for men and women and 2017 and 2018. In each case, `kerasformula` and `xvalPoly` are compared in terms of mean absolute error. (Note `kms` will standardize the outcome by default and so the test stats are on the same scale).

Men 2018 Competition
====================

``` r
Rx <- read.csv("Men_Rx_2018.csv")
```

``` r
colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
colnames(Rx) <- tolower(colnames(Rx))
colnames(Rx) <- gsub("x18", "open", colnames(Rx))

Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, open1percentile, overallpercentile)
Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
dim(Rx_complete)
```

    [1] 120216      5

``` r
head(Rx_complete)
```

      heightm weightkg age open1percentile overallpercentile
    1    1.70       86  28         99.9982          100.0000
    2    1.80       92  26         99.9941           99.9995
    3    1.79       86  25         99.9750           99.9989
    4    1.76       88  24         99.9899           99.9984
    5    1.80       88  28         99.9982           99.9979
    6    1.80       91  22         99.9156           99.9973

``` r
P <- ncol(model.matrix(overallpercentile ~ ., Rx_complete))

Rx_kms_out <- kms(overallpercentile ~ ., Rx_complete, 
                  layers = list(units = c(P, P, NA), 
                                activation = c("relu", "relu", "linear"), 
                                dropout = c(0.4, 0.3, NA), 
                                use_bias = TRUE, 
                                kernel_initializer = NULL, 
                                kernel_regularizer = "regularizer_l1_l2", 
                                bias_regularizer = "regularizer_l1_l2", 
                                activity_regularizer = "regularizer_l1_l2"),
                  seed=77777, Nepochs=5, 
                  validation_split = 0, pTraining = 0.9)
```

    ___________________________________________________________________________
    Layer (type)                     Output Shape                  Param #     
    ===========================================================================
    dense_1 (Dense)                  (None, 5)                     25          
    ___________________________________________________________________________
    dropout_1 (Dropout)              (None, 5)                     0           
    ___________________________________________________________________________
    dense_2 (Dense)                  (None, 5)                     30          
    ___________________________________________________________________________
    dropout_2 (Dropout)              (None, 5)                     0           
    ___________________________________________________________________________
    dense_3 (Dense)                  (None, 1)                     6           
    ===========================================================================
    Total params: 61
    Trainable params: 61
    Non-trainable params: 0
    ___________________________________________________________________________

``` r
Rx_complete_z <- as.data.frame(lapply(Rx_complete, kerasformula:::z))
xval.out <- xvalPoly(Rx_complete_z, maxDeg = 3, maxInteractDeg = 2)
```

    getPoly time in xvalPoly:  0.806 0.155 0.965 0 0 
    lm() time:  0.077 0.01 0.088 0 0 
    accuracy:  0.4026333 
    lm() time:  0.129 0.026 0.156 0 0 
    accuracy:  0.4012934 
    lm() time:  0.487 0.066 0.556 0 0 
    accuracy:  0.3996115 

``` r
xval.out
```

    [1] 0.4026333 0.4012934 0.3996115

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.1134672

Men 2017 Competition
====================

``` r
Rx <- read.csv("Men_Rx_2017.csv")
```

``` r
colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
colnames(Rx) <- tolower(colnames(Rx))
colnames(Rx) <- gsub("x17", "open", colnames(Rx))

Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, open1percentile, overallpercentile)
Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
dim(Rx_complete)
```

    [1] 41896     5

``` r
Rx_kms_out <- kms(overallpercentile ~ ., 
                  Rx_complete, seed=777, Nepochs=5, 
                  validation_split = 0, pTraining = 0.9)
```

    ___________________________________________________________________________
    Layer (type)                     Output Shape                  Param #     
    ===========================================================================
    dense_1 (Dense)                  (None, 256)                   1280        
    ___________________________________________________________________________
    dropout_1 (Dropout)              (None, 256)                   0           
    ___________________________________________________________________________
    dense_2 (Dense)                  (None, 128)                   32896       
    ___________________________________________________________________________
    dropout_2 (Dropout)              (None, 128)                   0           
    ___________________________________________________________________________
    dense_3 (Dense)                  (None, 1)                     129         
    ===========================================================================
    Total params: 34,305
    Trainable params: 34,305
    Non-trainable params: 0
    ___________________________________________________________________________

``` r
Rx_complete_z <- as.data.frame(lapply(Rx_complete, kerasformula:::z))
xval.out <- xvalPoly(Rx_complete_z, maxDeg = 3, maxInteractDeg = 1)
```

    getPoly time in xvalPoly:  0.218 0.056 0.275 0 0 
    lm() time:  0.025 0.004 0.029 0 0 
    accuracy:  0.4981946 
    lm() time:  0.039 0.008 0.047 0 0 
    accuracy:  0.4973153 
    lm() time:  0.09 0.018 0.109 0 0 
    accuracy:  0.4944461 

``` r
xval.out
```

    [1] 0.4981946 0.4973153 0.4944461

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.1349026

Women 2018 Competition
======================

``` r
Rx <- read.csv("Women_Rx_2018.csv")
```

``` r
colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
colnames(Rx) <- tolower(colnames(Rx))
colnames(Rx) <- gsub("x18", "open", colnames(Rx))

Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, open1percentile, overallpercentile)
Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
dim(Rx_complete)
```

    [1] 49853     5

``` r
Rx_kms_out <- kms(overallpercentile ~ ., 
                  Rx_complete, seed=777, Nepochs=5, 
                  validation_split = 0, pTraining = 0.9)
```

    ___________________________________________________________________________
    Layer (type)                     Output Shape                  Param #     
    ===========================================================================
    dense_1 (Dense)                  (None, 256)                   1280        
    ___________________________________________________________________________
    dropout_1 (Dropout)              (None, 256)                   0           
    ___________________________________________________________________________
    dense_2 (Dense)                  (None, 128)                   32896       
    ___________________________________________________________________________
    dropout_2 (Dropout)              (None, 128)                   0           
    ___________________________________________________________________________
    dense_3 (Dense)                  (None, 1)                     129         
    ===========================================================================
    Total params: 34,305
    Trainable params: 34,305
    Non-trainable params: 0
    ___________________________________________________________________________

``` r
Rx_complete_z <- as.data.frame(lapply(Rx_complete, kerasformula:::z))
```

``` r
xval.out <- xvalPoly(Rx_complete_z, maxDeg = 3, maxInteractDeg = 1)
```

    getPoly time in xvalPoly:  0.213 0.059 0.274 0 0 
    lm() time:  0.03 0.004 0.035 0 0 
    accuracy:  0.4366127 
    lm() time:  0.047 0.01 0.057 0 0 
    accuracy:  0.4347103 
    lm() time:  0.113 0.025 0.138 0 0 
    accuracy:  0.449277 

``` r
xval.out
```

    [1] 0.4366127 0.4347103 0.4492770

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.110857

Women 2017
==========

``` r
Rx <- read.csv("Women_Rx_2017.csv")
```

``` r
colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
colnames(Rx) <- tolower(colnames(Rx))
colnames(Rx) <- gsub("x17", "open", colnames(Rx))

Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, open1percentile, overallpercentile)
Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
dim(Rx_complete)
```

    [1] 13104     5

``` r
Rx_kms_out <- kms(overallpercentile ~ ., 
                  Rx_complete, seed=777, Nepochs=5, 
                  validation_split = 0, pTraining = 0.9)
```

    ___________________________________________________________________________
    Layer (type)                     Output Shape                  Param #     
    ===========================================================================
    dense_1 (Dense)                  (None, 256)                   1280        
    ___________________________________________________________________________
    dropout_1 (Dropout)              (None, 256)                   0           
    ___________________________________________________________________________
    dense_2 (Dense)                  (None, 128)                   32896       
    ___________________________________________________________________________
    dropout_2 (Dropout)              (None, 128)                   0           
    ___________________________________________________________________________
    dense_3 (Dense)                  (None, 1)                     129         
    ===========================================================================
    Total params: 34,305
    Trainable params: 34,305
    Non-trainable params: 0
    ___________________________________________________________________________

``` r
Rx_complete_z <- as.data.frame(lapply(Rx_complete, kerasformula:::z))
xval.out <- xvalPoly(Rx_complete_z, maxDeg = 3, maxInteractDeg = 1)
```

    getPoly time in xvalPoly:  0.048 0.008 0.058 0 0 
    lm() time:  0.007 0 0.007 0 0 
    accuracy:  0.5351013 
    lm() time:  0.012 0 0.012 0 0 
    accuracy:  0.5598741 
    lm() time:  0.03 0.006 0.036 0 0 
    accuracy:  0.8890716 

``` r
xval.out
```

    [1] 0.5351013 0.5598741 0.8890716

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.1536684
