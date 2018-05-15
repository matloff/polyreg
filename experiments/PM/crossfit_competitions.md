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
MAE_results <- matrix(nrow=2, ncol=4, dimnames=list(c("kms", "xvalPoly"), c("m2018", "m2017", "w2018", "w2017")))

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

    getPoly time in xvalPoly:  0.826 0.15 1 0 0 
    lm() time:  0.082 0.008 0.091 0 0 
    accuracy:  0.4026333 
    lm() time:  0.124 0.025 0.15 0 0 
    accuracy:  0.4012934 
    lm() time:  0.484 0.065 0.556 0 0 
    accuracy:  0.3996115 

``` r
xval.out
```

    [1] 0.4026333 0.4012934 0.3996115

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.1131168

``` r
MAE_results[1,1] <- Rx_kms_out$MAE_predictions
MAE_results[2,1] <- min(xval.out)
```

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
xval.out <- xvalPoly(Rx_complete_z, maxDeg = 3, maxInteractDeg = 1)
```

    getPoly time in xvalPoly:  0.231 0.061 0.298 0 0 
    lm() time:  0.027 0.003 0.032 0 0 
    accuracy:  0.4962781 
    lm() time:  0.042 0.009 0.052 0 0 
    accuracy:  0.4948952 
    lm() time:  0.097 0.02 0.117 0 0 
    accuracy:  0.5094317 

``` r
xval.out
```

    [1] 0.4962781 0.4948952 0.5094317

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.148042

``` r
MAE_results[1,2] <- Rx_kms_out$MAE_predictions
MAE_results[2,2] <- min(xval.out)
```

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
```

``` r
xval.out <- xvalPoly(Rx_complete_z, maxDeg = 3, maxInteractDeg = 1)
```

    getPoly time in xvalPoly:  0.225 0.063 0.3 0 0 
    lm() time:  0.032 0.004 0.036 0 0 
    accuracy:  0.436813 
    lm() time:  0.045 0.009 0.054 0 0 
    accuracy:  0.4345446 
    lm() time:  0.126 0.027 0.156 0 0 
    accuracy:  0.4339059 

``` r
xval.out
```

    [1] 0.4368130 0.4345446 0.4339059

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.122345

``` r
MAE_results[1,3] <- Rx_kms_out$MAE_predictions
MAE_results[2,3] <- min(xval.out)
```

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
xval.out <- xvalPoly(Rx_complete_z, maxDeg = 3, maxInteractDeg = 1)
```

    getPoly time in xvalPoly:  0.054 0.012 0.066 0 0 
    lm() time:  0.007 0.001 0.007 0 0 
    accuracy:  0.5422355 
    lm() time:  0.013 0.001 0.013 0 0 
    accuracy:  0.5667926 
    lm() time:  0.029 0.005 0.034 0 0 
    accuracy:  1.611625 

``` r
xval.out
```

    [1] 0.5422355 0.5667926 1.6116251

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.166333

``` r
MAE_results[1,4] <- Rx_kms_out$MAE_predictions
MAE_results[2,4] <- min(xval.out)
```

The results
===========

Out-of-sample mean absolute error (MAE) for `kms` vs. `xvalPoly` (for the latter, the lowest MAE for the models corresponding to the three polynomial degrees is selected).

``` r
MAE_results
```

                 m2018     m2017     w2018     w2017
    kms      0.1131168 0.1480420 0.1223450 0.1663330
    xvalPoly 0.3996115 0.4948952 0.4339059 0.5422355
