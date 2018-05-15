Crossfit Data
================

``` r
knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment="")
```

``` r
library(kerasformula)
```

Below find publically-available data from the Crossfit annual open, an amateur athletics competition consisting of five workouts each year. "Rx", as in "perscription", denotes the heaviest weights and most complex movements. In this analysis, we restrict the predictors to height, weight, age, region, and performance in the first of the five rounds. We also restrict analysis who did all five rounds "Rx" and reported that other data. The analysis is repeated for men and women and 2017 and 2018.

Men 2018 Competition
====================

``` r
Rx <- read.csv("Men_Rx_2018.csv")
colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
colnames(Rx) <- tolower(colnames(Rx))
colnames(Rx) <- gsub("x18", "open", colnames(Rx))

Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, regionname, open1percentile, overallpercentile)
Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
dim(Rx_complete)
```

    [1] 120216      6

``` r
head(Rx_complete)
```

      heightm weightkg age   regionname open1percentile overallpercentile
    1    1.70       86  28 Central East         99.9982          100.0000
    2    1.80       92  26  Canada East         99.9941           99.9995
    3    1.79       86  25 Europe South         99.9750           99.9989
    4    1.76       88  24  Canada East         99.9899           99.9984
    5    1.80       88  28  Canada East         99.9982           99.9979
    6    1.80       91  22  Canada East         99.9156           99.9973

``` r
Rx_kms_out <- kms(overallpercentile ~ ., 
                  Rx_complete, seed=777, Nepochs=5, 
                  validation_split = 0, pTraining = 0.9)
```

    ___________________________________________________________________________
    Layer (type)                     Output Shape                  Param #     
    ===========================================================================
    dense_1 (Dense)                  (None, 256)                   5632        
    ___________________________________________________________________________
    dropout_1 (Dropout)              (None, 256)                   0           
    ___________________________________________________________________________
    dense_2 (Dense)                  (None, 128)                   32896       
    ___________________________________________________________________________
    dropout_2 (Dropout)              (None, 128)                   0           
    ___________________________________________________________________________
    dense_3 (Dense)                  (None, 1)                     129         
    ===========================================================================
    Total params: 38,657
    Trainable params: 38,657
    Non-trainable params: 0
    ___________________________________________________________________________

``` r
Rx_kms_out$R2_predictions
```

              [,1]
    [1,] 0.6743121

``` r
Rx_kms_out$cor_kendals^2   # R^2 kendalls
```

              [,1]
    [1,] 0.4303534

``` r
Rx_kms_out$MSE_predictions
```

    [1] 0.0249767

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.1142787

Men 2017 Competition
====================

``` r
Rx <- read.csv("Men_Rx_2017.csv")
colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
colnames(Rx) <- tolower(colnames(Rx))
colnames(Rx) <- gsub("x17", "open", colnames(Rx))

Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, regionname, open1percentile, overallpercentile)
Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
dim(Rx_complete)
```

    [1] 41896     6

``` r
Rx_kms_out <- kms(overallpercentile ~ ., 
                  Rx_complete, seed=777, Nepochs=5, 
                  validation_split = 0, pTraining = 0.9)
```

    ___________________________________________________________________________
    Layer (type)                     Output Shape                  Param #     
    ===========================================================================
    dense_1 (Dense)                  (None, 256)                   5376        
    ___________________________________________________________________________
    dropout_1 (Dropout)              (None, 256)                   0           
    ___________________________________________________________________________
    dense_2 (Dense)                  (None, 128)                   32896       
    ___________________________________________________________________________
    dropout_2 (Dropout)              (None, 128)                   0           
    ___________________________________________________________________________
    dense_3 (Dense)                  (None, 1)                     129         
    ===========================================================================
    Total params: 38,401
    Trainable params: 38,401
    Non-trainable params: 0
    ___________________________________________________________________________

``` r
Rx_kms_out$R2_predictions
```

              [,1]
    [1,] 0.5659958

``` r
Rx_kms_out$cor_kendals^2   # R^2 kendalls
```

              [,1]
    [1,] 0.3267149

``` r
Rx_kms_out$MSE_predictions
```

    [1] 0.0311001

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.1357579

Women 2018 Competition
======================

``` r
Rx <- read.csv("Women_Rx_2018.csv")
colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
colnames(Rx) <- tolower(colnames(Rx))
colnames(Rx) <- gsub("x18", "open", colnames(Rx))

Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, regionname, open1percentile, overallpercentile)
Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
dim(Rx_complete)
```

    [1] 49853     6

``` r
Rx_kms_out <- kms(overallpercentile ~ ., 
                  Rx_complete, seed=777, Nepochs=5, 
                  validation_split = 0, pTraining = 0.9)
```

    ___________________________________________________________________________
    Layer (type)                     Output Shape                  Param #     
    ===========================================================================
    dense_1 (Dense)                  (None, 256)                   5632        
    ___________________________________________________________________________
    dropout_1 (Dropout)              (None, 256)                   0           
    ___________________________________________________________________________
    dense_2 (Dense)                  (None, 128)                   32896       
    ___________________________________________________________________________
    dropout_2 (Dropout)              (None, 128)                   0           
    ___________________________________________________________________________
    dense_3 (Dense)                  (None, 1)                     129         
    ===========================================================================
    Total params: 38,657
    Trainable params: 38,657
    Non-trainable params: 0
    ___________________________________________________________________________

``` r
Rx_kms_out$R2_predictions
```

              [,1]
    [1,] 0.5998907

``` r
Rx_kms_out$cor_kendals^2   # R^2 kendalls
```

              [,1]
    [1,] 0.4107823

``` r
Rx_kms_out$MSE_predictions
```

    [1] 0.02570468

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.1118306

Women 2017
==========

``` r
Rx <- read.csv("Women_Rx_2017.csv")
colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
colnames(Rx) <- tolower(colnames(Rx))
colnames(Rx) <- gsub("x17", "open", colnames(Rx))

Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, regionname, open1percentile, overallpercentile)
Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
dim(Rx_complete)
```

    [1] 13104     6

``` r
Rx_kms_out <- kms(overallpercentile ~ ., 
                  Rx_complete, seed=777, Nepochs=5, 
                  validation_split = 0, pTraining = 0.9)
```

    ___________________________________________________________________________
    Layer (type)                     Output Shape                  Param #     
    ===========================================================================
    dense_1 (Dense)                  (None, 256)                   5376        
    ___________________________________________________________________________
    dropout_1 (Dropout)              (None, 256)                   0           
    ___________________________________________________________________________
    dense_2 (Dense)                  (None, 128)                   32896       
    ___________________________________________________________________________
    dropout_2 (Dropout)              (None, 128)                   0           
    ___________________________________________________________________________
    dense_3 (Dense)                  (None, 1)                     129         
    ===========================================================================
    Total params: 38,401
    Trainable params: 38,401
    Non-trainable params: 0
    ___________________________________________________________________________

``` r
Rx_kms_out$R2_predictions
```

             [,1]
    [1,] 0.480236

``` r
Rx_kms_out$cor_kendals^2   # R^2 kendalls
```

              [,1]
    [1,] 0.2789829

``` r
Rx_kms_out$MSE_predictions
```

    [1] 0.03639114

``` r
Rx_kms_out$MAE_predictions
```

    [1] 0.1474672
