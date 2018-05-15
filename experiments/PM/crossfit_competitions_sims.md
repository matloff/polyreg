Crossfit Data
================

``` r
library(kerasformula)
library(polyreg)
library(ggplot2)
seed <- 12345
pTraining <- 0.9
```

Below find publically-available data from the Crossfit annual open, an amateur athletics competition consisting of five workouts each year. "Rx", as in "perscription", denotes the heaviest weights and most complex movements. In this analysis, we restrict the predictors to height, weight, age, region, and performance in the first of the five rounds. We also restrict analysis who did all five rounds "Rx" and reported that other data. The analysis is repeated for men and women and 2017 and 2018. In each case, `kerasformula` and `xvalPoly` are compared in terms of mean absolute error. (Note `kms` will standardize the outcome by default and so the test stats are on the same scale).

Men 2018 Competition
====================

``` r
competitions <- c("Men_Rx_2018", "Men_Rx_2017", "Women_Rx_2018", "Women_Rx_2017")

MAE_results <- matrix(nrow = 2*4*15, ncol=8)
colnames(MAE_results) <- c("MAE", "model", "competition", "openA", "openB", "seed", "N", "P")
MAE_results <- as.data.frame(MAE_results)

# which_opens
include_opens <- paste0("open", 1:5, "percentile")

r <- 1 # row of MAE_results
P <- 6

for(i in 1:4){
  
  Rx <- read.csv(paste0(competitions[i], ".csv"))   # slow step... data.table?
  colnames(Rx) <- gsub("[[:punct:]]", "", colnames(Rx))    # forgetmenot
  colnames(Rx) <- tolower(colnames(Rx))
  colnames(Rx) <- gsub("x18", "open", colnames(Rx))
  colnames(Rx) <- gsub("x17", "open", colnames(Rx))
  
  for(j in 1:4){
    for(k in (j+1):5){
      
      Rx_tmp <- dplyr::select(Rx, heightm, weightkg, age, overallpercentile)     
       Rx_tmp <- cbind(Rx_tmp, Rx[[include_opens[j]]], Rx[[include_opens[k]]])
       Rx_complete <- Rx_tmp[complete.cases(Rx_tmp), ]
       
       Rx_kms_out <- kms(overallpercentile ~ ., Rx_complete, 
                    layers = list(units = c(P, P, NA), 
                                activation = c("relu", "relu", "linear"), 
                                dropout = c(0.4, 0.3, NA), 
                                use_bias = TRUE, 
                                kernel_initializer = NULL, 
                                kernel_regularizer = "regularizer_l1_l2", 
                                bias_regularizer = "regularizer_l1_l2", 
                                activity_regularizer = "regularizer_l1_l2"),
                  seed=seed, Nepochs=5, 
                  validation_split = 0, pTraining = pTraining, verbose=0)
       Rx_complete_01 <- as.data.frame(lapply(Rx_complete, kerasformula:::zero_one))
       xval.out <- xvalPoly(Rx_complete_01, 
                          maxDeg = 3, maxInteractDeg = 2, 
                          nHoldout = floor(pTraining*nrow(Rx_complete)))
       
       MAE_results$MAE[r] <- Rx_kms_out$MAE_predictions
       MAE_results$MAE[r+1] <- min(xval.out)
       MAE_results$model[r] <- "kms"
       MAE_results$model[r+1] <- "polyreg"
       MAE_results$competition[r:(r+1)] <- competitions[i]
       MAE_results$openA[r:(r+1)] <- include_opens[[j]] # which rounds used to predict final
       MAE_results$openB[r:(r+1)] <- include_opens[[k]]
       MAE_results$seed[r:(r+1)] <- seed
       MAE_results$N[r:(r+1)] <- floor(pTraining*nrow(Rx_complete))
       MAE_results$P[r:(r+1)] <- Rx_kms_out$P
       r <- r + 2
       cat("finished simulation", r, "\n")
    }
  }  
  cat("finished", competitions[i], "\n")
}
```

The results
===========

Out-of-sample mean absolute error (MAE) for `kms` vs. `xvalPoly` (for the latter, the lowest MAE for the models corresponding to the three polynomial degrees is selected).

``` r
ggplot(MAE_results, aes(x = N, y  = MAE, color = model, pch=competition)) + geom_point()
```

![](crossfit_competitions_sims_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
MAE_results
```

               MAE   model   competition           openA           openB  seed
    1   0.09268332     kms   Men_Rx_2018 open1percentile open2percentile 12345
    2   0.09845982 polyreg   Men_Rx_2018 open1percentile open2percentile 12345
    3   0.06997361     kms   Men_Rx_2018 open1percentile open3percentile 12345
    4   0.10921087 polyreg   Men_Rx_2018 open1percentile open3percentile 12345
    5   0.06920743     kms   Men_Rx_2018 open1percentile open4percentile 12345
    6   0.11453900 polyreg   Men_Rx_2018 open1percentile open4percentile 12345
    7   0.05835935     kms   Men_Rx_2018 open1percentile open5percentile 12345
    8   0.07418255 polyreg   Men_Rx_2018 open1percentile open5percentile 12345
    9   0.07112396     kms   Men_Rx_2018 open2percentile open3percentile 12345
    10  0.10518538 polyreg   Men_Rx_2018 open2percentile open3percentile 12345
    11  0.06709697     kms   Men_Rx_2018 open2percentile open4percentile 12345
    12  0.11110825 polyreg   Men_Rx_2018 open2percentile open4percentile 12345
    13  0.06097750     kms   Men_Rx_2018 open2percentile open5percentile 12345
    14  0.07481714 polyreg   Men_Rx_2018 open2percentile open5percentile 12345
    15  0.07883466     kms   Men_Rx_2018 open3percentile open4percentile 12345
    16  0.11347619 polyreg   Men_Rx_2018 open3percentile open4percentile 12345
    17  0.06383485     kms   Men_Rx_2018 open3percentile open5percentile 12345
    18  0.07209419 polyreg   Men_Rx_2018 open3percentile open5percentile 12345
    19  0.05984704     kms   Men_Rx_2018 open4percentile open5percentile 12345
    20  0.07275715 polyreg   Men_Rx_2018 open4percentile open5percentile 12345
    21  0.09440996     kms   Men_Rx_2017 open1percentile open2percentile 12345
    22  0.10255934 polyreg   Men_Rx_2017 open1percentile open2percentile 12345
    23  0.07716992     kms   Men_Rx_2017 open1percentile open3percentile 12345
    24  0.11094495 polyreg   Men_Rx_2017 open1percentile open3percentile 12345
    25  0.07965579     kms   Men_Rx_2017 open1percentile open4percentile 12345
    26  0.09707374 polyreg   Men_Rx_2017 open1percentile open4percentile 12345
    27  0.07561122     kms   Men_Rx_2017 open1percentile open5percentile 12345
    28  0.10851431 polyreg   Men_Rx_2017 open1percentile open5percentile 12345
    29  0.09069661     kms   Men_Rx_2017 open2percentile open3percentile 12345
    30  0.11499086 polyreg   Men_Rx_2017 open2percentile open3percentile 12345
    31  0.07487309     kms   Men_Rx_2017 open2percentile open4percentile 12345
    32  0.09751620 polyreg   Men_Rx_2017 open2percentile open4percentile 12345
    33  0.07410852     kms   Men_Rx_2017 open2percentile open5percentile 12345
    34  0.10699325 polyreg   Men_Rx_2017 open2percentile open5percentile 12345
    35  0.07495597     kms   Men_Rx_2017 open3percentile open4percentile 12345
    36  0.09844750 polyreg   Men_Rx_2017 open3percentile open4percentile 12345
    37  0.07925495     kms   Men_Rx_2017 open3percentile open5percentile 12345
    38  0.10903266 polyreg   Men_Rx_2017 open3percentile open5percentile 12345
    39  0.19622631     kms   Men_Rx_2017 open4percentile open5percentile 12345
    40  0.10795536 polyreg   Men_Rx_2017 open4percentile open5percentile 12345
    41  0.09226367     kms Women_Rx_2018 open1percentile open2percentile 12345
    42  0.10478415 polyreg Women_Rx_2018 open1percentile open2percentile 12345
    43  0.08337254     kms Women_Rx_2018 open1percentile open3percentile 12345
    44  0.13284574 polyreg Women_Rx_2018 open1percentile open3percentile 12345
    45  0.08783822     kms Women_Rx_2018 open1percentile open4percentile 12345
    46  0.10346825 polyreg Women_Rx_2018 open1percentile open4percentile 12345
    47  0.06626291     kms Women_Rx_2018 open1percentile open5percentile 12345
    48  0.08821380 polyreg Women_Rx_2018 open1percentile open5percentile 12345
    49  0.08575813     kms Women_Rx_2018 open2percentile open3percentile 12345
    50  0.13622916 polyreg Women_Rx_2018 open2percentile open3percentile 12345
    51  0.07231336     kms Women_Rx_2018 open2percentile open4percentile 12345
    52  0.10293313 polyreg Women_Rx_2018 open2percentile open4percentile 12345
    53  0.06358446     kms Women_Rx_2018 open2percentile open5percentile 12345
    54  0.08954019 polyreg Women_Rx_2018 open2percentile open5percentile 12345
    55  0.09453945     kms Women_Rx_2018 open3percentile open4percentile 12345
    56  0.11282402 polyreg Women_Rx_2018 open3percentile open4percentile 12345
    57  0.06915858     kms Women_Rx_2018 open3percentile open5percentile 12345
    58  0.09492999 polyreg Women_Rx_2018 open3percentile open5percentile 12345
    59  0.07955127     kms Women_Rx_2018 open4percentile open5percentile 12345
    60  0.08170694 polyreg Women_Rx_2018 open4percentile open5percentile 12345
    61  0.11619483     kms Women_Rx_2017 open1percentile open2percentile 12345
    62  0.13836345 polyreg Women_Rx_2017 open1percentile open2percentile 12345
    63  0.08749113     kms Women_Rx_2017 open1percentile open3percentile 12345
    64  0.14863574 polyreg Women_Rx_2017 open1percentile open3percentile 12345
    65  0.11193769     kms Women_Rx_2017 open1percentile open4percentile 12345
    66  0.11604917 polyreg Women_Rx_2017 open1percentile open4percentile 12345
    67  0.08438283     kms Women_Rx_2017 open1percentile open5percentile 12345
    68  0.12949881 polyreg Women_Rx_2017 open1percentile open5percentile 12345
    69  0.16934012     kms Women_Rx_2017 open2percentile open3percentile 12345
    70  0.13458951 polyreg Women_Rx_2017 open2percentile open3percentile 12345
    71  0.10377900     kms Women_Rx_2017 open2percentile open4percentile 12345
    72  0.12155723 polyreg Women_Rx_2017 open2percentile open4percentile 12345
    73  0.10276844     kms Women_Rx_2017 open2percentile open5percentile 12345
    74  0.11885093 polyreg Women_Rx_2017 open2percentile open5percentile 12345
    75  0.11130030     kms Women_Rx_2017 open3percentile open4percentile 12345
    76  0.10909286 polyreg Women_Rx_2017 open3percentile open4percentile 12345
    77  0.11647202     kms Women_Rx_2017 open3percentile open5percentile 12345
    78  0.12075639 polyreg Women_Rx_2017 open3percentile open5percentile 12345
    79  0.08927585     kms Women_Rx_2017 open4percentile open5percentile 12345
    80  0.11732735 polyreg Women_Rx_2017 open4percentile open5percentile 12345
    81          NA    <NA>          <NA>            <NA>            <NA>    NA
    82          NA    <NA>          <NA>            <NA>            <NA>    NA
    83          NA    <NA>          <NA>            <NA>            <NA>    NA
    84          NA    <NA>          <NA>            <NA>            <NA>    NA
    85          NA    <NA>          <NA>            <NA>            <NA>    NA
    86          NA    <NA>          <NA>            <NA>            <NA>    NA
    87          NA    <NA>          <NA>            <NA>            <NA>    NA
    88          NA    <NA>          <NA>            <NA>            <NA>    NA
    89          NA    <NA>          <NA>            <NA>            <NA>    NA
    90          NA    <NA>          <NA>            <NA>            <NA>    NA
    91          NA    <NA>          <NA>            <NA>            <NA>    NA
    92          NA    <NA>          <NA>            <NA>            <NA>    NA
    93          NA    <NA>          <NA>            <NA>            <NA>    NA
    94          NA    <NA>          <NA>            <NA>            <NA>    NA
    95          NA    <NA>          <NA>            <NA>            <NA>    NA
    96          NA    <NA>          <NA>            <NA>            <NA>    NA
    97          NA    <NA>          <NA>            <NA>            <NA>    NA
    98          NA    <NA>          <NA>            <NA>            <NA>    NA
    99          NA    <NA>          <NA>            <NA>            <NA>    NA
    100         NA    <NA>          <NA>            <NA>            <NA>    NA
    101         NA    <NA>          <NA>            <NA>            <NA>    NA
    102         NA    <NA>          <NA>            <NA>            <NA>    NA
    103         NA    <NA>          <NA>            <NA>            <NA>    NA
    104         NA    <NA>          <NA>            <NA>            <NA>    NA
    105         NA    <NA>          <NA>            <NA>            <NA>    NA
    106         NA    <NA>          <NA>            <NA>            <NA>    NA
    107         NA    <NA>          <NA>            <NA>            <NA>    NA
    108         NA    <NA>          <NA>            <NA>            <NA>    NA
    109         NA    <NA>          <NA>            <NA>            <NA>    NA
    110         NA    <NA>          <NA>            <NA>            <NA>    NA
    111         NA    <NA>          <NA>            <NA>            <NA>    NA
    112         NA    <NA>          <NA>            <NA>            <NA>    NA
    113         NA    <NA>          <NA>            <NA>            <NA>    NA
    114         NA    <NA>          <NA>            <NA>            <NA>    NA
    115         NA    <NA>          <NA>            <NA>            <NA>    NA
    116         NA    <NA>          <NA>            <NA>            <NA>    NA
    117         NA    <NA>          <NA>            <NA>            <NA>    NA
    118         NA    <NA>          <NA>            <NA>            <NA>    NA
    119         NA    <NA>          <NA>            <NA>            <NA>    NA
    120         NA    <NA>          <NA>            <NA>            <NA>    NA
             N  P
    1   101460  5
    2   101460  5
    3    86850  5
    4    86850  5
    5    87443  5
    6    87443  5
    7    89497  5
    8    89497  5
    9    87278  5
    10   87278  5
    11   88128  5
    12   88128  5
    13   90699  5
    14   90699  5
    15   81505  5
    16   81505  5
    17   80512  5
    18   80512  5
    19   84092  5
    20   84092  5
    21   32724  5
    22   32724  5
    23   29630  5
    24   29630  5
    25   30063  5
    26   30063  5
    27   26327  5
    28   26327  5
    29   29136  5
    30   29136  5
    31   29061  5
    32   29061  5
    33   25708  5
    34   25708  5
    35   27620  5
    36   27620  5
    37   24727  5
    38   24727  5
    39   25623  5
    40   25623  5
    41   42352  5
    42   42352  5
    43   35620  5
    44   35620  5
    45   36221  5
    46   36221  5
    47   36585  5
    48   36585  5
    49   37214  5
    50   37214  5
    51   37767  5
    52   37767  5
    53   37851  5
    54   37851  5
    55   34420  5
    56   34420  5
    57   33838  5
    58   33838  5
    59   35382  5
    60   35382  5
    61    9558  5
    62    9558  5
    63    7354  5
    64    7354  5
    65    9373  5
    66    9373  5
    67    8597  5
    68    8597  5
    69    7222  5
    70    7222  5
    71    8725  5
    72    8725  5
    73    8069  5
    74    8069  5
    75    6966  5
    76    6966  5
    77    6574  5
    78    6574  5
    79    8380  5
    80    8380  5
    81      NA NA
    82      NA NA
    83      NA NA
    84      NA NA
    85      NA NA
    86      NA NA
    87      NA NA
    88      NA NA
    89      NA NA
    90      NA NA
    91      NA NA
    92      NA NA
    93      NA NA
    94      NA NA
    95      NA NA
    96      NA NA
    97      NA NA
    98      NA NA
    99      NA NA
    100     NA NA
    101     NA NA
    102     NA NA
    103     NA NA
    104     NA NA
    105     NA NA
    106     NA NA
    107     NA NA
    108     NA NA
    109     NA NA
    110     NA NA
    111     NA NA
    112     NA NA
    113     NA NA
    114     NA NA
    115     NA NA
    116     NA NA
    117     NA NA
    118     NA NA
    119     NA NA
    120     NA NA
