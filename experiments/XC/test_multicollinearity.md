Test Multicollinearity of Each layer output in Neural Networks
================

First model
-----------

``` r
library(car)
library(keras)
mnist <- dataset_mnist()
x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y

# reshape
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# rescale
x_train <- x_train / 255
x_test <- x_test / 255

y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)

## Defining the Model
model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 10, activation = 'relu', input_shape = c(784)) %>% 
  layer_dense(units = 10, activation = 'relu') %>%
  layer_dense(units = 10, activation = 'softmax')
summary(model)
```

    ## ___________________________________________________________________________
    ## Layer (type)                     Output Shape                  Param #     
    ## ===========================================================================
    ## dense_1 (Dense)                  (None, 10)                    7850        
    ## ___________________________________________________________________________
    ## dense_2 (Dense)                  (None, 10)                    110         
    ## ___________________________________________________________________________
    ## dense_3 (Dense)                  (None, 10)                    110         
    ## ===========================================================================
    ## Total params: 8,070
    ## Trainable params: 8,070
    ## Non-trainable params: 0
    ## ___________________________________________________________________________

``` r
## compile the model with appropriate loss function, optimizer, and metrics
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)
## 
history <- model %>% fit(
  x_train, y_train, 
  epochs = 30, batch_size = 128, 
  validation_split = 0.2
)


y <- apply(y_test, 1, which.max)
for (i in 1:3) {
  layer_model <- keras_model(inputs = model$input,
                             outputs = get_layer(model, index = i)$output)
  output <- predict(layer_model, x_test)
  df <- as.data.frame(cbind(y, output))
  vars <- paste(colnames(df)[-1], collapse = " + ")
  ff <- as.formula(paste("y", vars, sep = " ~ "))
  mod <- lm(ff, data=df)
  # check for perfect multicollinearity (exact linear relationship)
  ld.vars <- attributes(alias(mod)$Complete)$dimnames[[1]]
  if (!is.null(ld.vars)) {
    print("Perfect Multicollinearity occurs! Following variables are omitted:")
    print(ld.vars)
    formula.new <- as.formula(
      paste(
        paste(deparse(ff), collapse=""), 
        paste(ld.vars, collapse="-"),
        sep="-"
      )
    )
    mod <- lm(formula.new,data = df)
  }
  print("==============vif===============")
  vifs <- vif(mod)
  print(vifs)
  print("Percentage of VIF of coefficients that > 10:")
  print(mean(ifelse(vifs > 10, 1, 0)))
  print("Average VIF for all coefficients:")
  print(mean(vifs))
  print("================================")
}
```

    ## [1] "==============vif==============="
    ##       V2       V3       V4       V5       V6       V7       V8       V9 
    ## 1.446667 1.693511 2.486690 2.426699 1.799247 1.834349 2.518777 1.642510 
    ##      V10      V11 
    ## 2.156092 2.281227 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 0
    ## [1] "Average VIF for all coefficients:"
    ## [1] 2.028577
    ## [1] "================================"
    ## [1] "==============vif==============="
    ##       V2       V3       V4       V5       V6       V7       V8       V9 
    ## 1.810463 3.181093 3.951218 1.722111 2.045886 1.796799 2.865722 4.229225 
    ##      V10      V11 
    ## 2.280576 2.664334 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 0
    ## [1] "Average VIF for all coefficients:"
    ## [1] 2.654743
    ## [1] "================================"
    ## [1] "==============vif==============="
    ##           V2           V3           V4           V5           V6 
    ## 2.189511e+13 2.468096e+13 2.011231e+13 2.118642e+13 2.165339e+13 
    ##           V7           V8           V9          V10          V11 
    ## 1.780298e+13 2.071155e+13 2.094496e+13 1.887499e+13 2.055006e+13 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 1
    ## [1] "Average VIF for all coefficients:"
    ## [1] 2.084127e+13
    ## [1] "================================"

Second model
------------

``` r
## Defining the Model # 256, 128, 10
model2 <- keras_model_sequential() 
model2 %>% 
  layer_dense(units = 100, activation = 'relu', input_shape = c(784)) %>% 
  layer_dense(units = 50, activation = 'relu') %>%
  layer_dense(units = 10, activation = 'softmax')
summary(model2)
```

    ## ___________________________________________________________________________
    ## Layer (type)                     Output Shape                  Param #     
    ## ===========================================================================
    ## dense_4 (Dense)                  (None, 100)                   78500       
    ## ___________________________________________________________________________
    ## dense_5 (Dense)                  (None, 50)                    5050        
    ## ___________________________________________________________________________
    ## dense_6 (Dense)                  (None, 10)                    510         
    ## ===========================================================================
    ## Total params: 84,060
    ## Trainable params: 84,060
    ## Non-trainable params: 0
    ## ___________________________________________________________________________

``` r
## compile the model with appropriate loss function, optimizer, and metrics
model2 %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)
## 
history <- model2 %>% fit(
  x_train, y_train, 
  epochs = 30, batch_size = 128, 
  validation_split = 0.2
)


y <- apply(y_test, 1, which.max)
for (i in 1:3) {
  layer_model <- keras_model(inputs = model2$input,
                             outputs = get_layer(model2, index = i)$output)
  output <- predict(layer_model, x_test)
  df <- as.data.frame(cbind(y, output))
  vars <- paste(colnames(df)[-1], collapse = " + ")
  ff <- as.formula(paste("y", vars, sep = " ~ "))
  mod <- lm(ff, data=df)
  # check for perfect multicollinearity (exact linear relationship)
  ld.vars <- attributes(alias(mod)$Complete)$dimnames[[1]]
  if (!is.null(ld.vars)) {
    print("Perfect Multicollinearity occurs! Following variables are omitted:")
    print(ld.vars)
    formula.new <- as.formula(
      paste(
        paste(deparse(ff), collapse=""), 
        paste(ld.vars, collapse="-"),
        sep="-"
      )
    )
    mod <- lm(formula.new,data = df)
  }
  print("==============vif===============")
  vifs <- vif(mod)
  print(vifs)
  print("Percentage of VIF of coefficients that > 10:")
  print(mean(ifelse(vifs > 10, 1, 0)))
  print("Average VIF for all coefficients:")
  print(mean(vifs))
  print("================================")
}
```

    ## [1] "==============vif==============="
    ##        V2        V3        V4        V5        V6        V7        V8 
    ##  4.536417  3.230979  4.127490  4.305446  4.774784  6.402200 10.028756 
    ##        V9       V10       V11       V12       V13       V14       V15 
    ##  4.089588  6.682010  7.244165  5.753806 10.847889  7.439206  3.473305 
    ##       V16       V17       V18       V19       V20       V21       V22 
    ##  3.584776  4.637334  5.705306  3.815649  5.812197  5.995771  6.129901 
    ##       V23       V24       V25       V26       V27       V28       V29 
    ##  4.622227  3.660444  4.094138  7.264888  5.417153  3.671242  3.740504 
    ##       V30       V31       V32       V33       V34       V35       V36 
    ##  9.626484  5.000729  3.898053  4.720821  3.351841  3.283668  3.358210 
    ##       V37       V38       V39       V40       V41       V42       V43 
    ##  8.781175  7.102801  6.194946  7.406476  6.125610  3.822708  3.098106 
    ##       V44       V45       V46       V47       V48       V49       V50 
    ##  6.778395  5.326217  2.728986  4.967324  5.409926 10.946443  3.620753 
    ##       V51       V52       V53       V54       V55       V56       V57 
    ##  3.959439  7.677329  4.808583  3.215412  9.041172  3.239798  2.800004 
    ##       V58       V59       V60       V61       V62       V63       V64 
    ##  6.293416  2.805388  3.149473  4.500305  6.194416  7.280527  4.561650 
    ##       V65       V66       V67       V68       V69       V70       V71 
    ##  6.491281  7.039305  3.711323  4.182697  7.447860  6.252582  4.008049 
    ##       V72       V73       V74       V75       V76       V77       V78 
    ##  2.593865  2.720629  4.461064  9.655582  4.844066  2.055755  3.793045 
    ##       V79       V80       V81       V82       V83       V84       V85 
    ##  3.434356  2.744296  3.984772  5.931394 12.355629 13.281328  3.347514 
    ##       V86       V87       V88       V89       V90       V91       V92 
    ## 10.540669  4.531244  5.418851  3.134328  4.965125  6.217644  4.968757 
    ##       V93       V94       V95       V96       V97       V98       V99 
    ##  4.840326  3.877585  3.514533  5.411828  3.796020  3.818810  3.106423 
    ##      V100      V101 
    ##  6.532001  3.340744 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 0.06
    ## [1] "Average VIF for all coefficients:"
    ## [1] 5.285154
    ## [1] "================================"
    ## [1] "==============vif==============="
    ##        V2        V3        V4        V5        V6        V7        V8 
    ## 10.846036 13.772822 12.343875 21.366445  7.284953 16.098389 17.920573 
    ##        V9       V10       V11       V12       V13       V14       V15 
    ## 10.097489 11.798889  8.196397 12.289129 12.301875  8.306099  4.239906 
    ##       V16       V17       V18       V19       V20       V21       V22 
    ##  4.162183 12.229239  7.678081  7.181902  6.953941  6.355570  7.677137 
    ##       V23       V24       V25       V26       V27       V28       V29 
    ##  8.084094 11.726607  7.415067 12.352868  1.385079 18.456983 11.220826 
    ##       V30       V31       V32       V33       V34       V35       V36 
    ## 15.083065  9.804379 13.355918  6.860266 13.842871  6.523605 20.068790 
    ##       V37       V38       V39       V40       V41       V42       V43 
    ##  7.880678 12.842620 10.950680 10.634979  5.640742  5.077633  6.989395 
    ##       V44       V45       V46       V47       V48       V49       V50 
    ## 16.546772  7.591993 16.529296 11.034068  5.945894  5.400857 10.885525 
    ##       V51 
    ##  5.517079 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 0.52
    ## [1] "Average VIF for all coefficients:"
    ## [1] 10.29499
    ## [1] "================================"
    ## [1] "Perfect Multicollinearity occurs! Following variables are omitted:"
    ## [1] "V11"
    ## [1] "==============vif==============="
    ##       V2       V3       V4       V5       V6       V7       V8       V9 
    ## 1.781580 1.881558 1.795468 1.801628 1.769343 1.705527 1.750261 1.793727 
    ##      V10 
    ## 1.748892 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 0
    ## [1] "Average VIF for all coefficients:"
    ## [1] 1.780887
    ## [1] "================================"

Third model
-----------

``` r
## Defining the Model # 256, 128, 10
model3 <- keras_model_sequential() 
model3 %>% 
  layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>% 
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dense(units = 10, activation = 'softmax')
summary(model3)
```

    ## ___________________________________________________________________________
    ## Layer (type)                     Output Shape                  Param #     
    ## ===========================================================================
    ## dense_7 (Dense)                  (None, 256)                   200960      
    ## ___________________________________________________________________________
    ## dense_8 (Dense)                  (None, 128)                   32896       
    ## ___________________________________________________________________________
    ## dense_9 (Dense)                  (None, 10)                    1290        
    ## ===========================================================================
    ## Total params: 235,146
    ## Trainable params: 235,146
    ## Non-trainable params: 0
    ## ___________________________________________________________________________

``` r
## compile the model with appropriate loss function, optimizer, and metrics
model3 %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)
## 
history <- model3 %>% fit(
  x_train, y_train, 
  epochs = 30, batch_size = 128, 
  validation_split = 0.2
)


y <- apply(y_test, 1, which.max)
for (i in 1:3) {
  layer_model <- keras_model(inputs = model3$input,
                             outputs = get_layer(model3, index = i)$output)
  output <- predict(layer_model, x_test)
  df <- as.data.frame(cbind(y, output))
  vars <- paste(colnames(df)[-1], collapse = " + ")
  ff <- as.formula(paste("y", vars, sep = " ~ "))
  mod <- lm(ff, data=df)
  # check for perfect multicollinearity (exact linear relationship)
  ld.vars <- attributes(alias(mod)$Complete)$dimnames[[1]]
  if (!is.null(ld.vars)) {
    print("Perfect Multicollinearity occurs! Following variables are omitted:")
    print(ld.vars)
    formula.new <- as.formula(
      paste(
        paste(deparse(ff), collapse=""), 
        paste(ld.vars, collapse="-"),
        sep="-"
      )
    )
    mod <- lm(formula.new,data = df)
  }
  print("==============vif===============")
  vifs <- vif(mod)
  print(vifs)
  print("Percentage of VIF of coefficients that > 10:")
  print(mean(ifelse(vifs > 10, 1, 0)))
  print("Average VIF for all coefficients:")
  print(mean(vifs))
  print("================================")
}
```

    ## [1] "Perfect Multicollinearity occurs! Following variables are omitted:"
    ## [1] "V102"
    ## [1] "==============vif==============="
    ##        V2        V3        V4        V5        V6        V7        V8 
    ##  7.547033  6.513234  4.069649  7.915364  3.502342  5.465531 12.507120 
    ##        V9       V10       V11       V12       V13       V14       V15 
    ##  8.554493  8.444453  3.669072  5.598126  7.771989  5.661204  3.768564 
    ##       V16       V17       V18       V19       V20       V21       V22 
    ##  5.227363  3.087849  4.995560  6.265299  9.206464  5.021572 16.305529 
    ##       V23       V24       V25       V26       V27       V28       V29 
    ##  7.101228  7.960705  7.948761  6.534378  5.689506  9.645052  5.186820 
    ##       V30       V31       V32       V33       V34       V35       V36 
    ##  4.376705  9.769265  3.001924  6.910572  2.808622  4.103921  7.899001 
    ##       V37       V38       V39       V40       V41       V42       V43 
    ##  5.744066 10.116794  3.533803  5.202781 11.799059  4.965052  5.571559 
    ##       V44       V45       V46       V47       V48       V49       V50 
    ##  6.227966  5.826665  8.234880  3.437473  6.228147  5.410958  7.443128 
    ##       V51       V52       V53       V54       V55       V56       V57 
    ##  7.712523  6.796994  4.730901  9.495448  5.066630  5.695868 13.406341 
    ##       V58       V59       V60       V61       V62       V63       V64 
    ##  3.815074  4.067019  3.389051  6.602862  3.766733  5.177445 10.508699 
    ##       V65       V66       V67       V68       V69       V70       V71 
    ##  8.879725  7.731211  5.515310  7.628080  8.127954  8.339617  3.890565 
    ##       V72       V73       V74       V75       V76       V77       V78 
    ## 13.436513  6.301568  7.058254  4.695315  5.195511  6.997324  6.005192 
    ##       V79       V80       V81       V82       V83       V84       V85 
    ##  6.114230  7.201220  5.513868  3.959806  4.814599  8.192576  7.621969 
    ##       V86       V87       V88       V89       V90       V91       V92 
    ##  9.346218  7.499185  7.545895  3.120206  4.074442  4.813683  8.341466 
    ##       V93       V94       V95       V96       V97       V98       V99 
    ## 11.574216 12.871022  4.771665  4.630879  5.835404  8.049018  5.170221 
    ##      V100      V101      V103      V104      V105      V106      V107 
    ##  6.475836  4.440871  4.809551  2.845559  5.596763 15.423140  5.223963 
    ##      V108      V109      V110      V111      V112      V113      V114 
    ## 10.632637  6.110746  5.904639  7.848232 12.205245  8.037750  4.505498 
    ##      V115      V116      V117      V118      V119      V120      V121 
    ##  4.867907  6.290393  8.530108  3.723166  5.990476  2.310604  8.067069 
    ##      V122      V123      V124      V125      V126      V127      V128 
    ##  7.055878  8.184113  2.923863  5.421809  6.111627  8.422874  7.253709 
    ##      V129      V130      V131      V132      V133      V134      V135 
    ##  7.038257  7.922829  4.430675  7.031835  7.821509  5.687344  6.284072 
    ##      V136      V137      V138      V139      V140      V141      V142 
    ##  3.170910  3.438306 13.553102  6.660775  5.161331  3.712514  5.196036 
    ##      V143      V144      V145      V146      V147      V148      V149 
    ##  4.658211  4.099229  7.229060  7.812272  9.763516  6.444775  4.283143 
    ##      V150      V151      V152      V153      V154      V155      V156 
    ##  6.223838  2.825768  7.437129  4.245556 11.531765  7.366747  4.899586 
    ##      V157      V158      V159      V160      V161      V162      V163 
    ##  6.377462  9.711965  4.071035  6.962878  9.376104  4.928534  7.690135 
    ##      V164      V165      V166      V167      V168      V169      V170 
    ##  5.347472 11.295683  6.166759  7.929744  4.171401  8.522638  4.783075 
    ##      V171      V172      V173      V174      V175      V176      V177 
    ##  3.990237  6.163269  4.239249  7.695539  5.667172 10.880113  2.880086 
    ##      V178      V179      V180      V181      V182      V183      V184 
    ##  4.719246  6.736637 16.154700  4.181774  4.630446  7.545938  5.801075 
    ##      V185      V186      V187      V188      V189      V190      V191 
    ##  7.077265  6.348092  4.525605  3.803269 10.959970  4.083886  5.804828 
    ##      V192      V193      V194      V195      V196      V197      V198 
    ##  4.119895  7.004536 11.785827  7.769033  4.659598  6.903631  4.866840 
    ##      V199      V200      V201      V202      V203      V204      V205 
    ##  7.054053  7.516372  4.605310  7.725655  7.750682  3.191425  5.289897 
    ##      V206      V207      V208      V209      V210      V211      V212 
    ##  9.401132 10.868065  9.736338  5.301821 13.726375  4.019811  6.391113 
    ##      V213      V214      V215      V216      V217      V218      V219 
    ## 18.197771  8.164252  8.215721  4.795025  5.277172  3.954625  4.799211 
    ##      V220      V221      V222      V223      V224      V225      V226 
    ##  5.642396 11.162088  4.391683  3.720002  4.625917  4.963581  4.985495 
    ##      V227      V228      V229      V230      V231      V232      V233 
    ##  6.585450  7.296652  7.865876  9.217068  6.147536  7.674988  6.740688 
    ##      V234      V235      V236      V237      V238      V239      V240 
    ##  3.703391  4.200370  9.708912  5.460654  3.876806  6.803296  5.179139 
    ##      V241      V242      V243      V244      V245      V246      V247 
    ##  7.473416 10.379884  8.861926  4.290288  4.714435  4.934725  6.044634 
    ##      V248      V249      V250      V251      V252      V253      V254 
    ##  3.716100  3.203370  8.056910  4.655986  9.086173  5.092475  5.933353 
    ##      V255      V256      V257 
    ##  3.455133  6.037580 21.852243 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 0.09803922
    ## [1] "Average VIF for all coefficients:"
    ## [1] 6.643417
    ## [1] "================================"
    ## [1] "Perfect Multicollinearity occurs! Following variables are omitted:"
    ## [1] "V4"  "V56"
    ## [1] "==============vif==============="
    ##        V2        V3        V5        V6        V7        V8        V9 
    ## 20.559649 20.361071 10.270367  9.192690 19.256604 15.474982  4.580832 
    ##       V10       V11       V12       V13       V14       V15       V16 
    ## 26.299853 40.231137 37.004754 35.452537 14.179395 30.380496 16.011718 
    ##       V17       V18       V19       V20       V21       V22       V23 
    ## 14.662929 25.699650  5.289959 19.608441  7.869850 20.812494 15.032053 
    ##       V24       V25       V26       V27       V28       V29       V30 
    ## 13.396154 26.602284 18.647683 26.088508 18.417596 13.456200 25.196122 
    ##       V31       V32       V33       V34       V35       V36       V37 
    ## 17.795821 21.517206 12.967323 21.194181 25.458226 20.903212 10.099288 
    ##       V38       V39       V40       V41       V42       V43       V44 
    ## 14.046666 15.812183 12.320055 13.404284 16.423573 34.362176 11.660946 
    ##       V45       V46       V47       V48       V49       V50       V51 
    ## 26.089189 33.435067 22.578453 23.347751 27.439175 17.946961 16.370051 
    ##       V52       V53       V54       V55       V57       V58       V59 
    ## 29.758385 17.918800 26.565896 19.585444 19.729389 15.220741 14.764995 
    ##       V60       V61       V62       V63       V64       V65       V66 
    ## 21.395807 31.191659 32.688599  8.697531 26.634007 44.240415 12.866791 
    ##       V67       V68       V69       V70       V71       V72       V73 
    ## 15.944637 19.317399  4.118428 25.458769 16.577824 17.266123 13.350489 
    ##       V74       V75       V76       V77       V78       V79       V80 
    ## 17.723691 19.929432 25.240387 17.578826 31.002294 25.751870 16.979388 
    ##       V81       V82       V83       V84       V85       V86       V87 
    ## 10.477069 19.112987 14.095062 10.356505 15.124576  7.458330 33.823014 
    ##       V88       V89       V90       V91       V92       V93       V94 
    ## 12.374685 14.829383 16.670377 14.601876 13.784858 17.329695 19.187177 
    ##       V95       V96       V97       V98       V99      V100      V101 
    ## 16.656211 11.205040 10.106240  9.819206 25.587802 11.552029  8.222616 
    ##      V102      V103      V104      V105      V106      V107      V108 
    ## 27.298089 10.317012 24.852331  9.662904  8.671768 21.006508 17.691322 
    ##      V109      V110      V111      V112      V113      V114      V115 
    ## 12.266193 15.682918 17.829014 17.594165 27.946900 15.858070  8.066630 
    ##      V116      V117      V118      V119      V120      V121      V122 
    ## 26.756921 13.133092  8.221835 11.817421 12.762916 18.880186 20.977540 
    ##      V123      V124      V125      V126      V127      V128      V129 
    ##  8.821366 15.661299 34.263952 14.945081 33.484862 15.836441  9.308281 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 0.8809524
    ## [1] "Average VIF for all coefficients:"
    ## [1] 18.60869
    ## [1] "================================"
    ## [1] "Perfect Multicollinearity occurs! Following variables are omitted:"
    ## [1] "V11"
    ## [1] "==============vif==============="
    ##       V2       V3       V4       V5       V6       V7       V8       V9 
    ## 1.766018 1.861174 1.778313 1.781505 1.742850 1.691705 1.745373 1.766380 
    ##      V10 
    ## 1.745345 
    ## [1] "Percentage of VIF of coefficients that > 10:"
    ## [1] 0
    ## [1] "Average VIF for all coefficients:"
    ## [1] 1.764296
    ## [1] "================================"
