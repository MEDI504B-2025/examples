Deep learning
================

``` r
# Read clean data
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(caret)
```

    ## Loading required package: lattice
    ## 
    ## Attaching package: 'caret'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

``` r
library(pROC)
```

    ## Type 'citation("pROC")' for a citation.
    ## 
    ## Attaching package: 'pROC'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

``` r
# Build custom AUC function to extract AUC
# from the caret model object
eval_mod <- function(model, data) {
  pred <- predict(model, data)
  cm <- caret::confusionMatrix(pred, data$classes, positive="malignant")
  auc <- roc(data$classes,
             predict(model, data, type = "prob")[, "malignant"]) %>% auc()
  result <- c(cm$overall["Accuracy"],cm$byClass['Sensitivity'], cm$byClass['Specificity'], cm$byClass['F1'],AUC=auc)
  return(result)
}


bc_data <- readRDS("../EDA/bc_clean.RDS")
bc_data$classes <- as.factor(bc_data$classes)

set.seed(2024)
index <- caret::createDataPartition(bc_data$classes, p = 0.7, list = FALSE)
```

``` r
y = as.matrix(bc_data[,10])
y[which(y=="benign")] = 0
y[which(y=="malignant")] = 1
y = as.numeric(y)
x = as.numeric(as.matrix(bc_data[,1:9]))
x = matrix(as.numeric(x),ncol=9)

y_train <- y[index]
y_test <- y[-index]
x_train  <- x[index, ]
x_test <- x[-index, ]
```

``` r
library(deepnet)
set.seed(2024)

nn <- nn.train(x_train, y_train, hidden = c(10))
yy = nn.predict(nn, x_test)
print(head(yy))
```

    ##           [,1]
    ## [1,] 0.4666505
    ## [2,] 0.3203479
    ## [3,] 0.4602568
    ## [4,] 0.3169091
    ## [5,] 0.3960531
    ## [6,] 0.4044650

``` r
?nn.train
```

``` r
yhat = matrix(0,length(yy),1)
yhat[which(yy > mean(yy))] = 1
yhat[which(yy <= mean(yy))] = 0
cm = caret::confusionMatrix(factor(yhat),factor(y_test))
print(cm)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 114   0
    ##          1   5  69
    ##                                          
    ##                Accuracy : 0.9734         
    ##                  95% CI : (0.939, 0.9913)
    ##     No Information Rate : 0.633          
    ##     P-Value [Acc > NIR] : < 2e-16        
    ##                                          
    ##                   Kappa : 0.9436         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.07364        
    ##                                          
    ##             Sensitivity : 0.9580         
    ##             Specificity : 1.0000         
    ##          Pos Pred Value : 1.0000         
    ##          Neg Pred Value : 0.9324         
    ##              Prevalence : 0.6330         
    ##          Detection Rate : 0.6064         
    ##    Detection Prevalence : 0.6064         
    ##       Balanced Accuracy : 0.9790         
    ##                                          
    ##        'Positive' Class : 0              
    ## 
