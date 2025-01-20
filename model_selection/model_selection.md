Model Selection
================

# Import libraries

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
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

Build custom AUC function to extract AUC from the caret model object

``` r
eval_mod <- function(model, data) {
  pred <- predict(model, data)
  cm <- caret::confusionMatrix(pred, data$classes, positive="malignant")
  auc <- roc(data$classes,
             predict(model, data, type = "prob")[, "malignant"]) %>% auc()
  result <- c(cm$overall["Accuracy"],cm$byClass['Sensitivity'], cm$byClass['Specificity'], cm$byClass['F1'],AUC=auc)
  return(result)
}
```

# Read clean data

``` r
bc_data <- readRDS("../EDA/bc_clean.RDS")
bc_data$classes <- as.factor(bc_data$classes)
```

``` r
set.seed(2025)
index <- caret::createDataPartition(bc_data$classes, p = 0.7, list = FALSE)

train_data <- bc_data[index, ]
test_data  <- bc_data[-index, ]
```

Try doing backward selection using AIC first, then using BIC. Do you get
the same result?

``` r
M0 = glm(classes ~ 1, data = train_data, family = binomial)  # Null model
M1 = glm(classes ~ ., data = train_data, family= binomial)  # Full model
step.back.aic <- stats::step(M1, direction = "backward", trace = F, k = 2)
step.back.bic <-  stats::step(M1, direction = "backward", trace = F, k = log(dim(train_data)[1]))
summary(step.back.aic)
```

    ## 
    ## Call:
    ## glm(formula = classes ~ clump_thickness + marginal_adhesion + 
    ##     bare_nuclei + bland_chromatin + normal_nucleoli + mitosis, 
    ##     family = binomial, data = train_data)
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -12.1585     1.8981  -6.406 1.50e-10 ***
    ## clump_thickness     0.7900     0.1926   4.102 4.09e-05 ***
    ## marginal_adhesion   0.3478     0.1509   2.305  0.02117 *  
    ## bare_nuclei         0.5846     0.1351   4.328 1.50e-05 ***
    ## bland_chromatin     0.5789     0.1892   3.060  0.00222 ** 
    ## normal_nucleoli     0.3746     0.1385   2.705  0.00683 ** 
    ## mitosis             0.9787     0.3777   2.591  0.00957 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 578.843  on 440  degrees of freedom
    ## Residual deviance:  57.727  on 434  degrees of freedom
    ## AIC: 71.727
    ## 
    ## Number of Fisher Scoring iterations: 9

``` r
summary(step.back.bic)
```

    ## 
    ## Call:
    ## glm(formula = classes ~ clump_thickness + bare_nuclei + bland_chromatin + 
    ##     normal_nucleoli, family = binomial, data = train_data)
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -10.6522     1.4364  -7.416 1.21e-13 ***
    ## clump_thickness   0.9051     0.1764   5.131 2.88e-07 ***
    ## bare_nuclei       0.6289     0.1277   4.927 8.37e-07 ***
    ## bland_chromatin   0.5791     0.1719   3.369 0.000755 ***
    ## normal_nucleoli   0.3855     0.1265   3.047 0.002310 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 578.84  on 440  degrees of freedom
    ## Residual deviance:  69.02  on 436  degrees of freedom
    ## AIC: 79.02
    ## 
    ## Number of Fisher Scoring iterations: 8

Now Forward

``` r
step.fwd.aic <-  stats::step(M0, scope = list(lower = M0, upper = M1), direction = "forward", trace = FALSE, k = 2)
summary(step.fwd.aic)
```

    ## 
    ## Call:
    ## glm(formula = classes ~ uniformity_of_cell_size + bare_nuclei + 
    ##     clump_thickness + normal_nucleoli + bland_chromatin + mitosis + 
    ##     marginal_adhesion, family = binomial, data = train_data)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)             -11.7875     1.9322  -6.101 1.06e-09 ***
    ## uniformity_of_cell_size   0.1373     0.2026   0.678 0.498039    
    ## bare_nuclei               0.5634     0.1380   4.083 4.45e-05 ***
    ## clump_thickness           0.7389     0.2066   3.576 0.000349 ***
    ## normal_nucleoli           0.3381     0.1480   2.284 0.022356 *  
    ## bland_chromatin           0.5033     0.2220   2.268 0.023359 *  
    ## mitosis                   0.9395     0.3813   2.464 0.013745 *  
    ## marginal_adhesion         0.3240     0.1538   2.106 0.035173 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 578.843  on 440  degrees of freedom
    ## Residual deviance:  57.246  on 433  degrees of freedom
    ## AIC: 73.246
    ## 
    ## Number of Fisher Scoring iterations: 9

``` r
step.fwd.bic <-  stats::step(M0,scope = list(lower = M0, upper = M1), direction = "forward", trace = FALSE, k = log(dim(train_data)[1]))
summary(step.fwd.bic)
```

    ## 
    ## Call:
    ## glm(formula = classes ~ uniformity_of_cell_size + bare_nuclei + 
    ##     clump_thickness + normal_nucleoli, family = binomial, data = train_data)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -9.2420     1.2703  -7.275 3.46e-13 ***
    ## uniformity_of_cell_size   0.5505     0.1719   3.202  0.00137 ** 
    ## bare_nuclei               0.6126     0.1233   4.967 6.81e-07 ***
    ## clump_thickness           0.7213     0.1849   3.901 9.60e-05 ***
    ## normal_nucleoli           0.3645     0.1275   2.858  0.00426 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 578.843  on 440  degrees of freedom
    ## Residual deviance:  69.977  on 436  degrees of freedom
    ## AIC: 79.977
    ## 
    ## Number of Fisher Scoring iterations: 8

Exhaustive searches

``` r
library(bestglm)
```

    ## Loading required package: leaps

``` r
x_tr <- train_data 
x_tr$classes <- ifelse(train_data$classes=="malignant",1,0)

res.bestglm <-bestglm::bestglm(Xy = x_tr,
                               family = binomial,
                               IC = "AIC",     # Information criteria for
                               method = "exhaustive")
```

    ## Morgan-Tatar search since family is non-gaussian.

# Show top 5 models

``` r
res.bestglm$BestModels
```

    ##   clump_thickness uniformity_of_cell_size uniformity_of_cell_shape
    ## 1            TRUE                   FALSE                    FALSE
    ## 2            TRUE                   FALSE                    FALSE
    ## 3            TRUE                    TRUE                    FALSE
    ## 4            TRUE                    TRUE                    FALSE
    ## 5            TRUE                   FALSE                     TRUE
    ##   marginal_adhesion single_epithelial_cell_size bare_nuclei bland_chromatin
    ## 1              TRUE                       FALSE        TRUE            TRUE
    ## 2              TRUE                        TRUE        TRUE            TRUE
    ## 3              TRUE                       FALSE        TRUE            TRUE
    ## 4              TRUE                        TRUE        TRUE            TRUE
    ## 5              TRUE                       FALSE        TRUE            TRUE
    ##   normal_nucleoli mitosis Criterion
    ## 1            TRUE    TRUE  69.72746
    ## 2            TRUE    TRUE  70.14291
    ## 3            TRUE    TRUE  71.24577
    ## 4            TRUE    TRUE  71.44153
    ## 5            TRUE    TRUE  71.48178

``` r
summary(res.bestglm$BestModel)
```

    ## 
    ## Call:
    ## glm(formula = y ~ ., family = family, data = Xi, weights = weights)
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -12.1585     1.8981  -6.406 1.50e-10 ***
    ## clump_thickness     0.7900     0.1926   4.102 4.09e-05 ***
    ## marginal_adhesion   0.3478     0.1509   2.305  0.02117 *  
    ## bare_nuclei         0.5846     0.1351   4.328 1.50e-05 ***
    ## bland_chromatin     0.5789     0.1892   3.060  0.00222 ** 
    ## normal_nucleoli     0.3746     0.1385   2.705  0.00683 ** 
    ## mitosis             0.9787     0.3777   2.591  0.00957 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 578.843  on 440  degrees of freedom
    ## Residual deviance:  57.727  on 434  degrees of freedom
    ## AIC: 71.727
    ## 
    ## Number of Fisher Scoring iterations: 9

Same in caret

``` r
trainX <- train_data[,-10]
trainY <- train_data$classes

cctrl <- trainControl(method = "repeatedcv", 
                       number = 5, 
                       repeats = 3,  
                       savePredictions = TRUE,
                       summaryFunction = twoClassSummary,
                      classProbs = TRUE
)
```

``` r
set.seed(849)
caret_full <- train(classes~.,
                   data = train_data,
                   method = "glm", 
                   family = "binomial",
                   trControl = cctrl,
                   preProc = c("center", "scale"),
                   metric = "ROC",
                   trace = 0)
full <- eval_mod(caret_full, test_data)
```

    ## Setting levels: control = benign, case = malignant

    ## Setting direction: controls < cases

``` r
caret_fwd <- train(classes~.,
                             data = train_data,
                             method = "glmStepAIC", 
                             direction = "forward",
                             trControl = cctrl,
                             preProc = c("center", "scale"),
                             metric = "ROC",
                             trace = 0)
fwd <- eval_mod(caret_fwd, test_data)
```

    ## Setting levels: control = benign, case = malignant

    ## Setting direction: controls < cases

``` r
caret_back <- train(classes~.,
                   data = train_data,
                   method = "glmStepAIC", 
                   direction = "backward",
                   trControl = cctrl,
                   preProc = c("center", "scale"),
                   metric = "ROC",
                   trace = 0)
back <- eval_mod(caret_back, test_data)
```

``` r
?MASS::stepAIC
?stats::step
```

``` r
caret_both <- train(classes~.,
                    data = train_data,
                    method = "glmStepAIC", 
                    direction = "both",
                    trControl = cctrl,
                    preProc = c("center", "scale"),
                    metric = "ROC",
                    trace = 0)
both <- eval_mod(caret_both, test_data)
```

    ## Setting levels: control = benign, case = malignant

    ## Setting direction: controls < cases

``` r
rbind(full, fwd, back, both)
```

    ##       Accuracy Sensitivity Specificity        F1       AUC
    ## full 0.9787234   0.9565217   0.9915966 0.9705882 0.9923274
    ## fwd  0.9787234   0.9565217   0.9915966 0.9705882 0.9945195
    ## back 0.9787234   0.9565217   0.9915966 0.9705882 0.9939106
    ## both 0.9787234   0.9565217   0.9915966 0.9705882 0.9939106
