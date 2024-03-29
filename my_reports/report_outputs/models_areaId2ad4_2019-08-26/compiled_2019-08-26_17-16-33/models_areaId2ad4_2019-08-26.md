Models For Predicting Commuting Flow From areaId to admin4 units
================
Zhanhao Zhang
8/23/2019

## Introduction

We are going to explore a few models that can predict the commuting
flows from areaId to admin4 units. Attempted models are: gravity model,
disambiguation model, randomForest algorithm, gradient boosting
algorithm, and features mapping method.

### Libraries

``` r
library(here)
library(randomForest)
library(gbm)
library(data.table)

source(here("scripts/extract_models.R"))
```

### Required Data

To get these models work, we need a cleaned data that represents the
commuting flows from areaId to admin4. Required columns are listed below
(additional columns are optional and have no effect on the models
below):  
areaId: the areaId of the source.  
admin4Id: the admin4Id of the destination.  
admin2: the admin2 name of the destination.  
aid.pop: the population at the source areaId.  
ad4.pop: the population at the destination admin4.  
ad2.pop: the population at the destination admin2.  
d1: the geographic distance between source areaId and destination
admin2.  
d2: the geographic distance between source areaId and destination
admin4.  
w: the commuting flow from source areaId to destination admin4.

### Obtaining Data

``` r
get_data <- function(data_clean){
  dat <- data_clean[, c("areaId", "admin4Id", "admin2", "aid.pop", "ad4.pop", "d1", "d2", "w")]
  colnames(dat) <- c("areaId", "admin4Id", "admin2", "N1", "N2", "d1", "d", "w")
  dat <- dat[dat$N1 > 0,]
  idx <- c()
  ad2 <- sort(unique(dat$admin2))
  for(i in 1:nrow(dat)){
    idx <- c(idx, which(ad2 == dat[i,]$admin2))
  }
  dat <- data.frame(cbind(dat, data.frame(idx = idx)))
  N2 <- c()
  for(ad in ad2){
    pop <- unique(data_clean[data_clean$admin2==ad,]$ad2.pop)
    N2 <- c(N2, pop)
  }
  return(list(dat = dat, N2 = N2))
}

new_data_clean <- fread(here("data/clean/", "data_clean_pop2018.csv"), sep = ",",
                        header = T)
  
res <- get_data(new_data_clean)
dat <- data.frame(res$dat)
N2 <- res$N2

rel_gravity_dat_covs <- fread(here("data/clean", "rel_gravity_dat_covs.csv"),
                                   sep = ",", header = T)
rel_gravity_dat_covs <- data.frame(rel_gravity_dat_covs)
```

### Gravity Model

Gravity Model typically trains directly on the commuting data from
areaId to admin4 units in 2018. The covariates are: N1 (source
population), N2 (destination population), and d (geographic distance
between areaId and admin4). There are quite a few hyperparameters for
the gravity model:  
1\. poisson regression VS negative binomial regression  
2\. training on commuting data from 2015-2017 VS from 2018  
3\. the cutoff on geographic distance to split the data with two
separate models  
4\. exponential dependency VS power dependency on distance for data
below cutoff

Below is the super function to obtain the gravity model while users can
specify the above four
hyperparameters.

``` r
super_fit.2018 <- function(dat, model_name, use_prev = T, image_name = NULL,
                           cutoff = 0, base_dir = "./", exp_dependence_far = F,
                           exp_dependence_near = T, show_image = F){
  ###
  # Inputs:
  #   dat: the data to fit
  #   model_name: the name of model to use. Valid names are poisson and negbin.
  #   use_prev: whether or not to use the pretrained model on the dataset whose
  #             commuting flows are obtained in 2015-2017. DEFAULT = TRUE.
  #   image_name: the name of image that visualize the model performance. 
  #               DEFAULT = NULL.
  #   cutoff: the cutoff on geographic distance where we will use two models on
  #           data points from both sides of the cutoff. DEFAULT = 0, which means
  #           no cutoff. In the case of no cutoff, all data points are considered
  #           to be further than the cutoff, so you should only modify 
  #           exp_dependence_far if you want to choose power/exponential dependency
  #           on the geographic distance.
  #   base_dir: the directory to place the image that shows the model performance.
  #             DEFAULT = ./, which is the current directory.
  #   exp_dependence_far: whether to use expoential dependency or power dependency
  #                       on geographical distance for distance above the cutoff.
  #                       DEFAULT = FALSE, which means using power dependency.
  #   exp_dependence_near: whether to use expoential dependency or power dependency
  #                        on geographical distance for distance above the cutoff.
  #                        DEFAULT = TRUE, which means using exoonential dependency.
  #   show_image: whether to plot the image that visualize the model performance.
  #               DEFAULT = FALSE.
  #
  # Output:
  #   ssr: sum of residual squares of the current selected model on the given dataset.
  #   model_far: the fitted model for data points whose distances are above the cutoff.
  #   model_near: the fitted model for data points whose distances are below the cutoff. 
  ###
  
  if(model_name != "poisson" && model_name != "negbin"){
    cat("Usage: model_name = c(\"poisson\", \"negbin\")")
    return(NULL)
  }
  dat_near <- dat[dat$d < cutoff,]
  dat_far <- dat[dat$d >= cutoff,]
  
  pred_far <- NULL
  pred_near <- NULL
  
  if(use_prev){
    res <- super_model_fits(rel_gravity_dat_covs, model_name, image_name = NULL,
                            show_plots = F, cutoff = cutoff, use_time = F,
                            exp_dependence_far = exp_dependence_far,
                            exp_dependence_near = exp_dependence_near)
    pred_far <- predict(res$model_far, dat_far, type = "response")
    if(cutoff > 0){
      pred_near <- predict(res$model_near, dat_near, type = "response")
    }
  } else{
    
    formula_far <- w ~ log(N1) + log(N2) + d
    if(!exp_dependence_far){
      formula_far <- w ~ log(N1) + log(N2) + log(d)
    }
    formula_near <- w ~ log(N1) + log(N2) + d
    if(!exp_dependence_near){
      formula_near <- w ~ log(N1) + log(N2) + log(d)
    }
    
    model_far <- NULL
    model_near <- NULL
    if(model_name == "poisson"){
      model_far <- glm(data = dat_far, formula_far, family = "poisson")
      if(cutoff > 0){
        model_near <- glm(data = dat_near, formula_near, family = "poisson")
      }
    } else{
      model_far <- glm.nb(data = dat_far, formula_far)
      if(cutoff > 0){
        model_near <- glm.nb(data = dat_near, formula_near)
      }
    }
    pred_far <- model_far$fitted.values
    if(cutoff > 0){
      pred_near <- model_near$fitted.values
    }
  }
  
  dat_toPlot <- data.frame(N1 = c(dat_far$N1, dat_near$N1),
                           N2 = c(dat_far$N2, dat_near$N2),
                           d = c(dat_far$d, dat_near$d),
                           w = c(dat_far$w, dat_near$w),
                           pred = c(pred_far, pred_near))
  if(show_image){
    cor_plot_all(paste(image_name, "\n", sep=""), dat_toPlot, 
               paste(base_dir, image_name, sep = ""))
  }
  ssr <- sum((dat_toPlot$pred - dat_toPlot$w)^2)
  return(list(ssr = ssr, model_far = model_far, model_near = model_near))
}
```

The best models fit trained on dataset from 2015-2017 and from 2018 are
below:

Trained on 2015-2017 (exponential dependency on geographic distance for
data points whose distance below the cutoff, use negative binommial
regression, cutoff at 50km). SSR = 5662.416.
<img src="/Users/zhangji/Desktop/others/IHME/Macro/TravelingModel/traveling_fit/my_reports/data/images/experiments_areaId2ad4/plots/negbin_trainData=old_cutoff=50k_near=exp.png" title="Trained on 2015-2017" alt="Trained on 2015-2017" width="60%" style="display: block; margin: auto;" />
Trained on 2018 (exponential dependency on geographic distance for data
points whose distance below the cutoff, use poisson regression, cutoff
at 15km). SSR = 3635.204.
<img src="/Users/zhangji/Desktop/others/IHME/Macro/TravelingModel/traveling_fit/my_reports/data/images/experiments_areaId2ad4/plots/poisson_trainData=2018_cutoff=15k_near=exp.png" title="Trained on 2018" alt="Trained on 2018" width="60%" style="display: block; margin: auto;" />

The correlation plots for these two models
are:  
<img src="/Users/zhangji/Desktop/others/IHME/Macro/TravelingModel/traveling_fit/my_reports/data/images/experiments_areaId2ad4/plots/scatter_old.png" title="caption" alt="caption" width="50%" /><img src="/Users/zhangji/Desktop/others/IHME/Macro/TravelingModel/traveling_fit/my_reports/data/images/experiments_areaId2ad4/plots/scatter_new.png" title="caption" alt="caption" width="50%" />

### Disambiguation Model

Disambiguation Model is comprised of two parts: first map areaId to
admin2 units, then using the fact that we are leaving from a known
admin2 unit, calculate the probability of people traveling to each of
its admin4 units based on their population. The first part is fitted
using gravity model, either through commuting data from 2015-2017, or
through data from 2018. The second part is fitted based on pop^k, where
pop is the population at each admin4 unit. After several trials of
training using Particle Markov Chain Monte Carlo algorithm, I get a set
of converged k for each of the admin2 units: 0.2345718 for Baney,
0.5415359 for Luba, 0.4763921 for Malabo, and 2.4926470 for Riaba.

The correlation plots for the two parts of the disambiguation model
are:  
<img src="/Users/zhangji/Desktop/others/IHME/Macro/TravelingModel/traveling_fit/my_reports/data/images/experiments_areaId2ad2/plots/disambiguation_fir.png" title="caption" alt="caption" width="50%" /><img src="/Users/zhangji/Desktop/others/IHME/Macro/TravelingModel/traveling_fit/my_reports/data/images/experiments_areaId2ad2/plots/disambiguation_sec.png" title="caption" alt="caption" width="50%" />
Though the fit looks reasonable for now, this model works very poorly
for predicting commuting flow from areaId to admin4. One of the possible
reasons is that the probability of traveling to each admin4 unit is
proportional to pop^k (according to the assumption of the disambiguation
model), which smooth out the probability for each admin4 unit. In other
word, probabilities of admin4 units within the same admin2 area are
close to each other. However, in reality, there are usually very few
amount of people depart from each areaId, while there are a lot of
admin4 units within an admin2 area, so there is only a tiny amount of
commuting flow allocated to each admin4 unit. The correlation plot of
commuting flow from areaId to admin4 looks like the following (SSR =
8014.962):  
<img src="/Users/zhangji/Desktop/others/IHME/Macro/TravelingModel/traveling_fit/my_reports/data/images/experiments_areaId2ad2/plots/disambiguation_all.png" title="caption" alt="caption" width="50%" style="display: block; margin: auto;" />

### RandomForest Algorithm

``` r
rf <- randomForest(x = dat[, c("N1", "N2", "d")], y = dat$w, ntree = 5000)
ssr <- sum((rf$predicted - dat$w)^2)
plot(rf$predicted, dat$w, main = "randomForest performance",
     xlab = "predicted values", ylab = "original commuting flow",
     xlim = c(0, 20))
abline(a = 0, b = 1, col = "red")
```

![](models_areaId2ad4_2019-08-26_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
print(ssr)
```

    ## [1] 3805.317

### Gradient Boosting Algorithm

``` r
gb <- gbm(w ~ log(N1) + log(N2) + log(d), data = dat, n.trees = 20000,
          interaction.depth = 3, train.fraction = 0.8, cv.folds = 5)
```

    ## Distribution not specified, assuming gaussian ...

``` r
ssr <- sum((gb$fit - dat$w)^2)
plot(gb$fit, dat$w, main = "Gradient Boosting Performance", xlim = c(0, 20),
     xlab = "predicted values", ylab = "original commuting flows")
abline(a = 0, b = 1, col = "red")
```

![](models_areaId2ad4_2019-08-26_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
print(ssr)
```

    ## [1] 1269.072

You may notice that the Gradient Boosting Algorithm gives the best fit,
for it gives the lowest SSR. However, if we take a closer look by
splitting up the SSR for training set and validation set, we will
realize that the story is completely
different.

``` r
ssr_train <- gb$train.error[length(gb$train.error)] * nrow(dat) * gb$train.fraction
ssr_valid <- gb$valid.error[length(gb$valid.error)] * nrow(dat) * (1 - gb$train.fraction)
print(paste("SSR_TRAIN =", ssr_train, "SSR_VALID =", ssr_valid))
```

    ## [1] "SSR_TRAIN = 19.2333700442096 SSR_VALID = 1248.76027494515"

As you can see, the gradient boosting algorithm suffers from severe
overfitting.

### Features Mapping Method

The last approach, which is also subject to overfitting, is using the
idea of features mapping, where we map three features (source
population, destination population, and geographic distance) into
hundreds or even thousands of features, and then doing a regression on
those features for prediction. One thing to be aware of is that using
linear regression may give negative predictions, so one way to deal with
this is to use abs() or relu() to wrap the predicted values. Both
methods give similar results, where training accuracy is super high,
while the predicted values for the validation set can hardly fit along
the red line. Poisson regression has also be tried, but the sum of
residuals squares shoots up very quickly for the validation set as we
increase the number of features to get a better fit for the training
set. Thus, so far the best choice for the features mapping method is to
use linear regression and then a abs() or relu() as an activation
function.

``` r
feature_vector <- function(num_features_vec, dat_train, dat_valid, init_sd = 1){
  ssr_train_all <- c()
  ssr_valid_all <- c()
  ssr_total <- 10 ^ 6
  
  pre_train <- NULL
  pre_valid <- NULL
  
  for(num_features in num_features_vec){
    set.seed(123)
    transform_mat <- matrix(rnorm(3 * num_features, sd = init_sd), nrow = 3)
    
    train_orig_covs <- as.matrix(dat_train[, c("N1", "N2", "d")])
    train_features <- sin(train_orig_covs %*% transform_mat)
    train_dat_featured <- data.frame(cbind(data.frame(w = dat_train$w), train_features))
    valid_orig_covs <- as.matrix(dat_valid[, c("N1", "N2", "d")])
    valid_features <- sin(valid_orig_covs %*% transform_mat)
    valid_dat_featured <- data.frame(cbind(data.frame(w = dat_valid$w), valid_features))
    
    lm_fet <- lm(data = train_dat_featured, w ~ .)
    ssr_train <- sum((abs(lm_fet$fitted.values) - train_dat_featured$w)^2)
    pred_valid <- predict(lm_fet, valid_dat_featured, type = "response")
    ssr_valid <- sum((pred_valid - valid_dat_featured$w)^2)
    
    if(ssr_total > ssr_train + ssr_valid){
      ssr_total <- ssr_train + ssr_valid
      pre_train <- abs(lm_fet$fitted.values)
      pre_valid <- pred_valid
    }
    
    ssr_train_all <- c(ssr_train_all, ssr_train)
    ssr_valid_all <- c(ssr_valid_all, ssr_valid)
  }
  
  plot(num_features_vec, ssr_train_all, type = "l", lwd = 2, col = "blue",
       main = "Training and Validation SSR vs #vectors", xlab = "#vectors",
       ylab = "SSR", ylim = c(0, max(ssr_train_all, ssr_valid_all)))
  lines(num_features_vec, ssr_valid_all, lwd = 2, col = "red")
  lines(num_features_vec, ssr_train_all + ssr_valid_all, lwd = 2, col = "green")
  legend("topleft", c("train SSR", "valid SSR", "total SSR"), 
         col = c("blue", "red", "green"),
         bg = "gray", cex = 0.6, lwd = 2)
  return(list(pred_train = pre_train, pred_valid = pred_valid))
}

dat_tofit <- dat
colnames(dat_tofit) <- c("areaId", "admin4Id", "admin2", "N1", "N2", "d1",
                         "d", "w", "idx")

split_tofit <- sample(1:nrow(dat_tofit), 0.8 * nrow(dat_tofit))
dat_tofit_train <- dat_tofit[split_tofit,]
dat_tofit_valid <- dat_tofit[-split_tofit,]
fv_res <- feature_vector(c((1:40) * 20), dat_tofit_train, dat_tofit_valid,
                         init_sd = 1)
```

![](models_areaId2ad4_2019-08-26_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
par(mfrow = c(2, 1))
plot(fv_res$pred_train, dat_tofit_train$w, main = "training set",
     xlim = c(0, max(fv_res$pred_train, dat_tofit_train$w)),
     ylim = c(0, max(fv_res$pred_train, dat_tofit_train$w)))
abline(a = 0, b = 1, col = "red")
plot(fv_res$pred_valid, dat_tofit_valid$w, main = "validation set",
     xlim = c(0, max(fv_res$pred_valid, dat_tofit_valid$w)),
     ylim = c(0, max(fv_res$pred_valid, dat_tofit_valid$w)))
abline(a = 0, b = 1, col = "red")
```

![](models_areaId2ad4_2019-08-26_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
par(mfrow = c(1, 1))
```
