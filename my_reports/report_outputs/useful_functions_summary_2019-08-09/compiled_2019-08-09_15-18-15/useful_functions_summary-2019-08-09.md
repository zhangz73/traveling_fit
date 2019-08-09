Useful Functions Summary
================
Zhanhao Zhang
8/8/2019

## Data Transformation Helper Function

``` r
get_test_dat_mul <- function(test_dat_mul_with_w_filename){
  ###
  # get_test_dat_mul() returns the test_dat_mul to feed in multinomial
  # regression model for prediction
  # 
  # Input:
  #   test_dat_mul_with_w_filename: the filename of test_dat_mul_with_w that should
  #                                 contain everything useful for rel_dat_mul and
  #                                 test_dat_mul
  #
  # Output:
  #   The test_dat_mul to feed in multinomial regression model for prediction
  ###
  print("parsing data for prediction...")
  test_dat_mul_with_w <- fread(here("data/clean", test_dat_mul_with_w_filename), 
                               sep = ",", header = T)
  test_dat_mul_with_w <- data.frame(test_dat_mul_with_w)
  test_dat_mul <- test_dat_mul_with_w[,1:25]
  return(test_dat_mul)
}

get_rel_dat_mul <- function(test_dat_mul_with_w_filename){
  ###
  # get_rel_dat_mul() returns the rel_dat_mul to feed in multinomial
  # regression model for model fitting
  # 
  # Input:
  #   test_dat_mul_with_w_filename: the filename of test_dat_mul_with_w that should
  #                                 contain everything useful for rel_dat_mul and
  #                                 test_dat_mul
  #
  # Output:
  #   The rel_dat_mul to feed in multinomial regression model for fitting
  ###
  
  print("parsing data for model fitting...")
  test_dat_mul_with_w <- fread(here("data/clean", test_dat_mul_with_w_filename), 
                               sep = ",", header = T)
  test_dat_mul_with_w <- data.frame(test_dat_mul_with_w)
  centrals <- c("to", "ti_ban", "ti_mal", "ti_lub", "ti_ria", "ti_mok", "ti_ure")
  test_dat_mul <- test_dat_mul_with_w[,1:25]
  test_dat_mul <- data.frame(test_dat_mul)
  rel_dat_mul <- c()
  test_dat_mul_with_w <- data.table(test_dat_mul_with_w)
  for(c in centrals){
    tmp <- test_dat_mul_with_w[test_dat_mul_with_w[[c]] > 0,]
    dat <- tmp[rep(seq(1, nrow(tmp)), tmp[[c]])]
    lst_dest <- rep(c, nrow(dat))
    dat <- cbind(data.frame(dat[,1:25]), data.frame(dest = lst_dest))
    rel_dat_mul <- data.frame(rbind(data.frame(rel_dat_mul), data.frame(dat)))
  }
  return(rel_dat_mul)
}
```

## Take in data, and fit the multinomial regression model

``` r
mul_fit <- function(test_dat_mul_with_w_filename){
  ###
  # mul_fit() takes in the filename of cleaned data and return the fitted
  # multinomial regression model
  #
  # Input:
  #   test_dat_mul_with_w_filename: the filename of test_dat_mul_with_w that should
  #                                 contain everything useful for rel_dat_mul and
  #                                 test_dat_mul
  # Output:
  #   the fitted multinom object.
  ###

  rel_dat_mul <- get_rel_dat_mul(test_dat_mul_with_w_filename)
  print("model training...")
  mul_reg <- multinom(dest~. - w_total - areaId, data = rel_dat_mul, 
                        maxit = 5000)
  return(mul_reg)
}
```

## Take in fitted model, return probability coefficients

``` r
mul_predict <- function(coefficient_matrix, test_dat_mul_with_w_filename){
  ###
  # mul_predict() takes in a matrix of coefficients for multinomial regression
  # and return the predicted values as probabilities for rows in the test data
  #
  # Inputs:
  #   coefficient_matrix: a matrix of coefficients for multinomial regression,
  #                       which can be obtained via 
  #                         summary([multinom object])$coefficients or
  #                         the returned matrix from mul_sample() function
  #   test_dat_mul_with_w_filename: the filename of test_dat_mul_with_w that should
  #                                 contain everything useful for rel_dat_mul and
  #                                 test_dat_mul
  #
  # Output:
  #   a dataframe with predicted probabilities for each areaId
  ###
  
  # helper function for getting predicted values from multinomial
  # regression. The predicted values are probabilities.
  predict_mul <- function(coefs, dat){
    dat <- cbind(rep(1, nrow(dat)), dat)
    dat <- data.matrix(dat)
    rownames <- rownames(coefs)
    coefs <- matrix(coefs, ncol = ncol(coefs))
    init_pred <- coefs %*% t(dat)
    exp_pred <- rbind(rep(1, nrow(dat)), exp(init_pred))
    ones <- matrix(1, nrow = 1, ncol = nrow(exp_pred))
    denom <- ones %*% exp_pred
    exp_pred <- t(exp_pred)
    standardized_pred <- matrix(NA, nrow = nrow(exp_pred),
                                ncol = ncol(exp_pred))
    for(i in 1:ncol(exp_pred)){
      standardized_pred[,i] <- exp_pred[,i] / denom
    }
    
    colnames(standardized_pred) <- c("(Base)", rownames)
    return(standardized_pred)
  }

  test_data <- get_test_dat_mul(test_dat_mul_with_w_filename)
  print("model fitting...")
  test_data <- data.frame(test_data)
  dat <- test_data[,!(names(test_data) %in% c("areaId", "w_total"))]
  probs <- predict_mul(coefficient_matrix, dat)
  colnames(probs) <- c("p.ban", "p.lub", "p.mal", "p.mok", "p.ria", "p.ure", "p.off")
  result <- cbind(test_data[,names(test_data) %in% c("areaId")], probs)
  return(data.frame(result))
}
```

## Sample coefficients from multinomial regression model

``` r
mul_sample <- function(mul_reg){
  ###
  # mul_sample() takes in a multinom object and return a sample of coefficients
  #
  # Input:
  #   mul_reg: a multinom object
  #
  # Output:
  #   a matrix of coefficient sample who has the same format as the coefficient
  #   matrix contained in the summary([multinom object])
  ###
  
  summary_mul <- summary(mul_reg)
  coefs_tmp <- summary_mul$coefficients
  col_names <- colnames(coefs_tmp)
  coefs <- c()
  for(i in 1:nrow(coefs_tmp)){
    coefs <- c(coefs, coefs_tmp[i,])
  }
  sampled_mul <- mvrnorm(n = 1, 
                         mu = coefs,
                         Sigma = vcov(mul_reg))
  results <- matrix(as.matrix(sampled_mul), ncol = ncol(coefs_tmp), byrow = T)
  colnames(results) <- col_names
  rownames(results) <- rownames(coefs_tmp)
  return(results)
}
```
