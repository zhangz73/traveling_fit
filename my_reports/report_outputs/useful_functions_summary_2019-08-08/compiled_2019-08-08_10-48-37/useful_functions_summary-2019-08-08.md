Useful Functions Summary
================
Zhanhao Zhang
8/8/2019

## Take in data, and fit the multinomial regression model

``` r
mul_fit <- function(train_data_filename){
  ###
  # mul_fit() takes in the filename of cleaned data and return the fitted
  # multinomial regression model
  #
  # Input:
  #   train_data_filename: the filename of trained data, it must be in 
  #                        the ../data/clean directory
  #
  # Output:
  #   the fitted multinom object.
  ###
  train_data_dir <- paste("../data/clean/", train_data_filename, sep = "")
  train_data <- fread(train_data_dir, sep = ",", header = T)
  mul_reg <- multinom(dest~. - w_total - areaId, data = train_data, 
                        maxit = 5000)
  return(mul_reg)
}
```

## Take in fitted model, return probability coefficients

``` r
mul_predict <- function(coefficient_matrix, test_data_filename){
  ###
  # mul_predict() takes in a matrix of coefficients for multinomial regression
  # and return the predicted values as probabilities for rows in the test data
  #
  # Inputs:
  #   coefficient_matrix: a matrix of coefficients for multinomial regression,
  #                       which can be obtained via 
  #                         summary([multinom object])$coefficients or
  #                         the returned matrix from mul_sample() function
  #   test_data_filename: the filename for test data, must be in the
  #                       ../data/clean/ directory
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
  
  test_data_dir <- paste("../data/clean/", test_data_filename, sep = "")
  test_data <- fread(test_data_dir, sep = ",", header = T)
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
  results <- matrix(as.matrix(sampled_mul), ncol = ncol(coefs_tmp))
  colnames(results) <- col_names
  rownames(results) <- rownames(coefs_tmp)
  return(results)
}
```
