### Travel Model Fitting + Uncertainty ####
# Functions for performing draws from travel model fits
# 
# For now, we have multinomial model fit, based on population, distance, region of origin.
# 
###

# Load libraries ####
library(here)
library(data.table)
library(MASS)
library(nnet)

# Load and transform input data for prediction ####
get.data.for.prediction <- function(input.filename = "test_dat_mul_with_w.csv"){
  ###
  # get.data.for.prediction() returns the test_dat_mul to feed in multinomial
  # regression model for prediction
  # 
  # Input:
  #   input.filename: the filename of test_dat_mul_with_w that should
  #                   contain everything useful for rel_dat_mul and
  #                   test_dat_mul
  #                   default path = test_dat_mul_with_w.csv
  #
  # Output:
  #   The data to feed in multinomial regression model for prediction
  ###
  print("parsing data for prediction...")
  data <- fread(here("data/clean", input.filename), 
                               sep = ",", header = T)
  data <- data[,1:25]
  return(data)
}

get.input.data <- function(input.filename = "test_dat_mul_with_w.csv"){
  ###
  # get.input.data() returns the rel_dat_mul to feed in multinomial
  # regression model for model fitting
  # 
  # Input:
  #   input.filename: the filename of test_dat_mul_with_w that should
  #                   contain everything useful for rel_dat_mul and
  #                   test_dat_mul
  #                   default path = test_dat_mul_with_w.csv
  #
  # Output:
  #   The rel_dat_mul to feed in multinomial regression model for fitting
  ###
  print("parsing data for model fitting...")
  data.full <- fread(here("data/clean", input.filename), 
                               sep = ",", header = T)
  centrals <- c("to", "ti_ban", "ti_mal", "ti_lub", "ti_ria", "ti_mok", "ti_ure")
  input.data <- c()
  for(c in centrals){
    tmp <- data.full[data.full[[c]] > 0,]
    dat <- tmp[rep(seq(1, nrow(tmp)), tmp[[c]])]
    lst_dest <- rep(c, nrow(dat))
    dat <- cbind(data.frame(dat[,1:25]), data.frame(dest = lst_dest))
    input.data <- data.frame(rbind(data.frame(rel_dat_mul), data.frame(dat)))
  }
  return(input.data)
}

# Multinomial Model Fit ####
mul.fit <- function(input.data = NULL, filename = "test_dat_mul_with_w.csv"){
  ###
  # mul.fit() takes in the input.data from get.input.data() and returns
  #   a fit multinomial regression gravity model
  #
  # Input:
  #   input.data: data table formatted for a multinomial fit, 
  #               the output from from get.input.data() 
  #               
  #   filename:   If input.data == NULL, then we open and format the data at filename
  # Output: mul.reg
  #   the fitted multinomial model object
  ###
  if(is.null(input.data)){
    input.data <- get.input.data(input.filename = filename)
  }
  print("model training...")
  mul.reg <- multinom(dest~. - w_total - areaId, data = input.data, 
                      maxit = 5000)
  return(mul.reg)
}

# Predicted Values from Model Fit
predict.mul <- function(coefficient_matrix, data.pred = NULL, filename = "test_dat_mul_with_w.csv"){
  ###
  # mul_predict() takes in the matrix of coefficients from the multinomial regression 
  #  as well as the data set formatted for making predictions
  #  return the predicted values as probabilities for rows in the test data
  #
  # Inputs:
  #   coefficient_matrix: the matrix of coefficients for multinomial regression,
  #                       which can be obtained via coefficients()
  #
  #   data.pred: prediction data, the output from get.data.for.prediction()
  #   filname:   if data.pred == NULL, search for the data at this path location 
  #
  # Output:
  #   result: a data.table with predicted probabilities 
  #           of traveling from each areaId to each destination region
  ###
  
  # helper function for getting predicted values from multinomial ####
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
  
  # Load data from filename if none specified ####
  if (is.null(data.pred)){
    test_data <- get.data.for.prediction(filename)
  }
  print("model fitting...")
  test_data <- data.frame(test_data)
  dat <- test_data[,!(names(test_data) %in% c("areaId", "w_total"))]
  probs <- predict_mul(coefficient_matrix, dat)
  result <- cbind(test_data[,names(test_data) %in% c("areaId")], probs)
  colnames(result) <- c("areaId", "p.off", "p.ban", "p.lub", "p.mal", "p.mok", "p.ria", "p.ure")
  return(data.table(result))
}

# Perform Draw from Posterior Distribution ####
draw.sample.coefficients <- function(mul.reg){
  ###
  # draw.sample.coefficients() takes in a multinomial model object, such as the output from mul.fit(),
  #   performs a draw from the joint posterior distribution - assumed to be Gaussian -
  #   and returns a data table containing another coefficient matrix as in the output from predict.mul()
  #
  # Input:
  #   mul.reg: a multinomial model object
  #
  # Output:
  #   a matrix of model coefficients, with the same format is in the output of predict.mul()
  ###
  
  coefs_tmp <- coefficients(mul.reg)
  col_names <- colnames(coefs_tmp)
  coefs <- c()
  for(i in 1:nrow(coefs_tmp)){
    coefs <- c(coefs, coefs_tmp[i,])
  }
  sampled_mul <- mvrnorm(n = 1, 
                         mu = coefs,
                         Sigma = vcov(mul.reg))
  results <- matrix(as.matrix(sampled_mul), ncol = ncol(coefs_tmp), byrow = T)
  colnames(results) <- col_names
  rownames(results) <- rownames(coefs_tmp)
  return(results)
}

# Some test scripts ####
# # Read in data
# input.filename <- "test_dat_mul_with_w.csv"
# input.data <- get.input.data(input.filename = input.filename)
# # Fit data
# mul.model <- mul.fit(input.data = input.data)
# # predict data based off of fit model
# prediction.data <- get.data.for.prediction(input.filename = input.filename)
# mul.prediction <- predict.mul(coefficients(mul.model), prediction.data)
# # Perform draw
# mul.draw <- draw.sample.coefficients(mul.model)
# # Make predictions based off that draw => obtain probability distributions
# mul.draw.prediction <- predict.mul(mul.draw, prediction.data)