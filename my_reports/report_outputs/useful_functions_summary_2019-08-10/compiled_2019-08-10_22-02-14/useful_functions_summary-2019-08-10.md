Useful Functions Summary
================
Zhanhao Zhang
8/8/2019

## Script For Model Extraction

``` r
source(here("scripts", "extract_models.R"))
```

## Data Transformation Helper Function

``` r
get_test_dat_mul <- function(comprehensive_data_filename){
  ###
  # get_test_dat_mul() returns the test_dat_mul to feed in multinomial
  # regression model for prediction
  # 
  # Input:
  #   comprehensive_data_filename: the filename of comprehensive_data that should
  #                                 contain everything useful for rel_dat_mul and
  #                                 test_dat_mul
  #
  # Output:
  #   The test_dat_mul to feed in multinomial regression model for prediction
  ###
  print("parsing data for prediction...")
  comprehensive_data <- fread(here("data/clean", comprehensive_data_filename), 
                               sep = ",", header = T)
  comprehensive_data <- data.frame(comprehensive_data)
  test_dat_mul <- comprehensive_data[,1:25]
  return(test_dat_mul)
}

get_rel_dat_mul <- function(comprehensive_data_filename){
  ###
  # get_rel_dat_mul() returns the rel_dat_mul to feed in multinomial
  # regression model for model fitting
  # 
  # Input:
  #   comprehensive_data_filename: the filename of comprehensive_data that should
  #                                 contain everything useful for rel_dat_mul and
  #                                 test_dat_mul
  #
  # Output:
  #   The rel_dat_mul to feed in multinomial regression model for fitting
  ###
  
  print("parsing data for model fitting...")
  comprehensive_data <- fread(here("data/clean", comprehensive_data_filename), 
                               sep = ",", header = T)
  comprehensive_data <- data.frame(comprehensive_data)
  centrals <- c("to", "ti_ban", "ti_mal", "ti_lub", "ti_ria", "ti_mok", "ti_ure")
  test_dat_mul <- comprehensive_data[,1:25]
  test_dat_mul <- data.frame(test_dat_mul)
  rel_dat_mul <- c()
  comprehensive_data <- data.table(comprehensive_data)
  for(c in centrals){
    tmp <- comprehensive_data[comprehensive_data[[c]] > 0,]
    dat <- tmp[rep(seq(1, nrow(tmp)), tmp[[c]])]
    lst_dest <- rep(c, nrow(dat))
    dat <- cbind(data.frame(dat[,1:25]), data.frame(dest = lst_dest))
    rel_dat_mul <- data.frame(rbind(data.frame(rel_dat_mul), data.frame(dat)))
  }
  return(rel_dat_mul)
}

get_rel_gravity_dat_covs <- function(comprehensive_data_filename){
  ###
  # get_rel_gravity_dat_covs() is a function that retrieve useful data
  # for poisson regression and negative binomial regression
  #
  # Input:
  #   comprehensive_data_filename: the file that contains data for all useful
  #            covariates and response variable that will be used
  #            by the model.
  #
  # Output:
  #   return the dataframe that contains all relevant data to feed into
  #   poisson regression or negative binomial regression models.
  ###
  
  comprehensive_data <- fread(here("data/clean", comprehensive_data_filename), 
                               sep = ",", header = T)
  centrals <- c("to", "ti_ban", "ti_mal", "ti_lub", "ti_ria", "ti_mok", "ti_ure")
  d_names <- paste("d.", centrals, sep = "")
  t_names <- paste("t.", centrals, sep = "")
  N2_names <- paste("N2.", centrals, sep = "")
  
  all_names <- c(N2_names, d_names, t_names, centrals, "w_total")
  
  comprehensive_data <- data.frame(comprehensive_data)
  sub_df <- comprehensive_data[, !(names(comprehensive_data) %in% all_names)]
  sub_df <- sub_df[rep((1:nrow(comprehensive_data)), each=length(centrals)),]
  
  melt_helper <- function(names){
    sub <- melt(comprehensive_data[, (names(comprehensive_data) %in% names)])
    val <- sub$value
    variable <- sub$variable
    
    var_m <- matrix(variable, nrow = length(names), byrow = T)
    var_m <- as.character(matrix(var_m, nrow = 1))
    
    val_m <- matrix(val, nrow = length(names), byrow = T)
    val_m <- as.vector(matrix(val_m, nrow = 1))
    
    return(data.frame(variable = var_m, value = val_m))
  }
  
  sub_d <- melt_helper(d_names)$value
  sub_t <- melt_helper(t_names)$value
  sub_N2 <- melt_helper(N2_names)$value
  sub_w <- melt_helper(centrals)
  
  rel_gravity_dat_covs <- cbind(sub_df, data.frame(d = sub_d, t = sub_t, N2 = sub_N2,
                              w = sub_w$value, dest = sub_w$variable))
  rownames(rel_gravity_dat_covs) <- NULL
  return(rel_gravity_dat_covs)
}
```

## Take in data, and fit the multinomial regression model

``` r
mul_fit <- function(comprehensive_data_filename){
  ###
  # mul_fit() takes in the filename of cleaned data and return the fitted
  # multinomial regression model
  #
  # Input:
  #   comprehensive_data_filename: the filename of comprehensive_data that should
  #                                 contain everything useful for rel_dat_mul and
  #                                 test_dat_mul
  # Output:
  #   the fitted multinom object.
  ###

  rel_dat_mul <- get_rel_dat_mul(comprehensive_data_filename)
  print("model training...")
  mul_reg <- multinom(dest~. - w_total - areaId, data = rel_dat_mul, 
                        maxit = 5000)
  return(mul_reg)
}

poi_fit <- function(comprehensive_data_filename, 
                             model_name, cutoff = 0, use_time = F,
                             ad2_indicators = NULL,
                             use_distance_to_mal = FALSE,
                             exp_dependence_far = FALSE,
                             exp_dependence_near = TRUE){
  ###
  # poi_fit() is a powerful function that can generate
  # a model based on users choice.
  # Inputs:
  #   comprehensive_data_filename: the file that contains data for all useful
  #            covariates and response variable that will be used
  #            by the model.
  #   model_name: the name indicating type of a model. Available
  #               choices are: "negbin" for negative binomial regression,
  #               "poisson" for poisson regression,
  #               "linreg" for linear regression,
  #               "zinf" for zero-inflated negative binomial regression
  #   cutoff: the cutoff value for geographic distance or traveling 
  #           time. For traveling time, cutoff=50 for 50 minutes;
  #           for geographic distance, cutoff = 20*1000 for 20km
  #           DEFAULT = 0 (no cutoff).
  #   use_time: whether use traveling time for the distance covariates.
  #             TRUE for using traveling time, FALSE for using geographic
  #             distance. DEFAULT = FALSE
  #   ad2_indicators: a vector storing the names of ad2 indicators.
  #                   For instance, c("Peri", "Malabo", "Baney").
  #                   The name of these indicators should be consistent
  #                   with their corresponding column names.
  #                   DEFAULT = NULL (no ad2 indicators).
  #   use_distance_to_mal: whether to consider the distance from source
  #                        to Malabo for model fitting. Will use
  #                        traveling time to Malabo if use_time = TRUE,
  #                        geographic distance to Malabo othewise.
  #                        DEFAULT = FALSE.
  #   exp_dependence_far: whether to assume a exponential dependence on
  #                       distance on trips above the cutoff. TRUE for
  #                       exponential dependence, FALSE for power
  #                       dependence. DEFAULT = FALSE.
  #   exp_dependence_near: whether to assume a exponential dependence on
  #                       distance on trips below the cutoff. TRUE for
  #                       exponential dependence, FALSE for power
  #                       dependence. DEFAULT = TRUE
  #
  # Outputs:
  #   a list containing four variables:
  #     model_far: fitted model for trips above the cutoff
  #     model_near: fitted model for trips below the cutoff
  #     true_value_far: true values of the responsive variable for
  #                     trips above the cutoff
  #     true_value_near: true values of the responsive variable for
  #                     trips below the cutoff
  #   Note: for the case where cutoff=0 (no cutoff), the fitted
  #         model is stored in model_far, and the original values
  #         of the response variable are stored in true_value_far
  #         model_near and true_value_near will be substituted with
  #         a few basic values that won't interrupt the workflows
  #         for calculating total AIC and sum of residual squares.
  ###
  
  rel_gravity_dat_covs <- get_rel_gravity_dat_covs(comprehensive_data_filename)
  return(super_model_fits(rel_dat = rel_gravity_dat_covs, model_name = model_name, 
                          show_plots = F, use_time = use_time,
                          ad2_indicators = ad2_indicators,
                          use_distance_to_mal = use_distance_to_mal,
                          exp_dependence_far = exp_dependence_far,
                          exp_dependence_near = exp_dependence_near))
}
```

## Take in fitted model, return probability coefficients

``` r
mul_predict <- function(coefficient_matrix, comprehensive_data_filename){
  ###
  # mul_predict() takes in a matrix of coefficients for multinomial regression
  # and return the predicted values as probabilities for rows in the test data
  #
  # Inputs:
  #   coefficient_matrix: a matrix of coefficients for multinomial regression,
  #                       which can be obtained via 
  #                         summary([multinom object])$coefficients or
  #                         the returned matrix from mul_sample() function
  #   comprehensive_data_filename: the filename of comprehensive_data that should
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

  test_data <- get_test_dat_mul(comprehensive_data_filename)
  names <- unique(get_rel_dat_mul(comprehensive_data_filename)$dest)
  
  print("model fitting...")
  test_data <- data.frame(test_data)
  dat <- test_data[,!(names(test_data) %in% c("areaId", "w_total"))]
  probs <- predict_mul(coefficient_matrix, dat)
  names <- str_replace(names, "ti_", "p.")
  colnames(probs) <- str_replace(names, "to", "p.off")
  result <- cbind(test_data[,names(test_data) %in% c("areaId")], probs)
  return(data.frame(result))
}

poi_predict <- function(coefs, original_model, dat){
  ###
  # poi_predict() can takes the given coefficients and return the predicted
  # value for the given dataset, for poisson regression and negative binomial
  # regression
  #
  # Inputs:
  #   coefs: the given coefficients for poisson or negative binomial regression
  #   original_model: the glm() object of the original model returned by the
  #                   poi_fit()
  #   dat: the dataframe that contains data for all useful
  #            covariates and response variable that will be used
  #            by the model.
  #
  # Output:
  #   return the predicted values from a poisson or negative binomial regression
  #   model.
  #
  # Note:
  #   if you are using cutoff, please be sure that you only pass in the dataframe
  #   relevant to the current model. I.e, if cutoff = 20 * 1000 and you want to
  #   predict the values for distance_near, then please make sure the parameter
  #   dat doesn't contain any values with cutoff > 20 * 1000
  ###
  
  m <- original_model
  m$coefficients <- coefs
  return(predict(m, dat, type = "response"))
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
  
  coefs_tmp <- coefficients(mul_reg)
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

poi_sample <- function(model){
  ###
  # poi_sample() is a helper function for sampling the
  # coefficients of a model while taking the correlations
  # between each covariate into account. The sampling takes
  # place in a multivariate normal distribution.
  #
  # Input:
  #   model: the object of a lm() or glm(). Objects from
  #          multinom() are NOT accepted in this function.
  #
  # Output:
  #   a matrix containing one samples of coefficients of the given model.
  ###
  coefs <- model$coefficients
  cov_matrix <- vcov(model)
  return(mvrnorm(n = 1, mu = coefs, Sigma = cov_matrix))
}
```
