predict_allType <- function(coefs, dat, model_name){
  predict_lin <- function(coefs, dat){
    coefs <- matrix(coefs, nrow=1)
    pred <- coefs %*% t(dat)
    return(pred)
  }
  
  predict_poi <- function(coefs, dat){
    coefs <- matrix(coefs, nrow=1)
    pred <- coefs %*% t(dat)
    pred <- exp(pred)
    return(pred)
  }
  
  predict_mul <- function(coefs, dat){
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
  
  dat <- cbind(rep(1, nrow(dat)), dat)
  dat <- data.matrix(dat)
  if(model_name == "linreg"){
    return(predict_lin(coefs, dat))
  } else if(model_name == "poisson" || model_name == "negbin"){
    return(predict_poi(coefs, dat))
  } else if(model_name == "multinom"){
    return(predict_mul(coefs, dat))
  } else{
    print('Usage: predict_allType(coefs, dat, 
          model_name=c("linreg", "poisson", "negbin", "multinom"))')
    return(NULL)
  }
}

sample_coefs <- function(est, sd){
  sample <- matrix(NA, nrow = nrow(est), ncol = ncol(est))
  for(i in 1:nrow(est)){
    for(j in 1:ncol(est)){
      sample[i, j] <- rnorm(1, est[i, j], sd[i, j])
    }
  }
  rownames(sample) <- rownames(est)
  colnames(sample) <- colnames(est)
  return(sample)
}

predict_allType(nb$coefficients, small_df, "negbin")
predict_allType(s$coefficients, test_dat_mul[1:3,1:ncol(r)-1], "multinom")
