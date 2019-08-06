library(data.table)
library(MASS)
library(ggplot2)
library(reshape2)
library(stringr)
library(pscl)
library(nnet)
library(ggpubr)

# Correlation plots
cor_plot_single <- function(x, y, res, exp, x_name, y_name, title){
  ###
  # cor_plot_single() is a function for plotting a ratio plot
  #                   and a residual plot for a model with respect
  #                   to a given covariate
  # Inputs:
  #   x: the given covariate to be plotted against
  #   y: the y-axis, in most cases, it is w^{(D)} / w^{(M)}
  #   res: the residual, which is w^{(D)} - w^{(M)}
  #   exp: the expected value of the ratio plot, in most cases,
  #        it is a horizontal line y = 1
  #   x_name: the name for x-axis
  #   y_name: the name for y-axis
  #   title: the name of the title of the plots
  # Outputs:
  #   A list of ggplot() objects, containing one ratio plot and one
  #   residual plot
  ###
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  z <- which(y != 0)
  p2 <- ggplot() +
    geom_point(data = NULL, 
               aes(x = x[z], y = y[z]), col = "red") +
    ggtitle(paste(title, y_name, "vs", x_name)) +
    xlab(x_name) + ylab(y_name) + 
    geom_line(data = NULL, aes(x=x[z], y=exp[z]), col="black") +
    scale_y_log10(limits = c(0.01, 100)) + 
    theme(legend.position = "none")
  p2 <- p2 + scale_x_log10()
  
  p <- ggplot() + 
    geom_point(data=NULL, 
               aes(x=x, y=res, col="red")) +
    ggtitle(paste(title, "residual plot" , x_name)) +
    xlab(x_name) + ylab("residual plot") + 
    geom_hline(yintercept = 0, col="black") +
    theme(legend.position = "none") +
    scale_x_log10()
  
  return(list(ratio_plot = p2, residual_plot = p))
}

cor_plot_all <- function(title, dat, dir){
  ###
  # cor_plot_all() is a function to generate all plots for a model
  # 4 ratio plots and 4 residual plots, arranged in a 2x4 matrix
  # Inputs:
  #   title: describes the name and features of the model, which is
  #          to be used in the title of the plots and the filename
  #          of the generated image
  #   dat: a dataframe in the form of (w, d, N1, N2, t, pred), an
  #        auxilary data structure for the plotting
  ###
  p1 <- cor_plot_single(dat$d, dat$w / dat$pred, dat$w - dat$pred,
                        rep(1, length(dat$d)), "geographic distance",
                        "w^{(D)} / w^{(M)}", title)
  p2 <- cor_plot_single(dat$N1, dat$w / dat$pred, dat$w - dat$pred,
                        rep(1, length(dat$d)), "source population",
                        "w^{(D)} / w^{(M)}", title)
  p3 <- cor_plot_single(dat$N2, dat$w / dat$pred, dat$w - dat$pred,
                        rep(1, length(dat$d)), "destination population",
                        "w^{(D)} / w^{(M)}", title)
  p4 <- cor_plot_single(dat$t, dat$w / dat$pred, dat$w - dat$pred,
                        rep(1, length(dat$d)), "travel time",
                        "w^{(D)} / w^{(M)}", title)
  g1 <- ggarrange(p1$ratio_plot, p2$ratio_plot, p3$ratio_plot,    
                  p4$ratio_plot, ncol = 4, 
                  labels = c("A", "B", "C", "D"))
  g2 <- ggarrange(p1$residual_plot, p2$residual_plot, p3$residual_plot, 
                  p4$residual_plot, ncol = 4, 
                  labels = c("A", "B", "C", "D"))
  g <- ggarrange(g1, g2, nrow = 2)
  
  png(filename = paste(dir, ".png", sep=""), width = 1200, height = 600)
  print(g)
  dev.off()
}

super_model_fits <- function(rel_dat, model_name, image_name = NULL,
                             show_plots = TRUE, cutoff = 0, 
                             base_dir = "./", use_time = TRUE,
                             ad2_indicators = NULL,
                             use_distance_to_mal = FALSE,
                             exp_dependence_far = FALSE,
                             exp_dependence_near = TRUE){
  ###
  # super_model_fits() is a powerful function that can generate
  # a model based on users choice.
  # Inputs:
  #   rel_dat: the dataframe that contains data for all useful
  #            covariates and response variable that will be used
  #            by the model.
  #   model_name: the name indicating type of a model. Available
  #               choices are: "negbin" for negative binomial regression,
  #               "poisson" for poisson regression,
  #               "linreg" for linear regression,
  #               "zinf" for zero-inflated negative binomial regression
  #   image_name: the name of the image file for ratio plots and
  #               residual plots. image_name will also be used as 
  #               a title in the generated plots. 
  #               Do not need to specify if you do
  #               not wish to generate the plots.
  #   show_plots: whether to create and generate ratio plots and 
  #               residual plots. DEFAULT = TRUE.
  #   cutoff: the cutoff value for geographic distance or traveling 
  #           time. For traveling time, cutoff=50 for 50 minutes;
  #           for geographic distance, cutoff = 20*1000 for 20km
  #           DEFAULT = 0 (no cutoff).
  #   base_dir: the base directory to place the geneated plots.
  #             DEFAULT = "./", the current directory.
  #             Note: please use RELATIVE DIRECTORY!
  #   use_time: whether use traveling time for the distance covariates.
  #             TRUE for using traveling time, FALSE for using geographic
  #             distance. DEFAULT = TRUE.
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
  
  get_formula <- function(use_time, ad2_indicators,
                          use_distance_to_mal, exp_dependence){
    formula <- NULL
    labels <- c("log(N1)", "log(N2)")
    if(use_time){
      if(exp_dependence){
        labels <- c(labels, "t")
      } else{
        labels <- c(labels, "log(t)")
      }
    } else{
      if(exp_dependence){
        labels <- c(labels, "d")
      } else{
        labels <- c(labels, "log(d)")
      }
    }
    if(use_distance_to_mal){
      if(use_time){
        labels <- c(labels, "t_mal")
      } else{
        labels <- c(labels, "d_mal")
      }
    }
    if(!is.null(ad2_indicators)){
      labels <- c(labels, ad2_indicators)
    }
    formula <- reformulate(response = "w", termlabels = labels)
    return(formula)
  }
  
  get_model <- function(model_name, formula, data){
    model <- NULL
    if(model_name == "negbin"){
      model <- glm.nb(formula = formula, data = data)
    } else if(model_name == "poisson"){
      model <- glm(formula = formula, data = data, family = "poisson")
    } else if(model_name == "linreg"){
      model <- lm(formula = formula, data = data)
    } else if(model_name == "zinf"){
      model <- zeroinfl(formula = formula, data = data, dist = "negbin")
    } else{
      print("pick a model name from: negbin, poisson, linreg, zinf")
      return(NULL)
    }
    return(model)
  }
  
  rel_dat_copy <- rel_dat
  for(i in 1:ncol(rel_dat_copy)){
    rel_dat_copy[,i] <- unlist(rel_dat_copy[,i])
  }
  rel_gravity_far <- NULL
  rel_gravity_near <- NULL
  
  if(use_time){
    rel_gravity_far <- rel_dat[rel_dat_copy$t > cutoff,]
    rel_gravity_near <- rel_dat[rel_dat_copy$t <= cutoff,]
  } else{
    rel_gravity_far <- rel_dat[rel_dat$d > cutoff,]
    rel_gravity_near <- rel_dat[rel_dat$d <= cutoff,] 
  }
  
  formula_far <- get_formula(use_time, ad2_indicators,
                             use_distance_to_mal, exp_dependence_far)
  formula_near <- get_formula(use_time, ad2_indicators,
                              use_distance_to_mal, exp_dependence_near)
  
  model_far <- get_model(model_name, formula_far, rel_gravity_far)
  if(cutoff > 0){
    model_near <- get_model(model_name, formula_near, rel_gravity_near)
  }
  dir <- paste(base_dir, image_name, sep="")
  
  if(cutoff <= 0){
    
    dat <- cbind(rel_gravity_far, list(pred=model_far$fitted.values))
    for(i in 1:ncol(dat)){
      dat[,i] <- unlist(dat[,i])
    }
    if(show_plots){
      cor_plot_all(paste(image_name, "\n", sep=""), dat, dir)
    }
    df <- AIC(model_far) - AIC(model_far, k = 0)
    return(list(model_far=model_far, model_near=list(fitted.values=0, 
                                                     aic=df, loglik=0), 
                true_value_far=unlist(rel_gravity_far$w), 
                true_value_near=0))
  }
  
  dat_all <- cbind(rbind(rel_gravity_far, rel_gravity_near),
                   list(pred=c(model_far$fitted.values,
                               model_near$fitted.values)))
  
  for(i in 1:ncol(dat_all)){
    dat_all[,i] <- unlist(dat_all[,i])
  }
  if(show_plots){
    cor_plot_all(paste(image_name, "\n", sep=""), dat_all, dir)
  }
  return(list(model_far=model_far, model_near=model_near, 
              true_value_far=unlist(rel_gravity_far$w),
              true_value_near=unlist(rel_gravity_near$w)))
}