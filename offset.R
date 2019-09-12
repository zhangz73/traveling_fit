## Load function for extracting gravity models
source("extract_models.R")

### Load functions for correlation plots
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
  
  p_ratio <- ggplot(data = NULL, 
                    aes(x = x[z], y = y[z])) + geom_bin2d() + scale_x_log10() +
    scale_fill_gradient(low = "yellow", high = "red") +
    ggtitle(paste(title, y_name, "vs", x_name)) +
    xlab(x_name) + ylab(y_name) + 
    geom_line(data = NULL, aes(x=x[z], y=exp[z]), col="black") +
    scale_y_log10(limits = c(0.01, 100)) +
    theme(legend.position = "bottom")
  
  p <- ggplot() + 
    geom_point(data=NULL, 
               aes(x=x, y=res, col="red")) +
    ggtitle(paste(title, "residual plot" , x_name)) +
    xlab(x_name) + ylab("residual plot") + 
    geom_hline(yintercept = 0, col="black") +
    theme(legend.position = "none") +
    scale_x_log10()
  
  p_bin <- ggplot(data = NULL, aes(x = x, y = res)) + geom_bin2d() + scale_x_log10() +
    scale_fill_gradient(low = "yellow", high = "red") +
    ggtitle(paste(title, "heated residual plot" , x_name)) +
    xlab(x_name) + ylab("residual plot") + 
    geom_hline(yintercept = 0, col="black") +
    theme(legend.position = "bottom")
  
  return(list(ratio_plot = p2, heatmap_ratio = p_ratio,
              residual_plot = p, heatmap_res = p_bin))
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
  
  labels_names <- c("d", "N1", "N2")
  g1 <- ggarrange(p1$ratio_plot, p2$ratio_plot, p3$ratio_plot, ncol = 3, 
                  labels = labels_names)
  g2 <- ggarrange(p1$heatmap_ratio, p2$heatmap_ratio, p3$heatmap_ratio, 
                  ncol = 3, 
                  labels = labels_names)
  g3 <- ggarrange(p1$residual_plot, p2$residual_plot, p3$residual_plot, 
                  ncol = 3, 
                  labels = labels_names)
  g4 <- ggarrange(p1$heatmap_res, p2$heatmap_res, p3$heatmap_res, 
                  ncol = 3, 
                  labels = labels_names)
  g <- ggarrange(g1, g2, g3, g4, nrow = 4)
  
  png(filename = paste(dir, ".png", sep=""), width = 1000, height = 1300)
  print(g)
  dev.off()
}

## Load Data
rel_gravity_dat_covs <- read.csv("rel_gravity_dat_covs.csv")
BI_survey_data <- read.csv("BI_survey_data.csv")
n_sample <- data.frame(BI_survey_data[, c("areaId", "n")])
rel_gravity_dat_covs <- merge(rel_gravity_dat_covs, n_sample, by = "areaId")

## Extract models with/without offset on negative binomial cutoff at 20km
## Be sure to take out the rows where sample size = 0, otherwise weird
##  things can happen

smf_offset <- super_model_fits(rel_gravity_dat_covs[rel_gravity_dat_covs$n > 0,], 
                 "negbin", show_plots = F,
                 cutoff = 20 * 1000, use_time = F)
smf_nooffset <- super_model_fits(rel_gravity_dat_covs[rel_gravity_dat_covs$n > 0,], 
                               "negbin", show_plots = F,
                               cutoff = 20 * 1000, use_time = F, offset = F)

## Helper function to map the predicted flows to their corresponding
## covariates: d, N1, N2
get_total_w <- function(smf, rel_dat, cutoff){
  ###
  # Inputs:
  #   smf: returned object from the function super_model_fits
  #   rel_dat: same data fed into the super_model_fits
  #   cutoff: the cutoff value on distance for the super_model_fits
  #
  # Output:
  #   A dataframe with w, pred, d, N1 and N2, which can be passed into
  #   cor_plot_all() function to generate correlation plots
  ###
  
  dat_far <- rel_dat[rel_dat$d > cutoff,]
  dat_near <- rel_dat[rel_dat$d <= cutoff,]
  pred_far <- dat_far$N1 / dat_far$n * smf$model_far$fitted.values
  pred_near <- dat_near$N1 / dat_near$n * smf$model_near$fitted.values
  true_far <- dat_far$N1 / dat_far$n * smf$true_value_far
  true_near <- dat_near$N1 / dat_near$n * smf$true_value_near
  
  d <- c(dat_far$d, dat_near$d)
  N1 <- c(dat_far$N1, dat_near$N1)
  N2 <- c(dat_far$N2, dat_near$N2)
  
  return(data.frame(w = c(true_far, true_near), 
              pred = c(pred_far, pred_near), d = d, N1 = N1, N2 = N2))
} 

## Map the predicted commuting flows to their corresponding covariates
## and generate the correlation plots to visualize model performance
res_offset <- get_total_w(smf_offset, 
                          rel_gravity_dat_covs[rel_gravity_dat_covs$n > 0,], 
                          20 * 1000)
res_nooffset <- get_total_w(smf_nooffset, 
                            rel_gravity_dat_covs[rel_gravity_dat_covs$n > 0,], 
                            20 * 1000)
ssr1 <- sum((res_offset$w - res_offset$pred)^2)
ssr2 <- sum((res_nooffset$w - res_nooffset$pred)^2)

cor_plot_all("offset_negbin", res_offset, "./offset_negbin_heated")
cor_plot_all("nooffset_negbin", res_nooffset, "./nooffset_negbin_heated")


## Helper function to calculate the relative probabilities of 
## traveling to each of ad2 units for given areaIds as starting locations.
get_relative_probs <- function(res){
  ###
  # Input:
  #   res: the returned result from get_total_w()
  #
  # Output:
  #   return a dataframe containing probabilities of traveling to each
  #   ad2 unit from areaIds.
  ###
  
  pops <- BI_survey_data[, c("areaId", "pop", "ad2")]
  colnames(pops) <- c("areaId", "N1", "ad2")
  pops2 <- rel_gravity_dat_covs[, c("dest", "N2")]
  
  combined_df <- merge(merge(res, pops, by = "N1"), pops2, by = "N2")
  combined_df <- data.table(combined_df)
  
  result <- combined_df[, .(true_prob = sum(w), 
                            predict_prob = sum(pred),
                            N1 = mean(N1), N2 = mean(N2), d = mean(d)), 
                        by = .(areaId, ad2, dest)]
  total_w <- combined_df[, .(sw = sum(w)), by = areaId]
  total_pw <- combined_df[, .(pw = sum(pred)), by = areaId]
  
  result <- merge(merge(result, total_w, by = "areaId"), total_pw, by = "areaId")
  
  result <- data.frame(result)
  result$true_prob <- result$true_prob / result$sw
  result$predict_prob <- result$predict_prob / result$pw
  result <- result[result$sw > 0,]
  result <- result[, !colnames(result) %in% c("sw", "pw")]
  
  return(result)
}

## Calculate the relative probabilities for models with / without offset
## and generate correlation plots to visualize their performance
rel_probs_offset <- get_relative_probs(res_offset)
rel_probs_nooffset <- get_relative_probs(res_nooffset)
plot(rel_probs_offset$true_prob, rel_probs_offset$predict_prob, 
     xlim = c(0, 1), ylim = c(0, 1), xlab = "true probability",
     ylab = "predicted probability")
abline(a = 0, b = 1, col = "red")

colnames(rel_probs_offset) <- c("areaId", "ad2", "dest", "w", "pred", "N1", "N2", "d")
cor_plot_all("probability_offset", rel_probs_offset, "./probability_offset")
colnames(rel_probs_nooffset) <- c("areaId", "ad2", "dest", "w", "pred", "N1", "N2", "d")
cor_plot_all("probability_nooffset", rel_probs_nooffset, "./probability_nooffset")


