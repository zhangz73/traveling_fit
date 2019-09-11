source("extract_models.R")

rel_gravity_dat_covs <- read.csv("rel_gravity_dat_covs.csv")
BI_survey_data <- read.csv("BI_survey_data.csv")
n_sample <- data.frame(BI_survey_data[, c("areaId", "n")])
rel_gravity_dat_covs <- merge(rel_gravity_dat_covs, n_sample, by = "areaId")

smf_offset <- super_model_fits(rel_gravity_dat_covs[rel_gravity_dat_covs$n > 0,], 
                 "negbin", image_name = "offset",
                 cutoff = 20 * 1000, use_time = F)
smf_nooffset <- super_model_fits(rel_gravity_dat_covs[rel_gravity_dat_covs$n > 0,], 
                               "negbin", image_name = "no_offset",
                               cutoff = 20 * 1000, use_time = F, offset = F)
sum((smf_offset$model_far$fitted.values - smf_offset$true_value_far)^2) +
  sum((smf_offset$model_near$fitted.values - smf_offset$true_value_near)^2)
sum((smf_nooffset$model_far$fitted.values - smf_nooffset$true_value_far)^2) +
  sum((smf_nooffset$model_near$fitted.values - smf_nooffset$true_value_near)^2)

get_total_w <- function(smf, rel_dat, cutoff){
  dat_far <- rel_dat[rel_dat$d > cutoff,]
  dat_near <- rel_dat[rel_dat$d <= cutoff,]
  pred_far <- dat_far$N1 / dat_far$n * smf$model_far$fitted.values
  pred_near <- dat_near$N1 / dat_near$n * smf$model_near$fitted.values
  true_far <- dat_far$N1 / dat_far$n * smf$true_value_far
  true_near <- dat_near$N1 / dat_near$n * smf$true_value_near
  return(list(true_values = c(true_far, true_near), 
              pred_values = c(pred_far, pred_near)))
}
res_offset <- get_total_w(smf_offset, 
                          rel_gravity_dat_covs[rel_gravity_dat_covs$n > 0,], 
                          20 * 1000)
res_nooffset <- get_total_w(smf_nooffset, 
                            rel_gravity_dat_covs[rel_gravity_dat_covs$n > 0,], 
                            20 * 1000)

