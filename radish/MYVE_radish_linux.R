
library(corMLPE)
library(radish)
library(terra)
library(raster)
library(sf)
library(sp)

# Read in raster stack from disk, shouldn't need individual layers above
covariates <- stack("D:/Github/MYVE_LGA/data/covariates_stack_rescaled.grd")
names(covariates)

################ Spatial coordinates ###############
myve_coords <- as.data.frame(read.csv("~/Documents/radish/data/MYVECoordinates_nobaja.csv"))
myve_points <- SpatialPoints(myve_coords[,2:3])

################# Genetic distance convert to matrix
gendist_raw <- read.csv("~/Documents/radish/data/GeneticDist.csv")
gendist_raw$GenDisPC12 <- (gendist_raw$GenDisPC12 - min(gendist_raw$GenDisPC12, na.rm=TRUE)) / (max(gendist_raw$GenDisPC12,na.rm=TRUE - min(gendist_raw$GenDisPC12, na.rm=TRUE)))
ids <- sort(unique(c(gendist_raw$Indv1, gendist_raw$Indv2)))
gen_mat <- matrix(0, nrow = length(ids), ncol = length(ids),  dimnames = list(ids, ids))
gen_mat[cbind(gendist_raw$Indv1, gendist_raw$Indv2)] <- gendist_raw$GenDisPC12
gen_mat[cbind(gendist_raw$Indv2, gendist_raw$Indv1)] <- gendist_raw$GenDisPC12

# Reorder to match ordering of myve points
gen_mat <- gen_mat[myve_coords$Seq_ID, myve_coords$Seq_ID]
# Double check that ordering is correct
dimnames(gen_mat)[[1]] == myve_coords$Seq_ID

surface <- conductance_surface(covariates, myve_points, directions = 8)

# Loop through mlpe version of all univariate models, compare with aic
uni_models <- list()
for(idx in 1:length(surface$covariates)){
  uni_models[[idx]] <- radish(reformulate(surface$covariates[idx], response = "gen_mat"), surface, 
                      radish::loglinear_conductance, radish::mlpe, nu=105)
}

mlpe_models[2]$aic
model_aics <- vector()
model_lglk <- vector()
model_gradnorm <- vector()
model_coef <- vector()
model_ste <- vector()
model_zval <- vector()
model_pval <- vector()
model_alpha <- vector()
model_beta <- vector()
model_tau <- vector()
model_rho <- vector()
model_iters <- vector()

for(idx in 1:length(uni_models)){
  model_aics[idx] <- summary(uni_models[[idx]])$aic
  model_lglk[idx] <- summary(uni_models[[idx]])$loglik
  model_gradnorm[idx] <- summary(uni_models[[idx]])$gradnorm
  model_coef[idx] <- summary(uni_models[[idx]])$ztable[1]
  model_ste[idx] <- summary(uni_models[[idx]])$ztable[2]
  model_zval[idx] <- summary(uni_models[[idx]])$ztable[3]
  model_pval[idx] <- summary(uni_models[[idx]])$ztable[4]
  model_alpha[idx] <- summary(uni_models[[idx]])$phi[1]
  model_beta[idx] <- summary(uni_models[[idx]])$phi[2]
  model_tau[idx] <- summary(uni_models[[idx]])$phi[3]
  model_rho[idx] <- summary(uni_models[[idx]])$phi[4]
  model_iters[idx] <- summary(uni_models[[idx]])$iters
}


models_out <- data.frame(vars = names(covariates), aic = model_aics, loglike = model_lglk, gradnorm = model_gradnorm,
                         coefficient = model_coef, sterror = model_ste, zval = model_zval, pval = model_pval,
                         alpha = model_alpha, beta = model_beta, tau = model_tau, rho = model_rho, iterations = model_iters)
write.csv(models_out, "~/Documents/radish/Radish_univariate_models_out.csv")


# Read in raster stack from disk, shouldn't need individual layers above
#covariates_round2 <- stack("~/Documents/radish/data/covariates_stack_round2.grd")
# Now taking best model (TRI) and running each two-variable model
#twovar_models <- list()
#surface_round2 <- conductance_surface(covariates_round2, myve_points, directions = 8)
#covariates_round2 <- stack(list(covariates[[3]],covariates[[1]],covariates[[4]],covariates[[5]],
#                           covariates[[6]],covariates[[7]],covariates[[9]]))

#for(idx in 1:(length(surface_round2$covariates)-1)){
#  twovar_models[[idx]] <- radish(reformulate(termlabels = c(surface_round2$covariates[1],surface_round2$covariates[idx+1]), response = "gen_mat"), surface_round2, 
#                              radish::loglinear_conductance, radish::mlpe)
#  }


#for(idx in 1:length(twovar_models)){
#  model_aics[idx] <- summary(twovar_models[[idx]])$aic
#  model_lglk[idx] <- summary(twovar_models[[idx]])$loglik
#  model_gradnorm[idx] <- summary(twovar_models[[idx]])$gradnorm
#  model_coef[idx] <- summary(twovar_models[[idx]])$ztable[1]
#  model_ste[idx] <- summary(twovar_models[[idx]])$ztable[2]
#  model_zval[idx] <- summary(twovar_models[[idx]])$ztable[3]
#  model_pval[idx] <- summary(twovar_models[[idx]])$ztable[4]
#  model_alpha[idx] <- summary(twovar_models[[idx]])$phi[1]
#  model_beta[idx] <- summary(twovar_models[[idx]])$phi[2]
 # model_tau[idx] <- summary(twovar_models[[idx]])$phi[3]
#  model_rho[idx] <- summary(twovar_models[[idx]])$phi[4]
 # model_iters[idx] <- summary(twovar_models[[idx]])$iters
#}

#models_out <- data.frame(aic = model_aics, loglike = model_lglk, gradnorm = model_gradnorm,
#                         coefficient = model_coef, sterror = model_ste, zval = model_zval, pval = model_pval,
#                         alpha = model_alpha, beta = model_beta, tau = model_tau, rho = model_rho, iterations = model_iters)
#write.csv(models_out, "D:/UM/TXbats/LGA/Radish_twovariate_models_out.csv")


