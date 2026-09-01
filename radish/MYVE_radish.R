pak::pak("nspope/corMLPE")
pak::pak("nspope/radish")
library(corMLPE)
library(radish)
library(terra)
library(raster)
library(sf)
library(sp)

# Had to do this and add a line about C++ 14 to get this to work. 
usethis::edit_r_makevars()


data(melip)

setwd("D:/Github/MYVE_LGA/data")

env_layers <- list.files("./SpatialLayers", pattern=".tif$")

# scaling spatial covariates helps avoid numeric overflow
#covariates <- raster::stack(list(altitude = raster::scale(melip.altitude), 
 #                                forestcover = raster::scale(melip.forestcover)))

# Create raster stack - scaled
covariates1_scaled <- raster::stack(list(crop = raster::scale(raster(paste0("./SpatialLayers/",env_layers[1]))),
                                 dist2crop = raster::scale(raster(paste0("./SpatialLayers/",env_layers[2]))),
                                 dist2karst = raster::scale(raster(paste0("./SpatialLayers/",env_layers[3]))),
                                 dist2water = raster::scale(raster(paste0("./SpatialLayers/",env_layers[4]))),
                                 elevation = raster::scale(raster(paste0("./SpatialLayers/",env_layers[5]))),
                                 enm = raster::scale(raster(paste0("./SpatialLayers/",env_layers[6]))),
                                 forest = raster::scale(raster(paste0("./SpatialLayers/",env_layers[7]))),
                                 forest_gcam = raster::scale(raster(paste0("./SpatialLayers/",env_layers[8]))),
                                 ghfp1 = raster::scale(raster(paste0("./SpatialLayers/",env_layers[9]))),
                                 ghpf2 = raster::scale(raster(paste0("./SpatialLayers/",env_layers[10]))),
                                 karst = raster::scale(raster(paste0("./SpatialLayers/",env_layers[11]))),
                                 lights = raster::scale(raster(paste0("./SpatialLayers/",env_layers[12]))),
                                 slope = raster::scale(raster(paste0("./SpatialLayers/",env_layers[13]))),
                                 TPI = raster::scale(raster(paste0("./SpatialLayers/",env_layers[14]))),
                                 TRI = raster::scale(raster(paste0("./SpatialLayers/",env_layers[15]))),
                                 water = raster::scale(raster(paste0("./SpatialLayers/",env_layers[16])))
                                  ))

# Create raster stack - unscaled
covariates1 <- raster::stack(list(crop = raster(paste0("./SpatialLayers/",env_layers[1])),
                                  dist2crop = raster(paste0("./SpatialLayers/",env_layers[2])),
                                  dist2karst = raster(paste0("./SpatialLayers/",env_layers[3])),
                                  dist2water = raster(paste0("./SpatialLayers/",env_layers[4])),
                                  elevation = raster(paste0("./SpatialLayers/",env_layers[5])),
                                  enm = raster(paste0("./SpatialLayers/",env_layers[6])),
                                  forest = raster(paste0("./SpatialLayers/",env_layers[7])),
                                  forest_gcam = raster(paste0("./SpatialLayers/",env_layers[8])),
                                  ghfp1 = raster(paste0("./SpatialLayers/",env_layers[9])),
                                  ghpf2 = raster(paste0("./SpatialLayers/",env_layers[10])),
                                  karst = raster(paste0("./SpatialLayers/",env_layers[11])),
                                  lights = raster(paste0("./SpatialLayers/",env_layers[12])),
                                  slope = raster(paste0("./SpatialLayers/",env_layers[13])),
                                  TPI = raster(paste0("./SpatialLayers/",env_layers[14])),
                                  TRI = raster(paste0("./SpatialLayers/",env_layers[15])),
                                  water = raster(paste0("./SpatialLayers/",env_layers[16]))
))
# Have identical NoData mask on all layers
na_mask <- sum(covariates1_scaled)
covariates1 <- raster::mask(covariates1, na_mask)

# Write raster stack to disk
writeRaster(covariates1_scaled, filename = "covariates_stack_rescaled.grd", overwrite=TRUE)
writeRaster(covariates1, filename = "covariates_stack_raw.grd", overwrite=TRUE)
# Read in raster stack from disk, shouldn't need individual layers above
covariates <- stack("D:/Github/MYVE_LGA/data/covariates_stack_rescaled.grd")
names(covariates)

################ Spatial coordinates ###############
myve_coords <- as.data.frame(read.csv("D:/Github/MYVE_LGA/data/MYVECoordinates_nobaja.csv"))
head(myve_coords)
myve_points <- SpatialPoints(myve_coords[,2:3])

plot(covariates[["enm"]])
points(myve_points, pch = 19)

################# Genetic distance convert to matrix
gendist_raw <- read.csv("D:/Github/MYVE_LGA/data/GeneticDist.csv")
gendist_raw$GenDisPC12 <- (gendist_raw$GenDisPC12 - min(gendist_raw$GenDisPC12, na.rm=TRUE)) / (max(gendist_raw$GenDisPC12,na.rm=TRUE - min(gendist_raw$GenDisPC12, na.rm=TRUE)))
ids <- sort(unique(c(gendist_raw$Indv1, gendist_raw$Indv2)))
gen_mat <- matrix(0, nrow = length(ids), ncol = length(ids),  dimnames = list(ids, ids))
gen_mat[cbind(gendist_raw$Indv1, gendist_raw$Indv2)] <- gendist_raw$GenDisPC12
gen_mat[cbind(gendist_raw$Indv2, gendist_raw$Indv1)] <- gendist_raw$GenDisPC12

# Reorder to match ordering of myve points
gen_mat <- gen_mat[myve_coords$Seq_ID, myve_coords$Seq_ID]
# Double check that ordering is correct
dimnames(gen_mat)[[1]] == myve_coords$Seq_ID
dimnames(gen_mat) <- NULL
dimnames(gen_mat)[[1]][order(dimnames(gen_mat)[[1]])]

surface <- conductance_surface(covariates, myve_points, directions = 8)
surface_crop <- conductance_surface(stack(covariates[[1]]), myve_points, directions = 8)

fit_nnls_nate <- radish(melip.Fst ~ forestcover + altitude, surface, 
                   radish::loglinear_conductance, radish::leastsquares)
fit_nnls <- radish(gen_mat ~ enm, surface, 
                   radish::loglinear_conductance, radish::leastsquares)

summary(fit_nnls)

# refit with with a different measurement model that models
# dependence among pairwise measurements (radish::mlpe)
fit_mlpe_nate <- radish(melip.Fst ~ forestcover + altitude, surface, 
                   radish::loglinear_conductance, radish::mlpe)

fit_mlpe <- radish(gen_mat ~ enm, surface, 
                   radish::loglinear_conductance, radish::mlpe) # 5-10 minutes

summary(fit_mlpe)

# visualisation:
plot(fitted(fit_mlpe, "distance"), gen_mat, pch = 19,
     xlab = "Optimized resistance distance", ylab = "Genetic distance")

# visualise estimated conductance surface and asymptotic confidence intervals
fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)

plot(fitted_conductance[["est"]], 
     main = "Fitted conductance surface\n(enm)")
plot(fitted_conductance[["lower95"]], 
     main = "Fitted conductance surface\n(lower 95% CI)")
plot(fitted_conductance[["upper95"]], main = 
       "Fitted conductance surface\n(upper 95% CI)")

# visualise likelihood surface across grid (takes awhile)
theta <- as.matrix(expand.grid(forestcover=seq(-1,1,length.out=21), 
                               altitude=seq(-1,1,length.out=21)))
theta <- as.matrix(expand.grid(enm=seq(-1,1,length.out=21)))
grid <- radish_grid(theta, melip.Fst ~ forestcover + altitude, surface,
                    radish::loglinear_conductance, radish::mlpe)
grid <- radish_grid(theta, gen_mat ~ enm, surface,
                    radish::loglinear_conductance, radish::mlpe) #(5-10 minutes)

library(ggplot2)
ggplot(data.frame(loglik=grid$loglik, grid$theta), 
       aes(x=forestcover, y=altitude)) + 
  geom_tile(aes(fill=loglik)) + 
  geom_contour(aes(z=loglik), color="black") +
  annotate(geom = "point", colour = "red",
           x = coef(fit_mlpe)["forestcover"], 
           y = coef(fit_mlpe)["altitude"]) +
  theme_bw() +
  xlab(expression(theta[altitude])) +
  ylab(expression(theta[forestcover]))

ggplot(data.frame(loglik=grid$loglik, grid$theta), 
       aes(x=enm)) + 
  geom_tile(aes(fill=loglik)) + 
  geom_contour(aes(z=loglik), color="black") +
  theme_bw() +
  xlab(expression(theta[altitude])) +
  ylab(expression(theta[forestcover]))

# calculate resistance distances across grid
distances <- radish_distance(theta, ~forestcover + altitude, 
                             surface, radish::loglinear_conductance)

ibd <- which(theta[,1] == 0 & theta[,2] == 0)
plot(distances$distance[,,ibd], melip.Fst, pch = 19, 
     xlab = "Null resistance distance (IBD)", ylab = "Fst")

# model selection:
# fit a reduced model without "forestcover" covariate, and compare to 
# full model via a likelihood ratio test
fit_mlpe_reduced <- radish(melip.Fst ~ altitude, surface, 
                           radish::loglinear_conductance, radish::mlpe)
anova(fit_mlpe, fit_mlpe_reduced)

# test for an interaction
fit_mlpe_interaction <- radish(melip.Fst ~ forestcover * altitude, surface, 
                               radish::loglinear_conductance, radish::mlpe)
anova(fit_mlpe, fit_mlpe_interaction)

# test against null model of IBD
fit_mlpe_ibd <- radish(melip.Fst ~ 1, surface, 
                       radish::loglinear_conductance, radish::mlpe)
anova(fit_mlpe, fit_mlpe_ibd)

# categorical covariates:
# rasters of categorical covariates must have an associated RAT, see ?raster::ratify
# the names of levels are taken from the VALUE column of the RAT, if it exists (otherwise,
# the integer codes are used)
forestcover_class <- cut(raster::values(melip.forestcover), breaks = c(0, 1/3, 1/6, 1)) 
melip.forestcover_cat <- 
  raster::ratify(raster::setValues(melip.forestcover, as.numeric(forestcover_class)))

RAT <- levels(melip.forestcover_cat)[[1]]
RAT$VALUE <- levels(forestcover_class) #explicitly defines names in the RAT
levels(melip.forestcover_cat) <- RAT

covariates_cat <- raster::stack(list(forestcover = melip.forestcover_cat,
                                     altitude = melip.altitude)) 

surface_cat <- conductance_surface(covariates_cat, melip.coords, directions = 8)

fit_mlpe_cat <- radish(melip.Fst ~ forestcover + altitude, surface_cat, 
                       radish::loglinear_conductance, radish::mlpe)

# contrast coding is the default for R, and for this conductance model
# the (non-identifiable) intercept is omitted (e.g. only relative
# differences in conductance among levels are identifiable from the data)
summary(fit_mlpe_cat) 

# example of lower level interface:
# compute negative loglikelihood, gradient, Hessian for a given choice of
# of the conductance parameters theta, using a different measurement model
# (radish::generalized_wishart)
radish_algorithm(radish::loglinear_conductance(~forestcover + altitude, surface$x), 
                 radish::generalized_wishart, surface, 
                 pmax(melip.Fst, 0), nu = 1000, theta = c(-0.3, 0.3), 
                 gradient = TRUE, hessian = TRUE)$hessian
# numerical verification (not run)
#numDeriv::hessian(function(x)
#     radish_algorithm(radish::loglinear_conductance(~forestcover + altitude, surface$x), 
#                      radish::generalized_wishart, surface, 
#                      pmax(melip.Fst, 0), nu = 1000, theta = x)$objective,
#                  c(-0.3, 0.3))



# Converting ENM layer to resistance

enm_conductance <- raster("D:/Github/MYVE_LGA/data/SpatialLayers/ENM_MedianROC.tif")
plot(enm_conductance)
enm_resistance <- 1 - (enm_conductance/1000)
str(enm_resistance)
writeRaster(enm_resistance, "D:/Github/MYVE_LGA/data/SpatialLayers/ENM_MedianROC_resistance.tif")

# Cropping Flora's layers ###########################################
cropold <- raster("D:/UM/TXbats/LGA/FloraFiles/cropland_binary_NALCMS_azeq_5km_clip_mask.tif")
cropdistold <- raster("D:/UM/TXbats/LGA/FloraFiles/Dist2cropland_NALCMS_azeq_5km_clip_mask.tif")
lights <- raster("D:/UM/TXbats/LGA/FloraFiles/lights_NASA_NnadirCompSnowFree_azeq_5km_clip_mask.tif")
plot(crop)
plot(dist2crop)
plot(lights)
plot(covariates1[["elevation"]])

crop <- (cropold + covariates1[["elevation"]]) - covariates1[["elevation"]]
dist2crop <- (cropdistold + covariates1[["elevation"]]) - covariates1[["elevation"]]
lights <- (lights + covariates1[["elevation"]]) - covariates1[["elevation"]]

writeRaster(crop,"D:/Github/MYVE_LGA/data/cropland_binary_NALCMS_azeq_5km_clip_masked.tif")
writeRaster(dist2crop,"D:/Github/MYVE_LGA/data/Dist2cropland_NALCMS_azeq_5km_clip_masked.tif")
writeRaster(lights,"D:/Github/MYVE_LGA/data/lights_NASA_NnadirCompSnowFree_azeq_5km_clip_mask.tif")


########################### Candidate modes and testing runtimes

start_time <- Sys.time()
fit_mlpe1 <- radish(gen_mat ~ enm + TRI + ghfp1, surface, 
                   radish::loglinear_conductance, radish::mlpe) # 
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
fit_mlpe1 <- radish(gen_mat ~ enm + dist2karst + TRI + ghfp1 + crop + dist2water, surface, 
                    radish::loglinear_conductance, radish::leastsquares) # 
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
fit_mlpe1 <- radish(gen_mat ~ enm + dist2karst + TRI + ghfp1 + crop + dist2water, surface, 
                    radish::loglinear_conductance, radish::generalized_wishart) # 
end_time <- Sys.time()
print(end_time - start_time)

# IBD model

IBD <- radish(gen_mat ~ 1, surface, radish::loglinear_conductance, radish::mlpe)
IBD <- radish(gen_mat ~ 1, surface, radish::loglinear_conductance, radish::generalized_wishart, nu=1000)

ibd_aov <- anova(IBD, uni_models[[1]])
anova(IBD, uni_models[[16]])

# Loop through mlpe version of all univariate models, compare with aic
uni_models <- list()
for(idx in 1:length(surface$covariates)){
  uni_models[[idx]] <- radish(reformulate(surface$covariates[idx], response = "gen_mat"), surface, 
                      radish::loglinear_conductance, radish::mlpe)
}

sumtest <- summary(uni_models[[1]])

TRI <- radish(reformulate(surface$covariates[15], response = "gen_mat"), surface, 
                            radish::loglinear_conductance, radish::mlpe)
crop <- radish(reformulate(surface_crop$covariates[[1]], response = "gen_mat"), surface_crop, 
                            radish::loglinear_conductance, radish::generalized_wishart, nu=1000)

ones <- radish(reformulate(surface$covariates[[17]], response = "gen_mat"), surface, 
               radish::loglinear_conductance, radish::mlpe)

anova(IBD, crop)
anova(IBD, uni_models[[1]])

for(idx in 1:length(uni_models)){
print(uni_models[[idx]])
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
testsum <- summary(uni_models[[1]])
summary(uni_models[[14]])


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
write.csv(models_out, "D:/UM/TXbats/LGA/Radish_univariate_models_out4.csv")


covariates_round2 <- raster::stack(list(dist2karst = raster::scale(raster(paste0("./SpatialLayers/",env_layers[3]))),
                                        crop = raster::scale(raster(paste0("./SpatialLayers/",env_layers[1]))),
                                         dist2water = raster::scale(raster(paste0("./SpatialLayers/",env_layers[4]))),
                                         elevation = raster::scale(raster(paste0("./SpatialLayers/",env_layers[5]))),
                                         enm = raster::scale(raster(paste0("./SpatialLayers/",env_layers[6]))),
                                         forest = raster::scale(raster(paste0("./SpatialLayers/",env_layers[7]))),
                                         ghfp1 = raster::scale(raster(paste0("./SpatialLayers/",env_layers[9])))
                                        ))
# Have identical NoData mask on all layers
na_mask <- sum(covariates_round2)
covariates_round2 <- raster::mask(covariates_round2, na_mask)
writeRaster(covariates_round2, filename = "covariates_stack_round2.grd", overwrite=TRUE)

# Read in raster stack from disk, shouldn't need individual layers above
covariates_round2 <- stack("D:/Github/MYVE_LGA/data/covariates_stack_round2.grd")
# Now taking best model (TRI) and running each two-variable model
twovar_models <- list()
surface_round2 <- conductance_surface(covariates_round2, myve_points, directions = 8)
covariates_round2 <- stack(list(covariates[[3]],covariates[[1]],covariates[[4]],covariates[[5]],
                           covariates[[6]],covariates[[7]],covariates[[9]]))
covariates_round2 <- stack(list(dist2karst=covariates[[3]],crop=covariates[[1]],dist2water=covariates[[4]],elevation=covariates[[5]],
                                enm=covariates[[6]],forest=covariates[[7]],ghfp1=covariates[[9]]))
covariates_round2 <- covariates[[c(3,1,4,5,6,7,9)]]
for(idx in 1:(length(surface_round2$covariates)-1)){
  twovar_models[[idx]] <- radish(reformulate(termlabels = c(surface_round2$covariates[1],surface_round2$covariates[idx+1]), response = "gen_mat"), surface_round2, 
                              radish::loglinear_conductance, radish::mlpe)
  }


for(idx in 1:length(twovar_models)){
  model_aics[idx] <- summary(twovar_models[[idx]])$aic
  model_lglk[idx] <- summary(twovar_models[[idx]])$loglik
  model_gradnorm[idx] <- summary(twovar_models[[idx]])$gradnorm
  model_coef[idx] <- summary(twovar_models[[idx]])$ztable[1]
  model_ste[idx] <- summary(twovar_models[[idx]])$ztable[2]
  model_zval[idx] <- summary(twovar_models[[idx]])$ztable[3]
  model_pval[idx] <- summary(twovar_models[[idx]])$ztable[4]
  model_alpha[idx] <- summary(twovar_models[[idx]])$phi[1]
  model_beta[idx] <- summary(twovar_models[[idx]])$phi[2]
  model_tau[idx] <- summary(twovar_models[[idx]])$phi[3]
  model_rho[idx] <- summary(twovar_models[[idx]])$phi[4]
  model_iters[idx] <- summary(twovar_models[[idx]])$iters
}

models_out <- data.frame(aic = model_aics, loglike = model_lglk, gradnorm = model_gradnorm,
                         coefficient = model_coef, sterror = model_ste, zval = model_zval, pval = model_pval,
                         alpha = model_alpha, beta = model_beta, tau = model_tau, rho = model_rho, iterations = model_iters)
write.csv(models_out, "D:/UM/TXbats/LGA/Radish_twovariate_models_out.csv")


# Three variable models
covariates_round3 <- raster::stack(list(dist2karst = raster::scale(raster(paste0("./SpatialLayers/",env_layers[3]))),
                                        enm = raster::scale(raster(paste0("./SpatialLayers/",env_layers[6]))),
                                        crop = raster::scale(raster(paste0("./SpatialLayers/",env_layers[1]))),
                                        dist2water = raster::scale(raster(paste0("./SpatialLayers/",env_layers[4]))),
                                        elevation = raster::scale(raster(paste0("./SpatialLayers/",env_layers[5]))),
                                        forest = raster::scale(raster(paste0("./SpatialLayers/",env_layers[7]))),
                                        ghfp1 = raster::scale(raster(paste0("./SpatialLayers/",env_layers[9])))
))
# Have identical NoData mask on all layers
na_mask <- sum(covariates_round3)
covariates_round3 <- raster::mask(covariates_round3, na_mask)
writeRaster(covariates_round3, filename = "covariates_stack_round3.grd", overwrite=TRUE)

# Read in raster stack from disk, shouldn't need individual layers above
covariates_round3 <- stack("D:/Github/MYVE_LGA/data/covariates_stack_round3.grd")
# Now taking best model (TRI) and running each two-variable model
threevar_models <- list()
surface_round3 <- conductance_surface(covariates_round3, myve_points, directions = 8)

for(idx in 1:(length(surface_round3$covariates)-2)){
  threevar_models[[idx]] <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[idx+2]), response = "gen_mat"), surface_round3, 
                                 radish::loglinear_conductance, radish::mlpe)
}


for(idx in 1:length(threevar_models)){
  model_aics[idx] <- summary(threevar_models[[idx]])$aic
  model_lglk[idx] <- summary(threevar_models[[idx]])$loglik
  model_gradnorm[idx] <- summary(threevar_models[[idx]])$gradnorm
  model_coef[idx] <- summary(threevar_models[[idx]])$ztable[1]
  model_ste[idx] <- summary(threevar_models[[idx]])$ztable[2]
  model_zval[idx] <- summary(threevar_models[[idx]])$ztable[3]
  model_pval[idx] <- summary(threevar_models[[idx]])$ztable[4]
  model_alpha[idx] <- summary(threevar_models[[idx]])$phi[1]
  model_beta[idx] <- summary(threevar_models[[idx]])$phi[2]
  model_tau[idx] <- summary(threevar_models[[idx]])$phi[3]
  model_rho[idx] <- summary(threevar_models[[idx]])$phi[4]
  model_iters[idx] <- summary(threevar_models[[idx]])$iters
}

models_out <- data.frame(aic = model_aics, loglike = model_lglk, gradnorm = model_gradnorm,
                         coefficient = model_coef, sterror = model_ste, zval = model_zval, pval = model_pval,
                         alpha = model_alpha, beta = model_beta, tau = model_tau, rho = model_rho, iterations = model_iters)
write.csv(models_out, "D:/UM/TXbats/LGA/Radish_threevariate_models_out.csv")
summary(threevar_models[[1]])

# Four variable models
threeplusWater <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],surface_round3$covariates[4]), response = "gen_mat"), surface_round3, 
                                               radish::loglinear_conductance, radish::mlpe)
threepluselevation <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],surface_round3$covariates[5]), response = "gen_mat"), surface_round3, 
                         radish::loglinear_conductance, radish::mlpe)
threeplusforest <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],surface_round3$covariates[6]), response = "gen_mat"), surface_round3, 
                         radish::loglinear_conductance, radish::mlpe)
threeplusghfp <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],surface_round3$covariates[7]), response = "gen_mat"), surface_round3, 
                         radish::loglinear_conductance, radish::mlpe)

summary(threeplusWater) #AIC -25073
summary(threepluselevation) #AIC -25043
summary(threeplusforest) #AIC -25202
summary(threeplusghfp) #AIC -25124

# Five variable models
fourplusWater <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],
                                                    surface_round3$covariates[6], surface_round3$covariates[4]),response = "gen_mat"), surface_round3, 
                         radish::loglinear_conductance, radish::mlpe)
fourpluselevation <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],
                                                       surface_round3$covariates[6], surface_round3$covariates[5]),response = "gen_mat"), surface_round3, 
                             radish::loglinear_conductance, radish::mlpe)
fourplusghfp <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],
                                                  surface_round3$covariates[6], surface_round3$covariates[7]),response = "gen_mat"), surface_round3, 
                          radish::loglinear_conductance, radish::mlpe)

summary(fourplusWater) #AIC -25218
summary(fourpluselevation) #AIC -25200
summary(fourplusghfp) #AIC -25208


fivepluselevation <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],
                                                       surface_round3$covariates[6], surface_round3$covariates[4],surface_round3$covariates[5]),response = "gen_mat"), surface_round3, 
                            radish::loglinear_conductance, radish::mlpe)
fiveplusghfp <- radish(reformulate(termlabels = c(surface_round3$covariates[1],surface_round3$covariates[2],surface_round3$covariates[3],
                                                  surface_round3$covariates[6], surface_round3$covariates[4],surface_round3$covariates[7]),response = "gen_mat"), surface_round3, 
                       radish::loglinear_conductance, radish::mlpe)







########################## ResistanceGA testing ############################

library(devtools)
devtools::install_github("wpeterman/ResistanceGA2", build_vignettes = TRUE)

library(ResistanceGA2)











