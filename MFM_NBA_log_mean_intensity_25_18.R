######################
######################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack",
              "mvtnorm","dplyr","armspp","ggcorrplot","reshape2","ggpubr",
              "cowplot","doParallel","foreach","MixMatrix")
# "rjson","ggplot2","grid","inlabru","jpeg",
#  "RCurl","sp","INLA","Orcs","dplyr","gridExtra","ggsci"
# https://cran.r-project.org/web/packages/MixMatrix/vignettes/matrixnormal.html
ipak(packages)
######################
registerDoParallel(cores=50)
#setwd("/Users/fan/Documents/MFM_matrix/")
player_names_caps <- readRDS("./player_names_caps.rds")
load("./shotdata.Rdata")
load("./lgcp_fit.Rdata")
source("./MFM_MN.R")
#######################
dim(intensity.coef)
summary(coords) # y: 0, 50 (sideline to sideline); x: 0 - 46 (mid court line to baseline)
NameList
#######################
#######################
intensity_sum_array <- array(0, dim=c(25,18,nrow(intensity.coef)))
intensity_count_array <- array(0, dim=c(25,18,nrow(intensity.coef)))
for(k in 1:nrow(intensity.coef)){
  # temp_player_name <- shotDataf$PLAYER_NAME[l]
  for(l in 1:ncol(intensity.coef)){
    i <- ceiling(coords[l,2]/2)
    j <- ceiling(coords[l,1]/2)
    if(j > 18){
      next
    }
    intensity_sum_array[i,j,k] <- intensity_sum_array[i,j,k] + intensity.coef[k,l]
    intensity_count_array[i,j,k] <- intensity_count_array[i,j,k] + 1
  }
  # k <- which( (player_names_caps$name == temp_player_name) == TRUE ) # work with the k-th player in player_name_caps
  # i <- floor(shotDataf$LOC_X[l]/3)
  # j <- floor(shotDataf$LOC_Y[l]/3)
  # shotData_array[i,j,k] <- shotData_array[i,j,k] + 1
}
##########################
log_intensity_array <- log((intensity_sum_array/intensity_count_array))# log of the mean intensity
log_mean_intensity_array <- array(NA, dim=c(25,18,nrow(intensity.coef)))
for(k in 1:nrow(intensity.coef)){
  log_mean_intensity_array[,,k] <- log((intensity_sum_array[,,k]/intensity_count_array[,,k])/player_names_caps$caps[k])
}
# log_intensity_array <- log(intensity_sum_array/intensity_count_array) # log of the mean intensity
intensity_array <- intensity_sum_array/intensity_count_array # 
##########################
##########################
time1 <- Sys.time()
NBA_rlt <- foreach(i = 1:50) %dopar% {
  # MLE_shotData <- MLmatrixnorm(log_mean_intensity_array)
  temp_NBA_rlt <- CDMFM_new1_log(data=log_mean_intensity_array, 
                                 niterations=6000, 
                                 alpha=(1+25*2)/2, beta=diag(rep(0.5,dim(log_mean_intensity_array)[1])), 
                                 psi=(1+18*2)/2, rho=diag(rep(0.5,dim(log_mean_intensity_array)[2])), 
                                 GAMMA=3, 
                                 M0=(apply(log_mean_intensity_array, c(1,2), max) + apply(log_mean_intensity_array, c(1,2), min))/2, 
                                 Sigma0=diag( ((apply(log_mean_intensity_array, c(1), max)-apply(log_mean_intensity_array, c(1), min))/4)^2 ), 
                                 Omega0=diag( ((apply(log_mean_intensity_array, c(2), max)-apply(log_mean_intensity_array, c(2), min))/4)^2 ), 
                                 initNClusters=sample(2:12,size=1), 
                                 VN=VN, 
                                 MLE.initial=TRUE)
  getDahl(temp_NBA_rlt,burn=4000)
}
time2 <- Sys.time()
##########################
saveRDS(NBA_rlt, "NBA_log_mean_intensity_25_18_MLE_initial_6000iters_random_initclusters.rds") 
##########################