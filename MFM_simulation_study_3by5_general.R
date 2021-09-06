########################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack",
              "mvtnorm","dplyr","armspp","ggcorrplot","reshape2","ggpubr",
              "cowplot","doParallel","foreach","MixMatrix","abind",
              "mclust", "mcclust","fossil","kernlab", "mixtools",
              "npmr","tibble","doMC","clusterGeneration","ClusterR")
# "rjson","ggplot2","grid","inlabru","jpeg",
# "clusterGeneration"
#  "RCurl","sp","INLA","Orcs","dplyr","gridExtra","ggsci"
# https://cran.r-project.org/web/packages/MixMatrix/vignettes/matrixnormal.html
ipak(packages)
##########################
source("./MFM_MN.R")
#
M_list_3by5 <- readRDS("./M_3by5_list.rds")
#
MCMC.total <- 1000
MCMC.burnin <- ceiling(MCMC.total*0.30)
# MCMC.thin <- 10
# 
# registerDoParallel(cores=50)
# helper function AR(1) correlation matrix
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
###
# http://sherrytowers.com/2013/10/24/k-means-clustering/
# AIC for kmeans
kmeansAIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}
# BIC for kmeans
kmeansBIC <- function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + log(n)*m*k)
}
################
simulation_settings <- data.frame(expand.grid(signal_strength = 1, 
                                              noise_factor = c(0.5, 1, 1.5), 
                                              cluster_num = c(3),
                                              total_sample_size = c(500))) %>% rowid_to_column("setting_id")

################
iii = as.numeric(readline(prompt = "Enter the setting id: "))
################
signal_strength <- simulation_settings$signal_strength[iii]
noise_factor <- simulation_settings$noise_factor[iii]
cluster_num <- simulation_settings$cluster_num[iii]
total_sample_size <- simulation_settings$total_sample_size[iii]
setting_id <- simulation_settings$setting_id[iii]
###############
for(i in 1:cluster_num){
  assign(paste0("M",i), M_list_3by5[[i]])
}
# M2 : dumbell
# total_sample_size <- c(100,200,400)
num_replicates <- 100
# 
if(cluster_num == 3){
  cluster_proportion <- c(0.3,0.4,0.3)
} else if(cluster_num == 5){
  cluster_proportion <- rep(1/cluster_num, cluster_num)
}
# cluster_proportion <- c(0.25,0.25,0.25,0.25)
###
# registerDoParallel(cores = 50)
registerDoMC(20)
# sim_10by6 <- vector(mode="list", length = length(total_sample_size))
###
sim_3by5 <- foreach(j = 1:num_replicates, .errorhandling = "pass")%dopar%{
  # 
  set.seed(j+1234*setting_id)
  V <- (noise_factor^2)*clusterGeneration::rcorrmatrix(d = ncol(M1))
  U <- clusterGeneration::rcorrmatrix(d = nrow(M1))
  # U <- cov2cor(MCMCpack::rwish(nrow(M1)+1,diag(rep(1,nrow(M1)))))
  #
  if(cluster_num == 3){
    dat_sim_cluster1 <-  mniw::rMNorm(total_sample_size*cluster_proportion[1], 
                                      Lambda = M1*signal_strength, 
                                      SigmaR = U, SigmaC = V)
    dat_sim_cluster2 <-  mniw::rMNorm(total_sample_size*cluster_proportion[2], 
                                      Lambda = M2*signal_strength, 
                                      SigmaR = U, SigmaC = V)
    dat_sim_cluster3 <-  mniw::rMNorm(total_sample_size*cluster_proportion[3], 
                                      Lambda = M3*signal_strength, 
                                      SigmaR = U, SigmaC = V)
    # 
    dat_sim <- abind(dat_sim_cluster1, dat_sim_cluster2, dat_sim_cluster3,
                     along = 3)
    # 
    true_cluster_membership <- rep(c(1,2,3), c(total_sample_size*cluster_proportion[1],
                                               total_sample_size*cluster_proportion[2],
                                               total_sample_size*cluster_proportion[3]))
    
  } else if(cluster_num == 5){
    dat_sim_cluster1 <-  mniw::rMNorm(total_sample_size*cluster_proportion[1], 
                                      Lambda = M1*signal_strength, 
                                      SigmaR = U, SigmaC = V)
    dat_sim_cluster2 <-  mniw::rMNorm(total_sample_size*cluster_proportion[2], 
                                      Lambda = M2*signal_strength, 
                                      SigmaR = U, SigmaC = V)
    dat_sim_cluster3 <-  mniw::rMNorm(total_sample_size*cluster_proportion[3], 
                                      Lambda = M3*signal_strength, 
                                      SigmaR = U, SigmaC = V)
    dat_sim_cluster4 <-  mniw::rMNorm(total_sample_size*cluster_proportion[4], 
                                      Lambda = M4*signal_strength, 
                                      SigmaR = U, SigmaC = V)
    dat_sim_cluster5 <-  mniw::rMNorm(total_sample_size*cluster_proportion[5], 
                                      Lambda = M5*signal_strength, 
                                      SigmaR = U, SigmaC = V)
    
    dat_sim <- abind(dat_sim_cluster1, dat_sim_cluster2, dat_sim_cluster3,  dat_sim_cluster4,
                     dat_sim_cluster5,
                     along = 3)
    
    true_cluster_membership <- rep(c(1,2,3,4,5), c(total_sample_size*cluster_proportion[1],
                                                   total_sample_size*cluster_proportion[2],
                                                   total_sample_size*cluster_proportion[3],
                                                   total_sample_size*cluster_proportion[4],
                                                   total_sample_size*cluster_proportion[5]))
  }
  # combine them
  dat_sim_vector <- t(apply(dat_sim, c(3), c))
  
  # pairwise distance between matrices
  # nuclear norm, spectral norm
  dat_sim_matrix_dist_spectral <- matrix(0, nrow = total_sample_size, ncol = total_sample_size)
  dat_sim_matrix_dist_nuclear <- matrix(0, nrow = total_sample_size, ncol = total_sample_size)
  for(ii in 1:total_sample_size){
    for(jj in 1:total_sample_size){
      if(ii != jj){
        dat_sim_matrix_dist_spectral[ii,jj] <- norm(dat_sim[,,ii] - dat_sim[,,jj], type = "2")
        dat_sim_matrix_dist_nuclear[ii,jj] <- npmr::nuclear(dat_sim[,,ii] - dat_sim[,,jj])
      }
    }
  }
  # dat_sim <- array(NA, dim = c(nrow(M1),ncol(M1),total_sample_size))
  # run MFM analysis
  temp_time1 <- Sys.time()
  temp_MFM_rlt <- MFM_MxN_equal_cov(data=dat_sim, 
                                    niterations=MCMC.total, 
                                    alpha=(1+dim(dat_sim)[1])/2, beta=diag(rep(0.5,dim(dat_sim)[1])), 
                                    psi=(1+dim(dat_sim)[2])/2, rho=diag(rep(0.5,dim(dat_sim)[2])), 
                                    GAMMA=1, 
                                    M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2, 
                                    Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 ), 
                                    Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 ), 
                                    initNClusters=cluster_num*4, 
                                    MLE.initial=TRUE)
  # 
  # RI_MFM_trace <- rep(NA, MCMC.total)
  # ARI_MFM_trace <- rep(NA, MCMC.total)
  # for(ii in 1:MCMC.total){
  #   RI_MFM_trace[ii] <- fossil::rand.index(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout)
  #   ARI_MFM_trace[ii] <- fossil::adj.rand.index(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout)
  # }
  RI_MFM_trace <- sapply(1:MCMC.total, function(x) fossil::rand.index(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout))
  ARI_MFM_trace <- sapply(1:MCMC.total, function(x) mclust::adjustedRandIndex(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout))
  ###
  temp_time2 <- Sys.time()
  temp_DP_rlt <- DP_MxN_equal_cov(data=dat_sim, 
                                  niterations=MCMC.total, 
                                  alpha=(1+dim(dat_sim)[1])/2, beta=diag(rep(0.5,dim(dat_sim)[1])), 
                                  psi=(1+dim(dat_sim)[2])/2, rho=diag(rep(0.5,dim(dat_sim)[2])), 
                                  GAMMA=1, 
                                  M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2, 
                                  Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 ), 
                                  Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 ), 
                                  initNClusters=cluster_num*4, 
                                  MLE.initial=TRUE)
  #
  RI_DP_trace <- sapply(1:MCMC.total, function(x) fossil::rand.index(true_cluster_membership,temp_DP_rlt$Iterates[[ii]]$zout))
  ARI_DP_trace <- sapply(1:MCMC.total, function(x) mclust::adjustedRandIndex(true_cluster_membership,temp_DP_rlt$Iterates[[ii]]$zout))
  temp_time3 <- Sys.time()
  # do.call("c", lapply( ,function(x) fossil::rand.index(true_cluster_membership,x$zout) ) )
  temp_MFM_Dahl_rlt <- getDahl(temp_MFM_rlt,burn=MCMC.burnin)
  temp_DP_Dahl_rlt <- getDahl(temp_DP_rlt,burn=MCMC.burnin)
  ## k-means
  
  k_vec <- 1 : (total_sample_size - 1)
  temp_kmeans_list <- vector(mode="list", length=length(k_vec))
  temp_kmeans_aic <- rep(NA, length=length(k_vec))
  temp_kmeans_bic <- rep(NA, length=length(k_vec))
  ##
  for(k in 1:length(k_vec)){
    # kmeans each row is a data point
    temp_kmeans_list[[k]] <- kmeans(dat_sim_vector, 
                                    centers = k_vec[k])
    temp_kmeans_aic[k] <- kmeansAIC(temp_kmeans_list[[k]])
    temp_kmeans_bic[k] <- kmeansBIC(temp_kmeans_list[[k]])
  }
  # 
  temp_kmeans_rlt_aic <- temp_kmeans_list[which(temp_kmeans_aic==min(temp_kmeans_aic))]
  temp_kmeans_rlt_bic <- temp_kmeans_list[which(temp_kmeans_bic==min(temp_kmeans_bic))]
  temp_kmeans_rlt_oracle <- temp_kmeans_list[which(k_vec==cluster_num)]
  #
  if(max(temp_MFM_Dahl_rlt$zout) == total_sample_size){
    temp_kmeans_rlt_MFM <- vector(mode = "list", length = 1)
    temp_kmeans_rlt_MFM[[1]]$cluster <- 1:total_sample_size
  } else{
    temp_kmeans_rlt_MFM <- temp_kmeans_list[which(k_vec==max(temp_MFM_Dahl_rlt$zout))] 
  }

  temp_specc_rlt_oracle <- kernlab::specc(dat_sim_vector, centers=cluster_num)
  temp_specc_rlt_MFM <- try(kernlab::specc(dat_sim_vector, centers=max(temp_MFM_Dahl_rlt$zout)))
  # temp_time <- Sys.time()
  
  # finite mixtures of gaussian
  temp_mvnormalmixEM_rlt_oracle <- try(mvnormalmixEM(dat_sim_vector,
                                                     k = cluster_num))
  temp_mvnormalmixEM_rlt_MFM <- try(mvnormalmixEM(dat_sim_vector,
                                                  k = max(temp_MFM_Dahl_rlt$zout)))
  
  # matrix norm k-centroid clustering
  temp_kcentroid_rlt_nuclear_oracle <- Cluster_Medoids(dat_sim_matrix_dist_nuclear,
                                                       clusters = cluster_num)
  temp_kcentroid_rlt_spectral_oracle <- Cluster_Medoids(dat_sim_matrix_dist_spectral,
                                                        clusters = cluster_num)
  #
  temp_kcentroid_rlt_nuclear_MFM <- Cluster_Medoids(dat_sim_matrix_dist_nuclear,
                                                    clusters = cluster_num)
  temp_kcentroid_rlt_spectral_MFM <- Cluster_Medoids(dat_sim_matrix_dist_spectral,
                                                     clusters = cluster_num)
  #########################
  #########################
  # clustering performance
  ARI_MFM <-  mclust::adjustedRandIndex(true_cluster_membership,
                                        temp_MFM_Dahl_rlt$zout)
  ARI_DP <-  mclust::adjustedRandIndex(true_cluster_membership,
                                       temp_DP_Dahl_rlt$zout)
  #
  ARI_kmeans_aic <- mclust::adjustedRandIndex(true_cluster_membership,
                                              temp_kmeans_rlt_aic[[1]]$cluster)
  ARI_kmeans_bic <- mclust::adjustedRandIndex(true_cluster_membership,
                                              temp_kmeans_rlt_bic[[1]]$cluster)
  ARI_kmeans_oracle <- mclust::adjustedRandIndex(true_cluster_membership,
                                                 temp_kmeans_rlt_oracle[[1]]$cluster)
  
  ARI_kmeans_MFM <- mclust::adjustedRandIndex(true_cluster_membership,
                                              temp_kmeans_rlt_MFM[[1]]$cluster)
  
  ARI_specc_MFM <- ifelse(class(temp_specc_rlt_MFM) == "try-error", 
                          NA,
                          mclust::adjustedRandIndex(true_cluster_membership,
                                                    temp_specc_rlt_MFM@.Data))
  
  ARI_specc_oracle <- mclust::adjustedRandIndex(true_cluster_membership,
                                                temp_specc_rlt_oracle@.Data)
  #
  ARI_mvnormalmixEM_oracle <- ifelse(class(temp_mvnormalmixEM_rlt_oracle) == "try-error",
                                     NA,mclust::adjustedRandIndex(true_cluster_membership,
                                                                  apply(temp_mvnormalmixEM_rlt_oracle$posterior,
                                                                        1, function(x) which(x == max(x)))) )
  #  
  ARI_mvnormalmixEM_MFM <- ifelse(class(temp_mvnormalmixEM_rlt_MFM) == "try-error",
                                  NA,mclust::adjustedRandIndex(true_cluster_membership,
                                                               apply(temp_mvnormalmixEM_rlt_MFM$posterior,
                                                                     1, function(x) which(x == max(x)))) )
  # 
  ARI_kcentroid_nuclear_oracle <- mclust::adjustedRandIndex(true_cluster_membership,
                                                            temp_kcentroid_rlt_nuclear_oracle$clusters)
  
  ARI_kcentroid_nuclear_MFM <- mclust::adjustedRandIndex(true_cluster_membership,
                                                         temp_kcentroid_rlt_nuclear_MFM$clusters)
  
  ARI_kcentroid_spectral_oracle <- mclust::adjustedRandIndex(true_cluster_membership,
                                                             temp_kcentroid_rlt_spectral_oracle$clusters)
  
  ARI_kcentroid_spectral_MFM <- mclust::adjustedRandIndex(true_cluster_membership,
                                                          temp_kcentroid_rlt_spectral_MFM$clusters)
  
  # 
  #
  # RI_MFM <- fossil::rand.index(true_cluster_membership,
  #                              temp_MFM_Dahl_rlt$zout)
  # RI_kmeans_aic <- fossil::rand.index(true_cluster_membership,
  #                                     temp_kmeans_rlt_aic[[1]]$cluster)
  # RI_kmeans_bic <- fossil::rand.index(true_cluster_membership,
  #                                     temp_kmeans_rlt_bic[[1]]$cluster)
  # RI_kmeans_oracle <- fossil::rand.index(true_cluster_membership,
  #                                   temp_kmeans_rlt_oracle[[1]]$cluster)
  # RI_kmeans_MFM <- fossil::rand.index(true_cluster_membership,
  #                                     temp_kmeans_rlt_MFM[[1]]$cluster)
  # RI_specc_MFM <- fossil::rand.index(true_cluster_membership,
  #                                    temp_specc_rlt_MFM@.Data)
  # RI_specc_oracle <- fossil::rand.index(true_cluster_membership,
  #                                       temp_specc_rlt_oracle@.Data)
  
  # 
  Khat_MFM <- max(temp_MFM_Dahl_rlt$zout)
  Khat_DP <- max(temp_DP_Dahl_rlt$zout)
  Khat_kmeans_aic <- max(temp_kmeans_rlt_aic[[1]]$cluster)
  Khat_kmeans_bic <- max(temp_kmeans_rlt_bic[[1]]$cluster)
  Khat_kmeans_oracle <- cluster_num
  Khat_kmeans_MFM <- Khat_MFM
  Khat_specc_oracle <- cluster_num
  Khat_specc_MFM <- Khat_MFM 
  # RMSE V %x% U
  # VU_RMSE_MFM <- NA
  VU_RMSE_MFM <- sqrt(mean(((temp_MFM_Dahl_rlt$Vout %x% temp_MFM_Dahl_rlt$Uout) - (V %x% U))^2))
  VU_RMSE_DP <- sqrt(mean(((temp_DP_Dahl_rlt$Vout %x% temp_DP_Dahl_rlt$Uout) - (V %x% U))^2))
  # 
  M_RMSE_MFM <- sqrt(mean((temp_MFM_Dahl_rlt$Mout[,,temp_MFM_Dahl_rlt$zout] - dat_sim)^2))
  M_RMSE_DP <- sqrt(mean((temp_DP_Dahl_rlt$Mout[,,temp_DP_Dahl_rlt$zout] - dat_sim)^2))
  # M_RMSE_kmeans_aic <- sqrt(mean((temp_kmeans_rlt_aic[[1]]$centers[temp_kmeans_rlt_aic[[1]]$cluster,] - dat_sim_vector)^2))
  # M_RMSE_kmeans_bic <- sqrt(mean((temp_kmeans_rlt_bic[[1]]$centers[temp_kmeans_rlt_bic[[1]]$cluster,] - dat_sim_vector)^2))
  # M_RMSE_kmeans_3 <- sqrt(mean((temp_kmeans_rlt_oracle[[1]]$centers[temp_kmeans_rlt_oracle[[1]]$cluster,] - dat_sim_vector)^2))
  # M_RMSE_kmeans_MFM <- sqrt(mean((temp_kmeans_rlt_MFM[[1]]$centers[temp_kmeans_rlt_MFM[[1]]$cluster,] - dat_sim_vector)^2))
  # M_RMSE_specc_3 <- sqrt(mean((temp_specc_rlt_oracle@centers[temp_specc_rlt_oracle@.Data,] - dat_sim_vector)^2))
  # M_RMSE_specc_MFM <- sqrt(mean((temp_specc_rlt_MFM@centers[temp_specc_rlt_MFM@.Data,] - dat_sim_vector)^2))
  # 
  temp_rlt <- list(MFM_Dahl_rlt = temp_MFM_Dahl_rlt,
                   DP_Dahl_rlt = temp_DP_Dahl_rlt,
                   kmeans_rlt_aic = temp_kmeans_rlt_aic[[1]],
                   kmeans_rlt_bic = temp_kmeans_rlt_bic[[1]],
                   kmeans_rlt_oracle = temp_kmeans_rlt_oracle[[1]],
                   MFM_DP_time = difftime(temp_time3, temp_time2, units = "mins"),
                   MFM_Dahl_time = difftime(temp_time2,temp_time1,units = "mins"),
                   Khat_MFM=Khat_MFM,
                   Khat_DP=Khat_DP,
                   Khat_kmeans_aic=Khat_kmeans_aic,
                   Khat_kmeans_bic=Khat_kmeans_bic,
                   Khat_kmeans_oracle=Khat_kmeans_oracle,
                   Khat_kmeans_MFM=Khat_kmeans_MFM,
                   Khat_specc_oracle=Khat_specc_oracle,
                   Khat_specc_MFM=Khat_specc_MFM,
                   ARI_MFM=ARI_MFM,
                   ARI_DP=ARI_DP,
                   ARI_kmeans_aic=ARI_kmeans_aic,
                   ARI_kmeans_bic=ARI_kmeans_bic,
                   ARI_kmeans_oracle=ARI_kmeans_oracle,
                   ARI_kmeans_MFM=ARI_kmeans_MFM,
                   ARI_specc_oracle=ARI_specc_oracle,
                   ARI_specc_MFM=ARI_specc_MFM,
                   ARI_mvnormalmixEM_oracle = ARI_mvnormalmixEM_oracle,
                   ARI_mvnormalmixEM_MFM = ARI_mvnormalmixEM_MFM,
                   ARI_kcentroid_nuclear_oracle = ARI_kcentroid_nuclear_oracle,
                   ARI_kcentroid_nuclear_MFM = ARI_kcentroid_nuclear_MFM,
                   ARI_kcentroid_spectral_oracle = ARI_kcentroid_spectral_oracle,
                   ARI_kcentroid_spectral_MFM = ARI_kcentroid_spectral_MFM,
                   RI_MFM_trace=RI_MFM_trace,
                   ARI_MFM_trace=ARI_MFM_trace,
                   VU_RMSE_MFM=VU_RMSE_MFM,
                   VU_RMSE_DP=VU_RMSE_DP,
                   M_RMSE_MFM=M_RMSE_MFM,
                   M_RMSE_DP=M_RMSE_DP)
  temp_rlt
}
#####################################
# saveRDS(sim_10by6,
#         paste0("sim_10by6_",
#                signal_strength,
#                "_", noise_factor, "_", cluster_num, "clusters",".rds"))
saveRDS(sim_3by5,
        paste0("./simulation_results/sim_3by5_",setting_id,".rds"))
#####################################
#rm(sim_10by6)
#####################################