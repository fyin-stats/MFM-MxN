#######################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack",
              "mvtnorm","dplyr","armspp","ggcorrplot","reshape2","ggpubr",
              "cowplot","doParallel","foreach","MixMatrix","mclust","mcclust","fossil","abind","kernlab")
# "rjson","ggplot2","grid","inlabru","jpeg",
#  "RCurl","sp","INLA","Orcs","dplyr","gridExtra","ggsci"
# https://cran.r-project.org/web/packages/MixMatrix/vignettes/matrixnormal.html
ipak(packages)
#########################################################
registerDoParallel(cores=50)
source("MFM_MN.R")
### simulation study 17*12 based on the results from NBA data
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
#########################################################
# http://sherrytowers.com/2013/10/24/k-means-clustering/
########################################### AIC for kmeans
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
#setwd("/Users/fan/Documents/MFM_matrix/")
##
load("./Group1.Rdata")
load("./Group2.Rdata")
load("./Group3.Rdata")
load("./Group4.Rdata")
load("./Group5.Rdata")
load("./Group6.Rdata")
load("./Group7.Rdata")
load("./Group8.Rdata")
load("./Lambda.Rdata")
######################
M1_raw <- matrix(Lambda[[Group3[1]]]$mean,nrow=112,ncol=112)
M2_raw <- matrix(Lambda[[Group5[1]]]$mean,nrow=112,ncol=112)
M3_raw <- matrix(Lambda[[Group7[1]]]$mean,nrow=112,ncol=112)
########################
# M1 <- matrix(NA, nrow=25, ncol = 18)
# M2 <- matrix(NA, nrow=25, ncol = 18)
# M3 <- matrix(NA, nrow=25, ncol = 18)
########################
M1_sum <- matrix(0, nrow=25, ncol = 18)
M2_sum <- matrix(0, nrow=25, ncol = 18)
M3_sum <- matrix(0, nrow=25, ncol = 18)
M1_count <- matrix(0, nrow=25, ncol = 18)
M2_count <- matrix(0, nrow=25, ncol = 18)
M3_count <- matrix(0, nrow=25, ncol = 18)
#########################
for(k in 1:nrow(M1_raw)){
  for(l in 1:ncol(M1_raw)){
    i <- ceiling(k/4.5)
    j <- ceiling(l/5)
    if(j > 18){
      next
    }
    M1_sum[i,j] <- M1_sum[i,j] + M1_raw[k,l]
    M1_count[i,j] <- M1_count[i,j] + 1
    M2_sum[i,j] <- M2_sum[i,j] + M2_raw[k,l]
    M2_count[i,j] <- M2_count[i,j] + 1
    M3_sum[i,j] <- M3_sum[i,j] + M3_raw[k,l]
    M3_count[i,j] <- M3_count[i,j] + 1
  }
}
#########################
M1 <- log(M1_sum/M1_count)
M2 <- log(M2_sum/M2_count)
M3 <- log(M3_sum/M3_count)
#########################
rownames(M1) <- 1:25
colnames(M1) <- 1:18
M1_long <- melt(M1)
M1p <- ggplot(M1_long, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="white",high="blue") +
  labs(x="X", y="Y", title=paste("M")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

###
rownames(M2) <- 1:25
colnames(M2) <- 1:18
M2_long <- melt(M2)
M2p <- ggplot(M2_long, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="white",high="blue") +
  labs(x="X", y="Y", title=paste("M")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

###
rownames(M3) <- 1:25
colnames(M3) <- 1:18
M3_long <- melt(M3)
M3p <- ggplot(M3_long, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="white",high="blue") +
  labs(x="X", y="Y", title=paste("M")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))
###
# ggarrange(M1p, M2p, M3p)
###
jpeg("shape_25by18_simulation_M1.jpg")
image(x=1:25,y=1:18, z=M1, 
      useRaster=TRUE, 
      col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=18), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=25), labels=rownames(M1), lwd=0, pos=0)
dev.off()
###
jpeg("shape_25by18_simulation_M2.jpg")
image(x=1:25,y=1:18, z=M2, 
      useRaster=TRUE, 
      col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=18), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=25), labels=rownames(M1), lwd=0, pos=0)
dev.off()
###

# NBA_rlt_Dahl_log_intensity_25_18_50reps <- readRDS("NBA_log_intensity_25_18_MLE_initial_6000iters_random_initclusters.rds")
# NBA_rlt_Dahl_log_intensity <- NBA_rlt_Dahl_log_intensity_25_18_50reps[[5]]
# ###
# M1 <- NBA_rlt_Dahl_log_intensity$Mout[,,1]
# M2 <- NBA_rlt_Dahl_log_intensity$Mout[,,2]
# M3 <- NBA_rlt_Dahl_log_intensity$Mout[,,3]
# M4 <- NBA_rlt_Dahl_log_intensity$Mout[,,4]
###
# V <- NBA_rlt_Dahl_log_intensity$Vout
# # set.seed(0)
# U <- NBA_rlt_Dahl_log_intensity$Uout
### 0.2, 0.2, 0.3, 0.3
### sample size: 100, 200, 300
### number of replicates: 200
total_sample_size <- c(200)
rho <- c(0.9,0.6,0.3)
noise_sd <- c(1.5,1,0.5)
simulation_settings_df <- data.frame(expand.grid(rho, noise_sd))
colnames(simulation_settings_df) <- c("rho","noise_sd")
num_replicates <- c(50)
cluster_proportion <- c(0.3,0.4,0.3)
###
sim_25by18 <- vector(mode="list", length = nrow(simulation_settings_df))
###
for(i in 1:nrow(simulation_settings_df)){
  time1 <- Sys.time()
  temp_rho <- simulation_settings_df$rho[i]
  temp_noise_sd <- simulation_settings_df$noise_sd[i]
  # correlation matrix size: nrow(M1) * ncol(M1) by nrow(M1) * ncol(M1)
  V <- temp_noise_sd^2 * ar1_cor(ncol(M1), temp_rho)
  # temp_cov <- matrix(NA, nrow = nrow(M1)*ncol(M1), ncol = nrow(M1)*ncol(M1))
  # for(ii in 1:nrow(temp_cov)){
  #   for(jj in 1:ncol(temp_cov)){
  #     temp_ii_2d_index <- c(ii%%ncol(M1)+1, ii%/%ncol(M1)+1)
  #     temp_jj_2d_index <- c(jj%%ncol(M1)+1, jj%/%ncol(M1)+1)
  #     temp_cov[ii,jj] <- temp_noise_sd^2 * temp_rho^(abs(temp_ii_2d_index[1]-temp_jj_2d_index[1])+abs(temp_ii_2d_index[2]-temp_jj_2d_index[2]))
  #   }
  # }
  sim_25by18[[i]] <- foreach(j = 1:num_replicates, .errorhandling = "pass")%dopar%{
    #
    #
    set.seed(j)
    U <- cov2cor(MCMCpack::rwish(nrow(M1)+1,diag(rep(1,nrow(M1)))))
    dat_sim_cluster1 <- mniw::rMNorm(total_sample_size*cluster_proportion[1], 
                                      Lambda = M1, SigmaR = U, SigmaC = V)
    dat_sim_cluster2 <- mniw::rMNorm(total_sample_size*cluster_proportion[2], 
                                     Lambda = M2, SigmaR = U, SigmaC = V)
    dat_sim_cluster3 <- mniw::rMNorm(total_sample_size*cluster_proportion[3], 
                                     Lambda = M3, SigmaR = U, SigmaC = V)
    # dat_sim_cluster1 <- array(0, dim=c(nrow(M1),ncol(M1),total_sample_size*cluster_proportion[1]))
    # for(ii in 1:dim(dat_sim_cluster1)[3]){
    #   set.seed(1234*ii+j)
    #   dat_sim_cluster1[,,ii] <- 
    #   # dat_sim_cluster1[,,ii] <- M1 + matrix(mvtnorm::rmvnorm(1, mean = rep(0,nrow(temp_cov)), sigma = temp_cov),nrow=nrow(M1))
    # }
    # dat_sim_cluster2 <- array(0, dim=c(nrow(M2),ncol(M2),total_sample_size*cluster_proportion[2]))
    # for(ii in 1:dim(dat_sim_cluster2)[3]){
    #   set.seed(12345*ii+j)
    #   dat_sim_cluster2[,,ii] <- M2 + matrix(mvtnorm::rmvnorm(1, mean = rep(0,nrow(temp_cov)), sigma = temp_cov),nrow=nrow(M2))
    # }
    # dat_sim_cluster3 <- array(0, dim=c(nrow(M3),ncol(M3),total_sample_size*cluster_proportion[3]))
    # for(ii in 1:dim(dat_sim_cluster3)[3]){
    #   set.seed(123456*ii+j)
    #   dat_sim_cluster3[,,ii] <- M3 + matrix(mvtnorm::rmvnorm(1, mean = rep(0,nrow(temp_cov)), sigma = temp_cov),nrow=nrow(M3))
    #   # +matrix(rnorm(n=nrow(M3)*ncol(M3),mean=0,sd=0.1), nrow=nrow(M3), ncol=ncol(M3))
    # }
    # combine them
    dat_sim <- abind(dat_sim_cluster1, 
                     dat_sim_cluster2, 
                     dat_sim_cluster3,
                     along = 3)
    dat_sim_vector <- t(apply(dat_sim, c(3), c))
    true_cluster_membership <- rep(c(1,2,3), c(total_sample_size*cluster_proportion[1],
                                               total_sample_size*cluster_proportion[2],
                                               total_sample_size*cluster_proportion[3]))
    # dat_sim <- array(NA, dim = c(nrow(M1),ncol(M1),total_sample_size))
    # run MFM analysis
    temp_time1 <- Sys.time()
    temp_MFM_rlt <- CDMFM_new1_log(data=dat_sim, 
                                   niterations=1000, 
                                   alpha=(1+dim(dat_sim)[1])/2, beta=diag(rep(0.5,dim(dat_sim)[1])), 
                                   psi=(1+dim(dat_sim)[2])/2, rho=diag(rep(0.5,dim(dat_sim)[2])), 
                                   GAMMA=3, 
                                   M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2, 
                                   Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 ), 
                                   Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 ), 
                                   initNClusters=sample(2:10,size=1), 
                                   VN=VN, MLE.initial=TRUE)
    RI_MFM_trace <- rep(NA, 1000)
    ARI_MFM_trace <- rep(NA, 1000)
    for(ii in 1:1000){
      RI_MFM_trace[ii] <- fossil::rand.index(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout)
      ARI_MFM_trace[ii] <- fossil::adj.rand.index(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout)
    }
    # do.call("c", lapply( ,function(x) fossil::rand.index(true_cluster_membership,x$zout) ) )
    temp_MFM_Dahl_rlt <- getDahl(temp_MFM_rlt,burn=600)
    temp_time2 <- Sys.time()
    ## k-means
    k_vec <- 2:10
    temp_kmeans_list <- vector(mode="list", length=length(k_vec))
    temp_kmeans_aic <- rep(NA, length=length(k_vec))
    temp_kmeans_bic <- rep(NA, length=length(k_vec))
    for(k in 1:length(k_vec)){
      # kmeans each row is a data point
      temp_kmeans_list[[k]] <- kmeans(dat_sim_vector, centers = k_vec[k])
      temp_kmeans_aic[k] <- kmeansAIC(temp_kmeans_list[[k]])
      temp_kmeans_bic[k] <- kmeansBIC(temp_kmeans_list[[k]])
    }
    temp_kmeans_rlt_aic <- temp_kmeans_list[which(temp_kmeans_aic==min(temp_kmeans_aic))]
    temp_kmeans_rlt_bic <- temp_kmeans_list[which(temp_kmeans_bic==min(temp_kmeans_bic))]
    temp_kmeans_rlt_3 <- temp_kmeans_list[which(k_vec==3)]
    temp_kmeans_rlt_MFM <- temp_kmeans_list[which(k_vec==max(temp_MFM_Dahl_rlt$zout))]
    temp_time3 <- Sys.time()
    
    temp_specc_rlt_3 <- kernlab::specc(dat_sim_vector, centers=3)
    temp_specc_rlt_MFM <- kernlab::specc(dat_sim_vector, centers=max(temp_MFM_Dahl_rlt$zout))
    temp_time4 <- Sys.time()
    # clustering performance
    ARI_MFM <-  fossil::adj.rand.index(true_cluster_membership,
                                       temp_MFM_Dahl_rlt$zout)
    
    ARI_kmeans_aic <- fossil::adj.rand.index(true_cluster_membership,
                                             temp_kmeans_rlt_aic[[1]]$cluster)
    ARI_kmeans_bic <- fossil::adj.rand.index(true_cluster_membership,
                                             temp_kmeans_rlt_bic[[1]]$cluster)
    ARI_kmeans_3 <- fossil::adj.rand.index(true_cluster_membership,
                                           temp_kmeans_rlt_3[[1]]$cluster)
    ARI_kmeans_MFM <- fossil::adj.rand.index(true_cluster_membership,
                                             temp_kmeans_rlt_MFM[[1]]$cluster)
    ARI_specc_MFM <- fossil::adj.rand.index(true_cluster_membership,
                                        temp_specc_rlt_MFM@.Data)
    ARI_specc_3 <- fossil::adj.rand.index(true_cluster_membership,
                                            temp_specc_rlt_3@.Data)
    #
    RI_MFM <- fossil::rand.index(true_cluster_membership,
                                 temp_MFM_Dahl_rlt$zout)
    RI_kmeans_aic <- fossil::rand.index(true_cluster_membership,
                                        temp_kmeans_rlt_aic[[1]]$cluster)
    RI_kmeans_bic <- fossil::rand.index(true_cluster_membership,
                                        temp_kmeans_rlt_bic[[1]]$cluster)
    RI_kmeans_3 <- fossil::rand.index(true_cluster_membership,
                                      temp_kmeans_rlt_3[[1]]$cluster)
    RI_kmeans_MFM <- fossil::rand.index(true_cluster_membership,
                                             temp_kmeans_rlt_MFM[[1]]$cluster)
    RI_specc_MFM <- fossil::rand.index(true_cluster_membership,
                                            temp_specc_rlt_MFM@.Data)
    RI_specc_3 <- fossil::rand.index(true_cluster_membership,
                                          temp_specc_rlt_3@.Data)
    
    # 
    Khat_MFM <- max(temp_MFM_Dahl_rlt$zout)
    Khat_kmeans_aic <- max(temp_kmeans_rlt_aic[[1]]$cluster)
    Khat_kmeans_bic <- max(temp_kmeans_rlt_bic[[1]]$cluster)
    Khat_kmeans_3 <- 3
    Khat_kmeans_MFM <- Khat_MFM
    Khat_specc_3 <- 3
    Khat_specc_MFM <- Khat_MFM 
    # RMSE V %x% U
    # VU_RMSE_MFM <- NA
    VU_RMSE_MFM <- sqrt(mean(((temp_MFM_Dahl_rlt$Vout %x% temp_MFM_Dahl_rlt$Uout) - (V %x% U))^2))
    # 
    M_RMSE_MFM <- sqrt(mean((temp_MFM_Dahl_rlt$Mout[,,temp_MFM_Dahl_rlt$zout] - dat_sim)^2))
    M_RMSE_kmeans_aic <- sqrt(mean((temp_kmeans_rlt_aic[[1]]$centers[temp_kmeans_rlt_aic[[1]]$cluster,] - dat_sim_vector)^2))
    M_RMSE_kmeans_bic <- sqrt(mean((temp_kmeans_rlt_bic[[1]]$centers[temp_kmeans_rlt_bic[[1]]$cluster,] - dat_sim_vector)^2))
    M_RMSE_kmeans_3 <- sqrt(mean((temp_kmeans_rlt_3[[1]]$centers[temp_kmeans_rlt_3[[1]]$cluster,] - dat_sim_vector)^2))
    M_RMSE_kmeans_MFM <- sqrt(mean((temp_kmeans_rlt_MFM[[1]]$centers[temp_kmeans_rlt_MFM[[1]]$cluster,] - dat_sim_vector)^2))
    M_RMSE_specc_3 <- sqrt(mean((temp_specc_rlt_3@centers[temp_specc_rlt_3@.Data,] - dat_sim_vector)^2))
    M_RMSE_specc_MFM <- sqrt(mean((temp_specc_rlt_MFM@centers[temp_specc_rlt_MFM@.Data,] - dat_sim_vector)^2))
    # 
    temp_rlt <- list(MFM_Dahl_rlt = temp_MFM_Dahl_rlt,
                     kmeans_rlt_aic = temp_kmeans_rlt_aic[[1]],
                     kmeans_rlt_bic = temp_kmeans_rlt_bic[[1]],
                     kmeans_rlt_3 = temp_kmeans_rlt_3[[1]],
                     specc_time = temp_time4-temp_time3,
                     kmeans_time = temp_time3-temp_time2,
                     MFM_Dahl_time = temp_time2-temp_time1,
                     Khat_MFM=Khat_MFM,
                     Khat_kmeans_aic=Khat_kmeans_aic,
                     Khat_kmeans_bic=Khat_kmeans_bic,
                     Khat_kmeans_3=Khat_kmeans_3,
                     Khat_kmeans_MFM=Khat_kmeans_MFM,
                     Khat_specc_3=Khat_specc_3,
                     Khat_specc_MFM=Khat_specc_MFM,
                     ARI_MFM=ARI_MFM,
                     ARI_kmeans_aic=ARI_kmeans_aic,
                     ARI_kmeans_bic=ARI_kmeans_bic,
                     ARI_kmeans_3=ARI_kmeans_3,
                     ARI_kmeans_MFM=ARI_kmeans_MFM,
                     ARI_specc_3=ARI_specc_3,
                     ARI_specc_MFM=ARI_specc_MFM,
                     RI_MFM=RI_MFM,
                     RI_kmeans_aic=RI_kmeans_aic,
                     RI_kmeans_bic=RI_kmeans_bic,
                     RI_kmeans_3=RI_kmeans_3,
                     RI_kmeans_MFM=RI_kmeans_MFM,
                     RI_specc_3=RI_specc_3,
                     RI_specc_MFM=RI_specc_MFM,
                     RI_MFM_trace=RI_MFM_trace,
                     ARI_MFM_trace=ARI_MFM_trace,
                     VU_RMSE_MFM=VU_RMSE_MFM,
                     M_RMSE_MFM=M_RMSE_MFM,
                     M_RMSE_kmeans_aic=M_RMSE_kmeans_aic,
                     M_RMSE_kmeans_bic=M_RMSE_kmeans_bic,
                     M_RMSE_kmeans_3=M_RMSE_kmeans_3,
                     M_RMSE_kmeans_MFM=M_RMSE_kmeans_MFM,
                     M_RMSE_specc_3=M_RMSE_specc_3,
                     M_RMSE_specc_MFM=M_RMSE_specc_MFM)
    temp_rlt
  }
  time2 <- Sys.time()
  # sink("25by18_new_rho_kronecker_diffnoise_sd_status.txt")
  # print(temp_rho)
  # print(temp_noise_sd)
  # print(time2-time1)
  # cat("\n")
  # cat("\n")
  # sink()
}
####################
# MFM_simulation_study_25_18_new_kronecker_high_noise.R
####################
saveRDS(sim_25by18,
        paste0("sim_25by18_new_rho_kronecker_diffnoise_sd.rds"))
# sim_25by18_new_rho_kronecker_diffnoise_sd.rds