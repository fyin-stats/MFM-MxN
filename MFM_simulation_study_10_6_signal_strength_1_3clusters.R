###################################
############ MFM matrix normal simulation studies
###################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack",
              "mvtnorm","dplyr","armspp","ggcorrplot","reshape2","ggpubr",
              "cowplot","doParallel","foreach","MixMatrix","abind",
              "mclust", "mcclust","fossil","kernlab")
# "rjson","ggplot2","grid","inlabru","jpeg",
#  "RCurl","sp","INLA","Orcs","dplyr","gridExtra","ggsci"
# https://cran.r-project.org/web/packages/MixMatrix/vignettes/matrixnormal.html
ipak(packages)
#
source("MFM_MN.R")
registerDoParallel(cores=50)
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
###
# helper function low rank matrix 
# p number of rows
# q number of columns
# s: sparsity level (proportion of zero elements)
# R: rank
# low_rank_matrix <- function(p,q,s=0.5,R=min(p,q)){
#   prob <- sqrt(1 - (1-s)^(1/R))
#   M <- matrix(rbinom(n=p*q,size=1, prob = prob),
#                nrow=p, ncol=q)
#   return(M)
# }
####################################
### input -- sparsity level: s; rank of B: R; size of X: p1 and p2; sample size: n; beta: regression parameters for Z
### s = 
# s=0.5
# p1 = 10
# p2 = 5
# R = min(p1,p2)
# n = c(100,200,300)
#######
# prob <- sqrt(1 - (1-s)^(1/R) )
# B1 <- matrix(rbinom(n=p1*R,size=1, prob = prob),
#              nrow=p1, ncol=R) #p1 by R
# B2 <- matrix(rbinom(n=p2*R,size=2, prob=prob),
#              nrow=p2, ncol=R)
# 
# B1
# rownames(B1) <- 1:10
# colnames(B1) <- 1:6
# B1_long <- melt(B1)
# p1 <- ggplot(B1_long, aes(x = Var1, y = Var2)) +
#   geom_raster(aes(fill=value)) +
#   scale_fill_gradient(low="grey90", high="red") +
#   labs(x="X", y="Y", title=paste("")) +
#   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
#                      axis.text.y=element_text(size=9),
#                      plot.title=element_text(size=11))

############# simulation setting 1:
# Data: 10*6
# Cov1 : AR(1) with rho = 0.75
# Cov2 : a draw from standard wishart
# sparsity level fixed at 0.5
# Mean : 4 clusters: cross, T, dumbell, square
signal_strength = 1
noise_factor <- 1
# M1 : cross
M1 <- matrix(0,nrow=10,ncol=6)
M1[5:6,2:5] <- signal_strength
M1[2:9,3:4] <- signal_strength
rownames(M1) <- 1:10
colnames(M1) <- 1:6
M1_long <- melt(M1)
p1 <- ggplot(M1_long, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="white", high="black") +
  labs(x="", y="", title=paste("")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) +  guides(fill=FALSE)

# M2 : dumbell
M2 <- matrix(0,nrow=10,ncol=6)
M2[1,2:5] <- signal_strength
M2[2:9,3] <- signal_strength
M2[10,2:5] <- signal_strength
rownames(M2) <- 1:10
colnames(M2) <- 1:6
M2_long <- melt(M2)
p2 <- ggplot(M2_long, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="white", high="black") +
  labs(x="", y="", title=paste("")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) +  guides(fill=FALSE)
# square
M3 <- matrix(0,nrow=10,ncol=6)
M3[4:7,2:5] <- signal_strength
rownames(M3) <- 1:10
colnames(M3) <- 1:6
M3_long <- melt(M3)
p3 <- ggplot(M3_long, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="white", high="black") +
  labs(x="", y="", title=paste("")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) + guides(fill=FALSE)
#####
# random matrix with sparsity 0.25
# set.seed(0)
# M4 <- matrix(rbinom(60,1,prob=0.3)*signal_strength,nrow=10,ncol=6)
# rownames(M4) <- 1:10
# colnames(M4) <- 1:6
# M4_long <- melt(M4)
# p4 <- ggplot(M4_long, aes(x = Var1, y = Var2)) +
#   geom_raster(aes(fill=value)) +
#   scale_fill_gradient(low="white", high="red") +
#   labs(x="X", y="Y", title=paste("")) +
#   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
#                      axis.text.y=element_text(size=9),
#                      plot.title=element_text(size=11))
###
# jpeg("shape_10by6_simulation.jpg")
# ggarrange(p1,p2,p3,nrow=1,ncol=3)
# dev.off()
# 
# ###
# par(frow=c(1,3))
# https://stackoverflow.com/questions/27406508/r-axis-label-in-image
jpeg("shape_10by6_simulation_M1.jpg")
image(x=1:10,y=1:6, z=M1, useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
# axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
# axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
jpeg("shape_10by6_simulation_M2.jpg")
image(x=1:10,y=1:6, z=M2, useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
dev.off()
jpeg("shape_10by6_simulation_M3.jpg")
image(x=1:10,y=1:6, z=M3, useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
dev.off()
###
V <- noise_factor*ar1_cor(ncol(M1),0.9)
set.seed(0)
U <- cov2cor(MCMCpack::rwish(nrow(M1)+1,diag(rep(1,nrow(M1)))))
### 0.2, 0.2, 0.3, 0.3
### sample size: 100, 200, 300
### number of replicates: 200
total_sample_size <- c(100,200,400)
num_replicates <- c(100)
cluster_proportion <- c(0.3,0.3,0.4)
# cluster_proportion <- c(0.25,0.25,0.25,0.25)
###
sim_10by6 <- vector(mode="list", length = length(total_sample_size))
###
for(i in 1:length(total_sample_size)){
  sim_10by6[[i]] <- foreach(j = 1:num_replicates, .errorhandling = c("pass"))%dopar%{
    # 
    set.seed(j+1234*i)
    dat_sim_cluster1 <-  mniw::rMNorm(total_sample_size[i]*cluster_proportion[1], 
                                      Lambda = M1, 
                                      SigmaR = U, SigmaC = V)
    dat_sim_cluster2 <-  mniw::rMNorm(total_sample_size[i]*cluster_proportion[2], 
                                      Lambda = M2, 
                                      SigmaR = U, SigmaC = V)
    dat_sim_cluster3 <-  mniw::rMNorm(total_sample_size[i]*cluster_proportion[3], 
                                      Lambda = M3, 
                                      SigmaR = U, SigmaC = V)
    # combine them
    dat_sim <- abind(dat_sim_cluster1, dat_sim_cluster2, dat_sim_cluster3,
                     along = 3)
    dat_sim_vector <- t(apply(dat_sim, c(3), c))
    true_cluster_membership <- rep(c(1,2,3), c(total_sample_size[i]*cluster_proportion[1],
                                                 total_sample_size[i]*cluster_proportion[2],
                                                 total_sample_size[i]*cluster_proportion[3]))
    # dat_sim <- array(NA, dim = c(nrow(M1),ncol(M1),total_sample_size))
    # run MFM analysis
    temp_time1 <- Sys.time()
    temp_MFM_rlt <- CDMFM_new1_log(data=dat_sim, 
                                   niterations=1500, 
                                   alpha=(1+dim(dat_sim)[1])/2, beta=diag(rep(0.5,dim(dat_sim)[1])), 
                                   psi=(1+dim(dat_sim)[2])/2, rho=diag(rep(0.5,dim(dat_sim)[2])), 
                                   GAMMA=3, 
                                   M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2, 
                                   Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 ), 
                                   Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 ), 
                                   initNClusters=sample(2:10,size=1), 
                                   VN=VN, MLE.initial=TRUE)
    RI_MFM_trace <- rep(NA, 1500)
    ARI_MFM_trace <- rep(NA, 1500)
    for(ii in 1:1500){
      RI_MFM_trace[ii] <- fossil::rand.index(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout)
      ARI_MFM_trace[ii] <- fossil::adj.rand.index(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout)
    }
    # do.call("c", lapply( ,function(x) fossil::rand.index(true_cluster_membership,x$zout) ) )
    temp_MFM_Dahl_rlt <- getDahl(temp_MFM_rlt,burn=1000)
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
}
#####################################
saveRDS(sim_10by6,
        paste0("sim_10by6_",signal_strength,"_3clusters",".rds"))
#####################################
#####################################