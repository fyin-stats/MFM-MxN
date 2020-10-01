#####################
#####################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack",
              "mvtnorm","dplyr","armspp","ggcorrplot","reshape2","ggpubr",
              "cowplot","doParallel","foreach","MixMatrix","mclust","mcclust",
              "fossil","abind","xtable","spatstat","igraph","ggplot2")
# "rjson","ggplot2","grid","inlabru","jpeg",
#  "RCurl","sp","INLA","Orcs","dplyr","gridExtra","ggsci"
# https://cran.r-project.org/web/packages/MixMatrix/vignettes/matrixnormal.html
ipak(packages)
######################
######################
# setwd("/Users/fan/Documents/MFM_matrix/")
######################
sim_10by6_1_noise025 <- readRDS("./sim_10by6_1_noise025_3clusters.rds")
sim_10by6_1 <- readRDS("./sim_10by6_1_3clusters.rds")
# sim_10by6_1_noise025 <- readRDS("./sim_10by6_1_noise025_3clusters_copy2_augmented.rds")
# sim_10by6_1 <- readRDS("./sim_10by6_1_3clusters_augmented.rds")
#######################
signal_strength = 1
# noise_factor <- 1
# M1
M1 <- matrix(0,nrow=10,ncol=6)
M1[5:6,2:5] <- signal_strength
M1[2:9,3:4] <- signal_strength
rownames(M1) <- 1:10
colnames(M1) <- 1:6
# M2 : dumbell
M2 <- matrix(0,nrow=10,ncol=6)
M2[1,2:5] <- signal_strength
M2[2:9,3] <- signal_strength
M2[10,2:5] <- signal_strength
rownames(M2) <- 1:10
colnames(M2) <- 1:6
# M3 : square
M3 <- matrix(0,nrow=10,ncol=6)
M3[4:7,2:5] <- signal_strength
rownames(M3) <- 1:10
colnames(M3) <- 1:6
### Mean matrices
# sample size: 100
jpeg("shape_10by6_simulation_n100_recovered_M1.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[1]][[1]]$MFM_Dahl_rlt$Mout[,,1], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
#
jpeg("shape_10by6_simulation_n100_recovered_M3.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[1]][[1]]$MFM_Dahl_rlt$Mout[,,2], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
#
jpeg("shape_10by6_simulation_n100_recovered_M2.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[1]][[1]]$MFM_Dahl_rlt$Mout[,,3], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
#####################
# sample size: 200
jpeg("shape_10by6_simulation_n200_recovered_M1.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[2]][[1]]$MFM_Dahl_rlt$Mout[,,1], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
#
jpeg("shape_10by6_simulation_n200_recovered_M2.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[2]][[1]]$MFM_Dahl_rlt$Mout[,,2], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
#
jpeg("shape_10by6_simulation_n200_recovered_M3.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[2]][[1]]$MFM_Dahl_rlt$Mout[,,3], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
#
####################
## sample size: 400
jpeg("shape_10by6_simulation_n400_recovered_M2.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[3]][[1]]$MFM_Dahl_rlt$Mout[,,1], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
####################
jpeg("shape_10by6_simulation_n400_recovered_M3.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[3]][[1]]$MFM_Dahl_rlt$Mout[,,2], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
##
jpeg("shape_10by6_simulation_n400_recovered_M1.jpg")
image(x=1:10,y=1:6, z=sim_10by6_1[[3]][[1]]$MFM_Dahl_rlt$Mout[,,3], 
      useRaster=TRUE, col = grey(seq(1, 0, length = 256)),
      xlab="",ylab="")
axis(2, at=seq(0,1, length=6), labels=colnames(M1), lwd=0, pos=0)
axis(3, at=seq(0,1, length=10), labels=rownames(M1), lwd=0, pos=0)
dev.off()
##########
### Khat estimation by MFM
# high noise
Khat_MFM_n100 <- do.call("c", lapply(sim_10by6_1[[1]], function(x) max(x$Khat_MFM)))
Khat_MFM_n200 <- do.call("c", lapply(sim_10by6_1[[2]], function(x) max(x$Khat_MFM)))
Khat_MFM_n400 <- do.call("c", lapply(sim_10by6_1[[3]], function(x) max(x$Khat_MFM)))
###
xtable::xtable(rbind(table(factor(Khat_MFM_n100, levels=2:4)),
                     table(factor(Khat_MFM_n200,levels=2:4)),
                     table(factor(Khat_MFM_n400,levels=2:4))))

###########
# low noise
Khat_MFM_n100_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[1]], function(x) max(x$Khat_MFM)))
Khat_MFM_n200_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[2]], function(x) max(x$Khat_MFM)))
Khat_MFM_n400_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[3]], function(x) max(x$Khat_MFM)))
###
xtable::xtable(rbind(table(factor(Khat_MFM_n100_noise025, levels=2:4)),
                     table(factor(Khat_MFM_n200_noise025,levels=2:4)),
                     table(factor(Khat_MFM_n400_noise025,levels=2:4))))
########################
########################
# Rand index
########################
# high noise
RI_MFM_n100 <- do.call("c", lapply(sim_10by6_1[[1]], function(x) x$RI_MFM))
RI_MFM_n200 <- do.call("c", lapply(sim_10by6_1[[2]], function(x) x$RI_MFM))
RI_MFM_n400 <- do.call("c", lapply(sim_10by6_1[[3]], function(x) x$RI_MFM))
#
RI_kmeans_n100 <- do.call("c", lapply(sim_10by6_1[[1]], function(x) x$RI_kmeans_MFM))
RI_kmeans_n200 <- do.call("c", lapply(sim_10by6_1[[2]], function(x) x$RI_kmeans_MFM))
RI_kmeans_n400 <- do.call("c", lapply(sim_10by6_1[[3]], function(x) x$RI_kmeans_MFM))
# 
RI_spectral_n100 <- do.call("c", lapply(sim_10by6_1[[1]], function(x) x$RI_specc_MFM))
RI_spectral_n200 <- do.call("c", lapply(sim_10by6_1[[2]], function(x) x$RI_specc_MFM))
RI_spectral_n400 <- do.call("c", lapply(sim_10by6_1[[3]], function(x) x$RI_specc_MFM))

# low noise
RI_MFM_n100_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[1]], function(x) x$RI_MFM))
RI_MFM_n200_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[2]], function(x) x$RI_MFM))
RI_MFM_n400_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[3]], function(x) x$RI_MFM))
#
RI_kmeans_n100_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[1]], function(x) x$RI_kmeans_MFM))
RI_kmeans_n200_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[2]], function(x) x$RI_kmeans_MFM))
RI_kmeans_n400_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[3]], function(x) x$RI_kmeans_MFM))
# 
RI_spectral_n100_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[1]], function(x) x$RI_specc_MFM))
RI_spectral_n200_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[2]], function(x) x$RI_specc_MFM))
RI_spectral_n400_noise025 <- do.call("c", lapply(sim_10by6_1_noise025[[3]], function(x) x$RI_specc_MFM))

###################
#
RI_summary_table <- rbind(c(mean(RI_MFM_n100), mean(RI_kmeans_n100), mean(RI_spectral_n100), mean(RI_MFM_n100_noise025), mean(RI_kmeans_n100_noise025), mean(RI_spectral_n100_noise025)),
      c(mean(RI_MFM_n200), mean(RI_kmeans_n200), mean(RI_spectral_n200), mean(RI_MFM_n200_noise025), mean(RI_kmeans_n200_noise025), mean(RI_spectral_n200_noise025)),
      c(mean(RI_MFM_n400), mean(RI_kmeans_n400), mean(RI_spectral_n400), mean(RI_MFM_n400_noise025), mean(RI_kmeans_n400_noise025), mean(RI_spectral_n400_noise025)))
colnames(RI_summary_table) <- c("MFM","K-means","Spectral","MFM","K-means","Spectral")
rownames(RI_summary_table) <- c(100,200,400)

###########
#RMSE for V \otimes U
###########
# sim_10by6_1[[1]][[1]]$VU_RMSE_MFM
VU_RMSE_MFM_10by6_n100 <- do.call("c",lapply(sim_10by6_1[[1]], function(x) x$VU_RMSE_MFM))
VU_RMSE_MFM_10by6_n200 <- do.call("c",lapply(sim_10by6_1[[2]], function(x) x$VU_RMSE_MFM))
VU_RMSE_MFM_10by6_n400 <- do.call("c",lapply(sim_10by6_1[[3]], function(x) x$VU_RMSE_MFM))
VU_RMSE_MFM_10by6_df <- data.frame(RMSE = c(VU_RMSE_MFM_10by6_n100,VU_RMSE_MFM_10by6_n200,VU_RMSE_MFM_10by6_n400),
                                   n = factor(c(rep(100,100),rep(200,100),rep(400,100))), noise = "high noise")
# 
VU_RMSE_MFM_10by6_n100_noise025 <- do.call("c",lapply(sim_10by6_1_noise025[[1]], function(x) x$VU_RMSE_MFM))
VU_RMSE_MFM_10by6_n200_noise025 <- do.call("c",lapply(sim_10by6_1_noise025[[2]], function(x) x$VU_RMSE_MFM))
VU_RMSE_MFM_10by6_n400_noise025 <- do.call("c",lapply(sim_10by6_1_noise025[[3]], function(x) x$VU_RMSE_MFM))
VU_RMSE_MFM_10by6_noise025_df <- data.frame(RMSE = c(VU_RMSE_MFM_10by6_n100_noise025,VU_RMSE_MFM_10by6_n200_noise025,
                                            VU_RMSE_MFM_10by6_n400_noise025),
                                   n = factor(c(rep(100,100),rep(200,100),rep(400,100))), noise = "low noise")
# combine them together
VU_RMSE_MFM_10by6_all_df <- rbind(VU_RMSE_MFM_10by6_df, VU_RMSE_MFM_10by6_noise025_df)
# 
VU_RMSE_MFM_10by6_all_df %>% ggplot(aes(x=RMSE)) + geom_histogram() + facet_wrap(~n+noise,nrow=3) + ylab("") + theme(text = element_text(size=20),
                                                                                                                       axis.text.x = element_text(angle = 45, hjust=1),
                                                                                                                       strip.text = element_text(size = 20)) 

###########
## RI trajectory
sim_10by6_1_noise025 <- readRDS("./sim_10by6_1_noise025_3clusters.rds")
sim_10by6_1 <- readRDS("./sim_10by6_1_3clusters.rds")
###########
### n=100, 200, 400, large noise
###########
RI_trace_10by6_n100_df <- data.frame(RI=c(do.call("cbind",lapply(sim_10by6_1[[1]], function(x) x$RI_MFM_trace))),
                                     iter = rep(1:length(sim_10by6_1[[1]][[1]]$RI_MFM_trace), rep(length(sim_10by6_1[[1]]))), 
                                     noise_level = "1", n = 100)
RI_trace_10by6_n200_df <- data.frame(RI=c(do.call("cbind",lapply(sim_10by6_1[[2]], function(x) x$RI_MFM_trace))),
                                     iter = rep(1:length(sim_10by6_1[[2]][[1]]$RI_MFM_trace), rep(length(sim_10by6_1[[2]]))), 
                                     noise_level = "1", n = 200)
RI_trace_10by6_n400_df <- data.frame(RI=c(do.call("cbind",lapply(sim_10by6_1[[3]], function(x) x$RI_MFM_trace))),
                                     iter = rep(1:length(sim_10by6_1[[3]][[1]]$RI_MFM_trace), rep(length(sim_10by6_1[[3]]))), 
                                     noise_level = "1", n = 400)
##########
### n=100, 200, 400, small noise
##########
RI_trace_10by6_n100_noise025_df <- data.frame(RI=c(do.call("cbind",lapply(sim_10by6_1_noise025[[1]], function(x) x$RI_MFM_trace))),
                                     iter = rep(1:length(sim_10by6_1_noise025[[1]][[1]]$RI_MFM_trace), rep(length(sim_10by6_1_noise025[[1]]))), 
                                     noise_level = "0.25", n = 100)
RI_trace_10by6_n200_noise025_df <- data.frame(RI=c(do.call("cbind",lapply(sim_10by6_1_noise025[[2]], function(x) x$RI_MFM_trace))),
                                     iter = rep(1:length(sim_10by6_1_noise025[[2]][[1]]$RI_MFM_trace), rep(length(sim_10by6_1_noise025[[2]]))), 
                                     noise_level = "0.25", n = 200)
RI_trace_10by6_n400_noise025_df <- data.frame(RI=c(do.call("cbind",lapply(sim_10by6_1_noise025[[3]], function(x) x$RI_MFM_trace))),
                                     iter = rep(1:length(sim_10by6_1_noise025[[3]][[1]]$RI_MFM_trace), rep(length(sim_10by6_1_noise025[[3]]))), 
                                     noise_level = "0.25", n = 400)
##########
##########
RI_trace_10by6_df <- data.frame(rbind(RI_trace_10by6_n100_df,
                           RI_trace_10by6_n200_df,
                           RI_trace_10by6_n400_df,
                           RI_trace_10by6_n100_noise025_df,
                           RI_trace_10by6_n200_noise025_df,
                           RI_trace_10by6_n400_noise025_df))

RI_trace_10by6_df_summary <- RI_trace_10by6_df %>% group_by(noise_level, n, iter) %>% summarise(mean_RI = mean(RI))
# ri_trace_df_summary <- data.frame(ri_trace_df %>% group_by(method, setting, iter) %>% 
#                                     summarise(mean_RI = mean(RI)))
###
p_RI_trace <- ggplot(data=RI_trace_10by6_df, aes(x=iter, y=RI)) + geom_line(colour="darkgrey", size=0.1, alpha=0.2) + 
  facet_wrap(~noise_level+n, ncol=3) + geom_line(data=RI_trace_10by6_df_summary, 
                                                 aes(x=iter,y=mean_RI), colour="red", size=0.25) 

pdf("RI_trace_10by6.pdf")
print(p_RI_trace)
dev.off()
############################################################
############################################################

##################################################################################################################
##########################################
#### 25 by 18 results summary
##########################################
##############################
sim_25by18_new_rho_kronecker_high_noise_sd <- readRDS("sim_25by18_new_rho_kronecker_diffnoise_sd.rds")
###############################################
total_sample_size <- c(200)
rho <- c(0.9,0.6,0.3)
noise_sd <- c(1.5,1,0.5)
simulation_settings_df <- data.frame(expand.grid(rho, noise_sd))
colnames(simulation_settings_df) <- c("rho","noise_sd")
# simulation_settings_df <- data.frame(expand.grid(rho, noise_sd))
# colnames(simulation_settings_df) <- c("rho","noise_sd")
sim_25by18_new_rho_noise_sd <- sim_25by18_new_rho_kronecker_high_noise_sd
#################################################
simulation_settings_RI_df <- data.frame(RI_MFM=NA, RI_kmeans_MFM = NA, RI_spectral_MFM=NA)
for(i in 1:nrow(simulation_settings_df)){
  simulation_settings_RI_df <- rbind(simulation_settings_RI_df,
                                     data.frame(RI_MFM=mean(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[i]], function(x) x$RI_MFM))),
                                                RI_kmeans_MFM=mean(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[i]], function(x) x$RI_kmeans_MFM))),
                                                RI_spectral_MFM=mean(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[i]], function(x) x$RI_specc_MFM)))) )
}
simulation_settings_RI_df <- simulation_settings_RI_df[complete.cases(simulation_settings_RI_df),]
#
xtable(cbind(simulation_settings_df, simulation_settings_RI_df), digits = 3)
#########

##########
### Khat under different settings
##########
xtable::xtable(rbind(table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[1]], function(x) x$Khat_MFM)), levels=2:4)),
table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[2]], function(x) x$Khat_MFM)), levels=2:4)),
table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[3]], function(x) x$Khat_MFM)), levels=2:4)),
table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[4]], function(x) x$Khat_MFM)), levels=2:4)),
table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[5]], function(x) x$Khat_MFM)), levels=2:4)),
table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[6]], function(x) x$Khat_MFM)), levels=2:4)),
table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[7]], function(x) x$Khat_MFM)), levels=2:4)),
table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[8]], function(x) x$Khat_MFM)), levels=2:4)),
table(factor(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[9]], function(x) x$Khat_MFM)), levels=2:4))))

###########
## RI traceplot
RI_trace_25by18_df <- data.frame(RI=NA,iter=NA, rho=NA,noise_sd=NA)
############
for(i in 1:nrow(simulation_settings_df)){
  RI_trace_25by18_df <- rbind(RI_trace_25by18_df,
                           data.frame(RI=do.call("c", lapply(sim_25by18_new_rho_noise_sd[[i]], 
                                                             function(x) x$RI_MFM_trace)),
                                      iter=rep(1:length(sim_25by18_new_rho_noise_sd[[i]][[1]]$RI_MFM_trace), 
                                               length(do.call("c", lapply(sim_25by18_new_rho_noise_sd[[i]], 
                                                                   function(x) x$RI_MFM_trace)))/length(sim_25by18_new_rho_noise_sd[[i]][[1]]$RI_MFM_trace)),
                                      rho=simulation_settings_df$rho[i],
                                      noise_sd=simulation_settings_df$noise_sd[i]) 
                           )
}

##############
RI_trace_25by18_df <- RI_trace_25by18_df[complete.cases(RI_trace_25by18_df),]
##############
RI_trace_25by18_df_summary <- RI_trace_25by18_df %>% group_by(rho, noise_sd, iter) %>% summarise(mean_RI=mean(RI))
##############
p_RI_trace_25by18 <- ggplot(data=RI_trace_25by18_df, aes(x=iter, y=RI)) + geom_line(colour="darkgrey", size=0.1,alpha=0.2) + 
  facet_wrap(~noise_sd+rho, ncol=3) + geom_line(data=RI_trace_25by18_df_summary, 
                                                 aes(x=iter,y=mean_RI), colour="red", size=0.25) 

pdf("RI_trace_25by18.pdf")
print(p_RI_trace_25by18)
dev.off()
###############