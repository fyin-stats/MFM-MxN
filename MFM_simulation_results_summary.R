########################
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
#######
#######
library(stringi)
library(tidyr)
library(tibble)
library(xtable)
library(ggplot2)
# true intensity surface
M_list_3by5 <- readRDS("./M_3by5_list.rds")
# 
grid_3by5 <- expand.grid(x=1:3, y=1:5)
#
M_df_3by5 <- data.frame()
#
for(i in 1:length(M_list_3by5)){
  # 
  M_df_3by5 <- rbind(M_df_3by5,
                     data.frame(grid_3by5, 
                                value = c(M_list_3by5[[i]]), 
                                cluster_id = i) )
}
#######
M_df_3by5 %>% ggplot(aes(x = x, y = y, fill = value)) + 
  geom_tile() + facet_wrap(.~cluster_id) + theme(legend.position = "bottom",
                                                 axis.text=element_text(size=12),
                                                 axis.title=element_text(size=14,face="bold"),
                                                 strip.text = element_text(size = 12)) +
  scale_fill_distiller(palette='RdBu', name = "")
#######
# sim_3by5_1 <- readRDS("./simulation_results/sim_3by5_1.rds")
# simulation_settings_3by5 <- data.frame(expand.grid(signal_strength = 1, 
#                                               noise_factor = c(1), 
#                                               cluster_num = c(3),
#                                               total_sample_size = c(500))) %>% rowid_to_column("setting_id")
simulation_settings_3by5 <- data.frame(expand.grid(signal_strength = 1, 
                                              noise_factor = c(0.5, 1, 1.5), 
                                              cluster_num = c(3),
                                              total_sample_size = c(500))) %>% rowid_to_column("setting_id")

#
Khat_rlt_df_3by5 <- data.frame()
time_rlt_df_3by5 <- data.frame()
RMSE_VU_rlt_df_3by5 <- data.frame()
################
## 3 by 5
# 1:nrow(simulation_settings_3by5)
for(i in 2){
  #
  assign(paste0("sim_3by5_", i), readRDS( paste0("./simulation_results/sim_3by5_", i, ".rds")) )
  temp_rlt <- get(paste0("sim_3by5_", i))
  #
  temp_true_K <- simulation_settings_3by5$cluster_num[i]
  #
  # temp_rlt_df <- data.frame(Khat_DP = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
  #                           Khat_MFM = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
  #                           Khat_kmeans = do.call("c", lapply(temp_rlt, function(x) x$Khat_kmeans_MFM)),
  #                           Khat_specc = do.call("c", lapply(temp_rlt, function(x) x$Khat_specc_MFM)) )
  temp_rlt_df <- data.frame(ARI_DP = do.call("c", lapply(temp_rlt, function(x) x$ARI_DP)),
                            ARI_MFM = do.call("c", lapply(temp_rlt, function(x) x$ARI_MFM)),
                            ARI_kmeans = do.call("c", lapply(temp_rlt, function(x) x$ARI_kmeans_MFM)),
                            ARI_specc = do.call("c", lapply(temp_rlt, function(x) x$ARI_specc_MFM)),
                            ARI_mvnormalmixEM = do.call("c", lapply(temp_rlt, function(x) x$ARI_mvnormalmixEM_MFM)),
                            ARI_kcentroid_nuclear = do.call("c", lapply(temp_rlt, function(x) x$ARI_kcentroid_nuclear_MFM)),
                            ARI_kcentroid_spectral = do.call("c", lapply(temp_rlt, function(x) x$ARI_kcentroid_spectral_MFM)),
                            Khat_DP = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
                            Khat_MFM = do.call("c", lapply(temp_rlt, function(x) x$Khat_MFM)),
                            K_truth = temp_true_K,
                            setting_id = i) %>% rowid_to_column(var = "replicate_id")
  
  #
  Khat_rlt_df_3by5 <- rbind(Khat_rlt_df_3by5,
                             temp_rlt_df)
  
  #
  time_rlt_df_3by5 <- rbind(time_rlt_df_3by5,
                            data.frame(time_DP = do.call("c", lapply(temp_rlt, function(x) x$MFM_DP_time)),
                                       time_MFM = do.call("c", lapply(temp_rlt, function(x) x$MFM_Dahl_time)),
                                       setting_id = i) %>% rowid_to_column(var = "replicate_id") ) 
  
  #
  RMSE_VU_rlt_df_3by5 <- rbind(RMSE_VU_rlt_df_3by5,
                                data.frame(RMSE_VU_DP = do.call("c", lapply(temp_rlt, function(x) x$VU_RMSE_DP)),
                                           RMSE_VU_MFM = do.call("c", lapply(temp_rlt, function(x) x$VU_RMSE_MFM)),
                                           setting_id = i) %>% rowid_to_column(var = "replicate_id")) 
  
}
###################
###################
###################
###################
Khat_rlt_df_3by5 %>% group_by(setting_id) %>% summarise(n = n())
# RMSE VU
RMSE_VU_rlt_df_3by5 %>% dplyr::select(replicate_id, setting_id, contains("RMSE_VU")) %>% 
          pivot_longer(contains("RMSE_VU"), names_to = "method", values_to = "RMSE") %>% 
            inner_join(simulation_settings_3by5, by = c("setting_id" = "setting_id")) %>% 
  mutate(method = stri_sub(method,9,-1)) %>% ggplot(aes(x = RMSE, color = method)) +
  geom_density() + theme(legend.position = "bottom",
                         axis.text = element_text(size=12),
                         axis.title = element_text(size=14,face="bold"),
                         legend.text = element_text(size = 12)) 

# geom_histogram(position = "identity", alpha = 0.2)

###############
Khat_rlt_df_3by5 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
  pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
  group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
  inner_join(simulation_settings_3by5, by = c("setting_id" = "setting_id")) %>%
  mutate(method = stri_sub(method,6,-1)) %>% mutate(`rho factor` = factor(rho_factor),
                                                    `noise factor` = noise_factor,
                                                    `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
             y = Khat_accuracy, 
             color = method,
             shape = `rho factor`) ) +
  geom_point(size=2.5, alpha = 0.5, position = position_jitter(width = 0.1, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Khat accuracy") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom") 

# lapply(sim_3by5_1, 
#        function(x) x$MFM_Dahl_rlt)
# 
# ##
Khat_rlt_df_3by5 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% 
  group_by(setting_id, method) %>% summarise(mean_ARI = mean(ARI, na.rm = T)) %>% 
  mutate(method = stri_sub(method,5,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>% 
  unite("sample size, rho factor", total_sample_size:rho_factor, remove = FALSE,
        sep = ",") %>%
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>%
  ggplot(aes(x = `sample size, rho factor`, 
             y = mean_ARI, 
             shape = method,
             color = method) ) +
  geom_point(size=2.5, alpha = 0.5, position = position_jitter(width = 0.1, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Mean adjusted Rand Index") + xlab("Sample size, rho factor") +
  theme_bw() +
  theme(legend.position = "bottom") 

##
Khat_rlt_df_3by5 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% 
  group_by(setting_id, method) %>% summarise(mean_ARI = mean(ARI, na.rm = T),
                                             sd_ARI = sd(ARI, na.rm=T)) %>% 
  mutate(method = stri_sub(method,5,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>% 
  dplyr::select(method, mean_ARI, sd_ARI) %>%
  xtable(digits = 3)

# 3 by 5
time_rlt_df_3by5 %>% dplyr::select(replicate_id, setting_id, contains("time")) %>%
  pivot_longer(contains("time"), names_to = "method", values_to = "time") %>% 
  group_by(setting_id, method) %>% summarise(mean_time = mean(time)) %>% 
  mutate(method = stri_sub(method,6,-1)) %>%
  inner_join(simulation_settings_3by5, by = c("setting_id" = "setting_id")) %>%
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
             y = mean_time, 
             color = method,
             shape = `rho factor`) ) +
  geom_point(size=3.5, alpha = 0.75, position = position_jitter(width = 0.1, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Mean computation time (mins)") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14)) 

#
time_rlt_df_3by5 %>% dplyr::select(replicate_id, setting_id, contains("time")) %>%
  pivot_longer(contains("time"), names_to = "method", values_to = "time") %>% 
  group_by(setting_id, method) %>% summarise(mean_time = mean(time),
                                             median_time = median(time),
                                             SD_time = sd(time)) %>% 
  mutate(method = stri_sub(method,6,-1)) %>%
  inner_join(simulation_settings_3by5, by = c("setting_id" = "setting_id")) %>%
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>% dplyr::select(method, mean_time,
                                                           median_time, SD_time) %>%
  xtable(digits = 1)

# sim_3by5_1[[1]]$ARI_MFM
# ##
# do.call("c", lapply(sim_3by5_1, function(x) x$ARI_MFM)) %>% mean()
# # do.call("c", lapply(sim_3by5_1, function(x) x$ARI_MFM)) %>% mean()
# ## 
# do.call("c", lapply(sim_3by5_1, function(x) x$ARI_DP)) %>% mean()
# do.call("c", lapply(sim_3by5_1, function(x) x$Khat_DP)) %>% table()

########################################
########################################
## true values
M_list_10by6 <- readRDS("./M_10by6_list.rds")
grid_10by6 <- expand.grid(x=1:10, y=1:6)
#
M_df_10by6 <- data.frame()
#
for(i in 1:length(M_list_10by6)){
  # 
  M_df_10by6 <- rbind(M_df_10by6,
                     data.frame(grid_10by6, 
                                value = c(M_list_10by6[[i]]), 
                                cluster_id = i) )
}
####### 
M_df_10by6 %>% ggplot(aes(x = x, y = y, fill = value)) + 
  geom_tile() + facet_wrap(.~cluster_id) + 
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 12)) + 
  scale_fill_gradient(
    low = "white",
    high = "black"
  )
##
# Khat, ARI
# simulation_settings_10by6 <- data.frame(expand.grid(signal_strength = 1, noise_factor = c(1.5, 1, 0.5), cluster_num = c(3,5),
#                                               total_sample_size = c(200, 400, 100))) %>% rowid_to_column("setting_id")
simulation_settings_10by6 <- data.frame(expand.grid(signal_strength = 1, 
                                              noise_factor = c(1.5, 1, 0.5), 
                                              total_sample_size = c(100,200,300),
                                              rho_factor = c(0.9, 0.3),
                                              cluster_num = c(3,6))) %>% rowid_to_column("setting_id")

# data.frame(expand.grid(signal_strength = 1, 
#                        noise_factor = c(1.5, 1, 0.5), 
#                        total_sample_size = c(200, 500),
#                        rho_factor = c(0.9, 0.3),
#                        cluster_num = c(3,5))) %>% rowid_to_column("setting_id")

#
Khat_rlt_df_10by6 <- data.frame()
################
# computational time
time_rlt_df_10by6 <- data.frame()
RMSE_VU_rlt_df_10by6 <- data.frame()
## 10 by 6
for(i in c(1:nrow(simulation_settings_10by6))){
  #
  assign(paste0("sim_10by6_", i), readRDS( paste0("./simulation_results/sim_10by6_", i, ".rds")) )
  temp_rlt <- get(paste0("sim_10by6_", i))
  #
  temp_true_K <- simulation_settings_10by6$cluster_num[i]
  #
  # temp_rlt_df <- data.frame(Khat_DP = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
  #                           Khat_MFM = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
  #                           Khat_kmeans = do.call("c", lapply(temp_rlt, function(x) x$Khat_kmeans_MFM)),
  #                           Khat_specc = do.call("c", lapply(temp_rlt, function(x) x$Khat_specc_MFM)) )
  temp_rlt_df <- data.frame(ARI_DP = do.call("c", lapply(temp_rlt, function(x) x$ARI_DP)),
                            ARI_MFM = do.call("c", lapply(temp_rlt, function(x) x$ARI_MFM)),
                            ARI_kmeans = do.call("c", lapply(temp_rlt, function(x) x$ARI_kmeans_oracle)),
                            ARI_specc = do.call("c", lapply(temp_rlt, function(x) x$ARI_specc_oracle)),
                            ARI_kcentroid_nuclear = do.call("c", lapply(temp_rlt, function(x) x$ARI_kcentroid_nuclear_oracle)),
                            ARI_kcentroid_spectral = do.call("c", lapply(temp_rlt, function(x) x$ARI_kcentroid_spectral_oracle)),
                            Khat_DP = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
                            Khat_MFM = do.call("c", lapply(temp_rlt, function(x) x$Khat_MFM)),
                            K_truth = temp_true_K,
                            setting_id = i) %>% rowid_to_column(var = "replicate_id")
  
  #
  Khat_rlt_df_10by6 <- rbind(Khat_rlt_df_10by6,
                             temp_rlt_df)
  
  # 
  time_rlt_df_10by6 <- rbind(time_rlt_df_10by6,
                             data.frame(time_DP = do.call("c", lapply(temp_rlt, function(x) x$MFM_DP_time)),
                                        time_MFM = do.call("c", lapply(temp_rlt, function(x) x$MFM_Dahl_time)),
                                        setting_id = i) %>% rowid_to_column(var = "replicate_id") ) 
  
  # 
  RMSE_VU_rlt_df_10by6 <- rbind(RMSE_VU_rlt_df_10by6,
                                data.frame(RMSE_VU_DP = do.call("c", lapply(temp_rlt, function(x) x$VU_RMSE_DP)),
                                           RMSE_VU_MFM = do.call("c", lapply(temp_rlt, function(x) x$VU_RMSE_MFM)),
                                           setting_id = i) %>% rowid_to_column(var = "replicate_id")) 

}
########
Khat_rlt_df_10by6 %>% group_by(setting_id) %>% summarise(n = n()) %>% xtable()
########
########
# where are the NAs?
########
########
ARI_10by6 <- Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% pull(ARI) 
which(is.na(ARI_10by6) == TRUE)
sum(is.na(ARI_10by6))
# 7: 91
# 8: 77
# Khat accuracy
# Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
#   pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
#   group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
#   inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
#   mutate(method = stri_sub(method,6,-1)) %>% xtable()
###############
###############
# noise factor 1.5
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
  pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
  group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  mutate(method = stri_sub(method,6,-1)) %>% filter(noise_factor == 1.5)

# noise factor 1
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
  pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
  group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  mutate(method = stri_sub(method,6,-1)) %>% filter(noise_factor == 1)

# noise factor 0.5
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
  pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
  group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  mutate(method = stri_sub(method,6,-1)) %>% filter(noise_factor == 0.5)

###############
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
  pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
  group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  mutate(method = stri_sub(method,6,-1)) %>% mutate(`rho factor` = factor(rho_factor),
                                                    `noise factor` = noise_factor,
                                                    `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
                                                    y = Khat_accuracy, 
                                                    color = method,
                                                    shape = `rho factor`) ) +
  geom_point(size=4, alpha = 0.75, position = position_jitter(width = 0.25, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
                                    labeller = labeller(.rows = label_both, 
                                                        .cols = label_both)) + 
  ylab("Khat accuracy") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) 
# position = position_jitter(width = 0.1, height = 0.1)
# Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
#   pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
#   mutate(method = stri_sub(method,6,-1)) %>%
#   ggplot(aes(x = Khat)) + geom_histogram() + facet_grid(setting_id~method)

# ARI_kcentroid_nuclear_MFM
# ARI_kcentroid_spectral_MFM
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% 
  group_by(setting_id, method) %>% 
  summarise(mean_ARI = mean(ARI)) %>% 
  mutate(method = stri_sub(method,5,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>% 
  mutate_at(vars(signal_strength, noise_factor, cluster_num, total_sample_size, setting_id), funs(as.character)) %>%
  xtable(digits = 3) %>% print(include.rownames=FALSE)

## 
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% 
  group_by(setting_id, method) %>% summarise(mean_ARI = mean(ARI, na.rm = T)) %>% 
  mutate(method = stri_sub(method,5,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>% 
  unite("sample size, rho factor", total_sample_size:rho_factor, remove = FALSE,
        sep = ",") %>%
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>% 
  mutate(method = factor(method, levels = c("DP", "MFM", "kmeans", "specc", "kcentroid_spectral", "kcentroid_nuclear"),
                         labels = c("DP", "MFM", "K-means", "Spectral", "K-centroid spectral", "K-centroid nuclear"))) %>%
  ggplot(aes(x = `sample size, rho factor`, 
             y = mean_ARI, 
             shape = method,
             color = method) ) +
  geom_point(size=4, alpha = 0.5, position = position_jitter(width = 0.25, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Mean adjusted Rand Index") + xlab("Sample size, rho factor") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14)) 
####################
########## rho = 0.3
####################
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% 
  group_by(setting_id, method) %>% summarise(mean_ARI = mean(ARI, na.rm = T)) %>% 
  mutate(method = stri_sub(method,5,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>% filter(rho_factor == 0.3) %>% 
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
             y = mean_ARI, 
             shape = method,
             color = method) ) +
  geom_point(size=3.5, alpha = 0.5, position = position_jitter(width = 0.2, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Mean adjusted Rand Index") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom") + ggtitle("Rho factor = 0.3")

# rho = 0.9
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% 
  group_by(setting_id, method) %>% summarise(mean_ARI = mean(ARI, na.rm = T)) %>% 
  mutate(method = stri_sub(method,5,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>% filter(rho_factor == 0.9) %>% 
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
             y = mean_ARI, 
             shape = method,
             color = method) ) +
  geom_point(size=2.5, alpha = 0.5, position = position_jitter(width = 0.1, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Mean adjusted Rand Index") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom") + ggtitle("Rho factor = 0.9")

##################
# RMSE VU
RMSE_VU_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, contains("RMSE_VU")) %>%
  pivot_longer(contains("RMSE_VU"), names_to = "method", values_to = "RMSE_VU") %>% 
  group_by(setting_id, method) %>% summarise(mean_RMSE_VU = mean(RMSE_VU),
                                             median_RMSE_VU = median(RMSE_VU)) %>% 
  mutate(method = stri_sub(method,9,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id"))

# focus on a particular scenario
# rho = 0.3, cluster = 6
RMSE_VU_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, contains("RMSE_VU")) %>%
  pivot_longer(contains("RMSE_VU"), names_to = "method", values_to = "RMSE") %>% 
  mutate(method = stri_sub(method,9,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  filter(rho_factor == 0.9 & cluster_num == 6) %>% ggplot(aes(x = RMSE, color = method)) + 
  geom_density() + theme(legend.position = "bottom",
                         axis.text = element_text(size=12),
                         axis.title = element_text(size=14,face="bold"),
                         legend.text = element_text(size = 12)) + facet_wrap(.~total_sample_size,
                                                                             nrow = 1)


#####################
# computational time
#####################
time_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, contains("time")) %>%
  pivot_longer(contains("time"), names_to = "method", values_to = "time") %>% 
  group_by(setting_id, method) %>% summarise(mean_time = mean(time)) %>% 
  mutate(method = stri_sub(method,6,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
             y = mean_time, 
             color = method,
             shape = `rho factor`) ) +
  geom_point(size=3.5, alpha = 0.75, position = position_jitter(width = 0.1, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Mean computation time (mins)") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14)) 




# ## number of clusters
# do.call("c", lapply(sim_10by6_6, function(x) x$Khat_MFM)) %>% table()
# do.call("c", lapply(sim_10by6_5, function(x) x$Khat_DP)) %>% table()
# do.call("c", lapply(sim_10by6_8, function(x) x$Khat_MFM)) %>% table()
# do.call("c", lapply(sim_10by6_8, function(x) x$Khat_DP)) %>% table()
# #########################################
# #########################################
# ##
sim_10by6_28 <- readRDS("./simulation_results/sim_10by6_28.rds")
which(do.call("c", lapply(sim_10by6_28, function(x) "error" %in% class(x))) == TRUE)
####
do.call("c", lapply(sim_10by6_28, function(x) x$Khat_DP)) %>% table()
do.call("c", lapply(sim_10by6_28, function(x) x$Khat_MFM)) %>% table()
# # ##
# do.call("c", lapply(sim_10by6_1, function(x) x$Khat_MFM)) %>% table()
# do.call("c", lapply(sim_10by6_1, function(x) x$Khat_DP)) %>% table()
# # ##
# # 
# sim_10by6_8 <- readRDS("./simulation_results/sim_10by6_8.rds")
# ##
# do.call("c", lapply(sim_10by6_8, function(x) x$Khat_MFM)) %>% table()
# do.call("c", lapply(sim_10by6_8, function(x) x$Khat_DP)) %>% table()
# # 
# # 
# # ##
# # sim_10by6_4 <- readRDS("./simulation_results/sim_10by6_4.rds")
# # ##
# # do.call("c", lapply(sim_10by6_4, function(x) x$Khat_MFM)) %>% table()
# # do.call("c", lapply(sim_10by6_4, function(x) x$Khat_DP)) %>% table()

####
# 25 by 18
# sim_25by18_1 <- readRDS("./simulation_results/sim_25by18_1_new.rds")
# #####
# do.call("c", lapply(sim_25by18_1, function(x) max(x$MFM_Dahl_rlt$zout) )) %>% table()
simulation_settings_25by18 <- data.frame(expand.grid(signal_strength = 1, 
                                                     cluster_num = c(3),
                                                     noise_factor = c(1.5, 1, 0.5), 
                                                     total_sample_size = c(200),
                                                     rho_factor = c(0.9, 0.3))) %>% 
  rowid_to_column("setting_id")
#####
Khat_rlt_df_25by18 <- data.frame()
################
# computational time
time_rlt_df_25by18 <- data.frame()
RMSE_VU_rlt_df_25by18 <- data.frame()
## 10 by 6
for(i in c(1:nrow(simulation_settings_25by18))){
  #
  assign(paste0("sim_25by18_", i), readRDS( paste0("./simulation_results/sim_25by18_", i, "_new.rds")) )
  temp_rlt <- get(paste0("sim_25by18_", i))
  #
  temp_true_K <- simulation_settings_10by6$cluster_num[i]
  #
  # temp_rlt_df <- data.frame(Khat_DP = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
  #                           Khat_MFM = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
  #                           Khat_kmeans = do.call("c", lapply(temp_rlt, function(x) x$Khat_kmeans_MFM)),
  #                           Khat_specc = do.call("c", lapply(temp_rlt, function(x) x$Khat_specc_MFM)) )
  temp_rlt_df <- data.frame(ARI_DP = do.call("c", lapply(temp_rlt, function(x) x$ARI_DP)),
                            ARI_MFM = do.call("c", lapply(temp_rlt, function(x) x$ARI_MFM)),
                            ARI_kmeans = do.call("c", lapply(temp_rlt, function(x) x$ARI_kmeans_oracle)),
                            ARI_specc = do.call("c", lapply(temp_rlt, function(x) x$ARI_specc_oracle)),
                            ARI_kcentroid_nuclear = do.call("c", lapply(temp_rlt, function(x) x$ARI_kcentroid_nuclear_oracle)),
                            ARI_kcentroid_spectral = do.call("c", lapply(temp_rlt, function(x) x$ARI_kcentroid_spectral_oracle)),
                            Khat_DP = do.call("c", lapply(temp_rlt, function(x) x$Khat_DP)),
                            Khat_MFM = do.call("c", lapply(temp_rlt, function(x) x$Khat_MFM)),
                            K_truth = temp_true_K,
                            setting_id = i) %>% rowid_to_column(var = "replicate_id")
  
  #
  Khat_rlt_df_25by18 <- rbind(Khat_rlt_df_25by18,
                             temp_rlt_df)
  
  # 
  time_rlt_df_25by18 <- rbind(time_rlt_df_25by18,
                             data.frame(time_DP = do.call("c", lapply(temp_rlt, function(x) x$MFM_DP_time)),
                                        time_MFM = do.call("c", lapply(temp_rlt, function(x) x$MFM_Dahl_time)),
                                        setting_id = i) %>% rowid_to_column(var = "replicate_id") ) 
  
  # 
  RMSE_VU_rlt_df_25by18 <- rbind(RMSE_VU_rlt_df_25by18,
                                data.frame(RMSE_VU_DP = do.call("c", lapply(temp_rlt, function(x) x$VU_RMSE_DP)),
                                           RMSE_VU_MFM = do.call("c", lapply(temp_rlt, function(x) x$VU_RMSE_MFM)),
                                           setting_id = i) %>% rowid_to_column(var = "replicate_id")) 
  
}
########
#########################
Khat_rlt_df_25by18 %>% group_by(setting_id) %>% summarise(n = n()) %>% xtable()
#########################
#########################
#########################
#########################
Khat_rlt_df_25by18 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
  pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
  group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
  inner_join(simulation_settings_25by18, by = c("setting_id" = "setting_id")) %>%
  mutate(method = stri_sub(method,6,-1)) %>% mutate(`rho factor` = factor(rho_factor),
                                                    `noise factor` = noise_factor,
                                                    `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
             y = Khat_accuracy, 
             color = method,
             shape = `rho factor`) ) +
  geom_point(size=4, alpha = 0.75, position = position_jitter(width = 0.25, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Khat accuracy") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14)) 

#########################
RMSE_VU_rlt_df_25by18 %>% dplyr::select(replicate_id, setting_id, contains("RMSE_VU")) %>%
  pivot_longer(contains("RMSE_VU"), names_to = "method", values_to = "RMSE") %>% 
  mutate(method = stri_sub(method,9,-1)) %>%
  inner_join(simulation_settings_25by18, by = c("setting_id" = "setting_id")) %>%
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>% ggplot(aes(x = RMSE, color = method)) + 
  facet_grid(`noise factor`~`rho factor`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) +
  geom_density() + theme(legend.position = "bottom",
                         axis.text = element_text(size=12),
                         axis.title = element_text(size=14,face="bold"),
                         legend.text = element_text(size = 12)) 
#########################
#########################
Khat_rlt_df_25by18 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% 
  group_by(setting_id, method) %>% summarise(mean_ARI = mean(ARI, na.rm = T)) %>% 
  mutate(method = stri_sub(method,5,-1)) %>%
  inner_join(simulation_settings_25by18, by = c("setting_id" = "setting_id")) %>% 
  unite("sample size, rho factor", total_sample_size:rho_factor, remove = FALSE,
        sep = ",") %>%
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>% 
  mutate(method = factor(method, levels = c("DP", "MFM", "kmeans", "specc", "kcentroid_spectral", "kcentroid_nuclear"),
                         labels = c("DP", "MFM", "K-means", "Spectral", "K-centroid spectral", "K-centroid nuclear"))) %>%
  ggplot(aes(x = `sample size, rho factor`, 
             y = mean_ARI, 
             shape = method,
             color = method) ) +
  geom_point(size=4, alpha = 0.5, position = position_jitter(width = 0.25, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Mean adjusted Rand Index") + xlab("Sample size, rho factor") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14)) 
#########################
#########################
time_rlt_df_25by18 %>% dplyr::select(replicate_id, setting_id, contains("time")) %>%
  pivot_longer(contains("time"), names_to = "method", values_to = "time") %>% 
  group_by(setting_id, method) %>% summarise(mean_time = mean(time),
                                             median_time = median(time)) %>% 
  mutate(method = stri_sub(method,6,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  mutate(`rho factor` = factor(rho_factor),
         `noise factor` = noise_factor,
         `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
             y = mean_time, 
             color = method,
             shape = `rho factor`) ) +
  geom_point(size=3.5, alpha = 0.75, position = position_jitter(width = 0.1, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Mean computation time (mins)") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14)) 
#########################
#########################
# computational time
sim_10by6_29_GAMMA_2 <- readRDS("./simulation_results/sim_10by6_29_GAMMA_2.rds")
#
do.call("c", lapply(sim_10by6_29, function(x) x$Khat_DP)) %>% table()
do.call("c", lapply(sim_10by6_29_GAMMA_2, function(x) x$Khat_DP)) %>% table()
#
do.call("c", lapply(sim_10by6_29, function(x) x$Khat_MFM)) %>% table()
do.call("c", lapply(sim_10by6_29_GAMMA_2, function(x) x$Khat_MFM)) %>% table()
#
sim_10by6_35_GAMMA_2 <- readRDS("./simulation_results/sim_10by6_35_GAMMA_2.rds")
do.call("c", lapply(sim_10by6_35_GAMMA_2, function(x) x$Khat_MFM)) %>% table()

##########################
##########################
# 
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("ARI")) %>%
  pivot_longer(contains("ARI"), names_to = "method", values_to = "ARI") %>% 
  group_by(setting_id, method) %>% summarise(mean_ARI = mean(ARI, na.rm = T)) %>% 
  mutate(method = stri_sub(method,5,-1)) %>%
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id"))

#
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
  pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
  group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  mutate(method = stri_sub(method,6,-1)) %>% mutate(`rho factor` = factor(rho_factor),
                                                    `noise factor` = noise_factor,
                                                    `cluster number` = cluster_num) %>%
  ggplot(aes(x = factor(total_sample_size), 
             y = Khat_accuracy, 
             color = method,
             shape = `rho factor`) ) +
  geom_point(size=4, alpha = 0.75, position = position_jitter(width = 0.25, height=0)) + 
  facet_grid(`noise factor`~`cluster number`,
             labeller = labeller(.rows = label_both, 
                                 .cols = label_both)) + 
  ylab("Khat accuracy") + xlab("Sample size") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) 
###
Khat_rlt_df_10by6 %>% dplyr::select(replicate_id, setting_id, K_truth, contains("Khat")) %>%
  pivot_longer(contains("Khat"), names_to = "method", values_to = "Khat") %>% 
  group_by(setting_id, method) %>% summarise(Khat_accuracy = mean(Khat == K_truth)) %>% 
  inner_join(simulation_settings_10by6, by = c("setting_id" = "setting_id")) %>%
  mutate(method = stri_sub(method,6,-1)) %>% filter(setting_id == 29)