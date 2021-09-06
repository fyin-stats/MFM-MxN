##################
##################
player_names_caps <- readRDS("./player_names_caps.rds")
load("./shotdata.Rdata")
load("./lgcp_fit.Rdata")
source("./MFM_MN.R")
#######################
dim(intensity.coef)
summary(coords) # y: 0, 50 (sideline to sideline); x: 0 - 46 (mid court line to baseline)
NameList
NBA_log_intensity_post_results <- readRDS("./NBA_log_intensity_post_result.rds")
########################
NBA_log_intensity_post_results$Mout[,,1]
grid_df <- expand.grid(x=1:nrow(NBA_log_intensity_post_results$Mout[,,1]), 
                       y=1:ncol(NBA_log_intensity_post_results$Mout[,,1]))
#
NBA_log_intensity_Mout_df <- data.frame()
#
NBA_log_intensity_zout_list <- vector(mode="list",
                                      length = max(NBA_log_intensity_post_results$zout))
#
for(i in 1:max(NBA_log_intensity_post_results$zout)){
  # 
  NBA_log_intensity_Mout_df <- rbind(NBA_log_intensity_Mout_df,
                                     data.frame(grid_df, 
                                                M = 2 + c(NBA_log_intensity_post_results$Mout[,,i]),
                                                cluster = paste0("Cluster", i)) )
  # 
  NBA_log_intensity_zout_list[[i]] <- NameList[which(NBA_log_intensity_post_results$zout==i)]
}
# plot the resulting intensity matrices
NBA_log_intensity_Mout_df %>% ggplot(aes(x = x, y = y, fill = M)) + 
  geom_tile() + facet_wrap(.~cluster) + scale_fill_distiller(palette='RdBu', name = "log(intensity)") +
  theme(legend.text=element_text(size=15),
        legend.position = "bottom",
        strip.text = element_text(size = 12))
# name list
NBA_log_intensity_post_results$zout
# cluster 1
library(xtable)
matrix(NBA_log_intensity_zout_list[[1]], ncol=5) %>% xtable()
# cluster 2
matrix(NBA_log_intensity_zout_list[[2]], ncol=5) %>% xtable()
# cluster 3
matrix(NBA_log_intensity_zout_list[[3]], ncol=5) %>% xtable()
# cluster 4
matrix(NBA_log_intensity_zout_list[[4]], ncol=5) %>% xtable()
# cluster 5
matrix(NBA_log_intensity_zout_list[[5]], ncol=5) %>% xtable()
###########################
###########################
## covariance matrix visualization
###########################
NBA_log_intensity_25_18_MLE_initial_6000iters[[replicate_id]]$Uout %>% dim()
NBA_log_intensity_25_18_MLE_initial_6000iters[[replicate_id]]$Vout %>% dim()
###########################
grid_U <- expand.grid(x=1:dim(NBA_log_intensity_post_results$Uout)[1], 
                      y=1:dim(NBA_log_intensity_post_results$Uout)[2])
###########################
U_df <- data.frame(grid_U, value = c(NBA_log_intensity_post_results$Uout))
###########################
U_df %>% ggplot(aes(x = x, y = y, fill = value)) + 
  geom_tile() + theme(legend.position = "right",
                      
                      axis.text=element_text(size=12),
                      axis.title=element_text(size=14,face="bold"),
                      strip.text = element_text(size = 12)) +
  scale_fill_distiller(palette='RdBu', name = "")
############################
############################
grid_V <- expand.grid(x=1:dim(NBA_log_intensity_post_results$Vout)[1], 
                      y=1:dim(NBA_log_intensity_post_results$Vout)[2])
###########################
V_df <- data.frame(grid_V, value = c(NBA_log_intensity_post_results$Vout))
###########################
V_df %>% ggplot(aes(x = x, y = y, fill = value)) + 
  geom_tile() + theme(legend.position = "right",
                      legend.text=element_text(size=14),
                      axis.text=element_text(size=12),
                      axis.title=element_text(size=14,face="bold"),
                      strip.text = element_text(size = 12)) +
  scale_fill_distiller(palette='RdBu', name = "")
###########################
###########################
###########################