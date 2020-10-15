##############
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack",
              "mvtnorm","dplyr","armspp","ggcorrplot","reshape2","ggpubr",
              "cowplot","doParallel","foreach","MixMatrix","mclust","mcclust","fossil","plotly")
# "rjson","ggplot2","grid","inlabru","jpeg",
#  "RCurl","sp","INLA","Orcs","dplyr","gridExtra","ggsci"
# https://cran.r-project.org/web/packages/MixMatrix/vignettes/matrixnormal.html
ipak(packages)
###############
# post summary
# 25*28
setwd("/Users/fan/Documents/MFM_matrix/")
################
player_names_caps <- readRDS("./player_names_caps.rds")
load("./shotdata.Rdata")
load("./lgcp_fit.Rdata")
NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters <- readRDS("NBA_log_mean_intensity_25_18_MLE_initial_6000iters_random_initclusters.rds")
NBA_rlt_Dahl_log_intensity <- NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[4]]
do.call("c",lapply(NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters, function(x) max(x$zout))) %>% table()
which(do.call("c",lapply(NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters, function(x) max(x$zout))))
#
Khat <- max(NBA_rlt_Dahl_log_intensity$zout)
# NBA_M_10_10_all
NBA_M_25_18_all <- data.frame(Var1=NA, Var2=NA, value=NA, cluster=NA)
for(k in 1:Khat){
  temp_M <- NBA_rlt_Dahl_log_intensity$Mout[,,k]
  rownames(temp_M) <- 1:25
  colnames(temp_M) <- 1:18
  temp_M_long <- melt(temp_M)
  temp_M_long$cluster <- k
  NBA_M_25_18_all <- rbind(NBA_M_25_18_all, temp_M_long)
}
# remove the NA
NBA_M_25_18_all <- NBA_M_25_18_all[complete.cases(NBA_M_25_18_all),]
# visualize
ggplot(NBA_M_25_18_all, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill=exp(value))) +
  scale_fill_gradient(low="white",high="blue") +
  labs(x="X", y="Y", title=paste("M")) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) + facet_wrap(~cluster,
                                                                    nrow=2)
###############################
###############################
#####
mean_RI_25_18_vec <- rep(NA, length(NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters ))
for(i in 1:length(NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters)){
  mean_RI_25_18_vec[i] <- mean(do.call("c", lapply(NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters,function(x) fossil::adj.rand.index(NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[i]]$zout, x$zout)))[-i])
}
#####
# which(mean_RI_25_18_vec==max(mean_RI_25_18_vec))
# #####
# NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[45]]$zout %>% max()
# ######
# # 1 10 14 16 17 20 22 25 27 36 42 44 45
# NameList[NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[17]]$zout==2]

###### levelplot for covariance matrix estimates
# U
tool_matrix <- matrix(1, nrow = nrow(NBA_rlt_Dahl_log_intensity$Uout),ncol = ncol(NBA_rlt_Dahl_log_intensity$Uout))
#
X <- data.frame(cbind(which(tool_matrix==1, arr.ind = T),1:length(c(NBA_rlt_Dahl_log_intensity$Uout))))
colnames(X) <- c("x", "y", "linear_index")
Ufit_dat <- data.frame(y = X$y, x = X$x, value = c(NBA_rlt_Dahl_log_intensity$Uout))
#############
# ggsave(p, file = "../manuscript/plots/IntenFit.pdf", height = 5, width = 4)
# levelplot(value~x*y,Ufit_dat, 
#           col.regions = terrain.colors(100)[length(terrain.colors(100)):1], 
#           xlab = "X", ylab = "Y",
#           scales = list(x = list(at = seq(0, nrow(NBA_rlt_Dahl_log_intensity$Uout), length.out = 5)), 
#                         y = list(at = seq(0, ncol(NBA_rlt_Dahl_log_intensity$Uout), length.out = 5))))
#############
breaks <- seq(min(c(Ufit_dat$value)), max(c(Ufit_dat$value)), length.out = 100)
cols <- terrain.colors(100)[length(terrain.colors(100)):1]
# colorRampPalette(c("grey", "yellow", "green"))(length(breaks) - 1)
######
pU <- levelplot(value~x*y,Ufit_dat, 
                  col.regions = cols, at = breaks, 
                  xlab = "", ylab = "", main = "U",
                  scales = list(x = list(at = seq(0, nrow(NBA_rlt_Dahl_log_intensity$Uout), length.out = 6)), 
                                y = list(at = seq(0, ncol(NBA_rlt_Dahl_log_intensity$Uout), length.out = 6))))


###############
### V 
tool_matrix <- matrix(1, nrow = nrow(NBA_rlt_Dahl_log_intensity$Vout),ncol = ncol(NBA_rlt_Dahl_log_intensity$Vout))
#
X <- data.frame(cbind(which(tool_matrix==1, arr.ind = T),1:length(c(NBA_rlt_Dahl_log_intensity$Vout))))
colnames(X) <- c("x", "y", "linear_index")
Vfit_dat <- data.frame(y = X$y, x = X$x, value = c(NBA_rlt_Dahl_log_intensity$Vout))
#############
# ggsave(p, file = "../manuscript/plots/IntenFit.pdf", height = 5, width = 4)
# levelplot(value~x*y,Vfit_dat, 
#           col.regions = terrain.colors(100)[length(terrain.colors(100)):1], 
#           xlab = "X", ylab = "Y",
#           scales = list(x = list(at = seq(0, nrow(NBA_rlt_Dahl_log_intensity$Vout), length.out = 5)), 
#                         y = list(at = seq(0, ncol(NBA_rlt_Dahl_log_intensity$Vout), length.out = 5))))
#############
breaks <- seq(min(c(Vfit_dat$value)), max(c(Vfit_dat$value)), length.out = 100)
cols <- terrain.colors(100)[length(terrain.colors(100)):1]
# colorRampPalette(c("grey", "yellow", "green"))(length(breaks) - 1)
######
pV <- levelplot(value~x*y,Vfit_dat, 
                col.regions = cols, at = breaks, 
                xlab = "", ylab = "", main = "V",
                scales = list(x = list(at = seq(0, nrow(NBA_rlt_Dahl_log_intensity$Vout), length.out = 4)), 
                              y = list(at = seq(0, ncol(NBA_rlt_Dahl_log_intensity$Vout), length.out = 4))))

#### VxU
VxU <- NBA_rlt_Dahl_log_intensity$Vout %x% NBA_rlt_Dahl_log_intensity$Uout
####
tool_matrix <- matrix(1, nrow = nrow(VxU),ncol = ncol(VxU))
#
X <- data.frame(cbind(which(tool_matrix==1, arr.ind = T),1:length(c(VxU))))
colnames(X) <- c("x", "y", "linear_index")
VxUfit_dat <- data.frame(y = X$y, x = X$x, value = c(VxU))
#############
# ggsave(p, file = "../manuscript/plots/IntenFit.pdf", height = 5, width = 4)
# levelplot(value~x*y,Vfit_dat, 
#           col.regions = terrain.colors(100)[length(terrain.colors(100)):1], 
#           xlab = "X", ylab = "Y",
#           scales = list(x = list(at = seq(0, nrow(NBA_rlt_Dahl_log_intensity$Vout), length.out = 5)), 
#                         y = list(at = seq(0, ncol(NBA_rlt_Dahl_log_intensity$Vout), length.out = 5))))
#############
breaks <- seq(min(c(VxUfit_dat$value)), max(c(VxUfit_dat$value)), length.out = 100)
cols <- terrain.colors(100)[length(terrain.colors(100)):1]
# colorRampPalette(c("grey", "yellow", "green"))(length(breaks) - 1)
######
pVxU <- levelplot(value~x*y,VxUfit_dat, 
                col.regions = cols, at = breaks, 
                xlab = "", ylab = "", main = "V kronecker U",
                scales = list(x = list(at = seq(0, nrow(VxU), length.out = 4)), 
                              y = list(at = seq(0, ncol(VxU), length.out = 4))))

# #################################
# fossil::rand.index(NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[1]]$zout,
#                    NBA_rlt_Dahl_log_intensity_25_18_50reps[[3]]$zout)
# #####
# # player names
# NameList[NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[1]]$zout==1] %>% matrix(ncol=10) %>% xtable()
# 
# NameList[NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[1]]$zout==2]
# 
# NameList[NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[1]]$zout==3] 
# 
# NameList[NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters[[1]]$zout==3] %>% matrix(ncol=5) %>% xtable()
# 
# 
# 
# ######
# NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters_nclusters <- data.frame(count = c(do.call("c",lapply(NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters, function(x) max(x$zout))) %>% table()),
#                                                                       n_clusters = 2:6)
# 
# NBA_rlt_Dahl_log_mean_intensity_25_18_50reps_6000iters_nclusters %>% ggplot(aes(x=n_clusters,y=count)) + geom_bar(stat="identity") + xlab("Number of clusters") + ylab("Counts") + theme(text = element_text(size=20),
#                                                                                                                                                                                        axis.text.x = element_text(angle = 45, hjust=1),
#                                                                                                                                                                                        strip.text = element_text(size = 20)) 
#                                                                       