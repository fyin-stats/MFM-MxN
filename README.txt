# Readme file for the MFM-MxVN codes
# Please do not forget to set up the path as the local path of this folder before running the following scripts

# helper function
MFM_MN.R contains main functions for MFM-MxVN algorithm as well as the calculation of partition distribution

# simulation study - p=10, q=6
MFM_simulation_study_10_6_signal_strength_1_3clusters.R 
MFM_simulation_study_10_6_signal_strength_1_noise025_3clusters.R

# simulation study - p=25, q=18
MFM_simulation_study_25_18_new_kronecker_high_noise.R

# summary of simulation studies (need the output of MFM_simulation_study_10_6_signal_strength_1_3clusters.R, MFM_simulation_study_10_6_signal_strength_1_noise025_3clusters.R and MFM_simulation_study_25_18_new_kronecker_high_noise.R)
MFM_simulation_summary.R

# NBA data analysis
MFM_NBA_log_mean_intensity_25_18.R

# summary of NBA data analysis (need the output of MFM_NBA_log_mean_intensity_25_18.R)
MFM_NBA_post_summary.R



