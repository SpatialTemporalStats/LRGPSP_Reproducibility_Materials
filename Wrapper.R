################  Wrapper file for sequentially implementing each figure and table  ################
# This file outlines the reproducibility workflow of the article. 
# Please begin by downloading the entire repository as the file "LRGPSP_Reproducibility_Materials.zip" and extracting it as a folder named "LRGPSP_Reproducibility_Materials". 
# Next, please set your working directory to this folder. 
# The computation time reported below was recorded using (R 3.6.3) running on machine equipped with Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz and 125 GB RAM.



###### Load necessary R packages and functions
library(here)      # ! Please load this package after setting the working directory to folder "LRGPSP_Reproducibility_Materials"
library(assertthat)
library(autoFRK)
library(FRK)
library(geoR)
library(GpGp)
library(LatticeKrig)
library(fields)
library(mvtnorm)
library(SPlit)
library(support)
library(magrittr)
library(foreach)
library(doParallel)
library(purrr)
library(rlist)
library(sp)
library(ggplot2)
source(here("Rep_Points/k_DPP","helper-functions.R"))
source(here("Rep_Points/k_DPP","kdpp-sample.R"))
source(here("Rep_Points","SP-sample.R"))
source(here("Rep_Points","Rand-sample.R"))
source(here("Rep_Points","Grid-sample.R"))
source(here("S1_Gamma","helper-functions.R"))
source(here("S2_VersusK","helper-functions.R"))
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")



######   Figure 1 in Section 2.1 (and Figure S2 in Section S5)   #######################
# Figures 1 and S2 are used to illustrate the patterns of various kinds of rep-points, #
# including SPs, k-DPPs, random subsamples (Rands), and grid points (Grids).           #
# Code can be found in sub-repository "SP_Demo".                                       #
########################################################################################
# t1=proc.time()[[3]]
source(here("SP_Demo","SP_Demo.R"))
# t2=proc.time()[[3]]
# t2-t1=2344.043



######   Figure 2 and Table 1 in Section 3.1   #########################################
# Figure 2 and Table 1 show the values of smoothness parameter \gamma for Matern       #
# covariance function under various parameter settings.                                #
# Code and outputs can be found in sub-repository "S1_Gamma".                          #
########################################################################################
# Reproduce Figure 2 and Table 1 based on results from a single replicate, which takes about 2 hours.
source(here("S1_Gamma","Gamma_values_single.R"))
# Reproduce Figure 2 and Table 1 based on results from multiple replicates.
# source(here("S1_Gamma","Gamma_values.R"))



######   Figure 3 and Table 2 in Section 3.2   #########################################
# Figure 3 compares the rooted mean squared predictive errors (RMSPE) obtained by using#
# various kinds of rep-points under four scenarios. The details about these scenarios  # 
# are given in Section 3.2.                                                            #                                                                    #
# Table 2 shows how the energy distances of rep-points change as the size of rep-points#
# increases. Additionally, it compares energy distances of various kinds of rep-points.#
# Code and outputs can be found in sub-repository "S2_VersusK".                        #
########################################################################################
# Reproduce Figure 3 and Table 2 based on results from a single replicate. It takes about 16.5 hours.
source(here("S2_VersusK","Simu2.R"))
# Note that we provide three options for reproducing Figure 3 and Table 2. Readers could
# choose any of them by changing the code as shown in "Simu2.R".
# 1. As above, reproduce Figure 3 and Table 2 using outputs from a single replicate.
# 2. As done in the article, reproduce Figure 3 and Table 2 using outputs from 100 replicates. This option requires significantly more time.
# 3. Use our pre-generated outputs, which are provided in the sub-repository "S2_VersusK".



######   Figures 4 and 5 in Section 3.3 (and Figure S1 in Section S4)   ################
# Figures 4 and S1(a) show how log(MSPE)s of predictive process predictions with a     #
# sufficient number of rep-points change as the sample size n grows.                   #                                                                    #
# Figures 5 and S1(b) compare the smoothness parameters obtained by two approaches.    #
# Code and outputs can be found in sub-repository "S3_VersusN".                        #
########################################################################################
source(here("S3_VersusN","Simu3.R"))
# To reproduce these figures, we need outputs from 100 replicates. Obtaining outputs for 
# Figures 4 and 5 from a single replicate takes about 6.5 minutes. Obtaining outputs for
# Figure S1 from a single replicate takes about 11.9 minutes. If we use 4 cores in parallel 
# as given in "Simu3.R", this section will take about 7.7 hours.
# Readers could also uses our pre-generated outputs directly.



######   Figure 6 in Section 3.3        ################################################
# Figure 6 show how log(MSPE)s of predictive process predictions with a                #
# sufficient number of rep-points change versus the scales of nugget effect tau^2.     #                                                                    #
# Code and outputs can be found in sub-repository "S4_VersusTau".                      #
########################################################################################
source(here("S4_VersusTau","Simu4.R"))
# Obtaining outputs from a single replicate takes about 1 minute, which is short. 
# Therefore, we suggest readers reproduce Figure 6 using outputs from 100 replicates.
# This is what we did in the article and may take readers about 25 minutes if they use 4 cores in parallel.
# Additionally, we provide options to uses our pre-generated outputs directly.



######   Tables 3 and 4 in Section 4    ################################################
# Tables 3 and 4 compare the performance of various kinds of low-rank approximation    #
# methods using two real datasets.                                                     #
# Code and datasets can be found in sub-repository "Real_Data".                        #
########################################################################################
### To reproduce the results in Table 3
source(here("Real_Data","DP_LRK.R"))
# Obtaining outputs from a single replicate takes about 13.5 minutes. 

### To reproduce the results in Table 4
source(here("Real_Data","DTCO_LRK.R"))
# Obtaining outputs from a single replicate takes about 2.25 hours. 




