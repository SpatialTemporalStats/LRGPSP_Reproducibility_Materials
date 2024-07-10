# Reproducibility Instrucion for "Large-Scale Low-Rank Gaussian Process Prediction with Support Points"
This file documents the artifacts associated with the article (i.e., the data and code supporting the computational findings) and describes how to reproduce all figures and results. 

## Article Overview
Low-rank approximation is a popular strategy to tackle the “big $n$ problem” associated with large-scale Gaussian process regressions. Basis functions for developing low-rank structures are crucial and should be carefully specified. Predictive processes simplify the problem by inducing basis functions with a covariance function and a set of knots or representative points (rep-points). The existing literature suggests certain practical implementations of rep-points selection and covariance estimation; however, theoretical foundations explaining the influence of these two factors on predictive processes are lacking. In this paper, the asymptotic prediction performance of the predictive process and Gaussian process predictions are derived and the impacts of the selected rep-points and estimated covariance are studied. The use of support points (SPs) as knots, which best represent data locations, is advocated. Extensive simulation studies demonstrate the superiority of support points and verify our theoretical results. Real data of precipitation and ozone are used as examples, and the efficiency of our method over other widely used low-rank approximation methods is verified.

## Contents
#### Rep_Points
This sub-repository provides R scripts for selecting various kinds of rep-points, including those based on the fixed-size determinantal point process (k-DPPs), grid points (Grids), random samples (Rands), and support points (SPs).

- "k_DPP"
  - "kdpp-sample.R": functions for selecting k-DPPs
  - "helper-functions.R": helper functions for "kdpp-sample.R"
- "Grid-sample.R": functions for selecting Grids
- "Rand-sample.R": functions for selecting SPs
- "Rep_Points_README.md": more details about this sub-repository

#### SP_Demo
This sub-repository provides the R script for reproducing Figures 1 and S2, which illustrate the patterns of various kinds of rep-points under two location sets.

- "SP_Demo.R": R script for reproducing Figures 1 and S2

#### S1_Gamma
This sub-repository provides R scripts for reproducing Figure 2 and Table 1, which demonstrate values of smoothness parameter $\gamma$ for the Matern covariance function under various parameter settings. 

- "Gamma_values_single.R": R script for reproducing Figures 1 and S2 using results from a single replicate (less computational time)
- "Gamma_values_single.R": R script for reproducing Figures 1 and S2 using results from multiple replicates 
- "helper-functions.R": a helper function

#### S2_VersusK
This sub-repository includes R scripts and outputs for reproducing Figure 3 and Table 2. Figure 3 demonstrates the impact of rep-points on the predictive process by comparing the prediction performance of various types of representative points across four scenarios. Table 2 displays the variations in energy distances for rep-points with respect to their sizes.

- "Simu2.R": R script for reproducing Figure 3 and Table 2
- "helper-functions.R": a helper function
- "Setting1.csv", "Setting2.csv", "Setting3.csv", and "Setting4.csv": outputs from the simulated studies
- "Set1Sum.csv", "Set2Sum.csv", "Set3Sum.csv", and "Set4Sum.csv": summarized results for plotting Figure 3

#### S3_VersusN
This sub-repository includes R scripts and outputs for reproducing Figures 4, 5, and S1, which are used to verify Theorems 1 and 2. Figures 4 and S1(a) demonstrate the convergence rates of the predictive process predictions with a sufficient number of rep-points and an estimated covariance function. Figures 5 and S1(b) compare smoothness parameters $\gamma$ derived by different approaches.

- "Simu3.R": R script for reproducing Figures 4, 5, and S1
- "Setting1.csv" and "SettingS.csv": outputs from the simulated studies
- "Set1_GS.csv", "Set1_MSPE.csv", "SetS_GS.csv", and "SetS_MSPE.csv": summarized results for plotting Figures 4, 5, and S1

#### S4_VersusTau
This sub-repository includes R scripts and outputs for reproducing Figure 6, which demonstrates how the nugget effect influences predictive process predictions and is used to verify Theorem 3.

- "Simu4.R": R script for reproducing Figure 6
- "Tau.csv": outputs from the simulated studies

## Reproducibility Workflow
Please begin by downloading the entire repository as the file "LRGPSP_Reproducibility_Materials.zip" and extracting it as a folder named "LRGPSP_Reproducibility_Materials". Next, set your working directory to this folder. Then, load necessary R packages and functions. Finally, follow the "Wrapper.R" file to reproduce each figure and table sequentially. The computation time reported below was recorded using (R 3.6.3) running on machine equipped with Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz and 125 GB RAM.

#### Reproduce Figure 1 in Section 2.1 (and Figure S2 in Section S5)
Figures 1 and S2 demonstrate the patterns of various kinds of representative points under two location sets. Please refer to the "Wrapper.R" file for their reproducibility command. The total computational time is approximately 39.1 minutes, with the majority of the time spent obtaining k-DPPs.  For more detailed code and computational time, please refer to the sub-repository "SP_Demo". 

#### Reproduce Figure 2 and Table 1 in Section 3.1
Figure 2 and Table 1 illustrate the values of smoothness parameter $\gamma$ for the Matern covariance function under various parameter settings. Please refer to the "Wrapper.R" file for their reproducibility command. The total computational time is approximately 2 hours. For more detailed code, outputs, and computational time, please refer to the sub-repository "S1_Gamma". 

#### Reproduce Figure 3 and Table 2 in Section 3.2
Figure 3 and Table 2 demonstrate the influence of rep-points on the predictive process. Please refer to the "Wrapper.R" file for the commands to reproduce these results. There are three options for reproducibility. The first option, set as the default, reproduces Figure 3 and Table 2 using outputs from a single replicate, which takes approximately 16.5 hours. The second option, used in the article, reproduces the results using outputs from 100 replicates but requires significantly more time. The third option uses our pre-generated outputs directly. To select the second or third options, please manually adjust the settings in "Simu2.R". For more detailed code, outputs, and computational time, please refer to the sub-repository "S2_VersusK".

#### Reproduce Figures 4 and 5 in Section 3.3 (and Figure S1 in Section S4)
Figures 4, 5, and S1 are used to verify Theorems 1 and 2. Figures 4 and S1(a) show how the $\log$(MSPE)s of the predictive process predictions change as the sample size $n$ grows. Figure 5 and S1(b) display $\gamma$ obtained by two approaches. Please refer to the "Wrapper.R" file for the commands to reproduce these results. The outputs are based on 100 replicates. Generating outputs for Figures 4 and 5 from a single replicate takes approximately 6.5 minutes. For Figure S1, it takes about 11.9 minutes per replicate. Using 4 cores as specified in "Simu3.R", the entire process will take around 7.7 hours. Readers can adjust the number of cores used or opt to use our pre-generated outputs directly. For more detailed code, outputs, and computational time, please refer to the sub-repository "S3_VersusN".

## Data 

## Supporting Software Requirement
#### Version of primary software used
R version 3.6.3 (2020-02-29)

#### Libraries and dependencies used by the code

assertthat_0.2.1, autoFRK_1.4.3, doParallel_1.0.17, fields_13.3, foreach_1.5.2, FRK_2.0.5, geoR_1.9-2, ggplot2_3.3.6, GpGp_0.4.0, here_1.0.1, LatticeKrig_8.4, magrittr_2.0.3, mvtnorm_1.1-3, purrr_0.3.4, rlist_0.4.6.2, sp_1.4-6, Split_1.2, support_0.1.4





