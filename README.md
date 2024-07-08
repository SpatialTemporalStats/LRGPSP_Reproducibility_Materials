# Repeoducibility Instrucion for "Large-Scale Low-Rank Gaussian Process Prediction with Support Points"
This file documents the artifacts associated with the article (i.e., the data and code supporting the computational findings) and describes how to reproduce all figures and results. Please begin by downloading the entire repository as the file "LRGPSP_Reproducibility_Materials.zip" and extracting it as a folder named "LRGPSP_Reproducibility_Materials". Next, please set your working directory to this folder. The computation time reported below was recorded using (R 3.6.3) running on machine equipped with Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz and 125 GB RAM.

## Artical Overview
Low-rank approximation is a popular strategy to tackle the “big $n$ problem” associated with large-scale Gaussian process regressions. Basis functions for developing low-rank structures are crucial and should be carefully specified. Predictive processes simplify the problem by inducing basis functions with a covariance function and a set of knots or representative points (rep-points). The existing literature suggests certain practical implementations of rep-points selection and covariance estimation; however, theoretical foundations explaining the influence of these two factors on predictive processes are lacking. In this paper, the asymptotic prediction performance of the predictive process and Gaussian process predictions are derived and the impacts of the selected rep-points and estimated covariance are studied. The use of support points (SPs) as knots, which best represent data locations, is advocated. Extensive simulation studies demonstrate the superiority of support points and verify our theoretical results. Real data of precipitation and ozone are used as examples, and the efficiency of our method over other widely used low-rank approximation methods is verified.

## Contents
#### Rep_Points
This sub-repository provides R scripts for selecting various kinds of rep-points, including those based on the fixed-size determinantal point process (k-DPPs), grid points (Grids), random samples (Rands), and support points (SPs).

- "Rep_Points_README.md": details about this sub-repository
- "k_DPP"
  - "kdpp-sample.R": functions for selecting k-DPPs
  - "helper-functions.R": helper functions for "kdpp-sample.R"
- "Grid-sample.R": functions for selecting Grids
- "Rand-sample.R": functions for selecting SPs

#### SP_Demo
This sub-repository provides the R script for reproducing Figures 1 and S2, which illustrate the patterns of various kinds of rep-points under two location sets.

- "SP_Demo.R": R script for reproducing Figures 1 and S2

#### S1_Gamma
This sub-repository includes scripts for reproducing Figure 2 and Table 1, which demonstrate values of smoothness parameter $\gamma$ for the Matern covariance function under various parameter settings. 


## Reproducibility Workflow
Please begin by downloading the entire repository as the file "LRGPSP_Reproducibility_Materials.zip" and extracting it as a folder named "LRGPSP_Reproducibility_Materials". Next, set your working directory to this folder. Then, load necessary R packages and functions. Finally, follow the "Wrapper.R" file to reproduce each figure and table sequentially.

#### Figure 1 in Section 2.1 (and Figure S2 in Section S5)
Figures 1 and S2 demonstrate the patterns of various kinds of representative points under two location sets. For more detailed code and computational time, please refer to the sub-repository "SP_Demo". Additionally, refer to the "Wrapper.R" file for the reproducibility command. The total computational time for plotting Figures 1 and S2 is approximately 39.1 minutes, with the majority of the time spent obtaining k-DPPs.

#### Figure 2 and Table 1 in Section 3.1


## Data 

## Supporting Software Requirement
#### Version of primary software used
R version 3.6.3 (2020-02-29)

#### Libraries and dependencies used by the code

assertthat_0.2.1, autoFRK_1.4.3, doParallel_1.0.17, fields_13.3, foreach_1.5.2, FRK_2.0.5, geoR_1.9-2, ggplot2_3.3.6, GpGp_0.4.0, here_1.0.1, LatticeKrig_8.4, magrittr_2.0.3, mvtnorm_1.1-3, purrr_0.3.4, rlist_0.4.6.2, sp_1.4-6, Split_1.2, support_0.1.4





