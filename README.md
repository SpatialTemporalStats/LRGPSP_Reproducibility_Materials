# Reproducibility Instrucion for "Large-Scale Low-Rank Gaussian Process Prediction with Support Points"
This file documents the artifacts associated with the article (i.e., the data, code, and outputs supporting the computational findings) and describes how to reproduce all figures and results. 

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
This sub-repository provides R scripts for performing the simulated studies in Section 3.1 and reproducing Figure 2 and Table 1, which demonstrate values of smoothness parameter $\gamma$ for the Matern covariance function under various parameter settings. 

- "Gamma_values_single.R": R script for reproducing Figures 1 and S2 using results from a single replicate (less computational time)
- "Gamma_values.R": R script for performing the simulated studies in Section 3.1 and reproducing Figures 1 and S2 using results from multiple replicates 
- "helper-functions.R": a helper function

#### S2_VersusK
This sub-repository includes R scripts for performing the simulated studies in Section 3.2 and reproducing Figure 3 and Table 2, along with the corresponding outputs. Figure 3 demonstrates the impact of rep-points on the predictive process by comparing the prediction performance of various types of representative points across four scenarios. Table 2 displays the variations in energy distances for rep-points with respect to their sizes.

- "Simu2.R": R script for performing the simulated studies in Section 3.2 and reproducing Figure 3 and Table 2
- "helper-functions.R": a helper function
- "Setting1.csv", "Setting2.csv", "Setting3.csv", and "Setting4.csv": outputs from the simulated studies
- "Set1Sum.csv", "Set2Sum.csv", "Set3Sum.csv", and "Set4Sum.csv": summarized results for plotting Figure 3

#### S3_VersusN
This sub-repository includes R scripts for performing the simulated studies in Section 3.3 and reproducing Figures 4, 5, and S1, which are used to verify Theorems 1 and 2, along with the corresponding outputs. Figures 4 and S1(a) demonstrate the convergence rates of the predictive process predictions with a sufficient number of rep-points and an estimated covariance function. Figures 5 and S1(b) compare smoothness parameters $\gamma$ derived by different approaches.

- "Simu3.R": R script for performing the simulated studies in Section 3.3 and reproducing Figures 4, 5, and S1
- "Setting1.csv" and "SettingS.csv": outputs from the simulated studies
- "Set1_GS.csv", "Set1_MSPE.csv", "SetS_GS.csv", and "SetS_MSPE.csv": summarized results for plotting Figures 4, 5, and S1

#### S4_VersusTau
This sub-repository includes R scripts for performing the simulated studies in Section 3.3 and reproducing Figure 6, which demonstrates how the nugget effect influences predictive process predictions and is used to verify Theorem 3.

- "Simu4.R": R script for performing the simulated studies in Section 3.3 and reproducing Figure 6
- "Tau.csv": outputs from the simulated studies

#### Real_Data
This sub-repository includes R scripts for performing the real data examples in Section 4 and reproducing Tables 3 and 4, which compare the performance of various kinds of low-rank approximation methods.

- "anom1962.RData": the annual total precipitation anomalies data used in Section 4.1
- "Ozone_dat.csv": the level-2 total column ozone data in Section 4.2
- "DP_LRK.R": R script for performing the first real data example and reproducing Table 3
- "DTCO_LRK.R": R script for performing the second real data example and reproducing Table 4

  
## Reproducibility Workflow
Please begin by downloading the entire repository as the file "LRGPSP_Reproducibility_Materials.zip" and extracting it as a folder named "LRGPSP_Reproducibility_Materials". Next, set your working directory to this folder. Then, load necessary R packages and functions. Finally, follow the "Wrapper.R" file to reproduce each figure and table sequentially. The computation time reported below was recorded using (R 3.6.3) running on a machine equipped with Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz and 125 GB RAM.

#### Reproduce Figure 1 in Section 2.1 (and Figure S2 in Section S5)
Figures 1 and S2 demonstrate the patterns of various kinds of representative points under two location sets. Please refer to the "Wrapper.R" file for their reproducibility command. The total computational time is approximately 39.1 minutes, with the majority of the time spent obtaining k-DPPs.  For more detailed code and computational time, please refer to the sub-repository "SP_Demo". 

#### Reproduce Figure 2 and Table 1 in Section 3.1
Figure 2 and Table 1 illustrate the values of smoothness parameter $\gamma$ for the Matern covariance function under various parameter settings. Please refer to the "Wrapper.R" file for their reproducibility command. The total computational time is approximately 2 hours. For more detailed code, outputs, and computational time, please refer to the sub-repository "S1_Gamma". 

#### Reproduce Figure 3 and Table 2 in Section 3.2
Figure 3 and Table 2 demonstrate the influence of rep-points on the predictive process. Please refer to the "Wrapper.R" file for the commands to reproduce these results. There are three options for reproducibility. The first option, set as the default, reproduces Figure 3 and Table 2 using outputs from a single replicate, which takes approximately 16.5 hours. The second option, used in the article, reproduces the results using outputs from 100 replicates but requires significantly more time. The third option uses our pre-generated outputs directly. To select the second or third options, please manually adjust the settings in "Simu2.R". For more detailed code, outputs, and computational time, please refer to the sub-repository "S2_VersusK".

#### Reproduce Figures 4 and 5 in Section 3.3 (and Figure S1 in Section S4)
Figures 4, 5, and S1 are used to verify Theorems 1 and 2. Figures 4 and S1(a) show how the $\log$(MSPE)s of the predictive process predictions change as the sample size $n$ grows. Figure 5 and S1(b) display $\gamma$ obtained by two approaches. Please refer to the "Wrapper.R" file for the commands to reproduce these results. The outputs are based on 100 replicates. Generating outputs for Figures 4 and 5 from a single replicate takes approximately 6.5 minutes. For Figure S1, it takes about 11.9 minutes per replicate. Using 4 cores as specified in "Simu3.R", the entire process will take around 7.7 hours. Readers can adjust the number of cores used or opt to use our pre-generated outputs directly. For more detailed code, outputs, and computational time, please refer to the sub-repository "S3_VersusN".

#### Reproduce Figure 6 in Section 3.3
Figure 6 is used to verify Theorem 3, illustrating how the $\log$(MSPE)s of the predictive process predictions change as the scale of the nugget effect $\tau^2$ increases. For the commands to reproduce these results, please refer to the "Wrapper.R" file. Since it takes only about 1 minute to generate the outputs of a single replicate, we suggest reproducing Figure 6 with outputs from 100 replicates. Using 4 cores in parallel, as specified in "Simu4.R," this process takes about 25 minutes. Readers can adjust the number of cores used or opt to use our pre-generated outputs directly. For more detailed code, outputs, and computational time, please refer to the sub-repository "S4_VersusTau".

#### Reproduce Tables 3 and 4 in Section 4 (and Figure S3 in Section S5)
Tables 3 and 4 compare various low-rank approximation methods using two real datasets. For commands to reproduce these results, please refer to the "Wrapper.R" file. Generating outputs for a single replicate of Tables 3 and 4 takes about 13.5 minutes and 2.25 hours, respectively. Note that the computational times recorded here differ from those in the article due to variations in computational environments. However, these times still accurately reflect the comparative performance of the various methods. For the data, detailed code, and computational time, please refer to the sub-repository "Real_Data". 


## Data 
Two real datasets are utilized in this study to compare various popular low-rank approximation methods and underscore the advantages of the proposed approach, which uses a predictive process with an estimated covariance and support points. They are available as part of the paper’s supplementary material. 

#### Annual total precipitation anomalies data
The first dataset comprises annual total precipitation anomalies observed at 7352 weather stations in the United States in 1962, as discussed by [1] and [2]. It is available online at https://www.image.ucar.edu/Data/precip_tapering/. Readers can access it by loading the “anom1962.RData” file within the sub-repository “Real_Data”. 

#### Total column ozone data
The second dataset consists of 173,405 observations of the level-2 total column ozone for October 1, 1988, along with their locations. The entire total column ozone data were collected and preprocessed by NASA. The dataset used in this article is the same one discussed in [3] and [4]. It was kindly provided by one of the author of [4] and can be accessed by reading the “Ozone_dat.csv” in the sub-repository “Real_Data”. 

## Supporting Software Requirement
#### Version of primary software used
R version 3.6.3 (2020-02-29)

#### Libraries and dependencies used by the code

assertthat_0.2.1, autoFRK_1.4.3, doParallel_1.0.17, fields_13.3, foreach_1.5.2, FRK_2.0.5, geoR_1.9-2, ggplot2_3.3.6, GpGp_0.4.0, here_1.0.1, LatticeKrig_8.4, magrittr_2.0.3, mvtnorm_1.1-3, purrr_0.3.4, rlist_0.4.6.2, sp_1.4-6, Split_1.2, support_0.1.4

## Reference
[1] Kaufman, C. G., Schervish, M. J. and Nychka, D. W. (2008) ‘Covariance tapering for likelihood-based estimation in large spatial data sets’, Journal of the American Statistical Association, 103(484), pp. 1545–1555.

[2] Sang, H. and Huang, J.Z. (2012). 'A full scale approximation of covariance functions for large spatial data sets', Journal of the Royal Statistical Society: Series B (Statistical Methodology), 74(1), pp. 111–132.

[3] Cressie, N. and Johannesson, G. (2008), 'Fixed rank kriging for very large spatial data sets', Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(1), pp. 209-226.

[4] Meng, C., Zhang, X., Zhang, J., Zhong, W. and Ma, P. (2020), 'More efficient approximation of smoothing splines via space-filling basis selection', Biometrika, 107(3), pp. 723–735. 





