# Simulation Study
Code for implementing the simulation study in Section 5.

## Description
- **01_data_generation.R**
  generates data as described in Manuscript Section 5.1, under setting `A` (yielding approximately 80% death and 20% administrative censoring). Generated data are stored in sub-folder **./sim_data_A/**.

  To generate data under setting B (yielding approximately 50% death and 50% administrative censoring), modify the value of `setting` to `B` (lines 33-34) and run it. Generated data are stored in sub-folder **./sim_data_B/**.
* * * 
- **02_fit_jm_bivnorm_ng.R**
  estimates joint frailty models according to Ng et al. for data/scenarios generated under setting A. Results are stored in sub-folder **./sim_results_ng_A/**.
  When possible, this code runs in parallel using either the number of available cores minus 2, or up to a maximum of 40 cores (see `n.cores`). See "Execution time" below.
 
  To estimate Ng et al.'s models for data generated under setting B:
  1. First, modify the value of `setting` in *01_data_generation.R* to `B` and run it to generate data under setting B.
  2. Then, modify lines 35-36 in this file to set `setting` to `B` and run it.
   
  Results are stored in sub-folder **./sim_results_ng_B/**.
* * * 
- **02_fit_jmdf_gauss.R** and **02_fit_jmdf_unif.R**
  estimate JMDFs using Gaussian and Uniform initialization, respectively, for data/scenarios generated under setting A, initialization (ii) and threshold L = 1.5. Results are stored in sub-folder **./sim_results_A/**.
  When possible, these codes run in parallel using either the number of available cores minus 2, or up to a maximum of 20 cores (see `n.cores`). See "Execution time" below.

  To select different initialization (i or ii) and threshold L, modify the value of `folder` (lines 62-65: `I_L1`, `II_L1`, `II_L15`, or `II_L2`). 
  You must run the JMDF files separately for each initialization to get results corresponding to all choices.
   
  To estimate JMDFs for data generated under setting B:
  1. First, modify the value of `setting` in *01_data_generation.R* to `B` and run it to generate data under setting B;
  2. Then, modify lines 43-44 in this file to set `setting` to `B`;
  3. Select the desired initialization (i or ii) and threshold L (`folder`) and run it.
   Results are stored in sub-folder **./sim_results_B/**.
* * * 
- **03_sim_results_manuscript.R** reproduces all Tables and Figures in Section 5.2 of the Manuscript.
* * * 
- **04_sim_results_suppmatS1.R** reproduces all Tables and Figures in Supplementary Material S1.
   
  Before running this file, modify the value of `folder` in *02_fit_jmdf_gauss.R* and *02_fit_jmdf_unif.R* to `I_L1`, `II_L1`, and `II_L2`
  This selects the appropriate initialization (i or ii) and the threshold L. After making this change, run those files to obtain the results using different initializations. You must run the files separately for each initialization to get results  corresponding to all choices.
* * * 
- **04_sim_results_suppmatS2.R** reproduces all Tables and Figures in Supplementary Material S2.

  Before running this file:
  1. First, modify the value of `setting` in *01_data_generation.R* to `B` and run it to generate data under setting B;
  2. Then, modify the value of `setting` in *02_fit_jm_bivnorm_ng.R*, *02_fit_jmdf_gauss.R*, and *02_fit_jmdf_unif.R* to `B` and run those files  to obtain the results under setting B.
     For JMDF, select the desired initialization (i or ii) and the threshold L by modifying the value of `folder` (`I_L1`, `II_L1`, `II_L15`, or `II_L2`). You must run the JMDF files separately for each initialization to get results corresponding to all choices.
* * *   
- Sub-folder **./sim_functions/** contains some auxiliary functions to run the main files of the simulation study.   
- Sub-folders **./sim_data_*/** contain the generated data.
- Sub-folders **./sim_results_ng_*/** and **./sim_results_*/** contain the obtained results.

* * * 
## Execution time

| File | Model | Setting | Scenarios | Init, L | Execution time | Parallel cores |
|------|-------|---------| ----------|---------|----------------|----------------|
| 02_fit_jm_bivnorm_ng.R | JFM by Ng et al. | A | 1-9 | ${\theta_u^2= \theta_v^2}^{(0)}=0.1$, $\rho^{(0)}=0.1$ | 2.06 days | 40 |
| 02_fit_jmdf_gauss.R | JMDF, Gaussian | A | 1-9 | II, L = 1.5 | 25.05 min | 20 |
| 02_fit_jmdf_unif.R | JMDF, Unif | A | 1-9 | II, L = 1.5 | 20.55 min | 20 |


(Last update: November 28th, 2025)
