# JMdiscfrail
Code for implementing a joint model of recurrent and terminal events with discretely-distributed non-parametric frailty.


### Reference
Masci C, Spreafico M & Ieva F (2023). Joint modelling of recurrent and terminal events with discretely-distributed non-parametric frailty: application on re-hospitalizations and death in heart failure patients. Available soon as *arXiv* pre-print


### Data Availability
We cannot provide the original administrative data due to confidentiality.
We provide some fake datasets in order to allow researchers who want to replicate the same analysis to properly get how the code has to be run and how results are displayed and should be read.


## Description

- Files:
  - **main_Gauss.R**: Fit JMdiscfrail with Gaussian initialization. Display results.
  - **main_Unif.R**: Fit JMdiscfrail with Uniform initialization. Display results.
    
- Sub-folder **./functions/** contains some auxiliary functions to run the main files:
  - **JMdiscfrail.R**: Function that implements a joint model of recurrent and terminal events with discretely-distributed non-parametric frailty.
  - **stratified_baseline_surv.R**: Function that computes th stratified survival probability curves for recurrent and terminal processes from a joint model fitted using JMdiscfrail().
  - **log_lik_fun.R**: Functions for computing (classification) log-likelihood. Called in "JMdiscfrail.R".
  - **data_format.R**: Function for re-formatting dataset.
    
- Sub-folder **./data/** contains fake dataset along with their legends:
	- **data_legend.txt**: Variables legend of dataset 'fake_dataRD.Rdata'.
	- **fake_dataRD.Rdata**: Datasets related to 200 fake patients for recurrent (dataR) and terminal (dataD) events (see Supplementary Material A for further details).

## Software
- R software.
- Packages: data.table, MASS, mvtnorm, survival.
  
(Last update: October 24th, 2023)
