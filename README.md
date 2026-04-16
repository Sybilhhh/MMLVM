# MMLVM
Code for Biometrics paper 'Mixed membership latent variable model with unknown factors, factor loadings and number of extreme profiles'

# Authors
Yuyang HE, The Chinese University of Hong Kong 
Kai KANG, Sun Yat-sen University
Xinyuan SONG, The Chinese University of Hong Kong

The code is written in R and C++ linked via Rcpp, RcppArmadillo and RcppEigen

Prerequisites:
- R version 4.3.0

Requirement: 
libraries: 1. MASS 2. Rcpp; 3. RcppArmadillo; 4. progress; 5. bit; 6. psych; 7. lavaan; 8. LaplacesDemon

# 1. PPMI data analysis

Order selection and Bayesian estimation
File: PPMI_strudy/1_MMLVM/real.R

Run real.r with the current setting.
Output: Txt files such as 'alpha_mean.txt', 'b_mean.txt', and RData file to save working Rdata.

# 2. Simulation

File: simulation/1_sample_size/simulation.R
Run simulation.R with the current setting. It will first generate data for M = 3, Q = 3 using 'generate_data.r' and then estimate the order and parameters therein.

## 2.1 Sample size
Change 'N=500' in line 1 in generate_data.r to consider different sample size. 
Also, change 'P'(line 2),'Q'(line 5), 'R'(line 11), 'b_matrix'(line 29), 'zeta'(line 57) to consider different order and test order performance in simulation 2.
Output: Txt files such as 'alpha_mean.txt', 'b_mean.txt', and RData file to save working Rdata.

## 2.2 different priors considered in Sensitivity analysis 
Assign different priors for parameters. Run simulation.R with the current setting. It will do Bayesian estimation using priorII.

## 2.3 different distributions and measurement errors considered in Sensitivity analysis 
Assign different distributions for 'eps' (line 137) in file 'generate_data.r' to incorporate measurement errors in imaging data
Then run simulation.R

## 2.4 Trace plot
File: simulation/trace_plot/...
3 original MCMC train for order 'Q' and 'M' in Figure 2. Used to test convergence.

