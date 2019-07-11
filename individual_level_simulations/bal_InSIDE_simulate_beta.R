# Individual-level simulations with balanced pleiotropy and InSIDE assumption satisfied
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)

thetavec = c(0.2, 0, -0.2)
thetaUvec = c(0.3, 0.5)
Nvec = c(5e4, 8e4, 1e5)
prop_invalid_vec = c(0.3, 0.5, 0.7)

temp = as.integer(commandArgs(trailingOnly = TRUE)) # length=5
theta = thetavec[temp[1]] # True causal effect from X to Y
thetaU = thetaUx = thetaUvec[temp[2]] # Effect of the confounder on Y/X
N = Nvec[temp[3]] # Sample size for exposure X
prop_invalid = prop_invalid_vec[temp[4]] # Proportion of invalid IVs
repgrp = temp[5] # Simulation group number (1:8) - 200 simulations are partitioned into 8 groups for parallelization.

pthr = 5e-8 # p-value threshold for instrument selection
NxNy_ratio = 2 # Ratio of sample sizes for X and Y
M = 2e5 # Total number of independent SNPs representing the common variants in the genome

# Model parameters for effect size distribution
pi1=0.02*(1-prop_invalid); pi3=0.01
pi2=0.02*prop_invalid; 
sigma2x = sigma2y = 5e-5; sigma2u = 0; sigma2x_td = sigma2y_td = (5e-5)-thetaU*thetaUx*sigma2u

print(paste("N", N, "pthr", pthr, "pi1", pi1, "theta", theta, "thetaU", thetaU, "NxNy_ratio", NxNy_ratio))

# Sample size for X and Y
nx = N; ny = N/NxNy_ratio

# Due to large memory requirement of simulating genetic data for all the subjects at once,
# we simulate a batch of 500 subjects at a time, analyze each batch separately, 
# and use meta-analysis to get final estimate.
batch_size = 500 
betahat_all = vector("list", length = 25)

for (repind in 1:25){
    set.seed(89781*(repind+25*(repgrp-1)))
    
    # Generate indices of causal SNPs
    ind1 = sample(M, round(M*pi1))
    causalsnps = ind1
    ind2 = sample(setdiff(1:M,causalsnps), round(M*pi2))
    causalsnps = c(causalsnps,ind2)
    ind3 = sample(setdiff(1:M,causalsnps), round(M*pi3))
    causalsnps = c(causalsnps,ind3)
    
    # Simulate effect size
    gamma = phi = alpha = rep(0,M)
    
    gamma[ind1] = rnorm(length(ind1), mean = 0, sd = sqrt(sigma2x))
    gamma[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2x_td))
    alpha[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2y_td))
    alpha[ind3] = rnorm(length(ind3), mean = 0, sd = sqrt(sigma2y))
    
    betahat_x = betahat_y = rep(0,M)
    
    # Simulate data for exposure X and fit linear regression and meta-analysis to get summary statistics
    nbatch = nx/batch_size
    for (batch_ind in 1:nbatch){
        if (batch_ind%%5==0) print(batch_ind)
        # load(paste0("../geno/geno",batch_ind,".rda"))
        G = matrix(rbinom(batch_size*M, size=2, prob=0.3), ncol = M)
        G = (G-2*0.3)/sqrt(2*0.3*0.7)
        U = G %*% phi + rnorm(batch_size, mean = 0, sd = sqrt(1-M*pi2*sigma2u))
        X = G %*% gamma + thetaUx*U + rnorm(batch_size, mean=0, sd = sqrt(1-thetaUx^2-M*(pi1*sigma2x+pi2*sigma2x_td)))
        betahat_x = betahat_x + t(G)%*%X/batch_size
        rm(G,U,X)
    }
    betahat_x = betahat_x/nbatch
    
    # Simulate data for outcome Y and fit linear regression and meta-analysis to get summary statistics
    nbatch = ny/batch_size
    for (batch_ind in 1:nbatch){
        if (batch_ind%%5==0) print(batch_ind)
        # load(paste0("../geno/geno",200+batch_ind,".rda"))
        G = matrix(rbinom(batch_size*M, size=2, prob=0.3), ncol = M)
        G = (G-2*0.3)/sqrt(2*0.3*0.7)
        U = G %*% phi + rnorm(batch_size, mean = 0, sd = sqrt(1-M*pi2*sigma2u))
        X = G %*% gamma + thetaUx*U + rnorm(batch_size, mean=0, sd = sqrt(1-thetaUx^2-M*(pi1*sigma2x+pi2*sigma2x_td)))
        Y = G %*% alpha + theta*X + thetaU*U + rnorm(batch_size, mean=0, sd = sqrt(1-theta^2-thetaU^2-2*theta*thetaU*thetaUx-M*(pi3*sigma2y+pi2*sigma2y_td)))
        betahat_y = betahat_y + t(G)%*%Y/batch_size    
        rm(G,U,X,Y)
    }
    betahat_y = betahat_y/nbatch
    
    # Filter the SNPs that reach genome-wide significance in the study associated with X.
    ind_filter = which(2*pnorm(-sqrt(nx)*abs(betahat_x))<pthr)
    betahat.flt = cbind(betahat_x[ind_filter], betahat_y[ind_filter])
    colnames(betahat.flt) = c("betahat_x","betahat_y")
    
    print(paste("rep",repind))
    betahat_all[[repind]] = betahat.flt
    save(betahat_all, file = paste0("betahat_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))
}
