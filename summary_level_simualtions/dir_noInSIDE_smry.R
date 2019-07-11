# Summary level simulations with directional pleiotropy and InSIDE assumption violated
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(mr.raps)
library(MRPRESSO)
library(penalized)
source("MR_lasso.R")

thetavec = c(0.2, 0, -0.2)
thetaUvec = c(0.3, 0.5)
Nvec = c(5e4, 8e4, 1e5, 1.5e5, 2e5, 5e5, 1e6) # 1:7
prop_invalid_vec = c(0.3, 0.5, 0.7)

temp = as.integer(commandArgs(trailingOnly = TRUE))
theta = thetavec[temp[1]] # True causal effect from X to Y
thetaU = thetaUx = thetaUvec[temp[2]] # Effect of the confounder on Y/X
N = Nvec[temp[3]] # Sample size for exposure X
prop_invalid = prop_invalid_vec[temp[4]] # Proportion of invalid IVs

pthr = 5e-8 # p-value threshold for instrument selection
NxNy_ratio = 2 # Ratio of sample sizes for X and Y
M = 2e5 # Total number of independent SNPs representing the common variants in the genome

# Model parameters for effect size distribution
pi1=0.02*(1-prop_invalid); pi3=0.01
pi2=0.02*prop_invalid; 
sigma2x = sigma2y = 5e-5; sigma2u = 1e-4; sigma2x_td = sigma2y_td = (5e-5)-thetaU*thetaUx*sigma2u

print(paste("N", N, "pthr", pthr, "pi1", pi1, "theta", theta, "thetaU", thetaU, "prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio))

nx = N; ny = N/NxNy_ratio
est = matrix(NA, nrow = 200, ncol = 3+10*3)
mr_methods = c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS")
colnames(est) = c("numIV", "varX_expl","varY_expl", mr_methods, paste0(mr_methods,"_se"), paste0(mr_methods,"_time"))
# varX_expl: variance of X explained by IVs; varY_expl: variance of Y explained by IVs

for (repind in 1:200){
    set.seed(4365*repind)
    
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
    phi[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2u))
    alpha[ind2] = rnorm(length(ind2), mean = 0.005, sd = sqrt(sigma2y_td))
    alpha[ind3] = rnorm(length(ind3), mean = 0.005, sd = sqrt(sigma2y))
    
    # Generate summary statistics directly from summary-level model implied by individual-level model
    betax = gamma + thetaUx*phi
    betay = alpha + theta*betax + thetaU*phi
    betahat_x = betax + rnorm(M, mean = 0, sd = sqrt(1/nx))
    betahat_y = betay + rnorm(M, mean = 0, sd = sqrt(1/ny))
    
    # Filter the SNPs that reach genome-wide significance in the study associated with X.
    ind_filter = which(2*pnorm(-sqrt(nx)*abs(betahat_x))<pthr)
    numIV = length(ind_filter)
    est[repind,1] = numIV
    est[repind,2] = sum(betax[ind_filter]^2)
    est[repind,3] = sum(betay[ind_filter]^2)
    
    # MR analysis with 10 methods
    if (numIV>2){
        betahat_x.flt = betahat_x[ind_filter]
        betahat_y.flt = betahat_y[ind_filter]
        
        mr.obj = mr_input(bx = betahat_x.flt, bxse = rep(1/sqrt(nx), length(betahat_x.flt)),
                          by = betahat_y.flt, byse = rep(1/sqrt(ny), length(betahat_y.flt)))
        # 1. IVW
        T0 = proc.time()[3]
        res = mr_ivw(mr.obj)
        T1 = proc.time()[3]
        est[repind,3+1] = res$Estimate
        est[repind,3+11] = res$StdError
        est[repind,3+21] = T1-T0
        rm(res)
        
        # 2. median
        T0 = proc.time()[3]
        res = mr_median(mr.obj)
        T1 = proc.time()[3]
        est[repind,3+2] = res$Estimate
        est[repind,3+12] = res$StdError
        est[repind,3+22] = T1-T0
        rm(res)
        # 3. mode
        T0 = proc.time()[3]
        res = mr_mbe(mr.obj)
        T1 = proc.time()[3]
        est[repind,3+3] = res$Estimate
        est[repind,3+13] = res$StdError
        est[repind,3+23] = T1-T0
        rm(res)
        # 4.MR-PRESSO
        presso.df = data.frame(bx = betahat_x.flt, by = betahat_y.flt, bxse = rep(1/sqrt(nx),length(betahat_x.flt)), byse = rep(1/sqrt(ny), length(betahat_y.flt)))
        if (N<=2e5){
            T0 = proc.time()[3]
            res = mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "byse", SdExposure = "bxse", 
                            OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = presso.df, NbDistribution = 3000, SignifThreshold = 0.05)
            T1 = proc.time()[3]
            if (!is.na(res$`Main MR results`[2,"Causal Estimate"]) & !is.na(res$`Main MR results`[2,"Sd"])){
                est[repind,3+4] = res$`Main MR results`[2,"Causal Estimate"]
                est[repind,3+14] = res$`Main MR results`[2,"Sd"]
            } else{
                est[repind,3+4] = res$`Main MR results`[1,"Causal Estimate"]
                est[repind,3+14] = res$`Main MR results`[1,"Sd"]
            }
            est[repind,3+24] = T1-T0
            rm(res)
        }
        # 5. robust
        T0 = proc.time()[3]
        res = mr_ivw(mr.obj,"random", robust = TRUE)
        T1 = proc.time()[3]
        est[repind,3+5] = res$Estimate
        est[repind,3+15] = res$StdError
        est[repind,3+25] = T1-T0
        rm(res)
        # 6. MR-Lasso
        T0 = proc.time()[3]
        res = MR_lasso(presso.df$by, presso.df$bx, presso.df$byse)
        T1 = proc.time()[3]
        est[repind,3+6] = res$ThetaEstimate
        est[repind,3+16] = res$ThetaSE
        est[repind,3+26] = T1-T0
        rm(res)
        # 7. Egger
        T0 = proc.time()[3]
        res = mr_egger(mr.obj)
        T1 = proc.time()[3]
        est[repind,3+7] = res$Estimate
        est[repind,3+17] = res$StdError.Est
        est[repind,3+27] = T1-T0
        rm(res)
        # 8. contamination mixture
        T0 = proc.time()[3]
        res = mr_conmix(mr.obj)
        T1 = proc.time()[3]
        est[repind,3+8] = res$Estimate
        CIlength = res$CIUpper-res$CILower
        if (length(CIlength)>1) print(paste("Repind",repind,"conmix multimodal"))
        est[repind,3+18] = sum(CIlength)/1.96/2 ## Caution: this may be problematic
        est[repind,3+28] = T1-T0
        rm(res)
        # 9. MRMix
        theta_temp_vec = seq(-0.5,0.5,by=0.01)
        T0 = proc.time()[3]
        res = MRMix(betahat_x.flt, betahat_y.flt, sx=1/sqrt(nx), sy=1/sqrt(ny), theta_temp_vec, pi_init = 0.6, sigma_init = 1e-5)
        T1 = proc.time()[3]
        est[repind,3+9] = res$theta
        est[repind,3+19] = res$SE_theta
        est[repind,3+29] = T1-T0
        rm(res)
        # 10. MR-RAPS
        T0 = proc.time()[3]
        res = mr.raps.overdispersed.robust(presso.df$bx, presso.df$by, presso.df$bxse, presso.df$byse,
                                           loss.function = "huber", k = 1.345, initialization = c("l2"), 
                                           suppress.warning = FALSE, diagnostics = FALSE, niter = 20, tol = .Machine$double.eps^0.5)
        T1 = proc.time()[3]
        est[repind,3+10] = res$beta.hat
        est[repind,3+20] = res$beta.se
        est[repind,3+30] = T1-T0
        rm(res)
    }
    if (repind%%5==0){
        print(paste("Rep",repind,"numIV",numIV))
        save(est, file = paste0("est_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))
    }
}
save(est, file = paste0("est_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))

