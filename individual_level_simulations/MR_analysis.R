# Analyze summary statistics generated from individual level simulations
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(mr.raps)
library(MRPRESSO)
library(penalized)
source("/dcl01/chatterj/data/gqi/multitrait_effsize/code/MR_lasso.R")

thetavec = c(0.2, 0, -0.2)
Nvec = c(5e4, 8e4, 1e5)
prop_invalid_vec = c(0.3, 0.5, 0.7)

temp = as.integer(commandArgs(trailingOnly = TRUE))
theta = thetavec[temp[1]] # True causal effect from X to Y
thetaU = thetaUx = 0.3 # Effect of the confounder on Y/X
N = Nvec[temp[2]] # Sample size for exposure X
prop_invalid = prop_invalid_vec[temp[3]] # Proportion of invalid IVs

pthr = 5e-8
NxNy_ratio = 2
M = 2e5

print(paste("N", N, "pthr", pthr, "theta", theta, "thetaU", thetaU, "prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio))

nx = N; ny = N/NxNy_ratio
est = matrix(NA, nrow = 200, ncol = 3+10*3+2)
mr_methods = c("IVW", "median", "mode", "PRESSO", "Robust", "Lasso", "egger", "conmix", "MRMix", "RAPS")
colnames(est) = c("numIV", "varX_expl","varY_expl", mr_methods, paste0(mr_methods,"_se"), paste0(mr_methods,"_time"), "PRESSO_pval", "conmix_CIcover0")
# varX_expl: variance of X explained by IVs; varY_expl: variance of Y explained by IVs

load(paste0("betahat_200sim/betahat_200sim_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))
for (repind in 1:200){
    set.seed(8967*repind)
    numIV = nrow(betahat_200sim[[repind]])
    est[repind,1] = numIV
    
    if (numIV>2){
        betahat_x.flt = as.vector(betahat_200sim[[repind]][,1])
        betahat_y.flt = as.vector(betahat_200sim[[repind]][,2])
        
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
        if (nx<=2e5){
            T0 = proc.time()[3]
            res = try(mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "byse", SdExposure = "bxse", 
                                OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = presso.df, NbDistribution = 3000, SignifThreshold = 0.05))
            T1 = proc.time()[3]
            if (class(res)!="try-error"){
                if (!is.na(res$`Main MR results`[2,"Causal Estimate"]) & !is.na(res$`Main MR results`[2,"Sd"])){
                    est[repind,3+4] = res$`Main MR results`[2,"Causal Estimate"]
                    est[repind,3+14] = res$`Main MR results`[2,"Sd"]
                    est[repind,34] = res$`Main MR results`[2,"P-value"]
                } else{
                    est[repind,3+4] = res$`Main MR results`[1,"Causal Estimate"]
                    est[repind,3+14] = res$`Main MR results`[1,"Sd"]
                    est[repind,34] = res$`Main MR results`[1,"P-value"]
                }
                est[repind,3+24] = T1-T0
            } else{
                print(paste("Rep",repind,"numIV",numIV))
                print(res)
            }
            rm(res)
        }
        # 5. MR-Robust
        T0 = proc.time()[3]
        res = mr_ivw(mr.obj,"random", robust = TRUE)
        T1 = proc.time()[3]
        est[repind,3+5] = res$Estimate
        est[repind,3+15] = res$StdError
        est[repind,3+25] = T1-T0
        rm(res)
        # 6. MR-Lasso
        T0 = proc.time()[3]
        res = try(MR_lasso(presso.df$by, presso.df$bx, presso.df$byse))
        T1 = proc.time()[3]
        if (class(res)!="try-error"){
            est[repind,3+6] = res$ThetaEstimate
            est[repind,3+16] = res$ThetaSE
            est[repind,3+26] = T1-T0
        } else{
            print(paste("Rep",repind,"numIV",numIV))
            print(res)
        }
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
        est[repind,35] = ifelse(sum((res$CILower<=0)&(res$CIUpper>=0))>0, 1, 0)
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

