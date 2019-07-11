# This program runs individual level simulation
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)

thetavec = c(0.2, 0, -0.2)
thetaU = 0.3
Nvec = c(5e4, 8e4, 1e5) # 1:3
prop_invalid_vec = c(0.3, 0.5, 0.7)

for (theta in thetavec){
    for (N in Nvec){
        for (prop_invalid in prop_invalid_vec){
            betahat_200sim = list()
            for (repgrp in 1:8){
                load(paste0("betahat_all_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,"_repgrp",repgrp,".rda"))
                betahat_200sim = c(betahat_200sim, betahat_all)
                rm(betahat_all)
            }
            print(length(betahat_200sim))
            save(betahat_200sim, file = paste0("betahat_200sim/betahat_200sim_theta",theta,"_thetaU",thetaU,"_N",N,"_prop_invalid",prop_invalid,".rda"))
            rm(betahat_200sim)
        }
    }
}
