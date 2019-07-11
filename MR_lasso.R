# This function implements the MR-Lasso method
# Code is from paper
# Slob, Eric AW, and Stephen Burgess. "A Comparison Of Robust Mendelian Randomization Methods Using Summary Data." BioRxiv (2019): 577940.
MR_lasso<-function(betaYG,betaXG,sebetaYG){
    betaYGw = betaYG/sebetaYG # dividing the association estimates by sebetaYG is equivalent
    betaXGw = betaXG/sebetaYG # to weighting by sebetaYG^-2
    pleio = diag(rep(1, length(betaXG)))
    l1grid = c(seq(from=0.1, to=5, by=0.1), seq(from=5.2, to=10, by=0.2))
    # values of lambda for grid search
    l1grid_rse = NULL; l1grid_length = NULL; l1grid_beta = NULL; l1grid_se = NULL
    for (i in 1:length(l1grid)) {
        l1grid_which = which(attributes(penalized(betaYGw, pleio, betaXGw, lambda1=l1grid[i], trace=FALSE))$penalized==0) 
        l1grid_rse[i] = summary(lm(betaYG[l1grid_which]~betaXG[l1grid_which]-1, weights=sebetaYG[l1grid_which]^-2))$sigma 
        l1grid_length[i] = length(l1grid_which)
                                                                                             
        l1grid_beta[i] = lm(betaYG[l1grid_which]~betaXG[l1grid_which]-1, weights=sebetaYG[l1grid_which]^-2)$coef[1]
        l1grid_se[i] = summary(lm(betaYG[l1grid_which]~betaXG[l1grid_which]-1, 
                                  weights=sebetaYG[l1grid_which]^-2))$coef[1,2]/min(summary(lm(betaYG[l1grid_which]~betaXG[l1grid_which]-1, weights=sebetaYG[l1grid_which]^-2))$sigma, 1)
    }
    l1which_hetero = c(which(l1grid_rse[1:(length(l1grid)-1)]>1& diff(l1grid_rse)>qchisq(0.95, df=1)/l1grid_length[2:length(l1grid)]), length(l1grid))[1]
    # heterogeneity criterion for choosing lambda
    l1hetero_beta = l1grid_beta[l1which_hetero]
    l1hetero_se = l1grid_se[l1which_hetero] 
    list(ThetaEstimate=l1hetero_beta, ThetaSE=l1hetero_se )
}