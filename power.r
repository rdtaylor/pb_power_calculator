#!/usr/bin/Rscript
#
# power.r
#
# Simulates Plackett-Burman designs in order to determine a conservative
# confidence interval for power


require('FrF2', quietly=TRUE)
require('foreach', quietly=TRUE)
require('doMC', quietly=TRUE)


getConfigurations <- function() {

    n_cores_in_cpu <- 2

    n_factors <- 24
    # Use the maximum number of runs by default given n_factors
    n_runs <- n_factors + (4 - (n_factors %% 4))
    n_reps <- 1
    n_sims <- 10000
    design <- pb(nruns=n_runs, nfactors=n_factors, replications=n_reps)
    n <- nrow(design)

    alpha <- 0.05
    ybar <- 856
    sigma <- 19

    explanatory_vars <- colnames(design)
    effect_size <- 31
    effect_column <- 5

    configs <- Configurations(N_CORES_IN_CPU=n_cores_in_cpu,
                              N_FACTORS=n_factors,
                              N_RUNS=n_runs,
                              N_REPS=n_reps,
                              N_SIMS=n_sims,
                              ALPHA=alpha,
                              YBAR=ybar,
                              SIGMA=sigma,
                              EXPLANATORY_VARS=explanatory_vars,
                              EFFECT_SIZE=effect_size, 
                              EFFECT_COLUMN=effect_column,
                              DESIGN=design,
                              N=n)

    return(configs)
}


configClosure <- function(configs) {
    # Required by R to ensure immutability of the object
    configs <- configs

    getSetting <- function(element) {
        config <- getElement(configs, element)
        return(config)
    }
}


main <- function() {
    # Use parallel execution for the simulations
    # (registerDoMC must be done somewhere within scope of the foreach loop)
    registerDoMC(cores=configs("N_CORES_IN_CPU"))

    simulated_power <- runSims()
    conf_int <- calcConfInt(simulated_power)

    report <- getReport(conf_int)
    printReport(report)
}


runSims <- function() {
    effect_addition <- getEffectAdd()

    p_values <- foreach(i = seq(1, configs("N_SIMS")),
                        .combine='c',
                        .inorder=FALSE) %dopar% {
        runSim(effect_addition)
    }

    simulated_power <- getNullRejectedPct(p_values)
    return(simulated_power)
}


getEffectAdd <- function() {
    effect_addition <- rep(0, configs("N"))
    affected_rows <- which(configs("DESIGN")[, configs("EFFECT_COLUMN")] == 1)
    effect_addition[affected_rows] <- configs("EFFECT_SIZE")

    return(effect_addition)
}


runSim <- function(effect_addition) {
    resp <- getResponse(effect_addition)
    design_aov <- getANOVA(resp)
    pvalue <- getPFromANOVA(design_aov)

    return(pvalue)
}


getResponse <- function(effect_addition) {
    resp <- rnorm(n=configs("N"), mean=configs("YBAR"), sd=configs("SIGMA"))
    resp <- resp + effect_addition

    return(resp)
}


getANOVA <- function(resp) {
    design_formula <- paste0('resp ~ ',
                             paste(configs("EXPLANATORY_VARS"),
                             collapse=' + '))
    design_aov <- aov(as.formula(design_formula), data=configs("DESIGN"))

    return(design_aov)
}


getPFromANOVA <- function(design_aov) {
    aov_summary <- summary(design_aov)[[1]]
    pvalues <- aov_summary[[5]]
    pvalue <- pvalues[configs("EFFECT_COLUMN")]
    return(pvalue)
}


getNullRejectedPct <- function(p_values) {
    simulated_power <- mean(p_values <= configs("ALPHA"))
    return(simulated_power)
}


calcConfInt <- function(simulated_power) {
    Z <- qnorm(1 - configs("ALPHA")/2)
    SE <- sqrt(simulated_power * (1 - simulated_power) / configs("N_SIMS"))
    conf_int <- c(simulated_power - Z*SE, simulated_power + Z*SE)

    return(conf_int)
}


getReport <- function(conf_int) {
    conf_int_str <- paste0(conf_int, collapse=' ')

    report <- paste0((1 - configs("ALPHA")) * 100,
                     "% confidence interval for the power given ",
                     "alpha = ", configs("ALPHA"),
                     ", significant effects being ", configs("EFFECT_SIZE"),
                     " connections before errors, ", configs("N_REPS"),
                     " replications, and ", configs("N_SIMS"),
                     " simulations: \n",
                     conf_int_str, "\n")
    return(report)
}

printReport <- function(report) {
    cat(report)
}


# Needed to use FrF2::design objects in S4 objects as below
setOldClass("design") 

setClass("Configurations",
         representation(N_CORES_IN_CPU="numeric",
                        N_FACTORS="numeric",
                        N_RUNS="numeric",
                        N_REPS="numeric",
                        N_SIMS="numeric",
                        ALPHA="numeric",
                        YBAR="numeric",
                        SIGMA="numeric",
                        EXPLANATORY_VARS="character",
                        EFFECT_SIZE="numeric",
                        EFFECT_COLUMN="numeric",
                        DESIGN="design",
                        N="numeric"
                       ))

Configurations <- function(N_CORES_IN_CPU,
                           N_FACTORS,
                           N_RUNS,
                           N_REPS,
                           N_SIMS,
                           ALPHA,
                           YBAR,
                           SIGMA,
                           EXPLANATORY_VARS,
                           EFFECT_SIZE,
                           EFFECT_COLUMN,
                           DESIGN,
                           N) {

    new("Configurations", N_CORES_IN_CPU=N_CORES_IN_CPU,
                          N_FACTORS=N_FACTORS,
                          N_RUNS=N_RUNS,
                          N_REPS=N_REPS,
                          N_SIMS=N_SIMS,
                          ALPHA=ALPHA,
                          YBAR=YBAR,
                          EXPLANATORY_VARS=EXPLANATORY_VARS,
                          SIGMA=SIGMA,
                          EFFECT_SIZE=EFFECT_SIZE,
                          EFFECT_COLUMN=EFFECT_COLUMN,
                          DESIGN=DESIGN,
                          N=N)
}

# Start by instatiating the Configurations class and provide a global closure
# for easy access to configurations.
configs <- getConfigurations()
configs <- configClosure(configs)

main()
