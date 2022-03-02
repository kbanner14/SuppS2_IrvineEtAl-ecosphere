# README/MetadataS1

The _DataS1_ repository contains all functions (in _R_ directory), simulation results (in _RData_) and a \textsf{R} script file for excecuting simulations and recreating results presented in _Appendix-S2.R_. 

#### _R_ folder

Contents: 

- `setup_code.R` functions for simulating data under _2SppCt_, and provides relevant formats for using those data with the MLE-metric.  
- `mle_britzkeLRT_functions.R` functions for using the MLE-metric at the site and visit levels for the _MLESite_ and _Remove_ approaches, respectively. 
- `simulation_functions.R` Simulation wrapper functions for investigating parameter estimation of _2SppCt_, _Remove_, and _Naive_ approaches, and site-level decisions for _MLESite_. Includes NIMBLE model code for the _2SppCt_ model and a standard single-species occupancy model for fitting _Remove_ and _Naive_. 
- `z-decision-functions.R` includes functions for making site-level decisions about species presence for the Bayesian (_2SppCt_, _Remove_, _Naive_) and _MLESite_ approaches. 

#### _RData_ folder 

Contents: simulation results from all models and all scenarios with the following MCMC sampler settings. 

- Scenarios 1,2,3,4:  16 and 8 visits; niter = 10000, nburn = 5000, thin = 5, nchains = 3

- Scenario 5: required very long burnin and a more aggressive thin 

    - 8 visits: niter = 80000, nburn = 70000, thin = 10, nchains = 3
    - 16 visits: niter = 110000, nburn = 100000, thin = 10, nchains = 3 



