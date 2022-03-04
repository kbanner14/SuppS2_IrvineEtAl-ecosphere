# This R script goes through all code to run simulations (see warning)
# and recreate figures Appendix S2 using saved .rds files 
# from the simulations run by the authors.


# load packages necessary for data generation and 
# simulation functions 
library(nimble)
library(coda) # for effective size function
library(rstan) # for Rhat function
options(mc.cores = parallel::detectCores())
# load data generation and simulation functions
source("R/setup_code.R")
source("R/mle_britzkeLRT_functions.R")
source("R/simulation_functions.R") # includes NIMBLE Model code
source("R/z-decision-functions.R")

# WARNING: Careful tuning of the MCMC prior to running 
# simulations functions is required to asses convergence.
# An example of tuning is provided for S1 below 

# set arguments to S1
thetas <- matrix(c(.9,.1,.35,.65),
                nrow=2,
                ncol=2,byrow=TRUE)

# 55 sites, 16 visits, psi, lambda and alpha settings 
n.sites1 <-  5
n.sites2 <- 50
n.visits <- 16
psi.vec = c(.25,.75)
lambda.vec = c(0.3,10) 
alpha_vec = c(0.05, 0.1, 0.15, 0.2)

# simulate one data set
sim_dat <- sim_data(n.sites1 = n.sites1,
                    n.sites2 = n.sites2,
                    n.visits = n.visits,
                    psi.vec = psi.vec,
                    lambda.vec = lambda.vec, 
                    theta.mat = thetas)

n.sites.tot <- length(unique(sim_dat$MLE_Britzke_data$sites_ID))

spp1 <- matrix(sim_dat$MLE_Britzke_data$Spp1_AutoIDs,
               nrow =  n.sites.tot, ncol = n.visits, byrow=TRUE)
spp2 <- matrix(sim_dat$MLE_Britzke_data$Spp2_AutoIDs,
               nrow =  n.sites.tot, ncol = n.visits, byrow=TRUE)

# nimble model data
cd_data <- list(Spp1_AutoIDs = spp1,
                Spp2_AutoIDs = spp2,
                theta = thetas)

# nimble model constants
cd_constants <- list(n.sites =  n.sites.tot, n.revisits =  n.visits, n.spp = 2)

# initial values
cd_inits <- function(spp1, spp2) {
  naive_z1 <- as.numeric(apply(spp1, 1, sum) > 0)
  naive_z2 <- as.numeric(apply(spp2, 1, sum) > 0)
  Z_init <- cbind(naive_z1, naive_z2)
  psi_init <- runif(2)
  lambda_init <- rgamma(2, 2, .25)
  list(Z = Z_init,
       lambda = lambda_init,
       psi = psi_init)
}

# fit the model in one-step
cd_mcmc_simple <- nimbleMCMC(code = cd_nimble_ThetaKnown,
                             data = cd_data,
                             constants = cd_constants,
                             inits = cd_inits(spp1 = spp1, spp2 = spp2),
                             monitors = c("psi","lambda", "Z"),
                             thin = 5,
                             niter = 10000,
                             nburn = 5000, 
                             nchains = 3,
                             samplesAsCodaMCMC = TRUE)

# play around with thin, niter, nburn until convergence 
# looks decent on multiple iterations of fake data!
plot(cd_mcmc_simple[,111:114])
coda::traceplot(cd_mcmc_simple[,111])
coda::traceplot(cd_mcmc_simple[,112])

tmp <- summary(cd_mcmc_simple)
out <- as_tibble(cbind(tmp$statistics, tmp$quantiles)) %>%
  mutate(
    param = rownames(tmp$statistics),
    model = "2SppCt"
  ) %>%
  dplyr::select(param, model, everything()) %>%
  mutate(
    truth = c(sim_dat$MLE_Britzke_data$Z1, sim_dat$MLE_Britzke_data$Z2, lambda.vec, psi.vec),
    n_eff = coda::effectiveSize(cd_mcmc_simple)
  )

out$n_eff

mat <- matrix(NA, nrow = nrow(out), ncol = )
colnames(mat) <- c("Rhat")
for(a in 1:nrow(out)){
  tmp <- sapply(cd_mcmc_simple, function(x) x[,a])
  mat[a,] <- c(rstan::Rhat(tmp))
}
out$Rhat <- c(mat)     

out %>% filter(param == "lambda[1]"| param == "lambda[2]") %>% data.frame

# After tuning ALL scenarios, run the simulation. 
# WARNING: running simulations takes between 2-4 hours 
# on a laptop with 2.3 GHz 8-Core Intel Core i9 processor, 
# 16 GB 2667 MHz DDR4 RAM, and AMD Radeon Pro 
#  5500M 4 GB grahpics card.  Uncomment code below to run 
# simulations of your own.
# scenario 1
# start1 <- Sys.time()
# scen1_visits16 <- psi_cov_function(50, thetas = matrix(c(.9,.1,.35,.65),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 16,
#                              psi.vec = c(.25,.75),
#                              lambda.vec = c(0.3,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), 
#                              niter = 10000, nburn = 5000, thin = 5)
# end1 <- start1 - Sys.time()
# 
# # scenario 1; 8 visits
# start1b <- Sys.time()
# scen1_visits8 <- psi_cov_function(50, thetas = matrix(c(.65,.35,.1,.9),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 8,
#                              psi.vec = c(.25,.75),
#                              lambda.vec = c(0.3,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), 
#                              niter = 10000, nburn = 5000, thin = 5)
# end1b <- Sys.time() - start1b  
# 
# # scenario 2
# start2 <- Sys.time()
# scen2_visits16 <- psi_cov_function(50, thetas = matrix(c(.65,.35,.1,.9),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 16,
#                              psi.vec = c(.25,.75),
#                              lambda.vec = c(0.3,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), 
#                              niter = 10000, nburn = 5000, thin = 5)
# end2 <- Sys.time() - start2  
# 
# # S2 8 visits 
# start2b <- Sys.time()
# scen2_visits8 <- psi_cov_function(50, thetas = matrix(c(.65,.35,.1,.9),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 8,
#                              psi.vec = c(.25,.75),
#                              lambda.vec = c(0.3,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), 
#                              niter = 10000, nburn = 5000, thin = 5)
# end2b <- Sys.time() - start2b 
# 
# # scenario 3 
# start3 <- Sys.time()
# scen3_visits16 <- psi_cov_function(50, thetas = matrix(c(.65,.35,.4,.6),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 16,
#                              psi.vec = c(.25,.75),
#                              lambda.vec = c(0.3,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), niter = 10000, 
#                              nburn = 5000, thin = 5)
# end3 <- Sys.time() - start3 
# 
# start3b <- Sys.time()
# scen3_visits8 <- psi_cov_function(50, thetas = matrix(c(.65,.35,.4,.6),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 8,
#                              psi.vec = c(.25,.75),
#                              lambda.vec = c(0.3,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), niter = 10000, 
#                              nburn = 5000, thin = 5)
# end3b <- Sys.time() - start3b 
# 
# # scenario 4 
# start4 <- Sys.time()
# scen4_visits16 <- psi_cov_function(50, thetas = matrix(c(.9,.1,.35,.65),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 16,
#                              psi.vec = c(.5,.75),
#                              lambda.vec = c(0.3,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), 
#                              niter = 10000, nburn = 5000, thin = 5)
# end4 <- Sys.time() - start4
# 
# 
# # S4, 8 visits 
# start4b <- Sys.time()
# scen4_visits8 <- psi_cov_function(50, thetas = matrix(c(.9,.1,.35,.65),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 8,
#                              psi.vec = c(.5,.75),
#                              lambda.vec = c(0.3,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), 
#                              niter = 10000, nburn = 5000, thin = 5)
# end4b <- Sys.time() - start4b
# 
# # scenario 5
# 
# start5 <- Sys.time()
# scen5_visits16 <- psi_cov_function(50, thetas = matrix(c(.9,.1,.35,.65),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 16,
#                              psi.vec = c(.25,.75),
#                              lambda.vec = c(0.7,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2),
#                              niter = 110000, 
#                              nburn = 100000, thin = 10)
# end5 <- start5 - Sys.time()
# 
# start5b <- Sys.time()
# scen5_visits8 <- psi_cov_function(50, thetas = matrix(c(.9,.1,.35,.65),
#                                              nrow=2,
#                                              ncol=2,byrow=TRUE), 
#                              n.sites1 = 5,
#                              n.sites2 = 50,
#                              n.visits = 8,
#                              psi.vec = c(.25,.75),
#                              lambda.vec = c(0.7,10), 
#                              alpha_vec = c(0.05, 0.1, 0.15, 0.2), niter = 80000, 
#                              nburn = 70000, thin = 10)
# end5 <- start5 - Sys.time()
# 
########################
# Recreate Appendix S2 #
########################

# load tidyverse 
library(tidyverse)

# load functions for post processing simulation results
source("R/z-decision-functions.R")

# read in the data
# all .rds objects are in the RData directory, 
# this is located in the same directory this .Rmd file
# is in. 
# 1. set working directory to the Supplement folder 
# 2. read in data
scen1_v16 <- readRDS("RData/scen1_visits16.rds")
scen1_v8 <- readRDS("RData/scen1_visits8.rds")
scen2_v16 <- readRDS("RData/scen2_visits16.rds")
scen2_v8 <- readRDS("RData/scen2_visits8.rds")
scen3_v16 <- readRDS("RData/scen3_visits16.rds")
scen3_v8 <- readRDS("RData/scen3_visits8.rds")
scen4_v16 <- readRDS(file = "RData/scen4_visits16.rds")
scen4_v8 <- readRDS(file = "RData/scen4_visits8.rds")
scen5_v16 <- readRDS(file = "RData/scen5_visits16.rds")
scen5_v8 <- readRDS(file = "RData/scen5_visits8.rds")

# create one big data frame with all scenarios and visits
# create a list with all scenario-night combos. 
# this list will be used for data processing for all figs
sim_list <- list(scen1_v8, scen1_v16, 
                 scen2_v8, scen2_v16, 
                 scen3_v8, scen3_v16, 
                 scen4_v8, scen4_v16, 
                 scen5_v8, scen5_v16)

# convergence processing function. 
# filters out iterations that resulted in MCMC with 
# Rhat > 1.1
comp_check <- function(sim_results, rhat_cutoff = 1.1){
  out <- data.frame()
  # 2sppCt
  rhat <- sim_results$CtDetection %>% filter(param_clean != "Z") %>% 
    dplyr::select(Rhat)
  idx <- which(rhat > rhat_cutoff)
  ct_det <- sim_results$CtDetection %>% filter(param_clean != "Z")
  out <- rbind(out, ct_det[idx,c(1:3,7,11,15,17:20)])
  # naive 
  rhat_n <- sim_results$Naive %>% filter(param_clean != "Z") %>%
    dplyr::select(Rhat)
  idx_n <- which(rhat_n > rhat_cutoff)
  naive  <- sim_results$Naive %>% filter(param_clean != "Z")
  out <- rbind(out, naive[idx_n,c(1:3,7,11,15,17:20)])
  # remove 
  rhat_r <- sim_results$Naive %>% filter(param_clean != "Z") %>%
    dplyr::select(Rhat)
  idx_r <- which(rhat_r > rhat_cutoff)
  naive  <- sim_results$Naive %>% filter(param_clean != "Z")
  out <- rbind(out, naive[idx_r,c(1:3,7,11,15,17:20)])
  return(out)
}

# mcmc convergence check/summary function
converge_summ <- function(sim_list, rhat_cutoff = 1.1, names = 
                            c("N8_S1", "N16_S1", 
                              "N8_S4", "N16_S4",
                              "N8_S5", "N16_S5",
                              "N8_S2", "N16_S2",
                              "N8_S3", "N16_S3")){
  out <- dplyr::tibble()
  for(i in 1:length(sim_list)){
    df_summ <- comp_check(sim_results = sim_list[[i]],
                          rhat_cutoff = 1.1)
    df_summ$trt <- names[i]
    out <- dplyr::bind_rows(out, df_summ)
  }
  df_summ <- out %>% group_by(interaction(trt,param, model)) %>% 
    summarize(group = unique(trt), param = unique(param), 
              model = unique(model), 
              num_unconverged = length(unique(sim_iter)), 
              num_converged = 50 - num_unconverged) %>% 
    mutate(group = factor(group, levels = names))
  df_summ <- df_summ[order(df_summ$group), ]
  return(df_summ[,-1])
}


## ----param_est-setup-----------------------------------------------------
plot_estimation_df <- function(sim_output, 
                               psi.vec = c(0.25, 0.75), 
                               lambda.vec = c(0.3, 10), trt = "N8_S1", 
                               rhat_cutoff = 1.1){
  `%>%` <- magrittr::`%>%`
  model <- c("2sppCt", "Remove", "Naive")
  df_bysim <- data.frame()
  df_plot <- data.frame()
  df_all <- data.frame()
  for(i in 1:3){
    cap_summ <- sim_output[[i]] %>% 
      dplyr::filter(param_clean != "Z") %>%
      dplyr::filter(Rhat <= rhat_cutoff) %>% 
      dplyr::group_by(param) %>% 
      dplyr::summarise(cap_rate = mean(capture,na.rm=TRUE),
                       # avg_se = mean(se_psi_hat,na.rm=TRUE),
                       avg_est = mean(Mean,na.rm=TRUE),
                       avg_lb = mean(`2.5%`,na.rm=TRUE),
                       avg_ub = mean(`97.5%`,na.rm=TRUE),
                       avg_width = mean(cri_width,na.rm=TRUE),
                       sd_width = sd(cri_width,na.rm=TRUE), 
                       n_converged = length(unique(sim_iter))
      )
    cap_summ$model <- model[i]
    if(model[i] != "2sppCt"){
      cap_summ$dg <- c(1-exp(-lambda.vec),psi.vec)
    } else {
      cap_summ$dg <- c(lambda.vec,psi.vec)
    }
    cap_summ$trt <- trt 
    df_plot <- rbind(df_plot, cap_summ)
    ind_df <- sim_output[[i]] %>% 
      dplyr::filter(param_clean != "Z") %>% 
      dplyr::filter(Rhat <= rhat_cutoff) %>% 
      dplyr::select(sim_iter, model, param, Mean, `2.5%`, `25%`, `50%`, 
                    `75%`, `97.5%`, truth, capture, cri_width, Rhat)
    df_bysim <- rbind(df_bysim, ind_df)
    names(ind_df)[10] <- "dg"
    df_join <- left_join(cap_summ, ind_df, by = c("param", "dg"))
    df_all <- rbind(df_all, df_join)
  }
  out <- list('individual_sims' = df_bysim, 
              'overall_summary' = df_plot, 
              'overal_and_indiv' = df_all)
  return(out)
}

# overall average plot
# create data frame for plotting-- need dg values
# for all parameters
psi.vec <- c(0.25, 0.75)
psi_mat <- rbind(psi.vec, psi.vec, 
                 psi.vec, psi.vec, 
                 psi.vec, psi.vec,
                 c(0.5, 0.75), 
                 c(0.5, 0.75), 
                 psi.vec, psi.vec)
lambda_mat <- matrix(c(rep(c(0.3,10), 8), 0.7,10, 0.7,10), 
                     byrow = TRUE, ncol = 2)

# scenario-visit combinations (treatments)
scen_name <- c(rep(paste0("S", 1:5), each = 2))
visits <- c(rep(c("N8", "N16"), 5))
trt_vec <- paste(visits, scen_name, sep = "_")

# create one big data frame for plotting - 
# calls plot_estimation_df
df_all <- data.frame()
for(i in 1:length(sim_list)){
  df_temp <- plot_estimation_df(sim_list[[i]], psi.vec = psi_mat[i,], 
                                lambda.vec = lambda_mat[i, ], trt = trt_vec[i])
  df_all <- rbind(df_all, df_temp$overall_summary)
}

df_all$trt <- factor(df_all$trt, levels = c("N8_S1", "N16_S1", 
                                            "N8_S4", "N16_S4", 
                                            "N8_S5" ,"N16_S5", 
                                            "N8_S2", "N16_S2", 
                                            "N8_S3", "N16_S3"))
df_all$model <- factor(df_all$model, levels = c("2sppCt",
                                                "Remove",
                                                "Naive"))
df_all$visit <- substring(df_all$trt, first = 2, last = 3)
scen <- strsplit(as.character(df_all$trt), split = "_")
scen <- sapply(scen, function(x){x[2]})
df_all$scen <- scen

df_all2 <- data.frame()
for(i in 1:length(sim_list)){
  df_temp <- plot_estimation_df(sim_list[[i]], psi.vec = psi_mat[i,], 
                                lambda.vec = lambda_mat[i, ], 
                                trt = trt_vec[i])
  df_all2 <- rbind(df_all2, df_temp$overal_and_indiv)
}

df_all2$trt <- factor(df_all2$trt, levels = c("N8_S1", "N16_S1", 
                                              "N8_S4", "N16_S4", 
                                              "N8_S5" ,"N16_S5", 
                                              "N8_S2", "N16_S2", 
                                              "N8_S3", "N16_S3"))
df_all2$model <- factor(df_all2$model.x, levels = c("2sppCt",
                                                    "Remove",
                                                    "Naive"))
df_all2$visit <- substring(df_all2$trt, first = 2, last = 3)
scen <- strsplit(as.character(df_all2$trt), split = "_")
scen <- sapply(scen, function(x){x[2]})
df_all2$scen <- scen


## data frame to summarize parameter estimation for plotting
df_psi_sub2 <- df_all2 %>% 
  filter(param == "psi[1]") %>%
  filter(trt %in% c("N8_S1", "N8_S4", "N8_S5", 
                    "N8_S2", "N16_S2", "N8_S3", "N16_S3")) %>% 
  mutate(trt = factor(trt))
df_psi_sub <- df_all %>% 
  filter(param == "psi[1]") %>% 
  filter(trt %in% c("N8_S1", "N8_S4", "N8_S5", 
                    "N8_S2", "N16_S2", "N8_S3", "N16_S3")) %>% 
  mutate(trt = factor(trt))

df_p <- df_all %>% 
  filter(param == "p[1]") %>%
  filter(trt %in% c("N8_S1", "N8_S4", "N8_S5", 
                    "N8_S2", "N16_S2", "N8_S3", "N16_S3")) %>% 
  mutate(trt = factor(trt))

df_p_indiv <- df_all2 %>% filter(param == "p[1]") %>%
  filter(trt %in% c("N8_S1", "N8_S4", "N8_S5", 
                    "N8_S2", "N16_S2", "N8_S3", "N16_S3")) %>% 
  mutate(trt = factor(trt))

df_lambda_adj <- df_all2 %>% 
  filter(param == "lambda[1]") %>% 
  filter(trt %in% c("N8_S1", "N8_S4", "N8_S5", 
                    "N8_S2", "N16_S2", "N8_S3", "N16_S3")) %>% 
  mutate(trt = factor(trt)) %>% filter(cri_width <= 2)

df_lambda_summ <- df_lambda_adj %>% group_by(param, trt) %>% 
  summarize(cap_rate = mean(capture), 
            avg_est = mean(Mean), 
            avg_lb = mean(`2.5%`), 
            avg_ub = mean(`97.5%`), 
            avg_width = mean(cri_width), 
            sd_width = sd(cri_width),
            dg = unique(dg), 
            n_converged = length(Mean), 
            model = unique(model),
            trt = unique(trt), 
            visit = unique(visit), 
            scen = unique(scen))
df_lambda_summ <- df_lambda_summ[,c(1,3:8,10:11,9,2,12:13)]

df_psip_indiv <- rbind(df_psi_sub2, df_p_indiv, df_lambda_adj)
df_psip <- rbind(df_psi_sub, df_p, df_lambda_summ)
df_psip_indiv$param <- factor(df_psip_indiv$param, 
                              levels = c("psi[1]", "p[1]", "lambda[1]"))
levels(df_psip_indiv$param) <- c("psi", "p", "lambda")
df_psip$param <- factor(df_psip$param, 
                        levels = c("psi[1]", "p[1]", "lambda[1]"))
levels(df_psip$param) <- c("psi", "p", "lambda")


# make the figure for rare spp
ggplot(df_psip_indiv) + 
  geom_point(data = df_psip, aes(x = model, 
                                 y = avg_est, 
                                 color = cap_rate), 
             position = position_dodge(width = 0.2)) +
  geom_hline(data = df_psip,
             aes(yintercept = dg)) +
  geom_errorbar(aes(x = jitter(as.numeric(model)),
                    ymin = `2.5%`, ymax = `97.5%`,
                    color = cap_rate),
                lwd = 0.3,
                alpha = 0.3,
                width = 0) +
  geom_errorbar(data = df_psip, 
                aes(x = model, ymin = avg_lb, 
                    ymax = avg_ub, 
                    color = cap_rate), 
                size = 1, width = 0.4, 
                position = position_dodge(width = 0.2)) + 
  facet_grid(param ~ trt, labeller = label_parsed, scales = "free_y") +
  # geom_point(data = df_psi_sub, aes(x = model, y = dg), 
  # col = "magenta", size = 0.75) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_colour_gradient2(
    low = scales::muted("red"),
    mid = "lightgray",
    high = scales::muted("blue"),
    midpoint = 0.5,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour", 
    name = "Capture Rate"
  ) + labs(x = "Model", y = "", title = "Rare Species")


## parameter estimation for common species -- data frame 
# for plotting 
df_psi_sub2 <- df_all2 %>% 
  filter(param == "psi[2]") 

df_psi_sub <- df_all %>% 
  filter(param == "psi[2]")

df_p <- df_all %>% 
  filter(param == "p[2]")

df_p_indiv <- df_all2 %>% 
  filter(param == "p[2]") 

df_lambda_adj <- df_all2 %>% 
  filter(param == "lambda[2]") %>% 
  mutate(trt = factor(trt)) %>% filter(cri_width <= 2)

df_lambda_summ <- df_lambda_adj %>% group_by(param, trt) %>% 
  summarize(cap_rate = mean(capture), 
            avg_est = mean(Mean), 
            avg_lb = mean(`2.5%`), 
            avg_ub = mean(`97.5%`), 
            avg_width = mean(cri_width), 
            sd_width = sd(cri_width),
            dg = unique(dg), 
            n_converged = length(Mean), 
            model = unique(model),
            trt = unique(trt), 
            visit = unique(visit), 
            scen = unique(scen))
df_lambda_summ <- df_lambda_summ[,c(1,3:8,10:11,9,2,12:13)]

df_psip_indiv <- rbind(df_psi_sub2, df_p_indiv, df_lambda_adj)
df_psip <- rbind(df_psi_sub, df_p, df_lambda_summ)
df_psip_indiv$param <- factor(df_psip_indiv$param, 
                              levels = c("psi[2]", "p[2]", "lambda[2]"))
levels(df_psip_indiv$param) <- c("psi", "p", "lambda")
df_psip$param <- factor(df_psip$param, 
                        levels = c("psi[2]", "p[2]", "lambda[2]"))
levels(df_psip$param) <- c("psi", "p", "lambda")

# make the figure for rare spp

ggplot(df_psip_indiv) + 
  geom_point(data = df_psip, aes(x = model, 
                                 y = avg_est, 
                                 color = cap_rate), 
             position = position_dodge(width = 0.2)) +
  geom_hline(data = df_psip,
             aes(yintercept = dg)) +
  geom_errorbar(aes(x = jitter(as.numeric(model)),
                    ymin = `2.5%`, ymax = `97.5%`,
                    color = cap_rate),
                lwd = 0.3,
                alpha = 0.3,
                width = 0) +
  geom_errorbar(data = df_psip, 
                aes(x = model, ymin = avg_lb, 
                    ymax = avg_ub, 
                    color = cap_rate), 
                size = 1, width = 0.4, 
                position = position_dodge(width = 0.2)) + 
  facet_grid(param ~ trt, labeller = label_parsed, scales = "free_y") +
  # geom_point(data = df_psi_sub, aes(x = model, y = dg), 
  # col = "magenta", size = 0.75) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_colour_gradient2(
    low = scales::muted("red"),
    mid = "lightgray",
    high = scales::muted("blue"),
    midpoint = 0.5,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour", 
    name = "Capture Rate"
  ) + labs(x = "Model", y = "", title = "Common Species")



## summary of MCMC convergence
df_converge <- converge_summ(sim_list = list(scen1_v8, scen1_v16, 
                                             scen4_v8, scen4_v16, 
                                             scen5_v8, scen5_v16, 
                                             scen2_v8, scen2_v16, 
                                             scen3_v8, scen3_v16))
psi_converge <- df_converge[order(df_converge$param), ] 
psi_converge

## Setup for MLESite threshold comparison 
# set inputs for using simulation processing functions
# set up scenario and visit lables
scen_name <- c(rep(paste0("S", 1:5), each = 2))
visits <- c(rep(c("8N", "16N"), 5))

# create one big data frame from sim_list summaries
df_all <- data.frame()
for(i in 1:length(sim_list)){
  out <- suppressMessages(
    z_table_MLEsite(sim_results = sim_list[[i]]$MLESite, 
                    NA_combine = TRUE)
    
  )
  out$num_visit <- visits[i]
  out$scenario <- scen_name[i]
  out$trt <- with(out, interaction(num_visit, scenario))
  out$z_decision <- as.numeric(substring(out$z_decide.dg, 1,1))
  df_all <- rbind(df_all, out)
}

# relevel z_dg
df_all$z_dg <- factor(df_all$z_dg)
levels(df_all$z_dg) <- c("True Species Absent", "True Species Present")

# relevel z_decision
df_all$z_decision <- factor(df_all$z_decision)
levels(df_all$z_decision) <- c("Claim Species Absent", "Claim Species Present")

# relevel z_dg and z_decision 
df_all$z_dg <- relevel(df_all$z_dg, 
                       ref = "True Species Present")
df_all$z_decision <- relevel(df_all$z_decision, 
                             ref = "Claim Species Present")
# make alpha a factor
df_all$alpha <- factor(df_all$alpha)

# relevel to compare scenarios in text (left) to those not in text (right)
df_all$trt <- factor(df_all$trt, levels = c("8N.S1", "16N.S1", 
                                            "8N.S4", "16N.S4", 
                                            "8N.S5", "16N.S5", 
                                            "8N.S2", "16N.S2", 
                                            "8N.S3", "16N.S3"))


## df for plotting  
df_plot <- df_all 

# rare spp plot
mle_rare <- df_plot %>% filter(trt %in% c("8N.S1", "16N.S1", "8N.S4", 
                                          "16N.S4", "8N.S5", "16N.S5")) %>%
  filter(species == 1) %>%  
  ggplot(aes(y = cond_prob_z, x = trt, fill = alpha, col = alpha)) +
  geom_boxplot(varwidth = FALSE, notch = TRUE, alpha = 0.5, 
               outlier.colour = "white") + 
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  scale_fill_manual(values = c("#bdbdbd", "#dfc27d", "#bf812d", "#8c510a")) +
  scale_color_manual(values = c("#bdbdbd", "#dfc27d", "#bf812d","#8c510a")) + 
  xlab('Visit.Scenario') + 
  ylab('Calculated decision rate') + 
  ylim(0,1) + 
  facet_grid(z_decision ~ z_dg) +
  ggtitle('MLESite: Rare Species') + 
  geom_vline(xintercept = c(2.5, 4.5)) + 
  theme_bw(base_size = 16)

mle_common <- df_plot %>% filter(trt %in% c("8N.S1", "16N.S1", "8N.S4", 
                                            "16N.S4", "8N.S5", "16N.S5")) %>%
  filter(species == 2) %>%  
  ggplot(aes(y = cond_prob_z, x = trt, fill = alpha, col = alpha)) +
  geom_boxplot(varwidth = FALSE, notch = TRUE, alpha = 0.5, 
               outlier.colour = "white") + 
  geom_point(position = position_jitterdodge(), 
             alpha = 0.5) + 
  scale_fill_manual(values = c("#bdbdbd", "#dfc27d", "#bf812d", "#8c510a")) +
  scale_color_manual(values = c("#bdbdbd", "#dfc27d", "#bf812d", "#8c510a")) +
  xlab('Vist.Scenario') + 
  ylab('Calculated decision rate') + 
  ylim(0,1) + 
  facet_grid(z_decision ~ z_dg) +
  ggtitle('MLE-Site: Common Species') + 
  geom_vline(xintercept = c(2.5, 4.5)) + 
  theme_bw(base_size = 16) 

## ----MLESite-plots Fig S3
gridExtra::grid.arrange(mle_rare, mle_common)


## ----Bayes-decision_setup, 
# set inputs for using post processing functions
z_cutoffs <- c(0.05, 0.25, 0.5, 0.75, 0.95)

# one data frame for summarizing results
df_bayes <- data.frame()

for(i in 1:length(sim_list)){
  ctdet <- data.frame()
  for(j in 1:5){
    df_temp <- suppressMessages(
      z_table_nimble(sim_results = sim_list[[i]]$CtDetection, 
                     z_cutoff = z_cutoffs[j])
    )
    df_temp$z_cutoff <- z_cutoffs[j]
    ctdet <- rbind(ctdet, df_temp)
  }
  ctdet$num_visit <- visits[i]
  ctdet$scenario <- scen_name[i]
  ctdet$trt <- with(ctdet, interaction(num_visit, scenario))
  rem <- data.frame()
  for(j in 1:5){
    temp <- suppressMessages(
      z_table_nimble(sim_results = sim_list[[i]]$Remove, 
                     z_cutoff = z_cutoffs[j])
    )
    temp$z_cutoff <- z_cutoffs[j]
    rem <- rbind(rem, temp)
  }
  rem$num_visit <- visits[i]
  rem$scenario <- scen_name[i]
  rem$trt <- with(rem, interaction(num_visit, scenario))
  naive <- data.frame()
  for(j in 1:5){
    temp <- suppressMessages(
      z_table_nimble(sim_results = sim_list[[i]]$Naive, 
                     z_cutoff = z_cutoffs[j])
    )
    temp$z_cutoff <- z_cutoffs[j]
    naive <- rbind(naive, temp)
  }
  naive$num_visit <- visits[i]
  naive$scenario <- scen_name[i]
  naive$trt <- with(naive, interaction(num_visit, scenario))
  
  df_bayes <- rbind(df_bayes, ctdet, rem, naive)
}

df_bayes$z_dg <- factor(df_bayes$z_dg, 
                        labels = c("True Species Absent",
                                   "True Species Present"))

df_bayes$z_decision <- factor(df_bayes$z_decision, 
                              labels = c("Claim Species Absent", 
                                         "Claim Species Present"))

df_bayes$z_dg <- relevel(df_bayes$z_dg, 
                         ref = "True Species Present")
df_bayes$z_decision <- relevel(df_bayes$z_decision, 
                               ref = "Claim Species Present")

df_bayes$decision_rule <- paste0("z_cutoff = ", df_bayes$z_cutoff)  
df_bayes$z_cutoff <- factor(paste0("z_cutoff = ", df_bayes$z_cutoff))


## ----threshold-plots Figs S4/S5
# rare 2sppCt
rare2spp <- df_bayes %>% filter(trt %in% c("8N.S1", "16N.S1", "8N.S4", 
                                           "16N.S4", "8N.S5", "16N.S5")) %>%
  filter(species == 1) %>% filter(model == "2SppCt") %>%  
  ggplot(aes(y = cond_prob_z, x = trt, fill = z_cutoff, col = z_cutoff)) +
  geom_boxplot(varwidth = FALSE, notch = TRUE, alpha = 0.5, 
               outlier.colour = "white") +
  geom_point(position = position_jitterdodge(), 
             alpha = 0.5) + 
  scale_fill_manual(values = c("#a1d99b","#74c476","#41ab5d", 
                               "#238b45", "#00441b")) + 
  scale_color_manual(values = c("#a1d99b","#74c476","#41ab5d", 
                                "#238b45", "#00441b")) +
  xlab('Vist.Scenario') + 
  ylab('Calculated decision rate') + 
  ylim(0,1) + 
  facet_grid(z_decision ~ z_dg) +
  ggtitle('2SppCt: Rare Species') + 
  geom_vline(xintercept = c(2.5, 4.5)) + 
  theme_bw(base_size = 16) 

# rare remove 
rare_rem <- df_bayes %>% filter(trt %in% c("8N.S1", "16N.S1", "8N.S4", 
                                           "16N.S4", "8N.S5", "16N.S5")) %>%
  filter(species == 1) %>% filter(model == "Remove") %>%  
  ggplot(aes(y = cond_prob_z, x = trt, fill = z_cutoff, col = z_cutoff)) +
  geom_boxplot(varwidth = FALSE, notch = TRUE, alpha = 0.5, 
               outlier.colour = "white") +
  geom_point(position = position_jitterdodge(), 
             alpha = 0.5) + 
  scale_fill_manual(values = c("#bdbdbd", "#c7eae5", "#80cdc1", 
                               "#35978f", "#01665e")) +
  scale_color_manual(values = c("#bdbdbd", "#c7eae5", "#80cdc1", 
                                "#35978f", "#01665e")) +
  xlab('Vist.Scenario') + 
  ylab('Calculated decision rate') + 
  ylim(0,1) + 
  facet_grid(z_decision ~ z_dg) +
  ggtitle('Remove: Rare Species') + 
  geom_vline(xintercept = c(2.5, 4.5)) + 
  theme_bw(base_size = 16) 

# common 2sppct
common_2spp <- df_bayes %>% filter(trt %in% c("8N.S1", "16N.S1", "8N.S4", 
                                              "16N.S4", "8N.S5", "16N.S5")) %>%
  filter(species == 2) %>% filter(model == "2SppCt") %>%  
  ggplot(aes(y = cond_prob_z, x = trt, fill = z_cutoff, col = z_cutoff)) +
  geom_boxplot(varwidth = FALSE, notch = TRUE, 
               alpha = 0.5, 
               outlier.colour = "white") +
  geom_point(position = position_jitterdodge(), 
             alpha = 0.5) + 
  scale_fill_manual(values = c("#a1d99b","#74c476","#41ab5d", 
                               "#238b45", "#00441b")) +
  scale_color_manual(values = c("#a1d99b","#74c476","#41ab5d", 
                                "#238b45", "#00441b")) +
  xlab('Vist.Scenario') + 
  ylab('Calculated decision rate') + 
  ylim(0,1) + 
  facet_grid(z_decision ~ z_dg) +
  ggtitle('2SppCt: Common Species') + 
  geom_vline(xintercept = c(2.5, 4.5)) + 
  theme_bw(base_size = 16) 

# common remove 
common_rem <- df_bayes %>% filter(trt %in% c("8N.S1", "16N.S1", "8N.S4", 
                                             "16N.S4", "8N.S5", "16N.S5")) %>%
  filter(species == 2) %>%  filter(model == "Remove") %>% 
  ggplot(aes(y = cond_prob_z, x = trt, fill = z_cutoff, col = z_cutoff)) +
  geom_boxplot(varwidth = FALSE, notch = TRUE, 
               alpha = 0.5, outlier.colour = "white") +
  geom_point(position = position_jitterdodge(), 
             alpha = 0.5) + 
  scale_fill_manual(values = c("#bdbdbd", "#c7eae5", "#80cdc1", 
                               "#35978f", "#01665e")) +
  scale_color_manual(values = c("#bdbdbd", "#c7eae5", "#80cdc1", 
                                "#35978f", "#01665e")) +
  xlab('Vist.Scenario') + 
  ylab('Calculated decision rate') + 
  ylim(0,1) + 
  facet_grid(z_decision ~ z_dg) +
  ggtitle('Remove: Common Species') + 
  geom_vline(xintercept = c(2.5, 4.5)) + 
  theme_bw(base_size = 16) 


## ----Remove-threshold-rare Fig S4
gridExtra::grid.arrange(rare2spp, rare_rem)


## ----Remove-threshold-common Fig S5
gridExtra::grid.arrange(common_2spp, common_rem)


# recreate Tables S3 and S4
# with summary from "optimal thresholds" presented in text
# set inputs for using post processing functions
# set up scenario and visit labels
scen_name <- c(rep(paste0("S", 1:5), each = 2))
visits <- c(rep(c("8N", "16N"), 5))

# MLESite results df 
# create one big data frame from sim_list summaries
# For MLESite results
df_all <- data.frame()
for(i in 1:length(sim_list)){
  out <- suppressMessages(
    z_table_MLEsite(sim_results = sim_list[[i]]$MLESite, 
                    NA_combine = TRUE)
    
  )
  out$num_visit <- visits[i]
  out$scenario <- scen_name[i]
  out$trt <- with(out, interaction(num_visit, scenario))
  out$z_decision <- as.numeric(substring(out$z_decide.dg, 1,1))
  df_all <- rbind(df_all, out)
}

# relevel z_dg
df_all$z_dg <- factor(df_all$z_dg)
levels(df_all$z_dg) <- c("True Species Absent", "True Species Present")

# relevel z_decision
df_all$z_decision <- factor(df_all$z_decision)
levels(df_all$z_decision) <- c("Claim Species Absent", "Claim Species Present")
df_all$model <- "MLESite"

# 2Sppct and Remove results
# get one df to compare to MLESite results

# set z_cutoffs 
z_cutoffs <- c(0.05, 0.25, 0.75, 0.95)

# initialize df
df_bayes <- data.frame()

for(i in 1:length(sim_list)){
  ctdet <- data.frame()
  for(j in 1:4){
    df_temp <- suppressMessages(
      z_table_nimble(sim_results = sim_list[[i]]$CtDetection, 
                     z_cutoff = z_cutoffs[j])
    )
    df_temp$z_cutoff <- z_cutoffs[j]
    ctdet <- rbind(ctdet, df_temp)
  }
  ctdet$num_visit <- visits[i]
  ctdet$scenario <- scen_name[i]
  ctdet$trt <- with(ctdet, interaction(num_visit, scenario))
  rem <- data.frame()
  for(j in 1:4){
    temp <- suppressMessages(
      z_table_nimble(sim_results = sim_list[[i]]$Remove, 
                     z_cutoff = z_cutoffs[j])
    )
    temp$z_cutoff <- z_cutoffs[j]
    rem <- rbind(rem, temp)
  }
  rem$num_visit <- visits[i]
  rem$scenario <- scen_name[i]
  rem$trt <- with(rem, interaction(num_visit, scenario))
  df_bayes <- rbind(df_bayes, ctdet, rem)
}

df_bayes$z_dg <- factor(df_bayes$z_dg, 
                        labels = c("True Species Absent",
                                   "True Species Present"))

df_bayes$z_decision <- factor(df_bayes$z_decision, 
                              labels = c("Claim Species Absent", 
                                         "Claim Species Present"))
df_all$decision_rule <- paste0("alpha = ", df_all$alpha)
df_bayes$decision_rule <- paste0("z_cutoff = ", df_bayes$z_cutoff)  

levels(df_bayes$z_dg) <- c("True Species Absent", "True Species Present")

df_MLE <- df_all %>% dplyr::select(z_decide.dg, species, 
                                   cond_prob_z, z_dg, scenario, z_decision, decision_rule, 
                                   trt, model)
df_Ctdet <- df_bayes %>% filter(model == "2SppCt")
df_2SppCt <- df_Ctdet %>% dplyr::select(z_decide.dg, species, 
                                        cond_prob_z, z_dg, scenario, z_decision, decision_rule, 
                                        trt, model)
df_MLE2spp <- full_join(df_MLE, df_2SppCt)
df_MLE2spp$decision_rule <- factor(df_MLE2spp$decision_rule, 
                                   levels = c("alpha = 0.05", 
                                              "alpha = 0.1", 
                                              "alpha = 0.15", 
                                              "alpha = 0.2", 
                                              "z_cutoff = 0.95", 
                                              "z_cutoff = 0.75", 
                                              "z_cutoff = 0.25", 
                                              "z_cutoff = 0.05")) 

df_MLE2spp$z_dg <- relevel(df_MLE2spp$z_dg, ref = "True Species Present")
df_MLE2spp$z_decision <- relevel(df_MLE2spp$z_decision, ref = "Claim Species Present")

# MLEsite alpha = 0.1 S1,S4,S5, with 16N
df_MLEfinal <- df_MLE %>% filter(trt %in% c("16N.S1", "16N.S4", "16N.S5"),
                                 decision_rule == "alpha = 0.1") 
head(df_MLEfinal)
df_MLEfinal$model <- "MLESite"
df_bayesfinal <- df_bayes %>% filter(trt %in% 
                                       c("16N.S1", "16N.S4", "16N.S5"),
                                     decision_rule %in% c("z_cutoff = 0.25",
                                                          "z_cutoff = 0.75")
) %>% dplyr::select(z_decide.dg,
                    species,
                    cond_prob_z, 
                    z_dg, scenario,
                    z_decision, 
                    decision_rule, trt,
                    model)

df_compare <- full_join(df_MLEfinal, df_bayesfinal)
df_compare$decision_rule <- factor(df_compare$decision_rule, 
                                   levels = c("alpha = 0.1", 
                                              "z_cutoff = 0.75", 
                                              "z_cutoff = 0.25")) 

df_compare$z_dg <- relevel(df_compare$z_dg, 
                           ref = "True Species Present")
df_compare$z_decision <- relevel(df_compare$z_decision, 
                                 ref = "Claim Species Present")
df_compare$xlab <- factor(with(df_compare, interaction(trt, model)), 
                          levels = c("16N.S1.MLESite", "16N.S1.Remove", "16N.S1.2SppCt",
                                     "16N.S4.MLESite", "16N.S4.Remove","16N.S4.2SppCt", 
                                     "16N.S5.MLESite", "16N.S5.Remove","16N.S5.2SppCt"))

df_compare$xlab <- factor(df_compare$xlab, 
                          labels = c("S1.MLESite", "S1.Remove", "S1.2SppCt",
                                     "S4.MLESite", "S4.Remove","S4.2SppCt",
                                     "S5.MLESite", "S5.Remove","S5.2SppCt"))

# spp1 
df_plot <- df_compare %>% 
  filter(model == "MLESite" | 
           (z_dg == "True Species Absent" & 
              decision_rule == "z_cutoff = 0.75") | 
           (z_dg == "True Species Present" & 
              decision_rule == "z_cutoff = 0.25")) 


## ----table-results-------------------------------------------------------
results_summ <- df_plot %>%
  group_by(model, z_decide.dg, species, trt) %>%
  summarise(median = median(cond_prob_z),
            iqr = IQR(cond_prob_z), q1 = quantile(cond_prob_z, 0.25), 
            q3 = quantile(cond_prob_z, 0.75)) %>%
  filter(species == 1)

results_summ <- results_summ[order(results_summ$z_decide.dg, 
                                   results_summ$model), ]

names(results_summ) <- c("Approach", "Decision.Truth", "Species", 
                         "Visit.Scenario", "Median", "IQR", "Q1", "Q3")

results_summ[36:1, ]
