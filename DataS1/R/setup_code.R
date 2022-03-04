####Data S1 for Methods in Ecology and Evolution manuscript
##Modeling misclassification in multi-species acoustic data when
##estimating occupancy and relative activity
##by Wilson J Wright, Kathryn M Irvine, Emily S Almberg, and Andrea R Litt.
##This code reproduces the simulations described in the manuscript.

##Author: Wilson Wright
##R Version 3.5.1

## Edited by KM Irvine to strip it down for MLE vs. Bayes exploration 

##Description:
#Setup code for running a small simulation study.
#Data from two species are simulated and various false positive models
#are compared to explore how inferences could under this scenario.

#Function to simulate data, returns list formatted for Stan models.
#This is simulating two species activity with count detections and
#include misidentifications of individual detections.
sim_data <- function(n.sites1, n.sites2, n.visits,
                     psi.vec, lambda.vec, theta.mat){
  #n.sites1 = number of sites with confirmed visits denoted as 'a' sites
  #(assume first two visits for each of these sites is confirmed)
  #n.sites2 = number of sites without confirmed visits denoted as 'b' sites
  #n.visits = number of visits per site, assumed constant
  #psi.vec = vector of occupancy probabilities for species 1 and 2
  #lambda.vec = vector of call recording rates for species 1 and 2
  #theta.mat = confusion matrix for automatic call classifications
  #rows of theta.mat should sum to one.
  #diagonal elements are correct identification probabilities
  #first row corresponds to probabilities for calls truly from species 1

  stopifnot(n.sites1 > 0 && n.sites1 > 0)
  stopifnot(n.visits > 1 )

  n.total <- n.sites1 + n.sites2

  z1 <- rbinom(n.total, 1, psi.vec[1])
  z2 <- rbinom(n.total, 1, psi.vec[2])

  n.obs <- n.total * n.visits
  z1.rep <- rep(z1, each = n.visits)
  z2.rep <- rep(z2, each = n.visits)

  y1 <- z1.rep * rpois(n.obs, lambda.vec[1])
  y2 <- z2.rep * rpois(n.obs, lambda.vec[2])

  c1 <- c2 <- matrix(NA, ncol = 2, nrow = n.obs)
  c1.unconf <- c2.unconf <- numeric(n.obs)

  for(n in 1:n.obs){
    c1[n, ] <- rmultinom(1, y1[n], theta.mat[1, ])
    c2[n, ] <- rmultinom(1, y2[n], theta.mat[2, ])
    c1.unconf[n] <- c1[n, 1] + c2[n, 1]
    c2.unconf[n] <- c1[n, 2] + c2[n, 2]
  }

  #we just need c1 and c2  col sums for autoiIDs

  site_ID <- rep(c(1:n.total), each = n.visits)
  visit_ID <- rep(1:n.visits, n.total)



  MLE_data <-list(Z1=z1, Z2=z2,
                          sites_ID=site_ID, visit_ID = visit_ID,
                          Spp1_NumCalls=y1,Spp2_NumCalls=y2,
                          Spp1_AutoIDs=c1.unconf,Spp2_AutoIDs=c2.unconf)


  n.obs1 <- n.sites1 * n.visits  #confirmed obs
  n.obs2 <- n.sites2 * n.visits #unconfirmed obs
  n.conf <- 2
  n.unconf <- n.visits - n.conf
  keep.conf <- rep(c(rep(TRUE, n.conf), rep(FALSE, n.unconf)), n.sites1)
  c1.conf <- c1[1:n.obs1, ][keep.conf, ]
  c2.conf <- c2[1:n.obs1, ][keep.conf, ]
  c1.unconf <- c1[1:n.obs1, ][!keep.conf, ]
  c2.unconf <- c2[1:n.obs1, ][!keep.conf, ]
  c1.sum.a <- c1.unconf[, 1] + c2.unconf[, 1]
  c2.sum.a <- c1.unconf[, 2] + c2.unconf[, 2]

  c1.sum.b <- c1[(n.obs1 + 1):n.obs, 1] + c2[(n.obs1 + 1): n.obs, 1]
  c2.sum.b <- c1[(n.obs1 + 1):n.obs, 2] + c2[(n.obs1 + 1): n.obs, 2]

  sites.conf <- rep(c(1:n.sites1), each = n.visits)[keep.conf]
  sites.unconf.a <- rep(c(1:n.sites1), each = n.visits)[!keep.conf]
  sites.unconf.b <- rep(c((n.sites1 + 1):n.total), each = n.visits)


  Stan_model_data<-list(c1.conf = c1.conf, c2.conf = c2.conf,
                                    c1.sum.a = c1.sum.a, c2.sum.a = c2.sum.a,
                                    c1.sum.b = c1.sum.b, c2.sum.b = c2.sum.b,
                                    sites.conf = sites.conf, sites.unconf.a = sites.unconf.a,
                                    sites.unconf.b = sites.unconf.b)

  return(list("MLE_Britzke_data"= MLE_data, "Stan_model_data"=Stan_model_data))


  ##c1.conf is how calls confirmed to species 1 were classified for confirmed visits
  ##c2.conf is same thing for calls from species 2
  ##c1.sum.a is the unconfirmed calls classified TO species 1 from "a" sites
  ##c2.sum.a is same thing for calls classified as species 2
  ##c1.sum.b is counts of unconfirmed calls to species 1 from "b" sites
  ##c2.sum.b same thing for species 2
  ##all site variables just keep track of which sites correspond to each set of observations
}

#Function to put the needed data in a list for rstan for AmbigOnlyCts_2pp; count detection with only ambiguous cts

input_data1 <- function(sim_obj, n.visits){
  out <- list('n_sites_b' = length(unique(sim_obj$sites.unconf.b)),
              'n_visits_b' = n.visits,
              'n_unconf_b' = length(sim_obj$sites.unconf.b),
              'unconf1b' = sim_obj$c1.sum.b,
              'unconf2b' = sim_obj$c2.sum.b)
  return(out)
}


#function to store if there were any divergent transitions and/or
#max tree depth warnings for each model fit
save_warnings <- function(fit_obj){
  temp <- do.call(rbind, get_sampler_params(fit_obj, inc_warmup = FALSE))
  div <- mean(temp[, 'divergent__'] > 0)
  tree <- mean(temp[, 'treedepth__'] > 9)
  return(list(div = div, tree = tree))
}

#function to store the maximum rhat and minimum effective sample size.
save_convergence <- function(fit_sum){
  max.rhat <- max(fit_sum[, 'Rhat'])
  min.neff <- min(fit_sum[, 'n_eff'])
  return(list(rhat = max.rhat, neff = min.neff))
}

#Simulation wrapper function
#This simulates a dataset assuming the 2-spp count detection model.
#Then analyzes it with each of the four different models. Do this many times
#and save all of the model summaries.
#User can specify the number of simulated datasets and the parameters used
#to simulate the datasets.
sim_wrapper <- function(n.sim, n.sites1, n.sites2, n.visits,
                        psi1, lambda1, theta1,
                        psi2, lambda2, theta2){
  psi.vec <- c(psi1, psi2)
  lambda.vec <- c(lambda1, lambda2)
  theta.mat <- matrix(c(theta1, theta2), byrow = T, ncol =2)

  store.fit1 <- array(NA, c(8, 10, n.sim))


  dimnames(store.fit1)[[1]] <- c('psi1', 'psi2', 'lambda1', 'lambda2',
                                 'theta11', 'theta12', 'theta21', 'theta22')
  dimnames(store.fit1)[[2]] <- c('mean', 'se', 'sd', 'q2.5', 'q25', 'q50', 'q75', 'q97.5', 'n', 'r')
  dimnames(store.fit1)[[3]] <- c(1:n.sim)

  fit.div <- fit.tree <- fit.rhat <- fit.neff <- numeric(n.sim)

  for(n in 1:n.sim){
    temp.data <- sim_data(n.sites1, n.sites2, n.visits,
                          psi.vec, lambda.vec, theta.mat)
    model1.data <- input_data1(temp.data, n.visits)

    model1.fit <- sampling(model1, model1.data, iter = 800)


    model1.sum <- summary(model1.fit,
                          pars = c('psi1', 'psi2', 'lambda1', 'lambda2', 'theta1', 'theta2'),
                          use_cache = F)$summary

    model1.warn <- save_warnings(model1.fit)

    model1.conv <- save_convergence(model1.sum)

    store.fit1[,, n] <- model1.sum

    fit.div[n] <- model1.warn$div

    fit.tree[n] <- model1.warn$tree

    fit.rhat[n] <- model1.conv$rhat

    fit.neff[n] <- model1.conv$neff

  }

  return(list(store.fit1 = store.fit1,
              fit.div = fit.div,
              fit.tree = fit.tree,
              fit.rhat = fit.rhat,
              fit.neff = fit.neff))
}

##Function to help summarize results saved in the sim.wrapper
sim_summary <- function(sim_obj, psi1.true){
  keep1a <- sim_obj$fit.div == 0
  keep1b <- sim_obj$fit.tree == 0
  keep1c <- sim_obj$fit.rhat < 1.1
  keep1d <- sim_obj$fit.neff> 400


  keep1 <- apply(cbind(keep1a, keep1b, keep1c, keep1d),
                 1, all)

  psi1.means <- mean(sim_obj$store.fit1['psi1', 'mean', keep1])
  psi1.lwbnd <-mean(sim_obj$store.fit1['psi1', 'q2.5', keep1])
  psi1.upbnd <- mean(sim_obj$store.fit1['psi1', 'q97.5', keep1])
  psi1.cover <- mean(sim_obj$store.fit1['psi1', 'q2.5', keep1] < psi1.true &
                         sim_obj$store.fit1['psi1', 'q97.5', keep1] > psi1.true)

  out <- cbind(psi1.lwbnd, psi1.means, psi1.upbnd, psi1.cover)
  return(out)
}
