# mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# comma between r and simoutputFcn was giving us grief! 
# model code - count detection with known theta
cd_nimble_ThetaKnown <- nimbleCode({
  
  # Assuming theta is known, read in as data
  #--------------------------------------------------------------
  # priors: independent on lambda and psi
  #--------------------------------------------------------------
  #changing prior for Simulation 5 results only to see if helps?? from a shape =2 to shape =1
  for(k in 1:n.spp){
    lambda[k] ~ dgamma(shape = 2,
                       rate = 0.25)
    psi[k] ~ dbeta(1,1)
  }
  
  for(i in 1:n.sites){
    Z[i,1]~ dbinom(psi[1],1)
    Z[i,2]~ dbinom(psi[2],1)
    for(j in 1:n.revisits){
      Spp1_AutoIDs[i,j] ~ dpois(Z[i,1]*lambda[1]*theta[1,1] +
                                  Z[i,2]*lambda[2]*theta[2,1])
      Spp2_AutoIDs[i,j] ~ dpois(Z[i,1]*lambda[1]*theta[1,2] +
                                  Z[i,2]*lambda[2]*theta[2,2])
    }
  }
  
})# end model

# model code - independent vanilla occupancy models for 2 spp
cd_nimble_OccupancyModel <- nimbleCode({
  #priors
  for(k in 1:n.spp){
    psi[k] ~ dunif(0,1)
    p[k] ~ dunif(0,1)
  }
  for(i in 1:n.sites){
    Z[i,1]~ dbin(psi[1],1)
    Z[i,2]~ dbin(psi[2],1)
    
    Spp1[i] ~ dbin(p[1]*Z[i,1], size=n.revisits)
    Spp2[i] ~ dbin(p[2]*Z[i,2], size=n.revisits)
  }
  
})# end model

# thetas = matrix(c(.8,.2,.35,.65),nrow=2,=2,byrow=TRUE) 
# n.sites1 = 5
# n.sites2 = 50
# n.visits = 3
# psi.vec = c(.75,.25)
# lambda.vec = c(1,2)

# simulation function with defaults for scenario 1 from overleaf 
# doc
psi_cov_function <- function(n_iter, 
                             thetas = matrix(c(.9,.1,.35,.65),
                                             nrow=2,
                                             ncol=2,byrow=TRUE), 
                             n.sites1 = 5,
                             n.sites2 = 50,
                             n.visits = 16,
                             psi.vec = c(.25,.75),
                             lambda.vec = c(0.3,10), 
                             alpha_vec = c(0.05, 0.1, 0.15, 0.2), 
                             nburn = 2500, niter = 5000, thin = 5){
  
  summ_out_CtDet<-summ_out_remove<-summ_out_naive<-summ_out_mlesite<-as.list(1:n_iter)
  
  for(i in 1:n_iter){
    sim_dat <- sim_data(n.sites1 = n.sites1,
                        n.sites2 = n.sites2,
                        n.visits = n.visits,
                        psi.vec = psi.vec,
                        lambda.vec = lambda.vec, 
                        theta.mat = thetas)
    
    
    # fit nimble ct-detection model 
    n.sites.tot<-length(unique(sim_dat$MLE_Britzke_data$sites_ID))
    
    spp1 <- matrix(sim_dat$MLE_Britzke_data$Spp1_AutoIDs,
                   nrow =  n.sites.tot, ncol = n.visits, byrow=TRUE)
    spp2 <- matrix(sim_dat$MLE_Britzke_data$Spp2_AutoIDs,
                   nrow =  n.sites.tot, ncol = n.visits, byrow=TRUE)
    cd_data <- list(Spp1_AutoIDs = spp1,
                    Spp2_AutoIDs = spp2,
                    theta = thetas)
    
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
                                 thin = thin,
                                 niter = niter,
                                 nburn = nburn, 
                                 nchains = 3,
                                 samplesAsCodaMCMC = TRUE)
    
    # summary dataframe
    tmp <- summary(cd_mcmc_simple)
    model1_sum <- as_tibble(cbind(tmp$statistics, tmp$quantiles)) %>%
      mutate(
        param = rownames(tmp$statistics),
        model = "2SppCt",
        sim_iter = i
      ) %>%
      dplyr::select(param, model, sim_iter, everything()) %>%
      mutate(
        truth = c(sim_dat$MLE_Britzke_data$Z1, sim_dat$MLE_Britzke_data$Z2, lambda.vec, psi.vec),
        n_eff = effectiveSize(cd_mcmc_simple)
      )
    
    # add Rhat
    mat <- matrix(NA, nrow = nrow(model1_sum), ncol = )
    colnames(mat) <- c("Rhat")
    for(a in 1:nrow(model1_sum)){
      tmp <- sapply(cd_mcmc_simple, function(x) x[,a])
      mat[a,] <- c(Rhat(tmp))
    }
    model1_sum$Rhat <- c(mat)
    
    # capture
    model1_sum <- model1_sum %>%
      mutate(
        capture = ifelse(`2.5%` <= truth & truth <= `97.5%`, 1, 0),
        cri_width = `97.5%` - `2.5%`,
        mode = apply(do.call("rbind", cd_mcmc_simple), 2, function(x) getmode(x)),
        
        # NAs for latent Zs
        capture = ifelse(grepl("Z", param), NA, capture),
        cri_width = ifelse(grepl("Z", param), NA, cri_width),
        
        #NAs for other params
        mode = ifelse(!grepl("Z", param), NA, mode),
        
        # add species
        species = case_when(
          grepl("Z[[]", param) ~ strsplit(param, split = ", ") %>% 
            sapply(., function(x) x[2]) %>%
            strsplit(., split = "[]]") %>%
            unlist, 
          TRUE ~ strsplit(param, split = "[[]") %>%
            sapply(., function(x) x[2]) %>%
            strsplit(., split = "[]]") %>%
            unlist
        ),
        param_clean = strsplit(param, split = "[[]") %>% sapply(., function(x) x[1]),
        site = case_when(
          grepl("Z[[]", param) ~ strsplit(param, split = "[[]") %>%
            sapply(., function(x) x[2]) %>%
            strsplit(., split = ",") %>%
            sapply(., function(x) x[1])
        )
      ) %>%
      dplyr::select(model, sim_iter, param, param_clean, species, site, everything())
    
    summ_out_CtDet[[i]] <- model1_sum
    
    # Fit Naive Occupancy model + MLE remove with Occupancy
    # create df for britzke MLE functions just pulling off a
    test_df <- data.frame(do.call(cbind, 
                                  sim_dat$MLE_Britzke_data[c(3:4, 7:8)]))
    z_spp1 <- sim_dat$MLE_Britzke_data$Z1
    z_spp2 <- sim_dat$MLE_Britzke_data$Z2
    ambig_det <- test_df[,3:4]
    
    phi_mat <- t(thetas) # needed to define phi_mat
    mle_spp1 <- mle_closedform_2spp(ambig_det = ambig_det, 
                                    phi_mat = phi_mat,
                                    z_vec = z_spp1, 
                                    n_visit = n.visits, 
                                    spp_idx = 1, 
                                    adjusted = TRUE)
    # mle_spp1
    # implement for species 2 
    mle_spp2 <- mle_closedform_2spp(ambig_det = ambig_det, 
                                    phi_mat = phi_mat,
                                    z_vec = z_spp2, 
                                    n_visit = n.visits, 
                                    spp_idx = 2, 
                                    adjusted = TRUE)
    
    # create detection history matrix with Y_noFP MLE results
    dh1 <- mle_dh(mle_output = mle_spp1,
                  n_visit = n.visits, spp_idx = 1)
    dh2 <- mle_dh(mle_output = mle_spp2,
                  n_visit = n.visits, spp_idx = 2)
    dat_remove<-data.frame(Spp1 = rowSums(dh1$dh_remove), Spp2 = rowSums(dh2$dh_remove))
    dat_naive<-data.frame(Spp1 = rowSums(dh1$dh_naive), Spp2 = rowSums(dh2$dh_naive))
    cd_constants <- list(n.sites =  n.sites.tot, n.spp = 2, n.revisits = n.visits)
    n.sites.tot<-length(unique(sim_dat$MLE_Britzke_data$sites_ID))
    
    # initial values
    cd_inits <- function(spp1,spp2) {
      naive_z1 <- as.numeric(spp1 > 0)
      naive_z2 <- as.numeric(spp2 > 0)
      Z_init <- cbind(naive_z1, naive_z2)
      p_init <- runif(2)
      psi_init <- runif(2)
      list(Z = Z_init,
           p = p_init,
           psi = psi_init)
    }
    
    
    # fit the model in one-step
    cd_mcmc_naive <- nimbleMCMC(code = cd_nimble_OccupancyModel,
                                data = list(Spp1=dat_naive$Spp1,Spp2=dat_naive$Spp2),
                                constants = cd_constants,
                                inits = cd_inits(spp1=dat_naive$Spp1,spp2=dat_naive$Spp2),
                                monitors = c("psi", "p", "Z"),
                                thin = thin,
                                niter = niter,
                                nburn = nburn, nchains = 3,
                                samplesAsCodaMCMC = TRUE)
    
    # summary dataframe
    tmp <- summary(cd_mcmc_naive)
    model2_sum <- as_tibble(cbind(tmp$statistics, tmp$quantiles)) %>%
      mutate(
        param = rownames(tmp$statistics),
        model = "Naive",
        sim_iter = i
      ) %>%
      dplyr::select(param, model, sim_iter, everything()) %>%
      mutate(
        truth = c(sim_dat$MLE_Britzke_data$Z1, sim_dat$MLE_Britzke_data$Z2, 1-exp(-lambda.vec), psi.vec),
        n_eff = effectiveSize(cd_mcmc_naive)
      )
    
    # add Rhat
    mat <- matrix(NA, nrow = nrow(model2_sum), ncol = )
    colnames(mat) <- c("Rhat")
    for(a in 1:nrow(model2_sum)){
      tmp <- sapply(cd_mcmc_naive, function(x) x[,a])
      mat[a,] <- c(Rhat(tmp))
    }
    model2_sum$Rhat <- c(mat)
    
    # capture
    model2_sum <- model2_sum %>%
      mutate(
        capture = ifelse(`2.5%` <= truth & truth <= `97.5%`, 1, 0),
        cri_width = `97.5%` - `2.5%`,
        mode = apply(do.call("rbind", cd_mcmc_naive), 2, function(x) getmode(x)),
        
        # NAs for latent Zs
        capture = ifelse(grepl("Z", param), NA, capture),
        cri_width = ifelse(grepl("Z", param), NA, cri_width),
        
        #NAs for other params
        mode = ifelse(!grepl("Z", param), NA, mode),
        
        # add species
        species = case_when(
          grepl("Z[[]", param) ~ strsplit(param, split = ", ") %>% 
            sapply(., function(x) x[2]) %>%
            strsplit(., split = "[]]") %>%
            unlist, 
          TRUE ~ strsplit(param, split = "[[]") %>%
            sapply(., function(x) x[2]) %>%
            strsplit(., split = "[]]") %>%
            unlist
        ),
        param_clean = strsplit(param, split = "[[]") %>% sapply(., function(x) x[1]),
        site = case_when(
          grepl("Z[[]", param) ~ strsplit(param, split = "[[]") %>%
            sapply(., function(x) x[2]) %>%
            strsplit(., split = ",") %>%
            sapply(., function(x) x[1])
        )
      ) %>%
      dplyr::select(model, sim_iter, param, param_clean, species, site, everything())
    summ_out_naive[[i]] <- model2_sum
    
    
    cd_mcmc_remove <- nimbleMCMC(code = cd_nimble_OccupancyModel,
                                 data = list(Spp1=dat_remove$Spp1,Spp2=dat_remove$Spp2),
                                 constants = cd_constants,
                                 inits = cd_inits(spp1=dat_remove$Spp1,spp2=dat_remove$Spp2),
                                 monitors = c("psi","p", "Z"),
                                 thin = thin,
                                 niter = niter,
                                 nburn = nburn, nchains = 3,
                                 samplesAsCodaMCMC = TRUE)
    
    
    # summary dataframe
    tmp <- summary(cd_mcmc_remove)
    model3_sum <- as_tibble(cbind(tmp$statistics, tmp$quantiles)) %>%
      mutate(
        param = rownames(tmp$statistics),
        model = "Remove",
        sim_iter = i
      ) %>%
      dplyr::select(param, model, sim_iter, everything()) %>%
      mutate(
        truth = c(sim_dat$MLE_Britzke_data$Z1, sim_dat$MLE_Britzke_data$Z2, 1-exp(-lambda.vec), psi.vec),
        n_eff = effectiveSize(cd_mcmc_remove)
      )
    
    # add Rhat
    mat <- matrix(NA, nrow = nrow(model3_sum), ncol = )
    colnames(mat) <- c("Rhat")
    for(a in 1:nrow(model3_sum)){
      tmp <- sapply(cd_mcmc_remove, function(x) x[,a])
      mat[a,] <- c(Rhat(tmp))
    }
    model3_sum$Rhat <- c(mat)
    
    # capture
    model3_sum <- model3_sum %>%
      mutate(
        capture = ifelse(`2.5%` <= truth & truth <= `97.5%`, 1, 0),
        cri_width = `97.5%` - `2.5%`,
        mode = apply(do.call("rbind", cd_mcmc_remove), 2, function(x) getmode(x)),
        
        # NAs for latent Zs
        capture = ifelse(grepl("Z", param), NA, capture),
        cri_width = ifelse(grepl("Z", param), NA, cri_width),
        
        #NAs for other params
        mode = ifelse(!grepl("Z", param), NA, mode),
        
        # add species
        species = case_when(
          grepl("Z[[]", param) ~ strsplit(param, split = ", ") %>% 
            sapply(., function(x) x[2]) %>%
            strsplit(., split = "[]]") %>%
            unlist, 
          TRUE ~ strsplit(param, split = "[[]") %>%
            sapply(., function(x) x[2]) %>%
            strsplit(., split = "[]]") %>%
            unlist
        ),
        param_clean = strsplit(param, split = "[[]") %>% sapply(., function(x) x[1]),
        site = case_when(
          grepl("Z[[]", param) ~ strsplit(param, split = "[[]") %>%
            sapply(., function(x) x[2]) %>%
            strsplit(., split = ",") %>%
            sapply(., function(x) x[1])
        )
      ) %>%
      dplyr::select(model, sim_iter, param, param_clean, species, site, everything())
    summ_out_remove[[i]] <- model3_sum
    
    # model 4
    n_alpha <- length(alpha_vec)
    summ_alpha <- as.list(1:n_alpha)
    for(j in 1:n_alpha){
      spp1 <- mle_aggregated(
        ambig_det = ambig_det, 
        phi_mat = phi_mat,
        z_vec = z_spp1, 
        n_visit = n.visits, 
        spp_idx = 1, 
        adjusted = TRUE, 
        alpha = alpha_vec[j]
      )
      
      spp2 <- mle_aggregated(
        ambig_det = ambig_det, 
        phi_mat = phi_mat,
        z_vec = z_spp2, 
        n_visit = n.visits, 
        spp_idx = 2, 
        adjusted = TRUE,
        alpha = alpha_vec[j]
      ) 
      
      model4_sum <- bind_rows(
        spp1$summ_df, spp2$summ_df
      ) %>%
        mutate(
          sim_iter = i, 
          species = c(rep(1, nrow(spp1$summ_df)), rep(2, nrow(spp2$summ_df))),
          siteID = rep(seq(1:sum(n.sites1,n.sites2)), 2),
          model = "MLESite"
        ) %>%
        dplyr::select(
          sim_iter, species, model, everything()
        )
      model4_sum$alpha = alpha_vec[j]
      summ_alpha[[j]] <- model4_sum
    }
    
    summ_out_mlesite[[i]] <- do.call(rbind, summ_alpha)
    
    print(i)
    
    # output intermittently
    if(i %% 5 == 0){
      summ_df_ctdet <- do.call(rbind, summ_out_CtDet)
      summ_df_rem <- do.call(rbind, summ_out_remove)
      summ_df_naive <- do.call(rbind, summ_out_naive)
      summ_df_site <- do.call(rbind, summ_out_mlesite)
      
      tmp <- list("CtDetection"=summ_df_ctdet,"Remove"= summ_df_rem,"Naive"=summ_df_naive, "MLESite" = summ_df_site)
      saveRDS(tmp, file = paste0("partial_runs_", i, ".rds"))
    }
  }
  
  summ_df_ctdet <- do.call(rbind, summ_out_CtDet)
  summ_df_rem <- do.call(rbind, summ_out_remove)
  summ_df_naive <- do.call(rbind, summ_out_naive)
  summ_df_site <- do.call(rbind, summ_out_mlesite)
  
  # return(tab)
  return(list("CtDetection"=summ_df_ctdet,"Remove"= summ_df_rem,"Naive"=summ_df_naive, "MLESite" = summ_df_site))
}


MLESite_sim <- function(n_iter, 
                        thetas = matrix(c(.9,.1,.35,.65),
                                        nrow=2,
                                        ncol=2,byrow=TRUE), 
                        n.sites1 = 5,
                        n.sites2 = 50,
                        n.visits = 16,
                        psi.vec = c(.25,.75),
                        lambda.vec = c(0.3,10), 
                        alpha_vec = c(0.05, 0.1, 0.15, 0.2)){
    # set up storage
    n_alpha <- length(alpha_vec)
    summ_out_mlesite <- as.list(1:n_iter)
    # mle site for both spp
    for(i in 1:n_iter){
      sim_dat <- sim_data(n.sites1 = n.sites1,
                          n.sites2 = n.sites2,
                          n.visits = n.visits,
                          psi.vec = psi.vec,
                          lambda.vec = lambda.vec, 
                          theta.mat = thetas)
      test_df <- data.frame(do.call(cbind, 
                                    sim_dat$MLE_Britzke_data[c(3:4, 7:8)]))
      z_spp1 <- sim_dat$MLE_Britzke_data$Z1
      z_spp2 <- sim_dat$MLE_Britzke_data$Z2
      ambig_det <- test_df[,3:4]
      
      phi_mat <- t(thetas) # needed to define phi_mat
      # model 4
      summ_alpha <- as.list(1:n_alpha)
      for(j in 1:n_alpha){
        spp1 <- mle_aggregated(
          ambig_det = ambig_det, 
          phi_mat = phi_mat,
          z_vec = z_spp1, 
          n_visit = n.visits, 
          spp_idx = 1, 
          adjusted = TRUE, 
          alpha = alpha_vec[j]
        )
        
        spp2 <- mle_aggregated(
          ambig_det = ambig_det, 
          phi_mat = phi_mat,
          z_vec = z_spp2, 
          n_visit = n.visits, 
          spp_idx = 2, 
          adjusted = TRUE,
          alpha = alpha_vec[j]
        ) 
      
        model4_sum <- bind_rows(
          spp1$summ_df, spp2$summ_df
        ) %>%
          mutate(
            sim_iter = i, 
            species = c(rep(1, nrow(spp1$summ_df)), rep(2, nrow(spp2$summ_df))),
            siteID = rep(seq(1:sum(n.sites1,n.sites2)), 2),
            model = "MLESite"
          ) %>%
          dplyr::select(
            sim_iter, species, model, everything()
          )
        model4_sum$alpha = alpha_vec[j]
        summ_alpha[[j]] <- model4_sum
      }

    summ_out_mlesite[[i]] <- do.call(rbind, summ_alpha)
    print(i)
    }
    summ_df_site <- do.call(rbind, summ_out_mlesite)
  return(summ_df_site)
}

# mleSite_scen1 <- MLESite_sim(n_iter = 50) 

# check that alpha is working as it should...
# mleSite_scen1 %>% filter(siteID == 1 & species == 1)


