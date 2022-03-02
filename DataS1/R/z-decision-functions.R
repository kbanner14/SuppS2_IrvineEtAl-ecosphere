###############
## Functions ##
###############

# for nimble model output 

z_table_nimble <- function(sim_results, z_cutoff = 0.5){
  # filter to just look at Z states and estimates
  idx_z <- which(substring(sim_results$param, 1,1) == "Z")
  sim_results <- sim_results[idx_z, ]
  # create species column
  sim_results$species <- 
    sapply(strsplit(sim_results$param, ","), 
           function(x){
             return(as.numeric(substring(x[2], 2,2)))
           }
    )
  n_iter <- max(sim_results$sim_iter)
  n_spp <- max(sim_results$species)
  n_site <- nrow(sim_results)/(n_iter*n_spp)
  # use cutoff to create z-decision based on model, 
  # could use mode, but that's same as 0.5, i think..
  sim_results$z_decision <- ifelse(sim_results$Mean > 
                                     z_cutoff, 1, 0)
  # create decision table categories--- decision.truth
  sim_results$z_decide.dg <- interaction(sim_results$z_decision, 
                                         sim_results$truth)
  # check... 
  # car::some(sim_results[,c(1:4,19,20,13,21)]) 
  # summary data frame
  df_summ <- sim_results %>% group_by(species, 
                                      sim_iter, z_decide.dg) %>% 
    summarise(count = length(z_decide.dg), 
              z_decision = unique(z_decision),
              z_dg = unique(truth), 
              # prop = length(z_decide.dg)/n_site, 
              model = model[1])
  
  df_summ <- df_summ[order(df_summ$sim_iter), ]
  
  # count # z = 1 and z = 0 for each species for each sim_iter
  df_summ_z <- sim_results %>% group_by(species, sim_iter) %>% 
    summarise(z_1 = sum(truth), z_0 = n_site - sum(truth))
  
  # join with counts from 2x2 table
  df_out <- full_join(df_summ, df_summ_z)
  # find cond prob 
  df_out$cond_prob_z <- ifelse(df_out$z_dg == 1, df_out$count/df_out$z_1, 
                               df_out$count/df_out$z_0)
  return(df_out)
}

# for MLESite output 

z_table_MLEsite <- function(sim_results, NA_combine = FALSE, add_alpha_vec = NULL){
  # create table for each iteration 
  # head(sim_results)
  if(is.null(add_alpha_vec)){
  if(NA_combine == TRUE){
    sim_results$mle_decide <- replace_na(sim_results$Y_noFP_sl, 
                                         replace = 0)
    sim_results$z_decide.dg <- interaction(sim_results$mle_decide, 
                                           sim_results$Z_true)
    n_site <- nrow(sim_results)/(max(sim_results$sim_iter)*
                                   max(sim_results$species)*length(unique(sim_results$alpha)))
    df_summ <- sim_results %>% group_by(z_decide.dg, 
                                        species, 
                                        sim_iter, alpha) %>% 
      summarise(count = length(z_decide.dg), 
                # got rid of z_decision b/c makes duplicates
                z_dg = unique(Z_true), 
                model = "MLESite", 
                alpha = unique(alpha))
    df_summ <- df_summ[order(df_summ$sim_iter), ]
  } else {
    sim_results$z_decide.dg <- interaction(sim_results$Y_noFP_sl, 
                                           sim_results$Z_true)
    # car::some(sim_results)
    n_site <- nrow(sim_results)/(max(sim_results$sim_iter)*
                                   max(sim_results$species)*length(unique(sim_results$alpha)))
    df_summ <- sim_results %>% group_by(z_decide.dg, 
                                        species, 
                                        sim_iter, alpha) %>% 
      summarise(count = length(z_decide.dg), 
                z_decision = unique(Y_noFP_sl),
                z_dg = unique(Z_true), 
                # prop = length(z_decide.dg)/n_site, C
                model = "MLESite", 
                alpha = unique(alpha))
    df_summ <- df_summ[order(df_summ$sim_iter), ]
  } 
  } else {
    # add new Y_noFP_sl for more alpha vecs 
    # find out how many rows per alpha
    n_row <- nrow(sim_results %>% filter(alpha == 0.05))
    # order the output, and grab all of the results 
    sim_ordered <- sim_results[order(sim_results$alpha), ]
    mle_results <- cbind(sim_ordered[1:n_row, c(1:8)], Y_noFP_sl = NA, 
                               sim_ordered[1:n_row, c(10:12)])
    
    n_alpha_add <- length(add_alpha_vec)
    mle_add <- data.frame()
    for(i in 1:n_alpha_add){
      mle_results$alpha <- add_alpha_vec[i]
      mle_results$Y_noFP_sl <- ifelse(mle_results$lrt_pvalue_sl < add_alpha_vec[i], 
                                      1, 0)
      mle_add <- rbind(mle_add, mle_results)
    }
    sim_results <- rbind(sim_results, mle_add)
    if(NA_combine == TRUE){
      sim_results$mle_decide <- replace_na(sim_results$Y_noFP_sl, 
                                           replace = 0)
      sim_results$z_decide.dg <- interaction(sim_results$mle_decide, 
                                             sim_results$Z_true)
      n_site <- nrow(sim_results)/(max(sim_results$sim_iter)*
                                     max(sim_results$species)*length(unique(sim_results$alpha)))
      df_summ <- sim_results %>% group_by(z_decide.dg, 
                                          species, 
                                          sim_iter, alpha) %>% 
        summarise(count = length(z_decide.dg), 
                  # got rid of z_decision b/c makes duplicates
                  z_dg = unique(Z_true), 
                  model = "MLESite", 
                  alpha = unique(alpha))
      df_summ <- df_summ[order(df_summ$sim_iter), ]
    } else {
      sim_results$z_decide.dg <- interaction(sim_results$Y_noFP_sl, 
                                             sim_results$Z_true)
      # car::some(sim_results)
      n_site <- nrow(sim_results)/(max(sim_results$sim_iter)*
                                     max(sim_results$species)*length(unique(sim_results$alpha)))
      df_summ <- sim_results %>% group_by(z_decide.dg, 
                                          species, 
                                          sim_iter, alpha) %>% 
        summarise(count = length(z_decide.dg), 
                  z_decision = unique(Y_noFP_sl),
                  z_dg = unique(Z_true), 
                  # prop = length(z_decide.dg)/n_site, C
                  model = "MLESite", 
                  alpha = unique(alpha))
      df_summ <- df_summ[order(df_summ$sim_iter), ]
    } 

    } 
  
    
    # count # z = 1 and z = 0 for each species for each sim_iter
    df_summ_z <- sim_results %>% group_by(species, sim_iter, alpha) %>% 
      summarise(z_1 = sum(Z_true), z_0 = n_site - sum(Z_true))
    
    # join with counts from 2x2 table
    df_out <- full_join(df_summ, df_summ_z)
    # find cond prob 
    df_out$cond_prob_z <- ifelse(df_out$z_dg == 1, df_out$count/df_out$z_1, 
                                 df_out$count/df_out$z_0)
  return(df_out)
}

# create one big df with all models
table_process <- function(sim_run, z_cutoff = 0.5, NA_combine = FALSE, 
                          add_alpha_vec = NULL){
  # set up data frame to summarize decision results
  # create lists of tables for each model
  ct_det <- z_table_nimble(sim_results = sim_run$CtDetection, 
                           z_cutoff = z_cutoff)
  ct_det$alpha <- NA
  remove <- z_table_nimble(sim_results = sim_run$Remove, 
                           z_cutoff = z_cutoff)
  ct_det$alpha <- NA
  naive <- z_table_nimble(sim_results = sim_run$Naive, 
                          z_cutoff = z_cutoff)
  ct_det$alpha <- NA
  mle <- z_table_MLEsite(sim_results = sim_run$MLESite, 
                         NA_combine = NA_combine, add_alpha_vec = add_alpha_vec)
  df_out <- rbind(ct_det, remove, naive, mle)		
  df_out$model <- factor(df_out$model)
  df_out$species <- factor(df_out$species)
  df_out$alpha <- factor(df_out$alpha)
  return(df_out)
}

compare_threshold <- function(sim_run, z_cutoffs = c(0.5, 0.6, 0.7, 0.8), 
                              NA_combine = FALSE, add_alpha_vec = NULL) {
  df_compare <- data.frame()
  for(i in 1:length(z_cutoffs)){
    df <- table_process(sim_run = sim_run, z_cutoff = z_cutoffs[i], 
                        NA_combine = NA_combine, add_alpha_vec = add_alpha_vec)
    df$z_cutoff <- z_cutoffs[i]
    df_compare <- rbind(df_compare, df)
  }
  df_compare$z_cutoff <- factor(df_compare$z_cutoff)
  df_compare[which(df_compare$model == "MLESite"), "z_cutoff"] <- NA
  df_nimble <- df_compare %>% filter(model != "MLESite")
  df_mle <- df_compare %>% filter(model == "MLESite")
  out_list <- list("nimble_compare" = df_nimble, 
                   "MLESite_compare" = df_mle)
  return(out_list)
}

# #### development 
# scen1_v16 <- readRDS("../RData/rds-files_may12/scen1_visits16.rds")
# # set.seed(6121)
# # test <- car::some(scen1_v16$MLESite)
# # sim_results <- test
# sim_results <- scen1_v16$MLESite
# library(tidyverse)
# # check
# df_outNA <- z_table_MLEsite(sim_results, NA_combine = TRUE)
# df_out <- z_table_MLEsite(sim_results, NA_combine = FALSE)
# 
# for(i in 1:50){
#   test <- filter(df_out, species == 1, sim_iter == i, alpha == 0.05)  
#   print(paste("iteration", i, "; nrow = ", nrow(test)))
# }
# 
# filter(df_out, species == 1, sim_iter == 41, alpha == 0.05)
# filter(df_outNA, species == 1, sim_iter == 41, alpha == 0.05)
# # test table process
# sim_run <- scen1_v16
# 
# tp_NA <- table_process(sim_run = sim_run, NA_combine = TRUE)
# head(tp_NA)
# tp <- table_process(sim_run = sim_run, NA_combine = FALSE)
# head(tp)
# # looks good
# filter(tp_NA, species == 1, sim_iter == 41, 
#        alpha == 0.05 | is.na(alpha) == TRUE)
# filter(tp, species == 1, sim_iter == 41, 
#        alpha == 0.05 | is.na(alpha) == TRUE)
# 
# # now try compare_threshold 
# 
# ct_NA <- compare_threshold(sim_run = sim_run, NA_combine = TRUE)
# ct <- compare_threshold(sim_run = sim_run, NA_combine = FALSE)
# 
# filter(ct_NA$MLESite_compare, species == 1, sim_iter == 41, 
#        alpha == 0.05 | is.na(alpha) == TRUE)
# filter(ct$MLESite_compare, species == 1, sim_iter == 41, 
#        alpha == 0.05 | is.na(alpha) == TRUE)
# 
# match <- filter(ct_NA$nimble_compare, species == 1, sim_iter == 41, 
#        alpha == 0.05 | is.na(alpha) == TRUE) == filter(ct$nimble_compare, species == 1, sim_iter == 41, 
#        alpha == 0.05 | is.na(alpha) == TRUE) 
# sum(match == FALSE, na.rm = T)
# 
# # looks like everything is working! 
