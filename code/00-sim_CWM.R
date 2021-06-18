

# Functions for density and random deviates CWM ---------------------------

dCWM <- function(y, X, tau, beta, mu, Sigma, sigma, log = FALSE) {
  G <- ncol(mu)
  regr_log_density <-
    sapply(1:G, function(g)
      dnorm(
        x = y,
        mean = cbind(1, X) %*% beta[, g],
        sd = sqrt(sigma[g]),
        log = TRUE
      ))
  X_log_density <-
    sapply(1:G, function(g)
      mvtnorm::dmvnorm(
        x = X,
        mean = mu[, g],
        sigma = Sigma[, , g],
        log = TRUE
      ))
  log_comp_CWM <-
    sweep(
      (regr_log_density + X_log_density),
      MARGIN = 2,
      STATS = log(tau),
      FUN = "+"
    )
  res <- rowSums(exp(log_comp_CWM))
  if (log == TRUE) {
    res <- log(res)
  }
  res
}


rCWM <- function(n, tau, beta, mu, Sigma, sigma, ng_fixed=TRUE) {
  G <- ncol(mu)
  p <- nrow(mu)
  if(ng_fixed){
    n_g <- round(n*tau,0)
    class <- rep(1:G,n_g)
  } else{
    class <- sample(x = 1:G,size = n,replace = TRUE,prob = tau)
  }
  
  if (p == 1) {
    X <-
      as.matrix(sapply(1:n, function(i)
        MASS::mvrnorm(
          n = 1, mu = mu[, class[i]], Sigma = Sigma[, , class[i]]
        )))
    
  } else {
    X <-
      t(sapply(1:n, function(i)
        MASS::mvrnorm(
          n = 1, mu = mu[, class[i]], Sigma = Sigma[, , class[i]]
        )))
  }
  
    y <-
      sapply(1:n, function(i)
        rnorm(
          n = 1,
          mean = c(1, X[i, ]) %*% beta[,class[i]],
          sd = sqrt(sigma[class[i]])
        ))
  
  res <- data.frame(y,X)
  attr(res,"clust") <- class
  res
}

dCWM_KL <- function(yX, tau, beta, mu, Sigma, sigma, log = FALSE) {
  G <- ncol(mu)
  regr_log_density <-
    sapply(1:G, function(g)
      dnorm(
        x = yX[1],
        mean = c(1, yX[-1]) %*% beta[, g],
        sd = sqrt(sigma[g]),
        log = TRUE
      ))
  X_log_density <-
    sapply(1:G, function(g)
      mvtnorm::dmvnorm(
        x = yX[-1],
        mean = mu[, g],
        sigma = Sigma[, , g],
        log = TRUE
      ))
  log_comp_CWM <- regr_log_density + X_log_density + log(tau)
  res <- sum(exp(log_comp_CWM))
  if(log==TRUE){
    res <- log(res)
  }
  res
}


# Function for relabeling via max density ---------------------------------

estep_CWM <- function(y, X, tau, beta, mu, Sigma, sigma) {
  G <- ncol(mu)
  p <- nrow(mu)
  regr_log_density <-
    sapply(1:G, function(g)
      dnorm(
        x = y,
        mean = cbind(1, X) %*% beta[, g],
        sd = sqrt(sigma[g]),
        log = TRUE
      ))
  if(p>1){
  X_log_density <-
    sapply(1:G, function(g)
      mvtnorm::dmvnorm(
        x = X,
        mean = mu[, g],
        sigma = Sigma[, , g],
        log = TRUE
      ))
  } else {
    X_log_density <-
      sapply(1:G, function(g)
        stats::dnorm(
          x = X,
          mean = mu[, g],
          sd = sqrt(Sigma[, , g]),
          log = TRUE
        ))
  }
  
  log_comp_CWM <-
    sweep(
      (regr_log_density + X_log_density),
      MARGIN = 2,
      STATS = log(tau),
      FUN = "+"
    )
  
  z_max <- apply(log_comp_CWM, 1, max)
  log_density <- z_max + log(rowSums(exp(log_comp_CWM - z_max)))
  z <- exp(log_comp_CWM - log_density)
  z
}

cluster_relabeler_via_CWM_density <-
  function(X,
           group_representative_fake_obs,
           old_partition,
           new_partition,
           new_model,
           old_model) {
  # 
  flag_trimmed_units_old <- any(old_partition==0)
  flag_trimmed_units_new <- any(new_partition==0)
  
  K_old <- length(unique(old_partition))
  K_new <- length(unique(new_partition))
  
  if(flag_trimmed_units_old){
    K_old <- K_old-1
  }
  
  if(flag_trimmed_units_new){
    K_new <- K_new-1
  }
  
  tentative_relabeler <- mclust::map(estep_CWM(
    y = group_representative_fake_obs[,1],
    X = group_representative_fake_obs[,-1,drop=FALSE],
    tau = new_model$cw,
    beta = t(new_model$b),
    mu = t(new_model$center),
    Sigma = new_model$sigma,
    sigma = new_model$v
  ))
  
  if(any(tentative_relabeler==0)){ # when this happens it is due to spurious clusters, I substitute it with the missing label
    sporious_cluster_indicator <- setdiff(1:K_old,tentative_relabeler)
    tentative_relabeler[tentative_relabeler==0] <- sporious_cluster_indicator
  }
  
  if(any(duplicated(tentative_relabeler)) & K_new==K_old){ 
    sporious_cluster_indicator <- setdiff(1:K_old,tentative_relabeler) # when this happens it is due to spurious clusters, I substitute it with the missing label
    tentative_relabeler[duplicated(tentative_relabeler)] <- sporious_cluster_indicator
  }
  
  if(any(duplicated(tentative_relabeler)) & K_new>K_old){ 
    # FIXME need to recompute deepest obs for each
    warning("The clustering structure has been lost",call. = F)
  }
  # situations in which every deepest obs is in different group in the new partition
  if(K_new>K_old){
    new_group_labels <- setdiff(1:K_new,tentative_relabeler)
    tentative_relabeler <- c(tentative_relabeler,new_group_labels)
    
    new_y_representative_fake_obs <- sapply(new_group_labels,function(g)
      new_model$b[g,] %*% c(1,new_model$center[g,]))
    new_group_representative_fake_obs <- cbind(new_y_representative_fake_obs, new_model$center[new_group_labels,,drop=FALSE])
    group_representative_fake_obs <- rbind(group_representative_fake_obs, new_group_representative_fake_obs)
  }
  
  out <- list()
  out$relabeler <- tentative_relabeler
  
  if(K_new<K_old) {
    
    which_group_is_duplicated <- tentative_relabeler[duplicated(tentative_relabeler)]
    
    duplicated_y_representative_fake_obs <- sapply(which_group_is_duplicated,function(g)
      new_model$b[g,] %*% c(1,new_model$center[g,]))
    
    duplicated_group_representative_fake_obs <- cbind(duplicated_y_representative_fake_obs, new_model$center[which_group_is_duplicated,])
    
    positions_duplicated <- which(tentative_relabeler==which_group_is_duplicated) # pay attention, this works only if K_new == K_old-1
    
    NEW_group_representative_fake_obs <- group_representative_fake_obs[-positions_duplicated[2],]
    
    NEW_group_representative_fake_obs[positions_duplicated[1],] <- duplicated_group_representative_fake_obs
    
    NEW_tentative_relabeler <- mclust::map(estep_CWM( # this relabeler takes into account the fact that groups from the K_old groups in the old partition have been merged in K_new groups
      y = NEW_group_representative_fake_obs[,1],
      X = NEW_group_representative_fake_obs[,-1,drop=FALSE],
      tau = new_model$cw,
      beta = t(new_model$b),
      mu = t(new_model$center),
      Sigma = new_model$sigma,
      sigma = new_model$v
    ))
    
    new_partition_relabeled <- new_partition
    
    for (k in 1:K_old) {
      new_partition_relabeled[new_partition == k] <-
        tentative_relabeler[k]
    }
    
    out$new_partition_relabeled <- new_partition_relabeled
    out$relabeler <- NEW_tentative_relabeler
    group_representative_fake_obs <- NEW_group_representative_fake_obs
  } else if (!flag_trimmed_units_new) {
    out$new_partition_relabeled <-
      as.numeric(factor(new_partition, levels = out$relabeler))
  }
  
  out$group_representative_fake_obs <- group_representative_fake_obs
  out$trimmed_units_new <- flag_trimmed_units_new # If this is TRUE, `out$new partition relabeled` object should not be used, only `out$relabeler`
  out$trimmed_units_old <- flag_trimmed_units_old # If this is TRUE, `out$new partition relabeled` object should not be used, only `out$relabeler`
  out
}

# Functions for relabeling via depth --------------------------------------

deepest_groupwise_obs_detector <- function(X, partition){
  require(magrittr)
  N <- nrow(X)
  K <- length(unique(partition))
  obs_depth <- unlist(lapply(1:N, function(obs) depth::depth(u = X[obs,],
                                                             x = X[partition==partition[obs],,drop=FALSE],approx = TRUE)))
  # The following chain takes care of a lot of problems that may arise in computing the depth:
  # - for different groups, the same max depth value can be obtained
  # - within each class, more than one obs could have the same highest depth value
  # - tapply with which.max gives the id whithin class (not useful)
  
  data.frame(id=1:N,obs_depth,partition) %>% 
    dplyr::group_by(partition) %>% 
    dplyr::filter(obs_depth==max(obs_depth)) %>% 
    dplyr::distinct(partition,.keep_all = TRUE) %>% 
    dplyr::arrange(partition) %>% 
    dplyr::pull(id)
  
}


cluster_relabeler_via_depth <- function(X, deepest_obs, old_partition, new_partition){
  
  # 
  flag_trimmed_units_old <- any(old_partition==0)
  flag_trimmed_units_new <- any(new_partition==0)
  
  K_old <- length(unique(old_partition))
  K_new <- length(unique(new_partition))
  
  if(flag_trimmed_units_old){
    K_old <- K_old-1
  }
  
  if(flag_trimmed_units_new){
    K_new <- K_new-1
  }
  
  tentative_relabeler <- new_partition[deepest_obs]
  
  if(any(tentative_relabeler==0)){ # when this happens it is due to spurious clusters, I substitute it with the missing label
    sporious_cluster_indicator <- setdiff(1:K_old,tentative_relabeler)
    tentative_relabeler[tentative_relabeler==0] <- sporious_cluster_indicator
  }
  
  if(any(duplicated(tentative_relabeler)) & K_new==K_old){ 
    sporious_cluster_indicator <- setdiff(1:K_old,tentative_relabeler) # when this happens it is due to spurious clusters, I substitute it with the missing label
    tentative_relabeler[duplicated(tentative_relabeler)] <- sporious_cluster_indicator
  }
  
  if(any(duplicated(tentative_relabeler)) & K_new>K_old){ 
    # FIXME need to recompute deepest obs for each
    warning("The clustering structure has been lost",call. = F)
  }
  # situations in which every deepest obs is in different group in the new partition
  if(K_new>K_old){
    new_group_labels <- setdiff(1:K_new,tentative_relabeler)
    tentative_relabeler <- c(tentative_relabeler,new_group_labels)
    new_deepest_obs <- deepest_groupwise_obs_detector(X = X, partition = new_partition) # FIXME inefficient as it computes depth for all groups
    if(flag_trimmed_units_new){
      new_deepest_obs <- new_deepest_obs[-1] # remove the first deepest obs related to the trimmed units
    }
    deepest_obs <- c(deepest_obs, new_deepest_obs[new_group_labels])
  }
  
  out <- list()
  out$relabeler <- tentative_relabeler
  
  if(K_new<K_old) {
    
    # first I merge the components that were split in the old partition, and then I compute the relabeler and the new deepest obs
    
    old_partition_merged <- old_partition
    
    matched_old_groups_merged <- match(x = tentative_relabeler,table =  tentative_relabeler)
    
    
    for(k in 1:length(matched_old_groups_merged)){
      old_partition_merged[old_partition==k] <- matched_old_groups_merged[k]
    }
    
    if(flag_trimmed_units_old){
      old_partition_merged <- as.numeric(factor(old_partition_merged))-1 # minus one since the trimmed obs are denoted with 0
    }
    
    
    old_deepest_obs_merged <-
      deepest_groupwise_obs_detector(X = X, partition = old_partition_merged) 
    
    if (flag_trimmed_units_new) {
      old_deepest_obs_merged <-
        old_deepest_obs_merged[-1] # removing deepest obs related to the trimmed subset
    }
    
    
    tentative_relabeler <- new_partition[old_deepest_obs_merged] # this relabeler takes into account the fact that groups from the K_old groups in the old partition have been merged in K_new groups
    
    new_partition_relabeled <- new_partition
    
    
    
    for (k in 1:K_old) {
      new_partition_relabeled[new_partition == k] <-
        tentative_relabeler[k]
    }
    
    deepest_obs <- old_deepest_obs_merged
    out$new_partition_relabeled <- new_partition_relabeled
    out$relabeler <- tentative_relabeler
    
  } else if (!flag_trimmed_units_new) {
    out$new_partition_relabeled <-
      as.numeric(factor(new_partition, levels = out$relabeler))
  }
  
  out$deepest_obs <- deepest_obs
  out$trimmed_units_new <- flag_trimmed_units_new # If this is TRUE, `out$new partition relabeled` object should not be used, only `out$relabeler`
  out$trimmed_units_old <- flag_trimmed_units_old # If this is TRUE, `out$new partition relabeled` object should not be used, only `out$relabeler`
  out
}


# Function for building bivariate car-bike plot ---------------------------


plausible_solution_finder <- function(solutions_df, L, epsilon_ARI){ # L is the number of plausible solutions to retain
  
  criterion <- TRUE
  remaining_solutions_df <- solutions_df
  iter <- 1
  
  c_X_distinct <- unique(solutions_df$restr_factor_X)
  c_Y_distinct <- unique(solutions_df$restr_factor_Y)
  
  plausible_solutions_df <- # empty tibble
    tribble(
      ~ num_groups,
      ~ restr_factor_X,
      ~ restr_factor_Y,
      ~ model,
      ~ TBIC,
      ~ ARI_with_best,
      ~ stable_solution
    )
  
  while(criterion){
    
    best_model <- remaining_solutions_df %>%
      dplyr::filter(rank_solution == min(rank_solution))
    
    adjacent_c_X <- c_X_distinct[which(best_model$restr_factor_X==c_X_distinct)+c(-1,0,1)]
    adjacent_c_X <- adjacent_c_X[!is.na(adjacent_c_X)]
    
    adjacent_c_Y <- c_Y_distinct[which(best_model$restr_factor_Y==c_Y_distinct)+c(-1,0,1)]
    adjacent_c_Y <- adjacent_c_Y[!is.na(adjacent_c_Y)]
    
    # Keeping k equal to the one of the  best model, find the stable solutions
    # (ARI with respect to the best model is >= epsilon_ARI )
    
    similar_partitions_df <- remaining_solutions_df %>%
      dplyr::filter(num_groups == best_model$num_groups) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        ARI_with_best = mclust::adjustedRandIndex(x = model$assig, y = best_model$model[[1]]$assig),
        stable_solution = ifelse((ARI_with_best >=
                                    epsilon_ARI) &
                                   ((restr_factor_Y %in% adjacent_c_Y) &
                                      (restr_factor_X %in% adjacent_c_X)),
                                 # adjacent solutions
                                 TRUE,
                                 FALSE
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(stable_solution) %>% 
      dplyr::mutate(n_solution = iter, is_best = rank_solution == best_model$rank_solution)
    
    plausible_solutions_df <- bind_rows(plausible_solutions_df,filter(similar_partitions_df, is_best)) 
    plausible_solutions_df
    remaining_solutions_df = anti_join(
      remaining_solutions_df,
      similar_partitions_df)
    
    iter <- iter + 1
    
    if (iter > L | nrow(remaining_solutions_df) == 0) {
      criterion <- FALSE
    }
  }
  plausible_solutions_df
}

optimal_solution_finder <- function(df_plausible_solutions, epsilon_ARI=.8, threshold_n_sporious=10){
  
  I_subset <- 1:nrow(df_plausible_solutions)
  
  
  d_ARI <-
    sapply(1:nrow(df_plausible_solutions), function(i)
      purrr::map_dbl(
        1:nrow(df_plausible_solutions),
        ~ mclust::adjustedRandIndex(
          df_plausible_solutions$model[[i]]$assig,
          df_plausible_solutions$model[[.x]]$assig
        )
      ))
  
  diag(d_ARI) <- 0
  for (i in 1:nrow(df_plausible_solutions)) {
    if (i %in% I_subset) {
      spurious_solutions <- which(d_ARI[i, ] > epsilon_ARI)
      
      # Checking that the ith solution is not spurious by looking that each cluster size is > threshold_n_sporious
      if(any(table(df_plausible_solutions$model[[i]]$assig)[-1]<threshold_n_sporious)){ # the minus one eliminates the trimming subset
        I_subset <- setdiff(I_subset,i) # I eliminate that solution as it is spurious and I proceed with the search
        next
      }
      spurious_solutions[spurious_solutions>i]
      I_subset <- setdiff(I_subset, spurious_solutions)
    } else {
      next
    }
  }
  I_subset
}

solutions_ranking_identifier <-
  function(solutions_df,
           number_of_plausible_sol,
           epsilon_ARI,
           threshold_n_sporious) {
    
  
  # function that returns a list with 4 datasets, respectively containing the plausible solutions,
  # the optimal solutions, the spurious solutions and the augmented original dataset with the additional column for building the bivariate car plot in the CWRM context
  
  df_plausible_solutions <-
    plausible_solution_finder(solutions_df = solutions_df,
                              L = number_of_plausible_sol,
                              epsilon_ARI = epsilon_ARI)
  df_optimal_solutions <- 
    df_plausible_solutions %>% 
    slice(optimal_solution_finder(df_plausible_solutions = df_plausible_solutions, epsilon_ARI = epsilon_ARI, threshold_n_sporious=threshold_n_sporious)) %>% 
    rowid_to_column(var = "rank_best") %>% 
    select(-stable_solution, -n_solution,-is_best) 
  
  df_spurious_sol <-
    anti_join(df_plausible_solutions, df_optimal_solutions)
  
  df_4_plot_cars <-  solutions_df %>%
    mutate(
      is_stable = FALSE,
      is_best_interval = FALSE,
      n_sol = NA,
      is_sporious = ifelse(rank_solution %in% df_spurious_sol$rank_solution, TRUE, FALSE),
    ) %>% 
    left_join(select(df_optimal_solutions,1:6)) %>% 
    mutate(across(c(restr_factor_X, restr_factor_Y), .fns = factor)) 
  
  for(rank_best_mod in 1:nrow(df_optimal_solutions)) {
    
    current_optimal_solution <-
      slice(df_optimal_solutions, rank_best_mod)
    
    if (rank_best_mod < nrow(df_optimal_solutions)) {
      
      df_4_plot_cars <- df_4_plot_cars %>%
        rowwise() %>%
        mutate(
          ARI_with_best = mclust::adjustedRandIndex(x = model$assig,
                                                    y = current_optimal_solution$model[[1]]$assig),
          is_best_interval = ifelse(
            (rank_solution < slice(df_optimal_solutions, (rank_best_mod + 1))$rank_solution
            ) & (num_groups == current_optimal_solution$num_groups) & !(is_stable),
            TRUE,
            is_best_interval
          ),
          is_stable = ifelse((ARI_with_best >=
                                epsilon_ARI) &
                               num_groups == current_optimal_solution$num_groups,
                             TRUE,
                             is_stable
          ),
          n_sol = ifelse(
            num_groups == current_optimal_solution$num_groups &
              (is_best_interval | is_stable) & is.na(n_sol) & (!is_sporious),
            rank_best_mod,
            n_sol
          )
        ) %>% 
        ungroup()
    } else { # I do not compute best interval for the last optimal solution
      df_4_plot_cars <- df_4_plot_cars %>%
        rowwise() %>%
        mutate(
          ARI_with_best = mclust::adjustedRandIndex(x = model$assig,
                                                    y = current_optimal_solution$model[[1]]$assig),
          is_stable = ifelse((ARI_with_best >=
                                epsilon_ARI) &
                               num_groups == current_optimal_solution$num_groups,
                             TRUE,
                             is_stable
          ),
          n_sol = ifelse(
            num_groups == current_optimal_solution$num_groups &
              (is_best_interval | is_stable) & is.na(n_sol) & (!is_sporious),
            rank_best_mod,
            n_sol
          )
        ) %>% 
        ungroup()
    }
  }
  
  OUT <-
    list(
      plausible_solutions = df_plausible_solutions,
      optimal_solutions = df_optimal_solutions,
      spurious_solutions = df_spurious_sol,
      df_for_plotting = df_4_plot_cars
    )
  
  OUT
}

# Lighter version in which the optimals are immediately identified, without worrying about the spurious ----------


solutions_ranking_identifier_LIGHTER <-
  function(solutions_df,
           epsilon_ARI, number_of_optimal_sol) {
    
    
    # function that returns a list with 2 datasets, respectively containing 
    # the optimal solutions and the augmented original dataset with the additional column for building
    # the bivariate car plot in the CWRM context
    
    criterion <- TRUE
    remaining_solutions_df <- solutions_df
    iter <- 1
    
    plausible_solutions_df <- # empty tibble
      tribble(
        ~ num_groups,
        ~ restr_factor_X,
        ~ restr_factor_Y,
        ~ model,
        ~ TBIC,
        ~ ARI_with_best,
        ~ stable_solution
      )
    
    while(criterion){
      
      best_model <- remaining_solutions_df %>%
        dplyr::filter(rank_solution == min(rank_solution))
      
      # Find the similar solutions
      # (ARI with respect to the best model is >= epsilon_ARI )
      
      similar_partitions_df <- remaining_solutions_df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          ARI_with_best = mclust::adjustedRandIndex(x = model$assig, y = best_model$model[[1]]$assig),
          stable_solution = ifelse(ARI_with_best >=
                                      epsilon_ARI, 
                                   TRUE,
                                   FALSE
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::filter(stable_solution) %>% 
        dplyr::mutate(n_solution = iter, is_best = rank_solution == best_model$rank_solution)
      
      plausible_solutions_df <- bind_rows(plausible_solutions_df,filter(similar_partitions_df, is_best)) 
      remaining_solutions_df <-  anti_join(
        remaining_solutions_df,
        similar_partitions_df)
      
      iter <- iter + 1
      
      if (nrow(remaining_solutions_df) == 0 | (iter>number_of_optimal_sol)) {
        criterion <- FALSE
      }
    }
    
    df_optimal_solutions <- 
      plausible_solutions_df %>% 
      rowid_to_column(var = "rank_best") %>% 
      select(-stable_solution, -n_solution,-is_best) 
    
    df_4_plot_cars <-  solutions_df %>%
      mutate(
        is_stable = FALSE,
        is_best_interval = FALSE,
        n_sol = NA,
      ) %>% 
      left_join(select(df_optimal_solutions,1:6)) %>% 
      mutate(across(c(restr_factor_X, restr_factor_Y), .fns = factor)) 
    
    for(rank_best_mod in 1:nrow(df_optimal_solutions)) {
      
      current_optimal_solution <-
        slice(df_optimal_solutions, rank_best_mod)
      
      if (rank_best_mod < nrow(df_optimal_solutions)) {
        
        df_4_plot_cars <- df_4_plot_cars %>%
          rowwise() %>%
          mutate(
            ARI_with_best = mclust::adjustedRandIndex(x = model$assig,
                                                      y = current_optimal_solution$model[[1]]$assig),
            is_best_interval = ifelse(
              (rank_solution < slice(df_optimal_solutions, (rank_best_mod + 1))$rank_solution
              ) & (num_groups == current_optimal_solution$num_groups) & !(is_stable),
              TRUE,
              is_best_interval
            ),
            is_stable = ifelse((ARI_with_best >=
                                  epsilon_ARI) &
                                 num_groups == current_optimal_solution$num_groups,
                               TRUE,
                               is_stable
            ),
            n_sol = ifelse(
              num_groups == current_optimal_solution$num_groups &
                (is_best_interval | is_stable) & is.na(n_sol),
              rank_best_mod,
              n_sol
            )
          ) %>% 
          ungroup()
      } else { # I do not compute best interval for the last optimal solution
        df_4_plot_cars <- df_4_plot_cars %>%
          rowwise() %>%
          mutate(
            ARI_with_best = mclust::adjustedRandIndex(x = model$assig,
                                                      y = current_optimal_solution$model[[1]]$assig),
            is_stable = ifelse((ARI_with_best >=
                                  epsilon_ARI) &
                                 num_groups == current_optimal_solution$num_groups,
                               TRUE,
                               is_stable
            ),
            n_sol = ifelse(
              num_groups == current_optimal_solution$num_groups &
                (is_best_interval | is_stable) & is.na(n_sol),
              rank_best_mod,
              n_sol
            )
          ) %>% 
          ungroup()
      }
    }
    
    OUT <-
      list(
        optimal_solutions = df_optimal_solutions,
        df_for_plotting = df_4_plot_cars
      )
    
    OUT
  }
