library(tidyverse)
library(writexl)


# Full model output table -------------------------------------------------

load("~/Documents/02_Nebenprojekte/03_Illusory trait/Illusory_trait_analysis/results_cluster/01/riclpm_multivariate_results.Rdata")

# Function to create a table for one state component
create_param_table <- function(results, state_comp, trait_map) {
  
  trait_var <- trait_map[state_comp]
  
  # M1 = model 1 (index 1), M2 = model 6 (index 6, full model with 5 covariates)
  m1 <- results[[state_comp]][[1]]$pooled
  m2 <- results[[state_comp]][[6]]$pooled
  
  # Get all parameter names from M2 (superset of M1)
  all_params <- rownames(m2)
  
  # Build table
  tab <- data.frame(
    Parameter = all_params,
    M1_est = NA_real_,
    M1_se  = NA_real_,
    M1_p   = NA_real_,
    M2_est = m2[all_params, "est"],
    M2_se  = m2[all_params, "se"],
    M2_p   = m2[all_params, "p"]
  )
  
  # Fill in M1 values where available
  m1_params <- rownames(m1)
  tab[tab$Parameter %in% m1_params, "M1_est"] <- m1[m1_params[m1_params %in% all_params], "est"]
  tab[tab$Parameter %in% m1_params, "M1_se"]  <- m1[m1_params[m1_params %in% all_params], "se"]
  tab[tab$Parameter %in% m1_params, "M1_p"]   <- m1[m1_params[m1_params %in% all_params], "p"]
  
  # Round
  tab[, 2:7] <- round(tab[, 2:7], 3)
  
  return(tab)
}

# Generate tables for all dimensions
trait_map <- c(state_e = "trait_e", state_n = "trait_n", 
               state_o = "trait_o", state_a = "trait_a", state_c = "trait_c")

tables <- lapply(names(trait_map), function(sc) {
  create_param_table(results_all, sc, trait_map)
})
names(tables) <- names(trait_map)

# Print to check
print(tables$state_e)

create_label_map <- function(state_comp, other_constructs_map, trait_map) {
  
  constructs <- c(state_comp, other_constructs_map[[state_comp]])
  trait <- trait_map[state_comp]
  
  abbrev <- c(
    state_e = "E", state_n = "N", state_o = "O", state_a = "A", state_c = "C",
    happiness = "hap", sociability = "soc", positivity = "pos",
    negativity = "neg", mating = "mat", selfesteem = "se",
    duty = "dut", intellect = "int", adversity = "adv", deception = "dec",
    trait_e = "tr_E", trait_n = "tr_N", trait_o = "tr_O",
    trait_a = "tr_A", trait_c = "tr_C"
  )
  
  n_c <- length(constructs)
  clabs <- paste0("c", 1:n_c)
  cabbs <- abbrev[constructs]
  tr_abb <- abbrev[trait]
  
  map <- list()
  
  # Means
  for (k in 1:n_c) {
    map[[paste0("rim_", clabs[k])]] <- paste0("Mean[RI_", cabbs[k], "]")
  }
  map[[paste0("mt_", trait)]] <- paste0("Mean[", tr_abb, "]")
  
  # AR effects
  for (k in 1:n_c) {
    map[[paste0("ar_", clabs[k])]] <- paste0("AR[", cabbs[k], "]")
  }
  
  # CL effects: cl_cj_ck means j -> k
  for (k in 1:n_c) {
    for (j in 1:n_c) {
      if (k == j) next
      map[[paste0("cl_", clabs[j], "_", clabs[k])]] <-
        paste0("CL[", cabbs[j], "->", cabbs[k], "]")
    }
  }
  
  # First time point variances
  for (k in 1:n_c) {
    map[[paste0("v1_", clabs[k])]] <- paste0("Var_tp1[", cabbs[k], "]")
  }
  
  # First time point covariances
  for (k in 1:(n_c-1)) {
    for (j in (k+1):n_c) {
      map[[paste0("cov1_", clabs[k], clabs[j])]] <-
        paste0("Cov_tp1[", cabbs[k], ",", cabbs[j], "]")
    }
  }
  
  # Innovation variances
  for (k in 1:n_c) {
    map[[paste0("innov_", clabs[k])]] <- paste0("InnovVar[", cabbs[k], "]")
  }
  
  # Innovation covariances
  for (k in 1:(n_c-1)) {
    for (j in (k+1):n_c) {
      map[[paste0("icov_", clabs[k], clabs[j])]] <-
        paste0("InnovCov[", cabbs[k], ",", cabbs[j], "]")
    }
  }
  
  # RI variances
  for (k in 1:n_c) {
    map[[paste0("RIvar_", clabs[k])]] <- paste0("RIVar[", cabbs[k], "]")
  }
  
  # RI covariances
  for (k in 1:(n_c-1)) {
    for (j in (k+1):n_c) {
      map[[paste0("RIcov_", clabs[k], clabs[j])]] <-
        paste0("RICov[", cabbs[k], ",", cabbs[j], "]")
    }
  }
  
  # Trait variance
  map[[paste0("v_", trait)]] <- paste0("Var[", tr_abb, "]")
  
  # RI-trait covariances
  for (k in 1:n_c) {
    map[[paste0("c_", clabs[k], "_", trait)]] <-
      paste0("Cov[RI_", cabbs[k], ",", tr_abb, "]")
  }
  
  return(map)
}

# constructs in the same order as used 
other_constructs_map <- list(
  state_e = c("happiness", "sociability", "positivity", "negativity", "mating"),
  state_n = c("happiness", "selfesteem", "negativity", "positivity", "intellect"),
  state_o = c("selfesteem", "happiness", "negativity", "positivity", "deception"),
  state_a = c("happiness", "selfesteem", "negativity", "positivity", "adversity"),
  state_c = c("selfesteem", "happiness", "duty", "intellect", "adversity")
)

relabel_table <- function(tab, state_comp, other_constructs_map, trait_map) {
  label_map <- create_label_map(state_comp, other_constructs_map, trait_map)
  tab$Parameter_labeled <- sapply(tab$Parameter, function(p) {
    if (!is.null(label_map[[p]])) label_map[[p]] else p
  })
  tab <- tab %>% select(Parameter, Parameter_labeled, everything())
  return(tab)
}

reorder_table <- function(tab) {
  
  # Define parameter group order
  get_group <- function(p) {
    if (grepl("^RIVar", p)) return(1)
    if (grepl("^RICov", p)) return(2)
    if (grepl("^Var\\[tr", p)) return(3)
    if (grepl("^Cov\\[RI_", p)) return(4)
    if (grepl("^Mean\\[RI", p)) return(5)
    if (grepl("^Mean\\[tr", p)) return(6)
    if (grepl("^AR", p)) return(7)
    if (grepl("^CL", p)) return(8)
    if (grepl("^Var_tp1", p)) return(9)
    if (grepl("^Cov_tp1", p)) return(10)
    if (grepl("^InnovVar", p)) return(11)
    if (grepl("^InnovCov", p)) return(12)
    return(13)
  }
  
  tab$group <- sapply(tab$Parameter, get_group)
  tab <- tab %>% arrange(group) %>% select(-group)
  return(tab)
}

# Apply reordering
tables_labeled <- lapply(names(trait_map), function(sc) {
  tab <- relabel_table(tables[[sc]], sc, other_constructs_map, trait_map)
  tab <- tab %>%
    select(Parameter_labeled, M1_est, M1_se, M2_est, M2_se) %>%
    rename(Parameter = Parameter_labeled)
  reorder_table(tab)
})
names(tables_labeled) <- names(trait_map)

print(tables_labeled[["state_n"]])


# Export to Excel - one sheet per dimension
writexl::write_xlsx(tables_labeled, "full_output_tables_S1.xlsx")



# Fit indices, convergence, and nearPD tables -----------------------------
load("~/Documents/02_Nebenprojekte/03_Illusory trait/Illusory_trait_analysis/results_cluster/01/riclpm_multivariate_results.Rdata")
load("~/Documents/02_Nebenprojekte/03_Illusory trait/Illusory_trait_analysis/results_cluster/01/cov_checks_detrend.Rdata")
state_components <- c("state_e", "state_n", "state_o", "state_a", "state_c")

create_fit_table <- function(results_all, fit_measures_all) {
  do.call(rbind, lapply(names(results_all), function(sc) {
    do.call(rbind, lapply(c(0, 5), function(n_add) {
      model_label <- ifelse(n_add == 0, "M1", "M2")
      fm <- fit_measures_all[[sc]][[n_add + 1]]
      if (!is.null(fm)) {
        data.frame(
          Dimension = sc,
          Model = model_label,
          Chi_square = round(fm["chisq"], 2),
          df = round(fm["df"], 0),
          p = round(fm["pvalue"], 3),
          CFI = round(fm["cfi"], 3),
          TLI = round(fm["tli"], 3),
          RMSEA = round(fm["rmsea"], 3),
          SRMR = round(fm["srmr"], 3)
        )
      }
    }))
  }))
}

create_conv_table <- function(results_all) {
  do.call(rbind, lapply(names(results_all), function(sc) {
    do.call(rbind, lapply(c(0, 5), function(n_add) {
      model_label <- ifelse(n_add == 0, "M1", "M2")
      res <- results_all[[sc]][[n_add + 1]]
      if (!is.null(res)) {
        data.frame(
          Dimension = sc,
          Model = model_label,
          Converged = paste0(res$n_converged, " / 20")
        )
      }
    }))
  }))
}


create_nearPD_table <- function(cov_checks_detrend, state_components) {
  do.call(rbind, lapply(state_components, function(comp) {
    checks <- cov_checks_detrend[[comp]]
    n_applied <- sum(sapply(checks, function(x) x$nearPD_applied))
    max_changes <- sapply(checks, function(x) x$nearPD_max_change)
    min_eigs <- sapply(checks, function(x) x$min_eigenvalue)
    
    data.frame(
      Dimension = comp,
      nearPD_applied = paste0(n_applied, " / 20"),
      Max_abs_change = sprintf("%.2e", max(max_changes)),
      Smallest_eigenvalue_min = sprintf("%.2e", min(min_eigs)),
      Smallest_eigenvalue_max = sprintf("%.2e", max(min_eigs))
    )
  }))
}


fit_table_S1   <- create_fit_table(results_all, fit_measures_all)
conv_table_S1  <- create_conv_table(results_all)
nearPD_table_S1 <- create_nearPD_table(cov_checks_detrend, state_components)

write_xlsx(
  list(
    "Fit Indices" = fit_table_S1,
    "Convergence" = conv_table_S1,
    "nearPD" = nearPD_table_S1
  ),
  "diagnostic_tables_S1.xlsx"
)

