
################################################################################
#  04_RICLPM_MULTIVARIATE: RI-CLPM with cumulatively added constructs
#  Input: imp_detrend_list.Rdata + imp.Rdata (for sample size)
#
#  Per state component (e.g. state_e), models are fitted successively:
#    Model 0: state_e (univariate)
#    Model 1: state_e + construct 1
#    Model 2: state_e + construct 1 + construct 2
#    Model 3: state_e + construct 1 + construct 2 + construct 3
#    ...
#
#  Each model:
#    - Random intercept per construct
#    - Within-person AR(1) per construct (equal over time)
#    - Cross-lagged effects between all constructs (equal over time)
#    - Innovation variances per construct (equal over time)
#    - Innovation covariances between constructs (equal over time)
#
#  Fitting: On each of the m covariance matrices, then Rubin's Rules
################################################################################

library(lavaan)
library(Matrix)  # for nearPD

# ---- 1. Load ----------------------------------------------------------------  
sim <- "02"

folder.main.pre <- file.path(base_path, "output", sim)
folder.main     <- file.path(folder.main.pre)
if (!dir.exists(folder.main)) dir.create(folder.main, recursive = TRUE)

d2 <- read.csv(file.path(base_path, "data", "wide_dat_S2.csv"))

load( file.path( folder.main, "imp.Rdata" ) )
load( file.path( folder.main, "imp_detrend_list.Rdata" ) )

n_obs <- nrow(d2)
m_imp <- length(imp_detrend_list)
n_time <- 32

# ---- PARAMETER: Number of datasets to use ------------------------------------
m_use <- 20         # How many of the m_imp datasets to use?
# m_use = m_imp for all, e.g. 3 or 5 for quick tests
m_use <- min(m_use, m_imp)
m_indices <- 1:m_use  # Which datasets (the first m_use)

cat("n =", n_obs, ", m_total =", m_imp, ", m_use =", m_use, "\n\n")

# ---- 2. Construct assignment ------------------------------------------------
# Per state: which constructs are added cumulatively?
# Order determines the cumulative sequence of models.

state_components <- c("state_e", "state_n", "state_o", "state_a", "state_c")

other_constructs_map <- list(
  state_e = c("sociability", "happiness",  "tiredness", "typicality", "negativity"),
  state_n = c("stress", "happiness", "negativity", "duty", "typicality"),
  state_o = c("happiness", "sociability", "tiredness", "duty", "negativity"),
  state_a = c("sociability", "happiness", "tiredness", "hunger", "stress"),
  state_c = c("sociability", "duty", "tiredness", "happiness",  "negativity")
)

trait_map <- c(
  state_e = "trait_e",
  state_n = "trait_n",
  state_o = "trait_o",
  state_a = "trait_a",
  state_c = "trait_c"
)

# ---- 3. Model syntax generator (multivariate) --------------------------------
#
# Generates lavaan syntax for a multivariate RI-CLPM with
# any number of constructs and equality constraints over time.
generate_mv_riclpm_syntax <- function(construct_names, n_time, trait_names = NULL) {
  
  n_c <- length(construct_names)
  lines <- c()
  clabs <- paste0("c", 1:n_c)
  
  # ---- Random Intercepts ----
  lines <- c(lines, "# === Random Intercepts ===")
  for (k in 1:n_c) {
    loadings <- paste0("1*", construct_names[k], "_", 1:n_time, collapse = " + ")
    lines <- c(lines, paste0("RI_", clabs[k], " =~ ", loadings))
  }
  lines <- c(lines, "")
  
  # ---- Within-Person latent variables ----
  lines <- c(lines, "# === Within-Person Variables ===")
  for (k in 1:n_c) {
    for (t in 1:n_time) {
      lines <- c(lines, paste0("w", clabs[k], "_", t, " =~ 1*",
                               construct_names[k], "_", t))
    }
  }
  lines <- c(lines, "")
  
  # ---- Means: RI free, within = 0, manifest = 0 ----
  lines <- c(lines, "# === Means ===")
  for (k in 1:n_c) {
    lines <- c(lines, paste0("RI_", clabs[k], " ~ rim_", clabs[k], "*1"))
  }
  for (k in 1:n_c) {
    for (t in 1:n_time) {
      lines <- c(lines, paste0("w", clabs[k], "_", t, " ~ 0*1"))
    }
  }
  for (k in 1:n_c) {
    for (t in 1:n_time) {
      lines <- c(lines, paste0(construct_names[k], "_", t, " ~ 0*1"))
    }
  }
  if (!is.null(trait_names)) {
    for (tr in trait_names) {
      lines <- c(lines, paste0(tr, " ~ mt_", tr, "*1"))
    }
  }
  lines <- c(lines, "")
  
  # ---- AR effects with night lags ----
  lines <- c(lines, "# === AR Effects (equality constraints, night lag = 0) ===")
  for (k in 1:n_c) {
    ar_label <- paste0("ar_", clabs[k])
    for (t in 2:n_time) {
      lag_label <- ifelse((t - 1) %% 4 == 0, "0", ar_label)
      lines <- c(lines, paste0("w", clabs[k], "_", t, " ~ ", lag_label,
                               "*w", clabs[k], "_", t-1))
    }
  }
  lines <- c(lines, "")
  
  # ---- Cross-lagged effects with night lags ----
  if (n_c > 1) {
    lines <- c(lines, "# === Cross-Lagged Effects (equality constraints, night lag = 0) ===")
    for (k in 1:n_c) {
      for (j in 1:n_c) {
        if (k == j) next
        cl_label <- paste0("cl_", clabs[j], "_", clabs[k])
        for (t in 2:n_time) {
          lag_label <- ifelse((t - 1) %% 4 == 0, "0", cl_label)
          lines <- c(lines, paste0("w", clabs[k], "_", t, " ~ ", lag_label,
                                   "*w", clabs[j], "_", t-1))
        }
      }
    }
    lines <- c(lines, "")
  }
  
  # ---- Variance wave 1 (free per construct) ----
  lines <- c(lines, "# === Variance Wave 1 (free) ===")
  for (k in 1:n_c) {
    lines <- c(lines, paste0("w", clabs[k], "_1 ~~ v1_", clabs[k],
                             "*w", clabs[k], "_1"))
  }
  lines <- c(lines, "")
  
  # ---- Covariances wave 1 (free) ----
  if (n_c > 1) {
    lines <- c(lines, "# === Covariances Wave 1 (free) ===")
    for (k in 1:(n_c - 1)) {
      for (j in (k + 1):n_c) {
        lines <- c(lines, paste0("w", clabs[k], "_1 ~~ cov1_", clabs[k], clabs[j],
                                 "*w", clabs[j], "_1"))
      }
    }
    lines <- c(lines, "")
  }
  
  # ---- Innovation variances (equal over time, per construct) ----
  lines <- c(lines, "# === Innovation Variances (equality constraints) ===")
  for (k in 1:n_c) {
    innov_label <- paste0("innov_", clabs[k])
    for (t in 2:n_time) {
      lines <- c(lines, paste0("w", clabs[k], "_", t, " ~~ ", innov_label,
                               "*w", clabs[k], "_", t))
    }
  }
  lines <- c(lines, "")
  
  # ---- Innovation covariances (equal over time) ----
  if (n_c > 1) {
    lines <- c(lines, "# === Innovation Covariances (equality constraints) ===")
    for (k in 1:(n_c - 1)) {
      for (j in (k + 1):n_c) {
        icov_label <- paste0("icov_", clabs[k], clabs[j])
        for (t in 2:n_time) {
          lines <- c(lines, paste0("w", clabs[k], "_", t, " ~~ ", icov_label,
                                   "*w", clabs[j], "_", t))
        }
      }
    }
    lines <- c(lines, "")
  }
  
  # ---- RI variances and covariances ----
  lines <- c(lines, "# === RI Variances and Covariances ===")
  for (k in 1:n_c) {
    ri_var_label <- paste0("RIvar_", clabs[k])
    lines <- c(lines, paste0("RI_", clabs[k], " ~~ ", ri_var_label, "*RI_", clabs[k]))
  }
  if (n_c > 1) {
    for (k in 1:(n_c - 1)) {
      for (j in (k + 1):n_c) {
        ri_cov_label <- paste0("RIcov_", clabs[k], clabs[j])
        lines <- c(lines, paste0("RI_", clabs[k], " ~~ ", ri_cov_label, "*RI_", clabs[j]))
      }
    }
  }
  lines <- c(lines, "")
  
  # ---- Trait: variance + covariance with RIs ----
  if (!is.null(trait_names)) {
    lines <- c(lines, "# === Trait Variances and Covariances with RIs ===")
    for (tr in trait_names) {
      lines <- c(lines, paste0(tr, " ~~ v_", tr, "*", tr))
      for (k in 1:n_c) {
        lines <- c(lines, paste0("RI_", clabs[k], " ~~ c_", clabs[k], "_", tr, "*", tr))
      }
    }
    lines <- c(lines, "")
  }
  
  # ---- Observed residual variances = 0 ----
  lines <- c(lines, "# === Observed Residual Variances = 0 ===")
  for (k in 1:n_c) {
    for (t in 1:n_time) {
      lines <- c(lines, paste0(construct_names[k], "_", t, " ~~ 0*",
                               construct_names[k], "_", t))
    }
  }
  lines <- c(lines, "")
  
  # ---- RI orthogonal to within ----
  lines <- c(lines, "# === RI Orthogonal to Within ===")
  for (k in 1:n_c) {
    for (j in 1:n_c) {
      for (t in 1:n_time) {
        lines <- c(lines, paste0("RI_", clabs[k], " ~~ 0*w", clabs[j], "_", t))
      }
    }
  }
  
  paste(lines, collapse = "\n")
}

# ---- 4. Rubin's Rules -------------------------------------------------------
rubins_rules <- function(estimates, se_list) {
  m <- nrow(estimates)
  qbar <- colMeans(estimates)
  ubar <- colMeans(se_list^2)
  b <- apply(estimates, 2, var)
  t_var <- ubar + (1 + 1/m) * b
  se_pooled <- sqrt(t_var)
  z <- qbar / se_pooled
  p <- 2 * pnorm(-abs(z))
  data.frame(
    est = qbar, se = se_pooled, z = z, p = p,
    fmi = (b + b/m) / t_var,
    stringsAsFactors = FALSE
  )
}

# ---- 5. Fit models (cumulatively) -------------------------------------------
cat("=== RI-CLPM MULTIVARIATE (cumulative) ===\n\n")

t_global_start <- Sys.time()

results_all <- list()
fit_measures_all <- list()
neg_var_list <- list()

for (state_comp in state_components) {
  
  cat("============================================================\n")
  cat("STATE:", toupper(state_comp), "\n")
  cat("============================================================\n\n")
  
  other_csts <- other_constructs_map[[state_comp]]
  trait_var <- trait_map[state_comp]
  
  results_all[[state_comp]] <- list()
  fit_measures_all[[state_comp]] <- list()
  
  for (n_add in 0:length(other_csts)) {
    
    if (n_add == 0) {
      current_constructs <- state_comp
    } else {
      current_constructs <- c(state_comp, other_csts[1:n_add])
    }
    model_label <- paste(current_constructs, collapse = " + ")
    
    cat("--- Model", n_add, ":", model_label, "---\n"); flush.console()
    
    # Variable names (without trait)
    all_vars <- c()
    for (cst in current_constructs) {
      all_vars <- c(all_vars, paste0(cst, "_", 1:n_time))
    }
    
    # Variable names (with trait for covariance matrix)
    all_vars_with_trait <- c(all_vars, trait_var)
    
    # Model syntax (with trait)
    model_syntax <- generate_mv_riclpm_syntax(current_constructs, n_time,
                                              trait_names = trait_var)
    
    # Save syntax (for all state components)
    syntax_file <- file.path(folder.main,
                             paste0("08_lavaan_syntax_", state_comp, "_",
                                    length(current_constructs), "var.txt"))
    writeLines(c(
      paste("# RI-CLPM Syntax for:", model_label),
      paste("# State:", state_comp),
      paste("# Trait:", trait_var),
      paste("# Number of constructs:", length(current_constructs)),
      paste("# Generated at:", Sys.time()),
      paste("# Constructs:", paste(current_constructs, collapse = ", ")),
      paste("# Time points:", n_time),
      "",
      model_syntax
    ), syntax_file)
    cat("  -> Syntax saved:", syntax_file, "\n")
    
    estimates_mat <- NULL
    se_mat <- NULL
    std_estimates_mat <- NULL  
    std_se_mat <- NULL        
    fit_indices <- list()
    n_converged <- 0
    fit_times <- c()
    
    t_model_start <- Sys.time()
    
    for (m_idx in m_indices) {
      
      cat( paste0( m_idx, " " ) ); flush.console()
      
      t_fit_start <- Sys.time()
      
      # Compute covariance matrix + nearPD (incl. trait)
      cov_sub <- cov(imp_detrend_list[[m_idx]][, all_vars_with_trait], use = "everything")
      eigs <- eigen(cov_sub, only.values = TRUE)$values
      if (any(eigs <= 1e-10)) {
        cov_sub <- as.matrix(nearPD(cov_sub, corr = FALSE, keepDiag = TRUE,
                                    ensureSymmetry = TRUE)$mat)
        dimnames(cov_sub) <- list(all_vars_with_trait, all_vars_with_trait)
      }
      
      # Compute means vector (incl. trait)
      means_vec <- colMeans(imp_detrend_list[[m_idx]][, all_vars_with_trait], na.rm = TRUE)
      
      # lavaan fit
      fit <- tryCatch({
        lavaan(
          model = model_syntax,
          sample.cov = cov_sub,
          sample.mean = means_vec,
          sample.nobs = n_obs,
          estimator = "ML",
          meanstructure = TRUE,
          int.ov.free = FALSE
        )
      }, error = function(e) {
        cat("    m =", m_idx, ": ERROR -", e$message, "\n")
        NULL
      })
      
      t_fit_end <- Sys.time()
      fit_times <- c(fit_times, as.numeric(difftime(t_fit_end, t_fit_start, units = "secs")))
      
      if (is.null(fit) || !lavInspect(fit, "converged")) {
        if (!is.null(fit)) cat("    m =", m_idx, ": Not converged\n")
        next
      }
      
      n_converged <- n_converged + 1
      
      # Save fit object for first imputed dataset only
      if (m_idx == 1) {
        fit_file <- file.path(folder.main,
                              paste0("fit_", state_comp, "_",
                                     length(current_constructs), "var.Rdata"))
        save(fit, file = fit_file)
        cat("  -> Fit object saved:", fit_file, "\n")
      }
      
      
      pe <- parameterEstimates(fit)
      
      var_params <- pe[pe$op == "~~" & pe$lhs == pe$rhs, ]
      neg_vars <- var_params[var_params$est < 0, ]
      if (nrow(neg_vars) > 0) {
        neg_vars$state_comp <- state_comp
        neg_vars$model_label <- model_label
        neg_vars$n_constructs <- length(current_constructs)
        neg_vars$m_idx <- m_idx
        neg_var_list[[length(neg_var_list) + 1]] <- neg_vars
      }
      
      labeled <- pe[pe$label != "", ]
      labeled <- labeled[!duplicated(labeled$label), ]
      
      param_names <- labeled$label
      param_est <- labeled$est
      param_se <- labeled$se
      
      if (is.null(estimates_mat)) {
        estimates_mat <- matrix(NA, nrow = m_use, ncol = length(param_names),
                                dimnames = list(NULL, param_names))
        se_mat <- estimates_mat
      }
      
      estimates_mat[which(m_indices == m_idx), ] <- param_est
      se_mat[which(m_indices == m_idx), ] <- param_se
      
      std_sol <- tryCatch(
        standardizedSolution(fit),
        error = function(e) NULL
      )
      if (!is.null(std_sol)) {
        std_labeled <- std_sol[std_sol$label != "", ]
        std_labeled <- std_labeled[!duplicated(std_labeled$label), ]
        std_names <- std_labeled$label
        std_est   <- std_labeled$est.std
        std_se    <- std_labeled$se
        if (is.null(std_estimates_mat)) {
          std_estimates_mat <- matrix(NA, nrow = m_use, ncol = length(std_names),
                                      dimnames = list(NULL, std_names))
          std_se_mat <- std_estimates_mat
        }
        std_estimates_mat[which(m_indices == m_idx), ] <- std_est
        std_se_mat[which(m_indices == m_idx), ]        <- std_se
      }
      
      fm <- tryCatch(
        fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")),
        error = function(e) NULL
      )
      if (!is.null(fm)) fit_indices[[m_idx]] <- fm
    }
    
    t_model_end <- Sys.time()
    model_runtime <- as.numeric(difftime(t_model_end, t_model_start, units = "secs"))
    
    cat("\n  Converged:", n_converged, "/", m_use, "\n")
    cat("  Runtime: ", sprintf("%.1f", model_runtime), "s total,",
        sprintf("%.1f", mean(fit_times)), "s/dataset (Median:",
        sprintf("%.1f", median(fit_times)), "s)\n")
    
    if (!is.null(estimates_mat)) {
      valid <- complete.cases(estimates_mat)
      
      if (sum(valid) >= 2) {
        pooled <- rubins_rules(estimates_mat[valid, , drop = FALSE],
                               se_mat[valid, , drop = FALSE])
        rownames(pooled) <- colnames(estimates_mat)
        
        show_params <- grep("^(ar_|cl_)", rownames(pooled), value = TRUE)
        cat("\n  Dynamic parameters (within):\n")
        print(round(pooled[show_params, c("est", "se", "p")], 4))
        
        ri_params <- grep("^RI(var|cov)_", rownames(pooled), value = TRUE)
        if (length(ri_params) > 0) {
          cat("\n  Between-person (RI) variances/covariances:\n")
          print(round(pooled[ri_params, c("est", "se", "p")], 4))
        }
        
        results_all[[state_comp]][[n_add + 1]] <- list(
          constructs = current_constructs,
          label = model_label,
          pooled = pooled,
          n_converged = n_converged,
          runtime_total = model_runtime,
          runtime_per_dataset = mean(fit_times),
          n_constructs = length(current_constructs),
          n_variables = length(all_vars)
        )
      }
    }
    
    if (!is.null(std_estimates_mat)) {
      valid_std <- complete.cases(std_estimates_mat)
      if (sum(valid_std) >= 2) {
        pooled_std <- rubins_rules(std_estimates_mat[valid_std, , drop = FALSE],
                                   std_se_mat[valid_std, , drop = FALSE])
        rownames(pooled_std) <- colnames(std_estimates_mat)
        trait_cor_params <- grep(paste0("^c_c.*_", trait_map[state_comp], "$"),
                                 rownames(pooled_std), value = TRUE)
        if (length(trait_cor_params) > 0) {
          cat("\n  RI-Trait correlations (standardized):\n")
          print(round(pooled_std[trait_cor_params, c("est", "se", "p")], 4))
        }
        if (!is.null(results_all[[state_comp]][[n_add + 1]])) {
          results_all[[state_comp]][[n_add + 1]]$pooled_std <- pooled_std
        }
      }
    }
    
    if (length(fit_indices) > 0) {
      fm_mat <- do.call(rbind, fit_indices)
      fm_avg <- colMeans(fm_mat, na.rm = TRUE)
      cat("\n  Mean fit indices:\n")
      print(round(fm_avg, 4))
      fit_measures_all[[state_comp]][[n_add + 1]] <- fm_avg
    }
    
    cat("\n")
  }
}

if (length(neg_var_list) > 0) {
  neg_var_df <- do.call(rbind, neg_var_list)
  write.csv(neg_var_df,
            file.path(folder.main, "negative_variances.csv"),
            row.names = FALSE)
} else {
  write.csv(data.frame(message = "No negative variances found"),
            file.path(folder.main, "negative_variances.csv"),
            row.names = FALSE)
}

t_global_end <- Sys.time()
global_runtime <- difftime(t_global_end, t_global_start, units = "mins")
cat("\n=== TOTAL RUNTIME:", sprintf("%.1f", as.numeric(global_runtime)), "minutes ===\n\n")

# ---- 6. Summary: Fit comparison across cumulative models --------------------
cat("\n=== SUMMARY: FIT COMPARISON ===\n\n")

for (state_comp in state_components) {
  
  cat("---", toupper(state_comp), "---\n")
  cat(sprintf("%-4s %-45s %8s %8s %8s %8s %5s %8s\n",
              "Mod", "Constructs", "CFI", "TLI", "RMSEA", "SRMR", "Conv", "Runtime"))
  cat(paste(rep("-", 95), collapse = ""), "\n")
  
  other_csts <- other_constructs_map[[state_comp]]
  
  for (n_add in 0:length(other_csts)) {
    if (n_add == 0) {
      current <- state_comp
    } else {
      current <- c(state_comp, other_csts[1:n_add])
    }
    label <- paste(current, collapse = " + ")
    if (nchar(label) > 43) label <- paste0(substr(label, 1, 40), "...")
    
    fm <- if (n_add + 1 <= length(fit_measures_all[[state_comp]])) fit_measures_all[[state_comp]][[n_add + 1]] else NULL
    res <- if (n_add + 1 <= length(results_all[[state_comp]])) results_all[[state_comp]][[n_add + 1]] else NULL
    conv <- if (!is.null(res)) res$n_converged else 0
    
    if (!is.null(fm)) {
      rt <- if (!is.null(res) && !is.null(res$runtime_total)) sprintf("%6.1fs", res$runtime_total) else "NA"
      cat(sprintf("%-4d %-45s %8.3f %8.3f %8.3f %8.3f %3d/%d %8s\n",
                  n_add, label, fm["cfi"], fm["tli"], fm["rmsea"], fm["srmr"],
                  conv, m_use, rt))
    } else {
      cat(sprintf("%-4d %-45s %8s %8s %8s %8s %3d/%d %8s\n",
                  n_add, label, "NA", "NA", "NA", "NA", conv, m_use, "NA"))
    }
  }
  cat("\n")
}

# ---- 7. Summary: Cross-lagged effects ---------------------------------------
cat("\n=== SUMMARY: CROSS-LAGGED EFFECTS ===\n\n")

for (state_comp in state_components) {
  
  cat("---", toupper(state_comp), "---\n\n")
  
  other_csts <- other_constructs_map[[state_comp]]
  
  for (n_add in 0:length(other_csts)) {
    res <- if (n_add + 1 <= length(results_all[[state_comp]])) results_all[[state_comp]][[n_add + 1]] else NULL
    if (is.null(res)) next
    
    cat("  Model", n_add, "(", res$label, "):\n")
    
    pooled <- res$pooled
    cl_params <- grep("^cl_", rownames(pooled), value = TRUE)
    
    if (length(cl_params) > 0) {
      for (cl in cl_params) {
        est <- pooled[cl, "est"]
        se  <- pooled[cl, "se"]
        p   <- pooled[cl, "p"]
        sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                               ifelse(p < 0.05, "*", "")))
        cat(sprintf("    %-20s: %7.3f (SE=%6.3f, p=%6.4f) %s\n",
                    cl, est, se, p, sig))
      }
    }
    cat("\n")
  }
}

# ---- 8. Summary: Between-person (RI) variances + trait ----------------------
cat("\n=== SUMMARY: BETWEEN-PERSON (RI) VARIANCES/COVARIANCES + TRAIT ===\n\n")

for (state_comp in state_components) {
  
  cat("---", toupper(state_comp), "---\n\n")
  
  other_csts <- other_constructs_map[[state_comp]]
  trait_var <- trait_map[state_comp]
  
  for (n_add in 0:length(other_csts)) {
    res <- if (n_add + 1 <= length(results_all[[state_comp]])) results_all[[state_comp]][[n_add + 1]] else NULL
    if (is.null(res)) next
    
    cat("  Model", n_add, "(", res$label, "):\n")
    
    pooled <- res$pooled
    
    # RI variances and covariances
    ri_params <- grep("^RI(var|cov)_", rownames(pooled), value = TRUE)
    if (length(ri_params) > 0) {
      cat("    RI Variances/Covariances:\n")
      for (ri in ri_params) {
        est <- pooled[ri, "est"]
        se  <- pooled[ri, "se"]
        p   <- pooled[ri, "p"]
        sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                               ifelse(p < 0.05, "*", "")))
        cat(sprintf("      %-20s: %7.3f (SE=%6.3f, p=%6.4f) %s\n",
                    ri, est, se, p, sig))
      }
    }
    
    # Trait variance
    trait_var_param <- grep(paste0("^v_", trait_var, "$"), rownames(pooled), value = TRUE)
    if (length(trait_var_param) > 0) {
      cat("    Trait Variance:\n")
      for (tp in trait_var_param) {
        est <- pooled[tp, "est"]
        se  <- pooled[tp, "se"]
        p   <- pooled[tp, "p"]
        sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                               ifelse(p < 0.05, "*", "")))
        cat(sprintf("      %-20s: %7.3f (SE=%6.3f, p=%6.4f) %s\n",
                    tp, est, se, p, sig))
      }
    }
    
    # RI-Trait covariances (unstandardized)
    trait_cov_params <- grep(paste0("^c_c.*_", trait_var, "$"), rownames(pooled), value = TRUE)
    if (length(trait_cov_params) > 0) {
      cat("    RI-Trait Covariances (unstandardized):\n")
      for (tc in trait_cov_params) {
        est <- pooled[tc, "est"]
        se  <- pooled[tc, "se"]
        p   <- pooled[tc, "p"]
        sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                               ifelse(p < 0.05, "*", "")))
        cat(sprintf("      %-20s: %7.3f (SE=%6.3f, p=%6.4f) %s\n",
                    tc, est, se, p, sig))
      }
    }
    
    cat("\n")
  }
}

# ---- 9. Save ----------------------------------------------------------------
cat("Saving...\n")

save(results_all, fit_measures_all,
     file = file.path( folder.main, "riclpm_multivariate_results.Rdata" ))

cat(" -> riclpm_multivariate_results.Rdata\n")

# ---- 10. Report --------------------------------------------------------------
sink(file.path( folder.main, "07_riclpm_multivariate_report.txt" ))

cat("=== RI-CLPM MULTIVARIATE (cumulative models) ===\n\n")
cat("Equality Constraints: AR, Cross-Lagged, Innovation Variances/Covariances\n")
cat("Estimator: ML\n")
cat("n =", n_obs, ", m =", m_use, "\n\n")

for (state_comp in state_components) {
  
  cat("============================================================\n")
  cat("STATE:", toupper(state_comp), "\n")
  cat("============================================================\n\n")
  
  other_csts <- other_constructs_map[[state_comp]]
  
  for (n_add in 0:length(other_csts)) {
    res <- if (n_add + 1 <= length(results_all[[state_comp]])) results_all[[state_comp]][[n_add + 1]] else NULL
    if (is.null(res)) next
    
    cat("--- Model", n_add, ":", res$label, "---\n")
    cat("Converged:", res$n_converged, "/", m_use, "\n\n")
    
    cat("Pooled parameters (Rubin's Rules):\n")
    print(round(res$pooled, 4))
    
    if (!is.null(res$pooled_std)) {
      trait_var <- trait_map[state_comp]
      trait_cor_params <- grep(paste0("^c_c.*_", trait_var, "$"),
                               rownames(res$pooled_std), value = TRUE)
      if (length(trait_cor_params) > 0) {
        cat("\nRI-Trait correlations (standardized):\n")
        print(round(res$pooled_std[trait_cor_params, c("est", "se", "p")], 4))
      }
    }
    
    fm <- if (n_add + 1 <= length(fit_measures_all[[state_comp]])) fit_measures_all[[state_comp]][[n_add + 1]] else NULL
    if (!is.null(fm)) {
      cat("\nMean fit indices:\n")
      print(round(fm, 4))
    }
    
    cat("\n\n")
  }
}

sink()

cat(" -> 07_riclpm_multivariate_report.txt\n\n")
cat("=== DONE ===\n")