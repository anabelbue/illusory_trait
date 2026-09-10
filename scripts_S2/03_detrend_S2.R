
################################################################################
#  03_DETREND: Individual Detrending + nearPD + Covariance Matrix Checks
#  Input: imp.Rdata (from 01_imp.R)
#
#  Procedure:
#    1) Load imputed datasets
#    2) Per person and construct: remove trend (LOESS or linear)
#    3) Re-add individual means (preserve stable person differences)
#    4) Extract covariance matrices
#    5) nearPD correction (fix structural zero eigenvalues from detrending)
#    6) Check positive definiteness
#    7) Comparison: variances/covariances before vs after
#
#  NOTE: Person-wise detrending with identical time vector creates
#  structural zero eigenvalues in the covariance matrix (residuals exactly
#  orthogonal to trend basis functions). nearPD() finds the nearest
#  positive definite matrix with minimal change (Frobenius norm).
################################################################################

library(mice)
library(Matrix)   # for nearPD()

# ---- PARAMETERS --------------------------------------------------------------
detrend_method <- "loess"   # "linear" or "loess"
loess_span <- 3             # LOESS smoothing parameter
# 3 = very smooth (removes only global trends)
# 0.75 = moderate, 0.5 = flexible

# ---- 1. Load -----------------------------------------------------------------
# ============================================================
# USER CONFIGURATION - adjust this path to your setup
# ============================================================
base_path <- ""  # your workspace/repository path
# ============================================================

sim <- "02"

folder.main.pre <- file.path(base_path, "output", sim)
folder.main     <- file.path(folder.main.pre)
if (!dir.exists(folder.main)) dir.create(folder.main, recursive = TRUE)

d2 <- read.csv(file.path(base_path, "data", "wide_dat_S2.csv"))

load( file.path( folder.main, "imp.Rdata" ) )

cat("MICE object loaded: m =", imp$m, "\n\n")

constructs <- c("state_e", "state_n", "state_o", "state_c", "state_a", 
                "duty", "sociability", "negativity", "typicality", 
                "happiness", "stress", "hunger", "tiredness")
traits <- c("trait_e", "trait_n", "trait_o", "trait_a", "trait_c")
n_time <- 32

# ---- 2. Detrending function --------------------------------------------------
detrend_dataset <- function(d, constructs, n_time, method, span = 3) {
  
  d_detrend <- d
  
  for (cst in constructs) {
    
    vars <- paste0(cst, "_", 1:n_time)
    vars <- intersect(vars, names(d))
    if (length(vars) == 0) next
    
    mat <- as.matrix(d[, vars])
    
    for (i in 1:nrow(mat)) {
      y <- mat[i, ]
      if (sum(!is.na(y)) < 5) next
      
      person_mean <- mean(y, na.rm = TRUE)
      time_vec <- 1:length(y)
      
      if (method == "linear") {
        fit <- tryCatch(lm(y ~ time_vec), error = function(e) NULL)
        if (!is.null(fit)) {
          mat[i, ] <- (y - predict(fit)) + person_mean
        }
        
      } else if (method == "loess") {
        fit <- tryCatch(
          loess(y ~ time_vec, span = span, degree = 1),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          trend <- predict(fit, newdata = data.frame(time_vec = time_vec))
          mat[i, ] <- (y - trend) + person_mean
        }
      }
    }
    
    d_detrend[, vars] <- mat
  }
  
  return(d_detrend)
}

# ---- 3. Detrend all m datasets -----------------------------------------------
cat("Detrending", imp$m, "datasets (method:", detrend_method)
if (detrend_method == "loess") cat(", span =", loess_span)
cat(")...\n")

imp_detrend_list <- list()

for (m_idx in 1:imp$m) {
  d_imp <- complete(imp, action = m_idx)
  d_detrend <- detrend_dataset(d_imp, constructs, n_time,
                               method = detrend_method, span = loess_span)
  imp_detrend_list[[m_idx]] <- d_detrend
  cat("  Dataset", m_idx, "detrended\n")
}

cat("\n")

# ---- 4-6. Covariance matrices: EXTRACT, nearPD, CHECK -----------------------
cat("=== COVARIANCE MATRIX ANALYSIS (after detrending) ===\n\n")

state_components <- c("state_e", "state_n", "state_o", "state_a", "state_c")

trait_map <- c(state_e = "trait_e", state_n = "trait_n", state_o = "trait_o",
               state_a = "trait_a", state_c = "trait_c")

other_constructs_map <- list(
  state_e = c("sociability", "happiness",  "tiredness", "typicality", "negativity"),
  state_n = c("stress", "happiness", "negativity", "duty", "typicality"),
  state_o = c("happiness", "sociability", "tiredness", "duty", "negativity"),
  state_a = c("sociability", "happiness", "tiredness", "hunger", "stress"),
  state_c = c("sociability", "duty", "tiredness", "happiness",  "negativity")
)

cov_matrices_detrend_raw <- list()
cov_matrices_detrend     <- list()
cov_checks_detrend       <- list()
cov_matrices_before      <- list()  # Before matrices (same variable selection)

for (comp in state_components) {
  
  cat("Component:", toupper(comp), "\n")
  
  state_vars <- paste0(comp, "_", 1:32)
  trait_var <- trait_map[comp]
  other_constructs <- other_constructs_map[[comp]]
  other_vars <- c()
  for (oc in other_constructs) {
    other_vars <- c(other_vars, paste0(oc, "_", 1:32))
  }
  
  selected_vars <- c(state_vars, trait_var, other_vars)
  selected_vars <- intersect(selected_vars, names(d2))
  
  cat("  Additional constructs (", length(other_constructs), "):",
      paste(other_constructs, collapse = ", "), "\n")
  cat("  Total variables:", length(selected_vars),
      "(32 State + 1 Trait +", length(other_vars), "additional)\n")
  
  cov_matrices_detrend_raw[[comp]] <- list()
  cov_matrices_detrend[[comp]] <- list()
  cov_checks_detrend[[comp]] <- list()
  cov_matrices_before[[comp]] <- list()
  
  for (m_idx in 1:imp$m) {
    
    # Before covariance matrix (same variable selection, not detrended)
    d_imp <- complete(imp, action = m_idx)
    cov_matrices_before[[comp]][[m_idx]] <- cov(d_imp[, selected_vars], use = "everything")
    
    # After covariance matrix (detrended)
    cov_raw <- cov(imp_detrend_list[[m_idx]][, selected_vars], use = "everything")
    cov_matrices_detrend_raw[[comp]][[m_idx]] <- cov_raw
    
    # 5. nearPD if necessary
    eigs_raw <- eigen(cov_raw, only.values = TRUE)$values
    nearPD_applied <- any(eigs_raw <= 1e-10)
    nearPD_max_change <- 0
    
    if (nearPD_applied) {
      cov_fixed <- as.matrix(nearPD(cov_raw, corr = FALSE, keepDiag = TRUE,
                                    ensureSymmetry = TRUE)$mat)
      dimnames(cov_fixed) <- dimnames(cov_raw)
      nearPD_max_change <- max(abs(cov_fixed - cov_raw))
      cov_matrices_detrend[[comp]][[m_idx]] <- cov_fixed
      cat("  m =", m_idx, ": nearPD applied (max. change:", sprintf("%.2e", nearPD_max_change), ")\n")
    } else {
      cov_matrices_detrend[[comp]][[m_idx]] <- cov_raw
    }
    
    # 6. CHECK (after nearPD)
    cov_mat <- cov_matrices_detrend[[comp]][[m_idx]]
    eigenvals <- eigen(cov_mat, only.values = TRUE)$values
    
    cov_checks_detrend[[comp]][[m_idx]] <- list(
      is_positive_definite = all(eigenvals > 1e-10),
      min_eigenvalue = min(eigenvals),
      max_eigenvalue = max(eigenvals),
      condition_number = max(eigenvals) / max(min(eigenvals[eigenvals > 1e-10]), 1e-15),
      determinant = det(cov_mat),
      n_negative_eigenvalues = sum(eigenvals < 0),
      nearPD_applied = nearPD_applied,
      nearPD_max_change = nearPD_max_change
    )
  }
  
  # Summary
  pos_def_count <- sum(sapply(cov_checks_detrend[[comp]], function(x) x$is_positive_definite))
  cond_nums <- sapply(cov_checks_detrend[[comp]], function(x) x$condition_number)
  nearPD_count <- sum(sapply(cov_checks_detrend[[comp]], function(x) x$nearPD_applied))
  max_changes <- sapply(cov_checks_detrend[[comp]], function(x) x$nearPD_max_change)
  
  cat("  Positive Definite (after nearPD):", pos_def_count, "/", imp$m, "\n")
  cat("  nearPD applied:", nearPD_count, "/", imp$m, "matrices\n")
  if (nearPD_count > 0) {
    cat("  Max. change from nearPD:", sprintf("%.2e", max(max_changes)), "\n")
  }
  cat("  Condition Number - Min:", sprintf("%.2e", min(cond_nums)),
      ", Max:", sprintf("%.2e", max(cond_nums)), "\n\n")
}

# ---- 7. Comparison: variances/covariances before vs after -------------------
cat("=== COMPARISON BEFORE vs AFTER ===\n\n")

pdf(file.path( folder.main, "04_detrend_comparison.pdf" ), width = 16, height = 12)

for (comp in state_components) {
  
  avg_cov_before <- Reduce("+", cov_matrices_before[[comp]]) / imp$m
  avg_cov_after  <- Reduce("+", cov_matrices_detrend[[comp]]) / imp$m
  
  var_before <- diag(avg_cov_before)
  var_after  <- diag(avg_cov_after)
  
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  
  plot(var_before, var_after,
       xlab = "Variance (before detrending)", ylab = "Variance (after detrending)",
       main = "Variances", pch = 16, col = rgb(0, 0, 0, 0.3))
  abline(0, 1, col = "red", lwd = 2)
  
  var_reduction <- 100 * (1 - var_after / var_before)
  hist(var_reduction, breaks = 30,
       main = "Variance Reduction (%)",
       xlab = "% Reduction", col = "steelblue")
  abline(v = median(var_reduction, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
  
  cov_before_offdiag <- avg_cov_before[lower.tri(avg_cov_before)]
  cov_after_offdiag  <- avg_cov_after[lower.tri(avg_cov_after)]
  
  plot(cov_before_offdiag, cov_after_offdiag,
       xlab = "Covariance (before detrending)", ylab = "Covariance (after detrending)",
       main = "Covariances", pch = ".", col = rgb(0, 0, 0, 0.2))
  abline(0, 1, col = "red", lwd = 2)
  
  avg_cor_before <- cov2cor(avg_cov_before)
  avg_cor_after  <- cov2cor(avg_cov_after)
  
  cor_before_offdiag <- avg_cor_before[lower.tri(avg_cor_before)]
  cor_after_offdiag  <- avg_cor_after[lower.tri(avg_cor_after)]
  
  plot(cor_before_offdiag, cor_after_offdiag,
       xlab = "Correlation (before detrending)", ylab = "Correlation (after detrending)",
       main = "Correlations", pch = ".", col = rgb(0, 0, 0, 0.2))
  abline(0, 1, col = "red", lwd = 2)
  
  mtext(toupper(comp), outer = TRUE, cex = 1.5)
  
  cat(toupper(comp), ":\n")
  cat("  Median variance reduction:", sprintf("%.1f%%", median(var_reduction, na.rm = TRUE)), "\n")
  cat("  Median |correlation| before:", sprintf("%.3f", median(abs(cor_before_offdiag), na.rm = TRUE)), "\n")
  cat("  Median |correlation| after:", sprintf("%.3f", median(abs(cor_after_offdiag), na.rm = TRUE)), "\n\n")
}

dev.off()

cat(" -> 04_detrend_comparison.pdf\n\n")

# ---- 8. Save -----------------------------------------------------------------
cat("Saving...\n")

save(imp_detrend_list, file = file.path( folder.main, "imp_detrend_list.Rdata" ))
save(cov_matrices_before, file = file.path( folder.main, "cov_matrices_before.Rdata" ))
save(cov_matrices_detrend, file = file.path( folder.main, "cov_matrices_detrend.Rdata" ))
save(cov_matrices_detrend_raw, file = file.path( folder.main, "cov_matrices_detrend_raw.Rdata" ))
save(cov_checks_detrend, file = file.path( folder.main, "cov_checks_detrend.Rdata" ))

cat(" -> imp_detrend_list.Rdata (all detrended datasets)\n")
cat(" -> cov_matrices_before.Rdata (covariance matrices BEFORE detrending)\n")
cat(" -> cov_matrices_detrend.Rdata (covariance matrices after nearPD)\n")
cat(" -> cov_matrices_detrend_raw.Rdata (covariance matrices WITHOUT nearPD)\n")
cat(" -> cov_checks_detrend.Rdata (checks)\n\n")

# ---- 9. Report ---------------------------------------------------------------
sink(file.path( folder.main, "05_detrend_report.txt" ))

cat("=== DETRENDING REPORT ===\n\n")
cat("Method:", detrend_method, "\n")
if (detrend_method == "loess") cat("LOESS span:", loess_span, "\n")
cat("Datasets (m):", imp$m, "\n\n")

for (comp in state_components) {
  cat("---", toupper(comp), "---\n")
  
  pos_def <- sum(sapply(cov_checks_detrend[[comp]], function(x) x$is_positive_definite))
  cond_nums <- sapply(cov_checks_detrend[[comp]], function(x) x$condition_number)
  min_eigs <- sapply(cov_checks_detrend[[comp]], function(x) x$min_eigenvalue)
  nearPD_count <- sum(sapply(cov_checks_detrend[[comp]], function(x) x$nearPD_applied))
  max_changes <- sapply(cov_checks_detrend[[comp]], function(x) x$nearPD_max_change)
  
  cat("Positive Definite (after nearPD):", pos_def, "/", imp$m, "\n")
  cat("nearPD applied:", nearPD_count, "/", imp$m, "\n")
  if (nearPD_count > 0) {
    cat("Max. change from nearPD:", sprintf("%.2e", max(max_changes)), "\n")
  }
  cat("Condition Number (Min/Max):", sprintf("%.2e", min(cond_nums)), "/",
      sprintf("%.2e", max(cond_nums)), "\n")
  cat("Min Eigenvalue (Min/Max):", sprintf("%.2e", min(min_eigs)), "/",
      sprintf("%.2e", max(min_eigs)), "\n")
  
  avg_before <- Reduce("+", cov_matrices_before[[comp]]) / imp$m
  avg_after  <- Reduce("+", cov_matrices_detrend[[comp]]) / imp$m
  var_red <- 100 * (1 - diag(avg_after) / diag(avg_before))
  cat("Variance reduction: Median", sprintf("%.1f%%", median(var_red, na.rm = TRUE)),
      ", Range", sprintf("%.1f%% - %.1f%%", min(var_red, na.rm = TRUE), max(var_red, na.rm = TRUE)), "\n\n")
}

sink()

cat(" -> 05_detrend_report.txt\n\n")
cat("=== DONE ===\n")