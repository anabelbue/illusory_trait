################################################################################
#  DIAGNOSTICS: Convergence + Covariance Matrix Checks
#  Input: imp.Rdata (from 01_imputation.R)
################################################################################


library(mice)

# ============================================================
# USER CONFIGURATION - adjust this path to your setup
# ============================================================
base_path <- ""  # your workspace/repository path
# ============================================================
# ---- 1. Load -----------------------------------------------------------------
sim <- "02"

folder.main.pre <- file.path(base_path, "output", sim)
folder.main     <- file.path(folder.main.pre)
if (!dir.exists(folder.main)) dir.create(folder.main, recursive = TRUE)

d2 <- read.csv(file.path(base_path, "data", "wide_dat_S2.csv"))

load( file.path( folder.main, "imp.Rdata" ) )

cat("MICE object loaded: m =", imp$m, " datasets\n\n")

# ---- 2. Trace Plots (ALL variables) -----------------------------------------
cat("Generating trace plots for ALL variables...\n")

pdf(file.path( folder.main, "01_trace_plots_all.pdf" ), width = 14, height = 12)
plot(imp)
dev.off()

cat(" -> 01_trace_plots_all.pdf\n")
cat("   (Check whether chains have converged)\n\n")

# ---- 3. Density Plots with dataset numbers ----------------------------------
cat("Generating density plots with dataset labeling...\n")

pdf(file.path( folder.main, "02_density_plots_labeled.pdf" ), width = 14, height = 10)

vars_with_na <- names(d2)[colSums(is.na(d2)) > 0]
colors <- rainbow(imp$m)

for (v in vars_with_na) {
  tryCatch({
    
    # Pre-compute all densities
    all_dens <- list()
    for (m_idx in 1:imp$m) {
      d_imp <- complete(imp, action = m_idx)
      all_dens[[m_idx]] <- density(d_imp[[v]], na.rm = TRUE)
    }
    dens_true <- density(d2[[v]], na.rm = TRUE)
    
    # Compute X and Y ranges from ALL densities
    xmin <- min(c(sapply(all_dens, function(d) min(d$x)), min(dens_true$x)))
    xmax <- max(c(sapply(all_dens, function(d) max(d$x)), max(dens_true$x)))
    ymax <- max(c(sapply(all_dens, function(d) max(d$y)), max(dens_true$y)))
    
    # Plot setup
    plot(NULL, xlim = c(xmin, xmax),
         ylim = c(0, ymax * 1.1),
         main = paste("Observed vs Imputed:", v),
         xlab = v, ylab = "Density")
    
    # Draw imputed datasets (with numbers)
    for (m_idx in 1:imp$m) {
      dens <- all_dens[[m_idx]]
      lines(dens, col = colors[m_idx], lwd = 1)
      peak_idx <- which.max(dens$y)
      text(dens$x[peak_idx], dens$y[peak_idx], m_idx,
           col = colors[m_idx], cex = 0.7, font = 2)
    }
    
    # Observed data LAST (black, thick)
    lines(dens_true, col = "black", lwd = 4)
    
  }, error = function(e) {
    cat("Skipped:", v, "-", e$message, "\n")
  })
}

dev.off()

cat(" -> 02_density_plots_labeled.pdf\n")


# ---- 4. Covariance matrices: EXTRACT, CHECK, SAVE ---------------------------
cat("=== COVARIANCE MATRIX ANALYSIS ===\n\n")

state_components <- c("state_e", "state_n", "state_o", "state_a", "state_c")

trait_map <- c(state_e = "trait_e", state_n = "trait_n", state_o = "trait_o",
               state_a = "trait_a", state_c = "trait_c")

# Exactly 5 additional constructs per state component.
# This list is intentionally editable per state.
# Currently the same 5 constructs are entered for all;
# this can be adjusted per component later.
other_constructs_map <- list(
  state_e = c("sociability", "happiness",  "tiredness", "typicality", "negativity"),
  state_n = c("stress", "happiness", "negativity", "duty", "typicality"),
  state_o = c("happiness", "sociability", "tiredness", "duty", "negativity"),
  state_a = c("sociability", "happiness", "tiredness", "hunger", "stress"),
  state_c = c("sociability", "duty", "tiredness", "happiness",  "negativity")
)

# Storage for all covariance matrices and checks
cov_matrices <- list()
cov_checks <- list()

for (comp in state_components) {
  
  cat("Component:", toupper(comp), "\n")
  
  # Define variables
  state_vars <- paste0(comp, "_", 1:32)
  trait_var <- trait_map[comp]
  other_constructs <- other_constructs_map[[comp]]
  other_vars <- c()
  for (oc in other_constructs) {
    other_vars <- c(other_vars, paste0(oc, "_", 1:32))
  }
  
  selected_vars <- c(state_vars, trait_var, other_vars)
  selected_vars <- intersect(selected_vars, names(d2))
  
  cat("  Additional constructs:", paste(other_constructs, collapse = ", "), "\n")
  cat("  ", length(selected_vars), "variables (State + Trait + 5 additional constructs)\n")
  
  # ---- 4a. EXTRACT: Covariance matrices from all m datasets ----
  cov_matrices[[comp]] <- list()
  
  for (m_idx in 1:imp$m) {
    d_imp <- complete(imp, action = m_idx)
    cov_mat <- cov(d_imp[, selected_vars], use = "everything")
    cov_matrices[[comp]][[m_idx]] <- cov_mat
  }
  
  cat("  Extracted:", imp$m, "covariance matrices\n")
  
  # ---- 4b. CHECK: Positive definiteness + eigenvalues ----
  cov_checks[[comp]] <- list()
  
  for (m_idx in 1:imp$m) {
    cov_mat <- cov_matrices[[comp]][[m_idx]]
    eigenvals <- eigen(cov_mat, only.values = TRUE)$values
    
    is_pos_def <- all(eigenvals > 1e-10)
    cond_num <- max(eigenvals) / max(min(eigenvals[eigenvals > 1e-10]), 1e-15)
    
    cov_checks[[comp]][[m_idx]] <- list(
      is_positive_definite = is_pos_def,
      min_eigenvalue = min(eigenvals),
      max_eigenvalue = max(eigenvals),
      condition_number = cond_num,
      determinant = det(cov_mat),
      n_negative_eigenvalues = sum(eigenvals < 0)
    )
  }
  
  # Summary
  pos_def_count <- sum(sapply(cov_checks[[comp]], function(x) x$is_positive_definite))
  cat("  \u2713 Positive Definite:", pos_def_count, "/", imp$m, "\n")
  
  cond_nums <- sapply(cov_checks[[comp]], function(x) x$condition_number)
  cat("  Condition Number - Min:", sprintf("%.2e", min(cond_nums)), 
      ", Max:", sprintf("%.2e", max(cond_nums)), "\n\n")
}

# ---- 4c. SAVE: Covariance matrices + checks ----
cat("Saving covariance matrices and checks...\n")

save(cov_matrices, file = file.path( folder.main, "cov_matrices.Rdata" ) )
save(cov_checks, file = file.path( folder.main, "cov_checks.Rdata" ) )

cat(" -> cov_matrices.Rdata (all covariance matrices for all m datasets)\n")
cat(" -> cov_checks.Rdata (positive definiteness, eigenvalues, etc.)\n\n")

# ---- 5. Report ---------------------------------------------------------------
sink(file.path( folder.main, "03_diagnostics_report.txt" ) )

cat("=== MICE DIAGNOSTICS REPORT ===\n\n")
cat("Datasets (m):", imp$m, "\n")
cat("Iterations (maxit):", max(imp$iteration), "\n\n")

for (comp in state_components) {
  cat("---", toupper(comp), "---\n")
  
  pos_def_count <- sum(sapply(cov_checks[[comp]], function(x) x$is_positive_definite))
  cond_nums <- sapply(cov_checks[[comp]], function(x) x$condition_number)
  min_eigs <- sapply(cov_checks[[comp]], function(x) x$min_eigenvalue)
  
  cat("Positive Definite:", pos_def_count, "/", imp$m, "\n")
  cat("Condition Number (Min/Max):", sprintf("%.2e", min(cond_nums)), "/", 
      sprintf("%.2e", max(cond_nums)), "\n")
  cat("Min Eigenvalue (Min/Max):", sprintf("%.2e", min(min_eigs)), "/", 
      sprintf("%.2e", max(min_eigs)), "\n\n")
}

sink()

cat(" -> 03_diagnostics_report.txt\n\n")

cat("=== DONE ===\n")
cat("Saved for further analyses:\n")
cat("  - cov_matrices.Rdata\n")
cat("  - cov_checks.Rdata\n")