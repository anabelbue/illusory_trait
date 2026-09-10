library(mice)
library(dplyr)

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


cat("Dimensions:", nrow(d2), "rows x", ncol(d2), "columns\n")
cat("Overall proportion missing:", round(mean(is.na(d2)) * 100, 1), "%\n\n")

# Inspect missing pattern per construct
constructs <- c("state_e", "state_n", "state_o", "state_c", "state_a", 
                "duty", "sociability", "negativity", "typicality", 
                "happiness", "stress", "hunger", "tiredness")

cat("Proportion missing per construct:\n")
for (cst in constructs) {
  cols <- grep(paste0("^", cst, "_\\d+$"), names(d2), value = TRUE)
  pct  <- round(mean(is.na(d2[, cols])) * 100, 1)
  cat(sprintf("  %-15s: %5.1f%% missing\n", cst, pct))
}

# Check trait variables
cat("\nMissing in trait variables:\n")
traits_check <- c("trait_e", "trait_n", "trait_o", "trait_a", "trait_c")
for (tr in traits_check) {
  n_miss <- sum(is.na(d2[[tr]]))
  cat(sprintf("  %-10s: %d missing (%4.1f%%)\n", tr, n_miss,
              100 * n_miss / nrow(d2)))
}

# ---- 2. Strategy: Custom predictor matrix -----------------------------------
#

# Variable names
varnames   <- names(d2)
n_vars     <- length(varnames)
traits     <- c("trait_e", "trait_n", "trait_o", "trait_a", "trait_c")
window     <- 3   # +/- 3 time points as predictors

# Helper function: extract construct and time point from variable name
parse_varname <- function(v) {
  m <- regmatches(v, regexec("^(.+)_(\\d+)$", v))[[1]]
  if (length(m) == 3) return(list(construct = m[2], time = as.integer(m[3])))
  return(NULL)
}

# ---- 3. Build predictor matrix ----------------------------------------------
# Rows = variables to impute, columns = predictors
# 1 = use as predictor, 0 = do not use

pred <- matrix(0, nrow = n_vars, ncol = n_vars,
               dimnames = list(varnames, varnames))

for (i in seq_along(varnames)) {
  vi <- varnames[i]
  
  # id never as predictor, id never imputed
  if (vi == "id") next
  
  # Traits: other traits + some corresponding state time points as predictors
  if (vi %in% traits) {
    # Other traits
    pred[vi, setdiff(traits, vi)] <- 1
    # Corresponding state variables (5 evenly distributed time points)
    # e.g. trait_e -> state_e_1, state_e_9, state_e_17, state_e_24, state_e_32
    suffix <- sub("trait_", "state_", vi)
    state_preds <- paste0(suffix, "_", round(seq(1, 32, length.out = 5)))
    state_preds <- intersect(state_preds, varnames)
    pred[vi, state_preds] <- 1
    next
  }
  
  # Time-varying: parse
  parsed <- parse_varname(vi)
  if (is.null(parsed)) next
  
  cst <- parsed$construct
  tp  <- parsed$time
  
  # (a) All traits as predictors
  pred[vi, traits] <- 1
  
  # (b) Neighboring time points of the same construct (window)
  neighbor_tps <- max(1, tp - window):min(32, tp + window)
  neighbor_tps <- setdiff(neighbor_tps, tp)  # exclude itself
  neighbor_vars <- paste0(cst, "_", neighbor_tps)
  neighbor_vars <- intersect(neighbor_vars, varnames)
  pred[vi, neighbor_vars] <- 1
  
  # (c) Window +/- 1 of other constructs (cross-lagged structure)
  for (other_cst in setdiff(constructs, cst)) {
    for (lag in -1:1) {
      tp_lag <- tp + lag
      if (tp_lag >= 1 && tp_lag <= 32) {
        lag_var <- paste0(other_cst, "_", tp_lag)
        if (lag_var %in% varnames) {
          pred[vi, lag_var] <- 1
        }
      }
    }
  }
}

# Diagonal to 0 (variable cannot predict itself)
diag(pred) <- 0
# id never as predictor
pred[, "id"] <- 0

# Check: max. number of predictors per variable
pred_counts <- rowSums(pred)
cat("\nPredictors per variable:\n")
cat("  Min:", min(pred_counts[pred_counts > 0]), "\n")
cat("  Max:", max(pred_counts), "\n")
cat("  Median:", median(pred_counts[pred_counts > 0]), "\n")

# ---- 4. Choose imputation method --------------------------------------------
#
# Options when p >> n:
#   "pmm"  = Predictive Mean Matching (default, robust, preserves distribution)
#            -> works WITH the reduced predictor matrix
#   "cart" = Regression trees (alternative, does not need reduced matrix)
#   "rf"   = Random Forest (good but slow)
#
# Recommendation: pmm with the reduced predictor matrix above
#

meth <- make.method(d2)
meth["id"] <- ""   # do not impute id

# All set to PMM (default for numeric variables)
# If you prefer CART, change the next line:
# meth[meth != ""] <- "cart"

cat("\nMethods:\n")
print(table(meth))

# ---- 5. Run mice ------------------------------------------------------------
cat("\nStarting imputation (this may take a few minutes)...\n")
t_start <- Sys.time()

set.seed(1234)

imp <- mice(
  data            = d2,
  m               = 20,          # 20 imputed datasets (standard)
  method          = meth,
  predictorMatrix = pred,
  maxit           = 60,          # 60 iterations
  ridge           = 1e-04,       # Ridge penalty (default: 1e-05), stabilizes
  printFlag       = TRUE         # show progress
)

t_end <- Sys.time()
cat("\nRuntime:", round(difftime(t_end, t_start, units = "mins"), 1), "minutes\n")

# ---- Save -------------------------------------------------------------------
save(imp, d2, file = file.path(folder.main, paste0("imp.Rdata")) )