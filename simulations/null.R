# Clear the workspace
rm(list = ls())
# Set seed for reproducibility
set.seed(123)
##########################################################################################

 
run_sim <- function(rep_id, m, T.chan, gamma, alpha = 0.05 , dgp_type = "BB") {
  total_n <- m + m * T.chan
  # Define common base
  baseInfo <- list(
    n = total_n,
    basisType = "bspline",
    muInfo = list(type = "sin", a = 1),
    factor = 1
  )
  
  # Extend with DGP-specific fields
  if (dgp_type == "BB") {
    dataInfo <- c(baseInfo, list(nArgvals = 300))
    fd_data <- BB(dataInfo)
  } else if (dgp_type == "fIID") {
    dataInfo <- c(baseInfo, list(varType = "A", gaussian = TRUE))
    fd_data <- fIID(dataInfo)
  } else if (dgp_type == "fMA1") {
    dataInfo <- c(baseInfo, list(varType = "A", kappa = 0.7))
    fd_data <- fMA1(dataInfo)
  } else {
    stop(sprintf("Unsupported DGP type: %s", dgp_type))
  }
  
  # Ensure it's an fd object
  if (!inherits(fd_data, "fd")) {
    stop(sprintf("The object returned by %s is not of class 'fd'", dgp_type))
  }
  
  
  
  # Ensure the object is a functional data object
  if (!inherits(fd_data, "fd")) {
    stop(sprintf("The object returned by %s is not of class 'fd'", dgp_type))
  }
  
  # FPCA
  t_grid <- seq(0, 1, length.out = 301)
  values_all <- t(eval.fd(t_grid, fd_data))
  values_train <- values_all[1:m, ]
  values_monitor <- values_all[(m+1):(m + m * T.chan), ]
  Ly_train <- split(values_train, row(values_train))
  Lt_train <- replicate(length(Ly_train), t_grid, simplify = FALSE)
  fpca_train <- FPCA(Ly = Ly_train, Lt = Lt_train, optns = list(methodMuCovEst = 'smooth', FVEthreshold = 0.95, methodSelectK = 'FVE'))
  mean_train <- fpca_train$mu
  phi_train <- fpca_train$phi
  K <- ncol(phi_train)
  
  centered_monitor <- sweep(values_monitor, 2, mean_train)
  scores_monitor <- centered_monitor %*% phi_train * (1 / length(t_grid))
  scores_train <- fpca_train$xiEst
  scores_all <- rbind(scores_train, scores_monitor)
  
  # Prewhitening
  sample.omega <- var(scores_train)
  A_estimate <- ldl(sample.omega)$lower
  scores_whitened <- t(ginv(A_estimate) %*% t(scores_all))
  K <- min(ncol(scores_whitened), 8)
  scores_whitened <- scores_whitened[, 1:K, drop = FALSE]
  scores <- scores_whitened  # 🔧 Fix for CSMS input
  
  # Run tests
  ssms <- ssms.statistic.fpca(scores_all, m = m, T.chan = T.chan, alpha = alpha, gamma = gamma)
  rsms <- rsms.statistic.fpca(scores_whitened, m = m, T.chan = T.chan, gamm = gamma, alpha = alpha)
  csms <- csms.statistic.fpca(scores_all, m = m, T.chan = T.chan, gamm = gamma, alpha = alpha)
  
  return(data.frame(
    rep = rep_id,
    stat_rsms = rsms$statistic, reject_rsms = rsms$reject,
    stat_ssms = ssms$statistic, reject_ssms = ssms$reject,
    stat_csms = csms$statistic, reject_csms = csms$reject
  ))
}

################################################################################ 
#  set seed
 
set.seed(123)

# Set working directory if in RStudio
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Load required files and packages
source("genData.R")
source("Method.R")
library(pbmcapply)

# Detect all available cores
ncores <- parallel::detectCores() -3 

# Parameter options
m_vals <- c(100, 200, 500, 1000)
T_vals <- c(1, 2, 5, 10)
gamma_vals <- c(0, 0.15)
nsim <- 1000 
dgp_types <- c("BB", "fIID", "fMA1")

# Generate all combinations
param_grid <- expand.grid(m = m_vals, T.chan = T_vals, gamma = gamma_vals, dgp_type = dgp_types)

# Create output directory if missing
output_dir <- "output_null"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Main loop
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  m <- params$m
  T.chan <- params$T.chan
  gamma <- params$gamma
  dgp_type <- params$dgp_type
  
  cat(sprintf("\nRunning simulation: DGP = %s, m = %d, T = %d, gamma = %.2f\n", dgp_type, m, T.chan, gamma))
  
  # Parallel simulation with progress bar
  results <- pbmclapply(
    1:nsim,
    function(rep_id) {
      run_sim(rep_id = rep_id, m = m, T.chan = T.chan, gamma = gamma, alpha = 0.05, dgp_type = dgp_type)
    },
    mc.cores = ncores
  )
  
  # Combine results
  results_df <- do.call(rbind, results)
  
  # Save to file
  filename <- sprintf("monitoring_results_%s_m%d_T%d_gamma%s_nsim%d.csv",
                      dgp_type, m, T.chan, gsub("\\.", "p", as.character(gamma)), nsim)
  write.csv(results_df, file = file.path(output_dir, filename), row.names = FALSE)
}
