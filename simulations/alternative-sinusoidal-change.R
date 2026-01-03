# Clear the workspace
rm(list = ls())
# Set seed for reproducibility
set.seed(123)
##########################################################################################
run_sim_alt <- function(rep_id, m, T.chan, gamma, delta = 1, s.star = 200, alpha = 0.05, dgp_type = "BB") {
  total_n <- m + m * T.chan
  
  # Data generation
  baseInfo <- list(
    n = total_n,
    basisType = "bspline",
    muInfo = list(type = "sin", a = 1),
    factor = 1
  )
  
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
  
  # Inject smooth sinusoidal change after t* = m + s.star
  t_star <- m + s.star
  change_idx <- t_star:total_n
  num_change <- length(change_idx)
  
  # Define grid and basis
 
  t_grid <- seq(0, 1, length.out = 301)
  basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = 21)
  
  # For each time point after change, add sinusoidal drift
  for (i in seq_along(change_idx)) {
    rel_t <- (i - 1) / (m * T.chan)
    drift_val <- delta * sin(pi * rel_t)
    drift_vec <- rep(drift_val, length(t_grid))  # length 301
    
    # Convert to fd object (argvals = t_grid, not list)
    mu_fd <- Data2fd(argvals = t_grid, y = matrix(drift_vec, nrow = 301, ncol = 1), basisobj = basis)
    
    # Inject into fd_data
    fd_data$coefs[, change_idx[i]] <- fd_data$coefs[, change_idx[i]] + mu_fd$coefs[, 1]
  }
  
  
  # FPCA decomposition
  t_grid <- seq(0, 1, length.out = 301)
  values_all <- t(eval.fd(t_grid, fd_data))
  values_train <- values_all[1:m, ]
  values_monitor <- values_all[(m+1):total_n, ]
  
  Ly_train <- split(values_train, row(values_train))
  Lt_train <- replicate(length(Ly_train), t_grid, simplify = FALSE)
  
  fpca_train <- FPCA(
    Ly = Ly_train, Lt = Lt_train,
    optns = list(methodMuCovEst = 'smooth', FVEthreshold = 0.95, methodSelectK = 'FVE')
  )
  
  mean_train <- fpca_train$mu
  phi_train <- fpca_train$phi
  scores_train <- fpca_train$xiEst
  
  # Monitoring scores
  centered_monitor <- sweep(values_monitor, 2, mean_train)
  scores_monitor <- centered_monitor %*% phi_train * (1 / length(t_grid))
  scores_all <- rbind(scores_train, scores_monitor)
  
  # Prewhitening
  sample.omega <- var(scores_train)
  A_estimate <- ldl(sample.omega)$lower
  scores_whitened <- t(ginv(A_estimate) %*% t(scores_all))
  scores_whitened <- scores_whitened[, 1:min(ncol(scores_whitened), 8), drop = FALSE]
  
  # Run monitoring statistics
  ssms <- ssms.statistic.fpca.alt(scores_all, m = m, T.chan = T.chan, alpha = alpha, gamma = gamma)
  rsms <- rsms.statistic.fpca.alt(scores_whitened, m = m, T.chan = T.chan, gamm = gamma, alpha = alpha)
  csms <- csms.statistic.fpca.alt(scores_all, m = m, T.chan = T.chan, gamm = gamma, alpha = alpha)
  
  return(data.frame(
    rep = rep_id,
    stat_rsms = rsms$statistic, reject_rsms = rsms$reject, first_rsms = rsms$first_rejection,
    stat_ssms = ssms$statistic, reject_ssms = ssms$reject, first_ssms = ssms$first_rejection,
    stat_csms = csms$statistic, reject_csms = csms$reject, first_csms = csms$first_rejection
  ))
}

################################################################################ 
# Set seed
set.seed(123)

# Set working directory (optional for RStudio users)
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Load required files and packages
source("genData.R")
source("Method.R")
library(pbmcapply)

# Parallel cores
ncores <- parallel::detectCores()-3

# Parameters
 m <- 500
 T_vals <- c(1, 2)
 gamma_vals <- c(0, 0.15)
 
 delta_vals <- c(0.10,0.15,0.20,0.25,0.30,0.35)

 s_star_vals <- c(50, 200)

 

nsim <- 1000
dgp_types <- c("BB", "fIID", "fMA1")

# Create output directory
output_dir <- "output_sinusoidal_change"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Full grid
param_grid <- expand.grid(
  T.chan = T_vals,
  gamma = gamma_vals,
  delta = delta_vals,
  s.star = s_star_vals,
  dgp_type = dgp_types
)

# Run simulation loop
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ] 
  T.chan <- params$T.chan
  gamma <- params$gamma
  delta <- params$delta
  s.star <- params$s.star
  dgp_type <- params$dgp_type
  
  cat(sprintf("\nRunning: DGP = %s, m = %d, T = %d, gamma = %.2f, delta = %.2f, s.star = %d\n",
              dgp_type, m, T.chan, gamma, delta, s.star))
  
  results <- pbmclapply(
    1:nsim,
    function(rep_id) {
      run_sim_alt(rep_id = rep_id, m = m, T.chan = T.chan, gamma = gamma,
                  delta = delta, s.star = s.star, alpha = 0.05, dgp_type = dgp_type)
    },
    mc.cores = ncores
  )
  
  results_df <- do.call(rbind, results)
  
  filename <- sprintf("monitoring_results_%s_m%d_T%d_gamma%s_delta%s_sstar%d_nsim%d.csv",
                      dgp_type, m, T.chan,
                      gsub("\\.", "p", as.character(gamma)),
                      gsub("\\.", "p", as.character(delta)),
                      s.star, nsim)
  
  write.csv(results_df, file = file.path(output_dir, filename), row.names = FALSE)
}
