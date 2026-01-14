
required_packages <- c("fdapace", "parallel", "fda", "zoo", "sandwich", "MASS", "fastmatrix")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

library(fdapace)
library(parallel)
library(fda)
library(zoo)
library(MASS)
library(fastmatrix)


################################################################################
# Zhu et al. (2025) — Adjusted-Range Based KS Statistic
# Includes critical values for gamm ∈ {0, 0.15}, alpha ∈ {0.05, 0.10}

# ============================================================================== 
# Critical Value Lookup Table
# ============================================================================== 
get_critical_value_rsms <- function(T.chan, alpha = 0.05, gamm = 0, d) {
  # Validate input
  if (!T.chan %in% c(1, 2, 5, 10)) stop("T.chan must be one of 1, 2, 5, or 10.")
  if (!alpha %in% c(0.05, 0.1, 0.10)) stop("alpha must be 0.05 or 0.10.")
  if (!gamm %in% c(0, 0.15)) stop("gamm must be 0 or 0.15.")
  if (d < 1 || d > 8) stop("Only dimensions d = 1 to 8 are supported.")
  
  # Format alpha to fixed-point character
  alpha_key <- formatC(alpha, format = "f", digits = 2)
  t_key <- as.character(T.chan)
  
  # Column map: alpha → T.chan index
  col_map <- list(
    "0.05" = c(`1` = 1, `2` = 3, `5` = 5, `10` = 7),
    "0.10" = c(`1` = 2, `2` = 4, `5` = 6, `10` = 8)
  )
  
  # Validate mapping
  if (!alpha_key %in% names(col_map)) stop("Alpha key not found.")
  if (!t_key %in% names(col_map[[alpha_key]])) stop("T.chan key not found for this alpha.")
  
  # Retrieve index
  col_idx <- col_map[[alpha_key]][[t_key]]
  offset <- ifelse(gamm == 0, 0, 8)  # gamm = 0.15 is offset by 8 columns
  
  # Critical value matrix: rows = d (1:8), columns = 8 for gamm=0, 8 for gamm=0.15
  crit_mat <- matrix(c(
    # gamm=0 (columns 1-8)
    2.1, 1.5, 2.7, 2, 3.4, 2.5, 3.9, 2.8,      # d=1
    3.2, 2.4, 5.2, 4.1, 10.9, 8.3, 21.1, 15.9,  # d=2
    4, 3.3, 6.5, 5.1, 13.8, 10.8, 26.6, 21.1,   # d=3
    4.8, 3.9, 7.6, 6.2, 16.2, 12.8, 30.7, 24,   # d=4
    5.5, 4.5, 8.9, 7.3, 18.3, 14.8, 35.3, 28.6, # d=5
    6.1, 5, 9.8, 8.1, 20.8, 17.1, 40.1, 33,     # d=6
    6.7, 5.7, 10.9, 9.1, 22.6, 18.7, 43.1, 36,  # d=7
    7.4, 6.3, 11.9, 10.1, 24.9, 20.9, 46.4, 39.3,# d=8
    # gamm=0.15 (columns 9-16)
    2.7, 2, 3.3, 2.5, 3.9, 2.9, 4.3, 3.2,        # d=1
    4.8, 3.7, 8, 6.2, 19.2, 14.7, 44.6, 33.1,    # d=2
    5.8, 4.7, 9.6, 7.9, 23.8, 18.4, 53.8, 41.7,  # d=3
    6.9, 5.7, 11.5, 9.2, 28.1, 22.2, 63.4, 49.7, # d=4
    7.8, 6.5, 12.7, 10.4, 29.9, 24.5, 71.3, 55.3,# d=5
    8.6, 7.2, 13.8, 11.5, 33.8, 27.9, 77.9, 62.1,# d=6
    9.3, 7.9, 15.3, 12.9, 36.5, 30.3, 83.6, 67,  # d=7
    10.2, 8.6, 16.5, 13.9, 39.4, 32.8, 90.3, 73.3# d=8
  ), nrow = 8, byrow = TRUE)
  
  return(crit_mat[d, col_idx + offset])
}


# ============================================================================== 
# Adjusted-Range KS Test Statistic Function
# ============================================================================== 


rsms.statistic.fpca <- function(input.vec, m, T.chan, gamm = 0, alpha = 0.05) {
  
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  
  # Dimensions
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)  # Number of components (used as 'd')
  
  # Step 1: Demean using training sample
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Cumulative sums
  cumsum_train <- apply(input_demeaned[1:m, ], 2, cumsum)
  cumsum_monitor <- apply(input_demeaned[(m+1):(m + m*T.chan), ], 2, cumsum)
  
  # Step 3: Range-based variance normalization
  range_vals <- apply(cumsum_train, 2, function(x) max(x) - min(x))
  range_vals[range_vals < 1e-8] <- 1e-8  # Avoid division by zero
  V <- m * diag(1 / (range_vals^2))
  
  # Step 4: Compute scaled KS matrix
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  KS_mat <- cumsum_monitor %*% V %*% t(cumsum_monitor)
  KS_scaled <- KS_mat * scaling
  KS_stat <- max(KS_scaled)
  
  # Step 5: Compare to critical value
  crit_val <- get_critical_value_rsms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- KS_stat > crit_val
  
  return(list(statistic = KS_stat, critical_value = crit_val, reject = reject))
}

# ============================================================================== 
rsms.statistic.fpca.alt <- function(input.vec, m, T.chan, gamm = 0, alpha = 0.05) {
  input.vec <- as.matrix(input.vec)  # ensure matrix
  
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Handle vector case separately
  if (K == 1) {
    input_demeaned <- input.vec - mean(input.vec[1:m ])
    
    cumsum_train <- cumsum(input_demeaned[1:m])
    cumsum_monitor <- cumsum(input_demeaned[(m+1):(m + m*T.chan)])
    
    range_val <- max(cumsum_train) - min(cumsum_train)
    if (range_val < 1e-8) range_val <- 1e-8
    V_inv <- m / (range_val^2)
    
    time_idx <- 1:(m * T.chan)
    scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
      ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
    
    KS_scaled <- (cumsum_monitor^2) * V_inv * scaling
    KS_stat <- max(KS_scaled)
    
    crit_val <- get_critical_value_rsms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = 1)
    reject <- KS_stat > crit_val
    first_rejection <- if (reject) which(KS_scaled > crit_val)[1] else length(KS_scaled)
    
    return(list(
      statistic = KS_stat,
      critical_value = crit_val,
      reject = reject,
      first_rejection = first_rejection
    ))
  }
  
  # Multivariate case
  colmeans_m <- colMeans(input.vec[1:m, , drop = FALSE])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  cumsum_train <- apply(input_demeaned[1:m, , drop = FALSE], 2, cumsum)
  cumsum_monitor <- apply(input_demeaned[(m+1):(m + m*T.chan), , drop = FALSE], 2, cumsum)
  
  
  range_vals <- apply(cumsum_train, 2, function(x) max(x) - min(x))
  range_vals[range_vals < 1e-8] <- 1e-8
  V <- m * diag(1 / (range_vals^2))
  
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  KS_mat <- cumsum_monitor %*% V %*% t(cumsum_monitor)
  KS_scaled <- diag(KS_mat) * scaling
  KS_stat <- max(KS_scaled)
  
  crit_val <- get_critical_value_rsms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- KS_stat > crit_val
  first_rejection <- if (reject) which(KS_scaled > crit_val)[1] else length(KS_scaled)
  
  return(list(
    statistic = KS_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}

rsms.statistic.fpca.alt.old.version <- function(input.vec, m, T.chan, gamm = 0, alpha = 0.05) {
  
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  # Dimensions
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training data
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Cumulative sums
  cumsum_train <- apply(input_demeaned[1:m, ], 2, cumsum)
  cumsum_monitor <- apply(input_demeaned[(m+1):(m + m*T.chan), ], 2, cumsum)
  
  # Step 3: Range-based variance normalization
  range_vals <- apply(cumsum_train, 2, function(x) max(x) - min(x))
  range_vals[range_vals < 1e-8] <- 1e-8
  V <- m * diag(1 / (range_vals^2))
  
  # Step 4: Compute test sequence over monitoring period
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  # Apply test
  KS_mat <- cumsum_monitor %*% V %*% t(cumsum_monitor)
  KS_scaled <- diag(KS_mat) * scaling
  KS_stat <- max(KS_scaled)
  
  # Step 5: Rejection rule
  crit_val <- get_critical_value_rsms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- KS_stat > crit_val
  
  # Step 6: First time of rejection
  first_rejection <- if (reject) which(KS_scaled > crit_val)[1] else length(KS_scaled)
  
  return(list(
    statistic = KS_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}


################################################################################
# Chan et al.'s (2020) Self-normalized statistics 
# Includes critical values for gamm ∈ {0, 0.15}, alpha ∈ {0.05, 0.10}


# ============================================================================== 
# Critical Value Lookup Table
# ============================================================================== 


get_critical_value_ssms <- function(T.chan, alpha = 0.05, gamm = 0, d) {
  # Validate input
  if (!T.chan %in% c(1, 2, 5, 10)) stop("T.chan must be one of 1, 2, 5, or 10.")
  if (!alpha %in% c(0.05, 0.1, 0.10)) stop("alpha must be 0.05 or 0.10.")
  if (!gamm %in% c(0, 0.15)) stop("gamm must be 0 or 0.15.")
  if (d < 1 || d > 10) stop("Only dimensions d = 1 to 10 are supported.")
  
  # Format alpha for lookup
  alpha_key <- formatC(alpha, format = "f", digits = 2)
  t_key <- as.character(T.chan)
  
  col_map <- list(
    "0.05" = c(`1` = 1, `2` = 3, `5` = 5, `10` = 7),
    "0.10" = c(`1` = 2, `2` = 4, `5` = 6, `10` = 8)
  )
  
  if (!alpha_key %in% names(col_map)) stop("Alpha key not found.")
  if (!t_key %in% names(col_map[[alpha_key]])) stop("T.chan key not found for this alpha.")
  
  col_idx <- col_map[[alpha_key]][[t_key]]
  offset <- ifelse(gamm == 0, 0, 8)  # gamm=0.15 is offset by 8 columns
  
  # Critical value matrix: rows = d (1:10), columns = 8 for gamm=0, 8 for gamm=0.15
  crit_mat <- matrix(c(
    # gamm=0 (columns 1-8)
    34, 22.8, 44.1, 30, 55.5, 38.1, 59, 40.8,        # d=1
    68.7, 50.4, 93.3, 68.2, 113.4, 84.3, 126.8, 93.1,# d=2
    108.4, 83.7, 153.4, 115, 184.7, 140.5, 206.1, 155.9,# d=3
    165, 127.9, 219.2, 171.3, 265.5, 210.4, 293.8, 229.3,# d=4
    218.5, 175, 284.4, 226.9, 363.4, 287.1, 394.3, 319.1,# d=5
    279.6, 222.3, 375.2, 299.6, 457.9, 371.5, 515, 410.4,# d=6
    344, 282.5, 459.3, 373.5, 580.1, 466.9, 627.2, 509.3,# d=7
    425.5, 341.3, 559.8, 463.6, 694.9, 572.7, 774.4, 628.3,# d=8
    498.9, 414, 660.6, 548.2, 830, 692.2, 929.2, 768.5,   # d=9
    578.6, 486.9, 767.3, 646.8, 967.7, 824.6, 1079.5, 895.2,# d=10
    # gamm=0.15 (columns 9-16)
    43.6, 30.7, 53.1, 36, 60.9, 42.6, 67.8, 47.4,        # d=1
    89.7, 67.5, 111, 83.5, 128.7, 96.8, 140.5, 103.4,    # d=2
    143.9, 112.3, 185.9, 140.5, 204.8, 158.2, 219.6, 171.1, # d=3
    217.9, 171.7, 259.3, 207.9, 297.3, 234.7, 319.5, 252.3, # d=4
    290, 230.4, 345.9, 273.6, 404.5, 317.4, 425.7, 336.7,   # d=5
    358.6, 291.5, 436, 353.1, 508.4, 414, 546.3, 448.3,     # d=6
    445.9, 365.4, 542.4, 447.9, 641, 522.7, 664.9, 550.3,   # d=7
    530.4, 443.1, 647.1, 540.2, 774.7, 637.4, 827.6, 683.2, # d=8
    648.1, 532.1, 777.7, 647.9, 911.8, 768.5, 972.2, 810.6, # d=9
    742.6, 631, 906, 755.3, 1059.8, 898.9, 1136.9, 946.5    # d=10
  ), nrow = 10, byrow = TRUE)
  
  return(crit_mat[d, col_idx + offset])
}



#get_critical_value_ssms(T.chan = 5, alpha = 0.05, gamm = 0.15, d = 4)


#get_critical_value_ssms(T.chan = 2, alpha = 0.10, gamm = 0, d = 6)

# ============================================================================== 
ssms.statistic.fpca  <- function(input.vec, m, T.chan, alpha = 0.05, gamm = 0) { 
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training mean
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Cumulative sum over monitoring period (sn.k)
  sn.k <- apply(input_demeaned[(m+1):(m + m*T.chan), ], 2, cumsum)  # (m*T.chan) x K
  
  # Step 3: Compute denominator matrix D (based on training sample)
  cumsum_train <- apply(input_demeaned[1:m, ], 2, cumsum)
  D <- m^(-2) * t(cumsum_train) %*% cumsum_train
  
  # Step 4: Compute the M-statistic vector
  library(MASS)
  Mt.mat <- sn.k %*% ginv(D) %*% t(sn.k)  # size (m*T.chan) x (m*T.chan)
  
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2  * gamm)
  Mt.vec <- diag(Mt.mat) * scaling
  
  # Step 5: Final statistic
  Mt_stat <- max(Mt.vec)
  
  # Step 6: Compare to critical value from SSMS table
  crit_val <- get_critical_value_ssms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- Mt_stat > crit_val
  
  return(list(statistic = Mt_stat, critical_value = crit_val, reject = reject))
} 

# ==============================================================================  

ssms.statistic.fpca.alt <- function(input.vec, m, T.chan, alpha = 0.05,  gamm = 0) {
  input.vec <- as.matrix(input.vec)
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # ---- Univariate case ----
  if (K == 1) {
    input_demeaned <- input.vec - mean(input.vec[1:m])
    
    # Step 2: cumulative sums of monitoring sample
    sn.k <- cumsum(input_demeaned[(m+1):(m + m*T.chan), ])
    
    # Step 3: Denominator matrix D
    cumsum_train <- cumsum(input_demeaned[1:m])
    D <- m^(-2) * sum(cumsum_train^2)
    
    # Step 4: Compute test sequence
    time_idx <- 1:(m * T.chan)
    scaling <- m^(-1) * (1 + time_idx / m)^(-2 * gamm)
    Mt.vec <- (sn.k^2 / D) * scaling
    
    # Step 5–7: Compute and return
    Mt_stat <- max(Mt.vec)
    crit_val <- get_critical_value_ssms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = 1)
    reject <- Mt_stat > crit_val
    first_rejection <- if (reject) which(Mt.vec > crit_val)[1] else length(Mt.vec)
    
    return(list(
      statistic = Mt_stat,
      critical_value = crit_val,
      reject = reject,
      first_rejection = first_rejection
    ))
  }
  
  # ---- Multivariate case ----
  colmeans_m <- colMeans(input.vec[1:m, , drop = FALSE])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  sn.k <- apply(input_demeaned[(m+1):(m + m*T.chan), , drop = FALSE], 2, cumsum)
  cumsum_train <- apply(input_demeaned[1:m, , drop = FALSE], 2, cumsum)
  
  D <- m^(-2) * t(cumsum_train) %*% cumsum_train
  Mt.mat <- sn.k %*% MASS::ginv(D) %*% t(sn.k)
  
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2 * gamm)
  Mt.vec <- diag(Mt.mat) * scaling
  
  Mt_stat <- max(Mt.vec)
  crit_val <- get_critical_value_ssms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- Mt_stat > crit_val
  first_rejection <- if (reject) which(Mt.vec > crit_val)[1] else length(Mt.vec)
  
  return(list(
    statistic = Mt_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}



ssms.statistic.fpca.alt.old.version <- function(input.vec, m, T.chan, alpha = 0.05,  gamm = 0) { 
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training mean
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Cumulative sums of monitoring sample
  sn.k <- apply(input_demeaned[(m+1):(m + m*T.chan), ], 2, cumsum)
  
  # Step 3: Denominator matrix D
  cumsum_train <- apply(input_demeaned[1:m, ], 2, cumsum)
  D <- m^(-2) * t(cumsum_train) %*% cumsum_train
  
  # Step 4: Compute test sequence
  Mt.mat <- sn.k %*% MASS::ginv(D) %*% t(sn.k)
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2 * gamm)
  Mt.vec <- diag(Mt.mat) * scaling
  
  # Step 5: Max test statistic
  Mt_stat <- max(Mt.vec)
  
  # Step 6: Critical value and rejection
  crit_val <- get_critical_value_ssms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- Mt_stat > crit_val
  
  # Step 7: First rejection time
  first_rejection <- if (reject) which(Mt.vec > crit_val)[1] else length(Mt.vec)
  
  return(list(
    statistic = Mt_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}

# get_critical_value_rsms(T.chan = 5, alpha = 0.05, gamm = 0.15, d = 4)


# get_critical_value_rsms(T.chan = 2, alpha = 0.10, gamm = 0, d = 6)







################################################################################
# Standard HAC methods 
# Includes critical values for gamm ∈ {0, 0.15}, alpha ∈ {0.05, 0.10}
# ==============================================================================
# Critical Value Lookup Table for CSMS  
# ==============================================================================

get_critical_value_csms <- function(T.chan, alpha = 0.05, gamm = 0, d) {
  # Validate input
  if (!T.chan %in% c(1, 2, 5, 10)) stop("T.chan must be one of 1, 2, 5, or 10.")
  if (!alpha %in% c(0.05, 0.1, 0.10)) stop("alpha must be 0.05 or 0.10.")
  if (!gamm %in% c(0, 0.15)) stop("gamm must be 0 or 0.15.")
  if (d < 1 || d > 8) stop("Only dimensions d = 1 to 8 are supported.")
  
  # Normalize alpha and T.chan keys
  alpha_key <- formatC(alpha, format = "f", digits = 2)
  t_key <- as.character(T.chan)
  
  # Map alpha + T.chan → column index
  col_map <- list(
    "0.05" = c(`1` = 1, `2` = 3, `5` = 5, `10` = 7),
    "0.10" = c(`1` = 2, `2` = 4, `5` = 6, `10` = 8)
  )
  
  if (!alpha_key %in% names(col_map)) stop("Alpha key not found.")
  if (!t_key %in% names(col_map[[alpha_key]])) stop("T.chan key not found for this alpha.")
  
  col_idx <- col_map[[alpha_key]][[t_key]]
  offset <- ifelse(gamm == 0, 0, 8)  # gamm=0.15 offset by 8 columns
   
  # Critical value matrix: rows = d (1:8), columns = 8 for gamm=0, 8 for gamm=0.15
  crit_mat <- matrix(c(
    # gamm = 0 (columns 1--8): (T=1: 5%,10%), (T=2: 5%,10%), (T=5: 5%,10%), (T=10: 5%,10%)
    2.2, 1.7, 3.4, 2.7, 3.8, 3.0, 4.9, 3.6,      # d=1
    3.8, 3.1, 5.6, 4.5, 12.0, 9.8, 23.4, 19.4,   # d=2
    4.4, 3.8, 7.5, 6.1, 15.9, 13.5, 28.9, 24.2,  # d=3
    6.0, 5.2, 8.8, 7.4, 18.3, 15.7, 35.7, 29.7,  # d=4
    6.7, 5.6, 10.2, 8.6, 21.0, 18.0, 37.5, 32.6, # d=5
    7.2, 6.2, 11.3, 9.7, 23.6, 20.2, 42.9, 37.2, # d=6
    7.9, 7.1, 12.7, 11.1, 26.8, 23.2, 50.5, 43.9,# d=7
    8.8, 7.9, 13.7, 12.0, 29.7, 25.8, 54.8, 46.3,# d=8
    
    # gamm = 0.15 (columns 9--16) in the same order of (T, alpha)
    3.1, 2.5, 3.8, 3.1, 4.5, 3.7, 4.9, 3.8,      # d=1
    5.3, 4.3, 9.0, 7.1, 18.0, 14.9, 32.8, 26.9,  # d=2
    6.6, 5.5, 10.6, 9.0, 21.7, 18.1, 38.8, 32.5, # d=3
    7.7, 6.4, 12.3, 10.8, 26.6, 21.9, 47.6, 40.8,# d=4
    8.8, 7.8, 14.0, 11.6, 27.7, 24.4, 54.7, 44.7,# d=5
    10.0, 8.5, 15.3, 13.2, 31.0, 26.8, 62.4, 53.0,# d=6
    11.2, 9.8, 17.2, 15.1, 36.9, 31.1, 66.3, 55.9,# d=7
    11.9, 10.2, 18.7, 16.2, 38.6, 32.7, 70.8, 63.8 # d=8
  ), nrow = 8, byrow = TRUE)
  
  
  return(crit_mat[d, col_idx + offset])
}



# get_critical_value_csms(T.chan = 5, alpha = 0.05, gamm = 0.15, d = 4)


# get_critical_value_csms(T.chan = 2, alpha = 0.10, gamm = 0, d = 6)
# ==============================================================================

# ==============================================================================
# CSMS Test Statistic with HAC Normalization (based on training sample)
# ==============================================================================

csms.statistic.fpca <- function(input.vec, m, T.chan, alpha = 0.05, gamm = 0) {
  library(sandwich)  # For HAC estimator
  library(MASS)      # For ginv()
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training sample mean
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Estimate HAC long-run variance (Newey-West via dummy regression)
  training_sample <- input_demeaned[1:m, ]
  ts_sample <- ts(training_sample)
  lm_fit <- lm(ts_sample ~ 1)  # Intercept-only model for each PC
  hac_cov <- m * NeweyWest(lm_fit, prewhite = FALSE)
  
  if (any(is.na(hac_cov)) || any(is.nan(hac_cov))) {
    stop("HAC covariance estimation failed.")
  }
  
  # Step 3: Cumulative sum of monitoring period
  monitor_sample <- input_demeaned[(m+1):(m + m*T.chan), ]
  cumsum_monitor <- apply(monitor_sample, 2, cumsum)
  
  # Step 4: Build CSMS statistic
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  # Final CSMS statistic
  csms_values <- diag(cumsum_monitor %*% ginv(hac_cov) %*% t(cumsum_monitor)) * scaling
  csms_stat <- max(csms_values)
  
  # Step 5: Compare with critical value from CSMS table
  crit_val <- get_critical_value_csms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- csms_stat > crit_val
  
  return(list(statistic = csms_stat, critical_value = crit_val, reject = reject))
}

# ==============================================================================
csms.statistic.fpca.alt <- function(input.vec, m, T.chan, alpha = 0.05, gamm = 0) {
  library(sandwich)  # HAC estimator
  library(MASS)      # Generalized inverse
  
  input.vec <- as.matrix(input.vec)
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # === Univariate case ===
  if (K == 1) {
    input_demeaned <- input.vec - mean(input.vec[1:m])
    
    # Step 2: Estimate HAC variance using Newey-West (variance only)
    ts_sample <- ts(input_demeaned[1:m, 1])
    lm_fit <- lm(ts_sample ~ 1)
    hac_var <- m * as.numeric(NeweyWest(lm_fit, prewhite = FALSE))
    
    if (is.na(hac_var) || is.nan(hac_var) || hac_var < 1e-8) {
      stop("HAC variance estimation failed or is too small.")
    }
    
    # Step 3: Cumulative sum of monitoring period
    cumsum_monitor <- cumsum(input_demeaned[(m+1):(m + m*T.chan), 1])
    
    # Step 4: Scaling
    time_idx <- 1:(m * T.chan)
    scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
      ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
    
    csms_values <- (cumsum_monitor^2 / hac_var) * scaling
    csms_stat <- max(csms_values)
    
    # Step 5: Critical value
    crit_val <- get_critical_value_csms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = 1)
    reject <- csms_stat > crit_val
    first_rejection <- if (reject) which(csms_values > crit_val)[1] else length(csms_values)
    
    return(list(
      statistic = csms_stat,
      critical_value = crit_val,
      reject = reject,
      first_rejection = first_rejection
    ))
  }
  
  # === Multivariate case ===
  colmeans_m <- colMeans(input.vec[1:m, , drop = FALSE])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  training_sample <- input_demeaned[1:m, , drop = FALSE]
  ts_sample <- ts(training_sample)
  lm_fit <- lm(ts_sample ~ 1)
  hac_cov <- m * NeweyWest(lm_fit, prewhite = FALSE)
  
  if (any(is.na(hac_cov)) || any(is.nan(hac_cov))) {
    stop("HAC covariance estimation failed.")
  }
  
  monitor_sample <- input_demeaned[(m+1):(m + m*T.chan), , drop = FALSE]
  cumsum_monitor <- apply(monitor_sample, 2, cumsum)
  
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  csms_values <- diag(cumsum_monitor %*% MASS::ginv(hac_cov) %*% t(cumsum_monitor)) * scaling
  csms_stat <- max(csms_values)
  
  crit_val <- get_critical_value_csms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- csms_stat > crit_val
  first_rejection <- if (reject) which(csms_values > crit_val)[1] else length(csms_values)
  
  return(list(
    statistic = csms_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}

csms.statistic.fpca.alt.old.version <- function(input.vec, m, T.chan, alpha = 0.05, gamm = 0) {
  library(sandwich)  # For HAC estimator
  library(MASS)      # For ginv()
  input.vec <- as.matrix(input.vec)  # <-- Force matrix structure
  
  sample.size <- nrow(input.vec)
  K <- ncol(input.vec)
  
  # Step 1: Demean using training sample mean
  colmeans_m <- colMeans(input.vec[1:m, ])
  input_demeaned <- sweep(input.vec, 2, colmeans_m)
  
  # Step 2: Estimate HAC long-run variance (Newey-West via dummy regression)
  training_sample <- input_demeaned[1:m, ]
  ts_sample <- ts(training_sample)
  lm_fit <- lm(ts_sample ~ 1)  # Intercept-only model for each PC
  hac_cov <- m * NeweyWest(lm_fit, prewhite = FALSE)
  
  if (any(is.na(hac_cov)) || any(is.nan(hac_cov))) {
    stop("HAC covariance estimation failed.")
  }
  
  # Step 3: Cumulative sum of monitoring period
  monitor_sample <- input_demeaned[(m+1):(m + m*T.chan), ]
  cumsum_monitor <- apply(monitor_sample, 2, cumsum)
  
  # Step 4: Build CSMS statistic sequence
  time_idx <- 1:(m * T.chan)
  scaling <- m^(-1) * (1 + time_idx / m)^(-2) *
    ((time_idx / m) / (1 + time_idx / m))^(-2 * gamm)
  
  csms_values <- diag(cumsum_monitor %*% ginv(hac_cov) %*% t(cumsum_monitor)) * scaling
  csms_stat <- max(csms_values)
  
  # Step 5: Compare with critical value
  crit_val <- get_critical_value_csms(T.chan = T.chan, alpha = alpha, gamm = gamm, d = K)
  reject <- csms_stat > crit_val
  first_rejection <- if (reject) which(csms_values > crit_val)[1] else length(csms_values)
  
  return(list(
    statistic = csms_stat,
    critical_value = crit_val,
    reject = reject,
    first_rejection = first_rejection
  ))
}



