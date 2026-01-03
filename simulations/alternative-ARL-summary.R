# Load libraries
library(dplyr)
library(stringr)
library(ggplot2)
# Set working directory  

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    !is.null(rstudioapi::getActiveDocumentContext()$path)) {
  
  # Running in RStudio IDE
  script.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script.dir)
  
} else {
  # Running via Rscript or outside RStudio
  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("--file=", args, value = TRUE)
  if (length(file.arg) > 0) {
    script.path <- normalizePath(sub("--file=", "", file.arg))
    setwd(dirname(script.path))
  }
}


# Directories
input_dir <- "output_alt_abrupt_local_change"
plot_dir <- paste0(input_dir, "_ARL_figs")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Metadata extraction
extract_gamma <- function(name) {
  gamma_str <- str_extract(name, "gamma[0-9p]+")
  gamma_num <- as.numeric(str_replace(str_replace(gamma_str, "gamma", ""), "p", "."))
  return(gamma_num)
}

parse_filename <- function(name) {
  list(
    dgp = str_extract(name, "(?<=monitoring_results_)[a-zA-Z0-9]+"),
    m = as.numeric(str_extract(name, "(?<=_m)\\d+")),
    T.chan = as.numeric(str_extract(name, "(?<=_T)\\d+")),
    gamma = extract_gamma(name),
    delta = as.numeric(str_replace(str_extract(name, "delta[0-9p]+"), "delta", "") %>% str_replace("p", ".")),
    s.star = as.numeric(str_extract(name, "(?<=sstar)\\d+"))
  )
}

# Load and bind all CSVs
all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
all_results <- lapply(all_files, function(file) {
  df <- read.csv(file)
  df$file <- basename(file)
  return(df)
}) %>% bind_rows()

meta <- lapply(all_results$file, parse_filename)
meta_df <- bind_rows(meta)
all_results <- bind_cols(all_results, meta_df)

# Compute ARL
arl_data <- all_results %>%
  group_by(dgp, s.star, m, T.chan, gamma, delta) %>%
  summarise(
    ARL_R = mean(first_rsms, na.rm = TRUE),
    ARL_S = mean(first_ssms, na.rm = TRUE),
    ARL_C = mean(first_csms, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("ARL"),
    names_to = "Method",
    values_to = "ARL"
  ) %>%
  mutate(Method = recode(Method, ARL_R = "RSMS", ARL_S = "SSMS", ARL_C = "CSMS"))

# Filter for m = 500 and T = 1, 2
arl_data <- arl_data %>%
  filter(m == 500, T.chan %in% c(1, 2))

# Loop: each gamma as a separate figure
grouped <- split(arl_data, list(arl_data$dgp, arl_data$s.star, arl_data$m, arl_data$T.chan, arl_data$gamma), drop = TRUE)

for (key in names(grouped)) {
  df <- grouped[[key]]
  dgp_val <- unique(df$dgp)
  sstar_val <- unique(df$s.star)
  m_val <- unique(df$m)
  T_val <- unique(df$T.chan)
  gamma_val <- unique(df$gamma)
  
 
  # Compute y-axis limits:
  min_y <- min(arl_data$ARL[
    arl_data$dgp == dgp_val &
      arl_data$s.star == sstar_val &
      arl_data$m == m_val &
      arl_data$T.chan == T_val
  ], na.rm = TRUE)
  
  max_y <- m_val * T_val  # Upper limit: m * T
  
  plot <- ggplot(df, aes(x = delta, y = ARL, linetype = Method)) +
    geom_line(size = 1.1, color = "black") +
    scale_linetype_manual(
      values = c("RSMS" = "solid", "SSMS" = "dashed", "CSMS" = "dotdash")
    ) +
    labs(
      title = "",
      x = expression(Delta),
      y = "ARL",
      linetype = "Method"
    ) +
    coord_cartesian(ylim = c(min_y, max_y)) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 28),
      axis.title.y = element_text(size = 28),
      axis.text = element_text(size = 28), 
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20) # add more space at the bottom
    )
  
  
  
  folder_tag <- basename(input_dir)
  
  file_name <- paste0(
    "ARL_", folder_tag,
    "_dgp", dgp_val,
    "_k", sstar_val,
    "_m", m_val,
    "_T", T_val,
    "_gamma", formatC(gamma_val, format = "f", digits = 2),
    ".png"
  )
  
  ggsave(filename = file.path(plot_dir, file_name), plot = plot, width = 9, height = 6, dpi = 300)
}
