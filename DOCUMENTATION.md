
This repository contains R code for the simulation studies, critical value computation,
and empirical illustrations in the manuscript:

**“Online Monitoring of Functional Time Series with Adjusted-Range Self-Normalization: Applications to Financial Data”**

The code focuses on online change-point monitoring for (possibly dependent) functional
time series, using FPCA-based scores and HAC, Shao's (2010) self-normalization and the adjusted-range based self-normalization (Hong et al., 2024).

---

## Repository structure

- `critical-values/`
  - `RSMS-critical-values.R` : Monte Carlo critical values for RSMS
  - `SSMS-critical-values.R` : Monte Carlo critical values for SSMS
  - `CSMS-critical-values.R` : Monte Carlo critical values for CSMS
- `simulations/`
  - `Method.R` : core implementations (RSMS/SSMS/CSMS and utilities)
  - `genData.R` : data generating processes for functional time series
  - `null.R` : size experiments under the null
  - `alternative-*.R` : power / detection experiments under various alternatives
  - `alternative-ARL-summary.R` : summarise outputs (ARL / rejection rates)
- `empirical/`
  - `SPX.csv` : SPX/SP500 (not included due to copyright concern)
  - `sp500-plot` : construct functional data and plot mean curves
  - `sp500.R` : monitoring test on rates of return for SPX/SP500 data (requires `fdata`)
 
---

## Requirements

- R (>= 4.0 recommended)
- Common packages used in `simulations/Method.R` and `empirical/Method.R` include:
  `fdapace`, `fda`, `MASS`, `parallel`, `zoo`, `sandwich`, `fastmatrix`
- Simulation scripts also use: `pbmcapply` (parallel with progress bar)
- Critical value scripts use: `foreach`, `doSNOW`, `parallel`, `MASS`

Note: `Method.R` attempts to install missing packages automatically. If you prefer
manual installation, install the required packages first via `install.packages()`.

---

