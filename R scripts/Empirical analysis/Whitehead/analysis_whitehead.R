# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : analysis_whitehead.R                                       ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 20-11-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Main results and tables reported in the main text
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(cmdstanr)
library(dplyr)
library(kableExtra)
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: user-defined functions
# ─────────────────────────────────────────────────────────────────────────────

# Function for detect significant values and put it in bold
format_cell <- function(mean_val, lower, upper, latex = FALSE) {
  # 1. Verificar significancia
  is_sig <- (lower > 0) | (upper < 0)
  
  # 2. Formatear números a 2 decimales
  m_str <- sprintf("%.2f", mean_val)
  ci_str <- paste0("[", sprintf("%.2f", lower), ", ", sprintf("%.2f", upper), "]")
  
  # 3. Aplicar negrita CONDICIONAL (Detecta si es PDF/LaTeX o HTML)
  if (is_sig) {
    if (latex) {
      # Si estamos compilando a PDF (LaTeX)
      m_str <- paste0("\\textbf{", m_str, "}")
    } else {
      # Para todo lo demás (HTML, Viewer de RStudio, Word)
      m_str <- paste0("<b>", m_str, "</b>")
    }
  }
  
  return(list(mean = m_str, ci = ci_str))
}

# ELPDs differences
Delta_ELPD <- function(bfit, wfit) {
  D_ELPD <- sum(bfit$elpd_i - wfit$elpd_i)
  SE_D_ELPD <- sqrt(var(bfit$elpd_i - wfit$elpd_i) * bfit$N)
  c(D_ELPD, SE_D_ELPD)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: ELPD differences
# ─────────────────────────────────────────────────────────────────────────────

# Load experiments ELPDs
ELPDs_exps <- readRDS("Results/Rdata/ELPDs/whitehead_experiments_ELPDS.rds")
ELPDs_joint <- readRDS("Results/Rdata/ELPDs/whitehead_joint_ELPDS.rds")

# List wiht both ELPDs
ELPDs_list <- c(ELPDs_exps, ELPDs_joint)

# Suffix IDs for ELPDs per experiment and joint experiments
exps <- c("E1", "E2", "E3", "joint")

# Auxiliar function to compute ELPDs
get_comparison <- function(exp_suffix, target_model) {
  ref_val <- ELPDs_list[[paste0("ELPD_exgaussian_", exp_suffix)]]
  comp_val <- ELPDs_list[[paste0("ELPD_", target_model, "_", exp_suffix)]]
  res <- Delta_ELPD(ref_val, comp_val)
  return(round(res, 2))
}

# First column of the ELPDs data frame
ELPDs_table <- data.frame(Model = c("Ex-Gaussian", "Sh-lognormal", "Gaussian"))

# Iteramos sobre cada experimento para añadir sus dos columnas
for (e in exps) {
  # ex-gaussian vs shifted-lognormal
  sh_res <- get_comparison(e, "shlognormal")
  # ex-gaussian vs gaussian
  ga_res <- get_comparison(e, "gaussian")
  # Add this columns to our final table
  ELPDs_table[[paste0("ELPD_", e)]] <- c("0.00", as.character(sh_res[1]), as.character(ga_res[1]))
  ELPDs_table[[paste0("SE_", e)]]   <- c("---", as.character(sh_res[2]), as.character(ga_res[2]))
}

# Save ELPDs table
save(ELPDs_table, file = "Results/Rdata/Tables/whitehead_ELPDs_table.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: load fitted models
# ─────────────────────────────────────────────────────────────────────────────

# Load gaussian HFM fitted models
gaussian_E1  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E1.rds")
gaussian_E2  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E2.rds")
gaussian_E3  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E3.rds")
gaussian_all <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_joint.rds")

# Load ex-gaussian HFM fitted models
exgaussian_E1  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E1.rds")
exgaussian_E2  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E2.rds")
exgaussian_E3  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E3.rds")
exgaussian_all <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_joint.rds")

# Load ex-gaussian HFM fitted models
shlognormal_E1  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E1.rds")
shlognormal_E2  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E2.rds")
shlognormal_E3  <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E3.rds")
shlognormal_all <- readRDS("Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_joint.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: model-implied true correlations
# ─────────────────────────────────────────────────────────────────────────────

# Correlations
rho_ids <- c("Rho_beta[2,1]", "Rho_beta[3,1]", "Rho_beta[3,2]")

# Correlation results table
rho_table <- cbind(
  rbind(
    gaussian_E1$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975))),
    gaussian_E2$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975))),
    gaussian_E3$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975))),
    gaussian_all$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))
  ), 
  rbind(
    exgaussian_E1$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    exgaussian_E2$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    exgaussian_E3$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    exgaussian_all$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))[,-1]
  ),
  rbind(
    shlognormal_E1$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    shlognormal_E2$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    shlognormal_E3$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    shlognormal_all$summary(rho_ids, mean, quantile =~ quantile(.x, c(.025, .975)))[,-1]
  )
)

# Labels
experiments <- rep(c("Experiment 1", "Experiment 2", "Experiment 3", "All experiments"), each = 3)
pairs <- rep(c("Simon–Flanker", "Simon–Stroop", "Flanker–Stroop"), 4)

# Detect significant correlations for each model
g_vals <- mapply(format_cell, rho_table[,2], rho_table[,3], rho_table[,4], latex = FALSE)
eg_vals <- mapply(format_cell, rho_table[,5], rho_table[,6], rho_table[,7], latex = FALSE)
sl_vals <- mapply(format_cell, rho_table[,8], rho_table[,9], rho_table[,10], latex = FALSE)

# Final data frame for kable
df_final <- data.frame(
  Experiment = experiments,
  Pair = pairs,
  G_Mean = unlist(g_vals["mean",]),
  G_CI = unlist(g_vals["ci",]),
  EG_Mean = unlist(eg_vals["mean",]),
  EG_CI = unlist(eg_vals["ci",]),
  SL_Mean = unlist(sl_vals["mean",]),
  SL_CI = unlist(sl_vals["ci",])
)

# Kable table
kbl(df_final, 
    col.names = c("Experiment", "Pair", 
                  "$\\hat{\\rho}$", "95% CI", 
                  "$\\hat{\\rho}$", "95% CI", 
                  "$\\hat{\\rho}$", "95% CI"),
    escape = FALSE, 
    booktabs = TRUE,
    align = c("l", "l", "c", "c", "c", "c", "c", "c"), 
    caption = "Model-implied correlation estimates for the three HFMs across Whitehead et al. experiments") |> 
  kable_styling(full_width = F) |> 
  add_header_above(c(" " = 2, "Gaussian" = 2, "Ex-Gaussian" = 2, "Sh-lognormal" = 2)) |>
  collapse_rows(columns = 1, valign = "middle") |> 
  row_spec(0, bold = TRUE) |> 
  row_spec(c(3, 6, 9), hline_after = TRUE)

# Save this table for LaTeX table in the main text
g_vals <- mapply(format_cell, rho_table[,2], rho_table[,3], rho_table[,4], latex = TRUE)
eg_vals <- mapply(format_cell, rho_table[,5], rho_table[,6], rho_table[,7], latex = TRUE)
sl_vals <- mapply(format_cell, rho_table[,8], rho_table[,9], rho_table[,10], latex = TRUE)

# Final LaTeX data frame
whitehead_cors <- data.frame(
  Experiment = experiments,
  Pair = pairs,
  G_Mean = unlist(g_vals["mean",]),
  G_CI = unlist(g_vals["ci",]),
  EG_Mean = unlist(eg_vals["mean",]),
  EG_CI = unlist(eg_vals["ci",]),
  SL_Mean = unlist(sl_vals["mean",]),
  SL_CI = unlist(sl_vals["ci",])
)

save(whitehead_cors, file = "Results/Rdata/Tables/whitehead_cors.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: model-implied factor loadings
# ─────────────────────────────────────────────────────────────────────────────

# Correlation results table
lambda_table <- cbind(
  rbind(
    gaussian_E1$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975))),
    gaussian_E2$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975))),
    gaussian_E3$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975))),
    gaussian_all$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))
  ), 
  rbind(
    exgaussian_E1$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    exgaussian_E2$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    exgaussian_E3$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    exgaussian_all$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))[,-1]
  ),
  rbind(
    shlognormal_E1$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    shlognormal_E2$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    shlognormal_E3$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))[,-1],
    shlognormal_all$summary("Lambda_std", mean, quantile =~ quantile(.x, c(.025, .975)))[,-1]
  )
)

# Labels
experiments <- rep(c("Experiment 1", "Experiment 2", "Experiment 3", "All experiments"), each = 3)
pairs <- rep(c("Simon", "Flanker", "Stroop"), 4)

# Detect significant correlations for each model
g_vals <- mapply(format_cell, lambda_table[,2], lambda_table[,3], lambda_table[,4], latex = FALSE)
eg_vals <- mapply(format_cell, lambda_table[,5], lambda_table[,6], lambda_table[,7], latex = FALSE)
sl_vals <- mapply(format_cell, lambda_table[,8], lambda_table[,9], lambda_table[,10], latex = FALSE)

# Final data frame for kable
df_final <- data.frame(
  Experiment = experiments,
  Pair = pairs,
  G_Mean = unlist(g_vals["mean",]),
  G_CI = unlist(g_vals["ci",]),
  EG_Mean = unlist(eg_vals["mean",]),
  EG_CI = unlist(eg_vals["ci",]),
  SL_Mean = unlist(sl_vals["mean",]),
  SL_CI = unlist(sl_vals["ci",])
)

# Kable table
kbl(df_final, 
    col.names = c("Experiment", "Pair", 
                  "$\\hat{\\lambda}$", "95% CI", 
                  "$\\hat{\\lambda}$", "95% CI", 
                  "$\\hat{\\lambda}$", "95% CI"),
    escape = FALSE, 
    booktabs = TRUE,
    align = c("l", "l", "c", "c", "c", "c", "c", "c"), 
    caption = "Model-implied correlation estimates for the three HFMs across Whitehead et al. experiments") |> 
  kable_styling(full_width = F) |> 
  add_header_above(c(" " = 2, "Gaussian" = 2, "Ex-Gaussian" = 2, "Sh-lognormal" = 2)) |>
  collapse_rows(columns = 1, valign = "middle") |> 
  row_spec(0, bold = TRUE) |> 
  row_spec(c(3, 6, 9), hline_after = TRUE)

# Save this table for LaTeX table in the main text
g_vals <- mapply(format_cell, lambda_table[,2], lambda_table[,3], lambda_table[,4], latex = TRUE)
eg_vals <- mapply(format_cell, lambda_table[,5], lambda_table[,6], lambda_table[,7], latex = TRUE)
sl_vals <- mapply(format_cell, lambda_table[,8], lambda_table[,9], lambda_table[,10], latex = TRUE)

# Final LaTeX data frame
whitehead_loads <- data.frame(
  Experiment = experiments,
  Pair = pairs,
  G_Mean = unlist(g_vals["mean",]),
  G_CI = unlist(g_vals["ci",]),
  EG_Mean = unlist(eg_vals["mean",]),
  EG_CI = unlist(eg_vals["ci",]),
  SL_Mean = unlist(sl_vals["mean",]),
  SL_CI = unlist(sl_vals["ci",])
)

save(whitehead_loads, file = "Results/Rdata/Tables/whitehead_loads.rdata")

# ─────────────────────────────────────────────────────────────────────────────
