# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : analysis_viviani.R                                         ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 24-11-2025                                                 ║
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
library(posterior)
library(dplyr)
library(tidyr)
library(lavaan)
library(kableExtra)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Load empirical functions to rotate factor loadings
source("R scripts/R functions/empirical_functions.R")

# Function to compute experimental effects in linear mixed models with covariates
make_task_mean <- function(df_row, congruent = TRUE){
  ## Intercept
  mu <- df_row[["(Intercept)"]]
  
  ## Check if the first task it's not peripheal (the baseline)
  blk <- c("", "BlockPerifoveal", "BlockNavon",
           "BlockFigureGround", "BlockFlanker", "BlockSaliency")
  
  res <- numeric(6)
  for (j in 1:6){
    mu_j <- mu + if(j > 1) df_row[[ blk[j] ]] else 0
    if(!congruent){                                 
      mu_j <- mu_j + df_row[["CongC0"]] +
        if(j > 1) df_row[[ paste0(blk[j], ":CongC0") ]] else 0
    }
    res[j] <- mu_j
  }
  res
}

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

# Load fitted ex-gaussian HFM
exg_HFM  <- readRDS("Results/Stan/Viviani/Fitted models/viviani_exgaussian_GHFM.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Hierarchical models
# ─────────────────────────────────────────────────────────────────────────────

# Load fitted models
LMM_models <- readRDS("Results/Rdata/LMM/viviani_LMM_models.rds")

# LMM with covariates effects
fix_eff <- fixef(LMM_models$LMM_full)
rnd_eff <- ranef(LMM_models$LMM_full)$SSID |>               
  tibble::as_tibble(rownames = "SSID")

# Sum fixed plus random effects (only for coefs with random effects)
common_coef <- intersect(names(fix_eff), names(rnd_eff))
full_eff <- rnd_eff |>  
  mutate(across(all_of(common_coef), ~ . + fix_eff[cur_column()]))

# Add the fixed effects that don't appear in ranef (e.g., Trial_z, etc...)
only_fixed <- setdiff(names(fix_eff), names(full_eff))
for(nm in only_fixed){ full_eff[[nm]] <- fix_eff[[nm]] }

# Save task names to compute experimental effects
task_names <- c("Peripheral", "Perifoveal", "Navon",
                "FigureGround", "Flanker", "Saliency")

# Compute congruent expected value
C_scores <- full_eff %>%
  rowwise() %>%
  mutate(tmp = list(make_task_mean(cur_data(), TRUE))) %>%
  unnest_wider(tmp, names_sep = "_") %>%
  dplyr::select(starts_with("tmp_")) %>%
  setNames(task_names)

# Compute incongruent expected value
I_scores <- full_eff %>%
  rowwise() %>%
  mutate(tmp = list(make_task_mean(cur_data(), FALSE))) %>%
  unnest_wider(tmp, names_sep = "_") %>%
  dplyr::select(starts_with("tmp_")) %>%
  setNames(task_names)

# Compute experimental effect per subject and task
LMM_full_scores <- I_scores - C_scores

# Compute experimental effects without covariates (i.e., random slopes)
LMM_basic_scores <- coef(LMM_models$LMM_basic)$SSID[,-c(1:6)]
colnames(LMM_basic_scores) <- colnames(LMM_full_scores)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: ELPDs difference table
# ─────────────────────────────────────────────────────────────────────────────

# Load ELPDs difference table
viviani_ELPDs <- readRDS("Results/Rdata/ELPDs/viviani_HFM_ELPDs.rds")
ELPD_LMMs <- readRDS("Results/Rdata/ELPDs/viviani_LMM_ELPDs.rds")

# Data frame with ELPDs comparison
data.frame(
  model_1 = rep("exgaussian", 4),
  model_2 = c("shifted-lognormal", "LMM_full", "LMM_basic", "normal"),
  elpd_diff = c(
    sum(viviani_ELPDs$ELPD_exgaussian$elpd_i - viviani_ELPDs$ELPD_shlognormal$elpd_i),
    sum(viviani_ELPDs$ELPD_exgaussian$elpd_i - ELPD_LMMs$ELPD_LMM_full$elpd_i),
    sum(viviani_ELPDs$ELPD_exgaussian$elpd_i - ELPD_LMMs$ELPD_LMM_basic$elpd_i),
    sum(viviani_ELPDs$ELPD_exgaussian$elpd_i - viviani_ELPDs$ELPD_gaussian$elpd_i)
  ),
  elpd_diff_se = c(
    sd(viviani_ELPDs$ELPD_exgaussian$elpd_i - viviani_ELPDs$ELPD_shlognormal$elpd_i) * sqrt(length(viviani_ELPDs$ELPD_exgaussian$elpd_i)),
    sd(viviani_ELPDs$ELPD_exgaussian$elpd_i - ELPD_LMMs$ELPD_LMM_full$elpd_i) * sqrt(length(viviani_ELPDs$ELPD_exgaussian$elpd_i)),
    sd(viviani_ELPDs$ELPD_exgaussian$elpd_i - ELPD_LMMs$ELPD_LMM_basic$elpd_i) * sqrt(length(viviani_ELPDs$ELPD_exgaussian$elpd_i)),
    sd(viviani_ELPDs$ELPD_exgaussian$elpd_i - viviani_ELPDs$ELPD_gaussian$elpd_i) * sqrt(length(viviani_ELPDs$ELPD_exgaussian$elpd_i))
  ),
  #  Probability that the ex-gaussina model will outperform the other models when predicting new data
  P_better = round(c(1 - pnorm(0, 
                               mean = sum(viviani_ELPDs$ELPD_exgaussian$elpd_i - viviani_ELPDs$ELPD_shlognormal$elpd_i), 
                               sd = sd(viviani_ELPDs$ELPD_exgaussian$elpd_i - viviani_ELPDs$ELPD_shlognormal$elpd_i) * sqrt(length(viviani_ELPDs$ELPD_exgaussian$elpd_i))),
                     1 - pnorm(0, 
                               mean = sum(viviani_ELPDs$ELPD_exgaussian$elpd_i - ELPD_LMMs$ELPD_LMM_full$elpd_i), 
                               sd = sd(viviani_ELPDs$ELPD_exgaussian$elpd_i - ELPD_LMMs$ELPD_LMM_full$elpd_i) * sqrt(length(viviani_ELPDs$ELPD_exgaussian$elpd_i))),
                     1 - pnorm(0, 
                               mean = sum(viviani_ELPDs$ELPD_exgaussian$elpd_i - ELPD_LMMs$ELPD_LMM_basic$elpd_i), 
                               sd = sd(viviani_ELPDs$ELPD_exgaussian$elpd_i - ELPD_LMMs$ELPD_LMM_basic$elpd_i) * sqrt(length(viviani_ELPDs$ELPD_exgaussian$elpd_i))), 
                     1 - pnorm(0, 
                               mean = sum(viviani_ELPDs$ELPD_exgaussian$elpd_i - viviani_ELPDs$ELPD_gaussian$elpd_i), 
                               sd = sd(viviani_ELPDs$ELPD_exgaussian$elpd_i - viviani_ELPDs$ELPD_gaussian$elpd_i) * sqrt(length(viviani_ELPDs$ELPD_exgaussian$elpd_i)))), 3))

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: model-implied factor loadings
# ─────────────────────────────────────────────────────────────────────────────

# Frequentist factor loadings
LMM_full_EFA <- efa(data = apply(LMM_full_scores, 2, scale), nfactors = 2, 
                    estimator = "ml", rotation = "varimax", output = "lavaan")
LMM_basic_EFA <- efa(data = apply(LMM_basic_scores, 2, scale), nfactors = 2, 
                     estimator = "ml", rotation = "varimax", output = "lavaan")

# Rotate all factor loadings
exgaussian_loads <- MatchAlign(fit = exg_HFM, lambda_name = "Lambda_std", 
                               only_summary = FALSE)
aligned_exgauss_loads <- summarise_draws(exgaussian_loads$rotated_draws, mean, 
                                         quantile =~ quantile(.x, c(.025, .975)))

# Factor loadings results table
lambda_table <- cbind(
  aligned_exgauss_loads[,1], 
  standardizedsolution(LMM_basic_EFA)[1:12, c("est.std", "ci.lower", "ci.upper")],
  standardizedsolution(LMM_full_EFA)[1:12, c("est.std", "ci.lower", "ci.upper")], 
  aligned_exgauss_loads[,-1])

# Table labels
factors <- rep(c("Factor 1", "Factor 2"), each = 6)
tasks <- rep(c("Peripheral", "Perifoveal", "Navon", 
               "Figure-Ground", "Flanker", "Saliency"), 2)

# Detect significant correlations for each model
LMM_basic_vals  <- mapply(format_cell, lambda_table[,2], lambda_table[,3], lambda_table[,4], latex = FALSE)
LMM_full_vals <- mapply(format_cell, lambda_table[,5], lambda_table[,6], lambda_table[,7], latex = FALSE)
eg_vals <- mapply(format_cell, lambda_table[,8], lambda_table[,9], lambda_table[,10], latex = FALSE)

# Final data frame for kable
df_final <- data.frame(
  Factor = factors,
  Task = tasks,
  LMM_B_Mean = unlist(LMM_basic_vals["mean",]),
  LMM_B_CI = unlist(LMM_basic_vals["ci",]),
  LMM_F_Mean = unlist(LMM_full_vals["mean",]),
  LMM_F_CI = unlist(LMM_full_vals["ci",]),
  EG_Mean = unlist(eg_vals["mean",]),
  EG_CI = unlist(eg_vals["ci",])
)

# Kable table
kbl(df_final, 
    col.names = c("Factor", "Task", 
                  "$\\hat{\\lambda}$", "95% CI", 
                  "$\\hat{\\lambda}$", "95% CI", 
                  "$\\hat{\\lambda}$", "95% CI"),
    escape = FALSE, 
    booktabs = TRUE,
    align = c("l", "l", "c", "c", "c", "c", "c", "c"), 
    caption = "Model-implied factor loadings") |> 
  kable_styling(full_width = F) |> 
  add_header_above(c(" " = 2, "HM (no predictors)" = 2, "HM (with predictors)" = 2, "Ex-Gaussian" = 2)) |>
  collapse_rows(columns = 1, valign = "middle") |> 
  row_spec(0, bold = TRUE) |> 
  row_spec(6, hline_after = TRUE)

# Detect significant correlations for each model
LMM_basic_vals  <- mapply(format_cell, lambda_table[,2], lambda_table[,3], lambda_table[,4], latex = TRUE)
LMM_full_vals <- mapply(format_cell, lambda_table[,5], lambda_table[,6], lambda_table[,7], latex = TRUE)
eg_vals <- mapply(format_cell, lambda_table[,8], lambda_table[,9], lambda_table[,10], latex = TRUE)

# Final LaTeX data frame
viviani_loads <- data.frame(
  Factor = factors,
  Task = tasks,
  LMM_B_Mean = unlist(LMM_basic_vals["mean",]),
  LMM_B_CI = unlist(LMM_basic_vals["ci",]),
  LMM_F_Mean = unlist(LMM_full_vals["mean",]),
  LMM_F_CI = unlist(LMM_full_vals["ci",]),
  EG_Mean = unlist(eg_vals["mean",]),
  EG_CI = unlist(eg_vals["ci",])
)

save(viviani_loads, file = "Results/Rdata/Tables/viviani_loads.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Reliability of experimental effects
# ─────────────────────────────────────────────────────────────────────────────
# Hierarchical models split-half reliability values (estimates are medians)
# Retrieved from table S55 and S56 (https://osf.io/5sm9j/files/djszu)
viviani_reliab <- data.frame(
  Task = tasks[1:6],
  LMM_basic_est = c("0.750","0.769","0.708","0.904","0.933","0.684"),
  LMM_basic_CI  = c("[0.645, 0.828]", "[0.642, 0.854]", 
                    "[0.440, 0.830]", "[0.833, 0.947]",
                    "[0.872, 0.966]", "[0.400, 0.834]"),
  LMM_full_est = c("0.812","0.814","0.762","0.927","0.950","0.767"),
  LMM_full_CI  = c("[0.715, 0.876]", "[0.699, 0.884]", 
                   "[0.558, 0.864]", "[0.864, 0.966]",
                   "[0.907, 0.980]", "[0.570, 0.867]"),
  GenHFM_est = sprintf("%.3f", c(exg_HFM$summary("reliability", median)$median)),
  GenHFM_CI = paste0(
    "[", sprintf("%.3f", unlist(exg_HFM$summary("reliability", quantile =~ quantile(.x, .025))[,2])), 
    ", ",
    sprintf("%.3f", unlist(exg_HFM$summary("reliability", quantile =~ quantile(.x, .975))[,2])), 
    "]")
  )

# Final LaTeX table
kbl(viviani_reliab, 
    col.names = c("Task", 
                  "$\\hat{\\rho}$", "95% CI", 
                  "$\\hat{\\rho}$", "95% CI", 
                  "$\\hat{\\rho}$", "95% CI"),
    escape = FALSE, 
    booktabs = TRUE,
    align = c("l", "c", "c", "c", "c", "c", "c"),
    caption = "Reliability estimates of experimental effects using hierarchical models (HM) with and without predictors, and the Ex-Gaussian HFM</i>") |> 
  kable_styling(full_width = F) |> 
  add_header_above(c(" " = 1, 
                     "HM (no predictors)" = 2, 
                     "HM (with predictors)" = 2, 
                     "Ex-Gaussian HFM" = 2))

# Final LaTeX data frame
save(viviani_reliab, file = "Results/Rdata/Tables/viviani_reliab.rdata")

# ─────────────────────────────────────────────────────────────────────────────
