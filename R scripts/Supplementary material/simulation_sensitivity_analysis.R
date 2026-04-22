# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : simulation_sensitivity_analysis.R                          ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 31-03-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Analyze results from the Sensitivity Simulation Study
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(stringr)
library(afex)
library(emmeans)
library(effectsize)
library(kableExtra)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-defined functions
# ─────────────────────────────────────────────────────────────────────────────

# Fast function to generate data frame per measure (sensitivity version)
data_split_sens <- function(df_long, par, measure, statistic) {
  # Build filter: keep only the 3 HFM models
  idx <- df_long$measure == measure &
         df_long$statistic == statistic &
         df_long$model_pars %in% c("gaussian", "shlognormal", "exgaussian")
  
  # Handle NA parameter (model-level diagnostics like divergence, treedepth)
  if (is.na(par)) {
    idx <- idx & is.na(df_long$parameter)
  } else {
    idx <- idx & !is.na(df_long$parameter) & df_long$parameter == par
  }
  
  df <- as.data.frame(df_long[which(idx),])
  colnames(df) <- c("R", "L", "measure", "statistic", "parameter", "M", "value")
  for(j in 1:(ncol(df)-1)) {df[,j] <- as.factor(df[,j])}
  
  # 3 models per design cell → each = 3
  df$id <- rep(1:(nrow(df)/3), each = 3)
  return(df)
}

# Build emmeans table with MEASURES as header groups and MODELS as columns within
# Model order: G (Gaussian), SL (Shifted-Lognormal), EG (Ex-Gaussian)
build_measure_table <- function(anova_fits, measure_names, highlight_crits, format = "html") {
  
  models <- c("gaussian", "shlognormal", "exgaussian")  # G, SL, EG
  
  # Extract emmeans for each ANOVA fit
  emmeans_all <- lapply(anova_fits, function(fit) {
    em_R <- as.data.frame(emmeans(fit, specs = ~ R + M)) |>
      mutate(Factor_Name = "R", Level = as.character(R)) |>
      select(Factor_Name, Level, M, emmean)
    em_L <- as.data.frame(emmeans(fit, specs = ~ L + M)) |>
      mutate(Factor_Name = "L", Level = as.character(L)) |>
      select(Factor_Name, Level, M, emmean)
    bind_rows(em_R, em_L)
  })
  
  # Create base with Factor_Name and Level
  base <- emmeans_all[[1]] |>
    filter(M == models[1]) |>
    select(Factor_Name, Level) |>
    mutate(
      Level = factor(Level, levels = c("0.3", "0.5", "0.7")),
      Factor_Name = factor(Factor_Name, levels = c("R", "L"))
    ) |>
    arrange(Factor_Name, Level) |>
    mutate(Factor_Name = recode(Factor_Name,
                                "R" = "Reliability",
                                "L" = "True $\\Lambda$"))
  
  # For each MEASURE, add columns for each MODEL (G, SL, EG order)
  for (i in seq_along(measure_names)) {
    for (m in models) {
      vals <- emmeans_all[[i]] |>
        filter(M == m) |>
        mutate(
          Level = factor(Level, levels = c("0.3", "0.5", "0.7")),
          Factor_Name = factor(Factor_Name, levels = c("R", "L"))
        ) |>
        arrange(Factor_Name, Level) |>
        pull(emmean)
      base[[paste0(measure_names[i], "_", m)]] <- vals
    }
  }
  
  # Apply highlighting: for each measure, compare across models within each row
  for (i in seq_along(measure_names)) {
    col_names <- paste0(measure_names[i], "_", models)
    col_idx <- match(col_names, names(base))
    
    for (row in 1:nrow(base)) {
      raw_vals <- as.numeric(base[row, col_idx])
      
      best <- if (highlight_crits[i] == "min") {
        which.min(raw_vals)
      } else if (highlight_crits[i] == "max") {
        which.max(raw_vals)
      } else if (highlight_crits[i] == "0.95") {
        which.min(abs(raw_vals - 0.95))
      }
      
      for (j in seq_along(col_idx)) {
        val_str <- sprintf("%.3f", raw_vals[j])
        if (length(best) > 0 && j == best) {
          val_str <- if (format == "latex") paste0("\\textbf{", val_str, "}")
                     else paste0("<b>", val_str, "</b>")
        }
        base[row, col_idx[j]] <- val_str
      }
    }
  }
  
  base <- base |> mutate(across(everything(), as.character))
  return(base)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Load and prepare final data frame
# ─────────────────────────────────────────────────────────────────────────────

# Load aggregated simulation results
simres <- readRDS("Results/Rdata/Supplementary Material/sensitivity_simulation_results.rds")

df_long <- simres[, 1:(which(colnames(simres) == "REPLICATIONS")-1)] |>
  pivot_longer(
    cols = -c(I, M, reliab, std_lambda_j),
    names_to = "name",
    values_to = "value"
  ) |>
  select(-I, -M) |>
  mutate(
    measure = case_when(
      str_starts(name, "Tucker")       ~ "Tucker",
      str_starts(name, "gauss_first")  ~ "loo",
      str_starts(name, "gauss_best")   ~ "loo",
      str_starts(name, "loo_")         ~ "loo",
      TRUE ~ str_extract(
        name,
        paste(
          c("abias","arbias","stdbias","rmse","srmr","postsd","empirical",
            "rhat","essb","esst","ECR","power","divergence",
            "treedepth","ebfmi","time"),
          collapse = "|"
        )
      )
    ),
    fmatch = str_match(name, "(?:^|_)f([12])_(avg|disp)(?:_|$)"),
    f_idx = if_else(!is.na(fmatch[,1]), paste0("f", fmatch[,2]), NA_character_),
    statistic_from_f = if_else(!is.na(fmatch[,1]), fmatch[,3], NA_character_),
    statistic = coalesce(statistic_from_f,
                         if_else(str_detect(name, "(?:^|_)avg(?:_|$)"), "avg", NA_character_),
                         if_else(str_detect(name, "(?:^|_)disp(?:_|$)"), "disp", NA_character_)
    ),
    parameter = case_when(
      str_starts(name, "Tucker") ~ "global",
      TRUE ~ str_extract(
        name,
        paste(
          c("Rho","lambda","reliab",
            "mu_alpha","mu_theta","mu_sigma",
            "sd_alpha","sd_theta","sd_sigma",
            "elpd","both_weak","both_strong",
            "gauss_first","gauss_best",
            "weak","strong"),
          collapse = "|"
        )
      )
    ),
    model_pars = case_when(
      parameter %in% c("mu_alpha","mu_theta","mu_sigma",
                       "sd_alpha","sd_theta","sd_sigma") ~ "gaussian",
      str_detect(name, "_exgaussian")  ~ "exgaussian",
      str_detect(name, "_shlognormal") ~ "shlognormal",
      str_detect(name, "_gaussian")    ~ "gaussian",
      str_detect(name, "_mm")          ~ "mm",
      str_detect(name, "_sp")          ~ "spearman",
      str_detect(name, "_exg")         ~ "comparison",
      str_detect(name, "_shl")         ~ "comparison",
      str_starts(name, "gauss_best")   ~ "comparison",
      TRUE ~ NA_character_
    )
  )

# Tucker congruence: average f1 within each statistic (unidimensional → only f1)
tucker_df <- df_long |>
  filter(measure == "Tucker") |>
  group_by(reliab, std_lambda_j, measure, statistic, parameter, model_pars) |>
  summarise(
    value = mean(value[f_idx %in% c(NA, "f1")], na.rm = TRUE),
    .groups = "drop"
  )

# Others: standard averaging
others_df <- df_long |>
  filter(measure != "Tucker") |>
  group_by(reliab, std_lambda_j, measure, statistic, parameter, model_pars) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

df_long <- bind_rows(tucker_df, others_df)

# Summarise divided by 3000 (main study: 3 chains × 1000)
df_long <- df_long |>
  mutate(value = if_else(
    measure %in% c("divergence", "treedepth") & statistic %in% c("avg", "disp"),
    value * 3/2,
    value
  ))

# Add ESS% measure: proportion of independent samples
# Total draws = chains (2) × iter_sampling (1000) = 2000
ess_pct <- df_long |>
  filter(measure == "essb" & statistic == "avg") |>
  mutate(
    value = value / 2000 * 100,
    measure = "ess_pct"
  )
df_long <- bind_rows(df_long, ess_pct)

# Save long data frame
saveRDS(object = df_long, file = "Results/Rdata/Supplementary Material/sensitivity_simulation_longdata.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Convergence, efficiency and LOO model comparison
# ─────────────────────────────────────────────────────────────────────────────

df_long |> 
  filter(measure == "rhat") |> 
  summarise(rhat_avg = mean(value, na.rm = TRUE), .by = statistic) |> 
  as.data.frame()

df_long |> 
  filter(measure == "treedepth") |> 
  summarise(treedepth_avg = mean(value, na.rm = TRUE), 
            .by = c(statistic, model_pars)) |> 
  as.data.frame() |> na.omit()

df_long |> 
  filter(measure == "divergence") |> 
  summarise(divergence_avg = mean(value, na.rm = TRUE), 
            .by = c(statistic, model_pars)) |> 
  as.data.frame() |> na.omit()

df_long |> 
  filter(measure == "ebfmi") |> 
  summarise(ebfmi_avg = mean(value, na.rm = TRUE), 
            .by = c(statistic, model_pars)) |> 
  as.data.frame() |> na.omit()

df_long |>
  filter(measure == "loo" & !is.na(value)) |>
  as.data.frame()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: ANOVA meta-models
# ─────────────────────────────────────────────────────────────────────────────

# Correlations bias measures
ASRMR_Rho <- aov_car(value ~ R + L + Error(id/M), 
                     data = data_split_sens(df_long, measure = "srmr",
                                            statistic = "avg", par = "Rho"))
AAB_Rho <- aov_car(value ~ R + L + Error(id/M), 
                   data = data_split_sens(df_long, measure = "abias",
                                          statistic = "avg", par = "Rho"))
APSD_Rho <- aov_car(value ~ R + L + Error(id/M), 
                    data = data_split_sens(df_long, measure = "postsd",
                                           statistic = "avg", par = "Rho"))
AECR_Rho <- aov_car(value ~ R + L + Error(id/M), 
                    data = data_split_sens(df_long, measure = "ECR",
                                           statistic = "avg", par = "Rho"))

# Factor loadings bias measures
ARMSE_lam <- aov_car(value ~ R + L + Error(id/M), 
                     data = data_split_sens(df_long, measure = "rmse",
                                            statistic = "avg", par = "lambda"))
AAB_lam <- aov_car(value ~ R + L + Error(id/M), 
                   data = data_split_sens(df_long, measure = "abias",
                                          statistic = "avg", par = "lambda"))
APSD_lam <- aov_car(value ~ R + L + Error(id/M), 
                    data = data_split_sens(df_long, measure = "postsd",
                                           statistic = "avg", par = "lambda"))
AECR_lam <- aov_car(value ~ R + L + Error(id/M), 
                    data = data_split_sens(df_long, measure = "ECR",
                                           statistic = "avg", par = "lambda"))
Tucker_lam <- aov_car(value ~ R + L + Error(id/M), 
                      data = data_split_sens(df_long, measure = "Tucker",
                                             statistic = "avg", par = "global"))

# Diagnostics (correlation parameters)
Rhat_Rho <- aov_car(value ~ R + L + Error(id/M),
                    data = data_split_sens(df_long, measure = "rhat",
                                           statistic = "avg", par = "Rho"))
Div_mod  <- aov_car(value ~ R + L + Error(id/M),
                    data = data_split_sens(df_long, measure = "divergence",
                                           statistic = "avg", par = NA))
TD_mod   <- aov_car(value ~ R + L + Error(id/M),
                    data = data_split_sens(df_long, measure = "treedepth",
                                           statistic = "avg", par = NA))
ESS_Rho  <- aov_car(value ~ R + L + Error(id/M),
                    data = data_split_sens(df_long, measure = "ess_pct",
                                           statistic = "avg", par = "Rho"))

# Diagnostics (factor loading parameters)
Rhat_lam <- aov_car(value ~ R + L + Error(id/M),
                    data = data_split_sens(df_long, measure = "rhat",
                                           statistic = "avg", par = "lambda"))
ESS_lam  <- aov_car(value ~ R + L + Error(id/M),
                    data = data_split_sens(df_long, measure = "ess_pct",
                                           statistic = "avg", par = "lambda"))

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: ANOVA effect sizes (omega-squared)
# ─────────────────────────────────────────────────────────────────────────────

ANOVA_es <- data.frame(
  Predictor = omega_squared(ASRMR_Rho)[,1],
  ASRMR     = round(omega_squared(ASRMR_Rho)[,2], 2),
  AAB       = round(omega_squared(AAB_Rho)[,2], 2),
  APSD      = round(omega_squared(APSD_Rho)[,2], 2),
  AECR      = round(omega_squared(AECR_Rho)[,2], 2),
  ARMSE     = round(omega_squared(ARMSE_lam)[,2], 2),
  AAB       = round(omega_squared(AAB_lam)[,2], 2),
  APSD      = round(omega_squared(APSD_lam)[,2], 2),
  AECR      = round(omega_squared(AECR_lam)[,2], 2),
  Phi       = round(omega_squared(Tucker_lam)[,2], 2)
)

clean_colnames <- c("Predictor", "ASRMR", "AAB", "APSD", "AECR", 
                    "ARMSE", "AAB", "APSD", "AECR", "$\\phi$")

ANOVA_tbl <- ANOVA_es |> 
  mutate(across(where(is.numeric), function(x) {
    txt_val <- sprintf("%.2f", x)
    ifelse(x >= 0.14, paste0("<strong>", txt_val, "</strong>"), txt_val)
  }))

kbl(ANOVA_tbl, booktabs = TRUE, escape = FALSE,        
    col.names = clean_colnames, longtable = TRUE, linesep = "", 
    align = c("l", rep("c", 9)),
    caption = "Effect sizes ($\\omega^2_p$) for Correlations and Factor Loadings (Sensitivity)") %>%
  add_header_above(c(" " = 1, "Correlations" = 4, "Factor loadings" = 5), bold = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"), font_size = 10) %>% 
  footnote(general = c("Effect sizes >= .14 in bold. R = Reliability. L = True factor loadings. M = Estimation model."),
           general_title = "Note.", footnote_as_chunk = TRUE, 
           threeparttable = TRUE, escape = FALSE)

# Save LaTeX version
ANOVA_tbl_latex <- ANOVA_es |>
  mutate(across(where(is.numeric), function(x) {
    txt_val <- sprintf("%.2f", x)
    ifelse(x >= 0.14, paste0("\\textbf{", txt_val, "}"), txt_val)
  }))
save(ANOVA_tbl_latex, file = "Results/Rdata/Supplementary Material/ANOVA_table_sensitivity.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Emmeans table: correlations (bias measures)
# ─────────────────────────────────────────────────────────────────────────────

emmeans_rho <- build_measure_table(
  anova_fits = list(ASRMR_Rho, AAB_Rho, APSD_Rho, AECR_Rho),
  measure_names = c("ASRMR", "AAB", "APSD", "AECR"),
  highlight_crits = c("min", "min", "min", "0.95"),
  format = "html"
)

# Column names: G, SL, EG repeated per measure
colnames_rho <- c("", rep(c("G", "SL", "EG"), times = 4))
row_info_rho <- table(emmeans_rho$Factor_Name)[unique(emmeans_rho$Factor_Name)]

emmeans_rho[, 2:ncol(emmeans_rho)] %>%
  kable(format = "html", escape = FALSE,
        align = paste0("l", paste(rep("c", 12), collapse = "")),
        col.names = colnames_rho,
        caption = "Estimated Marginal Means: Correlations (Sensitivity)") %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "ASRMR" = 3, "AAB" = 3, "APSD" = 3, "AECR" = 3)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info_rho)) {
      x <- pack_rows(x, names(row_info_rho)[i], start_row, start_row + row_info_rho[i] - 1)
      start_row <- start_row + row_info_rho[i]
    }
    x
  }}() %>% 
  add_footnote(notation = "none", 
               label = c("Note: G = Gaussian HFM; SL = Shifted-Lognormal HFM; EG = Ex-Gaussian HFM.",
                         "Lowest ASRMR, AAB, APSD and ECR closest to .95 in bold."), 
               threeparttable = TRUE)

# Save LaTeX version
emmeans_rho_latex <- build_measure_table(
  anova_fits = list(ASRMR_Rho, AAB_Rho, APSD_Rho, AECR_Rho),
  measure_names = c("ASRMR", "AAB", "APSD", "AECR"),
  highlight_crits = c("min", "min", "min", "0.95"),
  format = "latex"
)
save(emmeans_rho_latex, file = "Results/Rdata/Supplementary Material/emmeans_rho_sensitivity.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Emmeans table: factor loadings (bias measures)
# ─────────────────────────────────────────────────────────────────────────────

emmeans_lam <- build_measure_table(
  anova_fits = list(ARMSE_lam, AAB_lam, APSD_lam, Tucker_lam, AECR_lam),
  measure_names = c("ARMSE", "AAB", "APSD", "Phi", "AECR"),
  highlight_crits = c("min", "min", "min", "max", "0.95"),
  format = "html"
)

colnames_lam <- c("", rep(c("G", "SL", "EG"), times = 5))
row_info_lam <- table(emmeans_lam$Factor_Name)[unique(emmeans_lam$Factor_Name)]

emmeans_lam[, -1] %>%
  kable(format = "html", escape = FALSE,
        align = paste0("l", paste(rep("c", 15), collapse = "")),
        col.names = colnames_lam,
        caption = "Estimated Marginal Means: Factor Loadings (Sensitivity)") %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "ARMSE" = 3, "AAB" = 3, "APSD" = 3, 
                     "$\\phi$" = 3, "AECR" = 3)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info_lam)) {
      x <- pack_rows(x, names(row_info_lam)[i], start_row, start_row + row_info_lam[i] - 1)
      start_row <- start_row + row_info_lam[i]
    }
    x
  }}() %>% 
  add_footnote(notation = "none", 
               label = c("Note: G = Gaussian HFM; SL = Shifted-Lognormal HFM; EG = Ex-Gaussian HFM.",
                         "Lowest ARMSE, AAB, APSD, highest phi, and ECR closest to .95 in bold."), 
               threeparttable = TRUE, escape = FALSE)

# Save LaTeX version
emmeans_lam_latex <- build_measure_table(
  anova_fits = list(ARMSE_lam, AAB_lam, APSD_lam, Tucker_lam, AECR_lam),
  measure_names = c("ARMSE", "AAB", "APSD", "Phi", "AECR"),
  highlight_crits = c("min", "min", "min", "max", "0.95"),
  format = "latex"
)
save(emmeans_lam_latex, file = "Results/Rdata/Supplementary Material/emmeans_lambda_sensitivity.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Emmeans table: diagnostics for correlations
# ─────────────────────────────────────────────────────────────────────────────

# Rhat and ESS% are Rho-specific; Div% and MaxTD% are model-level
emmeans_diag_rho <- build_measure_table(
  anova_fits = list(Rhat_Rho, Div_mod, TD_mod, ESS_Rho),
  measure_names = c("Rhat", "Div", "MaxTD", "ESS"),
  highlight_crits = c("min", "min", "min", "max"),
  format = "html"
)

colnames_diag_rho <- c("", rep(c("G", "SL", "EG"), times = 4))
row_info_diag_rho <- table(emmeans_diag_rho$Factor_Name)[unique(emmeans_diag_rho$Factor_Name)]

emmeans_diag_rho[, 2:ncol(emmeans_diag_rho)] %>%
  kable(format = "html", escape = FALSE,
        align = paste0("l", paste(rep("c", 12), collapse = "")),
        col.names = colnames_diag_rho,
        caption = "Estimated Marginal Means: Correlation Diagnostics (Sensitivity)") %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "$\\hat{R}$" = 3, "Div%" = 3, 
                     "MaxTD%" = 3, "ESS%" = 3)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info_diag_rho)) {
      x <- pack_rows(x, names(row_info_diag_rho)[i], start_row, start_row + row_info_diag_rho[i] - 1)
      start_row <- start_row + row_info_diag_rho[i]
    }
    x
  }}() %>% 
  add_footnote(notation = "none", 
               label = c("Note: G = Gaussian HFM; SL = Shifted-Lognormal HFM; EG = Ex-Gaussian HFM.",
                         "Rhat and ESS% for correlation parameters. Div% and MaxTD% are model-level.",
                         "ESS% = avg ESS / 2000 draws x 100."), 
               threeparttable = TRUE)

emmeans_diag_rho_latex <- build_measure_table(
  anova_fits = list(Rhat_Rho, Div_mod, TD_mod, ESS_Rho),
  measure_names = c("Rhat", "Div", "MaxTD", "ESS"),
  highlight_crits = c("min", "min", "min", "max"),
  format = "latex"
)
save(emmeans_diag_rho_latex, file = "Results/Rdata/Supplementary Material/emmeans_diag_rho_sensitivity.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11: Emmeans table: diagnostics for factor loadings
# ─────────────────────────────────────────────────────────────────────────────

# Rhat and ESS% are lambda-specific; Div% and MaxTD% are model-level (same as above)
emmeans_diag_lam <- build_measure_table(
  anova_fits = list(Rhat_lam, Div_mod, TD_mod, ESS_lam),
  measure_names = c("Rhat", "Div", "MaxTD", "ESS"),
  highlight_crits = c("min", "min", "min", "max"),
  format = "html"
)

colnames_diag_lam <- c("", rep(c("G", "SL", "EG"), times = 4))
row_info_diag_lam <- table(emmeans_diag_lam$Factor_Name)[unique(emmeans_diag_lam$Factor_Name)]

emmeans_diag_lam[, -1] %>%
  kable(format = "html", escape = FALSE,
        align = paste0("l", paste(rep("c", 12), collapse = "")),
        col.names = colnames_diag_lam,
        caption = "Estimated Marginal Means: Factor Loading Diagnostics (Sensitivity)") %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "$\\hat{R}$" = 3, "Div%" = 3, 
                     "MaxTD%" = 3, "ESS%" = 3)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info_diag_lam)) {
      x <- pack_rows(x, names(row_info_diag_lam)[i], start_row, start_row + row_info_diag_lam[i] - 1)
      start_row <- start_row + row_info_diag_lam[i]
    }
    x
  }}() %>% 
  add_footnote(notation = "none", 
               label = c("Note: G = Gaussian HFM; SL = Shifted-Lognormal HFM; EG = Ex-Gaussian HFM.",
                         "Rhat and ESS% for factor loading parameters. Div% and MaxTD% are model-level.",
                         "ESS% = avg ESS / 2000 draws x 100."), 
               threeparttable = TRUE)

emmeans_diag_lam_latex <- build_measure_table(
  anova_fits = list(Rhat_lam, Div_mod, TD_mod, ESS_lam),
  measure_names = c("Rhat", "Div", "MaxTD", "ESS"),
  highlight_crits = c("min", "min", "min", "max"),
  format = "latex"
)
save(emmeans_diag_lam_latex, file = "Results/Rdata/Supplementary Material/emmeans_diag_lam_sensitivity.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 12: Emmeans table: combined diagnostics (across Rho and lambda)
# ─────────────────────────────────────────────────────────────────────────────

# Average Rhat and ESS% across both Rho and lambda parameters
combined_diag <- df_long |>
  filter(measure %in% c("rhat", "ess_pct") & statistic == "avg" &
         parameter %in% c("Rho", "lambda") &
         model_pars %in% c("gaussian", "shlognormal", "exgaussian")) |>
  group_by(reliab, std_lambda_j, measure, statistic, model_pars) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(parameter = "combined")

df_long <- bind_rows(df_long, combined_diag)

# ANOVAs for combined diagnostics
Rhat_comb <- aov_car(value ~ R + L + Error(id/M),
                     data = data_split_sens(df_long, measure = "rhat",
                                            statistic = "avg", par = "combined"))
ESS_comb  <- aov_car(value ~ R + L + Error(id/M),
                     data = data_split_sens(df_long, measure = "ess_pct",
                                            statistic = "avg", par = "combined"))

# Build emmeans table: Rhat, Div%, MaxTD%, ESS% (combined across parameters)
emmeans_diag_comb <- build_measure_table(
  anova_fits = list(Rhat_comb, Div_mod, TD_mod, ESS_comb),
  measure_names = c("Rhat", "Div", "MaxTD", "ESS"),
  highlight_crits = c("min", "min", "min", "max"),
  format = "html"
)

colnames_diag_comb <- c("", rep(c("G", "SL", "EG"), times = 4))
row_info_diag_comb <- table(emmeans_diag_comb$Factor_Name)[unique(emmeans_diag_comb$Factor_Name)]

emmeans_diag_comb[, 2:ncol(emmeans_diag_comb)] %>%
  kable(format = "html", escape = FALSE,
        align = paste0("l", paste(rep("c", 12), collapse = "")),
        col.names = colnames_diag_comb,
        caption = "Estimated Marginal Means: Combined Diagnostics (Sensitivity)") %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "$\\hat{R}$" = 3, "Div%" = 3,
                     "MaxTD%" = 3, "ESS%" = 3)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info_diag_comb)) {
      x <- pack_rows(x, names(row_info_diag_comb)[i], start_row, start_row + row_info_diag_comb[i] - 1)
      start_row <- start_row + row_info_diag_comb[i]
    }
    x
  }}() %>%
  add_footnote(notation = "none",
               label = c("Note: G = Gaussian HFM; SL = Shifted-Lognormal HFM; EG = Ex-Gaussian HFM.",
                         "Rhat and ESS% averaged across correlation and factor loading parameters.",
                         "Div% and MaxTD% are model-level. ESS% = avg ESS / 2000 draws x 100."),
               threeparttable = TRUE)

# Save LaTeX version
emmeans_diag_comb_latex <- build_measure_table(
  anova_fits = list(Rhat_comb, Div_mod, TD_mod, ESS_comb),
  measure_names = c("Rhat", "Div", "MaxTD", "ESS"),
  highlight_crits = c("min", "min", "min", "max"),
  format = "latex"
)
save(emmeans_diag_comb_latex, file = "Results/Rdata/Supplementary Material/emmeans_diag_combined_sensitivity.rdata")

# ─────────────────────────────────────────────────────────────────────────────
