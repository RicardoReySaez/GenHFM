# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : simulation_analysis.R                                      ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 16-10-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# ANOVA meta-models. Notice that the variables are:
#   1. I: Sample size (100 or 200)
#   2. LF: Latent factors (unidimensional or 2-factors, both CFA or EFA)
#   3. DGM: Data-generation model (ex-gaussian or shifted-lognormal)
#   4. R: Reliability (0.3, 0.5 or 0.7)
#   5. L: Standardized factor loadings (0.3, 0.5 or 0.7)
#   6. M: Estimation method (aggregated effects, spearman, HFM or GenHFM)
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(dplyr)
library(tidyr)
library(stringr)
library(afex)
library(emmeans)
library(effectsize)
library(kableExtra)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: user defined functions
# ─────────────────────────────────────────────────────────────────────────────

# Fast function to generate data frame per measute
data_split <- function(df_long, par, measure, statistic) {
  # Split data frame
  df <- as.data.frame(df_long[which(df_long$measure == measure &
                                      df_long$parameter == par &
                                      df_long$statistic == statistic),])
  # Change column names
  colnames(df) <- c("DGM", "LF", "I", "R", "L", "measure", 
                    "statistic", "parameter", "M", "value")
  # Specify columns as factors
  for(j in 1:(ncol(df)-1)) {df[,j] <- as.factor(df[,j])}
  
  # Add an ID column for estimation method
  if(par == "Rho") {
    df$id <- rep(1:(nrow(df)/4), each = 4)
  } 
  if(par == "lambda" | par == "global") {
    df$id <- rep(1:(nrow(df)/2), each = 2)
  }
  # Return final data frame
  return(df)
}

# Marginal means per ANOVA model
marginal_means <- function(anova_fit, highlight_crit = "min", format = "latex") {
  
  # Emmeans per factor
  emmeans_list <- list(
    I   = emmeans(anova_fit, specs = ~ I + M),
    LF  = emmeans(anova_fit, specs = ~ LF + M),
    DGM = emmeans(anova_fit, specs = ~ DGM + M),
    R   = emmeans(anova_fit, specs = ~ R + M),
    L   = emmeans(anova_fit, specs = ~ L + M)
  )
  
  # Initial process
  table_data_raw <- do.call(rbind, emmeans_list) |> 
    as.data.frame() |> 
    mutate(across(c(I, LF, DGM, R, L), as.character)) |> 
    pivot_longer(cols = c(I, LF, DGM, R, L), names_to = "Factor_Name", values_to = "Level") |> 
    filter(Level != ".")
  
  # Check if this is for lambda or correlations
  orden_logico <- c("Naive", "Spearman", "gaussian", "skewed")
  niveles_presentes <- unique(as.character(table_data_raw$M))
  target_cols <- orden_logico[orden_logico %in% niveles_presentes]
  
  # Wide data frame
  table_data <- table_data_raw |> 
    mutate(M = factor(M, levels = target_cols)) |> 
    select(Factor_Name, Level, M, emmean) |> 
    pivot_wider(names_from = M, values_from = emmean) |> 
    mutate(Level = recode(Level, 
                          "2_CFA" = "CFA", 
                          "2_EFA" = "EFA", 
                          "shifted.lognormal" = "sh-lognormal")) |>
    mutate(Level = factor(Level, levels = c(
      "100", "200",                       
      "uni", "CFA", "EFA",                
      "exgaussian", "sh-lognormal",       
      "0.3", "0.5", "0.7"                 
    ))) |>
    mutate(Factor_Name = factor(Factor_Name, levels = c("I", "LF", "DGM", "R", "L"))) |>
    arrange(Factor_Name, Level) |>
    mutate(Factor_Name = recode(
      Factor_Name, 
      "I" = "Sample size", 
      "LF" = "Latent model", 
      "DGM" = "True model", 
      "R" = "Reliability", 
      "L" = "True $\\Lambda$")) |>
    select(Factor_Name, Level, all_of(target_cols)) |>
    mutate(across(all_of(target_cols), \(x) sprintf("%.3f", round(x, 3)))) |>
    mutate(across(everything(), as.character))
  
  # Control for bold letters
  mat_values <- as.matrix(table_data[, target_cols])
  
  # Function with several criteria
  find_target <- function(row_values) {
    num_vals <- as.numeric(row_values)
    if (all(is.na(num_vals))) return(NULL)
    
    if (highlight_crit == "min") return(which.min(num_vals))
    if (highlight_crit == "max") return(which.max(num_vals))
    if (highlight_crit == "0.95") return(which.min(abs(num_vals - 0.95)))
    return(NULL)
  }
  
  # Bold leters for latex or html
  for (i in 1:nrow(table_data)) {
    row_vals <- mat_values[i, ]
    target_idx <- find_target(row_vals)
    
    if (!is.null(target_idx)) {
      val_text <- table_data[i, target_cols[target_idx]]
      
      formatted_val <- if (format == "latex") {
        paste0("\\textbf{", val_text, "}")
      } else {
        paste0("<b>", val_text, "</b>")
      }
      
      table_data[i, target_cols[target_idx]] <- formatted_val
    }
  }
  
  return(table_data)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Prepare final data frame
# ─────────────────────────────────────────────────────────────────────────────

# Load aggregated simulation results
simres <- readRDS("Results/Rdata/Simulation study/aggregated_simulation_results.rds")

df_long <- simres[, 1:(which(colnames(simres) == "REPLICATIONS")-1)] |>
  pivot_longer(
    cols = -c(Model, M, I, reliab, std_lambda_j),
    names_to = "name",
    values_to = "value"
  ) |>
  mutate(
    # Identify measure (keep "Tuckey")
    measure = case_when(
      str_starts(name, "Tucker")      ~ "Tucker",
      str_starts(name, "skew_first")  ~ "skew_first",
      str_starts(name, "prob_skew")   ~ "prob_skew",
      TRUE ~ str_extract(
        name,
        paste(
          c("abias","arbias","stdbias","rmse","srmr","postsd","empirical",
            "rhat","essb","esst","ECR","power","loo","divergence",
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
      str_starts(name, "skew_first|prob_skew") ~ "global",
      TRUE ~ str_extract(
        name,
        paste(
          c("Rho","lambda","reliab",
            "mu_alpha","mu_theta","mu_sigma","mu_tau","mu_delta",
            "sd_alpha","sd_theta","sd_sigma","sd_tau","sd_delta",
            "crossload","Phi","elpd",
            "both_weak","both_strong","skew_weak","skew_strong"),
          collapse = "|"
        )
      )
    ),
    
    # Model variants from suffixes
    model_pars = case_when(
      parameter %in% c("mu_alpha","mu_theta","mu_sigma","mu_tau","mu_delta",
                       "sd_alpha","sd_theta","sd_sigma","sd_tau","sd_delta") ~ "skewed",
      str_detect(name, "_gauss") ~ "gaussian",
      str_detect(name, "_mm")    ~ "Naive",
      str_detect(name, "_sp")    ~ "Spearman",
      str_detect(name, "_skew")  ~ "skewed",
      str_starts(name, "skew_first|prob_skew") ~ "comparison",
      TRUE ~ NA_character_
    )
  )

# Tuckey: average f1/f2 within each statistic; if only f1 (or M == "uni"), we use f1
tuckey_df <- df_long |>
  filter(measure == "Tucker") |>
  group_by(Model, M, I, reliab, std_lambda_j, measure, statistic, parameter, model_pars) |>
  summarise(
    value = {
      if (first(M) == "uni") {
        mean(value[f_idx %in% c(NA, "f1")], na.rm = TRUE)
      } else {
        has_f1 <- any(f_idx == "f1", na.rm = TRUE)
        has_f2 <- any(f_idx == "f2", na.rm = TRUE)
        if (has_f1 && has_f2) {
          mean(value[f_idx %in% c("f1","f2")], na.rm = TRUE)
        } else if (has_f1) {
          mean(value[f_idx == "f1"], na.rm = TRUE)
        } else if (has_f2) {
          mean(value[f_idx == "f2"], na.rm = TRUE)
        } else {
          mean(value, na.rm = TRUE)
        }
      }
    },
    .groups = "drop"
  )

# Others: standard averaging (already separated by 'statistic' if present)
others_df <- df_long |>
  filter(measure != "Tucker") |>
  group_by(Model, M, I, reliab, std_lambda_j, measure, statistic, parameter, model_pars) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

df_long <- bind_rows(tuckey_df, others_df)

# Save this final df
saveRDS(object = df_long, file = "Results/Rdata/Simulation study/simulation_results_longdata.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: ANOVA meta-models for correlations and factor loadings
# ─────────────────────────────────────────────────────────────────────────────

# Correlations meta-models: Average SRMR
ASRMR_Rho <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                     data = data_split(df_long = df_long, measure = "srmr",
                                       statistic = "avg", par = "Rho"))

# Correlations meta-models: Average absolute bias
AAB_Rho <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                   data = data_split(df_long = df_long, measure = "abias",
                                     statistic = "avg", par = "Rho"))

# Correlations meta-models: Average uncertainty
APSD_Rho <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                     data = data_split(df_long = df_long, measure = "postsd",
                                       statistic = "avg", par = "Rho"))

# Correlations meta-models: Average Empirical Coverage Rate
AECR_Rho <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                    data = data_split(df_long = df_long, measure = "ECR",
                                      statistic = "avg", par = "Rho"))

# Factor loadings meta-models: Average SRMR
ARMSE_lam <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                     data = data_split(df_long = df_long, measure = "rmse",
                                       statistic = "avg", par = "lambda"))

# Factor loadings meta-models: Average absolute bias
AAB_lam <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                   data = data_split(df_long = df_long, measure = "abias",
                                     statistic = "avg", par = "lambda"))

# Factor loadings meta-models: Average uncertainty
APSD_lam <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                    data = data_split(df_long = df_long, measure = "postsd",
                                      statistic = "avg", par = "lambda"))

# Factor loadings meta-models: Average Empirical Coverage Rate
AECR_lam <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                    data = data_split(df_long = df_long, measure = "ECR",
                                      statistic = "avg", par = "lambda"))

# Factor loadings meta-models: Average Empirical Coverage Rate
Tucker_lam <- aov_car(value ~ (I + LF + DGM + R + L)^4 + Error(id/M), 
                      data = data_split(df_long = df_long, measure = "Tucker",
                                        statistic = "avg", par = "global"))

# ANOVA results based on effect sizes (omega-squared)
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

# Set new column names for ANOVA table
clean_colnames <- c("Predictor", "ASRMR", "AAB", "APSD", "AECR", 
                    "ARMSE", "AAB", "APSD", "AECR", "$\\phi$")

# Set in bold omega-squared values above or equal to 0.14
ANOVA_tbl <- ANOVA_es |> 
  mutate(across(where(is.numeric), function(x) {
    txt_val <- sprintf("%.2f", x)
    ifelse(x >= 0.14, 
           paste0("<strong>", txt_val, "</strong>"),
           txt_val)
  }))

# Final table results
kbl(ANOVA_tbl, 
    booktabs = TRUE, 
    escape = FALSE,        
    col.names = clean_colnames, 
    longtable = TRUE,
    linesep = "", 
    align = c("l", rep("c", 9)),
    caption = "Effect sizes ($\\omega^2_p$) for Correlations and Factor Loadings") %>%
  add_header_above(c(" " = 1, 
                     "Correlations" = 4, 
                     "Factor loadings" = 5), 
                   bold = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"), 
                font_size = 10) %>% 
  footnote(general = c("Effect sizes equal or higuer than .14 are marked in bold. I = Sample size (100 or 200). LF = Latent factors (unidimensional, two-factor CFA and two-factor CFA). DGM = Data-Generation Model (ex-Gaussian or Shifted-lognormal). R = Reliability (0.3, 0.5 or 0.7). L = Standardized factor loadings (0.3, 0.5 or 0.7). M = Estimation method (Gaussian HFM or True GenHFM). ASRMR = Average Standardized Root Mean Residuals. AAB = Average Absolute Bias. APSD = Average Posterior Standard Deviation (Average Standard Error in frequentist correlations). AECR = Average Empirical Coverage Rate. $\\phi$ = Tucker Congruence Index."),
           general_title = "Note.",
           footnote_as_chunk = TRUE, 
           threeparttable = TRUE,
           escape = FALSE)

# Prepare the same ANOVA table but just for LaTeX and save it
ANOVA_tbl <- ANOVA_es |>
  mutate(across(where(is.numeric), function(x) {
    txt_val <- sprintf("%.2f", x)
    ifelse(x >= 0.14,
           paste0("\\textbf{", txt_val, "}"),
           txt_val)
  }))

save(ANOVA_tbl, file = "Results/Rdata/Tables/ANOVA_table.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Estimated marginal means per model: true correlations
# ─────────────────────────────────────────────────────────────────────────────

# Correlations
marginal_means_rho <- cbind(
  marginal_means(anova_fit = ASRMR_Rho, highlight_crit = "min", format = "html"),
  marginal_means(anova_fit = AAB_Rho, highlight_crit = "min", format = "html")[-c(1:2)],
  marginal_means(anova_fit = APSD_Rho, highlight_crit = "min", format = "html")[-c(1:2)],
  marginal_means(anova_fit = AECR_Rho, highlight_crit = "0.95", format = "html")[-c(1:2)])

# 1. Definimos los nombres de las columnas para la tabla final
# Son 17 columnas en total (Level + 16 de datos)
colnames_final_rho <- c("", rep(c("A", "S", "G", "T"), times = 4))

# 2. Calculamos la agrupación de filas
# unique() asegura que mantenemos el orden de aparición de los factores
row_info_rho <- table(marginal_means_rho$Factor_Name)[unique(marginal_means_rho$Factor_Name)]

# 3. Generamos la tabla HTML
# Usamos el dataframe base pero seleccionamos por posición para evitar el error de nombres
marginal_means_rho[, 2:18] %>%  # Seleccionamos desde 'Level' hasta el final por posición
  kable(
    format = "html", 
    escape = FALSE,          
    align = "lcccccccccccccccc",
    col.names = colnames_final_rho,
    caption = "Table 3: Estimated Marginal Means for Fit and Uncertainty Indices"
  ) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(
    " " = 1, 
    "ASRMR" = 4, 
    "AAB" = 4, 
    "APSD/ASE" = 4, 
    "AECR" = 4
  )) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info_rho)) {
      x <- pack_rows(x, names(row_info_rho)[i], start_row, start_row + row_info_rho[i] - 1)
      start_row <- start_row + row_info_rho[i]
    }
    x
  }}() %>% 
  add_footnote(notation = "none", label = c("Note: A = Aggregated scores; S = Spearman correction for attenuation; G = Gaussian HFM; T = True GenHFM.", "The lowest SRMR, AAB, and APSD values, and the ECR value closest to .95, are highlighted in bold."), threeparttable = TRUE)

# Save Estimated Marginal Means for correlations in LaTeX format
marginal_means_rho_latex <- cbind(
  marginal_means(anova_fit = ASRMR_Rho, highlight_crit = "min", format = "latex"),
  marginal_means(anova_fit = AAB_Rho, highlight_crit = "min", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = APSD_Rho, highlight_crit = "min", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = AECR_Rho, highlight_crit = "0.95", format = "latex")[-c(1:2)])
save(marginal_means_rho_latex, file = "Results/Rdata/Tables/emmeans_rho.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Estimated marginal means per model: factor loadings
# ─────────────────────────────────────────────────────────────────────────────

# Correlations
marginal_means_lam <- cbind(
  marginal_means(anova_fit = ARMSE_lam, highlight_crit = "min", format = "html"),
  marginal_means(anova_fit = AAB_lam, highlight_crit = "min", format = "html")[-c(1:2)],
  marginal_means(anova_fit = APSD_lam, highlight_crit = "min", format = "html")[-c(1:2)],
  marginal_means(anova_fit = Tucker_lam, highlight_crit = "max", format = "html")[-c(1:2)],
  marginal_means(anova_fit = AECR_lam, highlight_crit = "0.95", format = "html")[-c(1:2)])

# 1. Definimos los nombres de las columnas para la tabla final
# Son 17 columnas en total (Level + 16 de datos)
colnames_final_lam <- c("", rep(c("Gaussian", "True"), times = 5))

# 2. Calculamos la agrupación de filas
# unique() asegura que mantenemos el orden de aparición de los factores
row_info_lam <- table(marginal_means_lam$Factor_Name)[unique(marginal_means_lam$Factor_Name)]

# 3. Generamos la tabla HTML
# Usamos el dataframe base pero seleccionamos por posición para evitar el error de nombres
marginal_means_lam[, -1] %>%  # Seleccionamos desde 'Level' hasta el final por posición
  kable(
    format = "html", 
    escape = FALSE,          
    align = "lcccccccc",
    col.names = colnames_final_lam,
    caption = "Estimated marginal means for factor loadings bias measures"
  ) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(
    " " = 1, 
    "ARMSE" = 2, 
    "AAB" = 2, 
    "APSD" = 2, 
    "Tucker's $\\phi$" = 2,
    "AECR" = 2
  )) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info_lam)) {
      x <- pack_rows(x, names(row_info_lam)[i], start_row, start_row + row_info_lam[i] - 1)
      start_row <- start_row + row_info_lam[i]
    }
    x
  }}() %>% 
  add_footnote(notation = "none", label = "Note: The lowest RMSE, AAB, and APSD values, the highest Tucker's $\\\\phi$ values, and the ECR value closest to .95 are in bold.", threeparttable = TRUE, escape = FALSE)

# Save Estimated Marginal Means for factor loadings in LaTeX format
marginal_means_lam_latex <- cbind(
  marginal_means(anova_fit = ARMSE_lam, highlight_crit = "min", format = "latex"),
  marginal_means(anova_fit = AAB_lam, highlight_crit = "min", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = APSD_lam, highlight_crit = "min", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = Tucker_lam, highlight_crit = "max", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = AECR_lam, highlight_crit = "0.95", format = "latex")[-c(1:2)])
save(marginal_means_lam_latex, file = "Results/Rdata/Tables/emmeans_lambda.rdata")

# ─────────────────────────────────────────────────────────────────────────────
