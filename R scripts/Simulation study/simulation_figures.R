# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : simulation_figures.R                                       ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 20-10-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Script to generate the figures used in the manuscript
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

# Libraries necessary for the script to function
library(SimDesign)
library(MASS)
library(dplyr)
library(ggplot2)
library(ggdist)
library(bayesplot)
library(patchwork)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Load simulation functions
source("R scripts/R functions/simulation_functions.R")

# Load simulation results
full_simres <- readRDS("Results/Rdata/Simulation study/raw_simulation_results.rds")

# Load metaparameters
load("Results/Rdata/Meta-parameters/meta_parameters.rdata")
fixed_objects <- list(meta_parameters = meta_parameters)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: compute results data frame
# ─────────────────────────────────────────────────────────────────────────────

# Simulation design
Design <- createDesign(Model = c("shifted.lognormal", "exgaussian"),
                       M = c("uni", "2_EFA","2_CFA"),
                       I = c(100, 200),
                       reliab = c(.3, .5, .7),
                       std_lambda_j = c(.3, .5, .7))

# Empty data frames per condition
temp_data_lists <- vector(mode = "list", length = nrow(Design))

# Generate results data frames
for(i in 1:nrow(Design)) {
  # Small progress bar
  cat(paste("\rCondition", i, "of", nrow(Design)))
  
  # Population parameters in this condition
  pop_params <- cond_parameters(condition = Design[i,], fixed_objects = fixed_objects)
  
  # Results in this condition
  results <- as.data.frame(merge(full_simres, Design[i,], by = names(Design[i,])))
  
  # ───────────────────────────────────── #
  #    bias in all parameter estimates    #
  # ───────────────────────────────────── #
  
  # Global summaries per parameter
  mu_alpha <- nonagg_summaries(results = results, parameter_name = "mu_alpha", parameter = pop_params$mu_alpha)
  mu_theta <- nonagg_summaries(results = results, parameter_name = "mu_theta", parameter = pop_params$mu_theta)
  mu_sigma <- nonagg_summaries(results = results, parameter_name = "mu_sigma", parameter = pop_params$mu_sigma)
  sd_alpha <- nonagg_summaries(results = results, parameter_name = "sd_alpha", parameter = pop_params$sd_alpha)
  sd_theta <- nonagg_summaries(results = results, parameter_name = "sd_theta", parameter = pop_params$sd_theta)
  sd_sigma <- nonagg_summaries(results = results, parameter_name = "sd_sigma", parameter = pop_params$sd_sigma)
  
  # Summaries for correlation matrices
  Rho_gauss <- nonagg_summaries(results = results, parameter_name = "Rho_gauss", parameter = pop_params$R_theta, Rho = TRUE)
  Rho_skew <- nonagg_summaries(results = results, parameter_name = "Rho_skew", parameter = pop_params$R_theta, Rho = TRUE)
  Rho_mm <- nonagg_summaries(results = results, parameter_name = "Rho_mm", parameter = pop_params$R_theta, method_of_moments = TRUE)
  Rho_sp <- nonagg_summaries(results = results, parameter_name = "Rho_sp", parameter = pop_params$R_theta, method_of_moments = TRUE)
  
  # Summaries for reliability
  reliab_gauss <- nonagg_summaries(results = results, parameter_name = "reliab_gauss", parameter = pop_params$reliab)
  reliab_skew  <- nonagg_summaries(results = results, parameter_name = "reliab_skew", parameter = pop_params$reliab)
  
  # Empty values (depends on latent structure)
  Phi_gauss <- nonagg_summaries(parameter_name = "Phi_gauss", empty = TRUE, Phi.lv = TRUE)
  Phi_skew <- nonagg_summaries(parameter_name = "Phi_skew", empty = TRUE, Phi.lv = TRUE)
  crossload_gauss <- nonagg_summaries(parameter_name = "crossload_gauss", empty = TRUE)
  crossload_skew <- nonagg_summaries(parameter_name = "crossload_skew", empty = TRUE)
  
  # Empty values (depends on specific model)
  mu_delta <- nonagg_summaries(parameter_name = "mu_delta", empty = TRUE)
  sd_delta <- nonagg_summaries(parameter_name = "sd_delta", empty = TRUE)
  mu_tau <- nonagg_summaries(parameter_name = "mu_tau", empty = TRUE)
  sd_tau <- nonagg_summaries(parameter_name = "sd_tau", empty = TRUE)
  
  # Unidimensional settings
  if(Design[i, "M"] == "uni"){
    # Gaussian and skewed lambdas
    lambda_gauss <- nonagg_summaries(results = results, 
                                     parameter_name = "lambda_gauss", 
                                     parameter = c(pop_params$lambda_std))
    lambda_skew <- nonagg_summaries(results = results, 
                                    parameter_name = "lambda_skew", 
                                    parameter = c(pop_params$lambda_std))
  }
  
  # Exploratory factor analysis
  if(Design[i, "M"] == "2_EFA"){
    # Main factor loadings
    lambda_gauss <- nonagg_summaries(results = results, 
                                     parameter_name = "lambda_gauss", 
                                     parameter = c(pop_params$lambda_std[1:3,1], 
                                                   pop_params$lambda_std[4:6,2]))
    lambda_skew  <- nonagg_summaries(results = results, 
                                     parameter_name = "lambda_skew", 
                                     parameter = c(pop_params$lambda_std[1:3,1], 
                                                   pop_params$lambda_std[4:6,2]))
    # Cross-loadings
    crossload_gauss <- nonagg_summaries(results = results, 
                                        parameter_name = "crossload_gauss", 
                                        parameter = c(pop_params$lambda_std[4:6,1], 
                                                      pop_params$lambda_std[1:3,2]))
    crossload_skew  <- nonagg_summaries(results = results, 
                                        parameter_name = "crossload_skew", 
                                        parameter = c(pop_params$lambda_std[4:6,1], 
                                                      pop_params$lambda_std[1:3,2]))
  }
  
  # Exploratory factor analysis
  if(Design[i, "M"] == "2_CFA"){
    # Factor loadings
    lambda_gauss <- nonagg_summaries(results = results, 
                                     parameter_name = "lambda_gauss", 
                                     parameter = c(pop_params$lambda_std[1:3,1], 
                                                   pop_params$lambda_std[4:6,2]))
    lambda_skew  <- nonagg_summaries(results = results, 
                                     parameter_name = "lambda_skew", 
                                     parameter = c(pop_params$lambda_std[1:3,1], 
                                                   pop_params$lambda_std[4:6,2]))
    # Common-factor correlation matrix
    Phi_gauss <- nonagg_summaries(results = results, 
                                  parameter_name = "Phi_gauss", 
                                  parameter = pop_params$Phi_cor[1,2], 
                                  Phi.lv = TRUE)
    Phi_skew <- nonagg_summaries(results = results, 
                                 parameter_name = "Phi_skew",
                                 parameter = pop_params$Phi_cor[1,2], 
                                 Phi.lv = TRUE)
  }
  
  # Extra parameters
  if(Design[i, "Model"] == "exgaussian"){
    # Save ex-gaussian parameters
    mu_tau <- nonagg_summaries(results = results, parameter_name = "mu_tau", parameter = pop_params$mu_tau)
    sd_tau <- nonagg_summaries(results = results, parameter_name = "sd_tau", parameter = pop_params$sd_tau)
  } else {
    # Save shifted-lognormal parameters
    mu_delta <- nonagg_summaries(results = results, parameter_name = "mu_delta", parameter = pop_params$mu_delta)
    sd_delta <- nonagg_summaries(results = results, parameter_name = "sd_delta", parameter = pop_params$sd_delta)
  }
  
  # Biases data frame per condition
  biases_i <- cbind(
    sim_ID = i,
    Design[i,], 
    rbind(
      # Population-level effects
      mu_alpha, mu_theta, mu_sigma, sd_alpha, sd_theta, sd_sigma,
      # Correlations
      Rho_gauss, Rho_skew, Rho_mm, Rho_sp,
      # Reliability
      reliab_gauss, reliab_skew,
      # Model-implied latent structure
      lambda_gauss, lambda_skew, crossload_gauss, crossload_skew, Phi_gauss, Phi_skew,
      # exgaussian and shifted-lognormal parameters
      mu_delta, sd_delta, mu_tau, sd_tau
    )
  )
  
  # Save final data frame
  temp_data_lists[[i]] <- biases_i
}

# Final long data frame
full_simres <- do.call(rbind, temp_data_lists)

# Aggregate results across the same parameter-type in each replica
full_simres_adj <- full_simres |> 
  group_by(sim_ID, Model, M, I, reliab, std_lambda_j, rep, parameter_set) |> 
  summarise(
    true_value    = mean(true_value,    na.rm = TRUE),
    est           = mean(est,           na.rm = TRUE),
    se            = mean(se,            na.rm = TRUE),
    bias          = mean(bias,          na.rm = TRUE),
    abias         = mean(abias,         na.rm = TRUE),
    rbias         = mean(rbias,         na.rm = TRUE),
    arbias_pct    = mean(arbias_pct,    na.rm = TRUE),
    squared_error = mean(squared_error, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    # Create a method variable
    method = case_when(
      parameter_set == "Rho_gauss"    ~ "normal",
      parameter_set == "Rho_mm"       ~ "naive",
      parameter_set == "Rho_skew"     ~ "skew",
      parameter_set == "Rho_sp"       ~ "spearman",
      parameter_set == "reliab_gauss" ~ "normal",
      parameter_set == "reliab_skew"  ~ "skew",
      grepl("_gauss$", parameter_set) ~ "normal",
      grepl("_skew$",  parameter_set) ~ "skew",
      TRUE ~ "skew"  
    ),
    # Remove sufix index
    parameter_set = sub("_(gauss|skew|mm|sp)$", "", parameter_set)
  ) %>%
  relocate(method, .after = parameter_set)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Plots – correlations by method and model
# ─────────────────────────────────────────────────────────────────────────────

# Colors for each model
pal_model <- c("Ex-Gaussian" = "#E08F2D",
               "Shifted-lognormal" = "#1E88E5")

rho_both_models_plot <-
  full_simres |>
  mutate(
    method = case_when(
      parameter_set == "Rho_gauss"    ~ "normal",
      parameter_set == "Rho_mm"       ~ "naive",
      parameter_set == "Rho_skew"     ~ "skew",
      parameter_set == "Rho_sp"       ~ "spearman",
      parameter_set == "reliab_gauss" ~ "normal",
      parameter_set == "reliab_skew"  ~ "skew",
      grepl("_gauss$", parameter_set) ~ "normal",
      grepl("_skew$",  parameter_set) ~ "skew",
      TRUE ~ "skew"
    ),
    parameter_set = sub("_(gauss|skew|mm|sp)$", "", parameter_set),
    Model = dplyr::recode(
      Model,
      "exgaussian"       = "Ex-Gaussian",
      "shifted.lognormal" = "Shifted-lognormal"
    ),
    method = factor(
      method,
      levels = c("naive", "spearman", "normal", "skew"),
      labels = c("Aggregated", "Spearman", "Gaussian HFM", "True GenHFM")
    )
  ) |>
  filter(
    parameter_set == "Rho",
    Model %in% c("Ex-Gaussian", "Shifted-lognormal")
  ) |>
  ggplot(aes(
    x     = est,
    y     = true_value,
    color = Model,
    shape = Model,
    linetype = Model
  )) +
  geom_abline(
    slope     = 1,
    intercept = 0,
    linetype  = "22",
    linewidth = 0.8,
    color     = "grey25"
  ) +
  stat_pointinterval(
    aes(group = interaction(true_value, Model)),
    orientation         = "horizontal",
    point_interval      = median_qi,
    .width              = c(0.80, 0.95),
    interval_size_range = c(0.5, 1),
    point_size          = 2.7,
    show_point          = TRUE,
    interval_alpha      = 0.45
  ) +
  facet_grid(
    method ~ reliab,
    labeller = labeller(
      method = label_value,
      reliab = function(x) paste0("Reliability = ", x)
    )
  ) +
  scale_color_manual(values = pal_model, name = "Model") +
  scale_shape_manual(values = c("Ex-Gaussian" = 16,
                                "Shifted-lognormal" = 17),
                     name = "Model") +
  scale_linetype_manual(values = c("Ex-Gaussian" = "solid",
                                   "Shifted-lognormal" = "dashed"),
                        name = "Model") +
  scale_x_continuous(
    breaks = seq(-0.3, 0.7, by = 0.15),
    labels = scales::number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  scale_y_continuous(
    breaks = seq(0.0, 0.6, by = 0.1),
    labels = scales::number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  coord_equal(
    clip = "on",
    xlim = c(-0.3, 0.7),
    ylim = c(0, 0.6)
  ) +
  labs(
    x = "Estimated correlation",
    y = "True correlation"
  ) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    plot.title.position = "panel",
    axis.title.x        = element_text(size = 28),
    axis.title.y        = element_text(size = 28),
    strip.text.y.right  = element_text(
      angle  = 270,
      vjust  = 0.5,
      hjust  = 0.5,
      margin = margin(l = 6, r = 6)
    ),
    strip.text          = element_text(face = "bold"),
    strip.background.y  = element_rect(fill = "grey92", color = NA),
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_line(linewidth = 0.4, color = "grey80"),
    panel.grid.major.y  = element_line(linewidth = 0.4, color = "grey80"),
    panel.spacing.x     = unit(10, "pt"),
    panel.spacing.y     = unit(10, "pt"),
    plot.title          = element_text(face = "bold", size = 28),
    legend.position     = "top",
    legend.text         = element_text(size = 22),
    legend.title        = element_text(size = 24, face = "bold")
  )

# Save the final figure
png(filename = "Figures/correlations_both_models.png",
    width = 4000, height = 4000, res = 300)
rho_both_models_plot
dev.off()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Plots – factor loadings by method and model
# ─────────────────────────────────────────────────────────────────────────────

lambda_plot_models <- full_simres |>
  filter(parameter_set %in% c("lambda_gauss", "lambda_skew")) |>
  mutate(
    method = case_when(
      parameter_set == "lambda_gauss" ~ "Gaussian HFM",
      parameter_set == "lambda_skew"  ~ "True GenHFM"
    ),
    method = factor(
      method,
      levels = c("Gaussian HFM", "True GenHFM")
    ),
    Model = factor(
      Model,
      levels = c("exgaussian", "shifted.lognormal"),
      labels = c("Ex-Gaussian", "Shifted-lognormal")
    )
  ) |>
  ggplot(aes(
    x        = est,
    y        = true_value,
    color    = Model,
    shape    = Model,
    linetype = Model
  )) +
  geom_abline(
    slope     = 1,
    intercept = 0,
    linetype  = "22",
    linewidth = 0.8,
    color     = "grey25"
  ) +  
  stat_pointinterval(
    aes(group = interaction(true_value, Model)),
    position            = "identity",
    orientation         = "horizontal",
    point_interval      = median_qi,
    .width              = c(0.80, 0.95),
    interval_size_range = c(0.5, 1),
    point_size          = 3.2,
    show_point          = TRUE,
    interval_alpha      = 0.45
  ) +
  facet_grid(
    method ~ reliab,
    labeller = labeller(
      method = label_value,
      reliab = function(x) paste0("Reliability = ", x)
    )
  ) +
  scale_color_manual(
    values = pal_model,
    name   = "Model"
  ) +
  scale_shape_manual(
    values = c("Ex-Gaussian" = 16,
               "Shifted-lognormal" = 17),
    name   = "Model"
  ) +
  scale_linetype_manual(
    values = c("Ex-Gaussian" = "solid",
               "Shifted-lognormal" = "dashed"),
    name   = "Model"
  ) +
  scale_x_continuous(
    breaks = seq(-0.3, 1.2, by = 0.2),
    labels = scales::number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  scale_y_continuous(
    breaks = seq(-0.3, 1.2, by = 0.2),
    labels = scales::number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  coord_equal(clip = "on", xlim = c(-.2, 0.9), ylim = c(-.2, 0.9)) +
  labs(
    x = "Estimated factor loading",
    y = "True factor loading"
  ) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    plot.title.position = "panel",
    axis.title.x        = element_text(size = 28),
    axis.title.y        = element_text(size = 28),
    strip.text.y.right = element_text(
      angle  = 270,
      vjust  = 0.5,
      hjust  = 0.5,
      margin = margin(l = 6, r = 6)
    ),
    strip.text         = element_text(face = "bold"),
    strip.background.y = element_rect(fill = "grey92", color = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.4, color = "grey80"),
    panel.grid.major.y = element_line(linewidth = 0.4, color = "grey80"),
    panel.spacing.x    = unit(10, "pt"),
    panel.spacing.y    = unit(10, "pt"),
    legend.position    = "top",
    legend.text        = element_text(size = 22),
    legend.title       = element_text(size = 24, face = "bold")
  )

# Save final plot
png(filename = "Figures/factor_loadings_simresults.png", 
    width = 4000, height = 4000, res = 300)
lambda_plot_models
dev.off()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Plots – ex-gaussian hierarchical component parameters
# ─────────────────────────────────────────────────────────────────────────────

exg_pars <- c("mu_alpha", "mu_theta", "mu_sigma", "mu_tau", 
              "sd_alpha", "sd_theta", "sd_sigma", "sd_tau")

hier_exg_dat <- full_simres_adj |> 
  filter(Model == "exgaussian") |> 
  filter(parameter_set %in% exg_pars) |> 
  mutate(parameter_set = if_else(
    parameter_set == "sd_theta",
    paste0("sd_theta_r", sprintf("%.0f", reliab * 100)),
    parameter_set
  ),
  par_lab = dplyr::recode(
    parameter_set,
    "mu_alpha"     = "mu[alpha]",
    "mu_sigma"     = "mu[sigma]",
    "mu_tau"       = "mu[tau]",
    "mu_theta"     = "mu[beta]",
    "sd_alpha"     = "sigma[alpha]",
    "sd_sigma"     = "sigma[sigma]",
    "sd_tau"       = "sigma[tau]",
    "sd_theta_r30" = "sigma[beta]*' ('*rho[xx]*' = 0.30)'",
    "sd_theta_r50" = "sigma[beta]*' ('*rho[xx]*' = 0.50)'",
    "sd_theta_r70" = "sigma[beta]*' ('*rho[xx]*' = 0.70)'"
  )
  )

true_vals_exg <- hier_exg_dat |> 
  distinct(par_lab, true_value)

exg_hier_plot <- ggplot(hier_exg_dat, aes(x = est)) +
  geom_histogram(bins = 30, color = "white", fill = "grey70") +
  geom_vline(
    data  = true_vals_exg,
    aes(xintercept = true_value),
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(
    ~ par_lab,
    scales  = "free",
    ncol    = 5,
    labeller = label_parsed
  ) +
  bayesplot::theme_default(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Estimated values",
    y = "Frequency",
  ) + scale_x_continuous(breaks = scales::pretty_breaks(4))


# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Plots – shifted-lognormal hierarchical component parameters
# ─────────────────────────────────────────────────────────────────────────────

shl_pars <- c("mu_alpha", "mu_theta", "mu_sigma", "mu_delta", 
              "sd_alpha", "sd_theta", "sd_sigma", "sd_delta")

hier_shl_dat <- full_simres_adj |> 
  filter(Model == "shifted.lognormal") |> 
  filter(parameter_set %in% shl_pars) |> 
  mutate(parameter_set = if_else(
    parameter_set == "sd_theta",
    paste0("sd_theta_r", sprintf("%.0f", reliab * 100)),
    parameter_set
  ),
  # Etiquetas bonitas para los facets (plotmath en formato texto)
  par_lab = dplyr::recode(
    parameter_set,
    "mu_alpha"     = "mu[alpha]",
    "mu_sigma"     = "mu[sigma]",
    "mu_delta"     = "mu[delta]",
    "mu_theta"     = "mu[beta]",
    "sd_alpha"     = "sigma[alpha]",
    "sd_sigma"     = "sigma[sigma]",
    "sd_delta"     = "sigma[delta]",
    "sd_theta_r30" = "sigma[beta]*' ('*rho[xx]*' = 0.30)'",
    "sd_theta_r50" = "sigma[beta]*' ('*rho[xx]*' = 0.50)'",
    "sd_theta_r70" = "sigma[beta]*' ('*rho[xx]*' = 0.70)'"
  )
  )

true_vals <- hier_shl_dat |> 
  distinct(par_lab, true_value)

shl_hier_plot <- ggplot(hier_shl_dat, aes(x = est)) +
  geom_histogram(bins = 30, color = "white", fill = "grey70") +
  geom_vline(
    data  = true_vals,
    aes(xintercept = true_value),
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(
    ~ par_lab,
    scales  = "free",
    ncol    = 5,
    labeller = label_parsed
  ) +
  bayesplot::theme_default(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Estimated values",
    y = "Frequency",
  ) + scale_x_continuous(breaks = scales::pretty_breaks(4))

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Plots – hierarchical component parameters plot
# ─────────────────────────────────────────────────────────────────────────────

# Bind both plots
combined_plot <- exg_hier_plot / shl_hier_plot +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

# Save final plot
png(filename = "Figures/hierarchical_pars_simresults.png", width = 5000, height = 3500, res = 300)
combined_plot
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
