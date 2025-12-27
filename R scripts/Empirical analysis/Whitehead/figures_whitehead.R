# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : figures_whitehead.R                                        ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 27-08-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Plots of the ex-gaussian and gaussian hierarchical factor models 
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(posterior)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions and fitted models
# ─────────────────────────────────────────────────────────────────────────────

# Empirical functions to generate posterior predictive model checks
source("R scripts/R functions/empirical_functions.R") 

# Function to extract and prepare factor loadings posterior draws
extract_lambdas <- function(fit, experiment, model, task_labels = c(
  "1" = "Simon",
  "2" = "Flanker",
  "3" = "Stroop"
)) {
  posterior::as_draws_df(fit$draws("Lambda_std")) |>
    tidyr::pivot_longer(
      cols = starts_with("Lambda_std"),
      names_to = "param",
      values_to = "lambda"
    ) |>
  dplyr::filter(stringr::str_detect(param, "Lambda_std\\[[0-9]+,1\\]")) |>
    dplyr::mutate(
      task_index = str_match(param, "Lambda_std\\[([0-9]+),1\\]")[, 2],
      task = recode(task_index, !!!task_labels),
      experiment = experiment,
      model = model
    ) |>
    dplyr::select(experiment, model, task, .draw, lambda)
}

# Load fitted models
gaussian_GHFM_fit_E1   <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E1.rds")
shlognorm_GHFM_fit_E1  <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E1.rds")
exgaussian_GHFM_fit_E1 <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E1.rds")
gaussian_GHFM_fit_E2   <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E2.rds")
shlognorm_GHFM_fit_E2  <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E2.rds")
exgaussian_GHFM_fit_E2 <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E2.rds")
gaussian_GHFM_fit_E3   <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E3.rds")
shlognorm_GHFM_fit_E3  <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E3.rds")
exgaussian_GHFM_fit_E3 <- readRDS(file = "Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Factor loadings posterior distribution
# ─────────────────────────────────────────────────────────────────────────────

# Set task labels
task_labels <- c(
  "1" = "Simon",
  "2" = "Flanker",
  "3" = "Stroop"
)

# Prepare the lambda data frame
lambda_df <- bind_rows(
  extract_lambdas(gaussian_GHFM_fit_E1,  "E1", "Gaussian"),
  extract_lambdas(exgaussian_GHFM_fit_E1,"E1", "Ex-Gaussian"),
  extract_lambdas(shlognorm_GHFM_fit_E1, "E1", "Sh-lognormal"),
  extract_lambdas(gaussian_GHFM_fit_E2,  "E2", "Gaussian"),
  extract_lambdas(exgaussian_GHFM_fit_E2,"E2", "Ex-Gaussian"),
  extract_lambdas(shlognorm_GHFM_fit_E2, "E2", "Sh-lognormal"),
  extract_lambdas(gaussian_GHFM_fit_E3,  "E3", "Gaussian"),
  extract_lambdas(exgaussian_GHFM_fit_E3,"E3", "Ex-Gaussian"),
  extract_lambdas(shlognorm_GHFM_fit_E3, "E3", "Sh-lognormal")
  ) |> 
  mutate(
    experiment = factor(experiment, levels = c("E1","E2","E3")) ,
    task       = factor(task,       levels = c("Simon","Flanker","Stroop"))
  )

# Color palette
colores_modelo <- c(
  "Gaussian"          = "grey40",
  "Ex-Gaussian"       = "#E08F2D",
  "Shifted-lognormal" = "#1E88E5"
)

# Experiment labels
lab_exp <- c(
  "E1" = "Experiment 1",
  "E2" = "Experiment 2",
  "E3" = "Experiment 3"
)

# Posterior distribution
plot_loads <- lambda_df |>
  mutate(
    model = dplyr::recode(
      model,
      "Gaussian"     = "Gaussian",
      "Ex-Gaussian"  = "Ex-Gaussian",
      "Sh-lognormal" = "Shifted-lognormal"
    ),
    model = factor(
      model,
      levels = c("Gaussian", "Ex-Gaussian", "Shifted-lognormal")
    )
  ) |>
  ggplot(aes(x = lambda, colour = model, fill = model)) +
  geom_density(alpha = 0.20, linewidth = 0.9, adjust = 1.1) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.6,
    color      = "firebrick4"
  ) +
  facet_grid(
    experiment ~ task,
    labeller = labeller(experiment = lab_exp)
  ) +
  scale_colour_manual(values = colores_modelo, name = "Model") +
  scale_fill_manual(values   = colores_modelo, name = "Model") +
  labs(
    x = "Standardized factor loadings",
    y = "Posterior density"
  ) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    legend.position    = "top",
    legend.title       = element_text(face = "bold", size = 24),
    legend.text        = element_text(size = 22),
    legend.key.width   = unit(1.4, "lines"),
    axis.title.x       = element_text(margin = margin(t = 8), size = 28),
    axis.title.y       = element_text(margin = margin(r = 8), size = 28),
    strip.text.x       = element_text(face = "bold", size = 18,
                                      margin = margin(t = 7, b = 7, l = 5, r = 5)),
    strip.text.y       = element_text(face = "bold", size = 18,
                                      margin = margin(t = 10, b = 10, l = 5, r = 5)),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    panel.border       = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(1.2, "lines"),
    panel.spacing.y    = unit(1.6, "lines")
  )

# Save this plot
png(filename = "Figures/whitehead_posterior_lambdas.png", width = 4000, height = 3000, res = 300)
plot_loads
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
