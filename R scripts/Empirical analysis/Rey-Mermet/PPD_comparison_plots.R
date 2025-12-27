# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : PPD_comparison.R                                           ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 16-11-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Plot posterior predictive distribution of observed RTs using Ex-Gaussian HFM
# and DDMs fitted in Rey-Mermet, Singmann & Oberauer, 2025
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(bayesplot)
library(ggplot2)
library(ggpubr)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User defined functions
# ─────────────────────────────────────────────────────────────────────────────

# Load empirical R functions
source("R scripts/R functions/empirical_functions.R")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Posterior Predictive Distribution: Experiment 1
# ─────────────────────────────────────────────────────────────────────────────

# Load Rey-Mermet, Singmann and Oberauer DDM fitted model (Exp. 1)
load("R scripts/Empirical analysis/Rey-Mermet/DDM_fitted_models/fit_wh_e1_1.rda")

# Just the first 100 posterior draws
pred_fit_wh_e1_1 <- pred_fit_wh_e1_1[1:100,]

# DDM posterior predictive density
color_scheme_set("pink")
ppd_wiener_exp1 <- ppc_dens_overlay(
  y = fit_wh_e1_1$data$rt,
  yrep = abs(pred_fit_wh_e1_1) # Negative RTs for lower bound predictions
)

# Remove posterior predictive distribution to save space
rm(pred_fit_wh_e1_1); gc()

# Posterior Predictive distribution: ex-gaussian model
exgaussian_fit_E1 <- readRDS("Results/Stan/Rey-Mermet/exgaussian_HFM_fit_E1.rds")

# Prepare Stan data (required for PPC)
load("R scripts/Empirical analysis/Rey-Mermet/datasets/e1fit.rdata")
sdata <- Stan.list.data(
  data = e1fit, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "congruency", 
  congruent_val = "congruent", 
  incongruent_val = "incongruent", 
  rt_var = "rt"
)

# Draw values from PPD
ppc_exg <- PPC_simdata(
  fit = exgaussian_fit_E1, 
  sdata = sdata, 
  n_draws = 100, 
  model = "exgaussian")

# Posterior Predictive distribution
color_scheme_set("brightblue")
ppd_exg_exp1 <- ppc_dens_overlay(
  y = sdata_to_longdf(sdata)$rt,
  yrep = t(ppc_exg[,grep("draw", colnames(ppc_exg))])
)

# Remove objects to save space
rm(fit_wh_e1_1)
rm(ppc_exg); gc()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Posterior Predictive Distribution: Experiment 2
# ─────────────────────────────────────────────────────────────────────────────

# Load Rey-Mermet, Singmann and Oberauer DDM fitted model (Exp. 2)
load("R scripts/Empirical analysis/Rey-Mermet/DDM_fitted_models/fit_wh_e2_1.rda")

# Just the first 100 posterior draws
pred_fit_wh_e2_1 <- pred_fit_wh_e2_1[1:100,]

# DDM posterior predictive density
color_scheme_set("pink")
ppd_wiener_exp2 <- ppc_dens_overlay(
  y = fit_wh_e2_1$data$rt,
  yrep = abs(pred_fit_wh_e2_1) # Negative RTs for lower bound predictions
)

# Remove posterior predictive distribution to save space
rm(pred_fit_wh_e2_1); gc()

# Posterior Predictive distribution: ex-gaussian model
exgaussian_fit_E2 <- readRDS("Results/Stan/Rey-Mermet/exgaussian_HFM_fit_E2.rds")

# Prepare Stan data (required for PPC)
load("R scripts/Empirical analysis/Rey-Mermet/datasets/e2.rdata")
sdata <- Stan.list.data(
  data = e2, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "congruency", 
  congruent_val = "congruent", 
  incongruent_val = "incongruent", 
  rt_var = "rt"
)

# Draw values from PPD
ppc_exg <- PPC_simdata(
  fit = exgaussian_fit_E2, 
  sdata = sdata, 
  n_draws = 100, 
  model = "exgaussian")

# Posterior Predictive distribution
color_scheme_set("brightblue")
ppd_exg_exp2 <- ppc_dens_overlay(
  y = sdata_to_longdf(sdata)$rt,
  yrep = t(ppc_exg[,grep("draw", colnames(ppc_exg))])
)

# Remove objects to save space
rm(fit_wh_e2_1)
rm(ppc_exg); gc()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Posterior Predictive Distribution: Experiment 3
# ─────────────────────────────────────────────────────────────────────────────

# Load Rey-Mermet, Singmann and Oberauer DDM fitted model (Exp. 3)
load("R scripts/Empirical analysis/Rey-Mermet/DDM_fitted_models/fit_wh_e3_1.rda")

# Just the first 100 posterior draws
pred_fit_wh_e3_1 <- pred_fit_wh_e3_1[1:100,]

# DDM posterior predictive density
color_scheme_set("pink")
ppd_wiener_exp3 <- ppc_dens_overlay(
  y = fit_wh_e3_1$data$rt,
  yrep = abs(pred_fit_wh_e3_1) # Negative RTs for lower bound predictions
)

# Remove posterior predictive distribution to save space
rm(pred_fit_wh_e3_1); gc()

# Posterior Predictive distribution: ex-gaussian model
exgaussian_fit_E3 <- readRDS("Results/Stan/Rey-Mermet/exgaussian_HFM_fit_E3.rds")

# Prepare Stan data (required for PPC)
load("R scripts/Empirical analysis/Rey-Mermet/datasets/e3.rdata")
sdata <- Stan.list.data(
  data = e3, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "congruency", 
  congruent_val = "congruent", 
  incongruent_val = "incongruent", 
  rt_var = "rt"
)

# Draw values from PPD
ppc_exg <- PPC_simdata(
  fit = exgaussian_fit_E3, 
  sdata = sdata, 
  n_draws = 100, 
  model = "exgaussian")

# Posterior Predictive distribution
color_scheme_set("brightblue")
ppd_exg_exp3 <- ppc_dens_overlay(
  y = sdata_to_longdf(sdata)$rt,
  yrep = t(ppc_exg[,grep("draw", colnames(ppc_exg))])
)

# Remove objects to save space
rm(fit_wh_e3_1)
rm(ppc_exg); gc()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: plot with all posterior distributions
# ─────────────────────────────────────────────────────────────────────────────

# Set a global theme
global_theme <- bayesplot::theme_default(base_size = 36) +
  theme(
    plot.title = element_text(hjust = 0, size = 38),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 36)
  )

# Plot limits
global_scales <- coord_cartesian(xlim = c(0, 1.8), ylim = c(0, 2.5))

# DDMs plot list
list_ddm <- list(
  ppd_wiener_exp1 + ggtitle("DDM: Experiment 1") + global_theme + global_scales,
  ppd_wiener_exp2 + ggtitle("DDM: Experiment 2") + global_theme + global_scales,
  ppd_wiener_exp3 + ggtitle("DDM: Experiment 3") + global_theme + global_scales
)

# Ex-Gaussian plot list
list_exg <- list(
  ppd_exg_exp1 + ggtitle("Ex-Gaussian: Experiment 1") + global_theme + global_scales,
  ppd_exg_exp2 + ggtitle("Ex-Gaussian: Experiment 2") + global_theme + global_scales,
  ppd_exg_exp3 + ggtitle("Ex-Gaussian: Experiment 3") + global_theme + global_scales
)

# First row of plots: DDMs
row_ddm <- ggarrange(plotlist = list_ddm, 
                     ncol = 3, nrow = 1, 
                     common.legend = TRUE, legend = "right")

# Second row of plots: Ex-Gaussian HFM
row_exg <- ggarrange(plotlist = list_exg, 
                     ncol = 3, nrow = 1, 
                     common.legend = TRUE, legend = "right")

# Joint both rows of plots
final_figure <- ggarrange(row_ddm, row_exg, ncol = 1, nrow = 2)

# Add X and Y labels and save the plot
png(filename = "Figures/PPD_comparison.png", 
    width = 7000, height = 4500, res = 300)
annotate_figure(final_figure,
                left = text_grob("Density", rot = 90, size = 40, family = "serif"),
                bottom = text_grob("Response Time", size = 40, family = "serif"))
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
