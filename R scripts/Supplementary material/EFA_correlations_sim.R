# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : EFA_correlations_sim.R                                     ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 06-08-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Small simulation so show that Factor models compute more precise
# correlations between measures
# 
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(SimDesign)
library(mvnfast)
library(dplyr)
library(tidyr)
library(ggplot2)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Simulation study
# ─────────────────────────────────────────────────────────────────────────────

# Simulation design
Design <- createDesign(N = c(50, 100, 500),
                       rho = seq(0, .95, .05))

# Number of items
items_5  <- list(items = 5)
items_10 <- list(items = 10)
items_20 <- list(items = 20)

# Data-generation function
Generate <- function(condition, fixed_objects) {
  R <- diag(1 - condition$rho, fixed_objects$items) + condition$rho
  dat <- mvnfast::rmvn(n = condition$N, mu = rep(0, fixed_objects$items), sigma = R)
  dat
}

# Dana analyse function
Analyse <- function(condition, dat, fixed_objects) {
  obs_cor <- cor(dat)[lower.tri(diag(fixed_objects$items))]
  efa_cor <- tcrossprod(factanal(x = dat, factors = 1)$loadings)[lower.tri(diag(fixed_objects$items))]
  res <- nc(observed = obs_cor, model = efa_cor)
  res
}

# Data summarise function
Summarise <- function(condition, results, fixed_objects) {
  obs_avg <- mean(colMeans(results[,grep("observed", colnames(results))]))
  mod_avg <- mean(colMeans(results[,grep("model", colnames(results))]))
  obs_SE  <- mean(colSDs(results[,grep("observed", colnames(results))]))
  mod_SE  <- mean(colSDs(results[,grep("model", colnames(results))]))
  sumres <- c(obs_avg = obs_avg, mod_avg = mod_avg, 
              obs_SE = obs_SE, mod_SE = mod_SE)
  sumres
}

# Run simulation study with 5 items
res_5 <- runSimulation(design=Design, replications=1e4L, generate=Generate, 
                     fixed_objects = items_5,
                     analyse=Analyse, summarise=Summarise, parallel = TRUE,
                     packages = c("mvnfast"),
                     ncores = 10)

# Run simulation study with 10 items
res_10 <- runSimulation(design=Design, replications=1e4L, generate=Generate, 
                        fixed_objects = items_10,
                        analyse=Analyse, summarise=Summarise, parallel = TRUE,
                        packages = c("mvnfast"),
                        ncores = 10)

# Run simulation study with 20 items
res_20 <- runSimulation(design=Design, replications=1e4L, generate=Generate, 
                       fixed_objects = items_20,
                       analyse=Analyse, summarise=Summarise, parallel = TRUE,
                       packages = c("mvnfast"),
                       ncores = 10)

# Combine suumarised results
equiv_pearson_CF <- cbind(items = c(rep(5, 60), rep(10, 60), rep(20, 60)), 
                                   rbind(res_5, res_10, res_20))

# Save summarized results
save(equiv_pearson_CF, file = "Results/Rdata/Supplementary material/equiv_pearson_CF.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Plot of results
# ─────────────────────────────────────────────────────────────────────────────

# Ratio between empirical standard deviations
df_ratio <- equiv_pearson_CF %>% 
  mutate(
    se_ratio = obs_SE / mod_SE,
    items    = factor(items, levels = c(5, 10, 20)),
    N        = factor(N)
  )

# Colors for each number of items
item_cols <- c(`5`  = "#CB8D78", `10` = "#9C584E", `20` = "#6E242B")

# Shapes for each item
item_shapes <- c(`5` = 16, `10` = 17, `20` = 15)

# Plot of ratios between empirical standard deviations
p_ratio <- ggplot(
  df_ratio,
  aes(x = rho, y = se_ratio,
      colour = items, shape = items, group = items)
) +
  geom_line(linewidth = 1.6) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = 2,
             linewidth  = 0.6, colour = "grey40") +
  facet_wrap(~ N, nrow = 1, labeller = labeller(N = \(x) paste0("N = ", x))) +
  scale_colour_manual(values = item_cols, name = "Number of items") +
  scale_shape_manual (values = item_shapes, name = "Number of items") +
  scale_x_continuous(
    breaks = seq(0, 0.9, 0.1),
    expand = expansion(mult = c(.02, .02))
  ) +
  scale_y_continuous(
    limits = c(0.75, NA),
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(.02, .05))
  ) +
  labs(
    x = bquote("True correlation" ~ (rho)),
    y = expression(SE[Pearson] / SE[Common-factor])
  ) +
  theme_bw(base_size = 24) %+replace%
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(1.2, "lines"),
    
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.text.x       = element_text(face = "bold"),
    
    axis.title.y       = element_text(angle = 90,
                                      vjust = .5, hjust = .5,
                                      margin = margin(r = 8)),
    axis.title.x       = element_text(margin = margin(t = 8)),
    
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    legend.key.width   = unit(.9, "cm")
  )

png(filename = "Figures/supplemental_SEratio.png", width = 5000, height = 3000, res = 300)
print(p_ratio)
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
