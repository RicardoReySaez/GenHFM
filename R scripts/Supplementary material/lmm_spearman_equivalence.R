# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : lmm_spearman_equivalence.R                                 ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 07-08-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Script to assess to what extent the hierarchical solution is better than
# the spearman correction for attenuation without informative priors
# (i.e., just using a frequentist hierarchical model)
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(SimDesign)
library(mvnfast)
library(lme4)
library(R2jags)
library(tidyverse)
library(ggplot2)
library(ggdist)
library(patchwork)
library(grid) 

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Rouder et al. (2023) simulation function
# retrieved from: https://github.com/PerceptionAndCognitionLab/ctx-inhibition/blob/public/papers/rev3/lib.R
sim2=function(vals,trueCor){
  I=vals$I #ppl
  L=vals$L # reps
  J=2
  alpha.var=diag(rep(vals$alpha.var,J))
  t.alpha=rmvnorm(I,rep(.8,J),alpha.var)
  t.theta=mvrnorm(I,
                  rep(.06,J),
                  matrix(ncol=2,
                         rep(vals$theta.var,4)*c(1,trueCor,trueCor,1)))
  t.s2=vals$s2
  K=2
  N=I*J*K*L
  sub=rep(1:I,each=J*K*L)
  task=rep(rep(1:J,each=K*L),I)
  cond=rep(rep(1:K,each=L),I*J)
  subtask=cbind(sub,task)
  t.cell=t.alpha[subtask]+(cond-1)*t.theta[subtask]
  rt=rnorm(N,t.cell,sqrt(t.s2))
  dat=data.frame(sub,task,cond,rt)
  out=list(dat=dat,t.theta=t.theta)
  return(out)
}

# Mehrvarz & Rouder (2025) bayesian hierarchical model in JAGS
# retrieved from: https://github.com/specl/ctx-pca-localize/blob/main/_logistics/mods.R
modJ.wishart = "
model{
  # Priors
  for (i in 1:I){
    alpha[i, 1:J] ~ dmnorm(nu, pPsi2)
    theta[i, 1:J] ~ dmnorm(mu, pSig2)
  }
  
  for (j in 1:J){
    nu[j] ~ dnorm(.8, pow(.3, -2))
    mu[j] ~ dnorm(.05, pow(.1, -2))
    pTau2[j] ~ dgamma(.5,.5)
  }
  
  pPsi2 ~ dwish(diagJ, J+1)
  pSig2 ~ dwish(diagJ*(.025^2), J+1)
  
  
  # Likelihood
  for (n in 1:N){
    center[n] = alpha[sub[n], task[n]] + (cond[n] - 1.5) * theta[sub[n], task[n]]
    y[n] ~ dnorm(center[n], pTau2[task[n]])
  }
}
"

# Rouder et al. (2023) spearman correction for attenuation function
# retrieved from: https://github.com/PerceptionAndCognitionLab/ctx-inhibition/blob/public/papers/rev3/lib.R
spearman=function(dat){
  mrt=tapply(dat$rt,list(dat$sub,dat$task,dat$cond),mean)
  sample.effect=mrt[,,2]-mrt[,,1]
  se <- function(x) sd(x)/sqrt(length(x))
  sert <- tapply(dat$rt, list(dat$sub,dat$task,dat$cond), se)
  # (squared) standard error of difference in means:
  se2.diff <- colMeans(sert[,,2]^2 + sert[,,1]^2)
  cov.diff <- cov(sample.effect) - diag(se2.diff)
  return(cov2cor(cov.diff))
}

# Rouder et al. (2023) simulation values
# Retrieved from: https://github.com/PerceptionAndCognitionLab/ctx-inhibition/blob/public/papers/rev3/p.Rmd
typical=list(
  I=200, #ppl
  L=100, # reps
  alpha.var=.2^2,
  theta.var=.025^2,
  s2=.200^2)

# Reliability assumed by Rouder et al. (2023) for all tasks
(.025^2/.2^2) / ((.025^2/.2^2) + 2 / 100) 

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: simulation study
# ─────────────────────────────────────────────────────────────────────────────

# Simulation conditions (true correlation value)
Design <- createDesign(trueCor = seq(.1, .9, .1))

# Fixed objects inside the simulation study
fixed_objects <- list(typical = typical, 
                      spearman = spearman, 
                      IW_mod = modJ.wishart)

# Data-generation function
Generate <- function(condition, fixed_objects) {
  dat <- sim2(vals = fixed_objects$typical, trueCor = condition$trueCor)$dat
  dat
}

# Analyse functions
Analyse <- function(condition, dat, fixed_objects) {
  # First analysis: spearman correction formula
  sp_cor <- fixed_objects$spearman(dat)[1,2]
  
  # Second analysis: Frequentist hierarchical model
  dat$task <- as.factor(dat$task)
  lme4_mod <- lme4::lmer(
    # Group-level alpha and theta values
    rt ~ 0 + task + cond:task + 
      # Subject-specific alpha values
      (0 + task | sub) + 
      # Subject-specific theta values  
      (0 + cond:task | sub), 
    data = dat,control = lme4::lmerControl(optimizer = "bobyqa"))
  # Check convergence
  if(!is.null(lme4_mod@optinfo$conv$lme4$code) | lme4::isSingular(lme4_mod)){
    stop("Convergence problem")
  }
  
  # Model-implied true correlation
  lme4_cor <- attr(lme4::VarCorr(lme4_mod)$sub.1, "correlation")[1,2]
  
  # Third analysis: Bayesian hierarchical model
  BHM_mod <- jags(
    data = list(
      "y" = dat$rt,
      "task" = as.numeric(dat$task),
      "sub" = dat$sub,
      "cond" = dat$cond,
      "I" = 200,
      "J" = 2,
      "diagJ" = diag(2),
      "N" = nrow(dat)
    ), parameters.to.save = c("pSig2"),
    n.chains = 2,
    n.iter = 2000, 
    n.burnin = 1000, 
    model.file = textConnection(fixed_objects$IW_mod))
  
  # Model-implied true correlation
  BHM_cor <- mean(apply(BHM_mod$BUGSoutput$sims.matrix[,2:5], 1, function(x) {
    cov2cor(solve(matrix(x, ncol=2)))[1,2]
  }))
  
  # Results: correlation estimates with each method
  res <- nc(sp_cor = sp_cor, lme4_cor = lme4_cor, BHM_cor = BHM_cor)
  res
}

# Summarise function
Summarise <- function(condition, results, fixed_objects) {
  ret <- c(avg_res = colMeans(results), SEs_res = colSDs(results))
  ret
}

# Run parallelized simulation study
sim_res <- runSimulation(
  design          =  Design, 
  generate        =  Generate, 
  analyse         =  Analyse, 
  summarise       =  Summarise,
  fixed_objects   =  fixed_objects,
  packages        =  c("R2jags", "lme4", "MASS"),
  replications    =  1000, 
  parallel        =  TRUE,
  ncores          =  10, 
  filename        =  "Results/Rdata/Supplementary material/lme4_spearman_equiv.rds")

sim_res <- readRDS("Results/Rdata/Supplementary material/lme4_spearman_equiv.rds")

# Extract full estimates
full_results <- SimExtract(sim_res, what = "results")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: plot simulation results
# ─────────────────────────────────────────────────────────────────────────────

df_long <- full_results %>% 
  pivot_longer(
    cols = c(sp_cor, lme4_cor, BHM_cor),
    names_to = "method",
    values_to = "estimate"
  ) %>% 
  mutate(
    method  = recode(method,
                     sp_cor   = "Spearman",
                     lme4_cor = "FHM",
                     BHM_cor  = "BHM"),
    method  = factor(method, levels = c("Spearman", "FHM", "BHM"), 
                     labels = c("Spearman", "Frequentist HM", 
                                "Bayesian HM")),
    trueCor = as.numeric(trueCor),
    rho_f  = factor(trueCor, levels = sort(unique(trueCor)))
  )

fill_cols <- c(
  Spearman = "#2A9D8F",
  `Frequentist HM` = "#E9C46A",
  `Bayesian HM` = "#7B4397"
)

line_cols <- c(
  Spearman = "#1F6F67",
  `Frequentist HM` = "#B99134",
  `Bayesian HM` = "#5A2F6E"
)

plots_lists <- vector("list", length = 9)

for (i in seq_along(plots_lists)) {
  rho_i <- i/10
  plots_lists[[i]] <- df_long |>
    filter(rho_f == as.character(rho_i)) |>
    ggplot(aes(fill = method, color = method, x = estimate)) +
    stat_slab(alpha = .25, adjust = 1.05, linewidth = .9) +
    stat_pointinterval(
      aes(y = 0, group = method),
      position = position_dodgejust(width = .20, justification = .10),
      point_size = 2.1, interval_size = 1.0,
      .width = c(.50, .80, .95)
    ) +
    geom_vline(xintercept = rho_i, linetype = "dotted", linewidth = .8, alpha = .8) +
    labs(title = bquote(rho == .(rho_i)), y = NULL, x = NULL) +
    scale_fill_manual(values = fill_cols,  name = NULL) +
    scale_colour_manual(values = line_cols, name = NULL) +
    scale_y_continuous(breaks = NULL) +
    bayesplot::theme_default(base_size = 20, base_family = "serif") %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 22, hjust = .5,
                                margin = margin(t = 0, r = 0, b = 12, l = 0))
    )
}

# Mosaico con alturas por fila
grid_3x3 <- wrap_plots(plots_lists, ncol = 3, byrow = TRUE) +
  plot_layout(heights = c(.9, 1.0, 1.2))

# Zona para la leyenda arriba del todo
legend_top <- patchwork::guide_area()

png(filename = "Figures/spearman_lme4_equivalence.png", width = 5000, height = 3500, res = 300)
(legend_top / grid_3x3) +
  plot_layout(heights = c(0.12, 1), guides = "collect") &
  theme(
    legend.position   = "top",
    legend.text       = element_text(size = 26),
    legend.title      = element_text(size = 26, face = "bold"),
    legend.key.size   = unit(1.2, "lines"),
    legend.box.margin = margin(2, 2, 8, 2),
    plot.margin = margin(10, 20, 10, 20)  # t, r, b, l (ajusta a tu gusto)
  )
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
