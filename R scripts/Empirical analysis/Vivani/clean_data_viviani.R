# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : clean_vivianni_data.R                                      ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 16-05-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# This script mimic Vivianni et al. (2024) data cleaning process
# Retrieved from: https://osf.io/5sm9j/files/osfstorage
# File is "ANALYSIS/Matlab/ManyStroopScript.m
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(openxlsx) # read Excel files
library(dplyr)    # data manipulation
library(tidyr)    # data manipulation
library(stringr)  # data manipulation
library(glmmTMB)  # Linear mixed models (with parallelized computation)
library(lme4)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-defined functions
# ─────────────────────────────────────────────────────────────────────────────

# Function to compute experimental effects with linear mixed model
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

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Clean data from Viviani et al
# ─────────────────────────────────────────────────────────────────────────────

# Download raw data from OSF
# download.file("https://osf.io/6qx98/download", mode = "wb",
#               destfile = "Data/Raw/Vivianni_2024_ManyStroopData.xlsx")

df <- read.xlsx(xlsxFile = "Data/Raw/Viviani_2024_ManyStroopData.xlsx") 

# Filter by experimental trial and remove block prefix
df <- df |> 
  filter(ExpTrial == 1) |> 
  mutate(Block = substr(Block, 5, nchar(Block)))

# Compute SSok flag per subkect
ss_ok <- df %>%
  group_by(SSID) %>%
  summarise(
    mean_iRT = mean(iRTok, na.rm = TRUE),
    mean_acc = mean(Correct, na.rm = TRUE)
  ) %>%
  mutate(
    zRT  = as.numeric(scale(mean_iRT)),
    ssOK = (abs(zRT) < 3) & (mean_acc > 0.75)
  ) %>%
  dplyr::select(SSID, ssOK)

# Join both data frames and filter by valid subjects
df_clean <- df %>%
  left_join(ss_ok, by = "SSID") %>%
  filter(
    !is.na(iRTok),
    !is.na(iRTpre),
    ssOK.x == 1,
    Correct == 1
  )

# Prepare variables for linear mixed models
df_clean <- df_clean %>%
  mutate(
    SubBlock_z = as.numeric(scale(SubBlockNum)),
    Trial_z    = as.numeric(scale(Trial)),
    iRTpre_z   = as.numeric(scale(iRTpre)),
    hRespC     = factor(RelevH, levels = c(-1, 1)),
    vRespC     = factor(RelevV, levels = c(-1, 1)),
    postERRC   = factor(PostERR, levels = c(0, 1)),
    CongC      = factor(CONG, levels = c(1, 0)),
    Block = factor(Block, levels = c("Peripheral", "Perifoveal",
                                     "Navon", "FigureGround",
                                     "Flanker", "Saliency"))
  )

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Estimate raw experimental effect per subject and task
# In the paper, these are GLM scores
# ─────────────────────────────────────────────────────────────────────────────

# Compute aggregated means for each subject in each condition and task
agg_means <-  tapply(df_clean$iRTok, list(df_clean$SSID, df_clean$Block, df_clean$CongC), mean)

# Finally, compute average experimental effect
GLM_scores_iRt <- agg_means[,,2] - agg_means[,,1]

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Estimate experimental effects with linear mixed model
# In the paper, these are LMM scores
# ─────────────────────────────────────────────────────────────────────────────

# Fit the full-linear mixed model
m_iRT <- glmmTMB(formula = iRTok ~ Trial_z*SubBlock_z + iRTpre_z + postERRC
                         + hRespC + vRespC + Block*CongC + (Block*CongC|SSID), 
                 data = df_clean, 
                 REML = FALSE, 
                 control = glmmTMBControl(
                   # Parallelized optimazion
                   parallel = 8, optimizer = optim,
                   optArgs = list(method = "L-BFGS-B"),
                   optCtrl=list(maxit=5e3)))

# Detect and exclude standardized outliers
zres <- as.numeric(scale(residuals(m_iRT)))
keep <- abs(zres) <= 3

# Refit the model without these outliers
m_iRT_trim <- update(m_iRT, data = df_clean[keep, ])

# Store fixed-and-random effects
fix_eff <- fixef(m_iRT_trim)$cond                       
rnd_eff <- ranef(m_iRT_trim)[[1]]$SSID |>               
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
LMM_scores_iRT <- I_scores - C_scores

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Estimate experimental effects with subject-specific OLS regression
# In the paper, these are RCA scores
# ─────────────────────────────────────────────────────────────────────────────

# Compute a linear mixed model for each subject
RCA_scores_iRT <- do.call(
  rbind,
  lapply(unique(df_clean$SSID), function(x) {
    # Fit OLS linear regression
    lm.fit <- lm(iRTok ~ Trial_z + iRTpre_z + postERRC + hRespC + vRespC + Block*CongC, 
                 data = df_clean[which(df_clean$SSID == x),])
    
    # Store experimental effects per tasl
    coefs <- data.frame(
      Peripheal    = coef(lm.fit)["CongC0"],
      Perifoveal   = coef(lm.fit)["CongC0"] + coef(lm.fit)["BlockPerifoveal:CongC0"],
      Navon        = coef(lm.fit)["CongC0"] + coef(lm.fit)["BlockNavon:CongC0"],
      FigureGround = coef(lm.fit)["CongC0"] + coef(lm.fit)["BlockFigureGround:CongC0"],
      Flanker      = coef(lm.fit)["CongC0"] + coef(lm.fit)["BlockFlanker:CongC0"],
      Saliency     = coef(lm.fit)["CongC0"] + coef(lm.fit)["BlockSaliency:CongC0"]
    )
    # Return experimental effects
    rownames(coefs) <- NULL
    return(coefs)
  })
)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Save processed data frames
# ─────────────────────────────────────────────────────────────────────────────

# Cleaned long data frame
vivianni_2024_df <- df_clean[, c("SSID", "Block", "CONG", "RT")]
# Scale to seconds and change column names
vivianni_2024_df$RT <- vivianni_2024_df$RT/1000
colnames(vivianni_2024_df) <- c("subject", "task", "condition", "RT")
# Bind experimental effects scores
vivianni_2024_effects <- cbind(
  GLM_scores_iRt, 
  LMM_scores_iRT, 
  LMMb_scores_iRT,
  RCA_scores_iRT
)
# Change column names
colnames(vivianni_2024_effects) <- paste0(task_names, rep(c("_GLM", "_LMM", "_RCA"), each = 6))

# Save data frames
save(vivianni_2024_df,      file = "Data/Processed/viviani_data.rdata")
save(vivianni_2024_effects, file = "Data/Processed/viviani_effects_data.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Structural Equation Modeling Analysis
# ─────────────────────────────────────────────────────────────────────────────

# Store column IDs
LMM_cols <- grep("LMM", colnames(vivianni_2024_effects))
GLM_cols <- grep("GLM", colnames(vivianni_2024_effects))
RCA_cols <- grep("RCA", colnames(vivianni_2024_effects))

# Latent dimensionalty assessment: EGA
EGAnet::EGA(vivianni_2024_effects[,LMM_cols])
EGAnet::EGA(vivianni_2024_effects[,GLM_cols])
EGAnet::EGA(vivianni_2024_effects[,RCA_cols])

# Multivariate normality assumption
psych::mardia(vivianni_2024_effects[,LMM_cols])
psych::mardia(vivianni_2024_effects[,GLM_cols])
psych::mardia(vivianni_2024_effects[,RCA_cols])

# Kaiser-Meyen-Olkin test
psych::KMO(vivianni_2024_effects[,LMM_cols])
psych::KMO(vivianni_2024_effects[,GLM_cols])
psych::KMO(vivianni_2024_effects[,RCA_cols])

# Bartlett test
psych::cortest.bartlett(vivianni_2024_effects[,LMM_cols])
psych::cortest.bartlett(vivianni_2024_effects[,GLM_cols])
psych::cortest.bartlett(vivianni_2024_effects[,RCA_cols])

# Exploratory Factor analysis
psych::fa(vivianni_2024_effects[,LMM_cols], nfactors = 2, 
          fm = "minres", rotate = "oblimin")
psych::fa(vivianni_2024_effects[,GLM_cols], nfactors = 2, 
          fm = "minres", rotate = "oblimin")
psych::fa(vivianni_2024_effects[,RCA_cols], nfactors = 2, 
          fm = "minres", rotate = "oblimin")

# Model fit comparison
lavaan::summary(lavaan::efa(data = vivianni_2024_effects[,LMM_cols], nfactors = 2, 
                            rotation = "oblimin", output = "lavaan", estimator = "MLR"),
                standardized = TRUE, fit.measures = TRUE)
lavaan::summary(lavaan::efa(data = vivianni_2024_effects[,GLM_cols], nfactors = 2, 
                            rotation = "oblimin", output = "lavaan", estimator = "MLR"),
                standardized = TRUE, fit.measures = TRUE)
lavaan::summary(lavaan::efa(data = vivianni_2024_effects[,RCA_cols], nfactors = 2, 
                            rotation = "oblimin", output = "lavaan", estimator = "MLR"),
                standardized = TRUE, fit.measures = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
