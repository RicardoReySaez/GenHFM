# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : clean_whitehead_data.R                                     ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 07-05-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# This script mimic Whitehead et al. (2020) data cleaning process
# Retrieved from: https://osf.io/7hp85/files/osfstorage
# File is "Analysis_ConflictCorr.Rmd
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Clean data from Experiment 1
# ─────────────────────────────────────────────────────────────────────────────

# Download raw data from OSF
# download.file("https://osf.io/c3x7j/download", 
#               destfile = "Data/Raw/whitehead_2020_Exp1_all.csv")

# Load raw data
data <- read.csv("Data/Raw/whitehead_2020_Exp1_all.csv", header = TRUE, 
                 na.strings=c(""))

# subjects to include
includesubs <- c(507,508,513,514,517,518,519,520,521,506,
                 515,523,524,525,526,527,528,529,530,532,
                 533,534,535,537,539,540,541,542,544,545,
                 546,547,548,549,551,553,555,559,560,561,
                 562,567,568,569,571,573,576,577,578,579,
                 583,584,585,586,587,588,591,592,593,594,
                 595,596,597,599,601,602,603,605,606,607,
                 608,609,610,611,613,614,615,616,617,619,
                 620,621,622,623,624,625,626,629,630,631,
                 633,634,635,637,639,640,641,642,643,645,
                 647,648,649,650,651,652,654,655,656,657,
                 658,659,660,661,662,663,664,665,666,667,
                 668,669,670,671,672,673,674,675,676,680,
                 678,679,681,683,685,687,689,690,692,693,
                 694,696,697,698,699,700,701,702,704,705,
                 706,708,710,711,712,713,715,716,718,719,
                 720,721,723,724,725,726,728,729,730,732,
                 733,734,736,737,738,739,740,741)

# Filter and clean Simon data
df.simon <- data %>% mutate(prevcon = lag(Congruency)) %>%
  mutate(StimSlideSimon.RT = as.integer(as.character(StimSlideSimon.RT)),
         StimSlideSimon.ACC = as.integer(as.character(StimSlideSimon.ACC)),
         BlockNum = as.integer(as.character(BlockNum))) %>%
  filter(Subject %in% includesubs & StimSlideSimon.RT != "" &
           (StimSlideSimon.RT > 200 & StimSlideSimon.RT < 3000) &
           StimSlideSimon.ACC == 1 & prevcon != 'NA' &
           BlockNum > 2) %>%
  mutate(RT = StimSlideSimon.RT, task = factor(1), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Filter and clean Flanker data
df.flanker <- data %>% mutate(prevcon = lag(Congruency)) %>%
  mutate(StimSlideFlanker.RT = as.integer(as.character(StimSlideFlanker.RT)),
         StimSlideFlanker.ACC = as.integer(as.character(StimSlideFlanker.ACC)),
         BlockNum = as.integer(as.character(BlockNum))) %>%
  filter(Subject %in% includesubs & StimSlideFlanker.RT != "" &
           (StimSlideFlanker.RT > 200 & StimSlideFlanker.RT < 3000) &
           StimSlideFlanker.ACC == 1 & prevcon != 'NA' &
           BlockNum > 2) %>%
  mutate(RT = StimSlideFlanker.RT, task = factor(2), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Filter and clean Stroop data
df.stroop <- data %>% mutate(prevcon = lag(Congruency)) %>%
  mutate(StimSlideStroop.RT = as.integer(as.character(StimSlideStroop.RT)),
         StimSlideStroop.RT = as.integer(as.character(StimSlideStroop.RT)),
         BlockNum = as.integer(as.character(BlockNum))) %>%
  filter(Subject %in% includesubs &  StimSlideStroop.RT != "" &
           (StimSlideStroop.RT > 200 & StimSlideStroop.RT < 3000) &
           StimSlideStroop.ACC == 1 & prevcon != 'NA' &
           BlockNum > 2) %>%
  mutate(RT = StimSlideStroop.RT, task = factor(3), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Create combined dataframe
data.cor.1 <- rbind(df.simon[,c("RT", "Congruency", "Subject", "task")],
                    df.flanker[,c("RT", "Congruency", "Subject", "task")],
                    df.stroop[,c("RT", "Congruency", "Subject", "task")])

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Clean data from Experiment 2
# ─────────────────────────────────────────────────────────────────────────────

# Download raw data from OSF: Flanker
# download.file("https://osf.io/gesr6/download", 
#               destfile = "Data/Raw/whitehead_2020_Exp2_flanker.csv")
# Download raw data from OSF: Simon
# download.file("https://osf.io/ut6nw/download", 
#               destfile = "Data/Raw/whitehead_2020_Exp2_simon.csv")
# Download raw data from OSF: Stroop
# download.file("https://osf.io/d57mg/download", 
#               destfile = "Data/Raw/whitehead_2020_Exp2_stroop.csv")

# Loading in the seperate raw datafiles
data.simon   <- read.csv("Data/Raw/whitehead_2020_Exp2_simon.csv", 
                         header = TRUE, na.strings=c(""))
data.stroop  <- read.csv("Data/Raw/whitehead_2020_Exp2_stroop.csv", 
                         header = TRUE, na.strings=c(""))
data.flanker <- read.csv("Data/Raw/whitehead_2020_Exp2_flanker.csv", 
                         header = TRUE, na.strings=c(""))


# Subjects to exclude
excludesubs <- c(115,116,126,140,148,153,160,175,
                 188,189,194,195,203,210,212,220,
                 229,233,237,239,243,250,253,297,
                 901,913,918,145,217,258,280,222)

# Filter and clean Simon data
df.simon <- data.simon %>% 
  mutate(StimSlideSimon.RT = as.integer(as.character(StimSlideSimon.RT))) %>%
  filter(!Subject %in% excludesubs &
           (StimSlideSimon.RT > 200 & StimSlideSimon.RT < 3000) &
           StimSlideSimon.ACC == 1 &
           BlockNum > 2) %>%
  mutate(RT = StimSlideSimon.RT, task = factor(1), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Filter and clean Flanker data
df.flanker <- data.flanker %>% 
  mutate(StimSlideFlanker.RT = as.integer(as.character(StimSlideFlanker.RT))) %>%
  filter(!Subject %in% excludesubs &
           (StimSlideFlanker.RT > 200 & StimSlideFlanker.RT < 3000) &
           StimSlideFlanker.ACC == 1 &
           BlockNum > 2) %>%
  mutate(RT = StimSlideFlanker.RT, task = factor(2), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Filter and clean Stroop data
df.stroop <- data.stroop %>% 
  mutate(StimSlideStroop.RT = as.integer(as.character(StimSlideStroop.RT)),
         BlockNum = as.integer(as.character(BlockNum)),
         Subject = as.integer(as.character(Subject))) %>%
  filter(!Subject %in% excludesubs & 
           (StimSlideStroop.RT > 200 & StimSlideStroop.RT < 3000) &
           StimSlideStroop.ACC == 1 &
           BlockNum > 2) %>%
  mutate(RT = StimSlideStroop.RT, task = factor(3), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Create combined dataframe
df.simon   <- select(df.simon, RT, Congruency, Subject, task)
df.flanker <- select(df.flanker, RT, Congruency, Subject, task)
df.stroop  <- select(df.stroop, RT, Congruency, Subject, task)

data.cor.2 <- rbind(df.simon,df.flanker,df.stroop)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Clean data from Experiment 3
# ─────────────────────────────────────────────────────────────────────────────

# Download raw data from OSF: Flanker
# download.file("https://osf.io/khc96/download",
#               destfile = "Data/Raw/whitehead_2020_Exp3_flanker.csv")
# Download raw data from OSF: Simon
# download.file("https://osf.io/kjyw7/download",
#               destfile = "Data/Raw/whitehead_2020_Exp3_simon.csv")
# Download raw data from OSF: Stroop
# download.file("https://osf.io/r4m9j/download",
#               destfile = "Data/Raw/whitehead_2020_Exp3_stroop.csv")

# Loading in the seperate raw datafiles
data.simon   <- read.csv("Data/Raw/whitehead_2020_Exp3_simon.csv", 
                         header = TRUE, na.strings=c(""))
data.stroop  <- read.csv("Data/Raw/whitehead_2020_Exp3_stroop.csv", 
                         header = TRUE, na.strings=c(""))
data.flanker <- read.csv("Data/Raw/whitehead_2020_Exp3_flanker.csv", 
                         header = TRUE, na.strings=c(""))


# Subjects to include
includesubs <- c(105,108,109,110,111,112,113,114,115,116,117,118,
                 119,120,121,122,123,124,125,126,127,128,129,130,
                 131,132,133,134,135,136,137,138,139,140,141,142,
                 143,144,145,146,147,148,149,150,151,152,153,154,
                 155,156,157,158,159,160,161,162,163,165,166,167,
                 168,169,170,171,172,173,174,175,176,177,178,179,
                 180,181,182,183,184,185,187,188,189,190,191,192,
                 193,194,195,196,197,198,199,200,201,202,203,204,
                 205,206,207,208,209,210,211,212,213,214,215,216,
                 217,218,219,220,221,222,223,224,225,226,228,229,
                 231,232,233,234,235,236,237,238,240,241,242,243,
                 244,246,247,248,249,250,251,253,254,255,256,257,
                 259,260,261,262,263,264,265,266,267,268,269,270,
                 271,273,274,275,277,278,279,280,281,282,283,284,
                 285,286,287,288,289,290,291,293,294,295,296,297,
                 299,301,303,305,306,307,308,309,311,312,313,314,
                 316,317,318,319,320,321,322,323,324,325,326,327,
                 328,330,331,332,333,334)

# Filter and clean Simon data
df.simon <- data.simon %>% 
  mutate(StimSlideSimon.RT = as.integer(as.character(StimSlideSimon.RT))) %>%
  filter(Subject %in% includesubs &
           (StimSlideSimon.RT > 200 & StimSlideSimon.RT < 3000) &
           StimSlideSimon.ACC == 1 &
           PracExp == "Exp")%>%
  mutate(RT = StimSlideSimon.RT, task = factor(1), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Filter and clean Flanker data
df.flanker <- data.flanker %>% 
  mutate(StimSlideFlanker.RT = as.integer(as.character(StimSlideFlanker.RT))) %>%
  filter(Subject %in% includesubs &
           (StimSlideFlanker.RT > 200 & StimSlideFlanker.RT < 3000) &
           StimSlideFlanker.ACC == 1 &
           PracExp == "Exp")%>%
  mutate(RT = StimSlideFlanker.RT, task = factor(2), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Filter and clean Stroop data
df.stroop <- data.stroop %>% 
  mutate(StimSlideStroop.RT = as.integer(as.character(StimSlideStroop.RT))) %>%
  filter(Subject %in% includesubs & 
           (StimSlideStroop.RT > 200 & StimSlideStroop.RT < 3000) &
           StimSlideStroop.ACC == 1 &
           PracExp == "Exp")%>%
  mutate(RT = StimSlideStroop.RT, task = factor(3), 
         Congruency = abs(as.numeric(as.character(Congruency))-1))

# Create combined dataframe
df.simon   <- select(df.simon, RT, Congruency, Subject, task)
df.flanker <- select(df.flanker, RT, Congruency, Subject, task)
df.stroop  <- select(df.stroop, RT, Congruency, Subject, task)

data.cor.3 <- rbind(df.simon,df.flanker,df.stroop)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Combine all experimental data frames
# ─────────────────────────────────────────────────────────────────────────────

# Add experiment and subject identificator in each data frame
data.cor.1.test <- data.cor.1 %>%
  select(Subject, Congruency, RT, task) %>% 
  mutate(experiment = factor(1), 
         Subject = factor(as.numeric(as.character(Subject))+0))
data.cor.2.test <- data.cor.2 %>% 
  select(Subject, Congruency, RT, task) %>% 
  mutate(experiment = factor(2), 
         Subject = factor(as.numeric(as.character(Subject))+1000))
data.cor.3.test <- data.cor.3 %>% 
  select(Subject, Congruency, RT, task) %>% 
  mutate(experiment = factor(3), 
         Subject = factor(as.numeric(as.character(Subject))+2000))

# Combine all data frames and filter data based on RTs
whitehead.data <- rbind(data.cor.1.test,data.cor.2.test,data.cor.3.test) %>% 
  filter(RT > 200, RT < 3000)

# Save cleaned data frame
save(whitehead.data, file = "Data/Processed/whitehead_data.rdata")

# ─────────────────────────────────────────────────────────────────────────────
