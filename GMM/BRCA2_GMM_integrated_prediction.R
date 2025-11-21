# Libraries
library(openxlsx)
library(ggplot2)
library(dplyr)
library(readr)
library(stats)
library(mclust)

# Functions
source("BRCA2/GMM_integration/functions_GMM_integrated.R")

# Set-up priorP
priorP <- 0.1

# Load data
data_df <- read.xlsx("BRCA2/GMM_integration/Fergus_Avengers_merged_data.xlsx")
data_df$AA_change <- ifelse(is.na(data_df$AA_change), paste0(data_df$RefAA, data_df$AApos, data_df$AltAA), data_df$AA_change) #Update AA change for variants in only one group
data_df$ClinVar <- ifelse(data_df$AA_change %in% c("T3013I", "K3015E"), "Benign", data_df$ClinVar) #Update 2 ClinVar vaitans

# Normalize per exon
data_df <- data_df %>%
  group_by(Exon) %>%
  mutate(post_pre_ratio_sns_loess_norm = per_exon_norm(cur_data())) %>%  
  ungroup()

# Global normalization
exon_medians <- data_df %>%
  group_by(Exon, Outcome) %>%
  summarise(Median_FS = median(post_pre_ratio_sns_loess_norm, na.rm = TRUE), .groups = "drop")

global_medians <- data_df %>%
  group_by(Outcome) %>%
  summarise(Global_Median_FS = median(post_pre_ratio_sns_loess_norm, na.rm = TRUE), .groups = "drop")

scaling_factors <- exon_medians %>%
  inner_join(global_medians, by = "Outcome") %>%
  mutate(Scaling_Factor = Global_Median_FS / Median_FS) %>%
  select(Exon, Outcome, Scaling_Factor)

data_df <- data_df %>%
  left_join(scaling_factors, by = c("Exon", "Outcome")) %>%
  mutate(Scaling_Factor = ifelse(is.na(Scaling_Factor), 1, Scaling_Factor),  # No scaling for other types
         final_norm = post_pre_ratio_sns_loess_norm * Scaling_Factor) %>%
  select(-Scaling_Factor)

data_df <- data_df[!is.na(data_df$final_norm),]

# Training data
training_df <- data_df[data_df$AA_change %in% c("R2502H", "D2665G", "I2675V", "V2728I", "A2770T", "S2835P", "A2912T", "S2922G", "T3013I", "K3015E", "K3059E", "N3124I"),]
training_df$lof <- ifelse(training_df$ClinVar == "Pathogenic", 1, 0)

# Run mclust
x_mclust.train <- MclustDA(data = training_df$final_norm, class = training_df$lof, G = 1:2, modelNames = c("V"), prop = c(1 - priorP, priorP))
x_mclust.predict <- predict.MclustDA(x_mclust.train, newdata = data_df$final_norm, newclass = training_df$lof, prop = c(1 - priorP, priorP), G = 1:2, modelNames = c("V"))
data_df$proba <- x_mclust.predict$z[,2]

# Classfication using prior-p
prior_tab <- read.xlsx("BRCA2/reference/scores/priorR.xlsx")

prior_vect <- as.vector(unlist(prior_tab[prior_tab$priorP == priorP, c(-1)]))

data_df$odd <- (data_df$proba * (1 - priorP)) / ((1 - data_df$proba) * priorP)

data_df <- data_df %>%
  mutate(odd.general.class = case_when(odd >= prior_vect[8] ~ "Pathogenic Strong",
                                       odd < prior_vect[8] & odd >= prior_vect[7] ~ "Pathogenic Strong",
                                       odd < prior_vect[7] & odd >= prior_vect[6] ~ "Pathogenic Moderate",
                                       odd < prior_vect[6] & odd >= prior_vect[5] ~ "Pathogenic Supporting",
                                       odd <= (prior_vect[1]) ~ "Benign Strong",
                                       odd > (prior_vect[1]) & odd <= (prior_vect[2]) ~ "Benign Strong",
                                       odd > (prior_vect[2]) & odd <= (prior_vect[3]) ~ "Benign Moderate",
                                       odd > (prior_vect[3]) & odd<= (prior_vect[4]) ~ "Benign Supporting",
                                       TRUE ~ "Uncertain"))

# Add "General" group for plotting and stats
data_df$group.general <- ifelse(data_df$odd.general.class %in% c("Pathogenic Strong", "Pathogenic Moderate", "Pathogenic Supporting"), "P", "VUS")
data_df$group.general <- ifelse(data_df$odd.general.class %in% c("Benign Strong", "Benign Moderate", "Benign Supporting"), "B", data_df$group.general)

# Save
write.xlsx(data_df, "BRCA2/GMM_integration/GMM_integrated_results.xlsx")


