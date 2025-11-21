##############################################
# Set-up
##############################################

# Libraries
library(openxlsx)
library(dplyr)
library(ggplot2)
library(scales)
library(limma)

# Functions
source("functions_GMM_integration.R")

##############################################
# Fergus' data
##############################################

# Load
fergus <- read.table("combined.raw.tsv", header = T)

# Add information
fergus$POS <- fergus$POS + 1
fergus$uPOS <- paste(fergus$POS, fergus$REF, fergus$ALT, sep = "_")
fergus$Splicing <- ifelse(fergus$SpliceAI_DS_AG > 0.2 | fergus$SpliceAI_DS_DG > 0.2, "Yes", "No")

# Compute D5 / lib
fergus <- fergus %>%
  
  mutate(across(starts_with("R") & ends_with("D5"), ~ . / fergus[[sub("D5", "lib", cur_column())]], .names = "{sub('_D5', '_PRE', col)}"))
# Compute D14 / lib
fergus <- fergus %>%
  mutate(across(starts_with("R") & ends_with("D14"), ~ . / fergus[[sub("D14", "lib", cur_column())]], .names = "{sub('_D14', '_POST', col)}"))

# Compute ratios
fergus <- fergus %>%
  mutate(across(matches("^R\\d+_POST$"), ~ . / fergus[[sub("_POST", "_PRE", cur_column())]], .names = "{sub('_POST', '_ratio', col)}"))

# Convert outcome
fergus <- fergus %>%
  mutate(Outcome = case_when(EventType == "Synonymous" ~ "Synonymous",
                             EventType == "StopGain" ~ "Non-sense",
                             EventType == "Missense" ~ "Non-synonymous",
                             TRUE ~ "Intronic"))

# Normalize
fergus <- fergus %>%
  mutate(expected_clinical = case_when(Outcome == "Synonymous" ~ "Functional",
                                       Outcome == "Non-sense" ~ "Non-functional",
                                       TRUE ~ "Unknown"))

# If splicing effect, remove from training
fergus$expected_clinical <- ifelse(fergus$Splicing == "Yes" & fergus$expected_clinical == "Functional", "Unknown", fergus$expected_clinical)

# By-experiment normalization
columns_to_normalize <- c("R1_ratio", "R2_ratio", "R3_ratio", "R4_ratio", "R5_ratio", "R6_ratio")

fergus <- fergus %>%
  group_by(exon) %>%
  do({
    temp_data <- .
    for (col in columns_to_normalize) {
      if (all(is.na(temp_data[[col]]))) {
        temp_data[[paste0(col, "_norm")]] <- NA
      } else {
        temp_data[[col]] <- as.numeric(temp_data[[col]])
        temp_data[[paste0(col, "_norm")]] <- normalization_per_exp(temp_data, col)
      }
    }
    temp_data
  }) %>%
  ungroup()

fergus$chr_pos <- as.numeric(gsub("_.*", "", fergus$uPOS))

fergus <- fergus %>%
  group_by(exon) %>%
  mutate(pos = chr_pos - min(chr_pos) + 1)

# Apply position biais for each experiment
columns_to_normalize <- c("R1_ratio_norm", "R2_ratio_norm", "R3_ratio_norm", "R4_ratio_norm", "R5_ratio_norm", "R6_ratio_norm")

fergus <- fergus %>%
  group_by(exon) %>%
  do({
    temp_data <- .
    for (col in columns_to_normalize) {
      if (all(is.na(temp_data[[col]]))) {
        temp_data[[paste0(col, "_sns")]] <- NA
      } else {
        temp_data[[paste0(col, "_sns")]] <- position_biais(temp_data, col)
      }
    }
    temp_data
  }) %>%
  ungroup()

fergus <- fergus %>%
  select(-expected_clinical, -chr_pos, -pos)

##############################################
# AVENGERS' data
##############################################

avengers <- read.xlsx("../results/BRCA2_all_counts.xlsx")

# Add splicing information
avengers$uPOS <- gsub(".*g\\.", "", avengers$g_nom)
avengers$uPOS <- gsub("(\\d+)([A-Z])>([A-Z])", "\\1_\\2_\\3", avengers$uPOS)
avengers$chr_pos <- as.numeric(gsub("_.*", "", avengers$uPOS))
splicing <- read.xlsx("../reference/all_exons_summary_spliceAI.xlsx")
splicing$uPOS <- gsub(".*g\\.", "", splicing$g.nom)
splicing$uPOS <- gsub("(\\d+)([A-Z])>([A-Z])", "\\1_\\2_\\3", splicing$uPOS)

splicing <- splicing %>%
  group_by(uPOS) %>%
  summarize(
    Affect_splicing = case_when(
      all(Affect_splicing == "Yes") ~ "Yes",
      all(Affect_splicing == "No") ~ "No",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )

avengers <- merge(avengers, splicing[, c("Affect_splicing", "uPOS")], by = "uPOS", all.x = T)

# Compute ratios
avengers <- avengers %>%
  group_by(Exon) %>%
  mutate(pos = chr_pos - min(chr_pos) + 1)

avengers$R7_PRE <- ifelse(avengers$D3_freq_RP1 > 1e-5, avengers$D3_freq_RP1, NA)
avengers$R8_PRE <- ifelse(avengers$D3_freq_RP2 > 1e-5, avengers$D3_freq_RP2, NA)
avengers$R7_POST <- ifelse(avengers$DMSO_D14_freq_RP1 > 0, avengers$DMSO_D14_freq_RP1, 0.00001)
avengers$R8_POST <- ifelse(avengers$DMSO_D14_freq_RP2 > 0, avengers$DMSO_D14_freq_RP2, 0.00001)
avengers$R7_ratio <- avengers$R7_POST / avengers$R7_PRE
avengers$R8_ratio <- ifelse(avengers$Experiment != "17.1", avengers$R8_POST / avengers$R8_PRE, avengers$R7_POST / avengers$R7_PRE)

# Add data for normalization and training
avengers <- avengers %>%
  mutate(expected_clinical = case_when((Outcome == "Synonymous") & !(ClinVar %in% c("Pathogenic", "Likely pathogenic")) ~ "Functional",
                                       (Outcome == "Non-sense") & !(ClinVar %in% c("Benign", "Likely benign")) ~ "Non-functional",
                                       TRUE ~ "Unknown"))

avengers$expected_clinical <- ifelse(avengers$Affect_splicing == "Yes" & avengers$expected_clinical == "Functional", "Unknown", avengers$expected_clinical)

# Normalize ratios
columns_to_normalize <- c("R7_ratio", "R8_ratio")

avengers <- avengers %>%
  group_by(Experiment) %>%
  do({
    temp_data <- .
    for (col in columns_to_normalize) {
      if (all(is.na(temp_data[[col]]))) {
        temp_data[[paste0(col, "_norm")]] <- NA
      } else {
        temp_data[[col]] <- as.numeric(temp_data[[col]])
        temp_data[[paste0(col, "_norm")]] <- normalization_per_exp(temp_data, col)
      }
    }
    temp_data
  }) %>%
  ungroup()

# Apply position biais for each experiment
columns_to_normalize <- c("R7_ratio_norm", "R8_ratio_norm")

avengers <- avengers %>%
  group_by(Exon) %>%
  do({
    temp_data <- .
    for (col in columns_to_normalize) {
      if (all(is.na(temp_data[[col]]))) {
        temp_data[[paste0(col, "_sns")]] <- NA
      } else {
        temp_data[[paste0(col, "_sns")]] <- position_biais(temp_data, col)
      }
    }
    temp_data
  }) %>%
  ungroup()

avengers <- avengers %>%
  select(-expected_clinical, -chr_pos, -pos)

##############################################
# Merge and process
##############################################

# Merge
merged <- merge(fergus, avengers, by = c("uPOS", "Outcome"), all = T)
merged <- merged[order(merged$uPOS),]
merged$Splicing <- ifelse(merged$Splicing == "Yes" | merged$Affect_splicing == "Yes", "Yes", "No")

# Get mean of replicates
merged$post_pre_ratio_sns_loess <- rowMeans(merged[, c("R1_ratio_norm_sns", "R2_ratio_norm_sns", 
                                                       "R3_ratio_norm_sns", "R4_ratio_norm_sns", 
                                                       "R5_ratio_norm_sns", "R6_ratio_norm_sns",
                                                       "R7_ratio_norm_sns", "R8_ratio_norm_sns")],
                                            na.rm = T)

# Add exon information
merged$Exon <- ifelse(is.na(merged$Exon), merged$exon, merged$Exon)
merged$Exon <- gsub("[^0-9.]", "", merged$Exon)

merged$exon <- NULL

# Add expected clinical
merged <- merged %>%
  mutate(expected_clinical = case_when(Outcome == "Synonymous" ~ "Functional",
                                       Outcome == "Non-sense" ~ "Non-functional",
                                       TRUE ~ "Unknown"))


merged$expected_clinical <- ifelse(merged$Splicing == "Yes" & merged$expected_clinical == "Functional", "Unknown", merged$expected_clinical)

# Clean
merged_clean <- merged[, c("uPOS", "AApos", "RefAA", "AltAA", "Exon", "Outcome", "expected_clinical", "Splicing", "c_nom", "g_nom", "AA_change", "dbSNP.ID", "ClinVar", "post_pre_ratio_sns_loess", "post_pre_ratio_sns_loess")]

# Save
write.xlsx(merged_clean, "Fergus_Avengers_merged_data.xlsx")
