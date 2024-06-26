###########################################################
#                                                         #
#     CT-TIME response monitoring              #
#                                                         #
###########################################################
# Load necessary packages
library(here)       # Generate paths
library(tidyverse)  # Data frame operations
library(caret)      # training models
library(gtsummary)  # Univariate analysis
library(survival)
library(survminer)
library(ggplot2)
library(tibble)
library(lmtest)
library(dplyr)
library(ggsankey)
###### Define paths and load data ######
bs_path <- here::here("..","1_data_processed", "raw_NANO_IMMUNOMIX_data.Rda") #input data
fu_path <- here::here("..","1_data_processed", "raw_NANO_IMMUNOMIX_timepoints.Rda") #input data
response_path <- here::here("..","0_data","IMMUNOMIX_NEW","Total_Immuno_Clinical_Database_20230927.csv") #input data
signature_path <- here::here("..","2_training","CT-TIME_signature.Rda") #input data
model_path <- here::here( "..","2_training","CT-TIME_model.rds") #input model
load(bs_path)
bs <- data_all[data_all$DataBatch=="IMMUNOMIX_Soft",]
bs$Timepoint <- "BS"
load(fu_path)
fu <- data_all[data_all$DataBatch=="IMMUNOMIX_Soft",]
fu <- subset(fu,Timepoint=="FU")
br <- data_all[data_all$DataBatch=="IMMUNOMIX_Soft",]
br <- subset(br,Timepoint=="BR")
response <- read.csv2(response_path)
load(signature_path)
glmnet_fit <- readRDS(model_path)
#########
bs[c(4:16,18:110)] <- scale(bs[c(4:16,18:110)])
fu[c(4:16,18:110)] <- scale(fu[c(4:16,18:110)])
br[c(4:16,18:110)] <- scale(br[c(4:16,18:110)])
############ Feature filtering ############
bs_features <- bs[,c(selected_feature_names,"Patient","LesionID","Timepoint","original_shape_VoxelVolume")]
bs_id <- data.frame(bs_features)
bs_id$CT_TIME <- predict(glmnet_fit,bs_id ,type = "prob")[,1]
bs_id$pred  <- predict(glmnet_fit,bs_id ,type = "raw")
bs <- bs_id[,c("Patient","LesionID","Timepoint","CT_TIME","pred","original_shape_VoxelVolume")]

fu_features <- fu[,c(selected_feature_names,"Patient","LesionID","Timepoint","original_shape_VoxelVolume")]
fu_id <- data.frame(fu_features)
fu_id$CT_TIME <- predict(glmnet_fit,fu_id ,type = "prob")[,1]
fu_id$pred  <- predict(glmnet_fit,fu_id ,type = "raw")
fu <- fu_id[,c("Patient","LesionID","Timepoint","CT_TIME","pred","original_shape_VoxelVolume")]

br_features <- br[,c(selected_feature_names,"Patient","LesionID","Timepoint","original_shape_VoxelVolume")]
br_id <- data.frame(br_features)
br_id$CT_TIME <- predict(glmnet_fit,br_id ,type = "prob")[,1]
br_id$pred  <- predict(glmnet_fit,br_id ,type = "raw")
br <- br_id[,c("Patient","LesionID","Timepoint","CT_TIME","pred","original_shape_VoxelVolume")]

#######
# df <- rbind(bs,fu,br)
# df <- merge(df,response[,c(4,7,22,24)], by.x="Patient",by.y="SAP",all=FALSE)
# df$clinical_benefit <- ifelse(df$PFS_months>=5,"Yes","No") 


all <- full_join(bs, fu, by = c("Patient", "LesionID"))
all <- full_join(all, br, by = c("Patient", "LesionID"))
all <- na.omit(all)

########## CT-TIME aggregation
aggregated_results <- all %>%
  group_by(Patient) %>%
  summarise(
    min_CT_TIME_BS = min(CT_TIME.x, na.rm = TRUE),
    max_CT_TIME_BS = max(CT_TIME.x, na.rm = TRUE),
    mean_CT_TIME_BS = mean(CT_TIME.x, na.rm = TRUE),
    min_CT_TIME_FU = min(CT_TIME.y, na.rm = TRUE),
    max_CT_TIME_FU = max(CT_TIME.y, na.rm = TRUE),
    mean_CT_TIME_FU = mean(CT_TIME.y, na.rm = TRUE),
    min_CT_TIME_BR = min(CT_TIME, na.rm = TRUE),
    max_CT_TIME_BR = max(CT_TIME, na.rm = TRUE),
    mean_CT_TIME_BR = mean(CT_TIME, na.rm = TRUE),
    sum_vol_BS = sum(original_shape_VoxelVolume.x, na.rm = TRUE),
    sum_vol_FU = sum(original_shape_VoxelVolume.y, na.rm = TRUE),
    sum_vol_BR = sum(original_shape_VoxelVolume, na.rm = TRUE)
  )

aggregated_results <- merge(aggregated_results,response[,c(1,4,7,14,22,24)], by.x="Patient",by.y="SAP")

# Create the new column based on the dictionary
aggregated_results$primary_tumor <- ifelse(
  aggregated_results$primary_tumor_localization %in% c("adrenal", "thymus", "thyroid","paraganglioma retroperitoneal"), "endocrine",
  ifelse(
    aggregated_results$primary_tumor_localization %in% c("pancreas", "liver", "biliary_tract", "cholangiocarcinoma", "hepatocarcinoma"), "hepatobiliary",
    ifelse(
      aggregated_results$primary_tumor_localization %in% c("colon", "esophagus", "gastric", "rectal", "colorectal", "UGE","stomach"), "gastrointestinal",
      ifelse(
        aggregated_results$primary_tumor_localization %in% c("kidney", "bladder", "prostate", "penis", "renal"), "genitourinary",
        ifelse(
          aggregated_results$primary_tumor_localization %in% c("cervix", "ovarian", "endometrial", "endometrium", "ovary", "ureteral", "urothelial"), "gynecologic",
          ifelse(
            aggregated_results$primary_tumor_localization %in% c("melanoma", "skin_uveal", "skin"), "skin",
            ifelse(
              aggregated_results$primary_tumor_localization %in% c("retroperitoneal", "bone", "sacral chordoma"), "sarcoma",
              ifelse(
                aggregated_results$primary_tumor_localization %in% c("breast TN", "breast/ovary", "breast"), "breast",
                ifelse(
                  aggregated_results$primary_tumor_localization %in% c("head_and_neck", "H&N"), "head&neck",
                  ifelse(
                    aggregated_results$primary_tumor_localization %in% c("mesothelioma", "pleura","lung"), "lung", NA
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

aggregated_results <- aggregated_results %>%
  distinct(Patient, .keep_all = TRUE)
names(aggregated_results)
aggregated_results$clinical_benefit <- ifelse(aggregated_results$PFS_months>=5,"Yes","No") 
aggregated_results$min_CT_TIME_change <- ((aggregated_results$min_CT_TIME_FU-aggregated_results$min_CT_TIME_BS))
aggregated_results$max_CT_TIME_change <- ((aggregated_results$max_CT_TIME_FU-aggregated_results$max_CT_TIME_BS))
aggregated_results$mean_CT_TIME_change <- ((aggregated_results$mean_CT_TIME_FU-aggregated_results$mean_CT_TIME_BS))
aggregated_results$VolumeChange <- 100*((aggregated_results$sum_vol_FU-aggregated_results$sum_vol_BS)/aggregated_results$sum_vol_BS)
aggregated_results$VolumetricResponse <- ifelse(100*((aggregated_results$sum_vol_FU-aggregated_results$sum_vol_BS)/aggregated_results$sum_vol_BS)>=0, "Progression","Regression")

# Earlier detection of progression?
vol_prog <- subset(aggregated_results,VolumetricResponse=="Progression")
exclude <- ifelse(vol_prog$Patient%in%c("10044748","11411793","16999701"),TRUE,FALSE)
vol_prog <- vol_prog[!exclude,]
vol_prog$min_CT_TIME_FU_dich <- ifelse(vol_prog$min_CT_TIME_FU>0.5,"inflamed","uninflamed")

# Create a boxplot for biomarker scores at follow-up for progressing patients
plot1 <- ggplot(vol_prog, aes(x = min_CT_TIME_FU_dich, y = PFS_months, fill = min_CT_TIME_FU_dich)) +
  geom_boxplot(color = "black") +  # Specify the outline color of the boxes
  geom_jitter(aes(color = primary_tumor), position = position_jitter(width = 0.3), size = 2, alpha = 0.5) +  # Add jitter with color
  scale_fill_manual(values = c("#b2182b", "#969696")) +  # Specify fill colors
  labs(title = "CT-TIME at Follow-up for Progressors",
       x = "Aggregated CT-TIME (dichotomized)",
       y = "Progression Free Survival (months)") +
  theme(legend.position = "bottom")  # Move legend to the bottom
# Perform t-test for biomarker scores
t_test_scores <- t.test(PFS_months ~ min_CT_TIME_FU_dich, data = vol_prog)
p_value_scores <- format.pval(t_test_scores$p.value, digits = 2)

# Add p-value to the plot
plot1 +
  annotate("text", x = 1.5, y = max(vol_prog$PFS_months),
           label = paste("p =", p_value_scores), hjust = 0, vjust = 1, color = "black")

# 
# 
# # Create a boxplot for biomarker scores at follow-up for progressing patients
# plot1 <- ggplot(vol_prog, aes(x = clinical_benefit, y = min_CT_TIME_FU, color = primary_tumor)) +
#   geom_boxplot(fill = "lightblue", color = "blue") +
#   geom_jitter(position = position_jitter(width = 0.3), size = 2, alpha = 0.5) +  # Add jitter
#   labs(title = "CT-TIME at Follow-up for Progressing Patients",
#        x = "Clinical Benefit (0: No Clinical Benefit, 1: Clinical Benefit)",
#        y = "Biomarker Score at Follow-up")+
#   theme(legend.position = "bottom")  # Move legend to the bottom
# 
# # Perform t-test for biomarker scores
# t_test_scores <- t.test(min_CT_TIME_FU ~ clinical_benefit, data = vol_prog)
# p_value_scores <- format.pval(t_test_scores$p.value, digits = 3)
# 
# # Add p-value to the plot
# plot1 +
#   annotate("text", x = 1.5, y = max(vol_prog$min_CT_TIME_FU),
#            label = paste("p =", p_value_scores), hjust = 0, vjust = 1, color = "black")

#
# # Create a boxplot for biomarker changes for progressing patients
# ggplot(vol_prog, aes(x = clinical_benefit, y = median_CT_TIME_change)) +
#   geom_boxplot(fill = "lightgreen", color = "darkgreen") +
#   labs(title = "Biomarker Changes for Progressing Patients",
#        x = "Clinical Benefit (0: No Clinical Benefit, 1: Clinical Benefit)",
#        y = "Biomarker Change from Baseline to Follow-up")

# Create a boxplot for biomarker scores at follow-up for progressing patients
plot1 <- ggplot(vol_prog, aes(x = trial_best_response, y = min_CT_TIME_FU, color = primary_tumor)) +
  geom_boxplot(fill = "lightblue", color = "blue") +
  geom_jitter(position = position_jitter(width = 0.3), size = 2, alpha = 0.5) +  # Add jitter
  labs(title = "Biomarker Scores at Follow-up for Progressing Patients",
       x = "Clinical Benefit (0: No Clinical Benefit, 1: Clinical Benefit)",
       y = "Biomarker Score at Follow-up")+
  theme(legend.position = "bottom")  # Move legend to the bottom

# Perform ANOVA for biomarker scores
anova_scores <- aov(min_CT_TIME_FU ~ trial_best_response, data = vol_prog)
p_value_scores <- summary(anova_scores)[["Pr(>F)"]][1]

# Add p-value to the plot
plot1 +
  annotate("text", x = length(levels(vol_prog$trial_best_response))/2, 
           y = max(vol_prog$min_CT_TIME_FU), 
           label = paste("p =", format.pval(p_value_scores, digits = 3)), hjust = 0, vjust = 1, color = "black")


aggregated_results$min_CT_TIME_BS_status <- ifelse(aggregated_results$min_CT_TIME_BS>0.5,"inflamed","uninflamed")
aggregated_results$min_CT_TIME_FU_status<- ifelse(aggregated_results$min_CT_TIME_FU>0.5,"inflamed","uninflamed")
aggregated_results$min_CT_TIME_BR_status <- ifelse(aggregated_results$min_CT_TIME_BR>0.5,"inflamed","uninflamed")

########## new sankey
# Step 1
df2 <- aggregated_results %>%
  make_long(min_CT_TIME_BS_status, min_CT_TIME_FU_status, min_CT_TIME_BR_status)

# Step 2
dagg <- df2 %>%
  dplyr::group_by(node,x) %>%
  tally()

# Step 3
df3 <- merge(df2, dagg, by.x = c('node','x'), by.y = c('node','x'), all.x = TRUE)

# Set colors
cols <- c("inflamed" = "#b2182b", "uninflamed" = "#969696")

# Create Sankey plot
pl <- ggplot(df3, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = factor(node),
  label = paste0(" n=", n)
)) +
  geom_sankey(
    flow.alpha = 0.5,
    color = "gray40",
    show.legend = TRUE,
    aes(fill = factor(node))
  ) +
  geom_sankey_label(
    size = 3,
    color = "white",
    fill = "gray40",
    hjust = -0.2
  ) +
  scale_fill_manual(values = cols) +  # Set node colors
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Sankey diagram of CT-TIME status")

pl


# Count changes in status
num_changes_inflamed_to_uninflamed <- sum(aggregated_results$min_CT_TIME_BS_status == "inflamed" & aggregated_results$min_CT_TIME_FU_status == "uninflamed")
num_changes_uninflamed_to_inflamed <- sum(aggregated_results$min_CT_TIME_BS_status == "uninflamed" & aggregated_results$min_CT_TIME_FU_status == "inflamed")
num_no_changes <- sum(aggregated_results$min_CT_TIME_BS_status == aggregated_results$min_CT_TIME_FU_status)


num_changes_inflamed_to_uninflamed_br <- sum(aggregated_results$min_CT_TIME_FU_status == "inflamed" & aggregated_results$min_CT_TIME_BR_status == "uninflamed")
num_changes_uninflamed_to_inflamed_br <- sum(aggregated_results$min_CT_TIME_FU_status == "uninflamed" & aggregated_results$min_CT_TIME_BR_status == "inflamed")
num_no_changes_br <- sum(aggregated_results$min_CT_TIME_FU_status == aggregated_results$min_CT_TIME_BR_status)

# Print the results
cat("Changes from inflamed to uninflamed:", num_changes_inflamed_to_uninflamed, "\n")
cat("Changes from uninflamed to inflamed:", num_changes_uninflamed_to_inflamed, "\n")
cat("No changes in status:", num_no_changes, "\n")


# Print the results
cat("Changes from inflamed to uninflamed BR:", num_changes_inflamed_to_uninflamed_br, "\n")
cat("Changes from uninflamed to inflamed BR:", num_changes_uninflamed_to_inflamed_br, "\n")
cat("No changes in status BR:", num_no_changes_br, "\n")

# Step 1
df2 <- aggregated_results %>%
  make_long(min_CT_TIME_BS_status, min_CT_TIME_FU_status, min_CT_TIME_BR_status)

# Step 2
dagg <- df2 %>%
  dplyr::group_by(node) %>%
  tally()

# Step 3
df3 <- merge(df2, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

# Set colors
cols <- c("inflamed" = "#b2182b", "uninflamed" = "#969696")

# Create Sankey plot
pl <- ggplot(df3, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = factor(node),
  label = ifelse(!is.na(n), paste0(node, " n=", n), NULL)
)) +
  geom_sankey(
    flow.alpha = 0.5,
    color = "gray40",
    show.legend = TRUE,
    aes(fill = factor(node))
  ) +
  geom_sankey_label(
    size = 3,
    color = "white",
    fill = "gray40",
    hjust = -0.2
  ) +
  scale_fill_manual(values = cols) +  # Set node colors
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank
  )
