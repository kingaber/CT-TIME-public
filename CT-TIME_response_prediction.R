###########################################################
#                                                         #
#     CT-TIME response prediction              #
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
###### Define paths and load data ######
data_path <- here::here("..","1_data_processed", "raw_NANO_IMMUNOMIX_data.Rda") #input data
response_path <- here::here("..","0_data","IMMUNOMIX_NEW","Total_Immuno_Clinical_Database_20231211.csv") #input data
signature_path <- here::here("..","2_training","CT-TIME_signature.Rda") #input data
model_path <- here::here( "..","2_training","CT-TIME_model.rds") #input model
load(data_path)
response <- read.csv2(response_path)
load(signature_path)
glmnet_fit <- readRDS(model_path)
###### Scale data and get outcome ######
### split cohorts
training <- data_all[data_all$DataBatch=="VHIO",]
#scale features
training[c(4:110)] <- scale(training[c(4:110)])

validation <- data_all[data_all$DataBatch=="IMMUNOMIX_Soft",]
#scale features
validation[c(4:110)] <- scale(validation[c(4:110)])

############ Feature filtering ############
training_features <- training[,c(selected_feature_names,"Patient","LesionID")]
training_id <- data.frame(training_features)

validation_features <- validation[,c(selected_feature_names,"Patient","LesionID")]
validation_id <- data.frame(validation_features)

######### 
training_id$CT_TIME <- predict(glmnet_fit,training_id ,type = "prob")[,1]
training_id$pred  <- predict(glmnet_fit,training_id ,type = "raw")

validation_id$CT_TIME <- predict(glmnet_fit,validation_id ,type = "prob")[,1]
validation_id$pred  <- predict(glmnet_fit,validation_id ,type = "raw")

######## Merge with Response Data #####
response <- subset(response, PFS_months > 0)
response$clinical_benefit <- ifelse(response$PFS_months>=5, "Yes","No")
response$clinical_benefit <- as.factor(response$clinical_benefit)

df <- merge(validation_id,response,by.x="Patient",by.y="SAP", all.y=TRUE)
# Remove rows with NA in the specific column
df <- df[complete.cases(df[, "original_shape_Elongation"]), ]

df_train <- merge(training_id,response,by.x="Patient",by.y="SAP", all.y=TRUE)
# Remove rows with NA in the specific column
df_train <- df_train[complete.cases(df_train[, "original_shape_Elongation"]), ]



# ######## Spatial heterogeneity
sub <- df[c("Patient","LesionID","CT_TIME","primary_tumor_localization", "pred","clinical_benefit","PFS_months","pfs_cens")]
sub <- sub[order(sub$CT_TIME),]
sub <- transform(sub,
                 ID = paste0(as.numeric(factor(Patient))))

####### cleanup locations #########
# Create the new column based on the dictionary
sub$primary_tumor <- ifelse(
  sub$primary_tumor_localization %in% c("adrenal", "thymus", "thyroid","paraganglioma retroperitoneal"), "endocrine",
  ifelse(
    sub$primary_tumor_localization %in% c("pancreas", "liver", "biliary_tract", "cholangiocarcinoma", "hepatocarcinoma"), "hepatobiliary",
    ifelse(
      sub$primary_tumor_localization %in% c("colon", "esophagus", "gastric", "rectal", "colorectal", "UGE","stomach"), "gastrointestinal",
      ifelse(
        sub$primary_tumor_localization %in% c("kidney", "bladder", "prostate", "penis", "renal"), "genitourinary",
        ifelse(
          sub$primary_tumor_localization %in% c("cervix", "ovarian", "endometrial", "endometrium", "ovary", "ureteral", "urothelial"), "gynecologic",
          ifelse(
            sub$primary_tumor_localization %in% c("melanoma", "skin_uveal", "skin"), "skin",
            ifelse(
              sub$primary_tumor_localization %in% c("retroperitoneal", "bone", "sacral chordoma"), "sarcoma",
              ifelse(
                sub$primary_tumor_localization %in% c("breast TN", "breast/ovary", "breast"), "breast",
                ifelse(
                  sub$primary_tumor_localization %in% c("head_and_neck", "H&N"), "head&neck",
                  ifelse(
                    sub$primary_tumor_localization %in% c("mesothelioma", "pleura","lung"), "lung", NA
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
# sub_resp <- subset(sub,clinical_benefit=="Yes")
# sub_prog <-subset(sub,clinical_benefit=="No")

# # lung location
# sub <- sub[grepl("LU|ung", sub$LesionID, ignore.case = TRUE), ]
# sub_resp <- subset(sub,clinical_benefit=="Yes")
# sub_prog <-subset(sub,clinical_benefit=="No")

# # lung tumors
# sub <- sub[grepl("ung", sub$primary_tumor_localization, ignore.case = TRUE), ]
# sub_resp <- subset(sub,clinical_benefit=="Yes")
# sub_prog <-subset(sub,clinical_benefit=="No")

# get lesion location
sub <- sub %>%
  mutate(
    lesion_location = case_when(
      grepl("^liver_", LesionID) ~ "liver",
      grepl("^LU_", LesionID) | grepl("^lung_", LesionID) ~ "lung",
      TRUE ~ "other"
    )
  )

unique(sub$Patient)
sub <- distinct(sub)

########## CT-TIME aggregation 
aggregated_results <- sub %>%
  group_by(Patient) %>%
  summarise(
    min_CT_TIME = min(CT_TIME, na.rm = TRUE),
    max_CT_TIME = max(CT_TIME, na.rm = TRUE),
    mean_CT_TIME = mean(CT_TIME, na.rm = TRUE),
    median_CT_TIME = median(CT_TIME, na.rm = TRUE)
  )

aggregated_results$mean_CT_TIME_status <- ifelse(aggregated_results$mean_CT_TIME>0.5,"inflamed","uninflamed")
aggregated_results$min_CT_TIME_status <- ifelse(aggregated_results$min_CT_TIME>0.5,"inflamed","uninflamed")
aggregated_results$max_CT_TIME_status <- ifelse(aggregated_results$max_CT_TIME>0.5,"inflamed","uninflamed")
aggregated_results$median_CT_TIME_status <- ifelse(aggregated_results$median_CT_TIME>0.5,"inflamed","uninflamed")


aggregated_results <- merge(aggregated_results,response, by.x="Patient",by.y="SAP")
names(aggregated_results)
aggregated_results <- distinct(aggregated_results)
unique(aggregated_results$Patient)
### some patients have 2 rows in response; I have only one image, so need to remove:
subset_duplicated <- subset(aggregated_results ,duplicated=="yes")
# subset_duplicated[c("Patient","min_CT_TIME", "trial_best_response","PFS_months")]
# rows_to_remove <- c(6,44,79,128,173,183,209,326)
aggregated_results <- aggregated_results[-c(6,44,79,128,173,183,209,326),]
#
aggregated_results$pfs_cens <- as.numeric(aggregated_results$pfs_cens)

cox_model_min <- coxph(Surv(PFS_months, pfs_cens) ~ min_CT_TIME, data = aggregated_results)
summary(cox_model_min)

cox_model_max <- coxph(Surv(PFS_months, pfs_cens) ~ max_CT_TIME, data = aggregated_results)
summary(cox_model_max)

cox_model_mean <- coxph(Surv(PFS_months, pfs_cens) ~ mean_CT_TIME, data = aggregated_results)
summary(cox_model_mean)

cox_model_median <- coxph(Surv(PFS_months, pfs_cens) ~ median_CT_TIME, data = aggregated_results)
summary(cox_model_median)

#### forest plot
# Extract hazard ratios and confidence intervals
hr_min <- exp(coef(cox_model_min))
ci_min <- exp(confint(cox_model_min))

hr_max <- exp(coef(cox_model_max))
ci_max <- exp(confint(cox_model_max))

hr_mean <- exp(coef(cox_model_mean))
ci_mean <- exp(confint(cox_model_mean))


hr_median <- exp(coef(cox_model_median))
ci_median <- exp(confint(cox_model_median))

# Create a data frame with the results
forest_data <- data.frame(
  Variable = c("min_CT_TIME", "max_CT_TIME", "mean_CT_TIME","median_CT_TIME"),
  HR = c(hr_min, hr_max, hr_mean,hr_median),
  Low = c(ci_min[1, 1], ci_max[1, 1], ci_mean[1, 1],ci_median[1, 1]),
  High = c(ci_min[1, 2], ci_max[1, 2], ci_mean[1, 2],ci_median[1, 2])
)

forest_plot <- ggplot(forest_data, aes(x = HR, y = Variable)) +
  geom_point(aes(color = Variable), size = 3) +
  geom_errorbarh(aes(xmin = Low, xmax = High), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.5, color = "black") +  # Add a thick line at HR = 1
  labs(x = "Hazard Ratio", y = "", title = "Forest Plot") +
  scale_x_continuous(limits = c(0, 3)) +  # Set x-axis limits
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend if not needed
  annotate("text", x = 2, y = nrow(forest_data) + 0.5, label = "Cox-PH for PFS", size = 3, color = "black")

# Show the plot
print(forest_plot)
# 
# # Create quartiles of the aggregated score
# aggregated_results$CT_TIME_quartile <- cut(aggregated_results$min_CT_TIME, quantile(aggregated_results$min_CT_TIME, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
# #aggregated_results$CT_TIME_dich <- cut(aggregated_results$min_CT_TIME, quantile(aggregated_results$min_CT_TIME, probs = c(0,  0.5, 1), na.rm = TRUE))
library(cutpointr)
#find optimal cutoff point
cp <- cutpointr( aggregated_results,  median_CT_TIME, clinical_benefit,
                 method = maximize_metric, metric = sum_sens_spec)
cp
#plot(cp)

# min_cp <- 0.656094
# max_cp <- 0.977487
# median_cp <- 0.529797

min_cp <- 0.66
max_cp <- 0.98
mean_cp <- 0.5

aggregated_results$min_CT_TIME_status <- ifelse(aggregated_results$min_CT_TIME>=min_cp, "inflamed","uninflamed")
aggregated_results$mean_CT_TIME_status <- ifelse(aggregated_results$mean_CT_TIME>=mean_cp, "inflamed","uninflamed")
aggregated_results$max_CT_TIME_status <- ifelse(aggregated_results$max_CT_TIME>=max_cp, "inflamed","uninflamed")
aggregated_results$median_CT_TIME_status <- ifelse(aggregated_results$median_CT_TIME>=mean_cp, "inflamed","uninflamed")
aggregated_results$prog_as_br <- ifelse(aggregated_results$trial_best_response=="PD","yes","no")
aggregated_results$or_as_br <- ifelse(aggregated_results$trial_best_response=="CR"|aggregated_results$trial_best_response=="PR","yes","no")

# Fit a Kaplan-Meier curve by quartiles
km_fit <- survfit(Surv(PFS_months, pfs_cens) ~ median_CT_TIME_status, data = aggregated_results)

# Plot the Kaplan-Meier curves
ggsurvplot(km_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c( "#b2182b", "#969696"),
           xlab = 'Progression free survival (months)',
           break.x.by = 5,
           xlim = c(0, 40) )


lung <- subset(aggregated_results,primary_tumor_localization=="lung")

km_fit_lung <- survfit(Surv(PFS_months, pfs_cens) ~ median_CT_TIME_status, data = lung)

# Plot the Kaplan-Meier curves
ggsurvplot(km_fit_lung,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c( "#b2182b", "#969696"),
           xlab = 'Progression free survival (months)',
           break.x.by = 5 ,
           xlim = c(0, 40))


aggregated_results$sex <- gsub("female ", "female", aggregated_results$sex)

table(aggregated_results$IO_combo)
aggregated_results$IO_combo <- gsub("combo ", "combo", aggregated_results$IO_combo)
aggregated_results$IO_combo <- gsub("comboimmuno", "combo_immuno", aggregated_results$IO_combo)
aggregated_results$IO_combo <- gsub("combo_immuno ", "combo_immuno", aggregated_results$IO_combo)

library(lubridate)
# Convert DOB and start_treatment_date to Date objects
betterDates <- as.Date(aggregated_results$DOB, "%d/%m/%y")
correctCentury <- as.Date(ifelse(betterDates > Sys.Date(), 
                                 format(betterDates, "19%y-%m-%d"), 
                                 format(betterDates)))

aggregated_results$age <- round(time_length(difftime(as.Date(aggregated_results$trial_start_date, format = "%d/%m/%y"),correctCentury), "years") ) # Calculate difference in years

# Create the new column based on the dictionary
sub$primary_tumor <- ifelse(
  sub$primary_tumor_localization %in% c("adrenal", "thymus", "thyroid","paraganglioma retroperitoneal"), "endocrine",
  ifelse(
    sub$primary_tumor_localization %in% c("pancreas", "liver", "biliary_tract", "cholangiocarcinoma", "hepatocarcinoma"), "hepatobiliary",
    ifelse(
      sub$primary_tumor_localization %in% c("colon", "esophagus", "gastric", "rectal", "colorectal", "UGE","stomach"), "gastrointestinal",
      ifelse(
        sub$primary_tumor_localization %in% c("kidney", "bladder", "prostate", "penis", "renal"), "genitourinary",
        ifelse(
          sub$primary_tumor_localization %in% c("cervix", "ovarian", "endometrial", "endometrium", "ovary", "ureteral", "urothelial"), "gynecologic",
          ifelse(
            sub$primary_tumor_localization %in% c("melanoma", "skin_uveal", "skin"), "skin",
            ifelse(
              sub$primary_tumor_localization %in% c("retroperitoneal", "bone", "sacral chordoma"), "sarcoma",
              ifelse(
                sub$primary_tumor_localization %in% c("breast TN", "breast/ovary", "breast"), "breast",
                ifelse(
                  sub$primary_tumor_localization %in% c("head_and_neck", "H&N"), "head&neck",
                  ifelse(
                    sub$primary_tumor_localization %in% c("mesothelioma", "pleura","lung"), "lung", NA
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

results_cont <- data.frame(Original = aggregated_results$clinical_benefit, Predicted = aggregated_results$min_CT_TIME)

# Create a boxplot
ggplot(results_cont, aes(x = Original, y = Predicted, fill = Original)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_fill_manual(values = c("#e74c3c","#3498db") )+
  labs(title = "Cohort 4: pan-cancer patients",
       x = "Clinical benefit at 5 months",
       y = "CT-TIME (Aggregated Score)")

results_cont <- data.frame(Original = lung$clinical_benefit, Predicted = lung$min_CT_TIME)

# Create a boxplot
ggplot(results_cont, aes(x = Original, y = Predicted, fill = Original)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_fill_manual(values = c("#e74c3c","#3498db") )+
  labs(title = "Cohort 4: lung cancer patients",
       x = "Clinical benefit at 5 months",
       y = "CT-TIME (Aggregated Score)")

library(pROC)
roc_min <- roc(response = aggregated_results$clinical_benefit, predictor = aggregated_results$min_CT_TIME, ci = TRUE, direction = "<")
roc_max <- roc(response = aggregated_results$clinical_benefit, predictor = aggregated_results$max_CT_TIME, ci = TRUE, direction = "<")
roc_mean <- roc(response = aggregated_results$clinical_benefit, predictor = aggregated_results$mean_CT_TIME, ci = TRUE, direction = "<")
roc_median <- roc(response = aggregated_results$clinical_benefit, predictor = aggregated_results$median_CT_TIME, ci = TRUE, direction = "<")

# Plot ROC curves
par(mfrow = c(1, 1))
plot(roc_min, col = "#01665E", lwd = 2, lty = 1, main = 'ROC: Clinical Benefit at 5 Months', xlab = "1 - Specificity", ylab = "Sensitivity")
lines(roc_max, col = "#01665E", lwd = 2, lty = 2)
lines(roc_mean, col = "#01665E", lwd = 1, lty = 1)
lines(roc_median, col = "#01665E", lwd = 1, lty = 2)

# Add legend
legend('bottomright', legend = c(
  paste('min_CT_TIME: ', format(round(roc_min$auc, 2), nsmall = 2), ' (', format(round(roc_min$ci[1], 2), nsmall = 2), '-', format(round(roc_min$ci[3], 2), nsmall = 2), ')', sep = ''),
  paste('max_CT_TIME: ', format(round(roc_max$auc, 2), nsmall = 2), ' (', format(round(roc_max$ci[1], 2), nsmall = 2), '-', format(round(roc_max$ci[3], 2), nsmall = 2), ')', sep = ''),
  paste('mean_CT_TIME: ', format(round(roc_mean$auc, 2), nsmall = 2), ' (', format(round(roc_mean$ci[1], 2), nsmall = 2), '-', format(round(roc_mean$ci[3], 2), nsmall = 2), ')', sep = ''),
  paste('median_CT_TIME: ', format(round(roc_median$auc, 2), nsmall = 2), ' (', format(round(roc_median$ci[1], 2), nsmall = 2), '-', format(round(roc_median$ci[3], 2), nsmall = 2), ')', sep = '')
  ), lty = c(1, 2, 1, 2), lwd = c(2, 2, 1, 1), col = c("#01665E", "#01665E", "#01665E", "#01665E"), cex = 0.6, bg = "white")

# Plot ROC curves for LUNG
lung <- subset(aggregated_results,primary_tumor_localization=="lung")

roc_min <- roc(response = lung$clinical_benefit, predictor = lung$min_CT_TIME, ci = TRUE, direction = "<")
roc_max <- roc(response = lung$clinical_benefit, predictor = lung$max_CT_TIME, ci = TRUE, direction = "<")
roc_mean <- roc(response = lung$clinical_benefit, predictor = lung$mean_CT_TIME, ci = TRUE, direction = "<")
roc_median <- roc(response = lung$clinical_benefit, predictor = lung$median_CT_TIME, ci = TRUE, direction = "<")


par(mfrow = c(1, 1))
plot(roc_min, col = "#01665E", lwd = 2, lty = 1, main = 'ROC: Clinical Benefit at 5 Months', xlab = "1 - Specificity", ylab = "Sensitivity")
lines(roc_max, col = "#01665E", lwd = 2, lty = 2)
lines(roc_mean, col = "#01665E", lwd = 1, lty = 1)
lines(roc_median, col = "#01665E", lwd = 1, lty = 2)

# Add legend
legend('bottomright', legend = c(
  paste('min_CT_TIME: ', format(round(roc_min$auc, 2), nsmall = 2), ' (', format(round(roc_min$ci[1], 2), nsmall = 2), '-', format(round(roc_min$ci[3], 2), nsmall = 2), ')', sep = ''),
  paste('max_CT_TIME: ', format(round(roc_max$auc, 2), nsmall = 2), ' (', format(round(roc_max$ci[1], 2), nsmall = 2), '-', format(round(roc_max$ci[3], 2), nsmall = 2), ')', sep = ''),
  paste('mean_CT_TIME: ', format(round(roc_mean$auc, 2), nsmall = 2), ' (', format(round(roc_mean$ci[1], 2), nsmall = 2), '-', format(round(roc_mean$ci[3], 2), nsmall = 2), ')', sep = ''),
  paste('median_CT_TIME: ', format(round(roc_median$auc, 2), nsmall = 2), ' (', format(round(roc_median$ci[1], 2), nsmall = 2), '-', format(round(roc_median$ci[3], 2), nsmall = 2), ')', sep = '')
  
  ), lty = c(1, 2, 1,2), lwd = c(2, 2, 1,1), col = c("#01665E","#01665E", "#01665E", "#01665E"), cex = 0.6, bg = "white")


cox_model_min <- coxph(Surv(PFS_months, pfs_cens) ~ min_CT_TIME, data = lung)
summary(cox_model_min)

cox_model_max <- coxph(Surv(PFS_months, pfs_cens) ~ max_CT_TIME, data = lung)
summary(cox_model_max)

cox_model_mean <- coxph(Surv(PFS_months, pfs_cens) ~ mean_CT_TIME, data = lung)
summary(cox_model_mean)

cox_model_median <- coxph(Surv(PFS_months, pfs_cens) ~ median_CT_TIME, data = lung)
summary(cox_model_median)

#### forest plot
# Extract hazard ratios and confidence intervals
hr_min <- exp(coef(cox_model_min))
ci_min <- exp(confint(cox_model_min))

hr_max <- exp(coef(cox_model_max))
ci_max <- exp(confint(cox_model_max))

hr_mean <- exp(coef(cox_model_mean))
ci_mean <- exp(confint(cox_model_mean))

hr_median <- exp(coef(cox_model_median))
ci_median <- exp(confint(cox_model_median))
# Create a data frame with the results
forest_data <- data.frame(
  Variable = c("min_CT_TIME", "max_CT_TIME", "mean_CT_TIME", "median_CT_TIME"),
  HR = c(hr_min, hr_max,hr_mean ,hr_median),
  Low = c(ci_min[1, 1], ci_max[1, 1], ci_mean[1, 1],ci_median[1, 1]),
  High = c(ci_min[1, 2], ci_max[1, 2], ci_mean[1, 2],ci_median[1,2])
)

forest_plot <- ggplot(forest_data, aes(x = HR, y = Variable)) +
  geom_point(aes(color = Variable), size = 3) +
  geom_errorbarh(aes(xmin = Low, xmax = High), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.5, color = "black") +  # Add a thick line at HR = 1
  labs(x = "Hazard Ratio", y = "", title = "Forest Plot") +
  scale_x_continuous(limits = c(0, 3)) +  # Set x-axis limits
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend if not needed
  annotate("text", x = 2, y = nrow(forest_data) + 0.5, label = "Cox-PH for PFS", size = 3, color = "black")

# Show the plot
print(forest_plot)

# # Create ROC curves pan cancer
# R_val <- pROC::roc(response = aggregated_results$clinical_benefit, predictor = aggregated_results$min_CT_TIME, ci = TRUE, direction="<", 
#                    xlab = "False Positive Rate", ylab = "True Positive Rate", 
#                    main = 'CT Radiomics Signature \n Model Logistic Regression', 
#                    cex.main = 1, col = "#01665E")
# plot(R_val, legacy.axes = TRUE, col ="#01665E", lwd = 3,lty = 1)
# legend('bottomright', legend = c(paste('Cohort 4: ',format(round(R_val$auc,2),nsmall = 2),
#                                        ' (',format(round(R_val$ci[1],2),nsmall = 2),'-',format(round(R_val$ci[3],2),nsmall = 2),')',sep ='')),
#        lty = c(1, 1), lwd = c(2.5, 2.5), col = c("#01665E", "#01665E"), cex = 0.6, bg = "white")
# 
# R_val$auc;R_val$ci
# # Area under the curve: 0.569
# # 95% CI: 0.5015-0.6366 (DeLong)
# 
# # Create ROC curves lung
# R_val <- pROC::roc(response = lung$clinical_benefit, predictor = lung$min_CT_TIME, ci = TRUE, direction="<", 
#                    xlab = "False Positive Rate", ylab = "True Positive Rate", 
#                    main = 'CT Radiomics Signature \n Model Logistic Regression - Training and Validation Set', 
#                    cex.main = 1, col = "#01665E")
# plot(R_val, legacy.axes = TRUE, col ="#01665E", lwd = 3,lty = 1)
# legend('bottomright', legend = c(paste('Cohort 4 (lung): ',format(round(R_val$auc,2),nsmall = 2),
#                                        ' (',format(round(R_val$ci[1],2),nsmall = 2),'-',format(round(R_val$ci[3],2),nsmall = 2),')',sep ='')),
#        lty = c(1, 1), lwd = c(2.5, 2.5), col = c("#01665E", "#01665E"), cex = 0.6, bg = "white")
# 
# 
# R_val$auc;R_val$ci
# # Area under the curve: 0.6875
# # 95% CI: 0.552-0.823 (DeLong)

########
###### finding cutoff
# library(cutpointr)
# #find optimal cutoff point
# cp <- cutpointr( aggregated_results,  min_CT_TIME, clinical_benefit, 
#                  method = maximize_metric, metric = sum_sens_spec)
# plot(cp)
# predict(cp, newdata = data.frame(training_selected_features))


aggregated_results$CT_TIME_dichot <- ifelse(aggregated_results$min_CT_TIME>=0.5, "Yes","No")
lung$CT_TIME_dichot <- ifelse(lung$min_CT_TIME>=0.5, "Yes","No")

# Calculate the confusion matrix
conf_matrix <- confusionMatrix(as.factor(aggregated_results$CT_TIME_dichot), reference = as.factor(aggregated_results$clinical_benefit))
conf_matrix

# Calculate the confusion matrix
conf_matrix <- confusionMatrix(as.factor(lung$CT_TIME_dichot), reference = as.factor(lung$clinical_benefit))
conf_matrix



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
