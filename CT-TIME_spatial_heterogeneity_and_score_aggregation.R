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
response_path <- here::here("..","0_data","IMMUNOMIX_NEW","Total_Immuno_Clinical_Database_20230927.csv") #input data
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
#sub_resp <- subset(sub,clinical_benefit=="Yes")
#sub_prog <-subset(sub,clinical_benefit=="No")

# # lung tumors
# sub <- sub[grepl("ung", sub$primary_tumor_localization, ignore.case = TRUE), ]
#sub_resp <- subset(sub,clinical_benefit=="Yes")
#sub_prog <-subset(sub,clinical_benefit=="No")

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

# 
# # Create the plot (spatial variability)
# plot <- ggplot(sub, aes(x = reorder(ID, CT_TIME), y = CT_TIME)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(data = sub,
#               aes(x = ID, y = CT_TIME, colour = pred), size = 0.3
#               , width = 0.1) +
#   labs(x = "Patients", y = "CT-TIME tumour score",
#        title = "") +
#   theme_minimal() +
#   theme(legend.position = "top",  # Adjust legend position as needed
#         axis.text.x = element_text(size = 5)) +
#   scale_color_manual(
#     values = c("#b2182b", "#969696"),
#     breaks = c("inflamed","uninflamed"),
#     labels = c("Inflamed","Uninflamed")
#   )
# # Show plot
# print(plot)


# # Create the plot (spatial variability)
# 
# plot <- ggplot(sub, aes(x = reorder(ID, CT_TIME), y = CT_TIME, color = clinical_benefit)) +
#   geom_boxplot(outlier.shape = NA, aes(group = clinical_benefit), color = "black") +  # Color boxplot outlines by clinical_benefit
#   geom_jitter(data = sub,
#               aes(x = ID, y = CT_TIME, colour = lesion_location), size = 1,
#               width = 0.6) +
#   labs(x = "Patient", y = "CT-TIME Tumor Score",
#        title = "") +
#   theme_minimal() +
#   theme(legend.position = "top",  # Adjust legend position as needed
#         axis.text.x = element_text(size = 5)) +
#   scale_color_manual(
#     values = c("#fb9a99", "#67a9cf", "#b2df8a"),
#     breaks = c("Liver", "Lung", "Other"),
#     labels = c("Liver", "Lung", "Other")
#   )
# 
# # Show plot
# print(plot)
# 
# 
# # Calculate coefficient of variation for each patient
# result <- sub %>%
#   group_by(Patient) %>%
#   summarise(cv = sd(CT_TIME) / mean(CT_TIME) * 100)
# 
# 
# # Classify based on 35% threshold and count the number of patients above it
# threshold <- 35
# above_threshold <- result %>%
#   filter(cv > threshold)
# 
# # Count the number of patients above the threshold
# num_patients_above_threshold <- nrow(above_threshold)
# 
# num_patients_with_multiple_observations <- sub %>%
#   group_by(Patient) %>%
#   summarise(
#     num_observations = n()
#   ) %>%
#   filter(num_observations > 1) %>%
#   nrow()
# 
# # Percentage
# num_patients_above_threshold/num_patients_with_multiple_observations*100


###### Median CT-TIME
# Calculate CT-TIME scores per patient
my_data <- sub %>%
  group_by(ID) %>%
  summarize(mean_CT_TIME = mean(CT_TIME),
            min_CT_TIME = min(CT_TIME),
            max_CT_TIME = max(CT_TIME),
            median_CT_TIME = median(CT_TIME),
            sd_CT_TIME = sd(CT_TIME)
            )
            

my_data <- distinct(my_data)
sub_data <- subset(my_data,sd_CT_TIME>0)
# low variability
sum(sub_data$sd_CT_TIME < 0.1)
#moderate variability
sum(sub_data$sd_CT_TIME < 0.3 & sub_data$sd_CT_TIME > 0.1)
#high variability
sum(sub_data$sd_CT_TIME > 0.3)

my_data$mean_CT_TIME_status <- ifelse(my_data$mean_CT_TIME>0.5,"inflamed","uninflamed")
my_data$min_CT_TIME_status <- ifelse(my_data$min_CT_TIME>0.5,"inflamed","uninflamed")
my_data$max_CT_TIME_status <- ifelse(my_data$max_CT_TIME>0.5,"inflamed","uninflamed")
my_data$median_CT_TIME_status <- ifelse(my_data$median_CT_TIME>0.5,"inflamed","uninflamed")

# Define colors
inflamed_color <- "#b2182b"
uninflamed_color <- "#969696"


# Create the plot (histogram)
histogram_plot <- ggplot(my_data, aes(x = mean_CT_TIME, fill = mean_CT_TIME_status)) +
  geom_histogram( color = "black", breaks = seq(0, 1, 0.05)) +
  scale_fill_manual(values = c(inflamed_color,uninflamed_color)) +
  labs(x = "Mean of CT-TIME", y = "Number of Patients",
       title = "") +
  theme_minimal() +
  theme(legend.position = "top",  # Adjust legend position as needed
        axis.text.x = element_text(size = 5))

# Show plot
print(histogram_plot)

# Create the plot (histogram)
histogram_plot_var <- ggplot(my_data, aes(x = min_CT_TIME, fill = min_CT_TIME_status)) +
  geom_histogram(color = "black", breaks = seq(0, 1, 0.05)) +
  scale_fill_manual(values = c(inflamed_color,  uninflamed_color)) +
  labs(x = "Minimum of CT-TIME", y = "Number of Patients",
       title = "") +
  theme_minimal() +
  theme(legend.position = "top",  # Adjust legend position as needed
        axis.text.x = element_text(size = 5))

# Show plot
print(histogram_plot_var)

# Create the plot (histogram)
histogram_plot <- ggplot(my_data, aes(x = median_CT_TIME, fill = median_CT_TIME_status)) +
  geom_histogram( color = "black", breaks = seq(0, 1, 0.05)) +
  scale_fill_manual(values = c(inflamed_color,uninflamed_color)) +
  labs(x = "Median of CT-TIME", y = "Number of Patients",
       title = "") +
  theme_minimal() +
  theme(legend.position = "top",  # Adjust legend position as needed
        axis.text.x = element_text(size = 5))

# Show plot
print(histogram_plot)

# Create the plot (histogram)
histogram_plot_var <- ggplot(my_data, aes(x = max_CT_TIME, fill = max_CT_TIME_status)) +
  geom_histogram(color = "black", breaks = seq(0, 1, 0.05)) +
  scale_fill_manual(values = c(inflamed_color,  uninflamed_color)) +
  labs(x = "Maximum of CT-TIME", y = "Number of Patients",
       title = "") +
  theme_minimal() +
  theme(legend.position = "top",  # Adjust legend position as needed
        axis.text.x = element_text(size = 5))

# Show plot
print(histogram_plot_var)

# Create the plot (histogram)
histogram_plot <- ggplot(my_data, aes(x = mean_CT_TIME, fill = mean_CT_TIME_status)) +
  geom_histogram( color = "black", breaks = seq(0, 1, 0.05)) +
  scale_fill_manual(values = c(inflamed_color,uninflamed_color)) +
  labs(x = "Mean of CT-TIME", y = "Number of Patients",
       title = "") +
  theme_minimal() +
  theme(legend.position = "top",  # Adjust legend position as needed
        axis.text.x = element_text(size = 5))

# Show plot
print(histogram_plot)
# 
# # Create the plot (histogram)
# histogram_plot_var <- ggplot(my_data, aes(x = iqr_CT_TIME, fill = iqr_CT_TIME_status)) +
#   geom_histogram(color = "black", breaks = seq(0, 1, 0.05)) +
#   scale_fill_manual(values = c(inflamed_color,  uninflamed_color)) +
#   labs(x = "Interquartile range of CT-TIME", y = "Number of Patients",
#        title = "") +
#   theme_minimal() +
#   theme(legend.position = "top",  # Adjust legend position as needed
#         axis.text.x = element_text(size = 5))
# 
# # Show plot
# print(histogram_plot_var)



# Create a histogram for the distribution of CT-TIME scores for all lesion locations
histogram_plot_all_locations <- ggplot(sub, aes(x = CT_TIME, fill = lesion_location)) +
  geom_histogram(binwidth = 0.1, position = "stack", alpha = 0.7) +
  labs(x = "CT-TIME Tumor Score", y = "Number of Tumors") +
  theme_minimal() +
  theme(legend.position = "top",  # Adjust legend position as needed
        axis.text.x = element_text(size = 8))

# Show plot
print(histogram_plot_all_locations)

# Create a histogram for the distribution of CT-TIME scores for primary tumor
histogram_plot_all_locations <- ggplot(sub, aes(x = CT_TIME, fill = primary_tumor)) +
  geom_histogram(binwidth = 0.1, position = "stack", alpha = 0.7) +
  labs(x = "CT-TIME Tumor Score", y = "Number of Tumors") +
  theme_minimal() +
  theme(legend.position = "top",  # Adjust legend position as needed
        axis.text.x = element_text(size = 8))

# Show plot
print(histogram_plot_all_locations)


# ###### check dependence
# # Assuming 'your_data_frame' is your data frame
# result_anova <- aov(CT_TIME ~ lesion_location * primary_tumor_localization, data = sub)
# 
# # Print ANOVA summary
# summary(result_anova)

# library(multcomp)
# library(tidyr)
# library(broom)
# 
# # Create a model
# model <- aov(CT_TIME ~ lesion_location * primary_tumor_localization, data = sub)
# 
# # Perform Tukey's HSD post-hoc test
# tukey_results <- TukeyHSD(model)
# 
# # Summarize the results
# tukey_summary <- tidy(tukey_results)
# 
# # Print the summary
# print(tukey_summary)

########## CT-TIME aggregation 
aggregated_results <- sub %>%
  group_by(Patient) %>%
  summarise(
    min_CT_TIME = min(CT_TIME, na.rm = TRUE),
    max_CT_TIME = max(CT_TIME, na.rm = TRUE),
    median_CT_TIME = median(CT_TIME, na.rm = TRUE)
  )

aggregated_results <- merge(aggregated_results,response, by.x="Patient",by.y="SAP")
names(aggregated_results)

#
survival_data <- with(aggregated_results, Surv(PFS_months, pfs_cens))
cox_model_min <- coxph(Surv(PFS_months, pfs_cens) ~ min_CT_TIME, data = aggregated_results)
summary(cox_model_min)

cox_model_max <- coxph(Surv(PFS_months, pfs_cens) ~ max_CT_TIME, data = aggregated_results)
summary(cox_model_max)

cox_model_median <- coxph(Surv(PFS_months, pfs_cens) ~ median_CT_TIME, data = aggregated_results)
summary(cox_model_median)


# Create quartiles of the aggregated score
aggregated_results$CT_TIME_quartile <- cut(aggregated_results$min_CT_TIME, quantile(aggregated_results$min_CT_TIME, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
#aggregated_results$CT_TIME_dich <- cut(aggregated_results$min_CT_TIME, quantile(aggregated_results$min_CT_TIME, probs = c(0,  0.5, 1), na.rm = TRUE))
aggregated_results$CT_TIME_dich <- ifelse(aggregated_results$min_CT_TIME>=0.5, "Inflamed","Uninflamed")
# Fit a Kaplan-Meier curve by quartiles
#km_fit <- survfit(Surv(PFS_months, pfs_cens) ~ CT_TIME_quartile, data = aggregated_results)
km_fit <- survfit(Surv(PFS_months, pfs_cens) ~ CT_TIME_dich, data = aggregated_results)

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
           break.x.by = 5 ,
           xlim = c(0, 40))

lung <- subset(aggregated_results,primary_tumor_localization=="lung")

km_fit_lung <- survfit(Surv(PFS_months, pfs_cens) ~ CT_TIME_dich, data = lung)

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
aggregated_results$primary_tumor <- ifelse(
  aggregated_results$primary_tumor_localization %in% c("adrenal", "thymus", "thyroid"), "endocrine",
  ifelse(
    aggregated_results$primary_tumor_localization %in% c("pancreas", "liver", "biliary_tract", "cholangiocarcinoma", "hepatocarcinoma"), "hepatobiliary",
    ifelse(
      aggregated_results$primary_tumor_localization %in% c("colon", "esophagus", "gastric", "rectal", "colorectal", "UGE"), "gastrointestinal",
      ifelse(
        aggregated_results$primary_tumor_localization %in% c("kidney", "bladder", "prostate", "penis", "renal"), "genitourinary",
        ifelse(
          aggregated_results$primary_tumor_localization %in% c("cervix", "ovarian", "endometrial", "endometrium", "ovary", "ureteral", "urothelial"), "gynecologic",
          ifelse(
            aggregated_results$primary_tumor_localization %in% c("melanoma", "skin_uveal"), "skin",
            ifelse(
              aggregated_results$primary_tumor_localization %in% c("retroperitoneal", "bone", "sacral chordoma"), "sarcoma",
              ifelse(
                aggregated_results$primary_tumor_localization %in% c("breast TN", "breast/ovary"), "breast",
                ifelse(
                  aggregated_results$primary_tumor_localization %in% c("head_and_neck", "H&N"), "head&neck",
                  ifelse(
                    aggregated_results$primary_tumor_localization %in% c("mesothelioma", "pleura"), "lung", NA
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
###### Univariate Analysis (target: PFS) ######
time_variable <- "PFS_months"
event_variable <- "pfs_cens"

#making formulas
univ_formulas <- sapply(names(aggregated_results[c(2:4)]),function(x)as.formula(paste("Surv(", time_variable, ",", event_variable, ") ~",x)))
#making a list of models
univ_models <- lapply(univ_formulas, function(x){coxph(x,data=aggregated_results)})
#extract data (here I've gone for HR and confint)
univ_results <- lapply(univ_models,function(x){return(exp(cbind(coef(x),confint(x))))})


df_univ_combined <- do.call(rbind, lapply(univ_results, function(result) {
  result_df <- as.data.frame(result)
  return(result_df)
}))

# Move row names to the first column
df_univ_combined <- rownames_to_column(df_univ_combined, var = "predictor")
df_univ_combined$model <- "Univariate"
# Rename the columns
colnames(df_univ_combined) <- c("predictor","estimate", "conf.low", "conf.high","cox_model")


######### UVA  plot ######
dotCOLS <- c("#bababa", "#4d4d4d")
barCOLS <- c("#bababa", "#4d4d4d")


ggplot(df_univ_combined, aes(x = estimate, xmax = conf.high, xmin = conf.low, y = predictor,color = cox_model,fill = cox_model)) +
  geom_point(size = 3, shape = 18, position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5), size = 1) +
  geom_vline(xintercept = 1, size = 1) +
  #geom_hline(yintercept = 0, size = 1) +
  scale_alpha_identity() +
  scale_fill_manual(values = barCOLS) +
  scale_color_manual(values = dotCOLS) +
  scale_x_continuous(name = "Hazard ratio", limits = c(0, 1.5)) +
  coord_cartesian(clip = "off") +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(0, "pt"),
    axis.line.x.bottom = element_line(size = 1),
    axis.text.y.left =  element_text(margin = margin(l = 20, unit = "pt")),
    strip.background = element_blank(),
    strip.text = element_blank()
  )
