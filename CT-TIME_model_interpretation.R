###########################################################
#                                                         #
#     CT-TIME model interpretation             #
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
library(corrplot)
###### Define paths and load data ######
data_path <- here::here("..","1_data_processed", "raw_NANO_IMMUNOMIX_data.Rda") #input data
ihc_path <- here::here("..","0_data","PREDICT","IHC_imp_2023-07-10_15_36.Rda") #input data
signature_path <- here::here("..","2_training","CT-TIME_signature.Rda") #input data
model_path <- here::here( "..","2_training","CT-TIME_model.rds") #input model
load(data_path)
load(ihc_path)
load(signature_path)
glmnet_fit <- readRDS(model_path)
###### Scale data and get outcome ######
### split cohorts
validation <- data_all[data_all$DataKernel=="Soft_predict",]
#scale features
validation[c(4:110)] <- scale(validation[c(4:110)])

############ Feature filtering ############

validation_features <- validation[,c(selected_feature_names,"Patient","LesionID")]
validation_id <- data.frame(validation_features)

######### 

validation_id$CT_TIME <- predict(glmnet_fit,validation_id ,type = "prob")[,1]
validation_id$pred  <- predict(glmnet_fit,validation_id ,type = "raw")

# Define a regular expression pattern for digits
pattern <- "\\d+"
# Use grep to find the numeric part in each string
validation_id$PatientNumber <- regmatches(validation_id$Patient, regexpr(pattern,validation_id$Patient))
ihc_imp$PatientNumber <- regmatches(ihc_imp$PredictID, regexpr(pattern,ihc_imp$PredictID))
ihc_imp <- subset(ihc_imp,Timepoint=="BL")
######## Merge with IHC #####
df <- merge(validation_id,ihc_imp,by="PatientNumber", all.y=TRUE)
# Remove rows with NA in the specific column
df <- df[complete.cases(df[, "original_shape_Elongation"]), ]

# Keep only biopsied lesions "_bio"
df_sub <- df[grepl("_bio", df[, "LesionID"]), ]

p <- ggboxplot(df_sub, x = "pred", y = "i_CD8",
               color = "pred", palette = "jco",add='jitter')
p + stat_compare_means()


# Scatter plot with correlation coefficient
#:::::::::::::::::::::::::::::::::::::::::::::::::
sp <- ggscatter(df_sub, x = "CT_TIME", y = "i_CD8",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#35978f", fill = "#35978f"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = , label.y = 8)
# icd163
sp <- ggscatter(df_sub, x = "CT_TIME", y = "i_CD163",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#969696", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = , label.y = 8)
# i cd3 tot
sp <- ggscatter(df_sub, x = "CT_TIME", y = "i_CD3_tot",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#969696", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = , label.y = 8)

# i cd3 tot
sp <- ggscatter(df_sub, x = "CT_TIME", y = "PDL1_T",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#969696", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = , label.y = 8)

######## Correlation of CT-TIME with IHC markers
cor_df <- df_sub[c("CT_TIME","Ki67","PDL1_T","CD31_D","i_CD163",
                     "i_CD8","i_FOX","i_CD3_tot")]

CorrMatrix <- cor(cor_df, method= "pearson")
testRes = cor.mtest(cor_df, conf.level = 0.95)
corrplot(CorrMatrix, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20')

######## Correlation of features with IHC markers
colnames(df_sub)[1:5] <- c("PatientNumber","shape_Elongation","firstorder_90Perc","glcm_MaxProb"  ,"glszm_SZNU")

cor_df <- df_sub[c("shape_Elongation","firstorder_90Perc","glcm_MaxProb"  ,"glszm_SZNU",
                   "Ki67","PDL1_T","CD31_D","i_CD163",
                   "i_CD8","i_FOX","i_CD3_tot")]

CorrMatrix <- cor(cor_df, method= "pearson")
testRes = cor.mtest(cor_df, conf.level = 0.95)
corrplot(CorrMatrix, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20')

