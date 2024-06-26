###########################################################
#                                                         #
#     CT-TIME signature                #
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
data_path <- here::here("..","1_data_processed", "harmonized_NANO_TCIA_data.Rda") #input data
robust_path <- here::here("..","0_data","NANOSTRING", "robustness.csv") #input data
signature_path <- here::here( "CT-TIME_signature.Rda") #output signature
load(data_path)
robust <- read.csv2(robust_path)
###### Scale data and get outcome ######
training <- data_cb[data_cb$DataBatch=="VHIO",]
#scale features
training[c(6:114)] <- scale(training[c(6:114)])

############ Feature filtering ############
robust_feature_names <- subset(robust,Concordance.Correlation.Coefficient>=0.75)$feature
robust_features <- training[,c(robust_feature_names)]
training_df <- data.frame(robust_features)
med_training <- median(training$GEP_score)
training_df$outcome <- ifelse(training$GEP_score>=med_training,"inflamed","uninflamed")
############ ElasticNet: hyperparameters tuning ############
# Step 1: Get alpha parameter using repeatedcv
set.seed(1234)
cv_5 = trainControl(method = "repeatedcv", number = 5,  repeats=10, returnResamp="all",
                    classProbs=TRUE, summaryFunction=twoClassSummary)

nano_elnet = train(
  outcome ~ ., data = training_df,
  method = "glmnet",metric = "ROC",
  trControl = cv_5
)

bestalpha = nano_elnet$bestTune$alpha #bestalpha 0.55

# Step 2: Get lambda 
set.seed(1234)
cv_5 = trainControl(method = "repeatedcv", number = 5, repeats=10, returnResamp="all",
                    classProbs=TRUE, summaryFunction=twoClassSummary)

nano_elnet = train(
  outcome ~ ., data = training_df,
  method = "glmnet",
  trControl = cv_5, metric = "ROC",
  tuneGrid = expand.grid(alpha = bestalpha,
                         lambda = seq(0.001,0.1,by = 0.001))
  
)

bestlambda = nano_elnet$bestTune$lambda #bestlambda 0.007


###### Elastic Net: feature selection (repeatedcv) ######
set.seed(1234)
cv_5 = trainControl(method = "repeatedcv", number = 5, repeats=10, returnResamp="all",
                    classProbs=TRUE, summaryFunction=twoClassSummary)

nano_fit = train(
  outcome ~ ., data = training_df,
  method = "glmnet",
  trControl = cv_5, metric = "ROC",
  tuneGrid = expand.grid(alpha = bestalpha,
                         lambda = bestlambda))

myCoefs1 <- coef(nano_fit$finalModel, nano_fit$bestTune$lambda) 
# Access the sparse matrix as a regular matrix
myCoefs1_dense <- as.matrix(myCoefs1)

# Extract the non-zero coefficients
non_zero_indices <- myCoefs1_dense[, 1] != 0
preselected_features_names  <- myCoefs1@Dimnames[[1]][non_zero_indices][-1] 

# Print the selected features
print(preselected_features_names)

preselected_features <- training_df[,c(preselected_features_names)]
training_df_2 <- data.frame(preselected_features)
heatmap(as.matrix(training_df_2))

# Remove redundant features
correlationMatrix <- cor(training_df_2)
redundant_feature_names <- findCorrelation(correlationMatrix, cutoff = .59,  names = TRUE, exact=TRUE) 
selected_features <- training_df_2[,-which(colnames(training_df_2) %in% redundant_feature_names)]
selected_feature_names <- colnames(selected_features) #
selected_feature_names 

training_selected_features <- data.frame(selected_features)
training_selected_features$outcome<- training_df$outcome
training_selected_features$outcome01 <- ifelse(training_selected_features$outcome=="inflamed",1,0)
training_selected_features$outcomeCONT <- training$GEP_score

#[1] "original_shape_Elongation"            "original_firstorder_90Percentile"    
#[3] "original_glcm_MaximumProbability"     "original_glszm_SizeZoneNonUniformity"

# GEP training signature 
# (1) pre-filter robust: 34 (robust features)
# (2) selected by Elastic Net: 13 (preselected features)
# (3) non-redundant: 4 (selected features)

save(selected_feature_names,file=signature_path)



# ###### Univariate Analysis (target: cont. GEP) ######
# 
# outcome_variable <- "outcomeCont"
# selected_formula <- as.formula(paste(outcome_variable, "~", paste(selected_feature_names , collapse = "+")))
# 
# tbl_uvregression(  training_selected_features[-c(5,6)] ,                       
#                    method = lm,
#                    y = outcomeCont,                               
#                    exponentiate = FALSE                     ## 
# )
# 
# ###### Multivariate Analysis ######
# mod2 <- lm(selected_formula, data= training_selected_features )
# tbl_regression(mod2)

