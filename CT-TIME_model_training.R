###########################################################
#                                                         #
#     CT-TIME model training              #
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
signature_path <- here::here("CT-TIME_signature.Rda") #input data
model_path <- here::here( "CT-TIME_model.rds") #output model
load(data_path)
load(signature_path)
###### Scale data and get outcome ######
### split cohorts
training <- data_cb[data_cb$DataBatch=="VHIO",]
#scale numeric features
training[c(6:114)] <- scale(training[c(6:114)])

############ Feature filtering ############
training_features <- training[,c(selected_feature_names)]
training_df <- data.frame(training_features)
med_training <- median(training$GEP_score)
training_df$outcome <- ifelse(training$GEP_score>=med_training,"inflamed","uninflamed")

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


############ Model Training: Logistic Regression ############
outcome_variable <- "outcome"
selected_formula <- as.formula(paste(outcome_variable, "~", paste(selected_feature_names , collapse = "+")))

#parameters from cv
bestalpha = 0.55
bestlambda = 0.007

set.seed(1234)
cv_5 = trainControl(method = "repeatedcv", number = 5, repeats=10, returnResamp="all",
                    classProbs=TRUE, summaryFunction=twoClassSummary)

glmnet_fit = train(
  selected_formula , data = training_df,
  method = "glmnet", 
  trControl = cv_5, metric = "ROC",
  tuneGrid = expand.grid(alpha = bestalpha,
                         lambda = bestlambda))


# Get logistic regression model coefficients (reference category: inflamed)
myCoef <- coef(glmnet_fit$finalModel, bestlambda) 

# Save final model
saveRDS(glmnet_fit, model_path)

############ Model cross-validation ############
set.seed(1234)
outcome_variable <- "outcome"
selected_formula_hl <- as.formula(paste(outcome_variable, "~", paste(selected_feature_names , collapse = "+")))

mean_AUC <- mean(glmnet_fit$resample$ROC) #0.78
sd_AUC <- sd(glmnet_fit$resample$ROC) #0.18

mean_Sens <- mean(glmnet_fit$resample$Sens) #0.70
sd_Sens <- sd(glmnet_fit$resample$Sens) #0.26
mean_Spec <- mean(glmnet_fit$resample$Spec) #0.68
sd_Spec <- sd(glmnet_fit$resample$Spec) #0.26
