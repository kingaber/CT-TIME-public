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
library(pROC)
library(PRROC)
library(SuperLearner)# ensemble
###### Define paths and load data ######
data_path <- here::here("..","1_data_processed", "harmonized_NANO_TCIA_data.Rda") #input data
signature_path <- here::here("CT-TIME_signature.Rda") #input data
logreg_path <- here::here( "CT-TIME_model_logreg.rds") #output model
dt_path <- here::here( "CT-TIME_model_dt.rds") #output model
rf_path <- here::here( "CT-TIME_model_rf.rds") #output model
svm_path <- here::here( "CT-TIME_model_svm.rds") #output model
knn_path <- here::here( "CT-TIME_model_knn.rds") #output model
nn_path <- here::here( "CT-TIME_model_nn.rds") #output model
ensemble_path <- here::here( "CT-TIME_model_ensemble.rds") #output model
load(data_path)
load(signature_path)
###### Scale data and get outcome ######
### split cohorts
training <- data_cb[data_cb$DataBatch=="VHIO",]
#scale features
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

############ Models training
##### ML balanced models #####
# Custom trainControl object for cross-validation
num_folds <- 5 
num_repeats <- 10
set.seed(1234)
custom_control <- trainControl(
  method = "repeatedcv",
  number = num_folds,
  repeats = num_repeats,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = TRUE
)


## Training and get metrics
models <- c("glm","rpart","rf","svmRadial","knn","nnet")
# Create an empty list to store results for each model
results_list <- list()
# Create an empty list to store results for each model
cv_results <- list()

for (model_name in models) {
  
  model <- train(
    selected_formula,  
    data = training_df,
    method = model_name,  
    trControl = custom_control
  )
  
  # Save the trained model
  saveRDS(model, file.path(getwd(), paste0(model_name, "_CT_TIME.rds")))
  
  # Make predictions on training and test sets
  train_predicted_prob <- predict(model, newdata = training_df, type = "prob")
  train_predicted_class <- predict(model, newdata = training_df, type = "raw")
  
  # Evaluate performance metrics on the test set
  conf_matrix <- confusionMatrix(train_predicted_class, reference = as.factor(training_df$outcome))
  sensitivity <- conf_matrix$byClass[1]
  specificity <- conf_matrix$byClass[2]
  precision <- conf_matrix$byClass[5]
  recall <- conf_matrix$byClass[6]
  
  # ROC curve and ROC-AUC
  outcome_numeric <- as.numeric(as.factor(training_df$outcome)) - 1 
  roc_obj <- roc(outcome_numeric, train_predicted_prob[[1]])
  roc_auc <- auc(roc_obj)
  
  # Precision-Recall curve and PR-AUC
  pr_obj <- pr.curve(outcome_numeric, train_predicted_prob[[1]])
  pr_auc <- pr_obj$auc.integral
  
  # Store metrics in the results list
  results_list[[model_name]] <- data.frame(
    Model = model_name,
    Precision = precision,
    Recall = recall,
    Specificity = specificity,
    ROC_AUC = roc_auc,
    PR_AUC = pr_auc
  )
  

  # Calculate mean and standard deviation for AUC, Sensitivity, Specificity
  mean_AUC <- mean(model$resample$ROC)
  sd_AUC <- sd(model$resample$ROC)
  
  mean_Sens <- mean(model$resample$Sens)
  sd_Sens <- sd(model$resample$Sens)
  
  mean_Spec <- mean(model$resample$Spec)
  sd_Spec <- sd(model$resample$Spec)
  
  # Store the results
  cv_results[[model_name]] <- data.frame(
    Model = model_name,
    Mean_AUC = mean_AUC,
    SD_AUC = sd_AUC,
    Mean_Sens = mean_Sens,
    SD_Sens = sd_Sens,
    Mean_Spec = mean_Spec,
    SD_Spec = sd_Spec
  )
  
}
  

# Combine results from the list into a single data frame
results_df <- do.call(rbind, results_list) #training
cv_results_df <- do.call(rbind, cv_results) #cv

# Print or further analyze the results
print(cv_results_df)


# Combine results from the list into a single data frame
cv_results_long <- tidyr::pivot_longer(cv_results_df, cols = starts_with("Mean_"), names_to = "Metric", values_to = "Value")

# Create a ggplot with error bars
ggplot(cv_results_long, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_errorbar(
    aes(ymin = Value - cv_results_long$SD_AUC, ymax = Value + cv_results_long$SD_AUC),
    position = position_dodge(0.9),
    width = 0.4,
    colour = "#636363",
    alpha = 0.9,
    size = 0.5
  ) +
  labs(title = "Models cross validation", x = "Model", y = "Value") +
  theme_minimal()



###### Ensemble model


##### Ensemble model ######
#listWrappers()
set.seed(1234)
train_data_sl <- training_df[c(1:4)]
train_labels_sl <- as.numeric(as.factor(training_df$outcome)) - 1 

sl = SuperLearner(Y = train_labels_sl, X =train_data_sl, family = binomial(),cvControl=list(V=5),
                  SL.library = c( "SL.glm" ,"SL.nnet" ,"SL.randomForest"    ))
predictions_train <- predict.SuperLearner(sl, newdata=train_data_sl,onlySL=TRUE)
predictions <- predictions_train

# Recode probabilities
conv.preds <- ifelse(predictions$pred>0.5,1,0)
# Create the confusion matrix
a <- as.factor(conv.preds)
b <- as.factor(train_labels_sl)
# Reorder the levels of 'a' to match the order of 'b'
a <- factor(a, levels = levels(b))
cm <- confusionMatrix(a, b)
# Return the confusion matrix
cm
# Calculate Sensitivity, Specificity, Precision, and Recall
sensitivity <- cm$byClass[1]
specificity <-  cm$byClass[2]
precision <- cm$byClass[5]
recall <- cm$byClass[6]
# ROC curve and ROC-AUC
roc_obj <- roc(b, as.numeric(a))
roc_auc <- auc(roc_obj)
# Precision-Recall curve and PR-AUC
pr_obj <- pr.curve(b, as.numeric(a))
pr_auc <- pr_obj$auc.integral

# Output the performance metrics
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("Specificity:", specificity, "\n")
cat("ROC-AUC:", roc_auc, "\n")
cat("PR-AUC:", pr_auc, "\n")

# Save the trained model
saveRDS(sl, file.path(getwd(), "Ensemble_CT_TIME.rds"))