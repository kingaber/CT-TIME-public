###########################################################
#                                                         #
#     CT-TIME model evaluation             #
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
data_path <- here::here("..","1_data_processed", "selected_features_NANO_TCIA_data.Rda") #input data
signature_path <- here::here("..","2_training","CT-TIME_signature.Rda") #input data
load(data_path)
load(signature_path)

###### Load and scale data and get outcome ######
### split cohorts
training <- data[data$Cohort=="Training",]
validation <- data[data$Cohort=="Validation",]

############ Feature filtering ############
training_id <- data.frame(training)
med_training <- median(training$GEP_score)
training_id$outcome <- ifelse(training$GEP_score>=med_training,"inflamed","uninflamed")

validation_id <- data.frame(validation)
med_validation <- median(validation$GEP_score)
validation_id$outcome <- ifelse(validation$GEP_score>=med_validation,"inflamed","uninflamed")

##### Evaluate models
models <- c("glm","rpart","rf","svmRadial","knn","nnet")
# Create an empty list to store results for each model
train_results_list <- list()
val_results_list <- list()

for (model_name in models) {
  # Load model
  model <- readRDS(file.path(here::here("..","2_training"), paste0(model_name, "_CT_TIME.rds")))
  # Make predictions on training and test sets
  train_predicted_prob <- predict(model, newdata = training_id, type = "prob")
  train_predicted_class <- predict(model, newdata = training_id, type = "raw")
  val_predicted_prob <- predict(model, newdata = validation_id, type = "prob")
  val_predicted_class <- predict(model, newdata = validation_id, type = "raw")
  
  # Evaluate performance metrics on the test set
  conf_matrix <- confusionMatrix(train_predicted_class, reference = as.factor(training_id$outcome))
  sensitivity <- conf_matrix$byClass[1]
  specificity <- conf_matrix$byClass[2]
  precision <- conf_matrix$byClass[5]
  recall <- conf_matrix$byClass[6]
  
  # ROC curve and ROC-AUC
  outcome_numeric <- as.numeric(as.factor(training_id$outcome)) - 1 
  roc_obj <- roc(outcome_numeric, train_predicted_prob[[1]])
  roc_auc <- auc(roc_obj)
  
  # Precision-Recall curve and PR-AUC
  pr_obj <- pr.curve(outcome_numeric, train_predicted_prob[[1]])
  pr_auc <- pr_obj$auc.integral
  
  # Store metrics in the results list
  train_results_list[[model_name]] <- data.frame(
    Model = model_name,
    Precision = precision,
    Recall = recall,
    Specificity = specificity,
    ROC_AUC = roc_auc,
    PR_AUC = pr_auc
  )
  
  # Evaluate performance metrics on the test set
  conf_matrix <- confusionMatrix(val_predicted_class, reference = as.factor(validation_id$outcome))
  sensitivity <- conf_matrix$byClass[1]
  specificity <- conf_matrix$byClass[2]
  precision <- conf_matrix$byClass[5]
  recall <- conf_matrix$byClass[6]
  
  # ROC curve and ROC-AUC
  outcome_numeric <- as.numeric(as.factor(validation_id$outcome)) - 1 
  roc_obj <- roc(outcome_numeric, val_predicted_prob[[1]])
  roc_auc <- auc(roc_obj)
  
  # Precision-Recall curve and PR-AUC
  pr_obj <- pr.curve(outcome_numeric, val_predicted_prob[[1]])
  pr_auc <- pr_obj$auc.integral
  
  # Store metrics in the results list
  val_results_list[[model_name]] <- data.frame(
    Model = model_name,
    Precision = precision,
    Recall = recall,
    Specificity = specificity,
    ROC_AUC = roc_auc,
    PR_AUC = pr_auc
  )
  
  
} 

# Combine results from the list into a single data frame
train_results_df <- do.call(rbind, train_results_list) 
validation_results_df <- do.call(rbind, val_results_list) 


# Convert relevant columns to character
train_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")] <- 
  lapply(train_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")], as.character)

# Reshape the data
train_results_long <- tidyr::pivot_longer(train_results_df, cols = c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC"), names_to = "Metric", values_to = "Value")

# Create a ggplot
ggplot(train_results_long, aes(x = Model, y = as.numeric(Value), fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  labs(title = "Models Training Results", x = "Model", y = "Value") +
  theme_minimal()


# Convert relevant columns to character
validation_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")] <- 
  lapply(validation_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")], as.character)

# Reshape the data
validation_results_long <- tidyr::pivot_longer(validation_results_df, cols = c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC"), names_to = "Metric", values_to = "Value")

# Create a ggplot
ggplot(validation_results_long, aes(x = Model, y = as.numeric(Value), fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  labs(title = "Models Validation Results", x = "Model", y = "Value") +
  theme_minimal()

### Best models ###
#rf, nnet, glm
# Load model
model <- readRDS(file.path(here::here("..","2_training"), "Ensemble_CT_TIME.rds"))
predictions_train <- predict.SuperLearner(model, newdata=training_id[c(1:4)],onlySL=TRUE)
predictions_val <- predict.SuperLearner(model, newdata=validation_id[c(1:4)],onlySL=TRUE)
### TRAIN
# Recode probabilities
conv.preds <- ifelse(predictions_train$pred>0.5,1,0)
# Create the confusion matrix
a <- as.factor(conv.preds)
b <- as.factor(as.numeric(as.factor(training_id$outcome)) - 1 )
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

en_train_results_list <- list()
en_train_results_list[["ensemble"]] <- data.frame(
  Model = "ensemble",
  Precision = precision,
  Recall = recall,
  Specificity = specificity,
  ROC_AUC = roc_auc,
  PR_AUC = pr_auc
)

# Convert each list to a data frame
en_train_results_df <- as.data.frame(en_train_results_list$ensemble)
train_results_df <- bind_rows(
  as.data.frame(train_results_list$glm),
  as.data.frame(train_results_list$rpart),
  as.data.frame(train_results_list$rf),
  as.data.frame(train_results_list$svmRadial),
  as.data.frame(train_results_list$knn),
  as.data.frame(train_results_list$nnet)
)

# Convert relevant columns to character for ensemble results
en_train_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")] <-
  lapply(en_train_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")], as.character)

# Convert relevant columns to character for individual model results
train_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")] <-
  lapply(train_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")], as.character)

# Combine the data frames
combined_df <- bind_rows(en_train_results_df, train_results_df)

# Reshape the data
combined_df_long <- tidyr::pivot_longer(combined_df, cols = c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC"), names_to = "Metric", values_to = "Value")

# Create a ggplot
ggplot(combined_df_long, aes(x = Model, y = as.numeric(Value), fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  labs(title = "Models Training Results", x = "Model", y = "Value") +
  theme_minimal()


### VAL
# Recode probabilities
conv.preds <- ifelse(predictions_val$pred>0.5,1,0)
# Create the confusion matrix
a <- as.factor(conv.preds)
b <- as.factor(as.numeric(as.factor(validation_id$outcome)) - 1 )
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

en_val_results_list <- list()
en_val_results_list[["ensemble"]] <- data.frame(
  Model = "ensemble",
  Precision = precision,
  Recall = recall,
  Specificity = specificity,
  ROC_AUC = roc_auc,
  PR_AUC = pr_auc
)

# Convert each list to a data frame
en_val_results_df <- as.data.frame(en_val_results_list$ensemble)

# Helper function to convert auc to character
convert_auc_to_character <- function(df) {
  df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")] <-
    lapply(df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")], as.numeric)
  return(df)
}
# Convert relevant columns to character for validation results
val_results_df <- bind_rows(
  convert_auc_to_character(as.data.frame(val_results_list$glm)),
  convert_auc_to_character(as.data.frame(val_results_list$rpart)),
  convert_auc_to_character(as.data.frame(val_results_list$rf)),
  convert_auc_to_character(as.data.frame(val_results_list$svmRadial)),
  convert_auc_to_character(as.data.frame(val_results_list$knn)),
  convert_auc_to_character(as.data.frame(val_results_list$nnet))
)

# Convert relevant columns to character
val_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")] <-
  lapply(val_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")], as.character)


# Convert relevant columns to character for ensemble results
en_val_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")] <-
  lapply(en_val_results_df[, c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC")], as.character)

# Combine the data frames
val_combined_df <- bind_rows(en_val_results_df,val_results_df)

# Reshape the data
combined_df_long <- tidyr::pivot_longer(val_combined_df, cols = c("Precision", "Recall", "Specificity", "ROC_AUC", "PR_AUC"), names_to = "Metric", values_to = "Value")

# Create a ggplot
ggplot(combined_df_long, aes(x = Model, y = as.numeric(Value), fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  labs(title = "Models Validation Results", x = "Model", y = "Value") +
  theme_minimal()

# Create ROC curves
R_val <- pROC::roc(response = validation_id$outcome, predictor = validation_id_new$pred_prob[,1], ci = TRUE, direction=">", 
                   xlab = "False Positive Rate", ylab = "True Positive Rate", 
                   main = 'CT Radiomics Signature \n Model Logistic Regression - Training and Validation Set', 
                   cex.main = 1, col = "#b2182b")
R_train <- pROC::roc(response = training_id$outcome, predictor = training_id$pred_prob[,1], ci = TRUE, direction=">", 
                     add = TRUE, col = "#969696")

# Plot both ROC curves on one plot
plot(R_val, legacy.axes = TRUE, col ="#c7eae5", lwd = 3,lty = 3)
plot(R_train, add = TRUE, col = "#80cdc1", lwd = 3 )
legend('bottomright', legend = c(paste('Cohort 2: ',format(round(R_val$auc,2),nsmall = 2),
                                       ' (',format(round(R_val$ci[1],2),nsmall = 2),'-',format(round(R_val$ci[3],2),nsmall = 2),')',sep =''),
                                 paste('Cohort 1: ',format(round(R_train$auc,2),nsmall = 2),
                                       ' (',format(round(R_train$ci[1],2),nsmall = 2),'-',format(round(R_train$ci[3],2),nsmall = 2),')',sep ='')),
       lty = c(1, 1), lwd = c(2.5, 2.5), col = c("#c7eae5", "#80cdc1"), cex = 0.8, bg = "white")

R_train$auc;R_train$ci

# Paper color palette:      
  # Inflamed Tumors: "#b2182b" (Red)
  # Uninflamed Tumors: "#878787" (Gray)
  # Responder: "#3498db" (Blue)
  # Non-Responder: "#e74c3c" (Red)
  # Section 1: "#c7eae5" (Light Blue)
  # Section 2: "#80cdc1" (Pale Green)
  # Section 3: "#35978f" (Medium Green)
  # Section 4: "#01665e" (Dark Green)

