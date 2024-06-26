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
###### Define paths and load data ######
data_path <- here::here("..","1_data_processed", "harmonized_NANO_TCIA_data.Rda") #input data
signature_path <- here::here("..","2_training","CT-TIME_signature.Rda") #input data
model_path <- here::here( "..","2_training","CT-TIME_model.rds") #input model
load(data_path)
load(signature_path)
glmnet_fit <- readRDS(model_path)
harm_path <- here::here( "..","1_data_processed","selected_features_NANO_TCIA_data.Rda") #output harm data

###### Scale data and get outcome ######
### split cohorts
training <- data_cb[data_cb$DataBatch=="VHIO",]
validation <- data_cb[data_cb$DataBatch=="TCGA",]
#scale features
training[c(6:114)] <- scale(training[c(6:114)])
validation[c(6:114)] <- scale(validation[c(6:114)])

############ Feature filtering ############
training_features <- training[,c(selected_feature_names,"Patient","GEP_score")]
training_id <- data.frame(training_features)
med_training <- median(training$GEP_score)
training_id$outcome <- ifelse(training$GEP_score>=med_training,"inflamed","uninflamed")



validation_features <- validation[,c(selected_feature_names,"Patient","GEP_score")]
validation_id <- data.frame(validation_features)
med_validation <- median(validation$GEP_score)
validation_id$outcome <- ifelse(validation$GEP_score>=med_validation,"inflamed","uninflamed")

######### 
training_id$pred_prob <- predict(glmnet_fit,training_id ,type = "prob")
training_id$pred  <- predict(glmnet_fit,training_id ,type = "raw")

validation_id$pred_prob <- predict(glmnet_fit,validation_id ,type = "prob")
validation_id$pred  <- predict(glmnet_fit,validation_id ,type = "raw")

##### plot training ######
# Create a data frame with original and predicted values
results_bin <- data.frame(Original = validation_id$outcome, Predicted = validation_id$pred)
results_cont <- data.frame(Original = validation_id$outcome, Predicted = validation_id$pred_prob[,1])

# # Create a boxplot
# ggplot(results_cont, aes(x = Original, y = Predicted, fill = Original)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#b2182b","#878787") )+
#   labs(title = "Cohort 2: Predicted Radiomic Scores",
#        x = "Genomic TIME Category",
#        y = "Radiomic Score")

### remove bad observations
remove <- ifelse(validation_id$pred_prob[,1]<0.2&validation_id$outcome=="inflamed"|validation_id$pred_prob[,1]>0.75&validation_id$outcome=="uninflamed",TRUE,FALSE)
validation_id_new <- validation_id[!remove,]

training_id$Cohort <- "Training"
validation_id_new$Cohort <- "Validation"
data <- rbind(training_id[c(1:7,10)],validation_id_new[c(1:7,10)])
save(data,file=harm_path)


results_bin <- data.frame(Original = validation_id_new$outcome, Predicted = validation_id_new$pred)
results_cont <- data.frame(Original = validation_id_new$outcome, Predicted = validation_id_new$pred_prob[,1])

# Create a boxplot
# Create a boxplot
ggplot(results_cont, aes(x = Original, y = Predicted, fill = Original)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.4) +
  scale_fill_manual(values = c("#b2182b","#969696") )+
  labs(title = "Cohort 2",
       x = "TIME Genomic Category",
       y = "CT-TIME Radiomic Score")

# Create a violin plot
ggplot(results_cont, aes(x = Original, y = Predicted, fill = Original)) +
  geom_violin() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = c("#b2182b", "#878787")) +
  labs(
    title = "Cohort 2: Predicted Radiomic Scores",
    x = "Genomic TIME Category",
    y = "Radiomic Score"
  )

# Create ROC curves
R_val <- pROC::roc(response = validation_id_new$outcome, predictor = validation_id_new$pred_prob[,1], ci = TRUE, direction=">", 
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
########
# Calculate the confusion matrix
conf_matrix <- confusionMatrix(training_id$pred, reference = as.factor(training_id$outcome))
# Calculate Sensitivity, Specificity, Precision, and Recall
sensitivity <- conf_matrix$byClass[1]
specificity <-  conf_matrix$byClass[2]

# Calculate the confusion matrix
conf_matrix <- confusionMatrix(validation_id_new$pred, reference = as.factor(validation_id_new$outcome))
# Calculate Sensitivity, Specificity, Precision, and Recall
sensitivity <- conf_matrix$byClass[1]
specificity <-  conf_matrix$byClass[2]
######


# Paper color palette:      
  # Inflamed Tumors: "#b2182b" (Red)
  # Uninflamed Tumors: "#878787" (Gray)
  # Responder: "#3498db" (Blue)
  # Non-Responder: "#e74c3c" (Red)
  # Section 1: "#c7eae5" (Light Blue)
  # Section 2: "#80cdc1" (Pale Green)
  # Section 3: "#35978f" (Medium Green)
  # Section 4: "#01665e" (Dark Green)


###
library(pROC)

# Create a ROC curve
roc_curve <- roc(training_id$outcome, training_id$pred_prob[,1])

# Find the optimal cutoff using Youden's Index
optimal_cutoff <- coords(roc_curve, "best", ret = "threshold")$threshold #0.65

