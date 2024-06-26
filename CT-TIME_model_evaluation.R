###########################################################
#                                                         #
#     CT-TINF model evaluation             #
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
model_path <- here::here( "CT-TIME_model.rds") #input model
load(data_path)
load(signature_path)
glmnet_fit <- readRDS(model_path)
###### Scale data and get outcome ######
### split cohorts
training <- data_cb[data_cb$DataBatch=="VHIO",]
#scale features
training[c(6:114)] <- scale(training[c(6:114)])

############ Feature filtering ############
training_features <- training[,c(selected_feature_names,"Patient","GEP_score")]
training_id <- data.frame(training_features)
med_training <- median(training$GEP_score)
training_id$outcome <- ifelse(training$GEP_score>=med_training,"inflamed","uninflamed")

training_id$pred_prob <- predict(glmnet_fit,training_id ,type = "prob")
training_id$pred  <- predict(glmnet_fit,training_id ,type = "raw")

##### plot training ######
# Create a data frame with original and predicted values
results_bin <- data.frame(Original = training_id$outcome, Predicted = training_id$pred)
results_cont <- data.frame(Original = training_id$outcome, Predicted = training_id$pred_prob[,1])

# Create a boxplot
ggplot(results_cont, aes(x = Original, y = Predicted, fill = Original)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.4) +
  scale_fill_manual(values = c("#b2182b","#969696") )+
  labs(title = "Cohort 1",
       x = "TIME Genomic Category",
       y = "CT-TIME Radiomic Score")

# Create a violin plot
ggplot(results_cont, aes(x = Original, y = Predicted, fill = Original)) +
  geom_violin() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = c("#b2182b", "#969696")) +
  labs(
    title = "Cohort 1: Predicted Radiomic Scores",
    x = "Genomic TIME Category",
    y = "Radiomic Score"
  )

# # ########### TODO: Survival analysis #####
# data_path2 <- here::here("data_processed", "training_clinical_data.csv") #input data
# data <- read.csv2(data_path2)
# #exclude patient 17376935
# data <- data[data$NHC!=17376935,]
# data$Patient <- as.character(data$NHC)
# 
# combined <- merge(training_id,data,by="Patient")
# combined$PFS_months <- (combined$PFS_days)/30.4
# fit <- survfit(Surv(PFS_months, progression) ~ GEP_dich, data = combined)
# print(fit)
# ggsurvplot(fit,
#            pval = TRUE, conf.int = TRUE,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c( "#b2182b", "#4d4d4d"),
#            xlab = 'Progression free survival (months)',
#            break.x.by = 5 )
# 
# ##### Cox models
# # Create a survival object
# surv_obj <- Surv(time = combined$PFS_months, event = combined$progression)
# 
# # Combine the dichotomized predictions and other selected features
# #cox_data <- combined[c("pred_dich")]
# cox_data <- combined[selected_feature_names] 
# # Create a Cox proportional hazards model
# cox_model <- coxph(surv_obj ~ ., data = cox_data)
# # Summarize the model
# summary(cox_model)
# 
# # #Calculate the concordance index
# # c_index <- survConcordance(cox_model, data = cox_data)$concordance
# # # Display the C-index
# # print(paste("Concordance Index (C-index):", round(c_index, 3)))
# # 
# # # Plot the survival curve
# # ggsurvplot(survfit(cox_model, newdata = cox_data), data = cox_data, risk.table = TRUE)
# 
# # Apply the Hosmer-Lemeshow test
# 
# 
# # ###### Univariate Analysis (target: CB 5 months) ######
# # combined$cb <- ifelse(combined$PFS_months>=5,1,0)
# # outcome_variable <- "cb"
# # selected_formula <- as.formula(paste(outcome_variable, "~", paste(names(combined[c(11:24,35)]), collapse = "+")))
# # 
# # combined[c(11:24,35)] <- scale(combined[c(11:24,35)])
# # tbl_uvregression(  combined[c(11:24,35,46)],                       
# #                    method = glm,
# #                    y = cb,                               
# #                    method.args = list(family = binomial),  ## define what type of glm want to run (logistic)
# #                    exponentiate = TRUE                     ## exponentiate to produce odds ratios (rather than log odds)
# # )
# # 
# # ###### Multivariate Analysis ######
# # mod2 <- glm(selected_formula, data= combined[c(11:24,35,46)] , 
# #             family=binomial(link = "logit"))
# # tbl_regression(mod2, exponentiate = TRUE)
# 
# 
# ###### Univariate Analysis (target: PFS) ######
# time_variable <- "PFS_months"
# event_variable <- "progression"
# 
# #making formulas
# univ_formulas <- sapply(names(combined[c(11:24,35)]),function(x)as.formula(paste("Surv(", time_variable, ",", event_variable, ") ~",x)))
# #making a list of models
# univ_models <- lapply(univ_formulas, function(x){coxph(x,data=combined)})
# #extract data (here I've gone for HR and confint)
# univ_results <- lapply(univ_models,function(x){return(exp(cbind(coef(x),confint(x))))})
# 
# 
# df_univ_combined <- do.call(rbind, lapply(univ_results, function(result) {
#   result_df <- as.data.frame(result)
#   result_df$model <- "Univariate"
#   return(result_df)
# }))
# 
# # Move row names to the first column
# df_univ_combined <- rownames_to_column(df_univ_combined, var = "label")
# 
# # Multivariate Cox Proportional Hazards Model
# selected_formula_cox_multi <- as.formula(paste("Surv(", time_variable, ",", event_variable, ") ~", paste(names(combined[c(11:24,35)]), collapse = "+")))
# cox_multi <- coxph(selected_formula_cox_multi, data = combined)
# 
# # Extracting coefficients and confidence intervals
# df_multi_combined <- exp(cbind(coef(cox_multi), confint(cox_multi)))
# df_multi_combined <- as.data.frame(df_multi_combined)
# df_multi_combined$model <- "Multivariate"  # Add a label for multivariate model
# # Move row names to the first column
# df_multi_combined <- rownames_to_column(df_multi_combined, var = "label")
# 
# 
# # Combine univariate and multivariate results
# df_combined <- rbind(df_univ_combined, df_multi_combined)
# df_combined$var_label <- "immune marker"
# 
# # Rename the columns
# colnames(df_combined) <- c("label","estimate", "conf.low", "conf.high","cox_model","var_label")
# 
# ######### UVA MVA plot ######
# dotCOLS <- c("#bababa", "#4d4d4d")
# barCOLS <- c("#bababa", "#4d4d4d")
# 
# 
# ggplot(df_combined, aes(x = estimate, xmax = conf.high, xmin = conf.low, y = label, color = cox_model, fill = cox_model)) +
#   geom_point(size = 3, shape = 18, position = position_dodge(width = 0.5)) +
#   geom_linerange(position = position_dodge(width = 0.5), size = 1) +
#   geom_vline(xintercept = 1, size = 1) +
#   #geom_hline(yintercept = 0, size = 1) +
#   facet_grid(var_label ~ ., scales = "free_y", space = "free_y") +
#   scale_alpha_identity() +
#   scale_fill_manual(values = barCOLS) +
#   scale_color_manual(values = dotCOLS) +
#   scale_y_discrete(name = "Immune Marker") +
#   scale_x_continuous(name = "Hazard ratio", limits = c(0.93, 1.08)) +
#   coord_cartesian(clip = "off") +
#   theme(
#     panel.background = element_blank(),
#     panel.spacing = unit(0, "pt"),
#     axis.line.x.bottom = element_line(size = 1),
#     axis.text.y.left =  element_text(margin = margin(l = 20, unit = "pt")),
#     strip.background = element_blank(), 
#     strip.text = element_blank()
#   )
# 
# # 
# # 
# # ########### simple LASSO (selected Tregs and macrophages)
# # library(glmnet)
# # set.seed(1234)
# # # Fit LASSO Cox regression
# # combined <- combined[combined$PFS_months > 0, ]
# # lasso_cox_model <- cv.glmnet(as.matrix(combined[c(11:24,35)]), Surv(combined$PFS_months, combined$progression), family = "cox", alpha = 1)
# # # Plot the cross-validated error as a function of lambda
# # plot(lasso_cox_model)
# # # Choose the best lambda
# # best_lambda <- lasso_cox_model$lambda.min
# # # Refit the model with the selected lambda
# # lasso_cox_final <- glmnet(as.matrix(combined[c(11:24,35)]), Surv(combined$PFS_months, combined$progression), family = "cox", alpha = 1, lambda = best_lambda)
# # # Display coefficients
# # coef(lasso_cox_final)
# # 
# # ##### Lasso Cox feature selection (T reg)---------------------------------------------
# # 
# # # Split data-----------------------------------------------------
# # num_vars <- names(combined[c(11:24,35)]) # Genomics
# # combined[c(11:24,35)] <- scale(combined[c(11:24,35)])
# # combined <- combined[combined$PFS_months > 0, ]
# # covariates <- c(num_vars,"PFS_months","progression")
# # features <- subset(combined, select = covariates )
# # 
# # # TRAIN TEST 80-20%
# # seed_n = 53 #56 svr and mean
# # set.seed(seed_n); rnorm(10)
# # split_data1 = partition(features, p = 0.8, cat_col = 'progression' )
# # train = split_data1[[1]]
# # test = split_data1[[2]]
# # 
# # set.seed(seed_n); rnorm(10)
# # n_seeds = sample(1:1500, 10, replace=FALSE)#
# # values = data.frame()
# # for (i in n_seeds) {
# #   
# #   # Feature selection: LASSO elastic net for cox model (most common features from different seeds) ------------------------
# #   set.seed(i)
# #   # model cv
# #   X <- na.omit(train[c('PFS_months', 'progression', num_vars)])
# #   names(X)[0:2] = c('time', 'status')
# #   x <- model.matrix(as.formula(paste('Surv(time, status) ~', paste(num_vars, collapse = '+'))), X)
# #   y = Surv(X$time, X$status)
# #   fit <- cv.glmnet(as.matrix(x), as.matrix(y), family = "cox")
# #   md1 <- glmnet(as.matrix(x), as.matrix(y), family = "cox", lambda = fit$lambda.min )
# #   fitted_vals <- predict(md1, newx = x)
# #   myCoefs1 <- coef(md1)
# #   vars_nonzero1 = myCoefs1@Dimnames[[1]][ which(myCoefs1 != 0 ) ] 
# #   
# #   set.seed(i)
# #   X2 <- na.omit(test[c('PFS_months', 'progression', num_vars)])
# #   names(X2)[0:2] = c('time', 'status')
# #   newx <- model.matrix(as.formula(paste('Surv(time, status) ~', paste( num_vars, collapse = '+'))), X2)
# #   newy = Surv(X$time, X$status)
# #   survival.results0 <- predict(md1, newx = newx)
# #   
# #   w.ROCtrain = risksetROC(Stime = X$time,  
# #                           status = X$status, 
# #                           marker = fitted_vals,
# #                           predict.time =7, plot = FALSE,
# #                           method = "Cox")
# #   
# #   w.ROCtest = risksetROC(Stime = X2$time,  
# #                          status = X2$status, 
# #                          marker = survival.results0,
# #                          predict.time = 7, plot = FALSE,
# #                          method = "Cox")
# #   
# #   
# #   print(paste0("N_seed = ", i,
# #                " - Train: ", round(w.ROCtrain$AUC,4),
# #                ", Test: ", round(w.ROCtest$AUC,4),
# #                ", Features: ", vars_nonzero1))
# #   
# #   sub <- cbind( i, round(w.ROCtrain$AUC,4), round(w.ROCtest$AUC,4),vars_nonzero1)   
# #   values <- rbind(values,sub)# Adding new variable to data
# #   
# #   
# # }
# # 
# # colnames(values) = c('N Seeds', 'Train',  'Test', 'Feature')
# # 
# # selected_frequency <- as.data.frame(table(values$Feature))
# # 
# # selected_frequency
# # 
# # most_frequent <- selected_frequency$Var1
# # #most_frequent <- as.character(selected_frequency$Var1[selected_frequency$Freq>median(selected_frequency$Freq)])
# # #most_frequent <- as.character(selected_frequency$Var1[selected_frequency$Freq>unname(quantile(selected_frequency$Freq,c(0.75)))])
# # most_frequent2 <- gsub("original_", "", most_frequent)
# # print(most_frequent2)


