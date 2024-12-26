library(survival)
library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)
library(rms)
source("bootstrap.R")
source("harrell_Uno.R")
source("calibrate_pac_fct.R")

#############################################################

# Function to select significant variables based on univariate Cox models
select_significant_vars <- function(training_data, set_number = 0, skip = FALSE) {
  covariates_set1 <- c("CEA_binary_20","KRAS","N_stage", "TBS_binary","extrahepatic.disease")
  covariates_set2 <- c("size_binary","KRAS","N_stage")
  covariates_set3 <- c("size_binary","number_binary","DFI.12", "N_stage", "CEA_binary_200")
  
  covariates <- if (set_number == 1) {
    covariates_set1
  } else if (set_number == 2) {
    covariates_set2
  } else if (set_number == 3) {
    covariates_set3
  } else if (set_number == 0) {
    c("age", "T_stage", "N_stage", "rightleft.Rec",  "cea_carcinoembryonic.antigen", "extrahepatic.disease", 
      "size", "number_liver_mets", "R0", "KRAS")
  }
  
  
  # If skip is TRUE, return the covariates as is
  if (skip) {
    return(covariates)
  }
  
  # Initialize an empty vector to store significant variables
  significant_vars <- c()
  
  # Loop through the covariates
  for (var in covariates) {
    # Dynamically create the formula
    formula <- as.formula(paste("Surv(time_to_event, died) ~", var))
    
    # Run the Cox proportional hazards model
    res.cox <- coxph(formula, data = training_data)
    
    # Extract the p-value from the summary
    p_value <- summary(res.cox)$coefficients[, "Pr(>|z|)"]
    
    # Check if the p-value is less than 0.05
    if (!is.na(p_value) && p_value < 0.05) {
      # If significant, store the variable name
      significant_vars <- c(significant_vars, var)
    }
  }
  
  # Return the significant variables
  return(significant_vars)
}

# Function to extract validation metrics from a validate object
extract_validation_metrics <- function(validation_obj, model_name) {
  # Define the metrics you want to extract
  required_metrics <- c("Dxy", "R2", "Slope", "Intercept", "Emax")
  
  # Initialize a list to store metrics
  metrics <- list(Model = model_name)
  
  # Loop through each required metric and extract its 'index' value
  for (metric in required_metrics) {
    if (metric %in% rownames(validation_obj)) {
      metrics[[metric]] <- validation_obj[metric, "index"]
    } else {
      metrics[[metric]] <- NA
    }
  }
  
  # Convert the list to a data frame
  return(as.data.frame(metrics))
}

oversample_matched_data <- function(train_data, number_times)
{### OVERSAMPLE THE MATCHED COHORT 24/8 #################
  # Replicate the data `number_times` and bind the rows together
  train_data <- do.call(rbind, replicate(number_times, train_data, simplify = FALSE))
  
}


#############################################################

### LOAD UNMATCHED TRAIN DATA #################################
train_filename <- paste0("traindata/unmatched_train_data.csv")
train_data_m1 <- read.csv(file = train_filename, stringsAsFactors = FALSE)
train_data_m1$size_binary <- ifelse(train_data_m1$size <= 5, 0, 1)

untreated_train_data_2a <- subset(train_data_m1, adjuvantchemo == 0)
treated_train_data_2b <- subset(train_data_m1, adjuvantchemo == 1)

### LOAD MATCHED TRAIN DATA #################################
train_data_equalized_3B <- read.csv(file = "traindata/matched_train_equalized_r1.csv", stringsAsFactors = FALSE)
train_data_equalized_3A <- read.csv(file = "traindata/matched_train_equalized_r1.csv", stringsAsFactors = FALSE)
#train_data_matched$size_binary <- ifelse(train_data_matched$size <= 5, 0, 1)
#train_data_equalized$size_binary <- ifelse(train_data_equalized$size <= 5, 0, 1)
train_data_equalized_r1 <- read.csv(file = "traindata/matched_train_equalized_r1.0_relaxed.csv", stringsAsFactors = FALSE)

untreated_train_data_3a <- subset(train_data_equalized_3A, adjuvantchemo == 0)
treated_train_data_3b <- subset(train_data_equalized_3B, adjuvantchemo == 1)
untreated_train_data_3a_r1 <- subset(train_data_equalized_r1, adjuvantchemo == 0)
treated_train_data_3b_r1 <- subset(train_data_equalized_r1, adjuvantchemo == 1)

### LOAD PSM DATA #################################
train_psm <- read.csv(file = "traindata/matchit_train_glm.csv", stringsAsFactors = FALSE)
untreated_psm <- subset(train_psm, adjuvantchemo == 0)
treated_psm <- subset(train_psm, adjuvantchemo == 1)


oversample_n <- 1
untreated_train_data_3a <- oversample_matched_data(untreated_train_data_3a,oversample_n)  # Duplicate all rows
treated_train_data_3b <- oversample_matched_data(treated_train_data_3b, oversample_n)  # Duplicate all rows
untreated_train_data_3a_r1 <- oversample_matched_data(untreated_train_data_3a_r1,oversample_n)  # Duplicate all rows
treated_train_data_3b_r1 <- oversample_matched_data(treated_train_data_3b_r1, oversample_n)  # Duplicate all rows

  
### SELECT COVARIATES #################################
set_number <-0
print("1 ...")
cov_1  <- select_significant_vars(train_data_m1,  set_number,TRUE)
print("2A ...")
cov_2a <-select_significant_vars(untreated_train_data_2a,  set_number,TRUE) 
print("2B ...")
cov_2b <- select_significant_vars(treated_train_data_2b,  set_number,TRUE)
print("3A ...")
cov_3a <-select_significant_vars(untreated_train_data_3a,  set_number,TRUE)
print("3B ... ")
cov_3b <- select_significant_vars(treated_train_data_3b,  set_number,TRUE)
print(cov_3b)

### FORMULAS #################################
formula_cox1 <- as.formula(paste("Surv(time_to_event, died) ~", paste(c(cov_1, "adjuvantchemo"), collapse = "+")))
formula_others <- as.formula(paste("Surv(time_to_event, died) ~", paste(cov_1, collapse = "+")))

### DO NOT FORGET TO TRAIN COX MODELS ON DUPLICATED DATA ################################
#cox1 <- coxph(as.formula(paste("Surv(time_to_event, died) ~", paste(c(cov_1, "adjuvantchemo"), collapse = "+"))), x = TRUE, data = train_data_m1)
cox1 <- cph(formula_cox1, x = TRUE,y = TRUE, 
            data = train_data_m1,surv = TRUE)
cox_untreated_2a <- cph(as.formula(paste("Surv(time_to_event, died) ~", paste(cov_2a, collapse = "+"))), x = TRUE,y = TRUE,
                        data = untreated_train_data_2a,surv = TRUE) # Model 2A
cox_treated_2b <- cph(as.formula(paste("Surv(time_to_event, died) ~", paste(cov_2b, collapse = "+"))), x = TRUE, y = TRUE,
                      data = treated_train_data_2b,surv = TRUE)    # Model 2B

cox_untreated_3a <- cph(as.formula(paste("Surv(time_to_event, died) ~", paste(cov_3a, collapse = "+"))), x = TRUE, y = TRUE,
                        data = untreated_train_data_3a,surv = TRUE) # Model 3A exact equalized
cox_untreated_3a_r1 <- cph(as.formula(paste("Surv(time_to_event, died) ~", paste(cov_3a, collapse = "+"))), x = TRUE, y = TRUE,
                           data = untreated_train_data_3a_r1,surv=TRUE )     # Model 3A relaxed equalized data

cox_treated_3b <- cph(as.formula(paste("Surv(time_to_event, died) ~", paste(cov_3b, collapse = "+"))), x = TRUE, y = TRUE,
                      data = treated_train_data_3b,surv = TRUE)     # Model 3B
cox_treated_3b_r1 <- cph(as.formula(paste("Surv(time_to_event, died) ~", paste(cov_3b, collapse = "+"))), x = TRUE, y = TRUE,
                           data = treated_train_data_3b_r1,surv=TRUE )     # Model 3B  r=1
cox_treated_psm <-cph(as.formula(paste("Surv(time_to_event, died) ~", paste(cov_3b, collapse = "+"))), x = TRUE, y = TRUE,
              data = treated_psm,surv=TRUE )     # PSM model
cox_untreated_psm <-cph(as.formula(paste("Surv(time_to_event, died) ~", paste(cov_3b, collapse = "+"))), x = TRUE, y = TRUE,
              data = untreated_psm,surv=TRUE )     # PSM model

################################

## Create the combined list of lists ################################

combined_list <- list(
  list(
    name = "1",
    model = cox1,
    data = train_data_m1,
    formula = formula_cox1
  ),
  list(
    name = "2A",
    model = cox_untreated_2a,
    data = untreated_train_data_2a,
    formula = formula_others
  ),
  list(
    name = "3A 1to1",
    model = cox_untreated_3a,
    data = untreated_train_data_3a,
    formula = formula_others
    
  ),
  list(
    name = "3A relaxed",
    model = cox_untreated_3a_r1,
    data = untreated_train_data_3a_r1,
    formula = formula_others
 
  ),
  list(
    name = "PSM A GLM",
    model = cox_untreated_psm,
    data = untreated_psm,
    formula = formula_others
  ),
  list(
    name = "2B",
    model = cox_treated_2b,
    data = treated_train_data_2b,
    formula = formula_others
  
  ),list(
    name = "3B 1to1",
    model = cox_treated_3b,
    data = treated_train_data_3b,
    formula = formula_others)
  ,list(
    name = "3B relaxed",
    model = cox_treated_3b_r1,
    data = treated_train_data_3b_r1,
    formula = formula_others),
  list(
    name = "PSM B GLM",
    model = cox_treated_psm,
    data = treated_psm,
    formula = formula_others
  )
)

extra_list <- list(
  list(
    name = "3A 1to1 dup",
    model = cox_untreated_3a,
    data = untreated_train_data_3a,
    formula = formula_others
    
  ),
  list(
    name = "3A relaxed dup",
    model = cox_untreated_3a_r1,
    data = untreated_train_data_3a_r1,
    formula = formula_others
    
  ),list(
    name = "3B 1to1 dup",
    model = cox_treated_3b,
    data = treated_train_data_3b,
    formula = formula_others)
  ,list(
    name = "3B relaxed dup",
    model = cox_treated_3b_r1,
    data = treated_train_data_3b_r1,
    formula = formula_others))
################################

source("bootstrap.R")
# BOOTSTRAP DISCRIMINATION INDICES
for (item in combined_list) {
  t_horizon<- 5
  if (grepl("PSM B", item$name)) { 
  #if (item$name == "1" ){
    print(paste("Validation:", item$name))
  
    discrim_res<-bootstrap_optimism(item$data, item$formula, item$name, n_bootstrap = 100)
    n_bootstrap <- 100
    #calib_res <- bootstrap_calibration_optimism(item$data, item$formula, item$name, n_bootstrap = n_bootstrap,t_horizon=t_horizon)
    }

}


# do 5-yr preds, range of death prob, HR table
for (item in combined_list) {
  #full_data <- read.csv("traindata/unmatched_train_data.csv", stringsAsFactors = FALSE)
  #predict.cox <- 1 - predictSurvProb(item$model,newdata=full_data,times=365*(1:5))
  ##column_name <- paste0("cox.5yr.", item$name)
  #full_data[[column_name]] <- ifelse(predict.cox[,5]==1,0.9999,predict.cox[,5])
  #column_name2 <- paste0("cox.5yr.", item$name,".cll")
  #full_data[[column_name2]] <- log(-log(1-full_data[[column_name]]))
  #write.csv(full_data, "traindata/unmatched_train_data.csv", row.names = FALSE)
  
  death_prob <- 1 - predictSurvProb(item$model, newdata = item$data, times = 365 * 5)
  range_death_prob <- range(death_prob)
  print(paste0("Model",item$name," range of predicted death prob:"))
  print(range_death_prob)
  
  #library(gtsummary)
  #library(gt)
 # library(webshot)
  #tbl <- tbl_regression(item$model, exponentiate = TRUE)
  
  # Convert to gt table
  #gt_tbl <- as_gt(tbl)
  
  # Save as PNG
  #gtsave(gt_tbl, paste0("HR_tables/Model",item$name, "_hr.png"))
}

# insample evaluation
source("harrell_Uno.R")
for (item in combined_list) {
  if (grepl("PSM", item$name)) {  
  ### UPDATE CONCORDANCE INSAMPLE METRICS ################################
  metrics_file <- "plots_and_metrics/insample_C_indices.csv"
  
    update_concordance(item$data, item$model, item$name, metrics_file)
  
  ### UPDATE CALIBRATION INSAMPLE METRICS ################################
  
  calib_file <- "plots_and_metrics/insample_calibration_metrics.csv"
  plot <- FALSE
  calibrate_pac(item$data, item$data,formula_others,calib_file,title_suffix = "")
  }
}

# Create an empty dataframe with required columns
res_C <- data.frame(
  Model = character(),
  Year = integer(),
  Harrell_C = numeric(),
  Harrell_C_2.5 = numeric(),
  Harrell_C_97.5 = numeric(),
  Uno_C = numeric(),
  Uno_C_2.5 = numeric(),
  Uno_C_97.5 = numeric(),
  Brier_Score = numeric(),
  Brier_R_squared = numeric(),
  stringsAsFactors = FALSE
)

# Save the empty dataframe to a CSV file
#write.csv(res_C, file = "plots_and_metrics/C_indices.csv", row.names = FALSE)

validate_df <- data.frame(
  Model = character(),
  Dxy = numeric(),
  R2 = numeric(),
  Slope = numeric(),
  D = numeric(),
  U = numeric(),
  Q = numeric(),
  G = numeric(),
  C_index = numeric(),
  stringsAsFactors = FALSE
)
#write.csv(validate_df, file = "plots_and_metrics/validate_rms.csv", row.names = FALSE)


