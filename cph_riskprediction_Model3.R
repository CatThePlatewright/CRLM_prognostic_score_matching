# Load required libraries


create_train_data_with_risk <- function(train_data_filename, set_number, seed) {
  library(survival)
  library(rms)
  library(predtools)
  library(magrittr)
  library(prodlim)
  library(riskRegression)
  library(pec)
  library(dplyr)
  library(ggplot2)
  library(predtools)
  library(pacman)
  pacman::p_load(survival,
                 Hmisc,
                 riskRegression)
  # Read the training data
  train_data_m3 <- read.csv(file = train_data_filename, stringsAsFactors = FALSE)
  train_data_m3$size_binary <- ifelse(train_data_m3$size <= 5, 0, 1)
  
  
  # Filter out untreated patients for MODEL 3 -- matching method
  untreated_train_data <- subset(train_data_m3, adjuvantchemo == 0)
  treated_train_data <- subset(train_data_m3, adjuvantchemo == 1)
  
  # Define covariate sets
  covariates_set1 <- c("CEA_binary_20", "KRAS", "N_stage", "TBS_binary", "extrahepatic.disease")
  covariates_set2 <- c("size_binary", "KRAS", "N_stage")
  covariates_set3 <- c("size_binary", "number_binary", "DFI.12", "N_stage", "CEA_binary_200")
  
  # Choose the covariate set based on the input
  covariates <- switch(set_number,
                       "1" = covariates_set1,
                       "2" = covariates_set2,
                       "3" = covariates_set3,
                       "0" = c("age", "T_stage", "N_stage", "rightleft.Rec", "cea_carcinoembryonic.antigen", "extrahepatic.disease", 
                         "size", "number_liver_mets", "R0", "KRAS"))
  
  # Create the Cox regression formula
  formula <- as.formula(paste("Surv(time_to_event, died) ~", paste(covariates, collapse = "+")))
  
  # Print the formula for verification
  print(formula)
  
  # Fit the Cox proportional hazards model
  cox1 <- coxph(formula, x = TRUE, data = untreated_train_data)
  
  # modified 9/20: train on treated data only!!
  cox1 <- coxph(formula, x = TRUE, data = treated_train_data)
  
  # Predict the probability of survival for all patients
  
  death_prob <- 1 - predictSurvProb(cox1, newdata = train_data_m3, times = 365 * 5)
  train_data_m3$dying_prob_5y <- ifelse(death_prob == 1, 0.9999, death_prob)
  # Print the range of predicted death probabilities
  range_death_prob <- range(death_prob)
  print("range of predicted death prob:")
  print(range_death_prob)
  # Save the datasets with predictions into a file
  
  # modified 9/20
  output_filename <- paste0("traindata/unmatched_with_risk_preds_treated.csv")
  write.csv(train_data_m3, file = output_filename, row.names = FALSE)
  
  # Return the coefficients for further analysis if needed
  return(coefficients(cox1))
}
# Define the train data filenames and set number
train_data_filename <- "traindata/unmatched_train_data"
set_number <- "0"

# Loop through seeds 100 to 124

create_train_data_with_risk(paste0(train_data_filename,".csv"),set_number, seed)


 
