#Cox multivariate:
# Install and load the Hmisc packageinstall.packages("Hmisc")
source("harrell_Uno.R")
library(Hmisc)
library(survival)
library(rms)
library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)
library(pacman)
library(survAUC)
library(caret)

pacman::p_load(
  rio,
  survival,
  rms,
  pec,
  survminer,
  riskRegression,
  timeROC,
  plotrix,
  splines,
  knitr,
  kableExtra,
  gtsummary,
  boot,
  tidyverse,
  rsample,
  gridExtra,
  dplyr,
  webshot
)

#### MODEL 1 #############
train_filename <- "traindata/unmatched_train_data_seed_100.csv"# "traindata/unmatched_train_data.csv"
test_filename <- "testdata/test_data_seed_100.csv" #"testdata/test_data.csv"
train_data <- read.csv(file =train_filename, stringsAsFactors = FALSE)
## CHANGE FILENAME for update_concordance() line 
test_data <- read.csv(file = test_filename, stringsAsFactors = FALSE)

train_data_m1<- train_data
test_data_m1 <- test_data

untreated_test_data <- subset(test_data_m1, adjuvantchemo == 0)
treated_test_data <- subset(test_data_m1, adjuvantchemo == 1)

set_number <-0 # if 0 then default full set of 10 variables
# cox model fitted on subsets of covariates 8/14
covariates_set1 <- c("CEA_binary_20","KRAS","N_stage", "TBS_binary","extrahepatic.disease")
covariates_set2 <- c("size_binary","KRAS","N_stage")
covariates_set3 <- c("size_binary","number_binary","DFI.12", "N_stage", "CEA_binary_200")

#set 3) node-positive primary, disease-free interval from primary to metastases <12 months, number of hepatic tumors >1, largest hepatic tumor >5 cm and carcinoembryonic antigen level >200 ng/ml
covariates_set3 <- c("size_binary","number_binary","DFI.12", "N_stage", "CEA_binary_200")
covariates <- if (set_number == 1) {
  covariates_set1
} else if (set_number == 2) {
  covariates_set2
} else if (set_number == 3) {
  covariates_set3
} else if (set_number == 0)  {
  c("age", "T_stage", "N_stage", "rightleft.Rec",  "cea_carcinoembryonic.antigen", "extrahepatic.disease", 
    "size", "number_liver_mets", "R0", "KRAS")
}

############ ASIDE: UNIVARIABLE TESTS TO DO MULTIVARIATE ################################
for (var in covariates) {
  # Dynamically create the formula
  formula <- as.formula(paste("Surv(time_to_event, died) ~", var))
  
  # Run the Cox proportional hazards model
  res.cox <- coxph(formula, data = train_data_m1)
  
  # Print the model results
  print(var)
  
  # Print the summary of the model
  print(summary(res.cox))
}

covariates <- c("age","extrahepatic.disease", "number_liver_mets", "R0")
#############################################################
formula <- as.formula(paste("Surv(time_to_event, died) ~", paste(c(covariates, "adjuvantchemo"), collapse = "+")))
cox1 <- coxph(formula, x=TRUE, data=train_data_m1)

predict.cox1A <- 1 - predictSurvProb(cox1,newdata=untreated_test_data,times=365*(1:5))
predict.cox1B <- 1 - predictSurvProb(cox1,newdata=treated_test_data,times=365*(1:5))

untreated_test_data$cox.5yr.1A <- predict.cox1A[,5]
untreated_test_data$cox.5yr.1A <- ifelse(untreated_test_data$cox.5yr.1A==1,0.9999,untreated_test_data$cox.5yr.1A)
untreated_test_data$cox.5yr.1A.cll <- log(-log(1-untreated_test_data$cox.5yr.1A))
treated_test_data$cox.5yr.1B <- predict.cox1B[,5]
treated_test_data$cox.5yr.1B <- ifelse(treated_test_data$cox.5yr.1B==1,0.9999,treated_test_data$cox.5yr.1B)
treated_test_data$cox.5yr.1B.cll <- log(-log(1-treated_test_data$cox.5yr.1B))

res_C = update_concordance(untreated_test_data,cox1,paste0("1A set",set_number),"plots_and_metrics_subset/new_C_indices.csv")
res_C = update_concordance(treated_test_data,cox1,paste0("1B set",set_number),"plots_and_metrics_subset/new_C_indices.csv")


#write.csv(untreated_test_data, "testdata/untreated_test_5yr_preds.csv", row.names = FALSE)
#write.csv(treated_test_data, "testdata/treated_test_5yr_preds.csv", row.names = FALSE)
