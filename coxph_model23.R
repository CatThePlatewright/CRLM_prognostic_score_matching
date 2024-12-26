source("harrell_Uno.R")
#Cox multivariate:
#Hosmer Lemeshow statistical test
#install.packages("survival")
library(survival)
library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)

train_filename <- "traindata/unmatched_train_data_seed_100.csv"# "traindata/unmatched_train_data.csv"
test_filename <- "testdata/test_data_seed_100.csv" #"testdata/test_data.csv"
train_data <- read.csv(file =train_filename, stringsAsFactors = FALSE)
## CHANGE FILENAME for update_concordance() line 
test_data <- read.csv(file = test_filename, stringsAsFactors = FALSE)

untreated_train_data2 <- subset(train_data, adjuvantchemo == 0)
treated_train_data2 <- subset(train_data, adjuvantchemo == 1)

### LOAD TEST DATA AND AUGMENT #################################
test_data <- read.csv(file = test_filename, stringsAsFactors = FALSE)
untreated_test_data <- subset(test_data, adjuvantchemo == 0)
treated_test_data <- subset(test_data, adjuvantchemo == 1)

set_number <-0
# cox multivariate MODEL 3 
train_data_m3a <- read.csv(paste0("traindata/matched_train_set0_seed100.csv"), stringsAsFactors = FALSE)
untreated_train_data <- subset(train_data_m3a, adjuvantchemo == 0)
train_data_m3a_r1 <- read.csv(paste0("traindata/matched_train_set0_r1_seed100.csv"), stringsAsFactors = FALSE)
untreated_train_data_r1 <- subset(train_data_m3a_r1, adjuvantchemo == 0)
train_data_m3b <- read.csv(paste0("traindata/matched_train_set0_seed100.csv"), stringsAsFactors = FALSE)
treated_train_data <- subset(train_data_m3b, adjuvantchemo == 1)

# train on subset 1 of covariates (8/15)
covariates_set1 <- c("CEA_binary_20","KRAS","N_stage", "TBS_binary","extrahepatic.disease")
covariates_set2 <- c("size_binary","KRAS","N_stage")
covariates_set3 <- c("size_binary","number_binary","DFI.12", "N_stage", "CEA_binary_200")
# Choose the covariate set based on the input
covariates <- if (set_number == 1) {
  covariates_set1
} else if (set_number == 2) {
  covariates_set2
} else if (set_number == 3) {
  covariates_set3
} else {
  c("age", "T_stage", "N_stage", "rightleft.Rec",  "cea_carcinoembryonic.antigen", "extrahepatic.disease", 
    "size", "number_liver_mets", "R0", "KRAS")
}

formula <- as.formula(paste("Surv(time_to_event, died) ~", paste(covariates, collapse = "+")))

print(formula)

cox_untreated2 <- coxph(formula, x=TRUE, data=untreated_train_data2) #Model 2A
cox_treated2 <- coxph(formula, x=TRUE, data=treated_train_data2)#Model 2B
cox_untreated <- coxph(formula, x=TRUE, data=untreated_train_data) #Model 3A r=0.6
cox_untreated_r1 <- coxph(formula, x=TRUE, data=untreated_train_data_r1) #Model 3A r=1
cox_treated <- coxph(formula, x=TRUE, data=treated_train_data) #Model 3B

# 4. predict 5-year death probabilities for test patients
predict.cox2A <- 1 - predictSurvProb(cox_untreated2,newdata=untreated_test_data,times=365*(1:5))
predict.cox2B <- 1 - predictSurvProb(cox_treated2,newdata=treated_test_data,times=365*(1:5))

predict.cox3Ar0.6 <- 1 - predictSurvProb(cox_untreated,newdata=untreated_test_data,times=365*(1:5))
predict.cox3Ar1 <- 1 - predictSurvProb(cox_untreated_r1,newdata=untreated_test_data,times=365*(1:5))
predict.cox3B <- 1 - predictSurvProb(cox_treated,newdata=treated_test_data,times=365*(1:5))


# Predicted probability of death within 1,2,3,4,5 years.

untreated_test_data$cox.5yr.2A <- predict.cox2A[,5]
untreated_test_data$cox.5yr.2A <- ifelse(untreated_test_data$cox.5yr.2A==1,0.9999,untreated_test_data$cox.5yr.2A)
untreated_test_data$cox.5yr.2A.cll <- log(-log(1-untreated_test_data$cox.5yr.2A))
untreated_test_data$cox.5yr.3Ar0.6 <- predict.cox3Ar0.6[,5]
untreated_test_data$cox.5yr.3Ar0.6 <- ifelse(untreated_test_data$cox.5yr.3Ar0.6==1,0.9999,untreated_test_data$cox.5yr.3Ar0.6)
untreated_test_data$cox.5yr.3Ar0.6.cll <- log(-log(1-untreated_test_data$cox.5yr.3Ar0.6))
untreated_test_data$cox.5yr.3Ar1 <- predict.cox3Ar1[,5]
untreated_test_data$cox.5yr.3Ar1 <- ifelse(untreated_test_data$cox.5yr.3Ar1==1,0.9999,untreated_test_data$cox.5yr.3Ar1)
untreated_test_data$cox.5yr.3Ar1.cll <- log(-log(1-untreated_test_data$cox.5yr.3Ar1))

treated_test_data$cox.5yr.2B <- predict.cox2B[,5]
treated_test_data$cox.5yr.2B <- ifelse(treated_test_data$cox.5yr.2B==1,0.9999,treated_test_data$cox.5yr.2B)
treated_test_data$cox.5yr.2B.cll <- log(-log(1-treated_test_data$cox.5yr.2B))
treated_test_data$cox.5yr.3B <- predict.cox3B[,5]
treated_test_data$cox.5yr.3B <- ifelse(treated_test_data$cox.5yr.3B==1,0.9999,treated_test_data$cox.5yr.3B)
treated_test_data$cox.5yr.3B.cll <- log(-log(1-treated_test_data$cox.5yr.3B))
#write.csv(untreated_test_data, "testdata/untreated_test_5yr_preds.csv", row.names = FALSE)
#write.csv(treated_test_data, "testdata/treated_test_5yr_preds.csv", row.names = FALSE)
library(gtsummary)
library(gt)
library(webshot)
# Save each table as a PNG image
#gtsave(as_gt(tbl_regression(cox_untreated2, exponentiate = TRUE)), "HR_tables/Model2A.png")
#gtsave(as_gt(tbl_regression(cox_treated2, exponentiate = TRUE)), "HR_tables/Model2B.png")
#gtsave(as_gt(tbl_regression(cox_untreated, exponentiate = TRUE)), "HR_tables/Model3A_r0.6.png")
#gtsave(as_gt(tbl_regression(cox_treated, exponentiate = TRUE)), "HR_tables/Model3B.png")
#gtsave(as_gt(tbl_regression(cox_untreated_r1, exponentiate = TRUE)), "HR_tables/Model3A_r1.png")


file_name = "plots_and_metrics_subset/new_C_indices.csv"

update_concordance(untreated_test_data,cox_untreated2,paste0("2A set",set_number),file_name)
update_concordance(untreated_test_data,cox_untreated,paste0("3A set",set_number," 1to1"),file_name)
update_concordance(untreated_test_data,cox_untreated_r1,paste0("3A set",set_number," r=1.0"),file_name)

update_concordance(treated_test_data,cox_treated2,paste0("2B set",set_number),file_name)
update_concordance(treated_test_data,cox_treated,paste0("3B set",set_number),file_name)
