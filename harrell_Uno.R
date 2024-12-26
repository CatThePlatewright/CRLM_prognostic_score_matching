library(survival)
library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               riskRegression,
               timeROC)

# Define the function
update_concordance <- function(data, cox_model, model_name, res_C_file) {
  # Read the initial dataframe from a CSV file

  res_C <- read.csv(res_C_file, stringsAsFactors = FALSE)
   # Calculate Harrell's C-index
  data$lp <- predict(cox_model, newdata = data)
  harrell_C_model <- concordance(Surv(time_to_event, died) ~ lp, data, reverse = TRUE)
  # Uno's C
  Uno_C_model <- concordance(Surv(time_to_event, died) ~ lp, data, reverse = TRUE, timewt = "n/G2")
  
  # Confidence intervals
  alpha <- 0.05
  
  harrell_CI_lower <- ifelse(harrell_C_model$concordance - qnorm(1 - alpha/2) * sqrt(harrell_C_model$var) < 0 &
                               abs(harrell_C_model$concordance - qnorm(1 - alpha/2) * sqrt(harrell_C_model$var)) < 0.00001,
                             0,
                             harrell_C_model$concordance - qnorm(1 - alpha/2) * sqrt(harrell_C_model$var))
  Uno_CI_lower <- ifelse(Uno_C_model$concordance - qnorm(1 - alpha/2) * sqrt(Uno_C_model$var) < 0 &
                           abs(Uno_C_model$concordance - qnorm(1 - alpha/2) * sqrt(Uno_C_model$var)) < 0.00001,
                         0,
                         Uno_C_model$concordance - qnorm(1 - alpha/2) * sqrt(Uno_C_model$var))
  
  
    # brier score
  b1 <- brier(cox_model, times = 365 * 5, newdata = data, ties = TRUE, detail = FALSE, timefix = TRUE, efron = FALSE)
  brier_score <- b1$brier  # Extract the Brier score
  # Print the structure of the object
  str(b1)
  
  # Print the summary of the object
  summary(b1)
  
  # List all the components in the object
  names(b1)
  res_C <- rbind(res_C, data.frame(
    Model = model_name, # Modify if necessary
    Harrell_C = harrell_C_model$concordance,
    Harrell_C_2.5 = harrell_CI_lower,
    Harrell_C_97.5 = harrell_C_model$concordance + qnorm(1 - alpha/2) * sqrt(harrell_C_model$var),
    Uno_C = Uno_C_model$concordance,
    Uno_C_2.5 = Uno_CI_lower,
    Uno_C_97.5 = Uno_C_model$concordance + qnorm(1 - alpha/2) * sqrt(Uno_C_model$var),
    Brier_Score = b1$brier))
  
  write.csv(res_C, res_C_file, row.names = FALSE)
    
  # Return the updated results dataframe
  return(res_C)
}


