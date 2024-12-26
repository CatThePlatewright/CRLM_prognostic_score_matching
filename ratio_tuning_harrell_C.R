library(survival)
library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)
validate_coxmodel <- function(model, n_sample=100)
{  set.seed(42)
validation_obj <- validate(model, B = n_sample)

# Extract relevant metrics from the validation object
dxy <- validation_obj['Dxy', 'index.corrected']
# Calculate the C-index
c_index <- (dxy / 2) + 0.5
return(c_index)
}
# Read the train data
file <- "traindata/matched_train_equalized_r"


# Generate range of ratio values from 0.3 to 2.1 in increments of 0.1
ratio_range <- seq(0.4, 4.0, by = 0.1)

# Create empty vectors to store C-index and Brier score values
harrell_C_untreated <- vector("numeric", length(ratio_range))
harrell_C_treated <- vector("numeric", length(ratio_range))


# Loop through each ratio and calculate performance metrics on the validation set
for (i in 1:length(ratio_range)) {
  ratio_number <- ratio_range[i]
  
  # Create file name with suffix
  file_name <- paste0(file, format(ratio_number, nsmall = 1, trim = TRUE), ".csv")
  train_data <- read.csv(file = file_name, stringsAsFactors = FALSE)
  
  # Split train data into training and validation sets
  #set.seed(123)
  #train_indices <- sample(nrow(train_data), 0.8 * nrow(train_data))
  #train_set <- train_data[train_indices, ]
  #val_set <- train_data[-train_indices, ]
  
  # Fit cox model on training set
  untreated_train_set <- subset(train_data, adjuvantchemo == 0)
  treated_train_set <- subset(train_data, adjuvantchemo == 1)
  
  cox_untreated <- cph(Surv(time_to_event, died) ~ age + T_stage + N_stage +
                               rightleft.Rec + cea_carcinoembryonic.antigen + extrahepatic.disease +
                               size + number_liver_mets + R0 + KRAS, x = TRUE, y=TRUE, surv=TRUE,data = untreated_train_set)
  
  cox_treated <- cph(Surv(time_to_event, died) ~ age + T_stage + N_stage +
                         rightleft.Rec + cea_carcinoembryonic.antigen + extrahepatic.disease +
                         size + number_liver_mets + R0 + KRAS, x = TRUE, y=TRUE, surv=TRUE,data = treated_train_set)
  # Calculate Harrell's C-index for validation set
  harrell_C_untreated[i] <- validate_coxmodel(cox_untreated)
  harrell_C_treated[i] <- validate_coxmodel(cox_treated)
}

print(harrell_C_untreated)
print(harrell_C_treated)
# Select the optimal ratio based on the performance metrics on the validation set
optimal_ratio_untreated <- ratio_range[which.max(harrell_C_untreated)]
optimal_ratio_treated <- ratio_range[which.max(harrell_C_treated)]
####### Plot harrell_C_treated ############################################################################
png(filename = "plots_and_metrics/bootstrap_C_index_vs_ratio.png", width = 800, height = 600)
plot(ratio_range, harrell_C_treated, type = "l", col = "blue", lwd = 2,
     xlab = "Ratio parameter r", ylab = "Bootstrapped Harrell's C index",
     main = "Bootstrapped Harrell's C index vs. ratio parameter")

# Add harrell_C_untreated to the same plot
lines(ratio_range, harrell_C_untreated, col = "red", lwd = 2)
abline(h = 0.607925002, col = "green", lwd = 2)  # Model 1 line
abline(h = 0.610707937, col = "orange", lwd = 2)  # Model 2A line
abline(h = 0.566546002, col = "cyan", lwd = 2)  # Model 2B line
# Close the PNG device to save the file
# Add a legend to the plot
legend("bottomright",legend = c("Model 1","Model 2A","Model 2B","Model 3A", "Model 3B"),
       col = c( "green","orange","cyan","red","blue"), lwd = 2, cex = 1.2)
dev.off()


print(paste("Optimal ratio A:", optimal_ratio_untreated, " B: ",optimal_ratio_treated))


# OLD: load the train data with optimal ratio
############################################################################
file_suffix <- paste("_ratio", format(optimal_ratio, nsmall = 1, trim = TRUE), ".csv", sep = "")
train_data <- read.csv(file = gsub(".csv", file_suffix, file), stringsAsFactors = FALSE)
# Fit the final model using the selected ratio on the entire train data
untreated_train_data <- subset(train_data, adjuvantchemo == 0)
res.cox_untreated_final <- coxph(Surv(time_to_event, died) ~ age + T_stage + N_stage +
                                   rightleft.Rec + cea_carcinoembryonic.antigen + extrahepatic.disease +
                                   size + number_liver_mets + R0 + KRAS, x = TRUE, data = untreated_train_data)


############################################################################
test_data <- read.csv(file = "testdata/test_data.csv", stringsAsFactors = FALSE)

# Evaluate the final model on the test data
untreated_test_data <- subset(test_data, adjuvantchemo == 0)
untreated_test_data$lp <- predict(res.cox_untreated_final, newdata = untreated_test_data)

# Calculate performance metrics on the test data
harrell_C_test <- concordance(Surv(time_to_event, died) ~ lp, untreated_test_data, reverse = TRUE)$concordance
# Calculate Harrell_C_2.5 and Harrell_C_97.5
harrell_CI_lower <- ifelse(harrell_C_test$concordance - qnorm(1 - alpha/2) * sqrt(harrell_C_test$var) < 0 &
                             abs(harrell_C_test$concordance - qnorm(1 - alpha/2) * sqrt(harrell_C_test$var)) < 0.00001,
                           0,
                           harrell_C_test$concordance - qnorm(1 - alpha/2) * sqrt(harrell_C_test$var))
harrell_C_2.5<- harrell_CI_lower
harrell_C_97.5 <- harrell_C_test$concordance + qnorm(1 - alpha/2) * sqrt(harrell_C_test$var)

Uno_CI_lower <- ifelse(Uno_C_test$concordance - qnorm(1 - alpha/2) * sqrt(Uno_C_test$var) < 0 &
                         abs(Uno_C_test$concordance - qnorm(1 - alpha/2) * sqrt(Uno_C_test$var)) < 0.00001,
                       0,
                       Uno_C_test$concordance - qnorm(1 - alpha/2) * sqrt(Uno_C_test$var))
Uno_C_2.5 <- Uno_CI_lower
Uno_C_97.5 <- Uno_C_test$concordance + qnorm(1 - alpha/2) * sqrt(Uno_C_test$var)
Brier_R_squared <- brier_test$rsquared
# Print the test scores
cat("Test scores:\n")
cat("Harrell's C-index on test data:", harrell_C_test, "\n")
cat("Uno's C-index on test data:", Uno_C_test, "\n")
cat("Brier score on test data:", brier_test$brier, "\n")




