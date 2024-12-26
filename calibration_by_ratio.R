# Section 1: Load Libraries ----
library(pacman)
pacman::p_load(survival, Hmisc, riskRegression, rms, dplyr)

# Section 2: Define Functions ----

# Function to calculate calibration metrics
calibrate_pac <- function(tutorial, tutorial1, title_suffix = "") {
  # Transform back into time in years
  tutorial$time_to_event <- tutorial$time_to_event / 365.25
  tutorial1$time_to_event <- tutorial1$time_to_event / 365.25
  
  # Training data
  if (grepl("^Model 1 [AB]$", title_suffix)) {
    cox1 <- coxph(Surv(time_to_event, died) ~ age + T_stage + N_stage +
                    rightleft.Rec + cea_carcinoembryonic.antigen + extrahepatic.disease +
                    size + number_liver_mets + R0 + KRAS + adjuvantchemo, x = TRUE, data = tutorial)
  } else {
    cox1 <- coxph(Surv(time_to_event, died) ~ age + T_stage + N_stage +
                    rightleft.Rec + cea_carcinoembryonic.antigen + extrahepatic.disease +
                    size + number_liver_mets + R0 + KRAS, x = TRUE, data = tutorial)
  }
  
  # Observed / Expected ratio
  t_horizon <- 5
  
  # Observed test
  obj <- summary(survfit(Surv(time_to_event, died) ~ 1, data = tutorial1), times = t_horizon)
  obs_t <- 1 - obj$surv # Death probability
  
  # Predicted risk
  tutorial1$pred <- riskRegression::predictRisk(cox1, newdata = tutorial1, times = t_horizon)
  # Expected
  exp_t <- mean(tutorial1$pred)
  
  # Calculate the OE values
  OE_t <- obs_t / exp_t
  
  # Calculate the confidence intervals for OE
  alpha <- 0.05
  OE_summary <- c(
    "Model" = title_suffix,
    "OE" = OE_t,
    "OE_2.5" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
    "OE_97.5" = OE_t * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  )
  
  # Calibration plot data (skipped plotting)
  tutorial1$pred.cll <- log(-log(1 - tutorial1$pred))
  vcal <- rms::cph(Surv(time_to_event, died) ~ rcs(pred.cll, 3), x = TRUE, y = TRUE, surv = TRUE, data = tutorial1) 
  
  dat_cal <- cbind.data.frame(
    "obs" = 1 - rms::survest(vcal, times = 5, newdata = tutorial1)$surv,
    "lower" = 1 - rms::survest(vcal, times = 5, newdata = tutorial1)$upper,
    "upper" = 1 - rms::survest(vcal, times = 5, newdata = tutorial1)$lower,
    "pred" = tutorial1$pred
  )
  dat_cal <- na.omit(dat_cal)
  dat_cal <- dat_cal[order(dat_cal$pred), ]
  
  absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)
  
  numsum_cph <- c(
    "ICI" = mean(absdiff_cph),
    setNames(quantile(absdiff_cph, c(0.5, 0.9), na.rm = TRUE), c("E50", "E90")),
    "Emax" = max(absdiff_cph)
  )
  
  # Calibration slope (fixed time point)
  tutorial1$lp <- predict(cox1, newdata = tutorial1)
  gval <- coxph(Surv(time_to_event, died) ~ lp, data = tutorial1)
  
  calslope_summary <- c(
    "calibration_slope" = gval$coef,
    "calslope_2.5"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
    "calslope_97.5" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
  )
  
  # Combine all the summaries into a data frame
  summary_row <- c(OE_summary, numsum_cph, calslope_summary)
  
  # Return the summary as a data frame
  return(as.data.frame(t(summary_row), stringsAsFactors = FALSE))
}

# Section 3: Data Preparation ----

# Read the test data
test_data <- read.csv(file = "testdata/test_data.csv", stringsAsFactors = FALSE)
untreated_test_data <- subset(test_data, adjuvantchemo == 0)

# Sample file name
file <- "traindata/matched_train_model3.csv"

# Generate range of ratio values from 1.0 to 2.1 in increments of 0.1
ratio_range <- seq(0.3, 2.1, by = 0.1)

# Initialize an empty data frame to store results
all_results <- data.frame()

# Section 4: Loop Through Ratios ----

# Loop through each ratio
for (i in 1:length(ratio_range)) {
  ratio_number <- ratio_range[i]
  
  # Create file name with suffix
  file_suffix <- paste("_ratio", format(ratio_number, nsmall = 1, trim = TRUE), ".csv", sep = "")
  
  # Read the CSV file corresponding to the current ratio
  train_data <- read.csv(file = gsub(".csv", file_suffix, file), stringsAsFactors = FALSE)
  
  # Filter out untreated patients
  untreated_train_data <- subset(train_data, adjuvantchemo == 0)
  
  # Calculate calibration metrics for untreated data
  untreated_results <- calibrate_pac(untreated_train_data, untreated_test_data, paste("Untreated", ratio_number))
  all_results <- rbind(all_results, untreated_results)
}

# Convert Ratio and Treatment columns
all_results <- all_results %>%
  mutate(Ratio = as.numeric(gsub(".* ([0-9.]+)", "\\1", Model)))

# Ensure necessary columns are numeric
numeric_columns <- c("OE", "ICI", "E50", "E90", "Emax", "calibration_slope.lp")
all_results[numeric_columns] <- lapply(all_results[numeric_columns], as.numeric)

# Section 5: Plotting Metrics ----
# Set up a 3x2 plotting layout
par(mfrow = c(3, 2), mar = c(3, 3, 2, 1))

plot_individual_metric <- function(metric_name, y_label, model1_ref, model2_ref) {
  untreated_values <- all_results[[metric_name]]
  
  # Determine appropriate y-axis limits based on the range of data and reference values
  y_min <- min(c(untreated_values, model1_ref, model2_ref), na.rm = TRUE)
  y_max <- max(c(untreated_values, model1_ref, model2_ref), na.rm = TRUE)
  
  # Set default y-axis limits
  ylim_range <- c(y_min, y_max)
  
  # Plot for untreated data
  plot(ratio_range, untreated_values, type = 'o', main = paste(y_label), 
       xlab = 'Ratio', ylab = y_label, ylim = ylim_range, cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.8)
  
  # Add lines for Model 1 A and Model 2 A reference values
  lines(c(0, max(ratio_range)), c(model1_ref, model1_ref), col = "red")
  lines(c(0, max(ratio_range)), c(model2_ref, model2_ref), col = "green")
}

# Reference values for Model 1 A and Model 2 A
model1_ref <- c(1.1219, 0.0816, 0.0932, 0.1180, 0.1328, 0.7780)
model2_ref <- c(1.1807, 0.1088, 0.1274, 0.1483, 0.1500, 0.5138)

# Plotting the metrics
plot_individual_metric("OE", "Observed/Expected Ratio", model1_ref[1], model2_ref[1])
plot_individual_metric("ICI", "Integrated Calibration Index", model1_ref[2], model2_ref[2])
plot_individual_metric("E50", "E50", model1_ref[3], model2_ref[3])
plot_individual_metric("E90", "E90", model1_ref[4], model2_ref[4])
plot_individual_metric("Emax", "Emax", model1_ref[5], model2_ref[5])
plot_individual_metric("calibration_slope.lp", "Calibration Slope", model1_ref[6], model2_ref[6])
