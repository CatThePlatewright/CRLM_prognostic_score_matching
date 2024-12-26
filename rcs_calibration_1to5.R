library(pacman)
pacman::p_load(survival,
               Hmisc,
               riskRegression)
library(rms)
library(survival)

calibrate_rcs <- function(tutorial, tutorial1, t_horizon, title_suffix = "", plot = FALSE) {
  
  # Read the existing OE table CSV file
  oe_table <- read.csv("plots_and_metrics/rcs_calibration_metrics.csv")
  
  # Transform back into time in years
  tutorial$time_to_event <- tutorial$time_to_event / 365.25
  tutorial1$time_to_event <- tutorial1$time_to_event / 365.25
  
  ################################################################################
  # Fit Cox PH model to model hazard of death. Compute OE.
  ################################################################################
  
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
  
  # Calculate mean Observed / Expected ratio OE
  # Observed risk of test set at t_horizon time
  obj <- summary(survfit(Surv(time_to_event, died) ~ 1, data = tutorial1), times = t_horizon)
  obs_t <- 1 - obj$surv # death probability
  
  # Predicted risk of test set
  tutorial1$pred <- riskRegression::predictRisk(cox1, newdata = tutorial1, times = t_horizon)
  
  # Expected
  exp_t <- mean(tutorial1$pred)
  
  # Calculate the OE values
  OE_t <- obs_t / exp_t
  # Calculate the confidence intervals for OE
  alpha <- 0.05
  OE_summary <- c( 
    "Model" = title_suffix,
    "Year" = t_horizon,
    "OE" = OE_t,
    "OE_2.5" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
    "OE_97.5" = OE_t * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  )
  
  ################################################################################
  # Calibration for predictions using RCS
  ################################################################################
  tutorial1$pred.cll <- log(-log(1 - tutorial1$pred))
  
  # Estimate actual risk
  vcal <- rms::cph(Surv(time_to_event, died) ~ rcs(pred.cll, 3),
                   x = T, y = T, surv = T, data = tutorial1) 
  
  dat_cal <- cbind.data.frame(
    "obs" = 1 - rms::survest(vcal, times = t_horizon, newdata = tutorial1)$surv,
    "lower" = 1 - rms::survest(vcal, times = t_horizon, newdata = tutorial1)$upper,
    "upper" = 1 - rms::survest(vcal, times = t_horizon, newdata = tutorial1)$lower,
    "pred" = tutorial1$pred
  )
  
  # Drop rows with NA values
  dat_cal <- na.omit(dat_cal)
  dat_cal <- dat_cal[order(dat_cal$pred), ]
  
  if (plot) {
    # Calibration plot
    par(xaxs = "i", yaxs = "i", las = 1)
    plot(
      dat_cal$pred, 
      dat_cal$obs,
      type = "l", 
      lty = 1, 
      xlim = c(0, 1),
      ylim = c(0, 1), 
      lwd = 2,
      xlab = "Predicted risk from developed model",
      ylab = "Predicted risk from refitted model", bty = "n"
    )
    lines(dat_cal$pred, 
          dat_cal$lower, 
          type = "l", 
          lty = 2, 
          lwd = 2)
    lines(dat_cal$pred, 
          dat_cal$upper,
          type = "l", 
          lty = 2, 
          lwd = 2)
    abline(0, 1, lwd = 2, lty = 2, col = 2)
    legend("bottomright",
           c("Ideal calibration",
             "Calibration curve based on secondary Cox model",
             "95% confidence interval"),
           col = c(2, 1, 1),
           lty = c(2, 1, 2),
           lwd = c(2, 2, 2),
           bty = "n",
           cex = 0.85)
    
    title(paste("Test data ", title_suffix, " Year ", t_horizon))
  }
  
  ################################################################################
  # ICI, E50, E90
  ################################################################################
  absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)
  
  numsum_cph <- c(
    "ICI" = mean(absdiff_cph),
    setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
    "Emax" = max(absdiff_cph)
  )
  
  ################################################################################
  # Calibration slope metrics
  ################################################################################
  
  tutorial1$lp <- predict(cox1, newdata = tutorial1)
  
  gval <- coxph(Surv(time_to_event, died) ~ lp, data = tutorial1)
  
  calslope_summary <- c(
    "calibration_slope" = gval$coef,
    "calslope_2.5"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
    "calslope_97.5" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
  )
  
  # Combine all the summaries into a data frame
  summary_row <- c(OE_summary, numsum_cph, calslope_summary)
  # Convert the vector to a data frame with one row
  summary_df <- as.data.frame(t(summary_row), stringsAsFactors = FALSE)
  # Append the new summary to the OE table
  oe_table <- rbind(oe_table, summary_df)
  
  # Write the updated table back to the CSV file
  write.csv(oe_table, "plots_and_metrics/rcs_calibration_metrics.csv", row.names = FALSE)
}


test_data <- read.csv("testdata/test_data.csv", stringsAsFactors = FALSE)
# Filter out untreated patients 
untreated_test_data <- subset(test_data, adjuvantchemo == 0)
# Filter out treated patients 
treated_test_data <- subset(test_data, adjuvantchemo == 1)

# load train data
train_data_m1 <- read.csv(file = "traindata/unmatched_train_data.csv", stringsAsFactors = FALSE)
untreated_train_data_m2 <- subset(train_data_m1, adjuvantchemo == 0)
treated_train_data_m2 <- subset(train_data_m1, adjuvantchemo == 1)

ratio <- "0.6"
train_data_m3 <- read.csv(file = paste("matchit_train.csv",sep = ""), stringsAsFactors = FALSE)
untreated_train_data_m3 <- subset(train_data_m3, adjuvantchemo == 0)
treated_train_data_m3 <- subset(train_data_m3, adjuvantchemo == 1)


models <- list(
  #train_data_m1,
  #untreated_train_data_m2,
  untreated_train_data_m3,
  #train_data_m1,
  #treated_train_data_m2,
  treated_train_data_m3
)

years <- 1:5
# Loop through each year
for (year in years) {
  cat("Year:", year, "\n")
  # Loop through each model 
  for (i in seq_along(models)) {
    train_data <- models[[i]]
    
    # Define the model label based on i
    #if (i <= 3) {
    #  model_label <- paste("Model", i, "A")
    #  test_data <-untreated_test_data
    #} else {
    #  model_label <- paste("Model", i - 3, "B")
    #  test_data <- treated_test_data
    #}
    if (i <= 1) {
      model_label <- paste("Model 3A MatchIT RF")
      test_data <-untreated_test_data
    } else {
      model_label <- paste("Model 3B MatchIT RF")
      test_data <- treated_test_data
    }
    # Save the plot as a PNG file
    #png(file = paste("plots_and_metrics/rcs_calibration_ratio/", model_label, "_year", year, ".png", sep = ""), width = 800, height = 600)
    calibrate_rcs(train_data, test_data, year, title_suffix = model_label)
    #dev.off()
    
  }
}

# Create an initial empty DataFrame for the OE table

oe_table <- data.frame(
  Model=character(),  
  Year = numeric(),
  OE = numeric(),
  OE_2.5 = numeric(),
  OE_97.5 = numeric(),
  ICI = numeric(),
  E50 = numeric(),
  E90 = numeric(),
  Emax = numeric(),
  calibration_slope = numeric(),
  calslope_2.5 = numeric(),
  calslope_97.5 = numeric()
)


# Write the empty oe_table to a CSV file
#write.csv(oe_table, "plots_and_metrics/rcs_calibration_metrics.csv", row.names = FALSE)


oe_table <- read.csv("plots_and_metrics/manual_bootstrap_calibration_1to5.csv")
library(dplyr)
library(stringr)
library(ggplot2)

# Filter the rows where "Model" column contains "A"
oe_table_A <- oe_table %>%
  filter(str_detect(Model, "A") | str_detect(Model, "1"))

# Filter the rows where "Model" column contains "B"
oe_table_B <- oe_table %>%
  filter(str_detect(Model, "B") | str_detect(Model, "1"))

# Plot for ICI for Models with "A"
ici_plot_A <- ggplot(oe_table_A, aes(x = Year, y = ICI, color = Model)) +
  geom_point() +
  geom_line() +
  labs(x = "Prediction time horizon (in years)", y = "ICI", color = "Model") +
  #ggtitle("ICI vs Prediction Year (untreated cohort)") +
  theme_bw() +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "black")

# Plot for ICI for Models with "B"
ici_plot_B <- ggplot(oe_table_B, aes(x = Year, y = ICI, color = Model)) +
  geom_point() +
  geom_line() +
  labs(x = "Prediction time horizon (in years)", y = "ICI", color = "Model") +
  #ggtitle("ICI vs Prediction Year (treated cohort)") +
  theme_bw() +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "black")

# Plot for OE ratio for Models with "A"
oe_plot_A <- ggplot(oe_table_A, aes(x = Year, y = OE_Ratio, color = Model)) +
  geom_point() +
  geom_line() +
  labs(x = "Prediction time horizon (in years)", y = "Observation-Estimation ratio", color = "Model") +
  #ggtitle("OE ratio vs Prediction Year (untreated cohort)") +
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black")

# Plot for OE ratio for Models with "B"
oe_plot_B <- ggplot(oe_table_B, aes(x = Year, y = OE_Ratio, color = Model)) +
  geom_point() +
  geom_line() +
  labs(x = "Prediction time horizon (in years)", y = "Observation-Estimation ratio", color = "Model") +
  #ggtitle("OE ratio vs Prediction Year (treated cohort)") +
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black")


# Create the plot
slope_plot <- ggplot(oe_table, aes(x = Year, y = calibration_slope.lp, color = Model)) +
  geom_point() +
  geom_line(aes(linetype = "calibration_slope")) +
  #geom_line(aes(y = calslope_2.5, linetype = "calslope_2.5"), data = oe_table) +
  #geom_line(aes(y = calslope_97.5, linetype = "calslope_97.5"), data = oe_table) +
  labs(x = "Year", y = "Calibration Slope", color = "Model") +
  ggtitle("Calibration Slope vs Prediction Year") +
  theme_bw() 
  #scale_linetype_manual(values = c("calibration_slope" = "solid", "calslope_2.5" = "dotted", "calslope_97.5" = "twodash"), name = "Line Legend")

# Save the plots
ggsave("plots_and_metrics/oe_plot_A.png", oe_plot_A, width = 8, height = 6, dpi = 300)
ggsave("plots_and_metrics/ici_plot_B.png", ici_plot_B, width = 8, height = 6, dpi = 300)
ggsave("plots_and_metrics/ici_plot_A.png", ici_plot_A, width = 8, height = 6, dpi = 300)
ggsave("plots_and_metrics/oe_plot_B.png", oe_plot_B, width = 8, height = 6, dpi = 300)

ggsave("plots_and_metrics/calslope_plot.png", slope_plot, width = 8, height = 6, dpi = 300)
