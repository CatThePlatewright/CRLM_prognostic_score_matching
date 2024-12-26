library(pacman)
pacman::p_load(survival,
               Hmisc,
               riskRegression)
library(rms)
library(survival)

calibrate_pac <- function(tutorial,tutorial1,title_suffix = "") {
  
  # transform back into time in years
  tutorial$time_to_event <- tutorial$time_to_event   / 365.25
  tutorial1$time_to_event <- tutorial1$time_to_event   / 365.25
  
  # training data 
  if (grepl("^Model 1 [AB]$", title_suffix)) {
    cox1 <- coxph(Surv(time_to_event, died) ~ age + T_stage + N_stage +
                    rightleft.Rec + cea_carcinoembryonic.antigen + extrahepatic.disease +
                    size + number_liver_mets + R0 + KRAS + adjuvantchemo, x = TRUE, data = tutorial)
  } else {
    cox1 <- coxph(Surv(time_to_event, died) ~ age + T_stage + N_stage +
                    rightleft.Rec + cea_carcinoembryonic.antigen + extrahepatic.disease +
                    size + number_liver_mets + R0 + KRAS, x = TRUE, data = tutorial)
  }# Observed / Expected ratio
  t_horizon <- 5
  
  # Observed test
  obj <- summary(survfit(
    Surv(time_to_event, died) ~ 1, 
    data = tutorial1),
    times = t_horizon)
  
  obs_t <- 1 - obj$surv # death probability
  
  # Predicted risk 
  tutorial1$pred <- riskRegression::predictRisk(cox1, 
                                                newdata = tutorial1,
                                                times = t_horizon)
  # Expected
  exp_t <- mean(tutorial1$pred)
  
  # Read the existing OE table CSV file
  oe_table <- read.csv("plots_and_metrics/rcs_calibration_metrics_5y.csv")
  
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
  
  
  # Calibration plot ----------------------------------
  tutorial1$pred.cll <- log(-log(1 - tutorial1$pred))
  
  
  # Estimate actual risk
  vcal <- rms::cph(Surv(time_to_event, died) ~ rcs(pred.cll, 3),
                   x = T,
                   y = T,
                   surv = T,
                   data = tutorial1
  ) 
  
  dat_cal <- cbind.data.frame(
    "obs" = 1 - rms::survest(vcal, 
                             times = 5, 
                             newdata = tutorial1)$surv,
    
    "lower" = 1 - rms::survest(vcal, 
                               times = 5, 
                               newdata = tutorial1)$upper,
    
    "upper" = 1 - rms::survest(vcal, 
                               times = 5, 
                               newdata = tutorial1)$lower,
    "pred" = tutorial1$pred
  )
  
  dat_cal <- dat_cal[order(dat_cal$pred), ]
  
  #dev.new()
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
  
  title(paste("Test data ", title_suffix))
  
  # Numerical measures
  absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)
  
  numsum_cph <- c(
    "ICI" = mean(absdiff_cph),
    setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
    "Emax" = max(absdiff_cph)
  )
  
  # calibration slope (fixed time point)-------------------------------------
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
  write.csv(oe_table, "plots_and_metrics/rcs_calibration_metrics_5y.csv", row.names = FALSE)
}


test_data <- read.csv("testdata/test_survival_tree_preds.csv", stringsAsFactors = FALSE)
# Filter out untreated patients 
untreated_test_data <- subset(test_data, adjuvantchemo == 0)
# Filter out treated patients 
treated_test_data <- subset(test_data, adjuvantchemo == 1)


models <- list(
  ost_model1, ost_model2, ost_model3
)
# Loop through each model
for (i in seq_along(models)) {
  # Define the model label based on i
  if (i <= 3) {
    model_label <- paste("Model", i, "A")
    test_data <- untreated_test_data
  } else {
    model_label <- paste("Model", i - 3, "B")
    test_data <- treated_test_data
  }
  
  # Save the plot as a PNG file
  png(file = paste("rcs_calibration_plot_5y_", model_label, ".png", sep = ""), width = 800, height = 600)
  calibrate_pac(train_data, test_data,title_suffix = model_label)
  dev.off()
  
}


# Create an initial empty DataFrame for the OE table

oe_table <- data.frame(
  Model=character(),
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
#write.csv(oe_table, "plots_and_metrics/rcs_calibration_metrics_5y.csv", row.names = FALSE)



