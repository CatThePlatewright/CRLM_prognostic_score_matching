library(pacman)
pacman::p_load(survival,
               Hmisc,
               riskRegression)
library(rms)
library(survival)

calibrate_pac <- function(tutorial, tutorial1, formula,file_path,title_suffix = "", set_number = 0, plot = FALSE) {
  # Read the CSV file
  oe_table <- read.csv(file_path)
  
  # Print the file path and covariate set to train on for confirmation
  print(paste("Training:", title_suffix))
  print(paste("training formula: ", formula))
  
  # Transform back into time in years
  tutorial$time_to_event <- tutorial$time_to_event / 365.25
  tutorial1$time_to_event <- tutorial1$time_to_event / 365.25
  
  cox1 <- coxph(formula, x = TRUE, data = tutorial)
  
  
  # Observed / Expected ratio
  t_horizon <- 5
  ################################
  ################################
  
  # Observed test
  obj <- summary(survfit(
    Surv(time_to_event, died) ~ 1, 
    data = tutorial1),
    times = t_horizon)
  
  obs_t <- 1 - obj$surv # Death probability
  
  
  # Predicted risk 
  tutorial1$pred <- riskRegression::predictRisk(cox1, newdata = tutorial1, times = t_horizon)
  tutorial1$pred <- tutorial1$pred <- ifelse(tutorial1$pred==1,0.9999,tutorial1$pred)
  
  # Check if there are any NA values in the predictions
  if (any(is.na(tutorial1$pred))) {
    warning("There are NA values in the predicted risks. Removing rows with NA predictions.")
    tutorial1 <- tutorial1[!is.na(tutorial1$pred), ]
    obs_t <- obs_t[!is.na(tutorial1$pred)]
  }
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
  
  ################################################################
  ################################################################
  
  tutorial1$pred.cll <- log(-log(1 - tutorial1$pred))
  
  # Estimate actual risk
  vcal <- rms::cph(Surv(time_to_event, died) ~ rcs(pred.cll, 3),
                   x = TRUE,
                   y = TRUE,
                   surv = TRUE,
                   data = tutorial1
  ) 
  
  dat_cal <- cbind.data.frame(
    "obs" = 1 - rms::survest(vcal, times = 5, newdata = tutorial1)$surv,
    "lower" = 1 - rms::survest(vcal, times = 5, newdata = tutorial1)$upper,
    "upper" = 1 - rms::survest(vcal, times = 5, newdata = tutorial1)$lower,
    "pred" = tutorial1$pred
  )
  # Drop rows with NA values
  dat_cal <- na.omit(dat_cal)
  dat_cal <- dat_cal[order(dat_cal$pred), ]
  # Calibration plot (only if plot argument is TRUE) ----------------------------------
  if (plot) {    
    # Save the plot as a PNG file
    png(file = paste("plots_and_metrics/calibration_plots/", model_label, ".png", sep = ""), width = 800, height = 600)
    
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
    dev.off()
  }
  
  # Numerical measures
  absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)
  
  numsum_cph <- c(
    "ICI" = mean(absdiff_cph),
    setNames(quantile(absdiff_cph, c(0.5, 0.9), na.rm = TRUE), c("E50", "E90")),
    "Emax" = max(absdiff_cph)
  )
  ################################################################
  ################################################################
  
  # Calibration slope (fixed time point)-------------------------------------
  tutorial1$lp <- predict(cox1, newdata = tutorial1)
  
  gval <- coxph(Surv(time_to_event, died) ~ lp, data = tutorial1)
  
  calslope_summary <- c(
    "calibration_slope" = gval$coef,
    "calslope_2.5" = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
    "calslope_97.5" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
  )
  
  # Combine all the summaries into a data frame
  summary_row <- c(OE_summary, numsum_cph, calslope_summary)
  # Convert the vector to a data frame with one row
  summary_df <- as.data.frame(t(summary_row), stringsAsFactors = FALSE)
  
  # Append the new summary to the OE table
  oe_table <- rbind(oe_table, summary_df)
  
  # Write the updated table back to the CSV file
  write.csv(oe_table, file_path, row.names = FALSE)
}




