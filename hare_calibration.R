##### R CODE FOR CONSTRUCTING CALIBRATION CURVES AND NUMERICAL

##### METRICS OF CALIBRATION USING HR



#install.packages("polspline")

library(survival)

library(rms)

library(randomForestSRC)

library(pec)

library(polspline)
# Define a global dataframe to store the metrics
metrics_df <- data.frame(Model = character(),
                         ICI = numeric(),
                         E50 = numeric(),
                         E90 = numeric(), 
                         Calibration_Slope = numeric(),
                         CI_2.5_percentile_Calibration_Slope = numeric(),
                         CI_97.5_percentile_Calibration_Slope =numeric(),
                         stringsAsFactors = FALSE,
                         row.names = NULL)  # Avoid automatic index column


test_data_m1 <- read.csv("testdata/test_data_model1.csv", stringsAsFactors = FALSE)
# Filter out untreated patients 
untreated_test_data_m1 <- subset(test_data_m1, adjuvantchemo == 0)
# Filter out treated patients 
treated_test_data_m1 <- subset(test_data_m1, adjuvantchemo == 1)

# Repeat for other models...
test_data_m2<- read.csv("testdata/test_data_m2.csv", stringsAsFactors = FALSE)
untreated_test_data_m2 <- subset(test_data_m2, adjuvantchemo == 0)
treated_test_data_m2 <- subset(test_data_m2, adjuvantchemo == 1)

test_data_m3<- read.csv("testdata/test_data_m3.csv", stringsAsFactors = FALSE)
untreated_test_data_m3 <- subset(test_data_m3, adjuvantchemo == 0)
treated_test_data_m3 <- subset(test_data_m3, adjuvantchemo == 1)

calibrate_survival_probabilities <- function(data, num_years, cox_col, cox_cll_col, title_suffix = "",x_pos = 0.1, y_pos = 0.5) {
  tutorial <- data
  print(as.matrix(tutorial[[cox_cll_col]]))
  # Calibration for predictions of SURVIVAL probabilities
  calibrate.cox <- hare(data = tutorial$time_to_event, delta = tutorial$died,
                        cov = as.matrix(tutorial[[cox_cll_col]]))

  # A sequence of DEATH probabilities ranging from the 1st to the 99th percentile of tutorial$cox.5yr with 100 equally spaced values.
  predict.grid.cox <- seq(quantile(tutorial[[cox_col]], probs = 0.01),
                          quantile(tutorial[[cox_col]], probs = 0.99), length = 100)
  print(predict.grid.cox)
  # map survival probabilities to the log cumulative hazard scale
  predict.grid.cox.cll <- log(-log(1 - predict.grid.cox))
  # Calibrating: refine the predicted survival probabilities based on the calibration model.
  predict.calibrate.cox <- phare(num_years * 365, predict.grid.cox.cll, calibrate.cox)
 # print(paste("predict.calibrate.cox: ",predict.calibrate.cox)) print(paste("length tutorial[[cox_col]]: ",length(tutorial[[cox_col]]))) print(paste("length predict calib cox: ",length(predict.calibrate.cox)))
 
  
   # Plot
  png(paste("calibration_plot_", num_years, "year", title_suffix, ".png", sep = ""))
  plot(predict.grid.cox, predict.calibrate.cox, type = "l", lty = 1, col = "red",
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Predicted probability of survival",
       ylab = "Observed probability of survival")
  title(paste(num_years, "year", "mortality risk", title_suffix))
  par(new = TRUE)
  plot(density(tutorial[[cox_col]]), axes = FALSE, xlab = NA, ylab = NA, main = "")
  axis(side = 4)
  dev.off()
  #calibration slope
  alpha <- 0.05
  # Calculate calibration slope and confidence intervals
  calslope_summary <- c(
    "calibration slope" = calibrate.cox$coef,
    "2.5 %"  = calibrate.cox$coef - qnorm(1 - alpha / 2) * sqrt(calibrate.cox$var),
    "97.5 %" = calibrate.cox$coef + qnorm(1 - alpha / 2) * sqrt(calibrate.cox$var)
  )
  
  
  # Calculate the metrics
  ICI <- round(mean(abs(tutorial[[cox_col]] - predict.calibrate.cox)), 2)
  E50 <- round(median(abs(tutorial[[cox_col]] - predict.calibrate.cox)), 2)
  E90 <- round(quantile(abs(tutorial[[cox_col]] - predict.calibrate.cox), probs = 0.9), 2)
  
  # Add the metrics to the global dataframe
  metrics_df <<- rbind(metrics_df, data.frame(Model = paste(title_suffix),
                                                            ICI = ICI,
                                                            E50 = E50,
                                                            E90 = E90,
                       Calibration_Slope = calslope_summary[1],
                       CI_2.5_percentile_Calibration_Slope = calslope_summary[2],
                       CI_97.5_percentile_Calibration_Slope = calslope_summary[3]))
  return(metrics_df)
}

#calibrate_survival_probabilities(untreated_test_data_m1, 1, "cox.1yr", "cox.1yr.cll", " Model 1 (untreated)")
#calibrate_survival_probabilities(untreated_test_data_m1, 2, "cox.2yr", "cox.2yr.cll", " Model 1 (untreated)")
#calibrate_survival_probabilities(untreated_test_data_m1, 3, "cox.3yr", "cox.3yr.cll", " Model 1 (untreated)")
#calibrate_survival_probabilities(untreated_test_data_m1, 4, "cox.4yr", "cox.4yr.cll", " Model 1 (untreated)")
print(calibrate_survival_probabilities(untreated_test_data_m1, 5, "cox.5yr", "cox.5yr.cll", " Model 1 (untreated)"))
#calibrate_survival_probabilities(untreated_test_data_m2, 1, "cox.1yr", "cox.1yr.cll", " Model 2A (untreated)")
#calibrate_survival_probabilities(untreated_test_data_m2, 2, "cox.2yr", "cox.2yr.cll", " Model 2A (untreated)")
#calibrate_survival_probabilities(untreated_test_data_m2, 3, "cox.3yr", "cox.3yr.cll", " Model 2A (untreated)")
#calibrate_survival_probabilities(untreated_test_data_m2, 4, "cox.4yr", "cox.4yr.cll", " Model 2A (untreated)")
calibrate_survival_probabilities(untreated_test_data_m2, 5, "cox.5yr", "cox.5yr.cll", " Model 2A (untreated)")

calibrate_survival_probabilities(untreated_test_data_m3, 1, "cox.1yr", "cox.1yr.cll", " Model 3A (untreated)")
calibrate_survival_probabilities(untreated_test_data_m3, 2, "cox.2yr", "cox.2yr.cll", " Model 3A (untreated)")
calibrate_survival_probabilities(untreated_test_data_m3, 3, "cox.3yr", "cox.3yr.cll", " Model 3A (untreated)")
calibrate_survival_probabilities(untreated_test_data_m3, 4, "cox.4yr", "cox.4yr.cll", " Model 3A (untreated)")
calibrate_survival_probabilities(untreated_test_data_m3, 5, "cox.5yr", "cox.5yr.cll", " Model 3A (untreated)")

calibrate_survival_probabilities(treated_test_data_m1, 1, "cox.1yr", "cox.1yr.cll", " Model 1 (treated)")
calibrate_survival_probabilities(treated_test_data_m1, 2, "cox.2yr", "cox.2yr.cll", " Model 1 (treated)")
calibrate_survival_probabilities(treated_test_data_m1, 3, "cox.3yr", "cox.3yr.cll", " Model 1 (treated)")
calibrate_survival_probabilities(treated_test_data_m1, 4, "cox.4yr", "cox.4yr.cll", " Model 1 (treated)")
calibrate_survival_probabilities(treated_test_data_m1, 5, "cox.5yr", "cox.5yr.cll", " Model 1 (treated)")

calibrate_survival_probabilities(treated_test_data_m2, 1, "cox.1yr", "cox.1yr.cll", " Model 2B (treated)")
calibrate_survival_probabilities(treated_test_data_m2, 2, "cox.2yr", "cox.2yr.cll", " Model 2B (treated)")
calibrate_survival_probabilities(treated_test_data_m2, 3, "cox.3yr", "cox.3yr.cll", " Model 2B (treated)")
calibrate_survival_probabilities(treated_test_data_m2, 4, "cox.4yr", "cox.4yr.cll", " Model 2B (treated)")
calibrate_survival_probabilities(treated_test_data_m2, 5, "cox.5yr", "cox.5yr.cll", " Model 2B (treated)")

calibrate_survival_probabilities(treated_test_data_m3, 1, "cox.1yr", "cox.1yr.cll", " Model 3B (treated)")
calibrate_survival_probabilities(treated_test_data_m3, 2, "cox.2yr", "cox.2yr.cll", " Model 3B (treated)")
calibrate_survival_probabilities(treated_test_data_m3, 3, "cox.3yr", "cox.3yr.cll", " Model 3B (treated)")
calibrate_survival_probabilities(treated_test_data_m3, 4, "cox.4yr", "cox.4yr.cll", " Model 3B (treated)")
calibrate_survival_probabilities(treated_test_data_m3, 5, "cox.5yr", "cox.5yr.cll", " Model 3B (treated)")



write.csv(metrics_df, file = "metrics_df_hare.csv", row.names = FALSE)
