library(rms)
library(Hmisc)
library(rms)

validate_coxmodel <- function(model,model_name,n, n_sample=100)
{ set.seed(42)
  
  reps <- 100; 
  dxy <- numeric(reps) 
  c_index <- numeric(reps)  
  slope<- numeric(reps)
  
  for(i in 1 : reps) {
    g <- update(model, subset=sample(1 : n, n, replace=TRUE))
    v <- validate(g, B=n_sample)
    dxy[i] <- v['Dxy', 'index.corrected']
    c_index[i] <- (dxy[i] / 2) + 0.5
    slope[i] <- v['Slope', 'index.corrected']
    
  }
#validation_obj <- validate(model, B = n_sample)
#dxy <- validation_obj['Dxy', 'index.corrected']
#r2 <- validation_obj['R2', 'index.corrected']
#slope <- validation_obj['Slope', 'index.corrected']
#d <- validation_obj['D', 'index.corrected']
#u <- validation_obj['U', 'index.corrected']
#q <- validation_obj['Q', 'index.corrected']
#g <- validation_obj['g', 'index.corrected']

# Calculate the C-index
#c_index <- (dxy / 2) + 0.5

return (list(dxy=mean(dxy),c_index=mean(c_index),CI_c_index = quantile(c_index, c(.025, .975)),
             slope=mean(slope), CI_slope = quantile(slope, c(.025, .975)) ))
}
# Function to calculate Harrell's C-index
calculate_harrell_c_index <- function(model, data) {
  # Predict risk scores
  data$lp <- predict(model, newdata = data, type = "lp")
  
  # Calculate Harrell's C-index using the rcorr.cens function from the rms package
  #c_index <- rcorr.cens(data$lp, Surv(data$time_to_event, data$died))["C Index"]
  c_index2 <- concordance(Surv(time_to_event, died) ~ lp, data, reverse = TRUE)$concordance
  return(c_index2)

}
# Function to calculate Uno's C-index
calculate_uno_c_index <- function(model, data) {
  data$lp <- predict(model, newdata = data, type = "lp")
  
  # Calculate Uno's C-index
  uno_c_index <- concordance(Surv(time_to_event, died) ~ lp, data, reverse = TRUE, timewt = "n/G2")$concordance
  
  return(uno_c_index)
}
calculate_brier <- function(model, data, t_horizon=5) {
  # Calculate calibration slope
  brier <- brier(model, times = 365 * t_horizon, 
               newdata = data, ties = TRUE, 
               detail = FALSE, timefix = TRUE, efron = FALSE)$brier
  
  return(brier)
}

# Functions to calculate calibration metrics
calculate_oe_ratio <- function(model, data,t_horizon=5) {

  surv_obj <- summary(survfit(Surv(time_to_event, died) ~ 1, data = data), times = t_horizon)
  obs_t <- 1 - surv_obj$surv # Death probability
  data$pred <- tryCatch({
    riskRegression::predictRisk(model, newdata = data, times = t_horizon)
  }, error = function(e) {
    warning("Error in predictRisk: ", e$message)
    rep(NA, nrow(data)) # Return NA values if there's an error
  })
  
  if (any(is.na(data$pred))) {
    warning("There are NA values in the predicted risks. Removing rows with NA predictions.")
    data <- data[!is.na(data$pred), ]
    obs_t <- obs_t[!is.na(data$pred)]
  }
  
  exp_t <- mean(data$pred)
  OE_t <- obs_t / exp_t
  return(OE_t)
}

calculate_ici <- function(model, data,t_horizon=5) {
  # Calculate ICI
  data$pred <-  riskRegression::predictRisk(model, newdata = data, times = t_horizon)
  #data$pred <- ifelse(data$pred==1,0.9999,data$pred)
  
  # Define a small epsilon to avoid log(0) issues
  epsilon <- 1e-10  # A small number to prevent log(0) and very large negative logs
  # Cap values close to 0 and 1
  data$pred <- pmin(pmax(data$pred, epsilon), 1 - epsilon)
  data$pred.cll <- log(-log(1 - data$pred))
  if (sum(is.infinite(data$pred.cll)) > 0 ){
  stop("Infinite values found in the dataset. Remove or handle them.")
}

  vcal <- rms::cph(Surv(time_to_event, died) ~ pred.cll, x = TRUE, y = TRUE, surv = TRUE, data = data)
  #vcal <- rms::cph(Surv(time_to_event, died) ~ rcs(pred.cll, 3), x = TRUE, y = TRUE, surv = TRUE, data = data)
  dat_cal <- cbind.data.frame(
    obs = 1 - rms::survest(vcal, times = 5, newdata = data)$surv,
    pred = data$pred
  )
  absdiff <- abs(dat_cal$pred - dat_cal$obs)
  ici <- mean(absdiff)
  emax <- max(absdiff)
  e50 <- quantile(absdiff, 0.5) 
  e90 <- quantile(absdiff, 0.9) 
  return(list(ici=ici,emax=emax,e50=e50, e90=e90))
}

calculate_cal_slope <- function(model, data) {
  alpha <- 0.05
  # Calculate calibration slope
  data$lp <- predict(model, newdata = data)
  gval <- coxph(Surv(time_to_event, died) ~ lp, data = data)
  cal_slope <- gval$coef
  slopeLB <- gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var)
  slopeUB <- gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)

  return(list(slope=cal_slope,slopeLB=slopeLB,slopeUB=slopeUB))
}



# Function to perform bootstrap process with additional metrics
bootstrap_calibration_optimism <- function(data, model_formula, model_name, n_bootstrap = 200,t_horizon=5) {
  n <- nrow(data)
  # Transform back into time in years
  data$time_to_event <- data$time_to_event / 365.25
  # Fit the model on the original dataset
  original_model <- cph(model_formula, x = TRUE, y = TRUE, data = data, surv = TRUE)
  
  
  # Calculate the original metrics on the training dataset
  #original_cal_slope <- calculate_cal_slope(original_model, data)
  original_oe_ratio <- calculate_oe_ratio(original_model, data,t_horizon)
  ici_obj <-calculate_ici(original_model, data,t_horizon)
  original_ici <- ici_obj$ici
  original_emax <- ici_obj$emax
  original_e50 <- ici_obj$e50
  original_e90 <- ici_obj$e90
  
  # Initialize vectors to store optimism values
  #optimism_cal_slope <- numeric(n_bootstrap)
  optimism_oe_ratio <- numeric(n_bootstrap)
  optimism_ici <- numeric(n_bootstrap)
  optimism_emax <- numeric(n_bootstrap)
  optimism_e50 <- numeric(n_bootstrap)
  optimism_e90 <- numeric(n_bootstrap)
  
  # Perform bootstrap resampling
  for (i in 1:n_bootstrap) {
    print(paste0("Model", model_name," bootstrap calibration no. ",i))
    # Generate a bootstrap sample
    bootstrap_sample <- data[sample(1:n, replace = TRUE), ]
        # Fit the model on the bootstrap sample
    bootstrap_model <- cph(model_formula, x = TRUE, y = TRUE, data = bootstrap_sample, surv = TRUE)
    
    # Calculate metrics on the bootstrap sample
    #bootstrap_cal_slope <- calculate_cal_slope(bootstrap_model, bootstrap_sample)
    bootstrap_oe_ratio <- calculate_oe_ratio(bootstrap_model, bootstrap_sample,t_horizon)
    ici_bootstrap <-calculate_ici(bootstrap_model, bootstrap_sample,t_horizon)
    bootstrap_ici <- ici_bootstrap$ici
    bootstrap_emax <- ici_bootstrap$emax
    bootstrap_e50 <- ici_bootstrap$e50
    bootstrap_e90 <- ici_bootstrap$e90
    # Calculate metrics on the original data using the bootstrap model
    #test_cal_slope <- calculate_cal_slope(bootstrap_model, data)
    test_oe_ratio <- calculate_oe_ratio(bootstrap_model, data,t_horizon)
    ici_test_obj <-calculate_ici(bootstrap_model, data,t_horizon)
    test_ici <- ici_test_obj$ici
    test_emax <- ici_test_obj$emax
    test_e50 <- ici_test_obj$e50
    test_e90 <- ici_test_obj$e90
    
    # Calculate optimism: bootstrap metric - test metric
    #optimism_cal_slope[i] <- bootstrap_cal_slope - test_cal_slope
    optimism_oe_ratio[i] <- bootstrap_oe_ratio - test_oe_ratio
    optimism_ici[i] <- bootstrap_ici - test_ici
    optimism_emax[i] <- bootstrap_emax - test_emax
    optimism_e50[i] <- bootstrap_e50 - test_e50
    optimism_e90[i] <- bootstrap_e90 - test_e90
    
  }
  print(optimism_ici)
  # Calculate average optimism
  #average_optimism_cal_slope <- mean(optimism_cal_slope)
  average_optimism_oe_ratio <- mean(optimism_oe_ratio)
  average_optimism_ici <- mean(optimism_ici)
  average_optimism_emax <- mean(optimism_emax)
  average_optimism_e50 <- mean(optimism_e50)
  average_optimism_e90 <- mean(optimism_e90)
  
  # Overfitting-corrected metrics
  #corrected_cal_slope <- original_cal_slope - average_optimism_cal_slope
  corrected_oe_ratio <- original_oe_ratio - average_optimism_oe_ratio
  corrected_ici <- original_ici - average_optimism_ici
  corrected_emax <- original_emax - average_optimism_emax
  corrected_e50 <- original_e50 - average_optimism_e50
  corrected_e90 <- original_e90 - average_optimism_e90
  print(original_ici)
  print(corrected_ici)
  rms_validate <- validate_coxmodel(item$model,item$name,n)
  
  # Save results to a CSV file
  results_file <- "plots_and_metrics/new_calib.csv"
  existing_results_df <- read.csv(results_file)
  
  new_results_df <- data.frame(
    Model = model_name,
    OE_Ratio = corrected_oe_ratio,    
    ICI = corrected_ici,
    E50 = corrected_e50,
    E90 = corrected_e90,
    Emax = corrected_emax,
    Calibration_Slope_RMS = rms_validate$slope,
    Slope_.025 = rms_validate$CI_slope[1],
    Slope_.975 = rms_validate$CI_slope[2]
  )
  
  final_results_df <- rbind(existing_results_df, new_results_df)
  write.csv(final_results_df, results_file, row.names = FALSE)

  return(new_results_df)
}

# Function to perform bootstrap process
bootstrap_optimism <- function(data, model_formula,model_name,n_bootstrap = 100) {
  n <- nrow(data) 
  # Fit the model on the original dataset
  original_model <- cph(model_formula, x = TRUE, y = TRUE, data = data, surv = TRUE)
  # evaluate the model on the training dataset 
  original_harrell_c_index <- calculate_harrell_c_index(original_model, data)
  original_uno_c_index <- calculate_uno_c_index(original_model, data)
  original_brier <- calculate_brier(original_model,data)
  
  # Initialize vectors to store optimism values
  optimism_harrell <- numeric(n_bootstrap)
  optimism_uno <- numeric(n_bootstrap)
  optimism_brier <- numeric(n_bootstrap)
  # Perform bootstrap resampling
  for (i in 1:n_bootstrap) {
    # Generate a bootstrap sample
    print(paste0("Model", model_name," bootstrap no. ",i))
    
    
    bootstrap_sample <- data[sample(1:n, replace = TRUE), ] # n << nrow(data) !!!
    
    # Fit the model on the bootstrap sample
    bootstrap_model <- cph(model_formula, x = TRUE, y = TRUE, data = bootstrap_sample, surv = TRUE)
    # Calculate C-indexes on the bootstrap sample
    bootstrap_harrell_c_index <- calculate_harrell_c_index(bootstrap_model, bootstrap_sample)
    bootstrap_uno_c_index <- calculate_uno_c_index(bootstrap_model, bootstrap_sample)
    bootstrap_brier <- calculate_brier(bootstrap_model,bootstrap_sample)
    
    # Calculate C-indexes on the original data using the bootstrap model
    test_harrell_c_index <- calculate_harrell_c_index(bootstrap_model, data)
    test_uno_c_index <- calculate_uno_c_index(bootstrap_model, data)
    test_brier <- calculate_brier(bootstrap_model,data)
    
    # Calculate optimism: bootstrap C-index - test C-index
    optimism_harrell[i] <- bootstrap_harrell_c_index - test_harrell_c_index
    optimism_uno[i] <- bootstrap_uno_c_index - test_uno_c_index
    optimism_brier[i] <- bootstrap_brier - test_brier
  }
  
  # Calculate average optimism
  average_optimism_harrell <- mean(optimism_harrell)
  average_optimism_uno <- mean(optimism_uno)
  average_optimism_brier <- mean(optimism_brier)
  # Overfitting-corrected C-index
  corrected_harrell_c_index <- original_harrell_c_index - average_optimism_harrell
  corrected_uno_c_index <- original_uno_c_index - average_optimism_uno
  uno_c_LB <-original_uno_c_index - quantile(optimism_uno, c(.025, .975))[1]
  uno_c_UB <-original_uno_c_index - quantile(optimism_uno, c(.025, .975))[2]

  corrected_brier <- original_brier - average_optimism_brier
  brier_LB <-original_brier - quantile(optimism_brier, c(.025, .975))[1]
  brier_UB <-original_brier - quantile(optimism_brier, c(.025, .975))[2]
  
  rms_validate <- validate_coxmodel(item$model,item$name, n)
  
    results_file <- "plots_and_metrics/new_discrimination.csv"
  existing_results_df <- read.csv(results_file)
  
  #  Combine results into a dataframe
  new_results_df <- data.frame( Model = model_name,
                                Dxy_RMS = rms_validate$dxy,
                               Harrell_C_RMS = rms_validate$c_index,
                                Harrell_C_.025 = rms_validate$CI_c_index[1],
                                Harrell_C_.975 = rms_validate$CI_c_index[2],
                                Uno_C = corrected_uno_c_index,
                                Uno_C_.025 = uno_c_UB,
                                Uno_C_.975 = uno_c_LB,
                                Brier_Score = corrected_brier,
                               Brier_.025 = brier_UB,
                               Brier_.975 = brier_LB)
  
  
  # Step 3: Append the new results to the existing dataframe
  updated_results_df <- rbind(existing_results_df, new_results_df)
  write.csv(updated_results_df, results_file, row.names = FALSE)
  # Return the results dataframe
  return(new_results_df)
}



# Define the column names that match the output of the bootstrap_optimism function

#
column_names <- c("Model", "OE_Ratio","ICI", "E50","E90", "Emax", "Calibration_Slope_RMS",
                  "Slope_.025","Slope_.975")

# Initialize an empty dataframe with the defined column names
results_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(results_df) <- column_names  # Set column names explicitly
#write.csv(results_df, "plots_and_metrics/new_calib.csv", row.names = FALSE)

#column_names <- c("Model", "Dxy_RMS","Harrell_C_RMS","Harrell_C_CI","Uno_C","Brier_Score")


