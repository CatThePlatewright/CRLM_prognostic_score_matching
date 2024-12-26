#############################################################

library(dplyr)

calculate_metrics <- function(metrics_file) {
  # Load the metrics data
  metrics_data <- read.csv(metrics_file, stringsAsFactors = FALSE)
  
  # Initialize a list to store the results
  results_list <- list()
  
  # Define the sets
  sets <- c("1A", "2A", "3A", "1B", "2B", "3B","3B r1")

  # Compute statistics for each set
  for (set in sets) {
    # Filter the data for the current set
    set_data <- metrics_data[grep(paste0("^", set, " "), metrics_data$Model),]
    
    # Calculate statistics for Harrell's C
    harrell_c_mean <- mean(set_data$Harrell_C, na.rm = TRUE)
    harrell_c_median <- median(set_data$Harrell_C, na.rm = TRUE)
    harrell_c_min <- min(set_data$Harrell_C, na.rm = TRUE)
    harrell_c_max <- max(set_data$Harrell_C, na.rm = TRUE)
    
    # Calculate statistics for Uno's C
    uno_c_mean <- mean(set_data$Uno_C, na.rm = TRUE)
    uno_c_median <- median(set_data$Uno_C, na.rm = TRUE)
    uno_c_min <- min(set_data$Uno_C, na.rm = TRUE)
    uno_c_max <- max(set_data$Uno_C, na.rm = TRUE)
    
    # Calculate statistics for Brier Score
    brier_score_mean <- mean(set_data$Brier_Score, na.rm = TRUE)
    brier_score_median <- median(set_data$Brier_Score, na.rm = TRUE)
    brier_score_min <- min(set_data$Brier_Score, na.rm = TRUE)
    brier_score_max <- max(set_data$Brier_Score, na.rm = TRUE)
    
    # Store the results in the list
    results_list[[set]] <- data.frame(
      Set = set,
      HarrellC_Mean = harrell_c_mean,
      HarrellC_Median = harrell_c_median,
      HarrellC_Min = harrell_c_min,
      HarrellC_Max = harrell_c_max,
      UnoC_Mean = uno_c_mean,
      UnoC_Median = uno_c_median,
      UnoC_Min = uno_c_min,
      UnoC_Max = uno_c_max,
      BrierScore_Mean = brier_score_mean,
      BrierScore_Median = brier_score_median,
      BrierScore_Min = brier_score_min,
      BrierScore_Max = brier_score_max
    )
  }
  
  # Combine the results into a single data frame
  results <- do.call(rbind, results_list)
  
  return(results)
}

# Call the function and save the results
metrics_file <- "plots_and_metrics/new_C_indices_experiments.csv"
metrics_results1 <- calculate_metrics(metrics_file)

# View the results
print(metrics_results1)
