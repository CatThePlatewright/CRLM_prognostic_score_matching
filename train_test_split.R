# Read the full imputed data
imputed_data <- read.csv(file = "traindata_main/CRLM_imputed.csv", stringsAsFactors = FALSE)
imputed_data$died <- ifelse(imputed_data$died, 1, 0) # Convert logical to binary (1/0)

# Convert months to days
days_per_month <- 30.4375  # Average number of days in a month
imputed_data$time_to_event <- imputed_data$time_to_event * days_per_month

# Loop through seeds 100 to 124
for (seed in 100:124) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Create an index vector with 70% of the data for training
  train_index <- sample(seq_len(nrow(imputed_data)), 0.7 * nrow(imputed_data))
  
  # Subset your data into training and testing sets
  train_data <- imputed_data[train_index, ]
  test_data <- imputed_data[-train_index, ]
  
  # Save the train and test data with the seed value in the filename
  train_filename <- paste0("traindata/unmatched_train_data_seed_", seed, ".csv")
  test_filename <- paste0("testdata/test_data_seed_", seed, ".csv")
  
  write.csv(train_data, file = train_filename, row.names = FALSE)
  write.csv(test_data, file = test_filename, row.names = FALSE)
}

cat("Data splits for seeds 100 to 124 have been saved.\n")




##### extra vars for subset experiments Aug 23 2024
test_data_extended <- test_data
train_data_extended <- train_data
test_data_extended$TBS_binary <-  ifelse(sqrt(test_data_extended$size^2 + test_data_extended$number_liver_mets^2) < 3, 0, 1)
test_data_extended$CEA_binary_200 <- ifelse(test_data_extended$cea_carcinoembryonic.antigen < 200, 0, 1)
test_data_extended$CEA_binary_20 <- ifelse(test_data_extended$cea_carcinoembryonic.antigen < 20, 0, 1)
test_data_extended$number_binary <- ifelse(test_data_extended$number_liver_mets <= 1, 0, 1)
test_data_extended$size_binary <- ifelse(test_data_extended$size < 5, 0, 1)
#write.csv(test_data_extended, file = "testdata_subset/test_data_extended.csv", row.names = FALSE)
train_data_extended$TBS_binary <-  ifelse(sqrt(train_data_extended$size^2 + train_data_extended$number_liver_mets^2) < 3, 0, 1)
train_data_extended$CEA_binary_200 <- ifelse(train_data_extended$cea_carcinoembryonic.antigen < 200, 0, 1)
train_data_extended$CEA_binary_20 <- ifelse(train_data_extended$cea_carcinoembryonic.antigen < 20, 0, 1)
train_data_extended$number_binary <- ifelse(train_data_extended$number_liver_mets <= 1, 0, 1)
train_data_extended$size_binary <- ifelse(train_data_extended$size < 5, 0, 1)
#write.csv(train_data_extended, file = "traindata_subset/unmatched_train_data_extended.csv", row.names = FALSE)
