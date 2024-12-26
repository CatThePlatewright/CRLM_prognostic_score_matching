source("calibrate_pac_fct.R")

run_calibration <- function(seed, set_number = 0) {
  print(paste("seed number: ",seed))
  ### LOAD TEST DATA #################################
  
  ######### !!!!!!!! to change to validation!!!!
  test_filename <- paste0("traindata/unmatched_train_data.csv") 
  ######### !!!!!!!! to change to validation!!!!
  test_data <- read.csv(file = test_filename, stringsAsFactors = FALSE)
  test_data$size_binary <- ifelse(test_data$size <= 5, 0, 1)
  
  ######### !!!!!!!! to change to validation!!!!
  untreated_test_data <- subset(test_data, adjuvantchemo == 0)
  treated_test_data <- subset(test_data, adjuvantchemo == 1)
  
  ### LOAD UNMATCHED TRAIN DATA #################################
  train_filename <- paste0("traindata/unmatched_train_data.csv")
  train_data_m1 <- read.csv(file = train_filename, stringsAsFactors = FALSE)
  train_data_m1$size_binary <- ifelse(train_data_m1$size <= 5, 0, 1)
  
  untreated_train_data_2a <- subset(train_data_m1, adjuvantchemo == 0)
  treated_train_data_2b <- subset(train_data_m1, adjuvantchemo == 1)
  
  ### LOAD MATCHED TRAIN DATA #################################
  train_filename_3a <- paste0("traindata/matched_train_set", set_number, "_r1_seed", seed, ".csv")
  train_filename_3b <- paste0("traindata/matched_train_set", set_number, "_seed", seed, ".csv")
  
  train_data_3a <- read.csv(file = "traindata/matched_train_ratio1.0set2.csv", stringsAsFactors = FALSE)
  #train_data_3a$size_binary <- ifelse(train_data_3a$size <= 5, 0, 1)
  
  untreated_train_data_3a <- subset(train_data_3a, adjuvantchemo == 0)
  
  train_data_3b <- read.csv(file = "traindata/matched_train_1_to_1set2.csv", stringsAsFactors = FALSE)
  #train_data_3b$size_binary <- ifelse(train_data_3b$size <= 5, 0, 1)
  
  treated_train_data_3b <- subset(train_data_3b, adjuvantchemo == 1)
  treated_train_data_3b_r1 <- subset(train_data_3a, adjuvantchemo == 1)
  

  models <- list(
    train_data_m1,
    untreated_train_data_2a,
    untreated_train_data_3a,
    train_data_m1,
    treated_train_data_2b,
    treated_train_data_3b#,treated_train_data_3b_r1
  )
  # Loop through each model
  for (i in seq_along(models)) {
    train_data <- models[[i]]
    
    # Define the model label based on i
    if (i < 3) {
      model_label <- paste0(i,"A subset")
      test_data <- untreated_test_data
    } 
    else if (i == 3)
    {
      model_label <- paste0("3A r1 matched on subset")
      test_data <- untreated_test_data
    }else if (i == 6)
    {
      model_label <- paste0("3B  matched on subset")
      test_data <- treated_test_data
    } else {
      model_label <- paste0(i - 3, "B subset")
      test_data <- treated_test_data
    }
    
    file_path <- "plots_and_metrics/calibration_metrics.csv"

    plot <- FALSE
    calibrate_pac(train_data, test_data,file_path,title_suffix = model_label,set_number=set_number,plot=plot)
    
  }
  
 }

# Loop through the seeds and call the function
run_calibration(seed,0)


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
#write.csv(oe_table, "plots_and_metrics/calibration_metrics.csv", row.names = FALSE)


