library(dplyr)
library(survival)
library(gtsummary)
library(survminer)

# Load the data
#data <- read.csv("traindata_main/original_data.csv")
#print(nrow(data))

# Drop rows with missing values in the 'adjuvantchemo' column
#data <- data[!is.na(data$adjuvantchemo), ]
#print(nrow(data))

# Remove rows where 'OS' column is equal to 0
#data <- data[data$OS != 0, ]
#print(nrow(data))
data <-read.csv("traindata/matched_train_equalized_r1.csv")
#data <-read.csv("traindata/unmatched_train_data.csv")

# Define categorical and continuous variables
categorical_vars <- c("gender","grade","rightleft.Rec", "extrahepatic.disease", "R0", "KRAS", "DFI.12", "T_stage", "N_stage")
categorical_vars <- c("gender","rightleft.Rec", "extrahepatic.disease", "R0", "KRAS", "DFI.12", "T_stage", "N_stage")

continuous_vars <- c("age", "size", "cea_carcinoembryonic.antigen", "number_liver_mets")



# Function to compute Mann-Whitney U test p-value for continuous variables
mann_whitney_test <- function(var) {
  p_val <- wilcox.test(data[[var]] ~ data$adjuvantchemo, data = data, na.rm = TRUE)$p.value
  return(p_val)
}

# Define function for Chi-squared test
chi_squared_test <- function(var) {
  tbl <- table(data[[var]], data$adjuvantchemo)
  chisq_result <- chisq.test(tbl, correct = FALSE)
  p_val <- chisq_result$p.value
  counts_percentages <- as.data.frame(tbl)
  
  # Calculate percentages
  percentages <- prop.table(tbl, 2) * 100
  counts_percentages$Percentages <- paste0(counts_percentages$Freq, " (", round(percentages[counts_percentages$Var2], 2), "%)")
  
  return(list(p_value = p_val, counts_percentages = counts_percentages))
}

# Apply function to categorical variables
results_list <- list()
for (cat_var in categorical_vars) {
  if (cat_var != "adjuvantchemo") {
    result <- chi_squared_test(cat_var)
    results_list[[cat_var]] <- list(
      Variable = cat_var,
      Test = "Chi-squared",
      p_value = result$p_value,
      Counts_Percentages = paste(result$counts_percentages$Percentages, collapse = ", ")
    )
  }
}

# Mann-Whitney U test and quantiles for continuous variables
for (cont_var in continuous_vars) {
  # Ensure 'adjuvantchemo' is binary and valid for grouping
  if(length(unique(data$adjuvantchemo)) == 2) {
    wilcox_result <- wilcox.test(data[[cont_var]] ~ data$adjuvantchemo)
    
    median_iqr <- tapply(data[[cont_var]], data$adjuvantchemo, function(x) {
      c(Median = median(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE))
    })
    
    median_iqr_treated <- median_iqr[[1]]
    median_iqr_untreated <- median_iqr[[2]]
    
    combined_median_iqr <- paste(
      paste(round(median_iqr_treated["Median"], 2), "(", round(median_iqr_treated["IQR"], 2), ")", sep=""),
      paste(round(median_iqr_untreated["Median"], 2), "(", round(median_iqr_untreated["IQR"], 2), ")", sep=""),
      sep = ", "
    )
    
    results_list[[cont_var]] <- list(
      Variable = cont_var,
      Test = "Mann-Whitney U",
      p_value = wilcox_result$p.value,
      Counts_Percentages = combined_median_iqr
    )
  }
}

# Combine results into a data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame, stringsAsFactors = FALSE))

# Assign column names
colnames(results_df) <- c("Variable", "Test", "p_value", "Counts_Percentages")

# View the results table
print(results_df)