# Copyright (c) [Year] [Your Name or Organization]
# This script is licensed under the MIT License.
# Purpose: Calculate the Dose Stability Index (DSI), adjusted for follow-up duration.

# Load necessary libraries
library(dplyr)    # For data manipulation
library(tidyr)    # For reshaping data
library(lubridate)# For date manipulation
library(here)     # For file path management
library(haven)    # For importing .sav files
library(labelled)
library(janitor)

# Step 1: Import the dataset
# Load the Warfarin dataset from the specified path
df <- read_sav(here("1_input", "Warfarin Data.sav"))

colnames(df)
# Define the columns for dates and weekly Warfarin doses
date_cols <- paste0("DATE$", sprintf("%02d", 1:24))  # Columns for dates (e.g., DATE$01 to DATE$24)
dose_cols <- paste0("WEEKLY_WARFARIN_DOSE$", sprintf("%02d", 1:24))  # Columns for weekly Warfarin doses

# Extract relevant columns
# Select only the necessary columns for DSI calculation
dsi_data <- df %>%
  select(WARFARIN_DATA_ENTRY_ID, all_of(c(date_cols, dose_cols)))

# Reshape the data into a long format for easier processing
# Convert wide-format data (multiple columns for dates and doses) into long format
dsi_long <- dsi_data %>%
  pivot_longer(
    cols = c(all_of(date_cols), all_of(dose_cols)),
    names_to = c(".value", "measurement"),
    names_pattern = "(.*)\\$(.*)"
  ) %>%
  filter(!is.na(DATE) & !is.na(WEEKLY_WARFARIN_DOSE)) %>%  # Remove rows with missing values
  mutate(
    DATE = as.Date(DATE, format = "%Y-%m-%d"),  # Convert DATE to Date type
    measurement = as.numeric(measurement)       # Convert measurement number to numeric
  ) %>%
  arrange(WARFARIN_DATA_ENTRY_ID, DATE)  # Sort by patient ID and date

# Filter out patients with fewer than 2 measurements
# Ensure that only patients with sufficient data are included in the analysis
dsi_long <- dsi_long %>%
  group_by(WARFARIN_DATA_ENTRY_ID) %>%
  filter(n() >= 2) %>%  # Keep only patients with 2 or more measurements
  ungroup()

# Display the extracted and reshaped data
# Preview the processed data
head(dsi_long)

# Function to calculate Dose Stability Index (DSI) for a single patient
# This function calculates the DSI based on the number of dose changes adjusted for follow-up duration
calculate_dsi <- function(patient_data) {
  # Calculate the number of dose changes
  dose_changes <- 0  # Counter for the number of dose changes
  current_dose <- patient_data$WEEKLY_WARFARIN_DOSE[1]  # Initialize with the first dose
  
  for (i in 2:nrow(patient_data)) {
    next_dose <- patient_data$WEEKLY_WARFARIN_DOSE[i]
    
    # Check if the dose has changed
    if (!is.na(current_dose) && !is.na(next_dose) && current_dose != next_dose) {
      dose_changes <- dose_changes + 1  # Increment the dose change counter
    }
    
    # Update the current dose
    current_dose <- next_dose
  }
  
  # Calculate the follow-up duration in months
  follow_up_duration <- as.numeric(
    difftime(patient_data$DATE[nrow(patient_data)], patient_data$DATE[1], units = "days") / 30.44
  )  # Convert days to months (average days per month = 30.44)
  
  # Calculate the rate of dose changes per month
  dose_change_rate <- dose_changes / follow_up_duration
  
  # Calculate the Dose Stability Index (DSI)
  dsi <- 1 / (1 + dose_change_rate)
  
  return(dsi)
}

# Calculate DSI for each patient
# Group the data by patient ID and calculate the DSI
dsi_results <- dsi_long %>%
  group_by(WARFARIN_DATA_ENTRY_ID) %>%
  summarise(
    DSI = calculate_dsi(cur_data())
  ) %>%
  ungroup()


dsi_results <- dsi_results %>% 
  clean_names() %>% 
  set_variable_labels(
    dsi  = "Dose stability index"
  )

# Export DSI results in .rds and .csv formats
# Save the results for further analysis or reporting
saveRDS(dsi_results, here("3_output", "dsi_results.rds"))  # Save as .rds
write.csv(dsi_results, here("3_output", "dsi_results.csv"), row.names = FALSE)  # Save as .csv

# Display the DSI results
print(dsi_results)
