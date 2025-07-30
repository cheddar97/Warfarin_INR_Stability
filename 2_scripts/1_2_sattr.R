# Copyright (c) [2025] [Samuel Tekle, MD]
# This script is licensed under License.
# Purpose: Calculate Stability-Adjusted Time in Therapeutic Range (saTTR) for patients on Warfarin therapy.
# saTTR is defined as TTR / (1 + SD[INR]), where SD[INR] is the standard deviation of INR values.

# Load necessary libraries
library(janitor)
library(dplyr)    # For data manipulation
library(tidyr)    # For reshaping data
library(lubridate)# For date manipulation
library(here)     # For file path management
library(haven)    # For importing .sav files
library(labelled)

# Step 1: Import the dataset
# Load the Warfarin dataset from the specified path
df <- read_sav(here("1_input", "Warfarin Data.sav"))

# Define the columns for dates and INR values
date_cols <- paste0("DATE$", sprintf("%02d", 1:24))  # Columns for dates (e.g., DATE$01 to DATE$24)
inr_cols <- paste0("INR$", sprintf("%02d", 1:24))    # Columns for INR values (e.g., INR$01 to INR$24)
goal_range <- "GOAL_THERAPEUTIC_RANGE"               # Column for therapeutic range

# Extract relevant columns
# Select only the necessary columns for TTR and saTTR calculation
ttr_data <- df %>%
  select(WARFARIN_DATA_ENTRY_ID, all_of(c(date_cols, inr_cols, goal_range, 
                                          "WARFARIN_DURATION")))

# Reshape the data into a long format for easier processing
# Convert wide-format data (multiple columns for dates and INR values) into long format
ttr_long <- ttr_data %>%
  pivot_longer(
    cols = c(all_of(date_cols), all_of(inr_cols)),
    names_to = c(".value", "measurement"),
    names_pattern = "(.*)\\$(.*)"
  ) %>%
  filter(!is.na(DATE) & !is.na(INR)) %>%  # Remove rows with missing dates or INR values
  mutate(
    DATE = as.Date(DATE, format = "%Y-%m-%d"),  # Convert DATE to Date type
    measurement = as.numeric(measurement)       # Convert measurement number to numeric
  ) %>%
  arrange(WARFARIN_DATA_ENTRY_ID, DATE)  # Sort by patient ID and date

# Define mapping for GOAL_THERAPEUTIC_RANGE
# Map the therapeutic range labels (e.g., "A", "B") to their corresponding min and max values
goal_mapping <- list("A" = c(2.0, 3.0), "B" = c(2.5, 3.5))

# Convert GOAL_THERAPEUTIC_RANGE values
# Extract the min and max values of the therapeutic range for each patient
ttr_long <- ttr_long %>%
  mutate(
    GOAL_MIN = sapply(GOAL_THERAPEUTIC_RANGE, function(x) goal_mapping[[x]][1]),
    GOAL_MAX = sapply(GOAL_THERAPEUTIC_RANGE, function(x) goal_mapping[[x]][2])
  )

# Check for NAs in GOAL_MIN and GOAL_MAX
# Issue a warning if any rows have invalid therapeutic range values
if (any(is.na(ttr_long$GOAL_MIN))) {
  warning("Some rows have invalid GOAL_THERAPEUTIC_RANGE values. Check the format.")
}

# Filter out individuals with fewer than 3 INR measurements
# Ensure that only patients with sufficient data are included in the analysis
ttr_long <- ttr_long %>%
  group_by(WARFARIN_DATA_ENTRY_ID) %>%
  filter(n() >= 3) %>%  # Keep only patients with 3 or more measurements
  ungroup()

# Display the extracted and reshaped data
# Preview the processed data
head(ttr_long)

# Function to calculate TTR for a single patient
# This function calculates the Time in Therapeutic Range (TTR) using the Rosendaal method
calculate_ttr_for_patient <- function(patient_data) {
  total_time_in_range <- 0
  total_observation_time <- 0
  
  # Loop through each interval between measurements
  for (i in 1:(nrow(patient_data) - 1)) {
    start_date <- patient_data$DATE[i]
    end_date <- patient_data$DATE[i + 1]
    start_inr <- patient_data$INR[i]
    end_inr <- patient_data$INR[i + 1]
    goal_min <- patient_data$GOAL_MIN[i]
    goal_max <- patient_data$GOAL_MAX[i]
    
    # Skip iteration if any value is NA or invalid
    if (is.na(start_inr) || is.na(end_inr) || is.na(goal_min) || is.na(goal_max)) {
      next
    }
    
    # Calculate the duration of the interval in days
    interval_duration <- as.numeric(difftime(end_date, start_date, units = "days"))
    
    # Linear interpolation of INR values over the interval
    time_in_range <- 0
    
    if (start_inr == end_inr) {
      # If INR is constant, check if it's within the therapeutic range
      if (start_inr >= goal_min && start_inr <= goal_max) {
        total_time_in_range <- total_time_in_range + interval_duration
      }
    } else {
      # Calculate the proportion of time spent in the therapeutic range
      if (start_inr < goal_min && end_inr < goal_min) {
        # Both INR values below the range
        time_in_range <- 0
      } else if (start_inr > goal_max && end_inr > goal_max) {
        # Both INR values above the range
        time_in_range <- 0
      } else {
        # Calculate the intersection points with the therapeutic range
        if (start_inr < goal_min) {
          # INR starts below the range and enters it
          time_to_enter <- (goal_min - start_inr) / (end_inr - start_inr) * interval_duration
          time_in_range <- interval_duration - time_to_enter
        } else if (start_inr > goal_max) {
          # INR starts above the range and enters it
          time_to_enter <- (start_inr - goal_max) / (start_inr - end_inr) * interval_duration
          time_in_range <- interval_duration - time_to_enter
        } else {
          # INR starts within the range
          if (end_inr < goal_min) {
            # INR exits the range below
            time_to_exit <- (start_inr - goal_min) / (start_inr - end_inr) * interval_duration
            time_in_range <- time_to_exit
          } else if (end_inr > goal_max) {
            # INR exits the range above
            time_to_exit <- (goal_max - start_inr) / (end_inr - start_inr) * interval_duration
            time_in_range <- time_to_exit
          } else {
            # INR stays within the range for the entire interval
            time_in_range <- interval_duration
          }
        }
      }
      total_time_in_range <- total_time_in_range + time_in_range
    }
    
    # Add the interval duration to the total observation time
    total_observation_time <- total_observation_time + interval_duration
  }
  
  # Calculate TTR as the proportion of time spent in the therapeutic range
  if (total_observation_time > 0) {
    ttr <- (total_time_in_range / total_observation_time) * 100
  } else {
    ttr <- NA  # Return NA if no valid observation time
  }
  return(ttr)
}

# Calculate TTR and SD[INR] for each patient
# Group the data by patient ID and calculate TTR and SD[INR] for each patient
ttr_results <- ttr_long %>%
  group_by(WARFARIN_DATA_ENTRY_ID) %>%
  summarise(
    TTR = calculate_ttr_for_patient(cur_data()),
    SD_INR = sd(INR, na.rm = TRUE)  # Calculate standard deviation of INR values
  ) %>%
  ungroup()

# Calculate Stability-Adjusted TTR (saTTR)
# saTTR is defined as TTR / (1 + SD[INR])
ttr_results <- ttr_results %>%
  mutate(saTTR = TTR / (1 + SD_INR))

ttr_results <- ttr_results %>% 
  clean_names() 

ttr_results <- ttr_results %>% 
  set_variable_labels(
    sa_ttr = "Stability adjusted time in therapuetic range",
    sd_inr= "Standard deviations in inr values",
    ttr = "Time in therapuetic range"
  )

# Export saTTR results in .rds and .csv formats
# Save the results for further analysis or reporting
saveRDS(ttr_results, here("3_output", "sattr_results.rds"))  # Save as .rds
write.csv(ttr_results, here("3_output", "sattr_results.csv"), row.names = FALSE)  # Save as .csv

# Display the saTTR results
print(ttr_results)

colnames(df)