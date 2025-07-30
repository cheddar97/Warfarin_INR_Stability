# Copyright (c) [Year] [Your Name or Organization]
# This script is licensed under the MIT License.
# Purpose: Calculate the INR Response Index, defined as the change in INR divided by the change in the weekly Warfarin dose.
# Patients with fewer than 3 measurements are excluded to avoid errors.
# Load necessary libraries
library(dplyr)
library(labelled)
library(tidyr)
library(lubridate)
library(here)
library(haven)

# Step 1: Import the dataset
df <- read_sav(here("1_input", "Warfarin Data.sav"))

# Define the columns for dates, INR values, and weekly Warfarin doses
date_cols <- paste0("DATE$", sprintf("%02d", 1:24))  
inr_cols <- paste0("INR$", sprintf("%02d", 1:24))    
dose_cols <- paste0("WEEKLY_WARFARIN_DOSE$", sprintf("%02d", 1:24))

# Extract relevant columns
response_data <- df %>%
  select(WARFARIN_DATA_ENTRY_ID, all_of(c(date_cols, inr_cols, dose_cols)))

# Reshape data into long format
response_long <- response_data %>%
  pivot_longer(
    cols = c(all_of(date_cols), all_of(inr_cols), all_of(dose_cols)),
    names_to = c(".value", "measurement"),
    names_pattern = "(.*)\\$(.*)"
  ) %>%
  filter(!is.na(DATE) & !is.na(INR) & !is.na(WEEKLY_WARFARIN_DOSE)) %>%
  mutate(
    DATE = as.Date(DATE, format = "%Y-%m-%d"),
    measurement = as.numeric(measurement)
  ) %>%
  arrange(WARFARIN_DATA_ENTRY_ID, DATE)

# Filter out patients with fewer than 3 measurements
response_long <- response_long %>%
  group_by(WARFARIN_DATA_ENTRY_ID) %>%
  filter(n() >= 3) %>%
  ungroup()

# Function to calculate INR Response Index per patient
calculate_inr_response_index <- function(patient_data) {
  iri_values <- numeric()
  dose_changes <- 0
  
  for (i in 1:(nrow(patient_data) - 1)) {
    start_inr <- patient_data$INR[i]
    end_inr <- patient_data$INR[i + 1]
    start_dose <- patient_data$WEEKLY_WARFARIN_DOSE[i]
    end_dose <- patient_data$WEEKLY_WARFARIN_DOSE[i + 1]
    
    # Skip iteration if any value is NA
    if (is.na(start_inr) || is.na(end_inr) || is.na(start_dose) || is.na(end_dose)) {
      iri_values <- c(iri_values, NA)
      next
    }
    
    # Calculate INR change and dose change
    delta_inr <- end_inr - start_inr
    delta_dose <- end_dose - start_dose
    
    # Count number of dose adjustments
    if (delta_dose != 0) {
      dose_changes <- dose_changes + 1
    }
    
    # Compute INR Response Index
    iri <- ifelse(delta_dose != 0, abs(delta_inr / delta_dose), NA)
    iri_values <- c(iri_values, iri)
  }
  
  # Compute single INR Response Index per patient
  mean_iri <- mean(iri_values, na.rm = TRUE)  # Mean INR response index
  sd_iri <- sd(iri_values, na.rm = TRUE)      # Standard deviation of IRI
  adj_iri <- mean_iri * (1 + (sd_iri / mean_iri)) * (1 + log1p(dose_changes))  # Adjust for variability & changes
  
  return(adj_iri)
}

# Calculate adjusted INR Response Index per patient
response_results <- response_long %>%
  group_by(WARFARIN_DATA_ENTRY_ID) %>%
  summarise(Adjusted_INR_Response_Index = calculate_inr_response_index(cur_data())) %>%
  ungroup()


response_results <- response_results %>% 
  clean_names() %>% 
  set_variable_labels(
    adjusted_inr_response_index = "Adjusted INR response index"
  )

# Save results
saveRDS(response_results, here("3_output", "adjusted_inr_response_index.rds"))
write.csv(response_results, here("3_output", "adjusted_inr_response_index.csv"), row.names = FALSE)

# Display results
print(response_results)
