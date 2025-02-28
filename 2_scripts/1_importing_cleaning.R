# Warfarin INR clinical research script=========================================
# 
#
#
# Written by : Samuel Tekle,MD
# Copyright (c) - 2025 Samuel Tekle
# 
# Started on : 6.2.2025
# Version control:
#
# Objectives: division of data set for convenient exploration and analysis

# Step 1: Load the necessary packages===========================================


# Load required packages using pacman


pacman::p_load(
  here,         # Simplifies file path management for reproducible workflows
  rio,          # Easy import/export of data in various formats (e.g., CSV, Excel)
  tidyverse,    # A collection of packages for data wrangling, visualization, and analysis
  janitor,      # Tools for cleaning and exploring messy data (e.g., clean_names())
  stringr,      # Functions for working with and manipulating text (strings)
  labelled,     # Adds and manages variable labels for better data understanding
  dplyr,        # Grammar for data manipulation (e.g., filter, mutate, summarize)
  purrr,        # Enhances functional programming for working with lists and vectors
  skimr,        # Generates quick and informative summaries of data sets
  forcats       # Simplifies working with categorical variables (factors)
)

# Step 2: Import the whole data set=============================================

df<- import(here("1_input", "Warfarin Data.sav"))

## Step 3: Subset the data into clinically relevant categories===================

### 1. Demographic Information =================================================


#### Step 1: Select relevant columns from the original data frame================
# This step creates a subset of the data with only the demographic columns we need for analysis.

demographic_data <- df %>% 
  select(
    WARFARIN_DATA_ENTRY_ID,  # Unique identifier for each case
    ZOBA,                    # Zonal address of the patient
    SUB_ZOBA,                # Sub-zonal address of the patient
    AGE,                     # Age of the patient
    GENDER,                  # Gender of the patient
    RELIGION,                # Religion of the patient
    ETHINICITY,              # Ethnicity of the patient
    MARIETAL_STATUS,         # Marital status of the patient
    EDUCATIONAL_STATUS,      # Educational level of the patient
    INCOME,                  # Monthly income of the patient
    OCCUPATION,              # Occupation of the patient
    FAMILY_SIZE              # Family size of the patient
  )

#### Step 2: Convert the data to a data frame (if it isn't already)=============
# This ensures that the data is in the correct format for further processing.
demographic_data <- as.data.frame(demographic_data)

####  Step 3: Recode all columns with labels automatically======================
# This step converts labeled variables into factors with proper labels.
# If a column has labels (e.g., 1 = "Male", 2 = "Female"), it will be converted to a factor.
demographic_data <- demographic_data %>%
  mutate(across(everything(), ~ 
                  if (!is.null(attributes(.)$labels)) {
                    factor(., 
                           levels = attributes(.)$labels, 
                           labels = names(attributes(.)$labels))
                  } else {
                    .
                  }
  ))

####  Step 4: Clean variable names and redefine Age as continuous===============
# This step cleans the column names (e.g., removes special characters, converts to lowercase)
# and ensures that the "age" column is treated as a numeric variable.
demographic_data <- demographic_data %>% 
  clean_names() %>%  # Clean column names (remove special characters, convert to lowercase, etc.)
  mutate(age = as.numeric(age))  # Convert age to numeric for continuous analysis

#### Step 5: Label the demographic dataset======================================
# This step adds descriptive labels to each variable for better understanding and documentation.
demographic_data <- demographic_data %>%
  set_variable_labels(
    warfarin_data_entry_id = "Unique identifier for each case",
    zoba = "Zonal Address of the patient",
    sub_zoba = "Sub-zonal Address of the patient",
    age = "Age (years)",
    gender = "Gender of the patient",
    religion = "Religion",
    ethinicity = "Ethnicity",
    marietal_status = "Marital status",
    educational_status = "Educational Level",
    income = "Monthly income",
    occupation = "Profession/Job",
    family_size = "Family size"
  )

# View the final cleaned and labeled demographic dataset
head(demographic_data)


### 2. Clinical History=========================================================

#### Step 1: Select relevant columns from the original data frame===============
# This step creates a subset of the data with only the columns we need for analysis.
clinical_history_data <- df %>% 
  select(
    WARFARIN_DATA_ENTRY_ID,  # Unique identifier for each case
    ANY_COMORBIDITIES,       # Presence of any comorbidities
    TYPE_OF_COMORBIDITIES,   # Type of comorbidities (multi-response variable)
    OTHER_COMORBIDITIES,     # Additional comorbidities not listed
    COMORBIDITY_1,           # First additional comorbidity
    COMORBIDITY_2,           # Second additional comorbidity
    ANY_SUBSTANCE_ABUSE,     # Presence of substance abuse
    TYPE                     # Type of substance abused
  )

#### Step 2: Define a function to clean multi-response variables================
# This function takes a data frame and a variable name, splits the multi-response variable into separate columns,
# and converts it into a binary format (1 for presence, 0 for absence).
clean_multi_response <- function(data, var_name) {
  # Extract labels from the variable (if any)
  labels <- attr(data[[var_name]], "labels")
  
  # Transform the data
  data_cleaned <- data %>%
    select(WARFARIN_DATA_ENTRY_ID, all_of(var_name)) %>%  # Keep the ID and the target variable
    mutate(SPLIT = strsplit(as.character(.data[[var_name]]), "")) %>%  # Split the multi-response string into individual values
    unnest(SPLIT) %>%  # Unnest the split values into separate rows
    mutate(
      SPLIT = factor(SPLIT, levels = labels, labels = names(labels)),  # Convert values to factors with proper labels
      Value = 1  # Add a binary indicator for presence
    ) %>%
    pivot_wider(
      names_from = SPLIT,  # Create new columns for each unique value
      values_from = Value,  # Fill with 1 for presence
      values_fill = list(Value = 0)  # Fill missing combinations with 0
    ) %>%
    select(-all_of(var_name))  # Remove the original variable after splitting
  
  return(data_cleaned)
}

#### Step 3: Apply the function to the multi-response variable "TYPE_OF_COMORBIDITIES"====

# This step cleans the multi-response variable and joins the cleaned data back to the main data frame.
cleaned_var_tc_data <- clean_multi_response(clinical_history_data, "TYPE_OF_COMORBIDITIES")

#### Step 4: Create a cleaned data frame with selected columns==================
cleaned_chx_data <- clinical_history_data %>% 
  select(
    WARFARIN_DATA_ENTRY_ID,  # Unique identifier
    ANY_COMORBIDITIES,       # Presence of any comorbidities
    OTHER_COMORBIDITIES,     # Additional comorbidities
    COMORBIDITY_1,           # First additional comorbidity
    COMORBIDITY_2,           # Second additional comorbidity
    ANY_SUBSTANCE_ABUSE,     # Substance abuse
    TYPE                     # Type of substance abused
  )

#### Step 5: Join the cleaned multi-response data back to the main data frame====
cleaned__chx_data <- left_join(cleaned_chx_data, cleaned_var_tc_data, by = "WARFARIN_DATA_ENTRY_ID")

#### Step 6: Clean column names ================================================
# Replace spaces, hyphens, and convert to lowercase for consistency.
colnames(cleaned__chx_data) <- gsub(" ", "_", colnames(cleaned__chx_data))  # Replace spaces with underscores
colnames(cleaned__chx_data) <- gsub("-", "_", colnames(cleaned__chx_data))  # Replace hyphens with underscores
colnames(cleaned__chx_data) <- tolower(colnames(cleaned__chx_data))  # Convert to lowercase

#### Step 7: Label all variable values==========================================
# Convert labeled variables to factors with proper labels.
cleaned__chx_data <- cleaned__chx_data %>%
  mutate(across(everything(), ~ 
                  if (!is.null(attributes(.)$labels)) {
                    factor(., levels = attributes(.)$labels, labels = names(attributes(.)$labels))
                  } else {
                    .
                  }
  ))

#### Step 8: Clean variable names and reclassify variable types=================
# Convert specific columns to factors for better analysis.
cleaned__chx_data <- cleaned__chx_data %>% 
  clean_names() %>%  # Clean column names (remove special characters, etc.)
  mutate(
    any_comorbidities = as.factor(any_comorbidities),
    other_comorbidities = as.factor(other_comorbidities),
    any_substance_abuse = as.factor(any_substance_abuse),
    ischemic_heart_disease = as.factor(ischemic_heart_disease),
    hypertension = as.factor(hypertension),
    cva = as.factor(cva),
    diabetes_mellitus = as.factor(diabetes_mellitus),
    copd_asthma = as.factor(copd_asthma),
    chronic_kidney_disease = as.factor(chronic_kidney_disease),
    chronic_liver_disease = as.factor(chronic_liver_disease),
    seizure_disorder = as.factor(seizure_disorder),
    hiv_aids = as.factor(hiv_aids)
  )

#### Step 9: Create a total comorbidities column ===============================
# Calculate the total number of comorbidities for each patient.
comorbidity_columns <- c(
  "ischemic_heart_disease", 
  "hypertension", 
  "cva", 
  "diabetes_mellitus", 
  "copd_asthma", 
  "chronic_kidney_disease", 
  "chronic_liver_disease", 
  "seizure_disorder", 
  "hiv_aids"
)

# Convert comorbidity columns to numeric (replace non-numeric values with 0)
cleaned__chx_data[comorbidity_columns] <- lapply(cleaned__chx_data[comorbidity_columns], function(x) {
  as.numeric(as.character(x))
})

# Replace NA values with 0 in comorbidity columns
cleaned__chx_data[comorbidity_columns][is.na(cleaned__chx_data[comorbidity_columns])] <- 0

# Calculate the total number of comorbidities for each row
cleaned__chx_data$total_comorbidities <- rowSums(cleaned__chx_data[comorbidity_columns], na.rm = TRUE)

# Add 1 to the total if "other_comorbidities" is "Yes"
cleaned__chx_data$total_comorbidities <- cleaned__chx_data$total_comorbidities + 
  ifelse(cleaned__chx_data$other_comorbidities == "Yes", 1, 0)

# Replace 1 with "Yes" and 0 with "No" for patients with comorbidities
cleaned__chx_data[cleaned__chx_data$total_comorbidities > 0, comorbidity_columns] <- 
  lapply(cleaned__chx_data[cleaned__chx_data$total_comorbidities > 0, comorbidity_columns], function(x) {
    ifelse(x == 1, "Yes", "No")
  })

# Set comorbidity values to NA for patients with no comorbidities
cleaned__chx_data[cleaned__chx_data$total_comorbidities == 0, comorbidity_columns] <- NA

#### Step 10: Label the dataset ================================================
# Add descriptive labels to each variable for better understanding.
cleaned__chx_data <- cleaned__chx_data %>%
  set_variable_labels(
    warfarin_data_entry_id = "Unique identifier for each case",
    any_comorbidities = "Presence of common chronic conditions",
    other_comorbidities = "Additional comorbidities",
    comorbidity_1 = "1st Additional comorbidity",
    comorbidity_2 = "2nd Additional comorbidity",
    any_substance_abuse = "Substance abuse",
    type = "Type of substance abused",
    ischemic_heart_disease = "Ischemic heart disease",
    hypertension = "Hypertension",
    cva = "Cerebrovascular accident or stroke",
    diabetes_mellitus = "Diabetes Mellitus",
    copd_asthma = "Chronic Obstructive Pulmonary Disorders or Asthma",
    chronic_kidney_disease = "Chronic Kidney disease",
    chronic_liver_disease = "Chronic Liver disease",
    seizure_disorder = "Seizure Disorder",
    hiv_aids = "HIV/AIDs",
    total_comorbidities = "Total number of comorbidities"
  )

# View the final cleaned and labeled dataset
head(cleaned__chx_data)


### 3. Warfarin Indication and Interventional History ==========================
####Step 1: Select relevant columns from the original data frame ===============
#This step creates a subset of the data with only the columns we need for analysis.

indication_data <- df %>% select(WARFARIN_DATA_ENTRY_ID, INDICATION, 
                                 Indication_cleaned, INTERVENTIONAL_STATUS, 
                                 TYPE_OF_INTERVENTION, VALVE_S_INTERVENED)


####Step 2: Create a cleaned data frame with selected columns ==================
#Start with the primary identifier and key variables for further processing.

indication_cleaned_data <- indication_data %>% select(WARFARIN_DATA_ENTRY_ID,
                                                     Indication_cleaned, 
                                                      INTERVENTIONAL_STATUS
                                                      )  # Start with the primary identifier

####Step 3: Define a function to clean multi-response variables ================
#This function takes a data frame and a variable name, splits the multi-response 
# variable into separate columns, and converts it into a binary format (1 for presence, 
#0 for absence).


# List of variables to clean
variables_to_clean <- c("TYPE_OF_INTERVENTION", 
                        "VALVE_S_INTERVENED")


# Function to clean multi-response variables
clean_multi_response <- function(data, var_name) {
  # Extract labels
  labels <- attr(data[[var_name]], "labels")
  
  # Transform data
  data_cleaned <- data %>%
    select(WARFARIN_DATA_ENTRY_ID, all_of(var_name)) %>%
    mutate(SPLIT = strsplit(as.character(.data[[var_name]]), "")) %>%
    unnest(SPLIT) %>%
    mutate(SPLIT = factor(SPLIT, 
                          levels = labels, 
                          labels = names(labels)),
           Value = 1) %>%  # Add a binary indicator for presence
    pivot_wider(names_from = SPLIT, 
                values_from = Value, 
                values_fill = list(Value = 0))  # Fill missing combinations with 0
  
  return(data_cleaned)
}

####Step 4: Apply the function to multi-response variables =====================
#List of variables to clean
# Loop through each variable to clean and merge
for (var in variables_to_clean) {
  cleaned_var_data <- clean_multi_response(indication_data, var)
  indication_cleaned_data <- left_join(indication_cleaned_data, cleaned_var_data, 
                            by = "WARFARIN_DATA_ENTRY_ID")
}

indication_cleaned_data <- indication_cleaned_data %>% 
  clean_names()

head(indication_cleaned_data)



#### Step 5 : Calculate the total_type_of_intervention =========================
# This counts the number of interventions performed for each patient based on the `type_of_intervention` column.
# NA values are treated as no intervention.
indication_cleaned_data <- indication_cleaned_data %>%
  mutate(
    total_type_of_intervention = ifelse(is.na(type_of_intervention), 0, 
                                        nchar(type_of_intervention))
  )

####Step 6: Calculate the total_number_of_valves_intervened=====================
# This sums up the number of valves intervened for each patient based on the 
# `valve_s_intervened` column. NA values are treated as no intervention.

indication_cleaned_data <- indication_cleaned_data %>%
  mutate(
    total_number_of_valves_intervened = ifelse(is.na(valve_s_intervened), 0, 
                                               nchar(valve_s_intervened))
  )


####Step 7: Convert binary columns to "Yes" and "No" ==========================
#List of columns to convert

binary_columns <- c(
  "valvular_repair", 
  "metalic_replacement", 
  "biologic_prosthesis", 
  "mitral_valve", 
  "aortic_valve", 
  "tricuspid_valve"
)

# Convert specified columns to "Yes" and "No", leaving NA as is
indication_cleaned_data <- indication_cleaned_data %>%
  mutate(across(all_of(binary_columns), ~ 
                  case_when(
                    . == 1 ~ "Yes",
                    . == 0 ~ "No",
                    TRUE ~ NA_character_  # Leave NA values unchanged
                  )
  )
)

####Step 8: Replace values in the type_of_intervention column ==================
#Convert codes to descriptive labels.
# Replace values in the type_of_intervention column
indication_cleaned_data <- indication_cleaned_data %>%
  mutate(
    type_of_intervention = case_when(
      type_of_intervention == "A" ~ "Valvular repair only",
      type_of_intervention == "AB" ~ "Both valvular repair and replacement",
      type_of_intervention == "B" ~ "Metallic replacement",
      type_of_intervention == "C" ~ "Biologic prosthesis replacement",
      TRUE ~ type_of_intervention  # Leave other values (including NA) unchanged
    )
  )


####Step 9: Replace values in the valve_s_intervened column ===================
#Convert codes to descriptive labels.

# Replace values in the valve_s_intervened column

indication_cleaned_data <- indication_cleaned_data %>%
  mutate(
    valve_s_intervened = case_when(
      valve_s_intervened == "A" ~ "MV only",
      valve_s_intervened == "AB" ~ "MV & AV",
      valve_s_intervened == "ABC" ~ "MV, AV, & TV",
      valve_s_intervened == "AC" ~ "MV & TV",
      valve_s_intervened == "B" ~ "AV only",
      TRUE ~ valve_s_intervened  # Leave other values (including NA) unchanged
    )
  )


#### Step 10: Label all variable values==========================================
# Convert labeled variables to factors with proper labels.

indication_cleaned_data <- indication_cleaned_data %>%
  mutate(across(everything(), ~ 
                  if (!is.null(attributes(.)$labels)) {
                    factor(., levels = attributes(.)$labels, 
                           labels = names(attributes(.)$labels))
                  } else {
                    .
                  }
  ))

colnames(indication_cleaned_data)

#### Step 11: Classify the variable type========================================
# Convert labeled variables to factors with proper labels.
# classification of the variable types 

indication_cleaned_data <- indication_cleaned_data %>% 
  mutate(
    indication_cleaned = as.factor(indication_cleaned),
    interventional_status = as.factor(interventional_status),
    type_of_intervention = as.factor(type_of_intervention),
    metalic_replacement = as.factor(metalic_replacement),
    biologic_prosthesis = as.factor(biologic_prosthesis),
    valve_s_intervened = as.factor(valve_s_intervened),
    mitral_valve = as.factor(mitral_valve),
    aortic_valve = as.factor(aortic_valve),
    tricuspid_valve = as.factor(tricuspid_valve),
    total_type_of_intervention = as.numeric(total_type_of_intervention),
    total_number_of_valves_intervened = as.numeric(total_number_of_valves_intervened)
      )


#### Step 12: Label the dataset ================================================
# Add descriptive labels to each variable for better understanding.
indication_cleaned_data <- indication_cleaned_data %>%
  set_variable_labels(
    warfarin_data_entry_id = "Unique identifier for each case",
    indication_cleaned = "Cleaned indication for intervention (categorical)",
    interventional_status = "Status of intervention (categorical)",
    type_of_intervention = "Type of intervention performed",
    valvular_repair = "Whether valvular repair was performed (Yes/No)",
    metalic_replacement = "Whether metallic replacement was performed (Yes/No)",
    biologic_prosthesis = "Whether biologic prosthesis was used (Yes/No)",
    valve_s_intervened = "Valves intervened during the procedure",
    mitral_valve = "Whether the mitral valve was intervened (Yes/No)",
    aortic_valve = "Whether the aortic valve was intervened (Yes/No)",
    tricuspid_valve = "Whether the tricuspid valve was intervened (Yes/No)",
    total_type_of_intervention = "Total number of types of interventions performed",
    total_number_of_valves_intervened = "Total number of valves intervened"
  )



### 4. Warfarin Treatment Details===============================================

treatment_details_data <- df %>% select(WARFARIN_DATA_ENTRY_ID, DD_warfarin_days, 
                                        WARFARIN_DURATION, GOAL_THERAPEUTIC_RANGE, 
                                        starts_with("WEEKLY_WARFARIN_DOSE")
                                        )

#### Step 1: Label all variable values==========================================
# Convert labeled variables to factors with proper labels.

treatment_details_data <- treatment_details_data %>%
  mutate(across(everything(), ~ 
                  if (!is.null(attributes(.)$labels)) {
                    factor(., levels = attributes(.)$labels, 
                           labels = names(attributes(.)$labels))
                  } else {
                    .
                  }
  ))

### Step 2: classify the variables nature=======================================

treatment_details_data <- treatment_details_data %>% 
  mutate(
    DD_warfarin_days = as.numeric(DD_warfarin_days),
    WARFARIN_DURATION = as.numeric(WARFARIN_DURATION),
    GOAL_THERAPEUTIC_RANGE = as.factor(GOAL_THERAPEUTIC_RANGE),
    across(starts_with("WEEKLY_WARFARIN_DOSE"), as.numeric)  # Corrected line
  ) %>% 
  clean_names()


####step 3: Add descriptive labels to each variable-============================

treatment_details_data <- treatment_details_data %>%
  set_variable_labels(
    warfarin_data_entry_id = "Unique identifier for each case",
    dd_warfarin_days = "Number of days on warfarin therapy",
    warfarin_duration = "Total duration of warfarin therapy (in days)",
    goal_therapeutic_range = "Target therapeutic range for warfarin",
    weekly_warfarin_dose_01 = "Weekly warfarin dose for week 1",
    weekly_warfarin_dose_02 = "Weekly warfarin dose for week 2",
    weekly_warfarin_dose_03 = "Weekly warfarin dose for week 3",
    weekly_warfarin_dose_04 = "Weekly warfarin dose for week 4",
    weekly_warfarin_dose_05 = "Weekly warfarin dose for week 5",
    weekly_warfarin_dose_06 = "Weekly warfarin dose for week 6",
    weekly_warfarin_dose_07 = "Weekly warfarin dose for week 7",
    weekly_warfarin_dose_08 = "Weekly warfarin dose for week 8",
    weekly_warfarin_dose_09 = "Weekly warfarin dose for week 9",
    weekly_warfarin_dose_10 = "Weekly warfarin dose for week 10",
    weekly_warfarin_dose_11 = "Weekly warfarin dose for week 11",
    weekly_warfarin_dose_12 = "Weekly warfarin dose for week 12",
    weekly_warfarin_dose_13 = "Weekly warfarin dose for week 13",
    weekly_warfarin_dose_14 = "Weekly warfarin dose for week 14",
    weekly_warfarin_dose_15 = "Weekly warfarin dose for week 15",
    weekly_warfarin_dose_16 = "Weekly warfarin dose for week 16",
    weekly_warfarin_dose_17 = "Weekly warfarin dose for week 17",
    weekly_warfarin_dose_18 = "Weekly warfarin dose for week 18",
    weekly_warfarin_dose_19 = "Weekly warfarin dose for week 19",
    weekly_warfarin_dose_20 = "Weekly warfarin dose for week 20",
    weekly_warfarin_dose_21 = "Weekly warfarin dose for week 21",
    weekly_warfarin_dose_22 = "Weekly warfarin dose for week 22",
    weekly_warfarin_dose_23 = "Weekly warfarin dose for week 23",
    weekly_warfarin_dose_24 = "Weekly warfarin dose for week 24"
  )


### 5. INR Monitoring ==========================================================

inr_monitoring_data <- df %>% select(WARFARIN_DATA_ENTRY_ID, starts_with("INR"), 
                                     TTR_percentage, TTR_GOOD_POOR, Days_in_TR)

#### Step 1: Label all variable values==========================================
# Convert labeled variables to factors with proper labels.

inr_monitoring_data <- inr_monitoring_data %>%
  mutate(across(everything(), ~ 
                  if (!is.null(attributes(.)$labels)) {
                    factor(., levels = attributes(.)$labels, 
                           labels = names(attributes(.)$labels))
                  } else {
                    .
                  }
  ))

#### Step 2: Label all variable values==========================================
# Convert labeled variables to factors with proper labels.


inr_monitoring_data <- inr_monitoring_data %>% 
  clean_names()

colnames(inr_monitoring_data)

####step 3: Add descriptive labels to each variable-============================
# Label the variables in the dataset
var_label(inr_monitoring_data) <- list(
  warfarin_data_entry_id = "Unique ID for warfarin data entry",
  inr_01 = "INR value at visit 1",
  inr_02 = "INR value at visit 2",
  inr_03 = "INR value at visit 3",
  inr_04 = "INR value at visit 4",
  inr_05 = "INR value at visit 5",
  inr_06 = "INR value at visit 6",
  inr_07 = "INR value at visit 7",
  inr_08 = "INR value at visit 8",
  inr_09 = "INR value at visit 9",
  inr_10 = "INR value at visit 10",
  inr_11 = "INR value at visit 11",
  inr_12 = "INR value at visit 12",
  inr_13 = "INR value at visit 13",
  inr_14 = "INR value at visit 14",
  inr_15 = "INR value at visit 15",
  inr_16 = "INR value at visit 16",
  inr_17 = "INR value at visit 17",
  inr_18 = "INR value at visit 18",
  inr_19 = "INR value at visit 19",
  inr_20 = "INR value at visit 20",
  inr_21 = "INR value at visit 21",
  inr_22 = "INR value at visit 22",
  inr_23 = "INR value at visit 23",
  inr_24 = "INR value at visit 24",
  ttr_percentage = "Time in therapeutic range (TTR) percentage",
  ttr_good_poor = "TTR classification (Good/Poor)",
  days_in_tr = "Total days within the therapeutic range"
)



### 6. Adverse Events===========================================================

adverse_events_data <- df %>% select(WARFARIN_DATA_ENTRY_ID, ANY_ADVERSE_EVENT, 
                                     TYPE_OF_ADVERSE_EVENT, BLEEDING_NATURE)


#### step 1: changing multiresponse variables to individual columns=============

clean_multi_response <- function(data, var_name) {
  # Extract labels from the variable (if any)
  labels <- attr(data[[var_name]], "labels")
  
  # Transform the data
  data_cleaned <- data %>%
    select(WARFARIN_DATA_ENTRY_ID, all_of(var_name)) %>%  # Keep the ID and the target variable
    mutate(SPLIT = strsplit(as.character(.data[[var_name]]), "")) %>%  # Split the multi-response string into individual values
    unnest(SPLIT) %>%  # Unnest the split values into separate rows
    mutate(
      SPLIT = factor(SPLIT, levels = labels, labels = names(labels)),  # Convert values to factors with proper labels
      Value = 1  # Add a binary indicator for presence
    ) %>%
    pivot_wider(
      names_from = SPLIT,  # Create new columns for each unique value
      values_from = Value,  # Fill with 1 for presence
      values_fill = list(Value = 0)  # Fill missing combinations with 0
    ) %>%
    select(-all_of(var_name))  # Remove the original variable after splitting
  
  return(data_cleaned)
}

#### Step 2: Apply the function to the multi-response variable "TYPE_OF_ADVERSE_EVENT, BLEEDING_NATURE"====

  
variables_to_clean <- c("TYPE_OF_ADVERSE_EVENT", "BLEEDING_NATURE") 

#List of variables to clean
# Loop through each variable to clean and merge
for (var in variables_to_clean) {
  cleaned_var_data <- clean_multi_response(adverse_events_data, var)
  adverserxn_cleaned_data <- left_join(adverse_events_data, cleaned_var_data, 
                                       by = "WARFARIN_DATA_ENTRY_ID")
}

#### step 3: changing 1 and zero to yes or no===================================
colnames(adverserxn_cleaned_data)
# Convert specified columns to "Yes" and "No", leaving NA as is
binary_columns <- c(
  "Gingival Bleeding", 
  "Nose Bleeding" ,"Ecchymosis", "Menorrhagia", 
  "Intra-articular bleeding",
  "GI Bleeding", "Hematuria", "Intracranial bleeding" 
)


adverserxn_cleaned_data <- adverserxn_cleaned_data %>%
  mutate(across(all_of(binary_columns), ~ 
                  case_when(
                    . == 1 ~ "Yes",
                    . == 0 ~ "No",
                    TRUE ~ NA_character_  # Leave NA values unchanged
                  )
  )
  )

#### Step 5 : Calculate the total_bleeding natures seen=========================
# This counts the number of interventions performed for each patient based on the `type_of_intervention` column.
# NA values are treated as no intervention.
adverserxn_cleaned_data <- adverserxn_cleaned_data %>%
  mutate(
    total_bleeding_nature = ifelse(is.na(BLEEDING_NATURE), 0, 
                                        nchar(BLEEDING_NATURE))
  )

#### Step 6: labeling and naming the vatiables=================================
colnames(adverserxn_cleaned_data)

adverserxn_cleaned_data <- adverserxn_cleaned_data %>% 
  clean_names()

# Add descriptive labels to each variable
adverserxn_cleaned_data <- adverserxn_cleaned_data %>%
  set_variable_labels(
    warfarin_data_entry_id = "Unique identifier for each case",
    any_adverse_event = "Presence of any adverse event (Yes/No)",
    type_of_adverse_event = "Type of adverse event experienced",
    bleeding_nature = "Nature of bleeding event",
    gingival_bleeding = "Gingival bleeding (Yes/No)",
    nose_bleeding = "Nose bleeding (Yes/No)",
    ecchymosis = "Ecchymosis (Yes/No)",
    menorrhagia = "Menorrhagia (Yes/No)",
    intra_articular_bleeding = "Intra-articular bleeding (Yes/No)",
    gi_bleeding = "Gastrointestinal (GI) bleeding (Yes/No)",
    hematuria = "Hematuria (Yes/No)",
    intracranial_bleeding = "Intracranial bleeding (Yes/No)",
    total_bleeding_nature = "Total number of bleeding events"
  )
#### step 7: Replace values in the type of adverse rxn column==================

adverserxn_cleaned_data <- adverserxn_cleaned_data %>%
  mutate(
    type_of_adverse_event = case_when(
      type_of_adverse_event == "A" ~ "Bleeding only",
      type_of_adverse_event == "AB" ~ "Both bleeding and VTE",
      type_of_adverse_event == "B" ~ "VTE only",
        TRUE ~ type_of_adverse_event  # Leave other values (including NA) unchanged
    )
  )


### 7. Co-Medications and Drug Interactions=====================================

co_medications_data <- df %>% select(WARFARIN_DATA_ENTRY_ID, ANY_CO_MEDICATIONS,
                                     LIST_OF_COMEDICATIONS, OTHER_MEDICATIONS, 
                                     WARFARIN_DRUG_INTERACTION)



####Step 1: Function to clean WARFARIN_DRUG_INTERACTION mulitple response=======

clean_warfarin_interaction <- function(data) {
  data_cleaned <- data %>%
    select(WARFARIN_DATA_ENTRY_ID, WARFARIN_DRUG_INTERACTION) %>%
    mutate(SPLIT = strsplit(as.character(WARFARIN_DRUG_INTERACTION), "")) %>%
    unnest(SPLIT) %>%
    group_by(WARFARIN_DATA_ENTRY_ID, SPLIT) %>%
    summarise(Count = n(), .groups = 'drop') %>%  # Correctly count occurrences for each patient
    pivot_wider(names_from = SPLIT, 
                values_from = Count, 
                values_fill = list(Count = 0)) %>%  # Fill missing combinations with 0
    rowwise() %>%
    mutate(TOTAL_INTERACTIONS = sum(c_across(c(A, B, C, D)), na.rm = TRUE),
           COUNT_A = ifelse(!is.na(A), A, 0),
           COUNT_B = ifelse(!is.na(B), B, 0),
           COUNT_C = ifelse(!is.na(C), C, 0),
           COUNT_D = ifelse(!is.na(D), D, 0)) %>%  # Count interactions for each severity level
    ungroup()
  
  return(data_cleaned)
}


# Clean WARFARIN_DRUG_INTERACTION
interaction_comx_data <- clean_warfarin_interaction(co_medications_data)

# Keep original column of WARFARIN_DRUG_INTERACTION
interaction_comx_data <- left_join( co_medications_data, interaction_comx_data,
                                   by = "WARFARIN_DATA_ENTRY_ID")


#### Step 1.1: Clean the multi-response list of comedications ==================

clean_multi_response <- function(data, var_name) {
  # Extract labels from the variable (if any)
  labels <- attr(data[[var_name]], "labels")
  
  # Transform the data
  data_cleaned <- data %>%
    select(WARFARIN_DATA_ENTRY_ID, all_of(var_name)) %>%  # Keep the ID and the target variable
    mutate(SPLIT = strsplit(as.character(.data[[var_name]]), "")) %>%  # Split the multi-response string into individual values
    unnest(SPLIT) %>%  # Unnest the split values into separate rows
    mutate(
      SPLIT = factor(SPLIT, levels = labels, labels = names(labels)),  # Convert values to factors with proper labels
      Value = 1  # Add a binary indicator for presence
    ) %>%
    pivot_wider(
      names_from = SPLIT,  # Create new columns for each unique value
      values_from = Value,  # Fill with 1 for presence
      values_fill = list(Value = 0)  # Fill missing combinations with 0
    ) %>%
    select(-all_of(var_name))  # Remove the original variable after splitting
  
  return(data_cleaned)
}

# Apply the function to the LIST_OF_COMEDICATIONS column
list_comx_data <- clean_multi_response(interaction_comx_data, "LIST_OF_COMEDICATIONS")

interaction_comx_data <- left_join(interaction_comx_data,list_comx_data, 
                                   by= "WARFARIN_DATA_ENTRY_ID" )

####step 2: cleaning and labeling and defining variables========================

interaction_comx_data <- interaction_comx_data %>% 
  clean_names() %>% 
  mutate(any_co_medications = case_when(
    any_co_medications == 1 ~  "Yes",
    any_co_medications == 2 ~  "No"
    ),
    any_co_medications = as.factor(any_co_medications)
  )



# Label the variables in the dataset
var_label(interaction_comx_data) <- list(
  warfarin_data_entry_id = "Unique ID for warfarin data entry",
  any_co_medications = "Indicator for co-medications (Yes/No)",
  list_of_comedications = "List of co-administered medications",
  other_medications = "Other medications taken by the patient",
  warfarin_drug_interaction = "Indicator for warfarin drug interactions",
  
  # Interaction types
  a = "Type A interaction (Minor or no clinical effect)",
  b = "Type B interaction (Unknown mechanism)",
  c = "Type C interaction (Monitor)",
  d = "Type D interaction (Avoid)",
  
  # Interaction counts
  total_interactions = "Total number of detected drug interactions",
  count_a = "Count of Type A interactions",
  count_b = "Count of Type B interactions",
  count_c = "Count of Type C interactions",
  count_d = "Count of Type D interactions",
  
  # Medications
  ramipril = "Ramipril",
  atrovastatin = "Atorvastatin ",
  bisoprolol = "Bisoprolol",
  furosamide = "Furosemide",
  atenolol = "Atenolol",
  digoxin = "Digoxin",
  omeprazole = "Omeprazole",
  isosorbide_dinitrate = "Isosorbide dinitrate",
  aspirin = "Aspirin",
  spironolactone = "Spironolactone",
  erythromycin = "Erythromycin",
  benzathine_benzylpenicillin = "Benzathine benzylpenicillin",
  dmp = "DMP",
  enalapril = "Enalapril",
  cardivolol = "Carvedilol",
  phenobarbitone = "Phenobarbital",
  amlodipin = "Amlodipine"
)

# Convert specified columns to "Yes" and "No", leaving NA as is
binary_columns <- c( "ramipril",                   
                      "atrovastatin",                "bisoprolol",                  "furosamide",                 
                      "atenolol",                    "digoxin",                     "omeprazole",                 
                     "isosorbide_dinitrate",        "aspirin",                     "spironolactone",             
                      "erythromycin",                "benzathine_benzylpenicillin", "dmp",                        
                      "enalapril",                   "cardivolol",                  "phenobarbitone",             
                      "amlodipin"                  

)


interaction_comx_data <- interaction_comx_data %>%
  mutate(across(all_of(binary_columns), ~ 
                  case_when(
                    . == 1 ~ "Yes",
                    . == 0 ~ "No",
                    TRUE ~ NA_character_  # Leave NA values unchanged
                  )
  )
  )

# remove reduntdant columns of the interaction counts

interaction_comx_data <- interaction_comx_data %>% 
  select(-c(a,b,c,d))

colnames(interaction_comx_data)

### 8. Warfarin Knowledge and Adherence=========================================

knowledge_adherence_data <- df %>% select(WARFARIN_DATA_ENTRY_ID, 
                                          WHAT_DOSE_OF_WARFARIN_WERE_YOU_T, 
                                          HOW_WOULD_YOU_KNOW_THE_MEDICATJ, 
                                          DO_YOU_KNOW_WHY_YOU_ARE_TAKING_W, 
                                          DO_REMEMBER_YOUR_LAST_INR, 
                                          WHEN_IS_YOUR_NEXT_BLOOD_TEST, 
                                          HAVE_YOU_MISSED_OR_TAKEN_EXCESSI)

####step 1: adding the total score of knowledge out of 6========================

head(knowledge_adherence_data)



# Create a new column that sums up all columns except WARFARIN_DATA_ENTRY_ID
knowledge_adherence_data <- knowledge_adherence_data %>%
  mutate(
    overall_knowledge = ifelse(
      rowSums(is.na(select(., -WARFARIN_DATA_ENTRY_ID))) > 0,  # Check if any column has NA
      NA,  # If any column has NA, set overall_knowledge to NA
      rowSums(select(., -WARFARIN_DATA_ENTRY_ID), na.rm = TRUE)  # Otherwise, calculate the sum
    )
  )
    
    # View the updated data frame
    head(knowledge_adherence_data)

colnames(knowledge_adherence_data)

#### step 2=====================================================================

knowledge_adherence_data <- knowledge_adherence_data %>% 
  set_variable_labels(
    overall_knowledge = "Sum of the component knowledge scores"
  )

## Save the cleaned subsets=====================================================
export(demographic_data, here("3_output", "demographic_data.rds"))
export(demographic_data, here("3_output", "demographic_data.csv"))

export(cleaned__chx_data, here("3_output", "cleaned__chx_data.rds"))
export(cleaned__chx_data, here("3_output", "cleaned__chx_data.csv"))

export(indication_cleaned_data, here("3_output", "indication_cleaned_data.rds"))
export(indication_cleaned_data, here("3_output", "indication_cleaned_data.csv"))


export(treatment_details_data, here("3_output", "treatment_details_data.rds"))
export(treatment_details_data, here("3_output", "treatment_details_data.csv"))


export(inr_monitoring_data, here("3_output", "inr_monitoring_data.rds"))
export(inr_monitoring_data, here("3_output", "inr_monitoring_data.csv"))



export(adverserxn_cleaned_data, here("3_output", "adverserxn_cleaned_data.rds"))
export(adverserxn_cleaned_data, here("3_output", "adverserxn_cleaned_data.csv"))


export(interaction_comx_data, here("3_output", "interaction_comx_data.rds"))
export(interaction_comx_data, here("3_output", "interaction_comx_data.csv"))

export(knowledge_adherence_data, here("3_output", "knowledge_adherence_data.rds"))
export(knowledge_adherence_data, here("3_output", "knowledge_adherence_data.csv"))




