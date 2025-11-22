####################################################################################
#
# Program name: LHS000501.R
#
# Description:
#   This R script processes the 2015 Colombia DHS Children’s Recode (KR file)
#   and generates a derived dataset used for the analysis of under-5 mortality.
#
# Version: R version 4.5.1
#
# Input:
#   - Raw DHS data file: COKR72FL.SAS7BDAT (Children’s Recode 2015)
#
# Output:
#   - Derived dataset: LHS000501.rda
#
# Required packages:
#   - haven (for reading DHS SAS files)
#   - dplyr (for data manipulation)
#
# Date created: 2025-11-22
#
####################################################################################

# Clear the workspace
rm(list = ls(all = TRUE))

# Load required packages
library(haven)
library(dplyr)
library(tidyr)

# Load the DHS Children’s Recode (KR) dataset Colombia (2015)
kr_2015 <- read_sas("data/raw/COKR72FL.SAS7BDAT")

# ------------------------------------------------------------------------------
# Derive variables
# ------------------------------------------------------------------------------
LHS00050_0 <- kr_2015 %>%
  rename(
    region = V101
  ) %>%
  mutate(
    # -------------------------------
    #  Design variables
    # -------------------------------
    strata = V022,
    psu = V001,
    weight = V005 / 1000000,
    normalized_weight = weight / mean(weight),

    # -------------------------------
    # Date of interview
    # -------------------------------
    date_interview = V008,

    # -------------------------------
    # Age at death
    # -------------------------------
    age_at_death = ifelse(is.na(B7), 999, B7),

    # -------------------------------
    # Age at interview
    # -------------------------------
    age_at_interview = date_interview - B3,

    # -------------------------------
    # Time variable
    # -------------------------------
    time = pmin(age_at_interview, age_at_death),

    # -------------------------------
    # Status variable
    # -------------------------------
    status = ifelse(age_at_death < age_at_interview, 1, 0),

    # -------------------------------
    # Child year of birth
    # -------------------------------
    yob = trunc(1900 + (B3 - 1) / 12),

    # -------------------------------
    # time interval (7 categories)
    # -------------------------------
    time_interval_c7 = case_when(
      time == 0 ~ 1,
      time %in% 1:5 ~ 2,
      time %in% 6:11 ~ 3,
      time %in% 12:23 ~ 4,
      time %in% 24:35 ~ 5,
      time %in% 36:47 ~ 6,
      time %in% 48:59 ~ 7
    ),

    # -------------------------------
    # Ethnicity (3 categories)
    # -------------------------------
    ethnicity_c3 = case_when(
      V131 == 1 ~ 1, # Indigenous
      V131 %in% c(3, 4, 5) ~ 2, # Afro-descendant
      V131 %in% c(2, 6) ~ 3 # Other
    ),
    # -------------------------------
    # ID variable
    # -------------------------------
    id = row_number()
  ) %>%
  # Keep only the requested variables
  select(
    id, region, strata, psu, normalized_weight, status, time, time_interval_c7, ethnicity_c3, yob
  ) %>%
  # Exclusions
  filter(yob < 2016 & time < 60) # keep only births before 2016 and under-5

# Function to create event variable
event_var <- function(df) {
  status_val <- unique(df$status)
  if (status_val == 0) {
    df$event <- rep(0, nrow(df))
  } else {
    df$event <- c(rep(0, nrow(df) - 1), 1)
  }
  df$age <- 1:nrow(df)
  return(df)
}

# Expanding the dataset to person-months format
LHS000501 <- LHS00050_0[rep(row_number(LHS00050_0), LHS00050_0$time_interval_c7), ] %>%
  group_by(id) %>%
  group_split() %>%
  lapply(., event_var) %>%
  bind_rows() %>%
  select(-time_interval_c7, -status, -time) %>%
  rename(status = event)

# Save the derived dataset
save(LHS000501, file = "data/derived/LHS000501.rda")
