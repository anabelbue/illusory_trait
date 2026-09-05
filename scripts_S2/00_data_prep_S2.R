# Pace data
library(tidyverse)
library(corrplot)
library(here)
library(psych)


#### load data, select variables, recode variables ####
pace_data = readRDS("/Users/anabelbue/Documents/01_Phd/05_PACE study/pace_analysis/data/pace_esm_data_s1.rds") %>% 
# rename variables
  dplyr::rename(
    "E1" = "soci",
    "E2" = "asse",
    "E3" = "acti",
    "A1" = "trus",
    "A2" = "kind",
    "A3" = "poli",
    "O1" = "crea",
    "O2" = "curi",
    "O3" = "arti",
    "C1" = "lazy",
    "C2" = "orga",
    "C3" = "reli",
    "N1" = "irri",
    "N2" = "rela",
    "N3" = "care",
    "happiness" = "happ", 
    "stress" = "stres", 
    "tiredness" = "tire", 
    "hunger" = "hung", 
    "sociability" = "sociality"
  )

# reverse code items
pace_data = pace_data |>
  mutate(
    C1 = 11 - C1,
    N2 = 11 - N2,
    N3 = 11 - N3
  )

# rename id and create composite scores 
pace_data = pace_data %>% 
  mutate(id = session, 
         state_e = rowMeans(select(., E1, E2, E3), na.rm = TRUE),
         state_n = rowMeans(select(., N1, N2, N3), na.rm = TRUE),
         state_o = rowMeans(select(., O1, O2, O3), na.rm = TRUE),
         state_c = rowMeans(select(., C1, C2, C3), na.rm = TRUE),
         state_a = rowMeans(select(., A1, A2, A3), na.rm = TRUE),
  )


vars <- c("state_e", "state_n", "state_o", "state_c", "state_a", 
          "duty", "sociability", "negativity", "typicality", 
          "happiness", "stress", "hunger", "tiredness")
          

# select variables of interest and order data by id and "created"
pace_data = pace_data |>
  select(id, condition, n_day, n_survey_per_day, n_sent_invitations, created, modified, ended, all_of(vars)) |>
  arrange(id, created)


#### data cleaning #### 
# remove empty surveys (which are duplicates)

# remove surveys that were not filled in (i.e., "modified" is NA)
pace_data = pace_data |>
  filter(!is.na(modified))

# check if there are any remaining duplicate surveys
pace_data |>
  group_by(id, created) |>
  filter(n() > 1)
# --> no

# remove rows where all items are missing
pace_data = pace_data |>
  filter(!if_all(all_of(vars), is.na))


# remove observations where the difference in time between starting and finishing the survey is greater than intended 
pace_data_r <- pace_data %>%
  mutate(
    timedifference = difftime(created, if_else(is.na(ended), modified, ended))
  ) %>%
  filter(
    is.na(timedifference) |
      (condition == "A" & abs(timedifference) < 2400) | # 40 minutes = 2400 seconds
      (condition == "B" & abs(timedifference) < 4200)  # 70 minutes = 4200 seconds
  )

nrow(pace_data) - nrow(pace_data_r) # 32 observations excluded 



# after day 4 it starts with day 3 again, so we could simply add a +1 to the day where it restarts again 

#id2: the not unique number 25 has the same day number as the first time (4)
# if this happens, I could check whether there are any distinct n_days in between (which would be day 5 here)
# and then assign the second number with the latter entry the following distinct number (6) 
# the difference between the newly assigned number and the following numbers then need to be added to all following days 
# if there is no day in between, then we have to create one, so simply add 1 plus one to that day and all following days
# the corrected n_sent_invitations colum should then be the n_day*n_survey_per_day 

# There may be two types: those where 
pace_data_corrected <- pace_data_r %>%
  group_by(id) %>%
  mutate(
    reset = n_day < lag(n_day, default = first(n_day)) |
      (n_day == lag(n_day, default = first(n_day)) & 
         as.numeric(difftime(created, lag(created, default = first(created)), units = "hours")) >= 12),
    offset = {
      offsets <- rep(0, n())
      n_day_temp <- n_day
      for (i in which(reset)) {
        last_day_before <- n_day_temp[i - 1]
        first_day_after <- n_day_temp[i]
        additional_offset <- last_day_before - first_day_after + 1
        offsets[i:n()] <- offsets[i:n()] + additional_offset
        n_day_temp[i:n()] <- n_day_temp[i:n()] + additional_offset
      }
      offsets
    },
    n_day_corrected = n_day + offset,
    modified = offset != 0
  ) %>%
  group_by(id, n_day_corrected) %>%
  mutate(
    n_survey_per_day_corrected = if_else(
      !reset & duplicated(n_survey_per_day),
      n_survey_per_day + cumsum(duplicated(n_survey_per_day)),
      n_survey_per_day
    ),
    corrected = !reset & duplicated(n_survey_per_day),
    has_zero = any(n_survey_per_day_corrected == 0),
    corrected = corrected | has_zero,
    n_survey_per_day_corrected = n_survey_per_day_corrected + if_else(has_zero, 1L, 0L)
  ) %>%
  ungroup() %>%
  mutate(
    n_sent_invitations_corrected = case_when(
      condition == "A" ~ ((n_day_corrected - 1) * 8) + n_survey_per_day_corrected,
      condition == "B" ~ ((n_day_corrected - 1) * 4) + n_survey_per_day_corrected
    )
  ) %>%
  select(-reset, -offset, -has_zero)


# (a) Check whether n_sent_invitations_corrected is unique per participant
pace_data_corrected %>%
  group_by(id) %>%
  filter(duplicated(n_sent_invitations_corrected)) %>%
  distinct(id, n_sent_invitations_corrected)

# (b) Check whether each combination of n_day_corrected and n_survey_per_day_corrected is unique per participant
pace_data_corrected %>%
  group_by(id, n_day_corrected, n_survey_per_day_corrected) %>%
  filter(n() > 1) %>%
  distinct(id, n_day_corrected, n_survey_per_day_corrected)



pace_data <- pace_data_corrected %>% 
  mutate(n_day = n_day_corrected, 
         n_survey_per_day = n_survey_per_day_corrected, 
         n_sent_invitations = n_sent_invitations_corrected) %>% 
  select(-n_survey_per_day_corrected, -n_day_corrected, -n_sent_invitations_corrected)



# create missing observations for surveys that were missed 
pace_data_filled <- pace_data %>%
  group_by(id) %>%
  mutate(n_surveys_max = if_else(condition == "A", 8L, 4L)) %>%
  group_by(id, n_day) %>%
  complete(
    n_survey_per_day = 1:first(n_surveys_max)
  ) %>%
  mutate(
    id = first(na.omit(id)),
    condition = first(na.omit(condition)),
    n_day = first(na.omit(n_day)),
    n_sent_invitations = case_when(
      condition == "A" ~ ((n_day - 1) * 8) + n_survey_per_day,
      condition == "B" ~ ((n_day - 1) * 4) + n_survey_per_day
    )
  ) %>%
  ungroup() %>%
  select(-n_surveys_max)




#Downsample condition A: split every day in two days that are added directly after each other 
pace_data_downsampled <- pace_data_filled %>%
  filter(condition == "A") %>%
  mutate(
    part = if_else(n_survey_per_day %in% c(1, 3, 5, 7), "odd", "even")
  ) %>%
  group_by(id) %>%
  mutate(
    max_day = max(n_day),
    n_day = if_else(
      part == "odd",
      n_day,                    # odd surveys: keep original day
      max_day + n_day           # even surveys: append after all original days
    ),
    n_survey_per_day = if_else(
      part == "odd",
      match(n_survey_per_day, c(1, 3, 5, 7)),
      match(n_survey_per_day, c(2, 4, 6, 8))
    ),
    n_sent_invitations = ((n_day - 1) * 4) + n_survey_per_day
  ) %>%
  ungroup() %>%
  select(-part, -max_day) %>%
  bind_rows(pace_data_filled %>% filter(condition == "B")) %>%
  arrange(id, n_day, n_survey_per_day)


# Find optimal 8-day window per participant
pace_data_final <- pace_data_downsampled %>%
  group_by(id) %>%
  group_modify(~{
    days <- unique(.x$n_day)
    
    if (length(days) < 8) return(.x %>% mutate(keep = FALSE, best_missing_rate = 1))
    
    best_start <- 1L
    best_missing_rate <- Inf
    
    for (i in 1:(length(days) - 7)) {
      window_days <- days[i:(i + 7)]
      window_data <- .x %>% filter(n_day %in% window_days)
      missing_rate <- mean(is.na(window_data$created))
      
      if (missing_rate < best_missing_rate) {
        best_missing_rate <- missing_rate
        best_start <- i
      }
    }
    
    best_window <- days[best_start:(best_start + 7)]
    
    .x %>%
      mutate(
        keep = n_day %in% best_window,
        best_missing_rate = best_missing_rate
      )
  }) %>%
  ungroup() %>%
  filter(keep, best_missing_rate <= 0.50) %>%
  group_by(id) %>%
  mutate(timepoint = row_number()) %>%
  ungroup() %>%
  select(-keep, -best_missing_rate, -n_day, -n_survey_per_day, -n_sent_invitations)

distinct(pace_data_final, id) # 970


# check how many non missing values are in the final sample 
pace_sample_check <- pace_data_final %>% 
  na.omit() 

nrow(pace_sample_check)

obs_per_person <- pace_sample_check %>% 
  group_by(id) %>% 
  summarize(n = n())

describe(obs_per_person$n)

#response rate
27.29/32

# Add baseline data ----------------------------------------------------------
baseline_data <-  read_csv("/Users/anabelbue/Documents/01_Phd/05_PACE study/pace_analysis/data/pace_baseline_data_s1.csv")
baseline_codebook <- readxl::read_excel("/Users/anabelbue/Documents/01_Phd/05_PACE study/pace_analysis/codebooks/baseline_codebook_pace.xlsx", 
                                        sheet =5)

recode_bfi_items <- baseline_codebook %>% 
  filter(startsWith(item_name, "BFI")) %>% 
  filter(recoded == 1) %>% 
  pull(item_name)

baseline_data[recode_bfi_items] <- 11- baseline_data[recode_bfi_items]

trait_data <- baseline_data %>% 
  mutate(id = session, 
         trait_e = rowMeans(select(., BFI_assertive_1, BFI_assertive_2, BFI_energy_1, BFI_energy_2, BFI_sociability_1, BFI_sociability_2), na.rm = TRUE),
         trait_n = rowMeans(select(., BFI_anxiety_1, BFI_anxiety_2, BFI_depression_1, BFI_depression_2, BFI_emotional_1, BFI_emotional_2), na.rm = TRUE),
         trait_o = rowMeans(select(., BFI_aesthetic_1, BFI_aesthetic_2,BFI_creative_1, BFI_creative_2, BFI_intellectual_1, BFI_intellectual_2), na.rm = TRUE),
         trait_a = rowMeans(select(., BFI_compassion_1, BFI_compassion_2, BFI_respect_1, BFI_respect_2, BFI_trust_1, BFI_trust_2), na.rm = TRUE),
         trait_c = rowMeans(select(., BFI_organization_1, BFI_organization_2, BFI_productive_1, BFI_productive_2, BFI_responsibility_1, BFI_responsibility_2), na.rm = TRUE)) %>% 
  select(id, gender, age, trait_e,  trait_n, trait_o, trait_a, trait_c)

# merge with esm data 
dat <- pace_data_final %>% 
  left_join(trait_data, by ="id")


# Find the covariates by inspecting the correlations ----------------------
dat <- dat %>%
  group_by(id) %>% 
  mutate(state_e_pm = mean(state_e, na.rm = TRUE),
         state_n_pm = mean(state_n, na.rm = TRUE), 
         state_o_pm = mean(state_o, na.rm = TRUE), 
         state_a_pm = mean(state_a, na.rm = TRUE), 
         state_c_pm = mean(state_c, na.rm = TRUE),
         happiness_pm = mean(happiness, na.rm = TRUE), 
         stress_pm = mean(stress, na.rm =TRUE), 
         tiredness_pm = mean(tiredness, na.rm = TRUE), 
         hunger_pm = mean(hunger, na.rm = TRUE), 
         duty_pm = mean(duty, na.rm =TRUE), 
         negativity_pm = mean(negativity, na.rm = TRUE),
         typicality_pm = mean(typicality, na.rm = TRUE),
         sociability_pm = mean(sociability, na.rm = TRUE)) %>%  
  ungroup() %>% 
  mutate(state_e_pmc = state_e - state_e_pm,
         state_n_pmc = state_n - state_n_pm, 
         state_o_pmc = state_o - state_o_pm, 
         state_a_pmc = state_a - state_a_pm, 
         state_c_pmc = state_c - state_c_pm,
         happiness_pmc = happiness - happiness_pm, 
         stress_pmc = stress - stress_pm, 
         tiredness_pmc = tiredness - tiredness_pm,
         hunger_pmc = hunger - hunger_pm, 
         duty_pmc = duty - duty_pm, 
         negativity_pmc = negativity - negativity_pm, 
         typicality_pmc = typicality - typicality_pm, 
         sociability_pmc = sociability -   sociability_pm)


# Within-person (WP) and between-person (BP) correlations
cor_table_lower_empty <- function(data ){
  table1 <- data
  table1[upper.tri(table1, diag=TRUE)] <-""
  table1 <- matrix(as.numeric(unlist(table1)),nrow=nrow(table1))
  table1 <- as.data.frame(table1)
}

# correlation person-means
wp_dat <- dat %>% select(contains("pmc"))
cor_wp<- round(cor(wp_dat, use ="complete.obs"),2)
corrplot(cor_wp,  type="lower", diag= TRUE, cl.cex = 1, tl.cex = 1, tl.col="black")
wp_cors_table <- cor_table_lower_empty(cor_wp)
names(wp_cors_table) <- names(wp_dat)
writexl::write_xlsx(wp_cors_table, "S2_wp_cor_table.xlsx")

# Plot BP correlation 
bp_dat <- dat %>% 
  select(matches("trait|pm$"), id) %>% 
  distinct(., id, .keep_all = TRUE) %>% 
  select(-id)
cor_bp<- round(cor(bp_dat, use ="complete.obs"),2)
corrplot(cor_bp,  type="lower", diag= TRUE, cl.cex = 1, tl.cex = 1, tl.col="black")
bp_cors_table <- cor_table_lower_empty(cor_bp)
names(bp_cors_table) <- names(bp_dat)
writexl::write_xlsx(bp_cors_table, "S2_bp_cor_table.xlsx")


# Create wide dataset -----------------------------------------------------
names(dat)

# Create wide dataset 
wide_dat <- dat %>%
  pivot_wider(
    id_cols = c(id, trait_e:trait_c), 
    names_from = timepoint,  
    values_from = c(state_e:tiredness)
  )
names(wide_dat)

write_csv(wide_dat, here("data", "wide_dat_S2.csv"))

#check response rates
response_rates <- wide_dat %>%
  select(starts_with("state_e_")) %>%
  summarise(across(everything(), ~ mean(!is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "survey", values_to = "response_rate") %>% 
  select(response_rate) %>% pull()


#check demographics
demo_dat <- dat %>% distinct(id, .keep_all = TRUE)
table(demo_dat$gender)
784/nrow(demo_dat)

psych::describe(demo_dat$age)

