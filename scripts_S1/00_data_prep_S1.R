# S1 Data prep & descriptive statistics 

## load packages
library(tidyverse)
library(here)
library(readr)
library(psych)
library(corrplot)
library(lubridate)
library(hms)
library(psych)

## load data 
state_dat <- read_csv(here("data", "Independent Effects L1 Data.csv"))
trait_dat <- read_csv(here("data", "Independent Effects L2 Data.csv"))
goal_dat <- read_csv(here("data", "Lifelogging_L1_ 7_13_15.csv"))
schedule <- readxl::read_excel(here("data", "text_schedule_S1.xlsx"))


# combine state and trait data 
c_dat <- merge(trait_dat, state_dat, by = "SID")
nrow(distinct(c_dat, SID)) #209

# add goal data
gc_dat <- merge(c_dat, goal_dat, by = "SID")
nrow(distinct(gc_dat, SID)) #138
## SID seems to be the participant id
nrow(distinct(state_dat, SID)) # 210
nrow(distinct(goal_dat, SID)) #138

### not all participants provided goal data = remove 

# remove non valid observations, i.e. that were completed more than an hour after receiving the invitation 
o_dat <- c_dat %>% filter(valid == 1)
nrow(c_dat) - nrow(o_dat) # 

# only keep particiapnts with at least 10 observations 
n_dat <- o_dat %>% 
         group_by(SID) %>% 
         mutate(n_obs_p = n()) %>% 
         filter(n_obs_p >=10) %>% 
          ungroup()
nrow(distinct(c_dat, SID)) - nrow(distinct(n_dat, SID))

## describe final sample
dat <- n_dat
nrow(distinct(dat, SID)) #204 participants 
nrow(dat)
describe(distinct(dat, SID, n_obs_p)$n_obs_p)



demo_dat <- dat %>% distinct(SID, .keep_all = TRUE)
table(demo_dat$Sex)
134/nrow(demo_dat)
psych::describe(demo_dat$Age)

# prepare the variables 
items_to_recode <- c("HEXACO53", "HEXACO35", "HEXACO41", "HEXACO59", # neuroticism
                     "HEXACO28", "HEXACO52", "HEXACO10","HEXACO46",  # extraversion
                     "HEXACO09", "HEXACO15", "HEXACO57", "HEXACO21",  # agreeableness
                     "HEXACO26", "HEXACO32", "HEXACO14", "HEXACO20", "HEXACO44", "HEXACO56", # conscientiousness
                     "HEXACO01", "HEXACO31", "HEXACO49", "HEXACO19", "HEXACO55") # openness

dat[items_to_recode] <- 6 - dat[items_to_recode]
dat <- dat %>%
  mutate(id = SID, 
         state_e = rowMeans(select(., ExpX, ExpDom), na.rm = TRUE),
         state_n = ExpE, 
         state_o = ExpO, 
         state_a = ExpA, 
         state_c = ExpC, 
         happiness = ExpHappy, 
         selfesteem= ExpFeel, 
         duty = ExpDuty, 
         intellect = ExpIntel, 
         adversity = ExpAdv, 
         mating = ExpMate, 
         positivity = ExpPos, 
         negativity = ExpNeg, 
         deception = ExpDec, 
         sociability = ExpSoc, 
         trait_e = rowMeans(select(., HEXACO04, HEXACO28, HEXACO52, HEXACO10, HEXACO34, HEXACO58, HEXACO16, HEXACO40, HEXACO22, HEXACO46), na.rm = TRUE),
         trait_n = rowMeans(select(., HEXACO05, HEXACO29, HEXACO53, HEXACO11, HEXACO35, HEXACO17, HEXACO41, HEXACO23, HEXACO47, HEXACO59), na.rm = TRUE),
         trait_o = rowMeans(select(., HEXACO01, HEXACO25, HEXACO07, HEXACO31, HEXACO13, HEXACO37, HEXACO49, HEXACO19, HEXACO43, HEXACO55), na.rm = TRUE), 
         trait_a = rowMeans(select(., HEXACO03, HEXACO27, HEXACO09, HEXACO33, HEXACO51, HEXACO15, HEXACO39, HEXACO57, HEXACO21, HEXACO45), na.rm = TRUE),
         trait_c = rowMeans(select(., HEXACO02, HEXACO26, HEXACO08, HEXACO32, HEXACO14, HEXACO38, HEXACO50, HEXACO20, HEXACO44, HEXACO56), na.rm = TRUE)) %>% 
  group_by(SID) %>% 
  mutate(state_e_pm = mean(state_e, na.rm = TRUE),
         state_n_pm = mean(state_n, na.rm = TRUE), 
         state_o_pm = mean(state_o, na.rm = TRUE), 
         state_a_pm = mean(state_a, na.rm = TRUE), 
         state_c_pm = mean(state_c, na.rm = TRUE),
         happiness_pm = mean(happiness, na.rm = TRUE), 
         selfesteem_pm = mean(selfesteem, na.rm =TRUE), 
         duty_pm = mean(duty, na.rm =TRUE), 
         intellect_pm = mean(intellect, na.rm = TRUE),
         adversity_pm = mean(adversity, na.rm = TRUE),
         mating_pm = mean(mating, na.rm = TRUE),
         positivity_pm = mean(positivity, na.rm = TRUE),
         negativity_pm = mean(negativity, na.rm = TRUE),
         deception_pm = mean(deception, na.rm = TRUE),
         sociability_pm = mean(sociability, na.rm = TRUE)) %>%  
  ungroup() %>% 
  mutate(state_e_pmc = state_e - state_e_pm,
         state_n_pmc = state_n - state_n_pm, 
         state_o_pmc = state_o - state_o_pm, 
         state_a_pmc = state_a - state_a_pm, 
         state_c_pmc = state_c - state_c_pm,
         happiness_pmc = happiness - happiness_pm, 
         selfesteem_pmc = selfesteem - selfesteem_pm, 
         duty_pmc = duty - duty_pm, 
         intellect_pmc = intellect - intellect_pm, 
         adversity_pmc = adversity - adversity_pm, 
         mating_pmc = mating - mating_pm, 
         positivity_pmc = positivity - positivity_pm, 
         negativity_pmc = negativity - negativity_pm, 
         deception_pmc = deception - deception_pm, 
         sociability_pmc = sociability - sociability_pm)



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
writexl::write_xlsx(wp_cors_table, "S1_wp_cor_table.xlsx")

# Plot BP correlation 
bp_dat <- dat %>% 
          select(matches("trait|pm$"), SID) %>% 
          distinct(., SID, .keep_all = TRUE) %>% 
          select(-SID)
cor_bp<- round(cor(bp_dat, use ="complete.obs"),2)
corrplot(cor_bp,  type="lower", diag= TRUE, cl.cex = 1, tl.cex = 1, tl.col="black")
bp_cors_table <- cor_table_lower_empty(cor_bp)
names(bp_cors_table) <- names(bp_dat)
writexl::write_xlsx(bp_cors_table, "S1_bp_cor_table.xlsx")

# Bring dataset in correct format -----------------------------------------
## every participant need to have exactly 56 rows, meaning that rows for missing surveys need to be added 
## this requires to identify the survey number of each response such that missings can be added in the right places 

prep_dat <- dat %>% select(id:trait_c, `Start Day`, StartDate)

weekday_levels <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")

# the surveys schedule contains the information at what time on which day the surveys are send out 
long_schedule <- schedule %>%
  mutate(across(Monday:Sunday, ~ as_hms(.))) %>% 
  pivot_longer(cols = Monday:Sunday, names_to = "day", values_to = "survey_time") %>%
  rename(interval = Interval)

# Process all participants
time_dat <- prep_dat %>%
  mutate(
    day = `Start Day`, 
    start_date = mdy_hm(StartDate),  # Use `mdy_hm()` since `StartDate` lacks seconds
    start_time = as_hms(start_date)  # Extract the time correctly
  ) %>% 
  arrange(id, start_date) %>%
  group_by(id) %>%
  mutate(
    first_day = match(first(day), weekday_levels),  # participants could start on different days, but the survey scheduled stayed alaways the same 
    day_num = match(day, weekday_levels),  # Convert weekday to numeric
    day_relative = (day_num - first_day) %% 7 + 1  # Ensure values are always 1–7 
  ) %>%
  left_join(long_schedule, by = "day") %>%
  mutate(
    time_diff = as.numeric(start_time) - as.numeric(survey_time),  # Difference in seconds
    time_diff = ifelse(time_diff < 0, NA, time_diff)  # Keep only earlier surveys
  ) %>%
  group_by(id, day, start_time) %>%
  slice_min(time_diff, with_ties = FALSE) %>%  # Select the closest earlier survey slot
  ungroup() %>%
  mutate(n_survey = (day_relative - 1) * 8 + interval) %>%
  arrange(id, start_date)

# Create full dataset with all 56 surveys per participant
full_dat <- expand.grid(
  id = unique(prep_dat$id),
  n_survey = 1:56
) %>%
  left_join(time_dat, by = c("id", "n_survey")) %>% 
  select(-c(`Start Day`:time_diff))

# Create wide dataset 
wide_dat <- full_dat %>%
  pivot_wider(
    id_cols = c(id, trait_e:trait_c), 
    names_from = n_survey,  
    values_from = c(state_e:sociability)
  )

# two rows for each participants were created since the traits had missing data due to the rows with missings
wide_dat <- wide_dat %>% filter(!is.na(trait_e)) # remove the unncessary rows 


# information on the distance between observations
gaps <- schedule %>%
  pivot_longer(cols = Monday:Sunday, names_to = "day", values_to = "time") %>%
  arrange(day, Interval) %>%
  group_by(day) %>%
  mutate(
    time_minutes = as.numeric(as_hms(time)) / 60,
    gap = time_minutes - lag(time_minutes)
  ) %>%
  filter(!is.na(gap))

# Summary across all days and gaps
gaps %>%
  ungroup() %>%
  summarise(
    mean_gap = mean(gap),
    sd_gap = sd(gap),
    min_gap = min(gap),
    max_gap = max(gap)
  )

# Find the 30 consecutive time points with the least missing data  --------
# First, how common are partially missing observations?

partly_missing <- dat %>% 
  filter(is.na(state_e) | is.na(state_n) | is.na(state_o) | is.na(state_a) | 
           is.na(state_c) | is.na(happiness) | is.na(selfesteem) | is.na(duty) | 
           is.na(intellect) | is.na(adversity) | is.na(mating) | is.na(positivity) | 
           is.na(negativity) | is.na(deception) | is.na(sociability)) %>% 
  rowwise() %>% 
  mutate(n_missing = sum(c_across(state_e:sociability) |> is.na())) %>% 
  ungroup() %>% 
  group_by(SID) %>% 
  summarise(
    n_obs = n(),
    avg_missing = mean(n_missing)
  ) %>% 
  filter(n_obs > 0)

nrow(partly_missing) # 130 participants have at least one partially missing observation
mean(partly_missing$n_obs) # on average, participants have 2.5 partially missing observations
mean(partly_missing$avg_missing) # on average, 1.2 items are missing 

# Since most surveys were completed either entirely or with only a very small number of missing items, it seems reasonable to use only a single variable of the dataset to inspect the completion rates
# find the response rate across all observations for state_e using the wide dataset
response_rates <- wide_dat %>%
  select(starts_with("state_e_")) %>%
  summarise(across(everything(), ~ mean(!is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "survey", values_to = "response_rate") %>% 
  select(response_rate) %>% pull()

# Now find the 32 consecutive time points with the least missing data for state_e
window_size <- 32
threshold_min <- 0.60   # set minimum acceptable completion rate

res <- do.call(rbind, lapply(1:(length(response_rates) - window_size + 1), function(i) {
  w <- response_rates[i:(i + window_size - 1)]
  w_sorted <- sort(w)
  
  data.frame(
    start = i,
    end = i + window_size - 1,
    mean = mean(w),
    min = min(w),
    second_min = w_sorted[2],
    sd = sd(w)
  )
}))

# Keep only windows whose minimum completion rate is above the threshold
res_ok <- subset(res, min >= threshold_min)

if (nrow(res_ok) == 0) {
  message("No 32-wave window meets the minimum threshold.")
} else {
  # Among eligible windows, choose the one with the highest mean
  best <- res_ok[which.max(res_ok$mean), ]
  print(best)
  
  # Optional: inspect the selected response rates
  best_window_values <- response_rates[best$start:best$end]
  print(best_window_values)
} 

names(wide_dat)
# Only keep the observations from the selected window from time 3 to 34
final_dat <- wide_dat %>% 
  select(
    id, trait_e:trait_c, 
    matches("^(state_e|state_n|state_o|state_a|state_c|happiness|selfesteem|duty|intellect|adversity|deception|sociability|positivity|negativity|mating)_(1[0-9]|2[0-9]|3[0-4]|[3-9])$")
  )

names(final_dat)
# # Change time ending such that it goes from 1 to 30 (instead of 3 to 32)
# cn <- colnames(final_dat)
# 
# # extract trailing number
# nums <- as.numeric(sub(".*_(\\d+)$", "\\1", cn))
# 
# # identify names that end in 3:32
# to_change <- !is.na(nums) & nums >= 3 & nums <= 32
# 
# # replace only those names
# cn[to_change] <- paste0(sub("_(\\d+)$", "", cn[to_change]), "_", nums[to_change] - 2)
# 
# colnames(final_dat) <- cn
# names(final_dat)

write_csv(final_dat, here("data", "wide_dat_S1.csv"))


### Taken from other script, check whether I need something later

coverage_extra <- lavInspect(extra_m1.fit_m, "coverage") %>% data.frame()
cor_table_lower_empty <- function(data ){
  table1 <- data
  table1[upper.tri(table1, diag=FALSE)] <-""
  table1 <- matrix(as.numeric(unlist(table1)),nrow=nrow(table1))
  # table1 <- as.data.frame(table1)
}
coverage <- cor_table_lower_empty(coverage_extra) %>% as.data.frame()
names(coverage) <- paste0("T",1:56)
writexl::write_xlsx(coverage, "coverage_S1.xlsx")
n_low_coverage <- coverage_extra[coverage_extra < .4] # 102 cases

completion_rate_item <- lavInspect(extra_m1_.fit, "coverage") %>% diag()

completion_rate_item [missings_per_item < .75] %>% length()


tdat <- dat %>% select(starts_with("state_e"))
n_missings <- sum(is.na(tdat))
total_obs <- 204*56
n_missings/total_obs
