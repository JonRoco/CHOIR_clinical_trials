#Load packages
library(survival)
library(tidyverse)

# Set Parameters and Seed:
set.seed(0)
n <- 185
n_trials <- 1e4
HR <- log(0.5)/log(0.65)

# Define Hazard Rates for Groups A and B.
# Both groups have the same hazard rate of 0.431 initially, which is used to simulate the survival times.
h1 <- -log(0.65)  #65% PFS
h2 <- -log(0.65) #65% PFSa


# Simulate exponetial survival times with hazard rate h1 and h2
timeA <- replicate(n_trials, rexp(n, rate = h1)) %>% as_tibble()
timeB <- replicate(n_trials, rexp(n, rate = h2)) %>% as_tibble()

surv_time <- bind_rows(timeA, timeB)

#' Simulate time on study. 
#' Assume total recruitment time = 3 and analysis will be done when the final patient has follow-up = 1.

tosA <- replicate(n_trials, runif(n, 1, 4)) %>% as_tibble()
tosB <- replicate(n_trials, runif(n, 1, 4)) %>% as_tibble()

study_time <- bind_rows(tosA,tosB)

# Determine Censoring (0) or Event (1):
#' If study_time > surv_time , then event = 1
#' If study_time <= surv_time, then event = 0
event_dat <- map2(study_time, surv_time,`>`) %>% 
  as_tibble() %>%
  mutate_if(is.logical,as.numeric)

# Calculate Observed Times (Censoring or Events):
time_dat <- map2(study_time, surv_time,`pmin`) %>% as_tibble() 
group <- rep(c("A","B"), each = n)

# Fit Cox Proportional Hazards Model and Calculate 95% CI
# Report the upper 95% CI of the HR estimate.
upperCI <- rep(0, n_trials)
j<- 0

for (i in colnames(time_dat)){
 j <- j+1
 x <- coxph(Surv(time_dat[[i]], event_dat[[i]]) ~ group) 
 upperCI[j] <- exp(confint(x))[2]
}

sum(upperCI < HR)/n_trials
