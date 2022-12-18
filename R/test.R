library(tidyr)
library(tidyverse)
library(lme4)
## Generate bootstrap samples
source("estimation_functions.R")

B <- 100
n <- nrow(epilepsy)

model <- run_model(epilepsy,"epilepsy")

bootstrap <- tibble(B = 1:B) %>%
  crossing(epilepsy) %>%
  group_by(B) %>%
  summarize(Index = sample(n, n, replace = TRUE),
            seizures = seizures[Index],
            age = age[Index],
            expind = expind[Index],
            treat = treat[Index],
            id = id[Index],
            .groups = "drop") %>%
  nest_by(B) %>% summarize(
    merge(enframe(run_model(data,"epilepsy")$beta,value='estimate'),
          enframe(run_model(data,"epilepsy")$test_stat,value='statistic')),
    .groups = "drop")

bootstrap %>%
  filter(name == "age") %>%
  ggplot(aes(x = estimate)) +
  geom_density()
## Compute bootstrap standard error and confidence interval
bootstrap %>%
  filter(name == "expind") %>%
  summarize(SD = sd(estimate),
            Lower95 = quantile(estimate, .025),
            Upper95 = quantile(estimate,.975),
            mean = mean(estimate))


t_value <- model$test_stat["expind:treat"]

bootstrap_1 <- tibble(B = 1:B) %>%
  crossing(epilepsy) %>%
  group_by(B) %>%
  summarize(Index1 = sample(n, n, replace = TRUE),
            Index2 = sample(n, n, replace = TRUE),
            seizures = seizures[Index1],
            age = age[Index1],
            expind = expind[Index1],
            treat = treat[Index2],
            id = id[Index1],
            .groups = "drop") %>%
  nest_by(B) %>%
  summarize(
    merge(enframe(run_model(data,"epilepsy")$beta,value='estimate'),
          enframe(run_model(data,"epilepsy")$test_stat,value='statistic')),
    .groups = "drop")

## Compute p-value
bootstrap_1 %>%
  filter(name == "expind:treat") %>%
  summarize(p = mean(abs(statistic) > abs(t_value)))
