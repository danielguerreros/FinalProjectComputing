library(tidyr)
library(tidyverse)
library(lme4)
## Generate bootstrap samples
source("estimation_functions.R")

B <- 100
n <- nrow(epilepsy)

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
  nest_by(B) %>% summarize(enframe((run_model(data, "epilepsy")$beta)), .groups = "drop")

## Compute bootstrap standard error and confidence interval
bootstrap %>%
  filter(name == "expind") %>%
  summarize(SD = sd(value),
            Lower95 = quantile(value, .025),
            Upper95 = quantile(value,.975))

