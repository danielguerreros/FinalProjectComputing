source("estimation_functions.R")
#' Confidence intervals for predictors
#'
#' @param data the data set used to fit the model
#' @param predictor the predictor whose confidence interval we want to calculate
#' @param n_boot the number of bootstraping samples to use
#'
#' @return a dataframe with the 95% confidence interval
#' @export
#'
#' @examples
#' ##### Epilepsy #####
#'
#' ## Load data
#' epilepsy <- read_csv("epilepsy.csv")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' confidenceInterval_age <- confidence_interval(epilepsy, "age",1000)
confidence_interval <- function(data, predictor,n_boot){
  if (!predictor %in% c("age","expind","treat","expind:treat"))
       stop("You must set predictor to one of: age,expind,treat,expind:treat")


  B <- n_boot
  n <- nrow(data)

  #Fit initial model
  model_0 = run_model(data,"epilepsy")

  bootstrap <- tibble(B = 1:B) %>%
    crossing(data) %>%
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
    filter(name == predictor) %>%
    summarize(Lower95 = quantile(estimate, .025),
              Upper95 = quantile(estimate,.975),
              stderror = sd(estimate),
              variable = predictor)

}

#' P Value of the predictor age
#'
#' @param data the data set used to fit the model
#' @param n_boot the number of bootstraping samples to use
#' @param t_value statistic of predictor when the model was fitted
#'
#' @return The pvalue of the hypothesis predictor age equal to 0
#' @export
#'
#' @examples
#' ##### Epilepsy #####
#'
#' ## Load data
#' epilepsy <- read_csv("epilepsy.csv")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' pValue_age <- pValue_age(epilepsy,1.16,1000)
pValue_age <- function(data,t_value,n_boot){

  B <- n_boot


  bootstrap <- tibble(B = 1:B) %>%
    crossing(data) %>%
    group_by(B) %>%
    summarize(Index1 = sample(n, n, replace = TRUE),
              Index2 = sample(n, n, replace = TRUE),
              seizures = seizures[Index1],
              age = age[Index2],
              expind = expind[Index1],
              treat = treat[Index1],
              id = id[Index1],
              .groups = "drop") %>%
    nest_by(B) %>%
    summarize(
      merge(enframe(run_model(data,"epilepsy")$beta,value='estimate'),
            enframe(run_model(data,"epilepsy")$test_stat,value='statistic')),
      .groups = "drop")

  ## Compute p-value
  bootstrap %>%
    filter(name == "age") %>%
    summarize(p = mean(abs(statistic) > abs(t_value)))
}

#' P Value of the predictor expind
#'
#' @param data the data set used to fit the model
#' @param n_boot the number of bootstraping samples to use
#' @param t_value statistic of predictor when the model was fitted
#'
#' @return The pvalue of the hypothesis predictor expind equal to 0
#' @export
#'
#' @examples
#' ##### Epilepsy #####
#'
#' ## Load data
#' epilepsy <- read_csv("epilepsy.csv")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' pValue_expind <- pValue(epilepsy,1.16,1000)
pValue_expind <- function(data,t_value,n_boot){

  B <- n_boot


  bootstrap <- tibble(B = 1:B) %>%
    crossing(data) %>%
    group_by(B) %>%
    summarize(Index1 = sample(n, n, replace = TRUE),
              Index2 = sample(n, n, replace = TRUE),
              seizures = seizures[Index1],
              age = age[Index1],
              expind = expind[Index2],
              treat = treat[Index1],
              id = id[Index1],
              .groups = "drop") %>%
    nest_by(B) %>%
    summarize(
      merge(enframe(run_model(data,"epilepsy")$beta,value='estimate'),
            enframe(run_model(data,"epilepsy")$test_stat,value='statistic')),
      .groups = "drop")

  ## Compute p-value
  bootstrap %>%
    filter(name == "expind") %>%
    summarize(p = mean(abs(statistic) > abs(t_value)))
}

#' P Value of the predictor treat
#'
#' @param data the data set used to fit the model
#' @param n_boot the number of bootstraping samples to use
#' @param t_value statistic of predictor when the model was fitted
#'
#' @return The pvalue of the hypothesis predictor treat equal to 0
#' @export
#'
#' @examples
#' ##### Epilepsy #####
#'
#' ## Load data
#' epilepsy <- read_csv("epilepsy.csv")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' pValue_treat <- pValue(epilepsy,1.16,1000)
pValue_treat <- function(data,t_value,n_boot){

  B <- n_boot


  bootstrap <- tibble(B = 1:B) %>%
    crossing(data) %>%
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
  bootstrap %>%
    filter(name == "treat") %>%
    summarize(p = mean(abs(statistic) > abs(t_value)))
}


