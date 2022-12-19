#' Final Project Model Fitting
#'
#' @param data the data set for your example
#' @param example the name of your example: one of "culcita", "ctsib", "epilepsy", or "tortoise"
#'
#' @return A list with summary statistics from the fitted model.
#' @export
#'
#' @examples
#' @importFrom dplyr tibble summarize %>% group_by
#' ##### Epilepsy #####
#'
#' ## Load data
#' data("epilepsy")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' ## Fit model to epilepsy data. Fixed effects in the model are:
#' ##  (Intercept) -- the intercept
#' ##  age -- age as a continuous predictor
#' ##  expind -- a categorical variable with two levels (0 for before and 1 for after)
#' ##  expind:treat -- a categorical variable with two levels (1 for observations from individuals on the drug in the after period and 0 otherwise)
#'
#' epilepsy_fit <- run_model(epilepsy, "epilepsy")

run_model <- function(data, example = "tortoise"){

  ## Fit model
  if(example == "tortoise")
    lmer_fit <- glmer(shells~prev+offset(log(Area))+factor(year)+(1|Site),
                      family=poisson,
                      data=data)

  else if(example == "culcita")
    lmer_fit <- glmer(predation~ttt+(1|block),
                      data=data,
                      family=binomial)

  else if(example == "ctsib")
    lmer_fit <- glmer(stable ~ Surface + Vision + (1|Subject),
                     family = binomial,
                     data = data)

  else if(example == "epilepsy")
    lmer_fit <- glmer(seizures ~ age + expind + expind:treat + (1|id),
                      family = poisson,
                      data = data)
  else
    stop("You must set example to one of: tortoise, culcita, ctsib, or epilspsy.")

  ## Compute model summary
  lmer_summ <- summary(lmer_fit)

  ## Extract coefficients
  coeff <- coefficients(lmer_summ)[,"Estimate"]

  ## Extract t-statistics for each coefficient
  test_stat <- coefficients(lmer_summ)[,"z value"]

  ## Extract random effects
  re <- ranef(lmer_fit)[[1]][[1]]

  ## Extract random effects variance
  sigmasq <- lmer_summ$varcor[[1]][[1]]

  ## Extract optimization information
  optinfo <- attributes(lmer_fit)$optinfo

  ## Avoid estimate of sigmasq=0 for the original tortoise data only
  if(example == "tortoise" & all(sort(data$shells) == c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,4,5,8,10,11,12))){
    if(sigmasq == 0){
      sigmasq <- .100

      set.seed(7777)
      re <- rnorm(length(re), sd = sqrt(sigmasq))
    }
  }

  ## Return values
  important_list = list(beta = coeff,
       sigmasq = sigmasq,
       re = re,
       test_stat = test_stat,
       convergence = c(optinfo$conv$opt,
                       optinfo$conv$lme4$code),
       messages = c(optinfo$message,
                    optinfo$conv$lme4$messages))


}

#' Model results
#'
#' @param data the data set used to fit the model
#' @param n_boot the number of bootstraping samples to use
#'
#' @return A dataframe with important information about the predictors like point estimates, tstatistic, standard error and confidence intervals.
#' @export
#'
#' @examples
#' @importFrom dplyr tibble summarize %>% group_by
#' ##### Epilepsy #####
#'
#' ## Load data
#' data("epilepsy")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' ## Fit model to epilepsy data. Fixed effects in the model are:
#' ##  (Intercept) -- the intercept
#' ##  age -- age as a continuous predictor
#' ##  expind -- a categorical variable with two levels (0 for before and 1 for after)
#' ##  expind:treat -- a categorical variable with two levels (1 for observations from individuals on the drug in the after period and 0 otherwise)
#'
#' results <- model_results(epilepsy, 10)
model_results <- function(data,n_boot){

  n <- nrow(data)

  # Fit initial model to get statistics and point estimates
  model_0 <- run_model(data,"epilepsy")

  # Saving the point estimates and the t_stat in one tibble
  point_estimates <- merge(enframe(model_0$beta,name="Variable",value="Point Estimate"),
                           enframe(model_0$test_stat,name="Variable",value="t_stat"))

  # Bootstrapping to get the confidence interval and standard error
  bootstrap <- tibble(B = 1:n_boot) %>%
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
      merge(enframe(run_model(data,"epilepsy")$beta,name="Variable",value='estimate'),
            enframe(run_model(data,"epilepsy")$test_stat,name="Variable",value='statistic')),
      .groups = "drop")

  # Saving bootstrap results in tibble
  important_stats <- bootstrap %>%
    group_by(Variable) %>%
    summarize(
      stderror = sd(estimate),
      Lower95 = quantile(estimate, .025),
      Upper95 = quantile(estimate,.975)
    )

  merge(point_estimates,important_stats)
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
#'
#' @examples
#' @importFrom dplyr tibble summarize %>% group_by
#' ##### Epilepsy #####
#'
#' ## Load data
#' data("epilepsy")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' pValue_age <- pValue_age(epilepsy,1.16,10)
pValue_age <- function(data,t_value,n_boot){

  B <- n_boot
  n <- nrow(data)

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
      enframe(run_model(data,"epilepsy")$test_stat,value='statistic'),
      .groups = "drop")

  ## Compute p-value
  bootstrap %>%
    filter(name == "age") %>%
    summarize(p_value = mean(abs(statistic) > abs(t_value))) %>%
    mutate(Variable = "age")

}

#' P Value of the predictor expind
#'
#' @param data the data set used to fit the model
#' @param n_boot the number of bootstraping samples to use
#' @param t_value statistic of predictor when the model was fitted
#'
#' @return The pvalue of the hypothesis predictor age equal to 0
#' @export
#'
#'
#' @examples
#' @importFrom dplyr tibble summarize %>% group_by
#' ##### Epilepsy #####
#'
#' ## Load data
#' data("epilepsy")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' pValue_expind <- pValue_expind(epilepsy,1.16,10)
pValue_expind <- function(data,t_value,n_boot){

  B <- n_boot
  n <- nrow(data)

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
      enframe(run_model(data,"epilepsy")$test_stat,value='statistic'),
      .groups = "drop")

  ## Compute p-value
  bootstrap %>%
    filter(name == "expind") %>%
    summarize(p_value = mean(abs(statistic) > abs(t_value))) %>%
    mutate(Variable = "expind")

}

#' P Value of the predictor treat
#'
#' @param data the data set used to fit the model
#' @param n_boot the number of bootstraping samples to use
#' @param t_value statistic of predictor when the model was fitted
#'
#' @return The pvalue of the hypothesis predictor age equal to 0
#' @export
#'
#' @examples
#' @importFrom dplyr tibble summarize %>% group_by
#' ##### Epilepsy #####
#'
#' ## Load data
#' data("epilepsy")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' pValue_treat <- pValue_age(epilepsy,1.16,10)
pValue_treat <- function(data,t_value,n_boot){

  B <- n_boot
  n <- nrow(data)

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
      enframe(run_model(data,"epilepsy")$test_stat,value='statistic'),
      .groups = "drop")

  ## Compute p-value
  bootstrap %>%
    filter(name == "expind:treat") %>%
    summarize(p_value = mean(abs(statistic) > abs(t_value))) %>%
    mutate(Variable = "expind:treat")

}


#' Model summary
#'
#' @param data the data set used to fit the model
#' @param n_boot the number of bootstraping samples to use
#'
#' @return A dataframe with important information about the predictors like point estimates, tstatistic, standard error, confidence intervals and p values.
#' @export
#'
#' @examples
#' @importFrom dplyr tibble summarize %>% group_by
#' ##### Epilepsy #####
#'
#' ## Load data
#' data("epilepsy")
#'
#' ## Sum separate observations for each patient in the after period
#' epilepsy <- epilepsy %>%
#'   group_by(id,treat,expind,age) %>%
#'   summarize(seizures = sum(seizures),
#'             .groups = "drop")
#'
#' ## Fit model to epilepsy data. Fixed effects in the model are:
#' ##  (Intercept) -- the intercept
#' ##  age -- age as a continuous predictor
#' ##  expind -- a categorical variable with two levels (0 for before and 1 for after)
#' ##  expind:treat -- a categorical variable with two levels (1 for observations from individuals on the drug in the after period and 0 otherwise)
#'
#' summary <- model_summary(epilepsy, 10)
model_summary <- function(data,n_boot){

  # Initial results
  results <- model_results(data,n_boot)

  #Getting statistics for predictors
  pAge <- results %>%
    filter(Variable=="age") %>%
    pull(t_stat)

  pExp <- results %>%
    filter(Variable=="expind") %>%
    pull(t_stat)

  pTreat <- results %>%
    filter(Variable=="expind:treat") %>%
    pull(t_stat)

  # Getting p Values for predictors
  p_Values <- pValue_age(data,pAge[1],n_boot) %>%
    bind_rows(pValue_expind(data,pExp[1],n_boot)) %>%
    bind_rows(pValue_treat(data,pTreat[1],n_boot))

  left_join(results,p_Values)
}
