
#' Title
#'
#' @param data
#' @param predictor
#' @param n_boot
#'
#' @return
#' @export
#'
#' @examples
confidence_interval <- function(data, predictor,n_boot){
  if (!data %in% c("age","expind","treat"))
       stop("You must set predictor to one of: age,expind,treat")
  B <- n_boot
  n <- nrow(data)

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
              Upper95 = quantile(estimate,.975))
}

#' Title
#'
#' @param data
#' @param n_boot
#'
#' @return
#' @export
#'
#' @examples
pValue <- function(data, n_boot){

  B <- n_boot
  model_0 <- run_model(epilepsy,"epilepsy")
  t_value <- model$test_stat["expind:treat"]


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
    filter(name == "expind:treat") %>%
    summarize(p = mean(abs(statistic) > abs(t_value)))
}
