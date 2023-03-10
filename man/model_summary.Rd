% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_functions.R
\name{model_summary}
\alias{model_summary}
\title{Model summary}
\usage{
model_summary(data, n_boot)
}
\arguments{
\item{data}{the data set used to fit the model}

\item{n_boot}{the number of bootstraping samples to use}
}
\value{
A dataframe with important information about the predictors like point estimates, tstatistic, standard error, confidence intervals and p values.
}
\description{
Model summary
}
