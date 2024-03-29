% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_0_model_modifications.R
\name{model_endo_exo_switch}
\alias{model_endo_exo_switch}
\title{Endo-Exo switch
Creates a new model base on an existing one with some endos and exos switched}
\usage{
model_endo_exo_switch(
  base_model,
  new_model_name,
  new_endo,
  new_exo,
  algo = TRUE,
  rcpp = FALSE,
  env = globalenv(),
  rcpp_path = getwd()
)
}
\arguments{
\item{base_model}{: the reference model to be used (obviously a thoR.model)}

\item{new_model_name}{: character string, must be different from base model name}

\item{new_endo}{: vector of  exogenous variables from base model to be switched to endo}

\item{new_exo}{: vector of  endogenous variables from base model to be switched to exo. length must match new_endo.}

\item{algo}{boolean. TRUE if using decomposition of model (recommended for large model). Default : TRUE}

\item{rcpp}{boolean. TRUE if using Rcpp for solver (recommended for large model). Default : FALSE}

\item{env}{Environment where to store the model object. Default : globalenv()}

\item{rcpp_path}{path to directory where to store the rcpp source files for the model. Default : working directory}
}
\value{
a thor model in the global environment
}
\description{
Endo-Exo switch
Creates a new model base on an existing one with some endos and exos switched
}
