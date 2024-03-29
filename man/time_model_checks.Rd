% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3_1_data_checks.R
\name{time_model_checks}
\alias{time_model_checks}
\title{Checking the viability of a variable as a time index
Condition : must have no duplicates, no NAs and sorted in ascending order}
\usage{
time_model_checks(database, index_time)
}
\arguments{
\item{database}{dataframe the database to check}

\item{index_time}{character string : name of the variable to use as time_index}
}
\value{
logical TRUE if the time_index chosen can be used for the solver
}
\description{
Checking the viability of a variable as a time index
Condition : must have no duplicates, no NAs and sorted in ascending order
}
\keyword{internal}
