% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_6_dico_builder.R
\name{build_dico}
\alias{build_dico}
\title{Building the variables map of the model}
\usage{
build_dico(equations_list_df, endos, exos, coeffs, p_endos, e_endos, h_endos)
}
\arguments{
\item{equations_list_df}{equations list dataframe}

\item{endos}{full endos vector}

\item{exos}{exo vector}

\item{coeffs}{coefficient vector}

\item{p_endos}{prologue endo vector}

\item{e_endos}{epilogue endo vector}

\item{h_endos}{heart endo vector}
}
\value{
list that gives info on each variables
}
\description{
Building the variables map of the model
}
\keyword{internal}
