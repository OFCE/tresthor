% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_3_quick_plot.R
\name{quick_plot}
\alias{quick_plot}
\title{Draw a line plot quickly}
\usage{
quick_plot(
  variables,
  database = data,
  start = NULL,
  end = NULL,
  index_time = "date",
  growth_rate = FALSE,
  title = NULL,
  colours = NULL
)
}
\arguments{
\item{variables}{variables to plot}

\item{database}{data.frame}

\item{start}{first date to plot Default: NULL (first observation)}

\item{end}{last date to plot. Default: NULL (last observation)}

\item{index_time}{Name of the column which contains the time info. Default: "date"}

\item{growth_rate}{TRUE if growth rates should be plotted instead of levels. Default : FALSE}

\item{title}{Title of the graph}

\item{colours}{vector of colours for the plotted line, must be at least as many colours as variables. extra colours will be ignored. If not enough colours are provided, the colours will revert to the defaul ones.}
}
\value{
ggplot graph
}
\description{
Draw a line plot quickly
}
