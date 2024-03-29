% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_thor_classes.R
\docType{class}
\name{thoR.equation-class}
\alias{thoR.equation-class}
\alias{thoR.equation}
\title{Single thoR equation.
Lone equations can be solved for a given endogenous,
and the coefficients may be estimated using various econometric methods if correctly parsed.}
\description{
Single thoR equation.
Lone equations can be solved for a given endogenous,
and the coefficients may be estimated using various econometric methods if correctly parsed.
}
\section{Slots}{

\describe{
\item{\code{equation_name}}{character string. Name of the equation.}

\item{\code{formula}}{character string. Formula of the equation as originally inputted.}

\item{\code{new_formula}}{character string. Formula that will be used in the solver.}

\item{\code{coefflist}}{character vector. List of coefficients.}

\item{\code{endogenous}}{character vector of length 1. endogenous variable of the equation.}

\item{\code{exogenous}}{character vector. exogenous variables of the equation.}

\item{\code{LHS}}{character string of the left hand side of the equation for use with econometric estimation functions.}

\item{\code{RHS}}{character string of the left hand side of the equation for use with econometric estimation functions.}

\item{\code{maxlag}}{numeric. maximum lag of the equation (useful for determining number of observations.)}

\item{\code{jacobian}}{character string. derivative of the function new formula.}

\item{\code{eq_f}}{function of evaluation of the new formula used in the solver.}

\item{\code{jac_f}}{function of evaluation of the jacobian used in the solver.}

\item{\code{econometric_safe}}{boolean. If the LHS and RHS can be used easily for econometric estimations(beta mode).}
}}

