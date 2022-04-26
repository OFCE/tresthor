#' Single thoR equation.
#' Lone equations can be solved for a given endogenous,
#' and the coefficients may be estimated using various econometric methods if correctly parsed.
#'
#' @slot equation_name character string. Name of the equation.
#' @slot formula character string. Formula of the equation as originally inputted.
#' @slot new_formula character string. Formula that will be used in the solver.
#' @slot coefflist character vector. List of coefficients.
#' @slot endogenous character vector of length 1. endogenous variable of the equation.
#' @slot exogenous character vector. exogenous variables of the equation.
#' @slot LHS character string of the left hand side of the equation for use with econometric estimation functions.
#' @slot RHS character string of the left hand side of the equation for use with econometric estimation functions.
#' @slot maxlag numeric. maximum lag of the equation (useful for determining number of observations.)
#' @slot jacobian character string. derivative of the function new formula.
#' @slot eq_f function of evaluation of the new formula used in the solver.
#' @slot jac_f function of evaluation of the jacobian used in the solver.
#' @slot econometric_safe boolean. If the LHS and RHS can be used easily for econometric estimations(beta mode).
#'
#' @export
thoR.equation<-setClass("thoR.equation",
                        slots = c(equation_name="character",
                                  formula="character",
                                  new_formula="character",
                                  coefflist="vector",
                                  endogenous="character",
                                  exogenous = "character",
                                  LHS="character",
                                  RHS="character",
                                  maxlag="numeric",
                                  jacobian="matrix",
                                  eq_f="function",
                                  jac_f="function",
                                  econometric_safe="logical"))

#' thoR models with multiple equations
#' Models are created using the create_model function
#' @export
thoR.model<-setClass("thoR.model",
           slots = c(model_name="character" ,
                     equation_list="data.frame",
                     var_map="list",
                     variables_matrix="data.frame",

                     endo_list="vector",
                     exo_list="vector",
                     coeff_list="vector",

                     prologue="logical",
                     heart="logical",
                     epilogue="logical",

                     prologue_jacobian="matrix",
                     heart_jacobian="matrix",
                     epilogue_jacobian="matrix",

                     prologue_endo="vector",
                     heart_endo="vector",
                     epilogue_endo="vector",

                     algo="logical",
                     rcpp="logical",
                     rcpp_source="character",

                     prologue_equations_f="function",
                     heart_equations_f="function",
                     epilogue_equations_f="function",
                     prologue_jacobian_f="function",
                     heart_jacobian_f="function",
                     epilogue_jacobian_f="function"
                     ))
