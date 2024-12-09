### Treating large models output from dynamo
library(tidyverse)
library(data.table)
library(tresthor)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(RcppXPtrUtils)
source("largemodel/translate_model_to_tresthor.R")

size = "large"
model <-translate_modelprg(modfile = str_c("largemodel/3me_",size,".prg"),
                   calibfile =  str_c("largemodel/calib_",size,".csv"),
                   base.year = 2019, last.year = 2050

                   )

tresthor::create_model(model_name = size ,
                       endogenous = model$endo,
                       exogenous = model$exo ,
                       coefficients = model$coef,
                       equations = model$equations,
                       rcpp = FALSE ,rcpp_path = "src/R_solver_files/" ,
                       no_var_map = TRUE )

# model$equations |> write.csv("largemodel/equations.csv")

create_model("opale",model_source = "inst/Opale/opale.txt",rcpp = TRUE,use.superlu = TRUE)
# opale@prologue_equations_f |> View()
## fonction d'évaluation des équations dans R

dat_base <- readRDS("inst/Opale/donnees_opale.rds")
coef <-   readRDS("inst/Opale/coefficients_opale.rds")
data <- add_coeffs(coef,database = dat_base,pos.coeff.name = 2, pos.coeff.value = 1)

source("R/equations_function.R")
sourceCpp("largemodel/computejacobian.cpp")
time = 80

pliop <- equations_r(t = 90, data)


mat_dat <- as.matrix
JacobianMatrix <- computeJacobian(equations_r, as.matrix(data_full) )

print(JacobianMatrix)


Rcpp::XPTr

f
