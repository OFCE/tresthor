### Treating large models output from dynamo
library(tidyverse)
library(data.table)
library(tresthor)
source("largemodel/translate_model_to_tresthor.R")
size = "large"

model <-translate_modelprg(modfile = str_c("largemodel/3me_",size,".prg"),
                   calibfile =  str_c("largemodel/calib_",size,".csv"),
                   base.year = 2019, last.year = 2050

                   )

tresthor::create_model(model_name = "large" ,
                       endogenous = model$endo,
                       exogenous = model$exo ,
                       coefficients = model$coef,
                       equations = model$equations,
                       rcpp = FALSE ,rcpp_path = "src/R_solver_files/" ,
                       no_var_map = TRUE )
