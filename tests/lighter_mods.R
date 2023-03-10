library(tidyverse)
library(tresthor)
library(Matrix)
library(Deriv)


### Create model light (more efficient version)
model_name = "thremee"
model_source = "tests/threeme_8x8_thor.txt"
algo = TRUE
rcpp = FALSE
rcpp_path = "tests"
endogenous = NULL
exogenous = NULL
coefficients = NULL
equations = NULL
env= globalenv()
no_var_map = FALSE

################################
#### 0.a  Arguments verification
################################
cat("Let's build the model :)
      \n
  Step 1: Running various checks on the inputs... \n\n")

assertthat::is.string(model_name)
if (!is.null(model_source)) {
  assertthat::is.readable(model_source)
  tresthor:::check_model_file(model_source)
  cat(paste0("The model will be created using the file ",basename(model_source),".\nIgnoring all other equations specifications."))

  manual_model_input = FALSE}

if (is.null(model_source)) {
  manual_model_input = TRUE

  if (!is.null(endogenous) & !is.null(equations)) {
    assertthat::assert_that(class(endogenous) == "character")
    assertthat::assert_that(class(equations) == "character")
    if (!is.null(exogenous)){assertthat::assert_that(class(exogenous) == "character")}
    if (!is.null(coefficients)){assertthat::assert_that(class(coefficients) == "character")}
  }else{
    stop(paste0("No model source is provided and no equation is specified.
Please either provide an external model file as model_source, or provide the equations and their specifications
(endogenous variables, and if needed exogenous variables and coefficients) manually in the appropriate arguments.
Refer to the documentation for more information."))
  }
}
assertthat::assert_that(class(algo) == "logical")
assertthat::assert_that(class(rcpp) == "logical")

if (rcpp){assertthat::is.dir(rcpp_path)}
################################
#### 0.b Setup
################################

options(stringsAsFactors=FALSE)
delta <- function(n,x, y=TRUE) NULL
drule[["delta"]] <- alist(x=1, y=NULL) # y is just a logical
drule[["abs"]] <- alist(x=ifelse(x==0, 0, sign(x)))

################################
#### 1.a Reading from outside source
################################
if ( manual_model_input == FALSE ){
  model_input <- readLines(model_source,warn = FALSE)
  eqlist <- unique(model_input[11:length(model_input)])
  endo <- tolower(unique(strsplit(model_input[2],split = ",")[[1]]))
  exo <- tolower(unique(strsplit(model_input[5],split = ",")[[1]]))
  coeff <-tolower(unique(strsplit(model_input[8],split = ",")[[1]]))
}
################################
#### 1.b Reading from manual input
################################
if ( manual_model_input == TRUE ){
  eqlist<-unique(equations)
  endo <- tolower(unique(endogenous))
  exo <- tolower(unique(exogenous))
  coeff <- tolower(unique(coefficients))
}

################################
#### 1.b Reading from manual input
################################
endo<-endo[order(endo)]
exo<-exo[order(exo)]
coeff<-coeff[order(coeff)]

################################
#### 1.c Checks for inputs
################################
eqlist   <- gsub("mylg\\(","lag\\(",eqlist)
###checks that all equations contain (=)
tresthor:::check_equation_input(eqlist)

###checks that all elements in variable vectors are viable
tresthor:::check_var_vector(endo,"Endogenous variables")
tresthor:::check_var_vector(exo,"Exogenous variables")
tresthor:::check_var_vector(coeff,"Coefficient variables")

tresthor:::check_variable_conflict(endo , exo)
tresthor:::check_variable_conflict(endo , coeff)
tresthor:::check_variable_conflict(exo , coeff)

################################
#### 2.a Equations indexing and splitting
################################
exo  <- tresthor:::is_in_formulas(exo,eqlist,"exogenous")
endo <- tresthor:::is_in_formulas(endo,eqlist,"endogenous")
coeff<- tresthor:::is_in_formulas(coeff,eqlist,"coefficients")
all_model_variables <- c(endo,exo, coeff)
all_model_variables <- all_model_variables[order(all_model_variables)]
##from here on, steps should be written with much minimalism in order to be used easily with model modifier functions
##A/creating indexes names, and LHS RHS
equations_list<-tresthor:::create_equations_list(eqlist)


if(tresthor:::parser_lag_delta_check(equations_list$equation)==FALSE){stop("Parser error on the lags and/or delta. Please check the model's formula specification.")}

memory  <- list(eqlist_raw = eqlist ,
                endo_raw = endo ,
                exo_raw = exo ,
                coeff_raw = coeff)

################################
#### 2.b generating the matrix that houses the contemporaneous endo
################################
cat("\n
  Step 2: Identifying the endogenous variables in the equations...
      \n")
##B/ generating eqns table

eqns <- tresthor:::table_contemporaneous_endos(formula_list = equations_list$formula,endogenous = endo,exogenous = exo, coefflist = coeff,equations_index = equations_list$id)


################################
#### 2.c checking parser and inputs
################################
tresthor:::check_eq_var_identification(endo = endo, eqns =eqns,names_of_equations = equations_list$name)

################################
#### 3.a decomposition
################################
if (algo == TRUE){
  cat("\n
  Step 3: Decomposing the model into blocks...
      \n")}else{
        cat("\n
  Step 3: Creating the unique block of the model...
      \n")
      }
## checks
if(algo == FALSE & nrow(eqns) > 75){cat("The model specified appears to be large. Consider using algo == TRUE to improve solver performances.")}
if(algo == TRUE  & nrow(eqns) < 50){cat("The model specified appears to be small. Consider using algo == FALSE to improve solver performances.")}

decomposition<-tresthor:::decomposing_model(endogenous_variables = endo,eq_var_matrix = eqns,decomposition = algo)
invisible(list2env(decomposition,environment()))

equations_list$part<-"tbd"
equations_list$part<-ifelse(equations_list$id %in% prologue_equations , "prologue",equations_list$part )
equations_list$part<-ifelse(equations_list$id %in% heart_equations , "heart",equations_list$part )
equations_list$part<-ifelse(equations_list$id %in% epilogue_equations , "epilogue",equations_list$part )

if("tbd" %in% equations_list$part){
  stop("Some equations were not identified in the decomposition.")
}else{cat("Decomposition algorithm successful. Below is the size of each block : \n")
  print(table(equations_list$part))}


