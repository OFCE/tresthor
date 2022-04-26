#' check if external model input is ok
#'
#' @param modelefile path to file to check
#' @importFrom purrr map
#' @importFrom purrr %>%
#' @return checks
#' @keywords internal
#'
check_model_file<-function(modelefile){

  scanner <- function(n){scan(modelefile,character(),sep = ",",skip= n-1,nlines = 1,quiet = TRUE)}

  list_scanned<-  c(1:10) %>% purrr::map(~scanner(.x))

   format_good <- grep(pattern = "endo",list_scanned[[1]]) * grep(pattern = "exo",list_scanned[[4]]) * grep(pattern = "coeff",list_scanned[[7]])

   input_good <- is.character(list_scanned[[2]]) * (is.character(list_scanned[[5]]) + is.null(list_scanned[[5]]) ) * (is.character(list_scanned[[8]]) + is.null(list_scanned[[8]]) )

   if (format_good == 0){stop("Something is wrong with the format of the model_source file. Refer to model_source_example() to see an example of a correct model input")}else{cat("Format of the model appears to be ok. \n")}
   if (input_good == 0){stop("Something is wrong with the inputs the model_source file. Refer to model_source_example() to see an example of a correct model input")}else{cat("Inputs of the model appear to be ok. \n")}
   }


#' Basic equation checks
#'
#' @param equation_vector character vector to check
#' @importFrom stringr str_count
#' @return check
#' @keywords internal
check_equation_input<- function(equation_vector){
   ev <- equation_vector
   check_ok <- TRUE
   wrongcharacters <- c("@|\\$|%|#|~|\\||;")


   if(prod(grepl(pattern = "=",x = ev))==0){
      cat("The following equations do not have the = sign: \n")
      print(ev[grepl(pattern = "=",x = ev)==FALSE])
      check_ok <- FALSE}

   if(sum( (stringr::str_count(ev,"=") > 1) ) > 0){
      cat("The following equations have too many equal signs: \n")
      print(ev[stringr::str_count(ev,"=") > 1])
      check_ok <- FALSE
   }

   if(sum( (stringr::str_count(ev,":") > 1) ) > 0){
      cat("Columns are used to give names to equations. The following equations have too many : signs: \n")
      print(ev[stringr::str_count(ev,":") > 1])
      check_ok <- FALSE
   }

   if(sum(grepl(pattern = wrongcharacters,x = ev))>0){
      cat("The following equations contain some invalid characters: \n")
      print(ev[grepl(pattern =wrongcharacters,x = ev)==TRUE])
      check_ok <- FALSE}

   if (check_ok == FALSE){stop("Please correct the input.")
      }else{
      cat("\nBasic equations checks OK \n")
   }
}


#' Checks for vectors
#'
#' @param variable_vector character vecteor of variable names to check
#' @param name object of the vector that is being checked
#'
#' @return check
#' @keywords internal
#'
check_var_vector<-function(variable_vector,name="variables"){
   vv<-variable_vector
   check_ok <- TRUE

   if (prod(grepl(pattern="^[a-z]",x = vv)) == 0){
      cat(paste0("Something is amiss with ",name,". Variables must start with a letter. \n"))
      print(vv[grepl(pattern="^[a-z]",x = vv)==FALSE])
      check_ok <- FALSE}

   if (sum(grepl(pattern="[^a-zA-Z0-9_]",x = vv)) > 0){
      cat(paste0("Something is amiss with ",name,". Only aplhanumeric characters and underscores are allowed for variable names. \n"))
      print(vv[grepl(pattern="[^a-zA-Z0-9_]",x = vv)==TRUE])
      check_ok <- FALSE}

   if (check_ok == FALSE){stop(paste0("Please correct the ",name," vector."))
      }else{
      cat(paste0(name," vector checks OK \n"))}
}

#' Checking if there is a variable conflict
#'
#' @param var1 First set of variables to check
#' @param var2 Second set of variables to check
#'
#' @return checks conflict
#' @keywords internal
check_variable_conflict<-function(var1,var2){
   if(length(intersect(var1,var2))>0){
      cat(paste0("The following variables appear in ",deparse(substitute(endo))," and ",deparse(substitute(exo)),":"))
      print(intersect(var1,var2))
      stop(paste0("Please check the composition of those vectors"))}
   }


#' Checking if model is well identified
#'
#' @param endo endogenous variables as declared
#' @param eqns equation identification matrix
#' @param names_of_equations vector with equation names
#'
#' @importFrom stats na.omit
#' @keywords internal
#' @return check if identification successful
check_eq_var_identification<-function(endo,eqns,names_of_equations){
   ok_checks <- TRUE
   vector_endo<-na.omit(unique(as.vector(as.matrix(eqns))))
   if(length(setdiff(endo,vector_endo))>0){
      cat("The following declared endogenous variables have not been detected in the equations at t=0. \n")
      print(setdiff(endo,vector_endo))
      ok_checks <- FALSE
   }

   if(length(setdiff(vector_endo,endo))>0){
      cat("The following variables appear in the equations and have been identified by the parser and are not endogenous. Perhaps they are exogenous variables. Or something is wrong in the syntax of the input. Please refer to user manual for further help. \n")
      print(setdiff(vector_endo,endo))
      ok_checks <- FALSE
   }

   if(length(endo) != nrow(eqns)){cat('The model needs to be correctly identified : as many equations as endogenous variables. \n')
      ok_checks <- FALSE
   }

   endot0_check<-apply(eqns,1,function(x)sum(is.na(x)))
   if(max(endot0_check)==ncol(eqns)){
      cat("The following equations have been found to have no endogenous variables in t=0")
      names(endot0_check)<-names_of_equations
      print(names(endot0_check[endot0_check==ncol(eqns)]))
      ok_checks <- FALSE
   }

   if (ok_checks==FALSE){stop("Please check the model input and try again.")}else{cat("The model and variables appear to be correctly identified. \n")}
}

###checking if variable appears in formula
#' which variables in the formula
#'
#' @param variables_to_test character vector of variables to test
#' @param formulas character vector of formulas
#' @param print_type display the variable types that are being tested
#' @keywords internal
#'
#' @return vector of variables that have been found in the formula
is_in_formulas <-function(variables_to_test, formulas, print_type =""){
   formulas <- gsub("^\\w+:","",formulas)
   variables_to_test<- tolower(variables_to_test)
   all_formulas_in_one <-tolower(paste(formulas,collapse = "-" ))

   variables_in_formula<-get_variables_from_string(all_formulas_in_one)
   absent_variables<-setdiff(variables_to_test,variables_in_formula)

   if(length(absent_variables)==0){
      cat(paste("All", print_type,"variables are present in the equations list. \n"))

      }else{
      cat(paste("The following", print_type,"variable(s) do(es) not appear in any equations, they should be dropped.\n"))
      print(absent_variables)}
   res <- setdiff(variables_to_test,absent_variables)
}

#' is the variable name acceptable
#'
#' @param text character string to test
#' @keywords internal
#'
#' @return boolean
acceptable_var_name<-function(text){
   check <- TRUE
   if(sum(grepl("^[a-z]",text))!=length(text)){
      cat('\nNames should start with a lower case letter only.\n')
      print(text[!grepl("^[a-z]",text)])
      check <- FALSE
   }
   if(sum(nchar(text)==0)>0){
      cat('\nNames should have at least one character.\n')
      check <- FALSE
   }

   treatment<- gsub("[a-z_0-9]","",text)
   if(sum(nchar(treatment)>0)>0){
      cat('\nInvalid characters found. Only lower case, alphanumerical characters and underscores are allowed.\n')
      print(text[nchar(treatment)>0])
      check <- FALSE
   }

   check
}

#' check if the lags and delta can be safely parser
#'
#' @param formula_list vector of formulas
#' @importFrom stringr str_count
#' @keywords internal
#' @return boolean
parser_lag_delta_check<-function(formula_list){
   check <- TRUE
   all_formulas_in_one <-tolower(paste(formula_list,collapse = "-" ))

   n_lags<-stringr::str_count(all_formulas_in_one,pattern = "(mylg|lag)\\(")
   correct_lags<- stringr::str_count(all_formulas_in_one,pattern = "(mylg|lag)\\([a-z](\\w+)?,")
   if(correct_lags != n_lags){cat("\nSome lagged variables are not correctly written in the equations. The parser might fail. Only one variable may be lagged at once :
This is ok : lag(my_variable42,...)
This is invalid : lag(my_variable42 + var2,...) ")
      check <- FALSE}

   n_deltas<-stringr::str_count(all_formulas_in_one,pattern = "delta\\(")
   correct_deltas<- stringr::str_count(all_formulas_in_one,pattern = "delta\\([0-9]+,")
   if(correct_deltas != n_deltas){cat("\nSome deltas are not correctly written in the equations. The parser might fail. Only one variable may be lagged at once :
This is ok : delta(42,...)
This is invalid :  delta(var42,...) ")
      check <- FALSE}

   check

}
