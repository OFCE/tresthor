##Single equation function

#' determining the maximum lag an equation with econometric format (newdiff instead of delta)
#' the formula must have been formatted to be evaluated by thor first
#'
#' @param new_formula_vector vector containing equations with the format that can be evaluated
#' @param default_error default value in case of error in algorithm. Default : NA
#'
#' @return maxlag value or NA if couldn't determine it.
#' @keywords internal
find_maxlag<-function(new_formula_vector,default_error = NA){
  if (sum(grepl(pattern = "lag\\(\\w+,-?\\(?[a-z]+",new_formula_vector))>0){
    cat("\n Some lags are variable dependant. maxlag cannot be computed, and is set to default.
        The equation cannot be used for econometric estimations in its current form.")
    maxlag <-default_error
  }else{
    functions_supported<-paste0(thor_functions_supported,collapse = "\\(|")
    data_test<-as.data.frame(cbind(lag=c(1:100),plopvar=1,maxlag=0))

    count_maxlag<-function(string){
    maxlag_function<-function(string){
      testlag<-string
      testlag<-gsub("lag","%",testlag)
      testlag<-gsub("newdiff","@",testlag)
      testlag<-gsub(paste0(functions_supported,"\\(|I\\("),"(",testlag)
      testlag<-gsub("[a-z](\\w+)?","data_test$plopvar",testlag)
      testlag<-gsub("%","lag",testlag)
      testlag<-gsub("@","newdiff",testlag)
      data_test[,"res"]=eval(parse(text=testlag))
      return<-sum(is.na(data_test[,"res"]))
    }
    safe_maxlag<-safely(maxlag_function,otherwise = NA)
    safe_maxlag(string)$result
    }
    maxlag_vec<-new_formula_vector %>% purrr::map_dbl(~count_maxlag(.x))
    if(sum(is.na(maxlag_vec))>0){
      cat("\n Some lags cannot be computed. maxlag is set to default.
        The equation cannot be used for econometric estimations in its current form.")
    maxlag <-default_error
    }else{maxlag<-max(maxlag_vec)}
  }
  maxlag
}


#' Prepares the RHS of the equation and reorder coefficient list for econometric estimation
#'
#' @param RHSformula right hand side of formula
#' @param LHSformula left hand side of formula
#' @param coeffs_variables list of coefficients
#' @param eq_vars equation variables
#'
#' @return list with modified RHS and coefficients in the right order
#' @keywords internal
reorder_coefficients <- function( LHSformula,RHSformula, coeffs_variables,eq_vars){
# LHSformula=LHS
# RHSformula=RHS
# coeffs_variables=coefflist
# eq_vars=eq_vars
  check_safe<-TRUE

  liste<-coeffs_variables
  liste<-liste[order(nchar(liste), liste)]
  RHStest<-RHSformula
  RHSformula_ec <- RHSformula
  for( i in seq_along(liste)){
    RHStest<-gsub( paste("-?","(\\s+)?",liste[i],"(\\s+)?","(\\*|\\+|\\-)","|(-|\\+)",liste[i],"$",sep = "") , paste("XX",i,"XX",sep="") , RHStest)
    RHSformula_ec<-gsub(paste(liste[i],"\\*\\(" ,sep = "" ),"I\\(",RHSformula_ec)
    RHSformula_ec<-gsub(paste(liste[i],"\\+" ,sep = "" ),"",RHSformula_ec)
    RHSformula_ec<-gsub(paste(liste[i],"\\*" ,sep = "" ),"",RHSformula_ec)
    RHSformula_ec<-gsub(paste("(-|\\+)",liste[i],"$",sep = "" ),"\\11",RHSformula_ec)
  }

RHS_seq<-str_extract_all(RHStest,"XX[0-9]+XX")
RHS_seqn<-as.numeric(vapply(RHS_seq[[1]],function(X)gsub("XX","",X),character(1)  ))
ordered_coeff<-liste[RHS_seqn]
if(length(ordered_coeff)!=length(coeffs_variables)){cat("\n Coefficients parser error, the formula may be difficult to use in econometric estimations.")
  check_safe<-FALSE}

testing_ols<-function(LHS=LHSformula,RHS=RHSformula_ec,eq_variables=eq_vars){
data_test<- as.data.frame(matrix(sample(1:100,size = length(eq_variables)*30,replace= TRUE),ncol=length(eq_variables)))
colnames(data_test)<-eq_variables
ols_test<- lm(paste(LHS,RHS,sep = "~"),data = data_test)
length(ols_test$coefficients)
}

safe_test<-safely(testing_ols,otherwise = 0)

n_coeff<-safe_test()$result

if(n_coeff %in% c(length(coeffs_variables),length(coeffs_variables)+1)){
  cat("\n RHS of equation appears to be ok for econometric estimations with the current coefficient list.")
}else{cat("\n RHS of equation does not appear to be ok for econometric estimations with the current coefficient list.")
  check_safe<-FALSE}
return<-list(check_safe,ordered_coeff,RHSformula_ec)
  }


#' Shows the formula with the coefficients' value in numeric instead of the coefficient variable names.
#'
#' @param formula a formula as a character string or a thor equation
#' @param coefflist vector with variables that could be coefficients. If the formula input is a thoR.equation model, this argument is ignored, and the coefficients are taken from the object parameters
#' @param round_digits Numeric. numbers of digits to round the value of the coefficient. Default : 4
#' @param quiet boolean. FALSE doesn't display the result. Default : TRUE
#' @param database database that contains the values of the formula coefficient. Name of coefficient must be as column names and values the same across all observations of the base. The first row will be used. Named vector should also work.
#'
#' @return character string
#' @export
#'
formula_with_coeffs  <-function(formula,coefflist=NULL,database,round_digits = 4,quiet=FALSE){
  options(scipen = 99)
    if(class(formula) != "character"){if(class(formula)!="thoR.equation"){stop("formula needs to be a character string of length 1 or a thoR.equation.")} }
    if(length(formula)!= 1){stop("formula needs to be a character string of length 1 ")}
    if(class(formula)== "character") { if(class(coefflist) != "character"){stop("If formula is not a thoR.equation object, the coefficients need to specificed as a character vector in coefflist.")}}

    if(class(formula) == "character"){
      detected_vars <- get_variables_from_string(formula)
      coefflist<-intersect(coefflist,detected_vars)
      string<-formula
   }
    if(class(formula) == "thoR.equation"){
      string <- formula@formula
      coefflist<- formula@coefflist
    }

    if(length(coefflist)==0){stop("no coefficients to replace.")}

    if(is.null(names(database))){stop("database needs to have names.")}else{
    plop_data_coeff<-as.data.frame(database[,coefflist])
    colnames(plop_data_coeff)<-coefflist
    naam<-names(plop_data_coeff)
    }

     if(prod(coefflist %in% naam)==0){stop("Missing coefficients in database.")}

    if(class(database)=="data.frame" | class(database)=="matrix"){
      if(sum(is.na(database[1,coefflist]))){stop("NAs in the first rows of the dataframe for coefficients. Coefficient should have the same values across all observations in the database.")}
      value_list<-as.numeric(as.vector(as.matrix(database[1,coefflist])))
    }
    if(class(database)=="numeric"){
      value_list<-as.vector(database)
    }

    names(value_list)<- naam

    coefflist<- coefflist[order(nchar(coefflist))]
    for(coeff in coefflist){
      value <- as.character(round(value_list[coeff],round_digits))
      string<-gsub(pattern =paste0("(?<![a-z0-9_])",coeff,"(?![a-z0-9_])"),replacement = value,string,perl = TRUE)
    }

 if(quiet==FALSE){cat(string)}
 return<-string
  }


#' @title Pull an equation from a model
#' @description This function returns a thoR.equation object created from an equation in a thoR.model.
#'
#' @param equation_name_or_id character string of length 1 with the equation name or id, must match the name or id used in the model.
#' @param model thoR.model object where to pull the equation from
#' @param endogenous character string of length 1 that is a variable of the equation and will be used as the equation's endogenous variable.
#' @param new_name character string to designate the new name of the model. if NULL, the name will be the same as the equation_name_or_id. Default : NULL.
#' @param default_maxlag_error if maxlag fails it will default to this value. Default : 5.
#' @param environment environment where the equation created must be loaded. Default : globalenv()
#' @return a thoR.equation from a thoR.model
#' @export
#'
equation_from_model<-function(equation_name_or_id,model,endogenous,new_name=NULL,default_maxlag_error = 5,environment = globalenv()){
  if(class(model)!="thoR.model"){stop("The model must be a thoR.model object.")}

  if(class(equation_name_or_id)!="character"){stop("Expecting a character string for equation_name_or_id.")}
  if(length(equation_name_or_id) !=1){stop("Only one equation at a time can be created.")}
  if(!is.null(new_name)){if(class(new_name)!="character"){stop("if a new_name is specified, it must be a character string")}}else{new_name <-equation_name_or_id}

  if(length(endogenous) !=1){stop("One exogenous only must be specified for the equation.")}
  if(class(endogenous)!= "character"){stop("Expecting a character string for endogenous.")}



  if(equation_name_or_id %in% model@equation_list$name){
    equation_name_or_id <- model@equation_list$id[which(model@equation_list$name==equation_name_or_id)]}
  if(!equation_name_or_id %in% model@equation_list$id){stop("equation_name_or_id was not found in the model.")}

  formula_m <- model@equation_list$equation[which(model@equation_list$id==equation_name_or_id)]

  formula_coeffs<- intersect(get_variables_from_string(formula_m),model@coeff_list)
  create_equation(equation_name = new_name,formula = formula_m,coefflist = formula_coeffs,endogenous = endogenous,env = environment,default_maxlag_error = default_maxlag_error )

}

