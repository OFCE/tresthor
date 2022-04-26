#' Information about variables in the model.
#'
#' @param model thoR.model of interest
#' @param variables character vector of variables to look up.
#'
#' @return Prints information about the variables : type and the equations where they appear.
#' @export
#'
var_info_model<-function(variables,model){
  if(class(model)!="thoR.model"){stop("The model must be a thoR.model object.")}
  if(class(variables)!="character"){stop("The variables must be specified in a string vector.")}

  model_var_list <- names(model@var_map)
  if(prod(variables %in% model_var_list) == 0 ){
    cat("\n The following variables were not found :\n")
    print(setdiff(variables,model_var_list))
    variables <- intersect(variables,model_var_list)
  }

 if(length(variables) == 0){stop("Variables are not in the model.")}

  var_info<-(purrr::set_names(variables) %>% purrr::map(~model@var_map[[.x]]))
  var = variables[1]

  print_info<-function(var){
    info_var<-var_info[[var]]
    equations_of_var<- model@equation_list[which(model@equation_list$name %in% info_var$equations),]
  cat(paste0("\n************************\n \n INFORMATION ON ",var, "\n\n TYPE:\n",paste(info_var$type,collapse = ", ")))
  cat(paste0("\n\n EQUATIONS:\n",paste(info_var$equations,collapse = ",")," \n"))
  cat("\n EQUATIONS FORMULAS:\n")
  cat(paste0(equations_of_var$id ," ",equations_of_var$name," : ",equations_of_var$equation, "\n \n"))
  cat("\n ")


  }
display<-purrr::set_names(variables) %>% purrr::map(~print_info(.x))

}


#' Information about equations in the model
#'
#' @param equations_name_or_id equations names or id in a character vector.
#' @param model thoR.model with the equations
#'
#' @return Prints information about the equations
#' @export
#'
equation_info_model <- function(equations_name_or_id,model){
  if(class(model)!="thoR.model"){stop("The model must be a thoR.model object.")}
  if(class(equations_name_or_id)!="character"){stop("The equations must be specified in a string vector.")}


  equations_list<-model@equation_list
  list_name_id <- unique(c(equations_list$id,equations_list$name))

  if(prod(equations_name_or_id %in% list_name_id) == 0 ){
    cat("\n The following equations were not found :\n")
    print(setdiff(equations_name_or_id,list_name_id))
    equations_name_or_id <- intersect(equations_name_or_id,list_name_id)
  }

  if(length(equations_name_or_id) == 0){stop("Equations specified are not in the model.")}


  equations_list_info<- equations_list[which(equations_list$id %in% equations_name_or_id |equations_list$name %in% equations_name_or_id ),c("id","name","equation","part")]

  equations_list_info$index <- c(1:nrow(equations_list_info))
  equations_id<-equations_list_info$id

  print_info<-function(equation_index){
    cat(paste0("\n************************\n\n EQUATION ID: ",equations_list_info$id[equation_index], " | NAME : ",equations_list_info$name[equation_index], " |  PART : ",equations_list_info$part[equation_index]))
    equation_vars <- get_variables_from_string(equations_list_info$equation[equation_index])
    endogenous<- intersect( equation_vars,model@endo_list)
    exogenous <- intersect( equation_vars,model@exo_list)
    coeff<- intersect( equation_vars,model@coeff_list)
    cat(paste0("\n\n ENDOGENOUS VARIABLES:\n",paste(endogenous,collapse = ", ")))
    cat(paste0("\n\n EXOGENOUS VARIABLES:\n",paste(exogenous,collapse = ", ")))
    cat(paste0("\n\n COEFFICIENT VARIABLES:\n",paste(coeff,collapse = ", ")))
    cat(paste0("\n\n FORMULAS:\n"),equations_list_info$equation[equation_index],"\n")
    # cat(paste0(equations_of_var$id ," ",equations_of_var$name," : ",equations_of_var$equation, "\n \n"))
    # cat("\n \n")
    }
    display<-equations_list_info$index %>% purrr::map(~print_info(.x))


}



