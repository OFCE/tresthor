###Build variables_info
#to improve over time

#' Building the variables map of the model
#'
#' @param equations_list_df equations list dataframe
#' @param endos full endos vector
#' @param exos exo vector
#' @param coeffs coefficient vector
#' @param p_endos prologue endo vector
#' @param e_endos epilogue endo vector
#' @param h_endos heart endo vector
#' @import purrr
#' @keywords internal
#' @return list that gives info on each variables
build_dico <-function(equations_list_df, endos, exos , coeffs,p_endos,e_endos,h_endos){

  var_vect <- c(endos,exos,coeffs)
  var_vect <-var_vect[order(var_vect)]

  # table_contemp <- table_contemporaneous_endos(formula_list = equations_list_df$formula,endogenous = var_vect,exogenous = c(),coefflist = c(),equations_index = equations_list_df$id )
  ## equation table
  equationvector<-equations_list_df$equation
  equationvector<-purrr::set_names(equationvector, equations_list_df$name)
  equation_map<- equationvector %>% map(~get_variables_from_string(.x) )
  var_vect <- set_names(var_vect)

  is_in_equation<- function(var){
  in_eq_true<-set_names(names(equationvector)) %>% map_lgl(~(var %in% equation_map[[.x]]))
  result<-list(names(in_eq_true[which(in_eq_true == TRUE)]))
  names(result)<-"equations"
  result}

  dico_equations <-var_vect %>% purrr::map(~is_in_equation(.x))



  type_of_variable <- function(variable){
    if(variable %in% p_endos){res<-c("endogenous","prologue")}
    if(variable %in% h_endos){res<-c("endogenous","heart")}
    if(variable %in% e_endos){res<-c("endogenous","epilogue")}

    if(variable %in% exos){res <- "exogenous"}
    if(variable %in% coeffs){res <- "coefficient"}
    result <- list(res)
    names(result)<-"type"
    result
  }

  dico_type <- set_names(var_vect) %>% purrr::map(~type_of_variable(.x))

 full_dico <- Map(c,dico_equations,dico_type)
  }
