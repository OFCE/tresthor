##checking data


#' Checking if all needed variables is in the data
#'
#' @param thor_model a thoR.model object or thoR.equation
#' @param database the database to check
#' @param quiet Display message if successful. Default : TRUE
#' @return boolean : TRUE if all ok, FALSE if there was a problem
#' @export
#'
data_model_checks <- function(thor_model,database,quiet= TRUE){
  if(class(thor_model)=="thoR.model"){
  vars<-names(thor_model@var_map)
  vars_data<-names(database)

  if(length(setdiff(vars,vars_data))==0){if (quiet == FALSE){cat("All model variables are in the database. \n")}
    check <-TRUE}else{check <-FALSE}

  if(length(setdiff(thor_model@endo_list,vars_data))>0){cat("The following endogenous variables are missing from the database. \n")
                                                        print(setdiff(thor_model@endo_list,vars_data))}
  if(length(setdiff(thor_model@exo_list,vars_data))>0){cat("The following exogenous variables are missing from the database. \n")
                                                        print(setdiff(thor_model@exo_list,vars_data))}
  if(length(setdiff(thor_model@coeff_list,vars_data))>0){cat("The following coefficient variables are missing from the database. \n")
                                                        print(setdiff(thor_model@coeff_list,vars_data))}
  }

  ##when model is equation
  if(class(thor_model)=="thoR.equation"){

    vars<-c(thor_model@coefflist,thor_model@endogenous,thor_model@exogenous)
    vars_data<-names(database)

    if(length(setdiff(vars,vars_data))==0){if (quiet == FALSE){cat("All model variables are in the database. \n")}
      check <-TRUE}else{check <-FALSE}

    if(length(setdiff(thor_model@endogenous,vars_data))>0){cat("The endogenous variable is missing from the database. \n")
      print(thor_model@endogenous)}
    if(length(setdiff(thor_model@exogenous,vars_data))>0){cat("The following exogenous variables are missing from the database. \n")
      print(setdiff(thor_model@exogenous,vars_data))}
    if(length(setdiff(thor_model@coefflist,vars_data))>0){cat("The following coefficient variables are missing from the database. \n")
      print(setdiff(thor_model@coefflist,vars_data))}

    }

  ###res
 check
  }


#' Checking the viability of a variable as a time index
#' Condition : must have no duplicates, no NAs and sorted in ascending order
#' @param database dataframe the database to check
#' @param index_time character string : name of the variable to use as time_index
#' @keywords internal
#'
#' @return logical TRUE if the time_index chosen can be used for the solver
time_model_checks <- function(database,index_time){
  check <- TRUE
  assertthat::is.string(index_time)
  #check that index_time is in database
  if(!index_time %in% names(database)){
  stop(paste0(index_time," is not a variable of the specified data base. \n"))
  }
  #check that index_time is a unique, non NA vector that is sorted in increasing order
  time_vec <- as.vector(database[,index_time])
  if (sum(duplicated(time_vec))>0){cat("The chosen index_time variable has some duplicated values. It is not a viable time index.\n")
    check <-FALSE}
  if (sum(is.na(time_vec))>0){cat("The chosen index_time variable has some NAs. It is not a viable time index.\n")
    check <-FALSE}
  if(identical(time_vec, time_vec[ordered(time_vec)])== FALSE){cat("The chosen index_time variable is not sorted in ascending order. It is not a viable time index.\n")
    check <-FALSE}

  check

}
#
