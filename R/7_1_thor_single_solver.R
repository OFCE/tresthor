# thoR.single.solver

#' thor solver
#'
#' @param equation A thoR.equation
#' @param first_period First period to solve. Same class as used in index_time. Must be in index_time.
#' @param last_period Last period to solve. Same class as used in index_time. Must be in index_time.
#' @param database data.frame database to be used. default is t_data
#' @param index_time name of column in the database which is meant to be use as the time vector
#' @param convergence_criteria initial convergence criteria of the model. Default : 10e^-16

#'
#' @return a data.frame with only the variables in the equation. The endogenous solved appears as endogenous.simul
#'
#' @export
#'
thor_equation_solver<-function(equation,
                      first_period,
                      last_period,
                      database = t_data,
                      index_time="date",
                      convergence_criteria=0.0000000000001){

  #############
  ##1.checks
  ##############
  model = equation
  assertthat::assert_that(convergence_criteria < 0.001)
  #check that data contains all variables used in model
  if(data_model_checks(model,database)==FALSE){stop("Missing variables in the database. The solver won't work.'")}
  if(time_model_checks(database,index_time)==FALSE){stop("The index_time variable is not suitable. The solver won't work'.")}
  #check that first and last period are part of index_time and that first <= last
  time_vec <- database[,index_time]
  if(!first_period %in% time_vec){stop("The specified first period is not found in the time index variable")}
  if(!last_period %in% time_vec){stop("The specified  last period is not found in the time index variable")}

  if(is.numeric(database[,index_time])){
    database[,index_time]<-as.character(database[,index_time])
    first_period <- as.character(first_period)
    last_period <- as.character(last_period)
    numeric_index_time = TRUE}else{numeric_index_time = FALSE}

  model_variables <-c(equation@coefflist,equation@endogenous,equation@exogenous)

  model_variables <- model_variables[order(model_variables)]
  t_data<-as.data.frame(database[,model_variables])
  rownames(t_data)<-time_vec
  #check for sufficient time_frame of data and that first and last period is ok
  # Can the solver run ant the first period?
  check_first_period <- time_solver_test_run(model,database,index_time,times = first_period)

  anchor_t <- which(database[,index_time]==first_period)
  final_t <- which(database[,index_time]==last_period)

  if(check_first_period == FALSE){
    #1 determine how early you can go

    max_t<- which(database[,index_time]==max(database[,index_time]))
    test_twelve<- c(min(max_t,anchor_t+12))
    diagnosis<-time_solver_test_run(model,database,index_time,times = database[c((anchor_t):test_twelve),index_time])
    if(sum(diagnosis)>0){cat(paste0("The solver cannot run from ",first_period,". It appears to be possible starting at ",min(names(diagnosis[TRUE])),". Try running it from there.") )}else{
      cat(paste0("The solver cannot run from ",first_period,", nor from the next 12 periods. Try running it from at least there, or check your data."))
    }
    stop("Solver cannot be ran from the first period specified.")
  }

  #sufficient data to run after?
  check_subdata<- TRUE

  if (!purrr::is_empty(model@exogenous) | !purrr::is_empty(model@coefflist) ){
  subdata <- database[c((anchor_t):final_t),c(index_time,model@exogenous,model@coefflist)]
  coeff_check <- (na_report_variables_times(subdata,subdata[,index_time],model@coefflist,index_time))
  if (length(coeff_check) == 0){cat(" ")}else{
    cat("The following coefficients are missing for at least one observation in the specified timeframe. \n")
    print(coeff_check)
    check_subdata<- FALSE
  }
  exo_check <- (na_report_variables_times(subdata,subdata[,index_time],model@exogenous,index_time))
  if (length(exo_check) == 0){cat(" ")}else{
    cat("The following exogenous variables are missing for at least one observation in the specified timeframe. \n")
    print(exo_check)
    check_subdata<- FALSE
  }
}
  if(check_subdata==FALSE){stop("Missing some variables to run the solver. Please check your data.")}

  #############
  ##2.solver
  #############
  ## extractions des infos nécessaires du modèle
  endogenous <-model@endogenous
  ## lists

  #################################################
  parameterconvergence<- convergence_criteria
  start=anchor_t
  end_period=final_t
  convergence=convergence_criteria
  #resetting endogenous variables that will be computed

  t_data[c(start:end_period),endogenous]<-NA

  ### Solver

  ###initialisation
    for (i in seq(start,end_period,1)) { #crochet solver
      timeref<-i
      t<-i
      timestart<-timeref-1
      if(timestart>0 & !is.na(t_data[timestart,endogenous])){t_data[timeref,endogenous]<-t_data[timestart,endogenous]}else
      {t_data[timeref,endogenous]<-1}


      ##### Solver starts here
      epsilon <- convergence

          x_n<-as.vector(t_data[timeref,endogenous]) %>% unlist()
          x_n1<-as.vector(1)

          if(x_n == x_n1){x_n1<- x_n + 1 }

          n=0 ##n number of iterations for convergence
          while ( (abs(x_n1-x_n)>epsilon) ){

            if (n>0) {
              x_n<-x_n1
              t_data[timeref,endogenous]<-x_n1 }

            Jacobian_n <- model@jac_f(t,t_data,model@jacobian)
            f_x_n      <- model@eq_f(t,t_data)

            inv_J_n<-as.matrix(solve(Jacobian_n, tol=10^-20, sparse=TRUE))

            x_n1<-x_n - (inv_J_n %*% f_x_n)
            names(x_n1) <- names(x_n)

            n<-n+1
            if(n>20){
              epsilon <- epsilon*10
              if(epsilon > parameterconvergence){parameterconvergence <- epsilon}}
          } ##convergence curl
          t_data[timeref,endogenous]<-x_n1
        }#time curls

    #crochet solver  R

  #############
  ##3.Results
  #############
  if(parameterconvergence>convergence_criteria){cat("\n The convergence criteria had to be increased to ",parameterconvergence," to solve the equation.\n")}

  rownames(t_data) <-rownames(database)
  database[,paste0(endogenous,".simul")]<-t_data[,endogenous]

  return <- database[,c(index_time,model_variables,paste0(endogenous,".simul"))]

}

