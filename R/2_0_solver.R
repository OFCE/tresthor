
#' thor solver
#'
#' @param model A thoR.model
#' @param first_period First period to solve. Same class as used in index_time. Must be in index_time.
#' @param last_period Last period to solve. Same class as used in index_time. Must be in index_time.
#' @param main_variable Main variable of interest for quick preview of the info.Must be a variable of the dataframe.
#' @param database data.frame database to be used. default is t_data
#' @param index_time name of column in the database which is meant to be use as the time vector
#' @param convergence_criteria initial convergence criteria of the model. Default : 10e^-09
#' @param dynamic_convergence numeric/integer. number of iterations before the solver increases the convergence criteria by a factor of 10. Default : 10
#' @param rcpp boolean. If the rcpp solver is to be used. Default : FALSE
#' @param main_variable_fx function to aply to the main variable displayed. Default : function(x)x
#' @param adv_diag : boolean : gives advanced informations aout the computations. Default : FALSE
#'
#' @return a dataframe that mirrors the database input, but with endogenous variables solved.
#' @import assertthat
#' @importFrom purrr is_empty
#' @export
#'
#'

thor_solver<-function(model,
                      first_period,
                      last_period,
                      main_variable = NULL,main_variable_fx=function(x)x,
                      database = t_data,
                      index_time="date",
                      convergence_criteria=0.000000001,
                      dynamic_convergence=10,
                      rcpp=FALSE,
                      adv_diag = FALSE){

  #############
  ##1.checks
  ##############
 model = model
  assertthat::assert_that(convergence_criteria < 0.01)
  #check that data contains all variables used in model
  if(data_model_checks(model,database)==FALSE){stop("Missing variables in the database. The solver won't work.'")}
  if(!is.null(main_variable)){if(!main_variable %in% names(database)){cat(paste0(main_variable, " doesn't seem to exist in the database. the argument will be ignored. \n"))
    main_variable <-NULL}}

  if(tresthor:::time_model_checks(database,index_time)==FALSE){stop("The index_time variable is not suitable. The solver won't work'.")}
  #check that first and last period are part of index_time and that first <= last
  time_vec <- database[,index_time]
  if(!first_period %in% time_vec){stop("The specified first period is not found in the time index variable")}
  if(!last_period %in% time_vec){stop("The specified  last period is not found in the time index variable")}


  if(is.numeric(database[,index_time])){
    database[,index_time]<-as.character(database[,index_time])
    first_period <- as.character(first_period)
    last_period <- as.character(last_period)
    numeric_index_time = TRUE}else{numeric_index_time = FALSE}
  # t_data<-

  model_variables <- names(model@var_map)
  model_variables <- model_variables[order(model_variables)]
  t_data<-as.data.frame(database[,model_variables])
  rownames(t_data)<-time_vec

  anchor_t <- which(database[,index_time]==first_period)
  final_t <- which(database[,index_time]==last_period)

  #check for sufficient time_frame of data and that first and last period is ok

  # Check that the period before first_period is full so as to initialize the model
  if(anchor_t == 1){stop("The first period cannot be the first observation the data base. The solver needs a previous full observation to initialize.")}else{
    ##check that the previous period has no na.
    check <-sum(is.na(t_data[anchor_t-1,model_variables]))
    if(check >0){stop(" The solver needs a previous full observation to initialize.")}
    }
  # Can the solver run ant the first period?
  check_first_period <- time_solver_test_run(model,database,index_time,times = first_period)



  if(check_first_period == FALSE){
    #1 determine how early you can go

    max_t<- which(database[,index_time]==max(database[,index_time]))
    test_twelve<- c(min(max_t,anchor_t+4))
    diagnosis<-time_solver_test_run(model,database,index_time,times = database[c((anchor_t):test_twelve),index_time])
    if(sum(diagnosis)>0){cat(paste0("The solver cannot run from ",first_period,". It appears to be possible starting at ",min(names(diagnosis[TRUE])),". Try running it from there.") )}else{
      cat(paste0("The solver cannot run from ",first_period,", nor from the next 12 periods. Try running it from at least there, or check your data."))
    }
    stop("Solver cannot be ran from the first period specified.")
    }

    #sufficient data to run after?
  check_subdata<- TRUE
  if (!purrr::is_empty(model@exo_list) | !purrr::is_empty(model@coeff_list) ){
  subdata <- database[c((anchor_t):final_t),c(index_time,model@exo_list,model@coeff_list)]
  coeff_check <- (na_report_variables_times(subdata,subdata[,index_time],model@coeff_list,index_time))
  if (length(coeff_check) == 0){cat(" \n")}else{
    cat("The following coefficients are missing for at least one observation in the specified timeframe. \n")
    print(coeff_check)
  check_subdata<- FALSE
  }
  exo_check <- (na_report_variables_times(subdata,subdata[,index_time],model@exo_list,index_time))
  if (length(exo_check) == 0){cat(" \n")}else{
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
  all_endo <-model@endo_list
  ## lists
  pro_list<-list( bool = model@prologue,
                  s_endo = model@prologue_endo,
                  jacobian  = model@prologue_jacobian,
                  jac_f= model@prologue_jacobian_f,
                  eq_f = model@prologue_equations_f,
                  s_name = "prologue",
                  s_criteria = convergence_criteria)

  hea_list<-list( bool = model@heart,
                  s_endo = model@heart_endo,
                  jacobian  = model@heart_jacobian,
                  jac_f= model@heart_jacobian_f,
                  eq_f = model@heart_equations_f,
                  s_name = "heart",
                  s_criteria = convergence_criteria)

  epi_list<-list( bool = model@epilogue,
                  s_endo = model@epilogue_endo,
                  jacobian  = model@epilogue_jacobian,
                  jac_f= model@epilogue_jacobian_f,
                  eq_f = model@epilogue_equations_f,
                  s_name = "epilogue",
                  s_criteria = convergence_criteria)
  stages_list<-list(pro_list,hea_list,epi_list)
  #################################################

  start=anchor_t
  end_period=final_t
  convergence=convergence_criteria
  #resetting endogenous variables that will be computed

  t_data[c(start:end_period),all_endo]<-NA

  ##for information
  if (!is.null(main_variable) == TRUE){
    cat(paste(main_variable,"before forecast :\n"))
    plop<-(database[,c(index_time,main_variable)])
    plop[,paste0("f_",main_variable)]<- main_variable_fx(plop[,main_variable])
    print(plop[c(max(1,anchor_t-8):anchor_t),c(index_time,paste0("f_",main_variable))])
  }

  #### Solver
 if (rcpp == FALSE){
  ###initialisation
  for (i in seq(start,end_period,1)) { #crochet solver
    timeref<-i
    timestart<-timeref-1
    t_data[timeref,all_endo]<-t_data[timestart,all_endo]
    t=i

    ##### Solver starts here
    for (i_stage in 1:3){
      list2env(stages_list[[i_stage]],envir =environment())

      if (stages_list[[i_stage]]$bool==TRUE){ ## solver for each stage sequentially

        epsilon <- c(rep(convergence,length(s_endo)))

        x_n<-as.vector(t_data[timeref,s_endo])
        x_n1<-as.vector(c(rep(1,length(s_endo))))

        if(sum(abs(x_n1-x_n)>c(rep(epsilon,length(x_n))))==0){x_n1<-as.vector(c(rep(2,length(s_endo))))}

        n=0 ##n number of iterations for convergence
        while ( sum(abs(x_n1-x_n)>c(rep(epsilon,length(x_n))))>0 ){

          if (n>0) {
            x_n<-x_n1

          if(length(s_endo)==1){t_data[timeref,s_endo]<-x_n1
                          }else{t_data[timeref,s_endo]<-x_n1[s_endo]}

          if(n > dynamic_convergence){
            epsilon<-epsilon*10
            n=1}
          }

          Jacobian_n <- jac_f(t,t_data,jacobian)
          f_x_n      <- eq_f(t,t_data)

          sparseJA <- as.matrix(Jacobian_n, "sparseMatrix")
          inv_J_n<-as.matrix(Matrix::solve(sparseJA, tol=10^-20, sparse=TRUE))
          # inv_J_n<-as.matrix(Matrix::solve(Jacobian_n, tol=10^-20))
          x_n1<-x_n - (inv_J_n %*% f_x_n)

          n<-n+1
        } ##convergence curl


        if(length(s_endo)==1){t_data[timeref,s_endo]<-x_n1
                        }else{t_data[timeref,s_endo]<-x_n1[s_endo]}

        if (epsilon[1] > stages_list[[i_stage]]$s_criteria){stages_list[[i_stage]]$s_criteria <- epsilon[1]}
      }
      }#stage curls


    cat(paste("\n ",database[timeref,index_time],"..."))

    }#time curls
   if(adv_diag == TRUE){
     cat("\n
        \n")
     if (model@prologue==TRUE){cat(paste0("\n Prologue max convergence criteria: ",stages_list[[1]]$s_criteria))}
     if (model@heart==TRUE){cat(paste0("\n Heart max convergence criteria: ",stages_list[[2]]$s_criteria))}
     if (model@epilogue==TRUE){cat(paste0("\n Epilogue max convergence criteria: ",stages_list[[3]]$s_criteria))}
   }#crochet solver  R
 }else{
  ##variables are stored in alphabetical order
   prologue_endo<-which(colnames(t_data) %in% model@prologue_endo) -1
   heart_endo<-which(colnames(t_data) %in% model@heart_endo)-1
   epilogue_endo<-which(colnames(t_data) %in% model@epilogue_endo)-1

  endo <-which(colnames(t_data) %in% model@endo_list)-1
  has_prologue <- model@prologue
  has_heart <- model@heart
  has_epilogue <- model@epilogue

  first_date<-anchor_t-1
  last_date<-final_t-1

  time_line <- as.character(database[,index_time])

  columns <- colnames(t_data)
  rows_t <- rownames(t_data)
  convergenceCriteria <- convergence_criteria
  t_data<-data.matrix(t_data)
  # sourceCpp("model_pomme_newton.cpp")
  Rcpp::sourceCpp(model@rcpp_source)
  t_data <- Rcpp_solver(t_data,
                        first_date, last_date,
                        convergenceCriteria,
                        endo, prologue_endo, heart_endo, epilogue_endo,
                        has_prologue, has_heart, has_epilogue, time_line)
  t_data <-as.data.frame(t_data)
  colnames(t_data)<-columns
  rownames(t_data)<-rows_t

}#crochet solver  Rcpp
  #############
  ##3.Results
  #############

 rownames(t_data) <-rownames(database)
 database[,all_endo]<-t_data[,all_endo]

 if (!is.null(main_variable) == TRUE){
   cat(paste(main_variable,"after forecast :\n"))
   plop<-(database[,c(index_time,main_variable)])
   plop[,paste0("f_",main_variable)]<- main_variable_fx(plop[,main_variable])
   print(plop[c(max(1,anchor_t-1):final_t),c(index_time,paste0("f_",main_variable))])
 }
if (numeric_index_time == TRUE){database[,index_time]<-as.numeric(database[,index_time])}
 return <- database

}


