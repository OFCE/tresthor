### check one

### variables report returns a table of variables with nas at given time periods
#' Report which variables have NAs for given periods
#'
#' @param database a dataframe containing all necessary variables
#' @param times periods in index_times to test
#' @param variables variables to check
#' @param index_time name of column in the database which is meant to be use as the time vector
#'
#' @return vector of variable names with NAs over the time period
#' @export
na_report_variables_times <- function(database,times,variables, index_time ="date"){
  #check time
  if(time_model_checks(database,index_time)==FALSE){stop("Problems with the index_time variable specified.")}

   database$TIME <- database[,index_time]
   if (all(times %in% database$TIME) == FALSE){stop( paste0('At least one of ',times, " was not found in ",index_time, "." ))}


  ####
    subdata <- database[which(database$TIME  %in% times),variables]
    res <- names(subdata)[unlist(lapply(subdata, function(x) any(is.na(x))))]

}

### Tests if the model can be solved at a given time (no missing data)

#' Tests if the solver has enough data to run
#' useful to check the need of lag datas
#' @param database a dataframe containing all necessary variables
#' @param model a thoR.model object or thoR.equation
#' @param index_time name of column in the database which is meant to be use as the time vector
#' @param times periods to test.
#'
#' @return boolean, TRUE if successful
#' @export
time_solver_test_run<-function(model , database , index_time = 'date' , times){
  times<-times[order(times)]
#check time
if(time_model_checks(database,index_time)==FALSE){stop("Problems with the index_time variable specified.")}
#check class model
assertthat::assert_that(class(model)=="thoR.equation"|class(model)=="thoR.model")
#check variables
if(data_model_checks(model,database)==FALSE){stop("Missing model variables in the database.")}

database$TIME <- database[,index_time]

if (all(times %in% database$TIME) == FALSE){stop( paste0('At least one of ',times, " was not found in ",index_time, "." ))}
  anchor_t<- which(database$TIME == times[1])




test_run<-function(time_t){
test_t <- which(database$TIME == time_t)

if(class(model)=="thoR.model"){
  if(anchor_t == 1){stop("The first period cannot be the first observation the data base. The solver needs a previous full observation to initialize.")}else{
    ##check that the previous period has no na.
    check <-sum(is.na(database[anchor_t-1,model@endo_list]))
    if(check >0){stop("The solver needs a previous full observation to initialize.")} }
  t_data<-database
  t_data[test_t,model@endo_list]<-t_data[test_t-1,model@endo_list]
  test_vect <- if(model@prologue == TRUE){sum(is.na( model@prologue_equations_f(test_t,t_data)) )}else{0} +
               if(model@heart == TRUE){sum(is.na( model@heart_equations_f(test_t,t_data)) )}else{0} +
               if(model@epilogue == TRUE){sum(is.na( model@epilogue_equations_f(test_t,t_data)) )}else{0}



  }

##when model is equation
if(class(model)=="thoR.equation"){
  if(anchor_t == 1){stop("The first period cannot be the first observation the data base. The solver needs a previous full observation to initialize.")}else{
    ##check that the previous period has no na.
    check <-sum(is.na(database[anchor_t-1,model@endogenous]))
    if(check >0){stop("The solver needs a previous full observation to initialize.")} }
  t_data<-database
  t_data[test_t,model@endogenous]<-t_data[test_t-1,model@endogenous]

  test_vect <- sum(is.na(model@eq_f(test_t,t_data)))

}

#res
if(test_vect > 0){cat(paste0("Insufficient data to solve the model at ",time_t,".\n"))
                  check <- FALSE}else{check <- TRUE}}

# time_t <- c("1981Q2","1981Q3","1981Q4","1982Q1","1982Q2","1982Q3","1982Q4","1983Q1","1983Q2","1983Q3")

checks<-purrr::map_lgl(set_names(times),~test_run(.x))

checks
}
