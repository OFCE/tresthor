

#' Simulate an equation
#'
#' @param thor_equation thoR.equation object
#' @param database data.frame with data for all variables in the equation.
#' @param start_sim  Start period to start the simulations.  Same class as used in index_time. Must be in index_time.
#' @param end_sim Last period to end the simulations.  Same class as used in index_time. Must be in index_time.
#' @param index_time Name of the column in the database which is meant to be use as the time vector.
#' @param residual_var character string, length 1. Variable of the equation that serves as the residual var.
#' @import dplyr
#' @import purrr
#' @return data.frame with the equation variable data. The "observed" variable represent the endogenous variable before simulation. The "simulated" variable is after estimation with the residual variable set 0. The residual is also estimated, as well as its contribution.
#' @export
#'
simulate_equation<-function(thor_equation, database, start_sim, end_sim, index_time="date", residual_var){
  if(class(thor_equation)!= "thoR.equation"){stop("thor_equation must be a thoR.equation object.")}


  if(!residual_var %in% thor_equation@exogenous){stop("the residual variable must be an exogenous variable of the thor_equation.")}

  if(data_model_checks(thor_equation,database)==FALSE){stop("Missing variables in the database. The solver won't work.'")}
  if(time_model_checks(database,index_time)==FALSE){stop("The index_time variable is not suitable. The solver won't work'.")}

  allvars <- c(index_time,thor_equation@endogenous,thor_equation@exogenous,thor_equation@coefflist)
  database <- database[,allvars]
  rownames(database)<- database[,index_time]
  database[,residual_var]<-0

  anchor_t <- which(database[,index_time]==start_sim)
  final_t  <- which(database[,index_time]==end_sim)

  simul_data <- thor_equation_solver(equation = thor_equation,first_period = start_sim,last_period = end_sim,database = database,index_time = index_time)%>% select(index_time,paste0(thor_equation@endogenous,".simul"))

colnames(database[,c(index_time,thor_equation@endogenous)])<- c(index_time,"observed")

residual<-purrr::quietly(create_equation)("residual",formula = thor_equation@formula,coefflist = thor_equation@coefflist,endogenous = residual_var,env=environment())$result
# switch data observ simul
 resid_data<- thor_equation_solver(residual,first_period = start_sim,last_period = end_sim,database = database,index_time = index_time) %>% select(index_time,paste0(residual_var,".simul"))

output_data<- database %>% left_join(simul_data,by = all_of(index_time)) %>% left_join(resid_data,by = all_of(index_time)) %>%
  rename(observed = thor_equation@endogenous ,
         simulated = paste0(thor_equation@endogenous,'.simul') ,
         residual = paste0(residual_var,".simul") ) %>%
  mutate(g_obs = observed/lag(observed) -1,
         g_sim = simulated/lag(simulated) -1,
         dlog_obs = log(observed) - log(lag(observed)),
         dlog_sim = log(simulated) -log(lag(simulated)),
         d_obs = observed - lag(observed),
         d_simul = simulated - lag(simulated),
         residual.contrib = residual - lag(residual)
          )

output_data<- output_data[(min(which(!is.na(output_data$observed)))):final_t,]
}


#' Dynamic contributions for ECM equations
#'
#' @param thor_equation thoR.equation object, must be an error correction model type of equation, otherwise results wont be accurate.
#' @param database data.frame with data for all variables in the equation.
#' @param start_sim  Start period to start the simulations.  Same class as used in index_time. Must be in index_time.
#' @param end_sim Last period to end the simulations.  Same class as used in index_time. Must be in index_time.
#' @param index_time Name of the column in the database which is meant to be use as the time vector.
#' @param residual_var character string, length 1. Variable of the equation that serves as the residual var.
#' @param regroup_these character vector of variables which contributions will be summed. Default: NULL
#' @param name_group character string to name the grouped variables contribution. Default: "group"
#' @param print_verif boolean. TRUE if verification of contributions and residuals should be printed. Default: TRUE
#' @return data.frame containing dynamic contributions as well as simulation results
#' @export
#'
dyn_contribs <- function(thor_equation,database, start_sim, end_sim, index_time="date", residual_var,regroup_these = NULL, name_group = "group",print_verif=TRUE){
  data_sim <- quietly(simulate_equation)(thor_equation , database ,start_sim, end_sim, index_time,residual_var)$result

  maxlag <- thor_equation@maxlag

 if (length(thor_equation@coefflist)>0){ formula_num <- formula_with_coeffs(thor_equation@formula,coefflist = thor_equation@coefflist, database = database,round_digits = 10,quiet = TRUE)
  }else{formula_num<-thor_equation@formula }

  deriv_table <- derivative_time_table(formula_num,maxlag)
  log.endo.true <- grepl("^log\\.",deriv_table[thor_equation@endogenous,1])

  contrib_var_list <- setdiff(rownames(deriv_table), c(residual_var,thor_equation@endogenous))

  endo <-thor_equation@endogenous
  exo <- contrib_var_list

  mainf_endo <- -1*(deriv_table[endo,3:(maxlag+2)])

  plop<-purrr::set_names(contrib_var_list) %>%  purrr::map(~calcul_contrib(.x,data = data_sim,deriv_table,mainf_endo,index_time))
  ##mise en forme des données

  plop_out<-list2DF(unlist(plop,recursive = FALSE))

  contrib<-plop_out[,c(1,grep(pattern = "\\.contrib",names(plop_out)))]
  colnames(contrib)<-gsub("^\\w+\\.(\\w+\\.contrib)$","\\1",colnames(contrib))
  colnames(contrib)<-gsub(paste0("^\\w+\\.(",index_time,")$"),"\\1",colnames(contrib))

  data_contrib <-  purrr::quietly(left_join)(data_sim,contrib)$result

if (print_verif){
  if(log.endo.true){
    verif_sim <- data_contrib$dlog_sim - rowSums(as.data.frame(data_contrib[,paste0(contrib_var_list,".contrib")]))
    verif_obs <- data_contrib$dlog_sim - data_contrib$dlog_obs + data_contrib$residual.contrib}else{
      verif_sim <- data_contrib$d_sim - rowSums(as.data.frame(data_contrib[,paste0(contrib_var_list,".contrib")]))
      verif_obs <- data_contrib$d_sim - data_contrib$d_obs + data_contrib$residual.contrib
    }
    cat("\n Verifications \nContributions: ")
    cat(round(verif_sim[max(length(verif_sim)-50,1):length(verif_sim)],6))
    cat("\n Verifications \nResiduals: ")
    cat(round(verif_sim[max(length(verif_obs)-50,1):length(verif_obs)],6))
    }
 if(!is.null(regroup_these)){
   #check that dummies are in the equation
    if(prod(regroup_these %in% thor_equation@exogenous)==0){stop("some of the variables to regroup specified were not found.")}else{
      group_contrib<-paste0(regroup_these,".contrib")
      data_contrib[,paste0(name_group,".contrib")]<-rowSums(data_contrib[,group_contrib])
      data_contrib[,group_contrib]<-NULL
    }
 }
  rownames(data_contrib)<-data_contrib[,index_time]

data_contrib
}

#' Quartely to yearly dynamic contributions
#'
#' @description Transforms  quarterly dynamic contribution data into yearly data, ready for graphics
#'
#' @param contrib_quarterly_database quarterly dynamic contributions data generated by dyn_contrib()
#' @param index_year year vector that matches the observations in the data
#' @param quarter_start numeric. Which is the first quarter of the database (1,2,3 or 4). Default: 1.
#'
#' @return annual database with yearly dynamic contribution
#' @import dplyr
#' @import purrr
#' @export
#'
yearly_contrib<-function(contrib_quarterly_database,index_year,quarter_start = 1){
  if(length(index_year) != nrow(contrib_quarterly_database)){stop()}
  variables = c(names(contrib_quarterly_database)[grep(pattern = ".contrib$",names(contrib_quarterly_database))],"dlog_obs")

  ann_list <- purrr::map(variables,~quarterly_to_annual(series = .x,database = contrib_quarterly_database,index_year = index_year,quarter_start = quarter_start))

  data_ann <- purrr::quietly(reduce)(ann_list,left_join)$result

  return<-data_ann %>% mutate_all(function(x){as.numeric(x)})
}



