

# Fonction qui retourne la formule txt de l'equation. eq_name = nom de l'equation, thor_model=objet modele
get_formula_txt <- function(eq_name,thor_model){
  equation_list <- thor_model@equation_list
  rownames(equation_list) <- equation_list$name
  formula_txt <- paste(equation_list[eq_name,"LHS"],equation_list[eq_name,"RHS"],sep="=")
  return(formula_txt)
}

# Fonction qui retourne la liste des variables d'une formule txt. thor_model=objet modele
get_variables_list <- function(formula_txt,thor_model){
  variables_list <- get_variables_from_string(formula_txt) %>%
    setdiff(.,thor_model@coeff_list)
  return(variables_list)
}

# Fonction qui retourne la liste des coefficients d'une formule txt. thor_model=objet modele DANS L'ORDRE d'apparition
get_coeff_list <- function(formula_txt,variables_list){
  coeff_list <- get_variables_from_string(formula_txt) %>%
    setdiff(.,variables_list)
  return(coeff_list)
}

# Fonction qui retourne un dataframe indiquant les variables qui doivent etre en log

#' @title Variables with or without logs
#' @param formula_txt plop
#' @param variables_list plop
#' @keywords internal
#' @importFrom gsubfn gsubfn
#' @return dataframe with variables that need to be transformed  by log
get_var_log_all <- function(formula_txt,variables_list){
  f<-gsub("\\\n","",formula_txt)
  f<-gsub("\\s+","",f)
  f<-gsub("(mylg|lag)\\((\\w+),-?([0-9]+)\\)","lagxxplopxx\\2xxplopxx\\3",f)
  f<-gsub("log\\((\\w+(xxplopxx\\w+)*)\\)","logxxplopxx\\1",f)
  f<-gsub("delta\\(([0-9]+),(\\w+(xxplopxx\\w+)*)\\)","(\\2-lagxxplopxx\\2xxplopxx\\1)",f)
  f<-gsub("lagxxplopxx(logxxplopxx)?lagxxplopxx(\\w+)xxplopxx([0-9]+)xxplopxx([0-9]+)","lagxxplopxx\\1\\2xxplopxx\\3+\\4",f)
  f<-gsub("logxxplopxxlagxxplopxx","lagxxplopxxlogxxplopxx",f)
  f<-gsub("lagxxplopxxlagxxplopxx","lagxxplopxx",f)
  f<-gsubfn::gsubfn("(xxplopxx\\w+)xxplopxx([0-9]+\\+[0-9]+)", function(n,m) paste0(n, "xxplopxx", eval(parse(text = m))), f)
  f<-gsub("=","-(",f)
  f<-gsub("$",")",f)

  var_log_all<-vector(length = length(variables_list))
  names(var_log_all)<-variables_list
  for (vari in variables_list){
    plop<-gsub(paste0(".*logxxplopxx",vari,'.*') ,paste0("log\\.",vari),f)
    if (grepl(paste0(".*logxxplopxx",vari,'.*'),f)==FALSE){plop<-vari}
    var_log_all[vari]<-plop
    rm(plop)
  }
  rm(vari)

  return(var_log_all)
}

# Fonction qui retourne une sub_database avec les variables en log si necessaire
get_sub_database <- function(database, var_log_all, variables_list,index_time="date"){
  var_log<-var_log_all[grepl("^log\\.\\w+",var_log_all)==TRUE]
  sub_database <- database[,c(index_time,variables_list)]
  for (vartolog in names(var_log) ){
    nom_var = var_log[vartolog]
    sub_database[,nom_var]<-log(sub_database[,vartolog])
  }
  rm(vartolog)
  return(sub_database)
}

# Fonction qui retourne les variables du LT + envoie dans la console les resultats d'estimation

#' get_var_lt
#'
#' @param formula_txt character string of formula
#' @param coeff_lt character vector of coefficients
#'
#' @importFrom gsubfn gsubfn
#' @keywords internal
get_var_lt <- function(formula_txt,coeff_lt){

  f<-gsub("\\\n","",formula_txt)
  f<-gsub("\\s+","",f)
  f<-gsub("(mylg|lag)\\((\\w+),-?([0-9]+)\\)","lagxxplopxx\\2xxplopxx\\3",f)
  f<-gsub("log\\((\\w+(xxplopxx\\w+)*)\\)","logxxplopxx\\1",f)
  f<-gsub("delta\\(([0-9]+),(\\w+(xxplopxx\\w+)*)\\)","(\\2-lagxxplopxx\\2xxplopxx\\1)",f)
  f<-gsub("lagxxplopxxlogxxplopxxlagxxplopxx(\\w+)xxplopxx([0-9]+)xxplopxx([0-9]+)","lagxxplopxxlogxxplopxx\\1xxplopxx\\2+\\3",f)
  f<-gsub("logxxplopxxlagxxplopxx","lagxxplopxxlogxxplopxx",f)
  f<-gsub("lagxxplopxxlagxxplopxx","lagxxplopxx",f)
  f<-gsubfn::gsubfn("(xxplopxx\\w+)xxplopxx([0-9]+\\+[0-9]+)", function(n,m) paste0(n, "xxplopxx", eval(parse(text = m))), f)
  f<-gsub("=","-(",f)
  f<-gsub("$",")",f)

  var_lt<-vector()
  for (clt in coeff_lt){
    # clt_isolated<-paste0("(?<![a-z0-9_])",clt,"(?![a-z0-9_])")
    # plop<-gsub(paste0(clt_isolated,"(\\w+)?.*"),"\\1",f,perl = TRUE)
    plop<-gsub(paste0(".*",clt,"(\\*(\\()?(\\w+))?.*"),"\\3",f)
    var_lt<-c(var_lt,plop)
    rm(plop)
  }

  rm(clt)
  names(var_lt)<-coeff_lt
  var_lt<-gsub("xxplopxx",".",var_lt)
  var_lt<-gsub("^lag\\.","",var_lt)
  var_lt<-gsub("\\.[0-9]+$","",var_lt)
  return(var_lt)
}

# Fonction qui retourne formula_txt avec les coefficients de LT estimes
#' Estimating long term part of ECM
#'
#' @param var_lt plop
#' @param var_log_all plop
#' @param endogenous_name plop
#' @param coeff_lt plop
#' @param sub_database plop
#' @param estim_start plop
#' @param estim_end plop
#' @param index_time plop
#' @import cointReg
#' @keywords internal
#'
#' @return data.frame with long-term coefficients values
get_coeff_lt_value <- function(var_lt,var_log_all,endogenous_name,coeff_lt,sub_database,estim_start,estim_end,index_time="date"){
  nclt<-length(coeff_lt)
  # nom des variables
  yt = as.character(var_log_all[which(names(var_log_all)==endogenous_name)])
  cst_lt <- as.character(names(var_lt[which(var_lt=="")]))
  xt = var_lt[which(var_lt!="")]

  # plage
  index_start = max(which(sub_database[,index_time]==estim_start) - 4 , 1)
  index_end = min(which(sub_database[,index_time]==estim_end) + 2 , nrow(sub_database) )

  # estimation
  LT<-cointReg::cointRegD(x=sub_database[c(index_start:index_end),xt],y=sub_database[c(index_start:index_end),yt],deter=c(rep(1,(index_end -index_start +1))), n.lead = 2 , n.lag = 2 ,kernel = "qs" , demeaning = TRUE , bandwidth = "and")
  LTcoeff<-LT$theta.all[1:nclt]

  #### Verfier le nom des coeffs LT et l'ordre : 1 constante, 2 coeff
  coeff_lt_value<-data.frame(coeff=c(cst_lt,names(xt)),value=LTcoeff)
  cat("\n \n Estimation du long-terme : \n")
  print(LT)
  return(coeff_lt_value)
}


# Fonction qui retourne le dataframe avec les coefficients de l'estimation +
# cree l'objet thoR.equation + envoie dans la console le detail de l'estimation
# Les cales doivent apparaitre à la fin de la formule !
# Aucun coefficient de doit apparaitre deux fois (sauf devant la cale)!
# Pas d'operations sur les coeffs dans la formule !
get_coeff_ct_value<-function(equation,coeff_ct,estim_start,estim_end,sub_database,index_time="date",const){
  # plage
  i_deb <- which(sub_database[,index_time]==estim_start)-equation@maxlag
  i_fin <- which(sub_database[,index_time]==estim_end)

  # formule
  formula_ols<-(paste(equation@LHS,equation@RHS,sep="~"))

  # donnees
  plop<-sub_database[i_deb:i_fin,]
  if(!const){formula_ols=paste0(formula_ols,"-1")}
  # regression MCO
  reg<-lm(eval(parse(text=formula_ols)), data=plop)
  cat("\n \n Estimation du court-terme : \n")
  print(summary(reg))
  # recuperer les coefficients
  my_coef <- coeff_ct
  my_l <- length(names(reg$coefficients))
  while (length(my_coef) < my_l) {
    my_coef <- c(my_coef,".")
  }
  names(reg$coefficients)<-my_coef
  coeff_ct_value<-data.frame(coeff=coeff_ct ,value=(stats::na.omit(reg$coefficients)))

  return(coeff_ct_value)
}

# Fonction qui fait l'estimation pour une equation
#' @title   Quickly estimate an OLS or ECM equation from a model
#'
#' @param thor_equation a thoR.equation object. Default : NULL
#' @param eq_name If no thor_equation is specified, name of the equation in the thor_model. Default : NULL
#' @param thor_model If no thor_equation is specified, the thoR.model where to get the equation eq_name. Default : NULL
#' @param database data.frame with the necessary data to estimate the equation. Default : NULL
#' @param index_time name of column in the database which is meant to be use as the time vector. Default : "date"
#' @param endogenous_name If no thor_equation is specified, endogenous variable of the equation, must be a variable of the equation. Default : NULL
#' @param coeff_lt list of long term coefficients. If NULL then the equation will only be estimated as OLS (one step). Default : NULL
#' @param estim_start First period of estimation, must be in index_time. Default : NULL
#' @param estim_end Last period of estimation, must be in index_time.. Default : NULL
#' @param const FALSE if the constant is not to be estimated. Default : TRUE
#'
#' @author Niamh Dunne
#' @import cointReg
#' @importFrom gsubfn gsubfn
#' @return database with the coefficients with all coefficients estimated
#' @export
#'
quick_estim <- function(thor_equation=NULL,eq_name=NULL,thor_model=NULL,endogenous_name=NULL,database,index_time="date",estim_start=NULL,estim_end=NULL,coeff_lt=NULL,const=TRUE){
  thor_equation_input <- FALSE
  if(is.null(thor_equation)){
    if(is.null(eq_name) | is.null(thor_model) | is.null(endogenous_name)){
      cat("\n If no thoR.equation is specified, please specify an equation from a thoR.model with an valid endogenous variable for this equation by specifying the arguments 'eq_name', 'thor_model' and 'endogenous_name'.")
      stop("No equation provided.")
    }else{
      ##check that this is valid
      if(inherits(thor_model,"thoR.model")==FALSE){stop("thor_model is not a thoR.model.")}
      if(!eq_name %in% thor_model@equation_list$name ){stop(paste(eq_name ,"cannot be found in the names of the thoR.model equations."))}
      formula_txt <- get_formula_txt(eq_name,thor_model)
      variables_list <- get_variables_list(formula_txt,thor_model)
      coeff_list <- get_coeff_list(formula_txt,variables_list)
    }

  }else{
    if(inherits(thor_equation,"thoR.equation")==FALSE){stop("thor_equation is not a thoR.equation.")}
      old_equation <- thor_equation
      thor_equation_input <- TRUE
      formula_txt <- thor_equation@formula
      formula_txt<-gsub("\\s+","",formula_txt)
      variables_list <- setdiff(get_variables_from_string(formula_txt), thor_equation@coefflist)
      eq_name<-thor_equation@equation_name
      endogenous_name <- thor_equation@endogenous
      coeff_list <-thor_equation@coefflist

  }

  cat(paste0("\n --------------- Estimating equation: ",eq_name," --------------- \n"))



  if(!is.null(coeff_lt)){if(prod(coeff_lt %in% coeff_list)==0){stop("Some coeff_lt were not found in the coefficients list of the equation.")}  }
  coeff_ct <- setdiff(coeff_list,coeff_lt)
  var_log_all <- get_var_log_all(formula_txt,variables_list)

  sub_database <- get_sub_database(database,var_log_all,variables_list,index_time)

  if (!is.null(coeff_lt)){
    var_lt <- get_var_lt(formula_txt,coeff_lt)
    coeff_lt_value <- get_coeff_lt_value(var_lt,var_log_all,endogenous_name,coeff_lt,sub_database,estim_start,estim_end,index_time)
    sub_database <- add_coeffs(listcoeff=coeff_lt_value,database=sub_database)
    database <- add_coeffs(listcoeff=coeff_lt_value,database=database)
    formula_ols <- formula_with_coeffs(formula_txt,coeff_lt,sub_database,round_digits=10)
  }else{
    formula_ols <- formula_txt
  }

  get_eq <- quietly(function(x) create_equation(x,
                                                formula = formula_ols,
                                                coefflist= coeff_list,
                                                endogenous= endogenous_name))

  coeff_ct_value <- get_coeff_ct_value(get_eq(eq_name)$result,coeff_ct,estim_start,estim_end,sub_database,index_time,const)
  database <- add_coeffs(listcoeff=coeff_ct_value,database=database)

  if(thor_equation_input){
    get_eq <- quietly(create_equation)(equation_name = old_equation@equation_name, formula = old_equation@formula,coefflist = old_equation@coefflist, endogenous = old_equation@endogenous  )
  }
  return(database)
}


# Fonction qui fait toutes les estimations par recursion, et renvoie le database avec tous les coefficients ajoutes
#' Quickly estimate multiple equations (OLS or ECM) from a model
#'
#' @param info_equations List (level 2) with the information about the equations. Each element (level 1) of the list must be a list named after an equation name of the model, that contains at least the following elements which correspond to the arguments of the `quick_estim()` function: endogenous_name ,coeff_lt, estim_start, estim_end and const.
#' @param thor_model thoR.model that contains the
#' @param database data.frame with the necessary data to estimate all the equations in the list. Default : NULL
#' @param index_time name of column in the database which is meant to be use as the time vector. Default : "date"
#' @author Niamh Dunne
#'
#' @return database with the coefficients with all coefficients estimated
#' @export
#' @importFrom purrr safely
#' @examples \dontrun{See the UK_exemple vignette for a comprehensive exemple.}
quick_estim_all <- function(info_equations,thor_model,database,index_time){
  n <- length(info_equations)
  if (n==0){
    return(database)
  }
  else{
    safe_quickestim<- purrr::safely(quick_estim,otherwise = NULL)

    res_quick_estim <- safe_quickestim(eq_name=names(info_equations)[1],
                            thor_model=thor_model,
                            database=database,
                            index_time=index_time,
                            endogenous_name=info_equations[[1]]$endogenous_name,
                            coeff_lt=info_equations[[1]]$coeff_lt,
                            estim_start=info_equations[[1]]$estim_start,
                            estim_end=info_equations[[1]]$estim_end,
                            const=info_equations[[1]]$const)

    if(is.null(res_quick_estim$error)){
      database<-res_quick_estim$result
    }else{
      cat(paste("Could not estimate equation:",names(info_equations)[1],"."))
      database<-database
    }
       return(quick_estim_all(info_equations[-1],
                           thor_model,
                           database,
                           index_time))


  }
}

