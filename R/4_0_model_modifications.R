################ENDO-EXO_SWITCH
###model_endo_exo_switch La fonction crée un nouveau modèle basé sur un modèle existant et fait passer les endogènes spécifiées en exogènes et exogènes spécifiées en exogènes, le même nombre d'endo et exo doit être renseigné.

#' Endo-Exo switch
#' Creates a new model base on an existing one with some endos and exos switched
#' @param base_model : the reference model to be used (obviously a thoR.model)
#' @param new_model_name  : character string, must be different from base model name
#' @param new_endo : vector of  exogenous variables from base model to be switched to endo
#' @param new_exo  : vector of  endogenous variables from base model to be switched to exo. length must match new_endo.
#' @param algo boolean. TRUE if using decomposition of model (recommended for large model). Default : TRUE
#' @param rcpp boolean. TRUE if using Rcpp for solver (recommended for large model). Default : FALSE
#' @param env Environment where to store the model object. Default : globalenv()
#' @param rcpp_path path to directory where to store the rcpp source files for the model. Default : working directory
#' @param use.superlu boolean. If SUPERLU library is installed, select TRUE to compile a rcpp model with superlu
#' @return a thor model in the global environment
#' @export
#' @import Deriv
#' @import assertthat
#' @import purrr
#' @import Rcpp
#' @import RcppArmadillo
#'
#'
model_endo_exo_switch<-function(base_model , new_model_name, new_endo, new_exo, algo = TRUE, rcpp = FALSE,env=globalenv(),rcpp_path = getwd() , use.superlu = FALSE) {
  new_endo<-unique(new_endo)
  new_exo<-unique(new_exo)
  ################################
  #### 0.a  New arguments
  ################################
  cat("Step 1: Running various checks on the inputs... \n")
  assertthat::is.string(new_model_name)
  if (new_model_name == base_model@model_name){stop("The new modified model must have a different name than the base model")}
  if (length(new_exo)!=length(new_exo)){stop("The number of endogenous and exogenous varibales switched must be equal")}
  if (prod(new_endo %in% base_model@exo_list)==0){stop("Some variables sepcified as new_endo are not exogenous variables of the base model")}
  if (prod(new_exo %in% base_model@endo_list)==0){stop("Some variables sepcified as new_exo are not endogenous variables of the base model")}
  if (rcpp){assertthat::is.dir(rcpp_path)}

###############################
#### 0.b Setup
################################

options(stringsAsFactors=FALSE)
delta <- function(n,x, y=TRUE) NULL
drule[["delta"]] <- alist(x=1, y=NULL) # y is just a logical


equations_list<-base_model@equation_list
################################
#### 1.a importing old model info
################################
endo  <- c(setdiff(base_model@endo_list,new_exo),new_endo)
exo   <- c(setdiff(base_model@exo_list,new_endo),new_exo)
coeff <- base_model@coeff_list
################################
#### 1.b Reordering variables
################################
endo<-endo[order(endo)]
exo<-exo[order(exo)]
coeff<-coeff[order(coeff)]

################################
#### 1.c Checks for inputs
################################
check_variable_conflict(endo , exo)
check_variable_conflict(endo , coeff)
check_variable_conflict(exo , coeff)

eqlist<-equations_list$equation
################################
#### 2.a Equations indexing and splitting
################################
exo  <- is_in_formulas(exo,eqlist,"exogenous")
endo <- is_in_formulas(endo,eqlist,"endogenous")
coeff<- is_in_formulas(coeff,eqlist,"coefficients")
all_model_variables <- c(endo,exo, coeff)
all_model_variables <- all_model_variables[order(all_model_variables)]
##from here on, steps should be written with much minimalism in order to be used easily with model modifier functions
##A/creating indexes names, and LHS RHS

#### 2.b generating the matrix that houses the contemporaneous endo
################################
cat("\n
  Step 2: Identifying the endogenous variables in the equations...
      \n")
##B/ generating eqns table
eqns <- table_contemporaneous_endos(formula_list = equations_list$formula,endogenous = endo,exogenous = exo, coefflist = coeff,equations_index = equations_list$id)

################################
#### 2.c checking parser and inputs
################################
check_eq_var_identification(endo = endo, eqns =eqns,names_of_equations = equations_list$name)

################################
#### 3.a decomposition
################################
if (algo == TRUE){
  cat("\n
  Step 3: Decomposing the model into blocks...
      \n")}else{
        cat("\n
  Step 3: Creating the unique block of the model...
      \n")
      }
## checks
if(algo == FALSE & nrow(eqns) > 75){cat("The model specified appears to be large. Consider using algo == TRUE to improve solver performances.")}
if(algo == TRUE  & nrow(eqns) < 50){cat("The model specified appears to be small. Consider using algo == FALSE to improve solver performances.")}

decomposition<-decomposing_model(endogenous_variables = endo,eq_var_matrix = eqns,decomposition = algo)
invisible(list2env(decomposition,environment()))

equations_list$part<-"tbd"
equations_list$part<-ifelse(equations_list$id %in% prologue_equations , "prologue",equations_list$part )
equations_list$part<-ifelse(equations_list$id %in% heart_equations , "heart",equations_list$part )
equations_list$part<-ifelse(equations_list$id %in% epilogue_equations , "epilogue",equations_list$part )

if("tbd" %in% equations_list$part){
  stop("Some equations were not identified in the decomposition.")
}else{cat("Decomposition algorithm successful. Below is the size of each block : \n")
  print(table(equations_list$part))}


################################
#### 4 Compute symbolic jacobians and the functions that go with
################################
cat("\n
   Step 4: Computing the blocks' symbolic jacobian matrices...
      \n")
eqns_as_list<-purrr::map(as.data.frame(t(eqns)),~unique(stats::na.omit(.x)))
cat("\n
   Step 5: Generating the model specific symbolic-to-numeric transition functions...
      \n")
prologue_jacobian <-matrix()
heart_jacobian <-matrix()
epilogue_jacobian <-matrix()

prologue_endo<-prologue_endo[order(prologue_endo)]
heart_endo<-heart_endo[order(heart_endo)]
epilogue_endo<-epilogue_endo[order(epilogue_endo)]

prologue_jacobian_f <-function(){print("none")}
heart_jacobian_f <-function(){print("none")}
epilogue_jacobian_f <-function(){print("none")}

prologue_equations_f <-function(){print("none")}
heart_equations_f <-function(){print("none")}
epilogue_equations_f <-function(){print("none")}

prologue_formulas <- equations_list[which(equations_list$id %in% prologue_equations),"new_formula" ]
heart_formulas <- equations_list[which(equations_list$id %in% heart_equations),"new_formula" ]
epilogue_formulas <- equations_list[which(equations_list$id %in% epilogue_equations),"new_formula" ]

dir.create("temp_paprfn",)

if(prologue == TRUE){
  prologue_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                    endo_vec = prologue_endo ,
                                                    equations_subset = prologue_equations,
                                                    eqns_vars_list = eqns_as_list)

  cat(get_functions_commands(list_of_equations = prologue_formulas,currentstep = "prologue_equations"),file=file.path("temp_paprfn","p_e_f.R"),sep="\n",append = TRUE)
  cat(get_jacobian_commands(jacobian = prologue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "prologue_jacobian"),file=file.path("temp_paprfn","p_j_f.R"),sep="\n",append = TRUE)
  source(file.path("temp_paprfn","p_e_f.R"),local = TRUE)
  source(file.path("temp_paprfn","p_j_f.R"),local = TRUE)
}

if(heart == TRUE){
  heart_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                 endo_vec = heart_endo ,
                                                 equations_subset = heart_equations,
                                                 eqns_vars_list = eqns_as_list)
  # heart_e_f_commands <- get_functions_commands(list_of_equations = heart_formulas,currentstep = "heart_equations")
  # heart_equations_f <- eval(parse(text=paste0(heart_e_f_commands,collapse = "\n")))
  #
  # heart_j_f_commands <- get_jacobian_commands(jacobian = heart_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian")
  # heart_jacobian_f <- eval(parse(text=paste0(heart_j_f_commands,collapse = "\n")))
  cat(get_functions_commands(list_of_equations = heart_formulas,currentstep = "heart_equations"),file=file.path("temp_paprfn","h_e_f.R"),sep="\n",append = TRUE)
  cat(get_jacobian_commands(jacobian = heart_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian"),file=file.path("temp_paprfn","h_j_f.R"),sep="\n",append = TRUE)
  source(file.path("temp_paprfn","h_e_f.R"),local = TRUE)
  source(file.path("temp_paprfn","h_j_f.R"),local = TRUE)
}

if(epilogue == TRUE){
  epilogue_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                    endo_vec = epilogue_endo ,
                                                    equations_subset = epilogue_equations,
                                                    eqns_vars_list = eqns_as_list)
  # epilogue_e_f_commands <- get_functions_commands(list_of_equations = epilogue_formulas,currentstep = "epilogue_equations")
  # epilogue_equations_f <- eval(parse(text=paste0(epilogue_e_f_commands,collapse = "\n")))
  #
  # epilogue_j_f_commands <- get_jacobian_commands(jacobian = epilogue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian")
  # epilogue_jacobian_f <- eval(parse(text=paste0(epilogue_j_f_commands,collapse = "\n")))
  cat(get_functions_commands(list_of_equations = epilogue_formulas,currentstep = "epilogue_equations"),file=file.path("temp_paprfn","e_e_f.R"),sep="\n",append = TRUE)
  cat(get_jacobian_commands(jacobian = epilogue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "epilogue_jacobian"),file=file.path("temp_paprfn","e_j_f.R"),sep="\n",append = TRUE)
  source(file.path("temp_paprfn","e_e_f.R"),local = TRUE)
  source(file.path("temp_paprfn","e_j_f.R"),local = TRUE)
}
################################
#### 5.b Rcpp files
################################
p_h_e_bool<-c(prologue,heart,epilogue)
p_h_e_jac <- list(prologue_jacobian , heart_jacobian , epilogue_jacobian )
rcpp_path<-gsub("/$","",rcpp_path)
rcpp_source <- "none"
if (rcpp == TRUE){
  cat("\n Creating the rcpp source files... \n")

  if(use.superlu){
    create_model_rcpp_source(
      model_name=new_model_name,p_h_e_bool,p_h_e_jac,all_model_vars = all_model_variables,rcpp_path = rcpp_path)
  }else{
    create_model_rcpp_source_nonsuperlu(
      model_name=new_model_name,p_h_e_bool,p_h_e_jac,all_model_vars = all_model_variables,rcpp_path = rcpp_path)
  }

  rcpp_source <-paste0(rcpp_path,"/",new_model_name,"_pomme_newton.cpp")
}
unlink("temp_paprfn",recursive = TRUE)
################################
#### 6. Create dico table
################################
cat("\n
  Step 6: Building the variables info map....
      \n")
var_map <- build_dico(equations_list_df = equations_list,endos = endo,exos = exo,coeffs = coeff,p_endos = prologue_endo,h_endos = heart_endo,e_endos = epilogue_endo)
################################
#### 7. Create the model object
################################


assign(new_model_name,thoR.model(model_name = new_model_name,
                                        equation_list = equations_list,
                                        var_map = var_map,
                                        variables_matrix = eqns,

                                        endo_list = tolower(endo) ,
                                        exo_list = tolower(exo),
                                        coeff_list = tolower(coeff),

                                        prologue = prologue ,
                                        heart = heart ,
                                        epilogue = epilogue,

                                        prologue_jacobian = prologue_jacobian ,
                                        heart_jacobian = heart_jacobian,
                                        epilogue_jacobian = epilogue_jacobian,

                                        prologue_endo = prologue_endo ,
                                        heart_endo = heart_endo ,
                                        epilogue_endo = epilogue_endo,

                                        algo=algo,
                                        rcpp = rcpp,
                                        rcpp_source=rcpp_source,

                                        prologue_equations_f=prologue_equations_f,
                                        heart_equations_f=heart_equations_f,
                                        epilogue_equations_f=epilogue_equations_f,

                                        prologue_jacobian_f=prologue_jacobian_f,
                                        heart_jacobian_f=heart_jacobian_f,
                                        epilogue_jacobian_f=epilogue_jacobian_f
),
envir = env )
if (rcpp == TRUE){
  cat("\n Compiling Rcpp functions... \n")
  Rcpp::sourceCpp(paste0(rcpp_path,"/",new_model_name,"_pomme_newton.cpp"))
}

cat("Model successfully built ! \n")
}##end of function

################REMOVE EQUATIONS

# model_remove_equation
#' Endo-Exo switch
#' Creates a new model base on an existing one with some endos and exos switched
#' @param base_model : the reference model to be used (obviously a thoR.model)
#' @param new_model_name  : character string, must be different from base model name
#' @param equations_to_remove : vector of equations to remove, selected by equation name or id
#' @param endos_to_remove  : vector of  endogenous variables from base model to remove. There must be as many as distinct equations, correspond to one equation each. If they still appear in the model, they will be considered as exogenous variables.
#' @param algo boolean. TRUE if using decomposition of model (recommended for large model). Default : TRUE
#' @param rcpp boolean. TRUE if using Rcpp for solver (recommended for large model). Default : FALSE
#' @param env Environment where to store the model object. Default : globalenv()
#' @param rcpp_path path to directory where to store the rcpp source files for the model. Default : working directory
#' @param use.superlu boolean. If SUPERLU library is installed, select TRUE to compile a rcpp model with superlu
#' @return a thor.model in the selected environment
#' @export
#' @import Deriv
#' @import assertthat
#' @import purrr
#' @import Rcpp
#' @import RcppArmadillo
#'
model_equations_remove<-function(base_model , new_model_name, equations_to_remove, endos_to_remove, algo = TRUE, rcpp = FALSE,env=globalenv(),rcpp_path = getwd() , use.superlu = FALSE) {
  equations_list_old<-base_model@equation_list
  ################################
  #### 0.a  New arguments
  ################################
  endos_to_remove<-unique(endos_to_remove)
  cat("Step 1: Running various checks on the inputs... \n")
  assertthat::is.string(new_model_name)

  if (new_model_name == base_model@model_name){stop("The new modified model must have a different name than the base model.")}

  ##check the equations id/name specified
  assertthat::assert_that(class(equations_to_remove)=="character")
  assertthat::assert_that(class(endos_to_remove)=="character")
  for (id in seq_along(equations_to_remove)){
     if (equations_to_remove[id] %in%equations_list_old$name){equations_to_remove[id]<- equations_list_old$id[which(equations_list_old$name==equations_to_remove[id])]}
  }
  if(sum(duplicated(equations_to_remove))>0){stop("Duplicated equations in the equations_to_remove vector.")}
  if(sum(!equations_to_remove %in%equations_list_old$id)>0){stop("Some equations specified were not found.")}
  ###check the endos
  if(sum(!endos_to_remove %in% base_model@endo_list)>0){stop("Some endogenous variables specified are not endogenous variables of the base model.")}
  if (length(endos_to_remove)!=length(equations_to_remove)){stop("The number of equations and endogenous variables removed must be equal.")}

  ##checking if each equation has an endo
  identification_matrix <-matrix(0,nrow=length(endos_to_remove),ncol = length(endos_to_remove))
  rownames(identification_matrix)<-equations_to_remove
  colnames(identification_matrix)<-endos_to_remove

  for (id in equations_to_remove){
    list_vars<-get_variables_from_string(equations_list_old$equation[which(equations_list_old$id==id)])
    endo_found<- intersect(list_vars,endos_to_remove)
    identification_matrix[id,endo_found]<-1
    }
  id_check <- c(rowSums(identification_matrix),colSums(identification_matrix))
  if(prod(id_check)==0){stop("Incoherent choices of endogenous variables to remove given the chosen equations to remove.")}

  if (rcpp){assertthat::is.dir(rcpp_path)}

  #removing the equations
  equations_list<-equations_list_old[-which(equations_list_old$id %in% equations_to_remove),]


  ###############################
  #### 0.b Setup
  ################################

  options(stringsAsFactors=FALSE)
  delta <- function(n,x, y=TRUE) NULL
  drule[["delta"]] <- alist(x=1, y=NULL) # y is just a logical

  ################################
  #### 1.a importing old model info
  ################################
  ##removed endos are temporarily assigned to exo variables. Model algorithms will remove unnecessary variables from the coef and exo list
  endo  <- c(setdiff(base_model@endo_list,endos_to_remove))
  exo   <- c(base_model@exo_list,endos_to_remove)
  coeff <- base_model@coeff_list
  ################################
  #### 1.b Reordering variables
  ################################
  endo<-endo[order(endo)]
  exo<-exo[order(exo)]
  coeff<-coeff[order(coeff)]

  ################################
  #### 1.c Checks for inputs
  ################################
  check_variable_conflict(endo , exo)
  check_variable_conflict(endo , coeff)
  check_variable_conflict(exo , coeff)

  eqlist<-equations_list$equation
  ################################
  #### 2.a Equations indexing and splitting
  ################################
  exo  <- is_in_formulas(exo,eqlist,"exogenous")
  endo <- is_in_formulas(endo,eqlist,"endogenous")
  coeff<- is_in_formulas(coeff,eqlist,"coefficients")
  all_model_variables <- c(endo,exo, coeff)
  all_model_variables <- all_model_variables[order(all_model_variables)]

  endo<-endo[order(endo)]
  exo<-exo[order(exo)]
  coeff<-coeff[order(coeff)]
  ##from here on, steps should be written with much minimalism in order to be used easily with model modifier functions
  ##A/creating indexes names, and LHS RHS

  #### 2.b generating the matrix that houses the contemporaneous endo
  ################################
  cat("\n
  Step 2: Identifying the endogenous variables in the equations...
      \n")
  ##B/ generating eqns table
  eqns <- table_contemporaneous_endos(formula_list = equations_list$formula,endogenous = endo,exogenous = exo, coefflist = coeff,equations_index = equations_list$id)

  ################################
  #### 2.c checking parser and inputs
  ################################
  check_eq_var_identification(endo = endo, eqns =eqns,names_of_equations = equations_list$name)

  ################################
  #### 3.a decomposition
  ################################
  if (algo == TRUE){
    cat("\n
  Step 3: Decomposing the model into blocks...
      \n")}else{
        cat("\n
  Step 3: Creating the unique block of the model...
      \n")
      }
  ## checks
  if(algo == FALSE & nrow(eqns) > 75){cat("The model specified appears to be large. Consider using algo == TRUE to improve solver performances.")}
  if(algo == TRUE  & nrow(eqns) < 50){cat("The model specified appears to be small. Consider using algo == FALSE to improve solver performances.")}

  decomposition<-decomposing_model(endogenous_variables = endo,eq_var_matrix = eqns,decomposition = algo)
  invisible(list2env(decomposition,environment()))

  equations_list$part<-"tbd"
  equations_list$part<-ifelse(equations_list$id %in% prologue_equations , "prologue",equations_list$part )
  equations_list$part<-ifelse(equations_list$id %in% heart_equations , "heart",equations_list$part )
  equations_list$part<-ifelse(equations_list$id %in% epilogue_equations , "epilogue",equations_list$part )

  if("tbd" %in% equations_list$part){
    stop("Some equations were not identified in the decomposition.")
  }else{cat("Decomposition algorithm successful. Below is the size of each block : \n")
    print(table(equations_list$part))}


  ################################
  #### 4 Compute symbolic jacobians and the functions that go with
  ################################
  cat("\n
   Step 4: Computing the blocks' symbolic jacobian matrices...
      \n")
  eqns_as_list<-purrr::map(as.data.frame(t(eqns)),~unique(stats::na.omit(.x)))
  cat("\n
   Step 5: Generating the model specific symbolic-to-numeric transition functions...
      \n")
  prologue_jacobian <-matrix()
  heart_jacobian <-matrix()
  epilogue_jacobian <-matrix()

  prologue_endo<-prologue_endo[order(prologue_endo)]
  heart_endo<-heart_endo[order(heart_endo)]
  epilogue_endo<-epilogue_endo[order(epilogue_endo)]

  prologue_jacobian_f <-function(){print("none")}
  heart_jacobian_f <-function(){print("none")}
  epilogue_jacobian_f <-function(){print("none")}

  prologue_equations_f <-function(){print("none")}
  heart_equations_f <-function(){print("none")}
  epilogue_equations_f <-function(){print("none")}

  prologue_formulas <- equations_list[which(equations_list$id %in% prologue_equations),"new_formula" ]
  heart_formulas <- equations_list[which(equations_list$id %in% heart_equations),"new_formula" ]
  epilogue_formulas <- equations_list[which(equations_list$id %in% epilogue_equations),"new_formula" ]

  dir.create("temp_paprfn",)

  if(prologue == TRUE){
    prologue_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                      endo_vec = prologue_endo ,
                                                      equations_subset = prologue_equations,
                                                      eqns_vars_list = eqns_as_list)

    cat(get_functions_commands(list_of_equations = prologue_formulas,currentstep = "prologue_equations"),file=file.path("temp_paprfn","p_e_f.R"),sep="\n",append = TRUE)
    cat(get_jacobian_commands(jacobian = prologue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "prologue_jacobian"),file=file.path("temp_paprfn","p_j_f.R"),sep="\n",append = TRUE)
    source(file.path("temp_paprfn","p_e_f.R"),local = TRUE)
    source(file.path("temp_paprfn","p_j_f.R"),local = TRUE)
  }

  if(heart == TRUE){
    heart_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                   endo_vec = heart_endo ,
                                                   equations_subset = heart_equations,
                                                   eqns_vars_list = eqns_as_list)
    # heart_e_f_commands <- get_functions_commands(list_of_equations = heart_formulas,currentstep = "heart_equations")
    # heart_equations_f <- eval(parse(text=paste0(heart_e_f_commands,collapse = "\n")))
    #
    # heart_j_f_commands <- get_jacobian_commands(jacobian = heart_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian")
    # heart_jacobian_f <- eval(parse(text=paste0(heart_j_f_commands,collapse = "\n")))
    cat(get_functions_commands(list_of_equations = heart_formulas,currentstep = "heart_equations"),file=file.path("temp_paprfn","h_e_f.R"),sep="\n",append = TRUE)
    cat(get_jacobian_commands(jacobian = heart_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian"),file=file.path("temp_paprfn","h_j_f.R"),sep="\n",append = TRUE)
    source(file.path("temp_paprfn","h_e_f.R"),local = TRUE)
    source(file.path("temp_paprfn","h_j_f.R"),local = TRUE)
  }

  if(epilogue == TRUE){
    epilogue_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                      endo_vec = epilogue_endo ,
                                                      equations_subset = epilogue_equations,
                                                      eqns_vars_list = eqns_as_list)
    # epilogue_e_f_commands <- get_functions_commands(list_of_equations = epilogue_formulas,currentstep = "epilogue_equations")
    # epilogue_equations_f <- eval(parse(text=paste0(epilogue_e_f_commands,collapse = "\n")))
    #
    # epilogue_j_f_commands <- get_jacobian_commands(jacobian = epilogue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian")
    # epilogue_jacobian_f <- eval(parse(text=paste0(epilogue_j_f_commands,collapse = "\n")))
    cat(get_functions_commands(list_of_equations = epilogue_formulas,currentstep = "epilogue_equations"),file=paste("temp_paprfn/e_e_f",".R",sep=""),sep="\n",append = TRUE)
    cat(get_jacobian_commands(jacobian = epilogue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "epilogue_jacobian"),file=paste("temp_paprfn/e_j_f",".R",sep=""),sep="\n",append = TRUE)
    source("temp_paprfn/e_e_f.R",local = TRUE)
    source("temp_paprfn/e_j_f.R",local = TRUE)
  }
  ################################
  #### 5.b Rcpp files
  ################################
  p_h_e_bool<-c(prologue,heart,epilogue)
  p_h_e_jac <- list(prologue_jacobian , heart_jacobian , epilogue_jacobian )
  rcpp_path<-gsub("/$","",rcpp_path)
  rcpp_source <- "none"
  if (rcpp == TRUE){
    cat("\n Creating the rcpp source files... \n")

    if(use.superlu){
      create_model_rcpp_source(
        model_name=new_model_name,p_h_e_bool,p_h_e_jac,all_model_vars = all_model_variables,rcpp_path = rcpp_path)
    }else{
      create_model_rcpp_source_nonsuperlu(
        model_name=new_model_name,p_h_e_bool,p_h_e_jac,all_model_vars = all_model_variables,rcpp_path = rcpp_path)
    }
    rcpp_source <-paste0(rcpp_path,"/",new_model_name,"_pomme_newton.cpp")
  }
  unlink("temp_paprfn",recursive = TRUE)
  ################################
  #### 6. Create dico table
  ################################
  cat("\n
  Step 6: Building the variables info map....
      \n")
  var_map <- build_dico(equations_list_df = equations_list,endos = endo,exos = exo,coeffs = coeff,p_endos = prologue_endo,h_endos = heart_endo,e_endos = epilogue_endo)
  ################################
  #### 7. Create the model object
  ################################


  assign(new_model_name,thoR.model(model_name = new_model_name,
                                              equation_list = equations_list,
                                              var_map = var_map,
                                              variables_matrix = eqns,

                                              endo_list = tolower(endo) ,
                                              exo_list = tolower(exo),
                                              coeff_list = tolower(coeff),

                                              prologue = prologue ,
                                              heart = heart ,
                                              epilogue = epilogue,

                                              prologue_jacobian = prologue_jacobian ,
                                              heart_jacobian = heart_jacobian,
                                              epilogue_jacobian = epilogue_jacobian,

                                              prologue_endo = prologue_endo ,
                                              heart_endo = heart_endo ,
                                              epilogue_endo = epilogue_endo,

                                              algo=algo,
                                              rcpp = rcpp,
                                              rcpp_source=rcpp_source,

                                              prologue_equations_f=prologue_equations_f,
                                              heart_equations_f=heart_equations_f,
                                              epilogue_equations_f=epilogue_equations_f,

                                              prologue_jacobian_f=prologue_jacobian_f,
                                              heart_jacobian_f=heart_jacobian_f,
                                              epilogue_jacobian_f=epilogue_jacobian_f
  ),
  envir = env )
  if (rcpp == TRUE){
    cat("\n Compiling Rcpp functions... \n")
    Rcpp::sourceCpp(paste0(rcpp_path,"/",new_model_name,"_pomme_newton.cpp"))
  }

  cat("Model successfully built ! \n")
}##end of function

################ADD EQUATIONS

#' Adds thor.equations to an existing model and builds the model.
#'
#' @param base_model : the reference model to be used (obviously a thoR.model)
#' @param new_model_name  : character string, must be different from base model name
#' @param thor_equations_add list of thoR equations to add.
#' @param algo boolean. TRUE if using decomposition of model (recommended for large model). Default : TRUE
#' @param rcpp boolean. TRUE if using Rcpp for solver (recommended for large model). Default : FALSE
#' @param env Environment where to store the model object. Default : globalenv()
#' @param rcpp_path path to directory where to store the rcpp source files for the model. Default : working directory
#' @param use.superlu boolean. If SUPERLU library is installed, select TRUE to compile a rcpp model with superlu
#' @return a thor.model in the selected environment
#' @export
#' @import Deriv
#' @import assertthat
#' @import purrr
#' @import Rcpp
#' @import RcppArmadillo
#'
model_equations_add<-function(base_model , new_model_name,
                              thor_equations_add= NULL,
                              algo = TRUE, rcpp = FALSE,env=globalenv(),rcpp_path = getwd() , use.superlu = FALSE ) {


  ###############################
  #### 0.b Setup
  ################################

  options(stringsAsFactors=FALSE)
  delta <- function(n,x, y=TRUE) NULL
  drule[["delta"]] <- alist(x=1, y=NULL) # y is just a logical

  equations_list_old<-base_model@equation_list
  endo_old <- base_model@endo_list
  exo_old <- base_model@exo_list
  coeff_old <- base_model@coeff_list

  ################################
  #### 1.a  New variables
  ################################
  thor_equations_add<-unique(thor_equations_add)
  if(prod(purrr::map_lgl(thor_equations_add,~class(.x)=="thoR.equation"))==0){stop("Equations added to the model must be of thoR.equation class.")}

  new_equations<-unique(purrr::map_chr(thor_equations_add,~.x@formula))


  new_endos<-unique(purrr::map_chr(thor_equations_add,~.x@endogenous))
  new_names<-unique(purrr::map_chr(thor_equations_add,~.x@equation_name))
  new_coeffs<-unique(unlist(purrr::map(thor_equations_add,~.x@coefflist)))


  new_formatted_formula<-tolower(gsub("\\s+|\\\n","",new_equations))
  new_formatted_formula<-unique(gsub("mylg\\(","lag\\(",new_formatted_formula))

  all_vars_new<-unique(unlist(purrr::map(new_formatted_formula,~get_variables_from_string(.x))))
  new_exos<-setdiff(setdiff(all_vars_new,new_endos),new_coeffs)

  if(sum(new_endos %in% endo_old) > 0){stop("New endogenous variables cannot be existing endogenous variables of the base model.")}
  if(length(new_endos) != length(thor_equations_add)){stop("There needs to be as many new endogenous variables as there are new equations.")}
  if(sum(new_endos %in% exo_old)>0){exo_old <- setdiff(exo_old, new_endos) }
  if(sum(new_endos %in% coeff_old)>0){coeff_old <- setdiff(coeff_old, new_endos) }
  if(sum(new_exos %in% endo_old)>0){new_exos <- setdiff(new_exos, endo_old) }
  if(sum(new_coeffs %in% endo_old)>0){new_coeffs <- setdiff(new_coeffs, endo_old) }
  if(sum(new_exos %in% exo_old)>0){new_exos <- setdiff(new_exos, exo_old) }
  if(sum(new_coeffs %in% coeff_old)>0){new_coeffs <- setdiff(new_coeffs, coeff_old) }

  ################################
  #### 1.b Compiling variables vectors
  ################################

  endo <- tolower(unique(c(endo_old,new_endos )))
  exo  <- tolower(unique(c(exo_old ,new_exos )))
  coeff <-tolower(unique(c(coeff_old ,new_coeffs )))

  endo<-endo[order(endo)]
  exo<-exo[order(exo)]
  coeff<-coeff[order(coeff)]



  ################################
  #### 1.c Compiling formulas
  ################################

  cat("Step 1: Combining new inputs and base model inputs... \n")
  assertthat::is.string(new_model_name)

  if (new_model_name == base_model@model_name){stop("The new modified model must have a different name than the base model.")}

  ##check the equations id/name specified
  if(length(new_names)!= length(new_formatted_formula)){stop("Some names or formulas of the new equations are the same.")}

  new_eqlist<-paste0(new_names,":",new_formatted_formula)
  new_old_eqlist<-ifelse(equations_list_old$name==equations_list_old$id,equations_list_old$equation,  paste0(equations_list_old$name,":",equations_list_old$equation))

  eqlist <- c(new_old_eqlist,new_eqlist)
  ###checks that all equations contain (=)
  check_equation_input(eqlist)

  ###checks that all elements in variable vectors are viable
  check_var_vector(endo,"Endogenous variables")
  check_var_vector(exo,"Exogenous variables")
  check_var_vector(coeff,"Coefficient variables")

  check_variable_conflict(endo , exo)
  check_variable_conflict(endo , coeff)
  check_variable_conflict(exo , coeff)

  if (rcpp){assertthat::is.dir(rcpp_path)}

  #Adding the equations to the list

  ################################
  #### 2.a Equations indexing and splitting
  ################################
  exo  <- is_in_formulas(exo,eqlist,"exogenous")
  endo <- is_in_formulas(endo,eqlist,"endogenous")
  coeff<- is_in_formulas(coeff,eqlist,"coefficients")
  all_model_variables <- c(endo,exo, coeff)
  all_model_variables <- all_model_variables[order(all_model_variables)]
  ##from here on, steps should be written with much minimalism in order to be used easily with model modifier functions
  ##A/creating indexes names, and LHS RHS
  equations_list<-create_equations_list(eqlist)
  if(parser_lag_delta_check(equations_list$equation)==FALSE){stop("Parser error on the lags and/or delta. Please check the model's formula specification")}

  memory  <- list(eqlist_raw = eqlist ,
                  endo_raw = endo ,
                  exo_raw = exo ,
                  coeff_raw = coeff)
  ################################
  #### 2.b generating the matrix that houses the contemporaneous endo
  ################################
  cat("\n
  Step 2: Identifying the endogenous variables in the equations...
      \n")
  ##B/ generating eqns table
  eqns <- table_contemporaneous_endos(formula_list = equations_list$formula,endogenous = endo,exogenous = exo, coefflist = coeff,equations_index = equations_list$id)

  ################################
  #### 2.c checking parser and inputs
  ################################
  check_eq_var_identification(endo = endo, eqns =eqns,names_of_equations = equations_list$name)

  ################################
  #### 3.a decomposition
  ################################
  if (algo == TRUE){
    cat("\n
  Step 3: Decomposing the model into blocks...
      \n")}else{
        cat("\n
  Step 3: Creating the unique block of the model...
      \n")
      }
  ## checks
  if(algo == FALSE & nrow(eqns) > 75){cat("The model specified appears to be large. Consider using algo == TRUE to improve solver performances.")}
  if(algo == TRUE  & nrow(eqns) < 50){cat("The model specified appears to be small. Consider using algo == FALSE to improve solver performances.")}

  decomposition<-decomposing_model(endogenous_variables = endo,eq_var_matrix = eqns,decomposition = algo)
  invisible(list2env(decomposition,environment()))

  equations_list$part<-"tbd"
  equations_list$part<-ifelse(equations_list$id %in% prologue_equations , "prologue",equations_list$part )
  equations_list$part<-ifelse(equations_list$id %in% heart_equations , "heart",equations_list$part )
  equations_list$part<-ifelse(equations_list$id %in% epilogue_equations , "epilogue",equations_list$part )

  if("tbd" %in% equations_list$part){
    stop("Some equations were not identified in the decomposition.")
  }else{cat("Decomposition algorithm successful. Below is the size of each block : \n")
    print(table(equations_list$part))}

  ################################
  #### 3.b modifying formulas for easy jacobian transformation
  ################################
  equations_list$new_formula <- formatting_formulas(equations_list$formula)

  ################################
  #### 4 Compute symbolic jacobians and the functions that go with
  ################################
  cat("\n
   Step 4: Computing the blocks' symbolic jacobian matrices...
      \n")
  eqns_as_list<-purrr::map(as.data.frame(t(eqns)),~unique(stats::na.omit(.x)))
  cat("\n
   Step 5: Generating the model specific symbolic-to-numeric transition functions...
      \n")
  prologue_jacobian <-matrix()
  heart_jacobian <-matrix()
  epilogue_jacobian <-matrix()

  prologue_endo<-prologue_endo[order(prologue_endo)]
  heart_endo<-heart_endo[order(heart_endo)]
  epilogue_endo<-epilogue_endo[order(epilogue_endo)]

  prologue_jacobian_f <-function(){print("none")}
  heart_jacobian_f <-function(){print("none")}
  epilogue_jacobian_f <-function(){print("none")}

  prologue_equations_f <-function(){print("none")}
  heart_equations_f <-function(){print("none")}
  epilogue_equations_f <-function(){print("none")}

  prologue_formulas <- equations_list[which(equations_list$id %in% prologue_equations),"new_formula" ]
  heart_formulas <- equations_list[which(equations_list$id %in% heart_equations),"new_formula" ]
  epilogue_formulas <- equations_list[which(equations_list$id %in% epilogue_equations),"new_formula" ]

  dir.create("temp_paprfn",)

  if(prologue == TRUE){
    prologue_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                      endo_vec = prologue_endo ,
                                                      equations_subset = prologue_equations,
                                                      eqns_vars_list = eqns_as_list)

    cat(get_functions_commands(list_of_equations = prologue_formulas,currentstep = "prologue_equations"),file=file.path("temp_paprfn","p_e_f.R"),sep="\n",append = TRUE)
    cat(get_jacobian_commands(jacobian = prologue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "prologue_jacobian"),file=file.path("temp_paprfn","p_j_f.R"),sep="\n",append = TRUE)
    source(file.path("temp_paprfn","p_e_f.R"),local = TRUE)
    source(file.path("temp_paprfn","p_j_f.R"),local = TRUE)
  }

  if(heart == TRUE){
    heart_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                   endo_vec = heart_endo ,
                                                   equations_subset = heart_equations,
                                                   eqns_vars_list = eqns_as_list)
    # heart_e_f_commands <- get_functions_commands(list_of_equations = heart_formulas,currentstep = "heart_equations")
    # heart_equations_f <- eval(parse(text=paste0(heart_e_f_commands,collapse = "\n")))
    #
    # heart_j_f_commands <- get_jacobian_commands(jacobian = heart_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian")
    # heart_jacobian_f <- eval(parse(text=paste0(heart_j_f_commands,collapse = "\n")))
    cat(get_functions_commands(list_of_equations = heart_formulas,currentstep = "heart_equations"),file=file.path("temp_paprfn","h_e_f.R"),sep="\n",append = TRUE)
    cat(get_jacobian_commands(jacobian = heart_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian"),file=file.path("temp_paprfn","h_j_f.R"),sep="\n",append = TRUE)
    source(file.path("temp_paprfn","h_e_f.R"),local = TRUE)
    source(file.path("temp_paprfn","h_j_f.R"),local = TRUE)
  }

  if(epilogue == TRUE){
    epilogue_jacobian <- symbolic_jacobian(equations_list_df = equations_list ,
                                                      endo_vec = epilogue_endo ,
                                                      equations_subset = epilogue_equations,
                                                      eqns_vars_list = eqns_as_list)
    # epilogue_e_f_commands <- get_functions_commands(list_of_equations = epilogue_formulas,currentstep = "epilogue_equations")
    # epilogue_equations_f <- eval(parse(text=paste0(epilogue_e_f_commands,collapse = "\n")))
    #
    # epilogue_j_f_commands <- get_jacobian_commands(jacobian = epilogue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "heart_jacobian")
    # epilogue_jacobian_f <- eval(parse(text=paste0(epilogue_j_f_commands,collapse = "\n")))
    cat(get_functions_commands(list_of_equations = epilogue_formulas,currentstep = "epilogue_equations"),file=file.path("temp_paprfn","e_e_f.R"),sep="\n",append = TRUE)
    cat(get_jacobian_commands(jacobian = epilogue_jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "epilogue_jacobian"),file=file.path("temp_paprfn","e_j_f.R"),sep="\n",append = TRUE)
    source(file.path("temp_paprfn","e_e_f.R"),local = TRUE)
    source(file.path("temp_paprfn","e_j_f.R"),local = TRUE)
  }
  ################################
  #### 5.b Rcpp files
  ################################
  p_h_e_bool<-c(prologue,heart,epilogue)
  p_h_e_jac <- list(prologue_jacobian , heart_jacobian , epilogue_jacobian )
  rcpp_path<-gsub("/$","",rcpp_path)
  rcpp_source <- "none"
  if (rcpp == TRUE){
    cat("\n Creating the rcpp source files... \n")

    if(use.superlu){
      create_model_rcpp_source(
        model_name=new_model_name,p_h_e_bool,p_h_e_jac,all_model_vars = all_model_variables,rcpp_path = rcpp_path)
    }else{
      create_model_rcpp_source_nonsuperlu(
        model_name=new_model_name,p_h_e_bool,p_h_e_jac,all_model_vars = all_model_variables,rcpp_path = rcpp_path)
    }


    rcpp_source <-paste0(rcpp_path,"/",new_model_name,"_pomme_newton.cpp")
  }
  unlink("temp_paprfn",recursive = TRUE)
  ################################
  #### 6. Create dico table
  ################################
  cat("\n
  Step 6: Building the variables info map....
      \n")
  var_map <- build_dico(equations_list_df = equations_list,endos = endo,exos = exo,coeffs = coeff,p_endos = prologue_endo,h_endos = heart_endo,e_endos = epilogue_endo)
  ################################
  #### 7. Create the model object
  ################################


  assign(new_model_name,thoR.model(model_name = new_model_name,
                                              equation_list = equations_list,
                                              var_map = var_map,
                                              variables_matrix = eqns,

                                              endo_list = tolower(endo) ,
                                              exo_list = tolower(exo),
                                              coeff_list = tolower(coeff),

                                              prologue = prologue ,
                                              heart = heart ,
                                              epilogue = epilogue,

                                              prologue_jacobian = prologue_jacobian ,
                                              heart_jacobian = heart_jacobian,
                                              epilogue_jacobian = epilogue_jacobian,

                                              prologue_endo = prologue_endo ,
                                              heart_endo = heart_endo ,
                                              epilogue_endo = epilogue_endo,

                                              algo=algo,
                                              rcpp = rcpp,
                                              rcpp_source=rcpp_source,

                                              prologue_equations_f=prologue_equations_f,
                                              heart_equations_f=heart_equations_f,
                                              epilogue_equations_f=epilogue_equations_f,

                                              prologue_jacobian_f=prologue_jacobian_f,
                                              heart_jacobian_f=heart_jacobian_f,
                                              epilogue_jacobian_f=epilogue_jacobian_f
  ),
  envir = env )
  if (rcpp == TRUE){
    cat("\n Compiling Rcpp functions... \n")
    Rcpp::sourceCpp(paste0(rcpp_path,"/",new_model_name,"_pomme_newton.cpp"))
  }

  cat("Model successfully built ! \n")
}##end of function
