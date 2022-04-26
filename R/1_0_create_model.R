# #New Create Model
# model_name = "model"
# model_source = "essais/opale_prevision_lag.txt"
# algo = TRUE
# rcpp = TRUE
# endogenous = NULL
# exogenous = NULL
# coefficients = NULL
# equations = NULL
# env= globalenv()
# rcpp_path = getwd()
#
#' @title Create a thoR.model object
#' @description This function creates a thoR.model object by reading a model in text file or specified from a manual input in the function arguments.
#'
#' @param model_name character. Name of the model object. Default: 'prev'
#' @param model_source path. path to the .txt file that contains the model info. Default: 'NULL'
#' @param algo boolean. TRUE if using decomposition of model (recommended for large model). Default : TRUE
#' @param rcpp boolean. TRUE if using Rcpp for solver (recommended for large model). Default: FALSE
#' @param endogenous character vector of endogenous variables (only read if no input is provided in model_source). Default: NULL
#' @param exogenous character vector of exogenous variables (only read if no input is provided in model_source, optionnal). Default: NULL
#' @param coefficients character vector of coefficients (only read if no input is provided in model_source, optionnal). Default: NULL
#' @param equations character vector containing equations (only read if no input is provided in model_source, optionnal). Default: NULL
#' @param env Environment where to store the model object. Default: globalenv()
#' @param rcpp_path path to directory where to store the rcpp source files for the model. Default: working directory
#'
#' @return A thoR.model created in the environment.
#' @examples
#' \dontrun{
#'  create_model("My_Model",model_source_example(TRUE))
#'  create_model("Opale",system.file("models", "opale.txt", package = "tresthor"),rcpp=TRUE)
#'  }
#'
#' @rdname create_model
#' @export
#' @import Deriv
#' @import assertthat
#' @import purrr
#' @import Rcpp
create_model <- function(model_name = "model",
                         model_source = NULL,
                         algo = TRUE,
                         rcpp = FALSE,
                         rcpp_path = getwd(),
                         endogenous = NULL,
                         exogenous = NULL,
                         coefficients = NULL,
                         equations = NULL, env= globalenv()) {
  ################################
  #### 0.a  Arguments verification
  ################################
  cat("Let's build the model :)
      \n
  Step 1: Running various checks on the inputs... \n\n")

  assertthat::is.string(model_name)
  if (!is.null(model_source)) {
    assertthat::is.readable(model_source)
    check_model_file(model_source)
    cat(paste0("The model will be created using the file ",basename(model_source),".\nIgnoring all other equations specifications."))

    manual_model_input = FALSE}

  if (is.null(model_source)) {
    manual_model_input = TRUE

    if (!is.null(endogenous) & !is.null(equations)) {
      assertthat::assert_that(class(endogenous) == "character")
      assertthat::assert_that(class(equations) == "character")
      if (!is.null(exogenous)){assertthat::assert_that(class(exogenous) == "character")}
      if (!is.null(coefficients)){assertthat::assert_that(class(coefficients) == "character")}
      }else{
        stop(paste0("No model source is provided and no equation is specified.
Please either provide an external model file as model_source, or provide the equations and their specifications
(endogenous variables, and if needed exogenous variables and coefficients) manually in the appropriate arguments.
Refer to the documentation for more information."))
      }
  }
  assertthat::assert_that(class(algo) == "logical")
  assertthat::assert_that(class(rcpp) == "logical")

  if (rcpp){assertthat::is.dir(rcpp_path)}
  ################################
  #### 0.b Setup
  ################################

  options(stringsAsFactors=FALSE)
  delta <- function(n,x, y=TRUE) NULL
  drule[["delta"]] <- alist(x=1, y=NULL) # y is just a logical


  ################################
  #### 1.a Reading from outside source
  ################################
  if ( manual_model_input == FALSE ){
    model_input <- readLines(model_source,warn = FALSE)
    eqlist <- unique(model_input[11:length(model_input)])
    endo <- tolower(unique(strsplit(model_input[2],split = ",")[[1]]))
    exo <- tolower(unique(strsplit(model_input[5],split = ",")[[1]]))
    coeff <-tolower(unique(strsplit(model_input[8],split = ",")[[1]]))
    }
  ################################
  #### 1.b Reading from manual input
  ################################
  if ( manual_model_input == TRUE ){
    eqlist<-unique(equations)
    endo <- tolower(unique(endogenous))
    exo <- tolower(unique(exogenous))
    coeff <- tolower(unique(coefficients))
  }

  ################################
  #### 1.b Reading from manual input
  ################################
  endo<-endo[order(endo)]
  exo<-exo[order(exo)]
  coeff<-coeff[order(coeff)]

  ################################
  #### 1.c Checks for inputs
  ################################
  eqlist   <- gsub("mylg\\(","lag\\(",eqlist)
  ###checks that all equations contain (=)
  check_equation_input(eqlist)

  ###checks that all elements in variable vectors are viable
  check_var_vector(endo,"Endogenous variables")
  check_var_vector(exo,"Exogenous variables")
  check_var_vector(coeff,"Coefficient variables")

  check_variable_conflict(endo , exo)
  check_variable_conflict(endo , coeff)
  check_variable_conflict(exo , coeff)
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
  if(parser_lag_delta_check(equations_list$equation)==FALSE){stop("Parser error on the lags and/or delta. Please check the model's formula specification.")}

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
    rcpp_path<-normalizePath(rcpp_path)
    rcpp_source <- "none"
  if (rcpp == TRUE){
    cat("\n Creating the rcpp source files... \n")
    create_model_rcpp_source(model_name,p_h_e_bool,p_h_e_jac,all_model_vars = all_model_variables,rcpp_path = rcpp_path)
    rcpp_source <- file.path(rcpp_path,paste0(model_name,"_pomme_newton.cpp"))
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


  assign(model_name,thoR.model(model_name = model_name,
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
    Rcpp::sourceCpp(file.path(rcpp_path,paste0(model_name,"_pomme_newton.cpp")))
  }

 cat(" Model successfully built ! \n")
  }##end of function
