
#' Creating a thoR.equation
#'
#' @param equation_name character string to name the equations
#' @param formula character string of full formula
#' @param coefflist coefficient variables for information. Create equation will try to create a RHS and LHS that is compatible with econometric estimations
#' @param endogenous the endogenous variable of the equation (only one)
#' @param env environment where the the equation should be loaded
#' @param default_maxlag_error number of lag to default to if the maximum lag of the equation cannot be determined
#'
#' @return thoR.equation object
#' @export
#'
create_equation<-function(equation_name="my_equation",
                          formula,
                          coefflist=NULL,
                          endogenous=NULL,
                          env = globalenv(),
                          default_maxlag_error = 5

                          ){

  options(stringsAsFactors=FALSE)
  delta <- function(n,x, y=TRUE) NULL
  drule[["delta"]] <- alist(x=1, y=NULL) # y is just a logical
  ################################
  #### 0.a  Arguments verification
  ################################
  old_formula <- formula
  formula   <- tolower(gsub("\\s+","",formula))
  formula   <- gsub("mylg\\(","lag\\(",formula)
  coefflist   <- tolower(gsub("\\s+","",coefflist))
  endogenous   <- tolower(gsub("\\s+","",endogenous))

  assertthat::is.string(equation_name)
  if(acceptable_var_name(equation_name)==FALSE){stop("Invalid characters found in 'equation_name'.")}

  assertthat::assert_that(class(formula) == "character")
  if(parser_lag_delta_check(formula)==FALSE){stop("Parser error on the lags and/or delta. Please check the equation's formula specification.")}

  assertthat::assert_that(class(endogenous) == "character")
  if(acceptable_var_name(endogenous)==FALSE){stop("Invalid characters found in 'endogenous'.")}

  if (!is.null(coefflist)){assertthat::assert_that(class(coefflist) == "character")
                            if(acceptable_var_name(coefflist)==FALSE){stop("Invalid characters found in 'coefflist'.")}}

  if(length(equation_name)!=1){stop("Only one name can be specified.")}
  if(length(formula)!=1){stop("Only one equation can be specified.")}
  if(length(endogenous)!=1){stop("Only one endogenous variables can be specified.")}
  ################################
  #### 1  Variables check and extraction
  ################################
  econometric_safe <- TRUE
  eq_vars   <- get_variables_from_string(formula)
  eq_vars   <- unique(eq_vars[order(eq_vars)])
  coefflist <- unique(coefflist[order(coefflist)])

  wrongcharacters <- c("@|\\$|%|#|~|\\||;|:")

  if(acceptable_var_name(eq_vars)==FALSE){stop("Invalid characters found in the equation variables.")}
  if(sum(grepl(wrongcharacters,eq_vars))>0){stop("Invalid characters found in the formula")}
  if(stringr::str_count(formula,"=") != 1) {stop("The formula can only have one = sign")}
  if(length(setdiff(coefflist,eq_vars))>0){cat("The following coefficients were not found in the formula:")
                                           print(setdiff(coefflist,eq_vars))
                                           coefflist<-intersect(coefflist,eq_vars)}
  if(length(setdiff(endogenous,eq_vars))>0){stop(paste("The endogenous variable specified was not found in the formula."))}

  coefflist <- unique(coefflist[order(coefflist)])
  cat("\nEquation variables identified :\n")
  cat(eq_vars)
  cat("\n Exogenous variables :\n")
  cat(setdiff(eq_vars,c(coefflist,endogenous)))
  exogenous<- setdiff(eq_vars,c(coefflist,endogenous))
  exogenous<- exogenous[order(exogenous)]

  ################################
  #### 3  creating the formula
  ################################
  LHS<-(strsplit(formula,split="=")[[1]])[1]
  RHS<-(strsplit(formula,split="=")[[1]])[2]

  new_formula <- paste(LHS,RHS, sep = "-(")
  new_formula <- paste(new_formula, ")", sep = "")
  new_formula <- formatting_formulas(new_formula)

  LHS<-gsub("delta\\(1","newdiff\\(1" ,LHS)
  RHS<-gsub("delta\\(1","newdiff\\(1"  ,RHS)
  LHS<-gsub("(mylg|lag)\\((\\w+),-?([0-9]+)\\)","lag\\(\\2,\\3\\)"  ,LHS)
  RHS<-gsub("(mylg|lag)\\((\\w+),-?([0-9]+)\\)","lag\\(\\2,\\3\\)"  ,RHS)
  ################################
  #### 3  computing derivative
  ################################
  jacobian<-matrix(0,nrow = 1,ncol = 1)
  rownames(jacobian)<-equation_name
  colnames(jacobian)<-endogenous

  my_formula<- parse(text=new_formula)
  my_derivative <-Deriv(new_formula,endogenous)
  jacobian[1,1]<-paste(my_derivative)

  ################################
  #### 4  creating the functions
  ################################
  dir.create("temp_paprfn")
  cat(get_functions_commands(list_of_equations = new_formula,currentstep = "equation"),file=paste("temp_paprfn/eq_f",".R",sep=""),sep="\n",append = TRUE)
  cat(get_jacobian_commands(jacobian = jacobian,database_name = "t_data",output_matrix_name = "Jacobian_n",info.matrix = "jacobian",currentstep = "jacobian"),file=paste("temp_paprfn/j_f",".R",sep=""),sep="\n",append = TRUE)

  source("temp_paprfn/eq_f.R",local = TRUE)
  source("temp_paprfn/j_f.R",local = TRUE)
  unlink("temp_paprfn",recursive = TRUE)

  ################################
  #### 5 formatting for estimations later
  ################################

  coeff_in_LHS <- intersect(get_variables_from_string(LHS), coefflist)
  if(length(coeff_in_LHS)>0){econometric_safe <- FALSE
  cat("\n Some coefficient variables were found on the left-hand side of the equations. Coefficients will not reordered, beware if you are trying to estimate the coefficients econometrically.")
  ordered_coeff<-coefflist}else{
  #### 5.A Reordering coeffs in the order they appear if the equation is wrote for econometrics, otherwise cannot be used
  #### works if all coeffs are in the RHS
    econometric_checks<-reorder_coefficients(LHS,RHS,coefflist,eq_vars)
    if(econometric_checks[[1]]==FALSE){
      ordered_coeff<-coefflist
      econometric_safe<- FALSE
      RHS <-econometric_checks[[3]]
    }else{
      ordered_coeff<-econometric_checks[[2]]
      RHS <-econometric_checks[[3]]
    }
}
  ###Finding the lag maximum for the equation
  tesmaxlag<-find_maxlag(c(RHS,LHS))
  if(is.na(tesmaxlag) ){
    econometric_safe <- FALSE
    maxlag<-default_maxlag_error
  }else{
    maxlag<-tesmaxlag
  }
  ###finding max lag per equation
  ##If there are lags that are variable dependant, the equation is not safe for econometric estimations and maxlag cannot be computer, is set to 10.

  #########################################################################


  assign(equation_name,thoR.equation(equation_name=equation_name,
                                     formula   = old_formula,
                                     new_formula = new_formula,
                                     coefflist =ordered_coeff,
                                     endogenous=endogenous,
                                     exogenous = exogenous ,
                                     LHS=LHS,
                                     RHS=RHS,
                                     maxlag=maxlag,
                                     jacobian = jacobian,
                                     eq_f=equation_f,
                                     jac_f=jacobian_f,
                                     econometric_safe = econometric_safe),  envir = env )

}

