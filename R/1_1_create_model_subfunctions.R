#' Getting the old style list2DF
#'
#' @param xlist a list
#' @return dataframe with NAs when the length of the elements differ
#' @keywords internal
list_2DF_withNA <- function(xlist){
  maximus <- max(length(xlist))
  res <- xlist %>% map(~if(length(.x != maximus)){.x <- c(.x, rep(NA,(maximus - length(.x)) ))})
}


#' Generating equation index id
#'
#' @param equation_list a vector
#' @return character vector with eq_###
#' @keywords internal
generate_equations_index <- function(equation_list){
  n_eqlist <- length(equation_list)

  num_vec <- c(1:n_eqlist)
  size_char <- nchar(as.character(n_eqlist))

  index_char <- formatC(num_vec, width = size_char, format = "d", flag = "0")
  equations_index<-paste("eq",index_char,sep = "_")

  }


#' Create the equations_list table from an eqlist
#'
#' @param eq_list vector of equations as raw input
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom splitstackshape cSplit
#' @return data.frame equations lists
#' @keywords internal
create_equations_list<-function(eq_list){
  eqlist<-tolower(gsub("\\s+|\\\n","",eq_list) )


  indexing <- as.data.frame(splitstackshape::cSplit(as.data.frame(eqlist),1,sep = ":" ,type.convert = FALSE))
  if (ncol(indexing)==1){indexing[,2]<-NA_character_}

  names(indexing) <- c("equation_name","equation")
  indexing$name_mod <- !is.na(indexing$equation)
  indexing$id <- generate_equations_index(eqlist)
  indexing$name <- ifelse(indexing$name_mod == TRUE , indexing$equation_name,indexing$id )
  indexing$equation <- ifelse(indexing$name_mod == FALSE , indexing$equation_name,indexing$equation )
  # check for conflicts in name duplicates and if names match incorrect id
  id <- indexing$id
  name <- indexing$name

  ##duplicates correction
  if (sum(duplicated(name))>0){
    cat("Some equations names are duplicates. They will be altered. To avoid this, modify the names in the model input.\n")
    name <- make.unique(name,sep = "_")
  }


  ##conflict corretion
  name_in_id <- name %in% id
  name_is_id <- name == id
  name_conflict_id <- name_in_id + name_is_id == 1

  if(sum(name_conflict_id)>0){
    cat("The following equation names clash with existing ids, they will revert to id : \n")
    print(name[which(name_conflict_id==1)])
    name[which(name_conflict_id==1)] <- id[which(name_conflict_id==1)]
  }

  check_var_vector(name)
  indexing$name <- name
  indexing <- indexing %>% dplyr::select("id","name","equation")
  equations_list <- as.data.frame(splitstackshape::cSplit(indexing,"equation",sep = "=" ,type.convert = FALSE,drop = FALSE)) %>% dplyr::rename(LHS = "equation_1", RHS = "equation_2")
  equations_list$formula<-paste(equations_list$LHS,equations_list$RHS, sep = "-(")
  equations_list$formula<-paste(equations_list$formula, ")", sep = "")
  rownames(equations_list)<- equations_list$id
  return<-equations_list
}


#' Creating the table of contemporenous endos
#'
#' @param formula_list vector containing all formulas
#' @param endogenous vector of endogenous variables
#' @param exogenous vector of exogenous variableq
#' @param coefflist vector of coefficient
#' @param equations_index vector of equation index
#' @param functionsthor vector of supported functions : defaults to thor_functions_supported
#'
#' @return data.frame with all contemporaneous endos per equation
#' @importFrom splitstackshape cSplit
#' @importFrom purrr %>%
#' @importFrom purrr map_dbl
#' @importFrom purrr map2_df
#' @importFrom methods new
#' @keywords internal
table_contemporaneous_endos<-function(formula_list,endogenous,exogenous, coefflist,equations_index,functionsthor = thor_functions_supported){
  endo  <- endogenous
  exo   <- exogenous
  coeff <- coefflist


  function_pattern<- c(functionsthor,toupper(functionsthor)) %>% paste0('\\(') %>% paste(collapse = "|")
  function_pattern<-paste0("-|\\*|/|\\^|\\(|\\)|",function_pattern)


  eqns<-gsub("(mylg|lag)\\(\\w+,(-)?([0-9]|\\(\\w+\\+[0-9]\\)|\\(\\w+\\))\\)","LAGFLAG",formula_list)
  eqns<-gsub("delta\\(([0-9]+),","delta\\1\\(",eqns)
  eqns<-gsub("delta[0-9]+","+",eqns)

  eqns<-gsub(function_pattern,"+",eqns)
  eqns<-gsub("[0-9]+\\.[0-9]+e","+",eqns)
  eqns<-gsub("\\++","+",eqns)
  eqns<-gsub("^\\+","",eqns)

  #### 3.B Creating the matrix that lists all the variables per equation
  eqns <-splitstackshape::cSplit(as.matrix(eqns),1, sep= "+", drop=TRUE , type.convert=FALSE  )
  eqns<-as.data.frame(apply(eqns,1:2,function(x)gsub("^[0-9]+\\.?[0-9]*$|LAGFLAG","",x)),stringsAsFactors=FALSE,row.names = equations_index)
  eqns<-as.data.frame(apply(eqns,1:2,function(x)gsub('\\s+', "",x)),stringsAsFactors=FALSE,row.names = equations_index)
  ## eqns contains only contemporaneous variables
  eqns[eqns==""]=NA

  variables_matrix<-as.matrix(eqns)
  ## keep only endogenous variables
  remove_those<- c(exo, coeff)
  if(length(remove_those)>0){
       variables_matrix[variables_matrix %in% remove_those] <- NA}

  #remove unnecessary columns
  eqns<-Filter(function(x) !all(is.na(x)), as.data.frame(variables_matrix))
  simple_eqns<-apply(eqns, 1, paste, collapse="+")
  simple_eqns<-gsub("\\+NA","",simple_eqns)
  simple_eqns<-gsub("NA\\+","",simple_eqns)
  eqns <- as.data.frame(splitstackshape::cSplit(as.data.frame(simple_eqns),c(1) ,sep= "+", drop=TRUE , type.convert=FALSE  ))
  eqns_list <- strsplit(simple_eqns,"\\+") %>% purrr::map(~unique(.x))

  na_needed<- function(vector,target_length){
    size <-length(vector)
    na_needed <- target_length - size
  }
  nas_vec <- eqns_list %>% purrr::map_dbl(~na_needed(.x,target_length = max(lengths(eqns_list)) ) )

  fill_with_na<-function(vector,n_nas){
    if (n_nas > 0){k <- length(vector) - n_nas + 1
                  vector[k:length(vector)] <- NA_character_}
    vector
  }


  eqns_bis <- list_2DF_withNA(eqns_list)
  eqns_ter <- as.data.frame(t(purrr::map2_df(eqns_bis,nas_vec,~fill_with_na(.x,.y))))



  return(eqns_ter)
}


#' Formatting formula
#' rearranges logs and lags to fit  derivators and such.
#'
#' @param formula_list vector of formula to be modified
#'
#' @return character vector of reformatted formulas easy for computation of jacobians
#' @keywords internal
formatting_formulas <- function(formula_list){
  source_formula <- formula_list

  f_mod <- gsub("(mylg|lag)\\((\\w+),-?0\\)|(mylg|lag)\\((LOG\\(\\w+\\)),-?0\\)","\\2",source_formula)
  f_mod <- gsub("(mylg|lag)\\((LOG\\(\\w+\\)),-?0\\)","\\2",f_mod)
  f_mod <- gsub("(mylg|lag)\\((log\\(\\w+\\)),-?0\\)","\\2",f_mod)
  f_mod <- gsub("(mylg|lag)\\((\\w+),-?(\\w+)\\)","lag\\.\\2\\.\\3",f_mod)
  f_mod <- gsub("(mylg|lag)\\((\\w+),-?\\((\\w+)\\)\\)","lag\\.\\2\\.\\3",f_mod)
  f_mod <- gsub("(mylg|lag)\\((\\w+),-?\\((\\w+)\\+(\\w+)\\)\\)","lag\\.\\2\\.\\3\\.\\4",f_mod)
  f_mod <- gsub("(mylg|lag)\\((LOG)\\((\\w+)\\),-?(\\w+)\\)","lag\\.\\2\\.\\3\\.\\4",f_mod)
  f_mod <- gsub("(mylg|lag)\\((LOG)\\((\\w+)\\),-?\\((\\w+)\\)\\)","lag\\.\\2\\.\\3\\.\\4",f_mod)
  f_mod <- gsub("(mylg|lag)\\((LOG)\\((\\w+)\\),-?\\((\\w+)\\+(\\w+)\\)\\)","lag\\.\\2\\.\\3\\.\\4.\\4",f_mod)
  f_mod <- gsub("LOG","log",f_mod)
  f_mod <- gsub("lag\\.log\\.(\\w+)\\.(\\w+)\\.(\\w+)","log\\(lag\\.\\1\\.\\2\\.\\3\\)",f_mod)
  f_mod <- gsub("lag\\.log\\.(\\w+)\\.(\\w+)","log\\(lag\\.\\1\\.\\2\\)",f_mod)

  return(f_mod)
}


#' Determine all variables present in a formula
#'
#' @param string character string of the formula
#'
#' @return character vector with all variables detected in the formula.
#' @export
get_variables_from_string <- function(string){
  string<-gsub("\\s+","",string)
  function_pattern<- c(thor_functions_supported,toupper(thor_functions_supported),"delta","lag","mylg") %>% paste0('\\(') %>% paste(collapse = "|")
  symbol_pattern <- "(\\+|\\*|,|-|/|\\^|\\(|\\)|=)|\\\\"

  transfo1<-gsub(function_pattern,"@",string)
  transfo2<-gsub(symbol_pattern,"@",transfo1)
  transfo3<-gsub("@+",",",transfo2)
  transfo4<-gsub(",[0-9]+(\\.[0-9]+)?",",",transfo3)
  transfo4<-gsub("^[0-9]+(\\.[0-9]+)?","",transfo4)
  transfo5<-gsub("^,|,$","",transfo4)
  plop<-  unique(as.vector(strsplit(transfo5,",",))[[1]])
  plop<- plop[!plop%in%""]

  plop
}

