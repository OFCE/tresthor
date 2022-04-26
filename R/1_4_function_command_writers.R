#' Generate the lines of the function for evaluating equations
#' stored as a vector for easy manipulation
#' @param list_of_equations vector of equations
#' @param database_name database of the solver (t_data)
#' @param output_vector_name vector name as used in solver f_x_n
#' @param currentstep will be used in the function name
#' @importFrom purrr %>%
#' @return character vector with the lines to build the function
#' @keywords internal
get_functions_commands<-function(list_of_equations,
                                 database_name ="t_data",
                                 output_vector_name="f_x_n",
                                 currentstep="prologue_equations"){

  vectorsize <- length(list_of_equations)
  plop_col<-vector()

  ###Line 1 should read : current_step_f<-function(){
  line1<- paste(currentstep,"_f<-" ,sep="")
  line1a<- gsub("<-","<-function\\(t=1,t_data\\)\\{",line1)
  plop_col<-c(plop_col,line1a)

  ###Line 2 creates the zero matrix  should read :  Jacobian_n<-matrix(0,ncol = 201,nrow = 201)
  line2 <- paste(output_vector_name,"<- cBLERGrepBLERG0,",vectorsize,"GRELB","GRELB",sep = "")
  line2a <- gsub("BLERG","\\(",line2)
  line2b <- gsub("GRELB","\\)",line2a)
  plop_col<-c(plop_col,line2b)

  ###
  for ( i in 1:length(list_of_equations) )
  { functions_supported<-paste0(thor_functions_supported,"\\(") %>% paste(collapse="|")
  function_patterns<- paste0("(",functions_supported,"|lag\\.|delta\\(|newdiff\\(|\\.e[0-9]\\()")

  old_expression<-tolower(list_of_equations[i])

  expression_functions<-gsub(function_patterns,"\\U\\1",old_expression,perl = TRUE)
  ##change delta to newdiff
  delta_expression<-gsub("DELTA\\(([0-9]+),","NEWDIFF\\(\\1,",expression_functions,perl = TRUE)
  ##deal with lags
  lag_expression<-gsub("LAG\\.(\\w+)\\.(\\w+)\\.(\\w+)","LAG\\(\\1,\\U\\2\\+\\U\\3\\)",delta_expression,perl = TRUE)
  lag_expression<-gsub("LAG\\.(\\w+)\\.(\\w+)","LAG\\(\\1,\\U\\2)",lag_expression,perl = TRUE)
  lag_expression_bis<-gsub("(LAG\\(\\w+,)(\\w+)\\+(\\w+)\\)",paste("\\1",toupper(database_name),"\\[T,\\'\\2\\'\\]\\+",toupper(database_name),"\\[T,\\'\\3\\'\\]\\)",sep=""),lag_expression,perl = TRUE)
  lag_expression_bis<-gsub("(LAG\\(\\w+,)(\\w+)\\)",paste("\\1",toupper(database_name),"\\[T,\\'\\2\\'\\]\\)",sep=""),lag_expression_bis,perl = TRUE)
  lag_expression_ter<-gsub(paste(toupper(database_name),"\\[T,\\'([0-9]+)'\\]",sep=""),"\\1",lag_expression_bis,perl = TRUE)
  ##replace all variables with database_name[,variable_name]
  new_expression1<-gsub("([a-z]([a-z0-9_]+)?)",paste(database_name,"\\[,\\'\\1\\'\\]",sep=""),lag_expression_ter)
  new_expression2<-tolower(new_expression1)
  ##add the [t] at the end and big parenthis for the whole expression
  new_expression3<- if(grepl(database_name, new_expression2) ==  TRUE ){paste("(",new_expression2,")[t]",sep="")}else{new_expression2}
  ##lowercase the expression
  new_expression4<-gsub("\\.t_data\\[,\\'e([0-9])\\'\\]","e\\1",new_expression3)
  final_expression<-paste(output_vector_name,"[",i,"]<-",new_expression4,sep="")
  plop_col<-c(plop_col,final_expression)

  }

  finalline <- paste("return <- ", output_vector_name, "BLERG",sep="")
  finallineb<- gsub("BLERG","\\}",finalline)
  plop_col<-c(plop_col,finallineb)
  return<-plop_col
}


#' Generate the command lines of function to evaluate the jacobian
#'
#' @param jacobian jacobian symbolic matrix
#' @param database_name database of the solver (t_data)
#' @param output_matrix_name output matrix name as used in solver
#' @param info.matrix where to source the matrix
#' @param currentstep will be used in the function name
#' @keywords internal
#'
#' @return character vector with the lines to build the function
get_jacobian_commands<-function(jacobian,database_name,output_matrix_name="output",info.matrix="jacobian",currentstep="prologue_jacobian"){

  matrixsize = ncol(jacobian)
  plop_col<-vector()
  ###Line 1 should read : currentstep_f<-function(){
  line1<- paste(currentstep,"_f<-" ,sep="")
  line1a<- gsub("<-","<-function\\(t=1,t_data,jacobian\\)\\{",line1)
  # cat(line1a,file=paste(name_of_file,".R",sep=""),sep="\n",append = TRUE)
  plop_col<-c(plop_col,line1a)
  ###Line 2 creates the zero matrix  should read :  Jacobian_n<-matrix(0,ncol = 201,nrow = 201)
  line2 <- paste(output_matrix_name,"<- matrixBLERG0,ncol = ",matrixsize,",nrow = ",matrixsize,"GRELB",sep = "")
  line2a <- gsub("BLERG","\\(",line2)
  line2b <- gsub("GRELB","\\)",line2a)
  # cat(line2b,file=paste(name_of_file,".R",sep=""),sep="\n",append = TRUE)
  plop_col<-c(plop_col,line2b)
  ###

  ###cases where jacobian is 1 or -1 are so frequent so lets group those
  final_expression<-paste(output_matrix_name,"which",info.matrix,"==1 <-1",sep="")
  fe1<-gsub("which","\\[which\\(",final_expression)
  fe2<-gsub("==1","==\\'1\\'\\)\\]",fe1)
  plop_col<-c(plop_col,fe2)
  # cat(fe2,file=paste(name_of_file,".R",sep=""),sep="\n",append = TRUE)

  final_expression<-paste(output_matrix_name,"which",info.matrix,"==-1 <- -1",sep="")
  fe1<-gsub("which","\\[which\\(",final_expression)
  fe2<-gsub("==-1","==\\'-1\\'\\)\\]",fe1)
  plop_col<-c(plop_col,fe2)
  # cat(fe2,file=paste(name_of_file,".R",sep=""),sep="\n",append = TRUE)

  ### write all other cases
  for ( i in 1:ncol(jacobian) )
  { for ( j in 1:ncol(jacobian) )
  { if (jacobian[i,j]%in% c("0","1","-1")==FALSE) {

    functions_supported<-paste0(thor_functions_supported,"\\(") %>% paste(collapse="|")
    function_patterns<- paste0("(",functions_supported,"|lag\\.|delta\\(|newdiff\\(|\\.e[0-9]\\()")

    old_expression<-tolower(jacobian[i,j])

    expression_functions<-gsub(function_patterns,"\\U\\1",old_expression,perl = TRUE)
    ##change delta to newdiff
    delta_expression<-gsub("DELTA\\(([0-9]+),","NEWDIFF\\(\\1,",expression_functions,perl = TRUE)
    ##deal with lags
    lag_expression<-gsub("LAG\\.(\\w+)\\.(\\w+)\\.(\\w+)","LAG\\(\\1,\\U\\2\\+\\U\\3\\)",delta_expression,perl = TRUE)
    lag_expression<-gsub("LAG\\.(\\w+)\\.(\\w+)","LAG\\(\\1,\\U\\2)",lag_expression,perl = TRUE)
    lag_expression_bis<-gsub("(LAG\\(\\w+,)(\\w+)\\+(\\w+)\\)",paste("\\1",toupper(database_name),"\\[T,\\'\\2\\'\\]\\+",toupper(database_name),"\\[T,\\'\\3\\'\\]\\)",sep=""),lag_expression,perl = TRUE)
    lag_expression_bis<-gsub("(LAG\\(\\w+,)(\\w+)\\)",paste("\\1",toupper(database_name),"\\[T,\\'\\2\\'\\]\\)",sep=""),lag_expression_bis,perl = TRUE)
    lag_expression_ter<-gsub(paste(toupper(database_name),"\\[T,\\'([0-9]+)'\\]",sep=""),"\\1",lag_expression_bis,perl = TRUE)
    ##replace all variables with database_name[,variable_name]
    new_expression1<-gsub("([a-z]([a-z0-9_]+)?)",paste(database_name,"\\[,\\'\\1\\'\\]",sep=""),lag_expression_ter)
    new_expression2<-tolower(new_expression1)
    ##add the [t] at the end and big parenthis for the whole expression
    new_expression3<- if(grepl(database_name, new_expression2) ==  TRUE ){paste("(",new_expression2,")[t]",sep="")}else{new_expression2}
    ##lowercase the expression
    new_expression4<-gsub("\\.t_data\\[,\\'e([0-9])\\'\\]","e\\1",new_expression3)
    final_expression<-paste(output_matrix_name,"[",i,",",j,"]<-",new_expression4,sep="")
    plop_col<-c(plop_col,final_expression)
    # cat( final_expression,file=paste(name_of_file,".R",sep=""),sep="\n",append = TRUE)
  }
  }
  }
  ##last line returns the jacobian and closes the functions
  finalline <- paste("return <- ", output_matrix_name, "BLERG",sep="")
  finallineb<- gsub("BLERG","\\}",finalline)
  # cat(finallineb,file=paste(name_of_file,".R",sep=""),sep="\n",append = TRUE)
  plop_col<-c(plop_col,finallineb)
  return<-plop_col
}
# prologue_function<-get_functions_commands(
#                 list_of_equations = equations_list[which(equations_list$part=="prologue"), "new_formula"],
#                 database_name ="t_data",
#                 output_vector_name="f_x_n",
#                 currentstep="prologue_equations")



# eval(parse(text=paste0(prologue_function,collapse = "\n")))
