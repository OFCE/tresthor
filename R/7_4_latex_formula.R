

#' Transform a formula into LaTex code for display
#'
#' @param equation a formula as a character string or a thoR.equation
#' @param coeff_list a vector of coeeficients that could be in the equation. If equation is a thoR.equation object, this argument is skipped. If NULL coefficients that are still presnet in the formula will not be formatted.
#'
#' @return a character string written into LaTex to use with Rmarkdown or LaTex
#' @export
#' @import assertthat
#' @examples
#' my_eq <- "delta(1,log(var_1)) = cst + cf_1 *delta(1,log(var_2)) - 0.5 * lag(var_1,1)"
#' cat(formula_latex(my_eq,c("cst","cf1")))
formula_latex<-function(equation,coeff_list = NULL){
  if(inherits(equation,'thoR.equation',TRUE)){
    coeffs <- equation@coefflist
    formula <- equation@formula
  }else{
    assertthat::is.string(equation)
    if(!is.null(coeff_list)){coeffs<-coeff_list
    }else{
      coeffs <- c()
    }
    formula <- equation
  }

  formula<-tolower(formula)
  variables <- setdiff(get_variables_from_string(formula),coeffs)

  str_var <- paste0( "((?<![a-z0-9_])",paste(variables,collapse = "(?![a-z0-9_])|(?<![a-z0-9_])"),"(?![a-z0-9_]))")
  str_coeff <- paste0( "((?<![a-z0-9_])",paste(coeffs,collapse = "(?![a-z0-9_])|(?<![a-z0-9_])"),"(?![a-z0-9_]))")
  str_funs <- paste0("(",paste(thor_functions_supported,collapse = "\\(|"),"\\()")


  ### to do : replace variable names by pure text bold
  ### replace coeffs by pure text
  formula <-gsub("\\+\\+","+" ,formula,perl = TRUE )
  formula <-gsub("\\+-","-" ,formula,perl = TRUE )
  formula <-gsub("-\\+","-" ,formula,perl = TRUE )
  formula <-gsub("--","+" ,formula,perl = TRUE )
  formula <-gsub("\\\n","\\\\\\\\" ,formula,perl = TRUE )
  formula<-gsub("\\s+","",formula,perl= TRUE)
  formula <-gsub("(lag|mylg)\\((\\w+),-?(.)\\)"," \\\\textbf{\\U\\2}_{t-\\3} " ,formula,perl = TRUE )

  if(length(variables)>0){
    formula <-gsub(str_var," \\\\textbf{\\U\\1}_{t} " ,formula,perl = TRUE )
  }
  if(length(coeffs)>0){
    formula <-gsub(str_coeff," \\\\textit{\\U\\1} " ,formula,perl = TRUE )
  }
  formula <-gsub(str_funs," \\\\mathit{\\U\\1}" ,formula,perl = TRUE )
  formula <-tolower(formula)

  formula <-gsub("delta\\(([0-9]+),"," \\\\Delta_{\\1}(" ,formula,perl = TRUE )
  formula <-gsub("^\\s+","" ,formula,perl = TRUE )
  formula <-gsub("\\s+$","" ,formula,perl = TRUE )

}
