##computing symbolic jacobian


#' Computing the symbolic jacobian
#'
#' @param equations_list_df dataframe with the equation infos and formulas.
#' @param eqns_vars_list equations variable matrix as a list.
#' @param endo_vec vector of the endogenous variables by which to compute the partial derivative.
#' @param equations_subset vector with the ids of the equations.
#' @param id_col column in equations_list_df that contains the id. Default : "id"
#' @param formula_col column in equations_list_df that contains the formula. Default : "new_formula"
#'
#' @importFrom purrr map
#' @import Deriv
#' @return character, named matrix of jacobians
#' @keywords internal
symbolic_jacobian<- function(equations_list_df, eqns_vars_list, endo_vec, equations_subset, id_col = "id",formula_col = "new_formula"){
   delta <- function(n,x, y=TRUE) NULL
   drule[["delta"]] <- alist(x=1, y=NULL) # y is just a logical
  if(length(equations_subset) != length(endo_vec)){stop(paste0("The length of the endogenous vector (",endo_vec, ") does not match the number of equations."))}

   eqns_subset <- equations_subset[order(equations_subset)]
   equations_list_df$id_bis <- equations_list_df[,id_col]
   formulas <- equations_list_df[which(equations_list_df$id_bis %in% equations_subset)  ,c("id_bis",formula_col)]
   eqns_var_info <- eqns_vars_list[which(names(eqns_vars_list) %in% equations_subset )]

   if(nrow(formulas) != length(equations_subset)){stop("Some equations of the subset were not found in the equation_list data.frame.")}
   if(length(eqns_var_info) != length(equations_subset)){stop("Some equations of the subset were not found in the eqns_var list.")}

   if (identical(formulas$id_bis, names(eqns_var_info)) == FALSE ){
     cat("Warning : the order of equations ids in the formula table do not match the order of the names of the equations in thelist. They will be reordered. \n")
     formulas <- formulas[order(formulas$id_bis,decreasing = FALSE),]
     eqns_var_info <- eqns_var_info[order(names(eqns_var_info),decreasing = FALSE)]
   }

   jacobian<-matrix(0,nrow = length(eqns_var_info),ncol = length(endo_vec))

    rownames(jacobian)<-eqns_subset
    colnames(jacobian)<-endo_vec
    jacobian<-cbind(jacobian, rowname=rownames(jacobian) )

  partial_derivatives_of_equation<-function(equation_name,jacobian){
    my_f <- parse(text=formulas$new_formula[which(formulas$id_bis ==equation_name)] )
    my_varlist <- intersect(as.vector(eqns_var_info[[equation_name]]),endo_vec )
    line_in_jacobian<- jacobian[equation_name,]
    for ( j in 1:length(my_varlist)){
      h<-my_varlist[j]
      my_derivative <-Deriv(my_f,h)
      line_in_jacobian[h]<-paste(my_derivative)
    }
    return(line_in_jacobian)
  }


   new_jacobian<-purrr::map_df(equations_subset,~partial_derivatives_of_equation(.x,jacobian))
   new_jacobian<-as.data.frame(new_jacobian) ##removing tibble to be able to give names

   new_jacobian<-new_jacobian[order(new_jacobian$rowname),]
   rownames(new_jacobian)<-new_jacobian$rowname
   new_jacobian$rowname<- NULL

   jacobian_output<-as.matrix(new_jacobian)

   jacobian_output
}


