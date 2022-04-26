
#' Quickly solve a static equation
#'
#' @param formula equation formula, must not have any symbolic variable other than the endogenous variable. Default: "5*x=210"
#' @param endogenous Variable to solve as it appears in the equation. Default: "x"
#' @param init numeric initial value for the solver to start the solving algorithm. Default : 1
#' @param quiet TRUE to not display the result in the console. Default : FALSE
#'
#' @return a numeric value, a possible solution of the equation
#' @export
#' @import assertthat
#' @import purrr
#' @import dplyr
#' @examples
#' quick_solve("log(x)+x=10","x")
#' quick_solve("3*x+4=10","x")
quick_solve <- function(formula="5*x=210", endogenous="x", init=1,quiet = FALSE){
  assertthat::assert_that(is.numeric(init))
  assertthat::assert_that(is.character(endogenous))

  res<-purrr::quietly(create_equation)("quick_solve_eq",
                  formula = formula,
                  endogenous = endogenous,env = environment())$result
  if(quick_solve_eq@maxlag != 0){stop("quick_solve only works for static equations.")}
  if(length(quick_solve_eq@exogenous) != 0){stop("quick_solve only works equations without exogenous variables.")}
  if(length(quick_solve_eq@coefflist) != 0){stop("quick_solve only works equations without coefficients.")}

  data <- data.frame(obs = c(0,1), ## obs servira d'indicateur de temps
                        x = c(init,NA))
  names(data)<-c("obs",endogenous)
  safe_solver <- purrr::safely(thor_equation_solver)
  solution <- purrr::quietly(safe_solver)(quick_solve_eq,first_period = 1,last_period = 1,
                          database = data , index_time = "obs")$result

  if(!is.null(solution$error)){
    cat("\nCouldn't solve the equation. Maybe try another init value. Or the equation cannot be solved. \nAdditionally : \n")
    cat(solution$error$message,"\n")
    stop("Couldn't solve the equation.")
    }

 res<- solution$result %>% filter(obs == 1) %>% select(paste0(quick_solve_eq@endogenous,".simul"))
 if(is.na(res[1,1])){cat("Found NA.")}
 message <- paste("\n",endogenous ,"=",res[1,1],"\n")
 if(quiet==FALSE){cat(message)}

 return <- res[1,1]

}
