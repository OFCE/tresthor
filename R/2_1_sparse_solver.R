## Solver front end for models built with `create_model_sparse()`.
##
## The R side does nothing per-iteration: it prepares the data matrix once,
## hands it to the generated C++ and copies the endogenous columns back.

#' @title Solve a sparse thoR model
#'
#' @param model a `thoR.model` built by `create_model_sparse()`
#' @param first_period first period to solve, as found in `index_time`
#' @param last_period last period to solve
#' @param database data.frame holding the data
#' @param index_time name of the time column in `database`
#' @param convergence_criteria relative step tolerance (`rtol`). A block has
#'   converged when every endogenous variable satisfies
#'   `|dx| <= rtol*|x| + atol`. Relative rather than absolute because macro
#'   models mix variables of very different magnitude. Default 1e-10.
#' @param atol absolute step tolerance, which governs variables sitting at or
#'   near zero. Default 1e-8.
#' @param max_iter maximum Newton iterations per block and period. Default 100.
#' @param damping boolean. Backtrack the Newton step when it would increase the
#'   residual. Default TRUE.
#' @param verbose boolean. Print progress per period. Default TRUE.
#' @param diagnostics boolean. Return iteration counts and residuals as
#'   attributes of the result. Default FALSE.
#'
#' @return the input data.frame with the endogenous variables solved.
#' @export
thor_solver_sparse <- function(model,
                               first_period,
                               last_period,
                               database,
                               index_time = "date",
                               convergence_criteria = 1e-10,
                               atol = 1e-8,
                               max_iter = 100L,
                               damping = TRUE,
                               verbose = TRUE,
                               diagnostics = FALSE) {

  stopifnot(convergence_criteria < 0.01, atol >= 0, max_iter >= 1L)

  time_vec <- database[[index_time]]
  numeric_index_time <- is.numeric(time_vec)
  key <- as.character(time_vec)

  anchor_t <- match(as.character(first_period), key)
  final_t  <- match(as.character(last_period),  key)
  if (is.na(anchor_t)) stop("The first period is not found in '", index_time, "'.")
  if (is.na(final_t))  stop("The last period is not found in '", index_time, "'.")
  if (final_t < anchor_t) stop("The last period comes before the first period.")
  if (anchor_t == 1L) {
    stop("The first period cannot be the first observation: the solver needs a previous full observation to initialise.")
  }

  ## Columns of the data matrix are the model variables in alphabetical order.
  ## The generated C++ indexes into exactly this ordering.
  model_variables <- sort(c(model@exo_list, model@endo_list, model@coeff_list))
  missing <- setdiff(model_variables, names(database))
  if (length(missing)) {
    stop("Missing variables in the database: ",
         paste(utils::head(missing, 10), collapse = ", "),
         if (length(missing) > 10) sprintf(" (and %d more)", length(missing) - 10) else "")
  }

  M <- as.matrix(database[, model_variables, drop = FALSE])
  storage.mode(M) <- "double"

  ## The period before the first solved one has to be complete: it seeds the
  ## Newton start value and supplies every lagged term.
  na_prev <- model_variables[is.na(M[anchor_t - 1L, ])]
  if (length(na_prev)) {
    stop("The observation before the first period has missing values, so the solver cannot initialise: ",
         paste(utils::head(na_prev, 10), collapse = ", "),
         if (length(na_prev) > 10) sprintf(" (and %d more)", length(na_prev) - 10) else "")
  }

  ## Exogenous variables and coefficients must be present over the whole span.
  fixed <- c(model@exo_list, model@coeff_list)
  if (length(fixed)) {
    span <- M[anchor_t:final_t, fixed, drop = FALSE]
    bad <- fixed[apply(is.na(span), 2L, any)]
    if (length(bad)) {
      stop("Exogenous variables or coefficients are missing over the solved period: ",
           paste(utils::head(bad, 10), collapse = ", "),
           if (length(bad) > 10) sprintf(" (and %d more)", length(bad) - 10) else "")
    }
  }

  if (!exists("sparse_solver", mode = "function")) {
    compile_model_cpp(model@rcpp_source)
  }

  t_run <- system.time({
    out <- sparse_solver(M,
                         as.integer(anchor_t - 1L),   # 0-based rows for C++
                         as.integer(final_t - 1L),
                         convergence_criteria,
                         atol,
                         as.integer(max_iter),
                         isTRUE(damping),
                         isTRUE(verbose))
  })

  solved <- out$data
  colnames(solved) <- model_variables
  database[, model@endo_list] <- solved[, model@endo_list, drop = FALSE]

  if (verbose) {
    cat("Solved ", final_t - anchor_t + 1L, " periods in ",
        round(t_run[["elapsed"]], 2), " s (",
        sum(out$iterations), " Newton iterations, worst scaled step ",
        format(max(out$convergence), digits = 3), " (converges at 1), max |residual| ",
        format(max(out$residuals), digits = 3), ")\n", sep = "")
  }

  if (numeric_index_time) database[[index_time]] <- as.numeric(database[[index_time]])

  if (diagnostics) {
    attr(database, "iterations") <- stats::setNames(out$iterations, key[anchor_t:final_t])
    attr(database, "residuals")  <- stats::setNames(out$residuals,  key[anchor_t:final_t])
    attr(database, "convergence") <- stats::setNames(out$convergence, key[anchor_t:final_t])
    attr(database, "elapsed")    <- t_run[["elapsed"]]
  }
  database
}
