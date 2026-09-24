## Sparse model builder.
##
## Same front end as `create_model()` (parsing, checks, block decomposition),
## but from the symbolic jacobian onwards nothing dense is ever built: the
## jacobian is kept as triplets and handed straight to the C++ generator.

#' @title Create a thoR.model with a sparse solver
#' @description Builds a model whose jacobians are stored and generated in
#'   sparse form. Intended for large models (thousands of equations), where the
#'   dense representation used by `create_model()` is not viable: on a model
#'   with n equations the dense path allocates an n x n matrix per Newton
#'   iteration and loops over n^2 cells at build time, while real macro models
#'   have well under 1% of those entries filled.
#'
#' @param model_name character. Name of the model object.
#' @param model_source path to the .txt model file.
#' @param algo boolean. TRUE to decompose the model into prologue/heart/epilogue.
#' @param rcpp_path directory where the generated C++ source is written.
#' @param endogenous,exogenous,coefficients,equations manual model input, used
#'   only when `model_source` is NULL.
#' @param env environment in which to create the model object.
#' @param no_var_map boolean. TRUE to skip building the variable map.
#' @param compile boolean. TRUE to compile the generated C++ immediately.
#' @param cache directory in which to cache the compiled object between R
#'   sessions, FALSE to disable, or NULL (the default) for
#'   `tresthor_cache_dir()`. Rebuilding a model whose equations have not
#'   changed then costs no compilation at all.
#'
#' @return A `thoR.model`, created in `env` and also returned invisibly.
#' @import Deriv
#' @import assertthat
#' @export
create_model_sparse <- function(model_name = "model",
                                model_source = NULL,
                                algo = TRUE,
                                rcpp_path = getwd(),
                                endogenous = NULL,
                                exogenous = NULL,
                                coefficients = NULL,
                                equations = NULL,
                                env = globalenv(),
                                no_var_map = TRUE,
                                compile = TRUE,
                                cache = NULL) {

  t_start <- Sys.time()
  cat("Building model '", model_name, "' (sparse)\n\n", sep = "")

  options(stringsAsFactors = FALSE)
  delta <- function(n, x, y = TRUE) NULL
  drule[["delta"]] <- alist(x = 1, y = NULL)
  ## quote() rather than alist(), which is equivalent here but keeps R's code
  ## checker from reading the rule's `x` as an undefined global.
  drule[["abs"]]   <- list(x = quote(ifelse(x == 0, 0, sign(x))))

  ################################
  #### 1. Read and check the model
  ################################
  cat("Step 1: reading and checking the model...\n")
  if (!is.null(model_source)) {
    assertthat::is.readable(model_source)
    check_model_file(model_source)
    model_input <- readLines(model_source, warn = FALSE)
    eqlist <- unique(model_input[11:length(model_input)])
    endo   <- tolower(unique(strsplit(model_input[2], split = ",")[[1]]))
    exo    <- tolower(unique(strsplit(model_input[5], split = ",")[[1]]))
    coeff  <- tolower(unique(strsplit(model_input[8], split = ",")[[1]]))
  } else {
    if (is.null(endogenous) || is.null(equations)) {
      stop("No model source provided and no equations specified.")
    }
    eqlist <- unique(equations)
    endo   <- tolower(unique(endogenous))
    exo    <- tolower(unique(exogenous))
    coeff  <- tolower(unique(coefficients))
  }

  endo  <- sort(endo);  exo <- sort(exo);  coeff <- sort(coeff)
  eqlist <- gsub("mylg\\(", "lag\\(", eqlist)

  check_equation_input(eqlist)
  check_var_vector(endo,  "Endogenous variables")
  check_var_vector(exo,   "Exogenous variables")
  check_var_vector(coeff, "Coefficient variables")
  check_variable_conflict(endo, exo)
  check_variable_conflict(endo, coeff)
  check_variable_conflict(exo, coeff)

  exo   <- is_in_formulas(exo,   eqlist, "exogenous")
  endo  <- is_in_formulas(endo,  eqlist, "endogenous")
  coeff <- is_in_formulas(coeff, eqlist, "coefficients")
  all_model_variables <- sort(c(endo, exo, coeff))

  equations_list <- create_equations_list(eqlist)
  if (parser_lag_delta_check(equations_list$equation) == FALSE) {
    stop("Parser error on the lags and/or delta. Please check the model's formulas.")
  }

  ################################
  #### 2. Endogenous variables per equation
  ################################
  cat("Step 2: identifying the endogenous variables in the equations...\n")
  eqns <- table_contemporaneous_endos(formula_list = equations_list$formula,
                                      endogenous = endo, exogenous = exo,
                                      coefflist = coeff,
                                      equations_index = equations_list$id)
  check_eq_var_identification(endo = endo, eqns = eqns,
                              names_of_equations = equations_list$name)

  ################################
  #### 3. Decomposition
  ################################
  cat("Step 3: decomposing the model into blocks...\n")
  decomposition <- decomposing_model(endogenous_variables = endo,
                                     eq_var_matrix = eqns,
                                     decomposition = algo)

  ## Unpacked explicitly rather than with list2env(), so that the block
  ## variables are visible to readers and to R's code checker.
  prologue <- decomposition$prologue
  heart    <- decomposition$heart
  epilogue <- decomposition$epilogue
  prologue_endo <- sort(decomposition$prologue_endo)
  heart_endo    <- sort(decomposition$heart_endo)
  epilogue_endo <- sort(decomposition$epilogue_endo)
  prologue_equations <- decomposition$prologue_equations
  heart_equations    <- decomposition$heart_equations
  epilogue_equations <- decomposition$epilogue_equations

  equations_list$part <- "tbd"
  equations_list$part[equations_list$id %in% prologue_equations] <- "prologue"
  equations_list$part[equations_list$id %in% heart_equations]    <- "heart"
  equations_list$part[equations_list$id %in% epilogue_equations] <- "epilogue"
  if ("tbd" %in% equations_list$part) {
    stop("Some equations were not identified in the decomposition.")
  }
  print(table(equations_list$part))

  equations_list$new_formula <- formatting_formulas(equations_list$formula)

  ################################
  #### 4. Sparse symbolic jacobians
  ################################
  cat("\nStep 4: computing the sparse symbolic jacobians...\n")
  eqns_as_list <- purrr::map(as.data.frame(t(eqns)), ~unique(stats::na.omit(.x)))

  build_block <- function(flag, block_endo, block_eqs, label) {
    if (!isTRUE(flag) || length(block_eqs) == 0L) return(NULL)
    j <- symbolic_jacobian_sparse(equations_list_df = equations_list,
                                  eqns_vars_list = eqns_as_list,
                                  endo_vec = block_endo,
                                  equations_subset = block_eqs)
    cat("   ", label, ": "); print(j)
    ## residual formulas, in the same order as the jacobian rows
    pos <- match(j$equations, as.character(equations_list$id))
    list(jac = j, formulas = equations_list$new_formula[pos])
  }

  blocks <- list(
    prologue = build_block(prologue, prologue_endo, prologue_equations, "prologue"),
    heart    = build_block(heart,    heart_endo,    heart_equations,    "heart"),
    epilogue = build_block(epilogue, epilogue_endo, epilogue_equations, "epilogue")
  )

  ################################
  #### 5. C++ generation
  ################################
  cat("\nStep 5: generating the sparse C++ solver...\n")
  rcpp_source <- create_model_rcpp_sparse(model_name, blocks,
                                          all_model_variables, rcpp_path)
  cat("   ", if (isTRUE(attr(rcpp_source, "unchanged"))) "unchanged: " else "written to ",
      rcpp_source, " (", round(file.info(rcpp_source)$size / 1024), " KB)\n", sep = "")

  ################################
  #### 6. Model object
  ################################
  var_map <- if (no_var_map) {
    list(no_var_map = TRUE,
         message = "To use var_map functionalities, rebuild with no_var_map = FALSE")
  } else {
    build_dico(equations_list_df = equations_list, endos = endo, exos = exo,
               coeffs = coeff, p_endos = prologue_endo, h_endos = heart_endo,
               e_endos = epilogue_endo)
  }

  none <- function() print("none")
  model <- thoR.model(
    model_name = model_name,
    equation_list = equations_list,
    var_map = var_map,
    variables_matrix = eqns,
    endo_list = tolower(endo), exo_list = tolower(exo), coeff_list = tolower(coeff),
    prologue = prologue, heart = heart, epilogue = epilogue,
    ## the dense slots stay empty: the sparse jacobians live in `sparse_jacobians`
    prologue_jacobian = matrix(), heart_jacobian = matrix(), epilogue_jacobian = matrix(),
    prologue_endo = prologue_endo, heart_endo = heart_endo, epilogue_endo = epilogue_endo,
    algo = algo, rcpp = TRUE, rcpp_source = rcpp_source,
    prologue_equations_f = none, heart_equations_f = none, epilogue_equations_f = none,
    prologue_jacobian_f = none, heart_jacobian_f = none, epilogue_jacobian_f = none)

  attr(model, "sparse_jacobians") <- lapply(blocks, function(b) b$jac)
  attr(model, "sparse") <- TRUE

  assign(model_name, model, envir = env)

  if (compile) {
    cat("\nStep 6: compiling...\n")
    el <- system.time(compile_model_cpp(rcpp_source, cache = cache))[["elapsed"]]
    cat("   compiled in ", round(el, 1), " s",
        if (el < 1) "  (from cache)" else "", "\n", sep = "")
  }

  cat("\nModel built in ",
      round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), " s\n", sep = "")
  invisible(get(model_name, envir = env))
}
