## ---------------------------------------------------------------------------
## Solving ThreeME from its .txt model file
##
## A worked example of the sparse path: read a model from a text file, build
## it, solve it over a horizon, and check the answer.
##
##   Rscript tests/solve_threeme.R          # 4x4, the quicker one
##   Rscript tests/solve_threeme.R 8x8
##
## Run from the root of the package.
## ---------------------------------------------------------------------------

library(tresthor)

## ---- 1. What to run -------------------------------------------------------

classification <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(classification)) classification <- "4x4"

model_file <- file.path("tests", paste0("threeme_", classification, "_thor.txt"))
data_file  <- file.path("tests", paste0("data3me_", classification, ".rds"))

first_period <- 2016
last_period  <- 2050

## Where the generated C++ goes. Keep it somewhere stable rather than a temp
## directory: the model object refers to this file, so a saved model can be
## reloaded and re-solved later without rebuilding.
build_dir <- file.path("tests", "build")
dir.create(build_dir, showWarnings = FALSE, recursive = TRUE)


## ---- 2. Build the model from the .txt file --------------------------------
##
## The file has a fixed layout: endogenous variables on line 2, exogenous on
## line 5, coefficients on line 8, then one equation per line from line 11.
##
## create_model_sparse() parses it, splits the model into prologue / heart /
## epilogue, differentiates every equation symbolically, writes a C++ solver
## specialised to this model, and compiles it.

cat("\n== Building ThreeME", classification, "==\n")

model <- create_model_sparse(
  model_name   = paste0("threeme", classification),
  model_source = model_file,
  rcpp_path    = build_dir,
  no_var_map   = TRUE     # skip the variable map; set FALSE if you want var_map
)

## `model` is also assigned in the global environment under `model_name`.


## ---- 3. Load the data -----------------------------------------------------
##
## One row per period, one column per variable. Every exogenous variable and
## coefficient must be present over the whole solved horizon, and the period
## *before* first_period must be complete: it seeds the solver and supplies
## the lagged terms.

data_3me <- readRDS(data_file)

cat("\ndata: ", nrow(data_3me), " periods x ", ncol(data_3me), " variables",
    " (", min(data_3me$year), "-", max(data_3me$year), ")\n", sep = "")


## ---- 4. Solve -------------------------------------------------------------
##
## Returns the same data.frame with the endogenous variables filled in over
## [first_period, last_period]. diagnostics = TRUE attaches per-period
## iteration counts, residuals and convergence measures.

cat("\n== Solving", first_period, "to", last_period, "==\n")

result <- thor_solver_sparse(
  model                = model,
  first_period         = first_period,
  last_period          = last_period,
  database             = data_3me,
  index_time           = "year",
  convergence_criteria = 1e-10,   # relative: |dx| <= rtol*|x| + atol
  atol                 = 1e-8,    # absolute floor, for variables near zero
  max_iter             = 100,
  verbose              = TRUE,
  diagnostics          = TRUE
)


## ---- 5. Check the answer --------------------------------------------------
##
## Convergence is reported as a scaled step: <= 1 means every variable met the
## tolerance. Worth checking the residuals too -- that is the model's own
## equations evaluated at the solution, so it catches a "converged" answer
## that does not actually satisfy the model.

cat("\n== Diagnostics ==\n")
cat("Newton iterations per period (total ", sum(attr(result, "iterations")), "):\n", sep = "")
print(attr(result, "iterations"))

cat("\nworst scaled step over all periods:",
    format(max(attr(result, "convergence")), digits = 3), " (converged if <= 1)\n")

resid <- model_residuals(model, result,
                         periods = first_period:last_period, index_time = "year")
cat("max |residual| by block:\n")
print(apply(resid, 2, max))


## ---- 6. Look at the results -----------------------------------------------

vars <- intersect(c("gdp", "pib", "ch", "u", "co2", "ems_co2"), names(result))
if (length(vars)) {
  cat("\n== A few series ==\n")
  print(result[result$year %in% c(2016, 2020, 2030, 2040, 2050), c("year", vars)],
        row.names = FALSE)
}


## ---- 7. Save --------------------------------------------------------------

saveRDS(result, file.path(build_dir, paste0("threeme", classification, "_solved.rds")))
saveRDS(model,  file.path(build_dir, paste0("threeme", classification, "_model.rds")))

cat("\nSaved to", build_dir, "\n")

## To re-solve later without rebuilding, in a fresh session:
##
##   library(tresthor)
##   model <- readRDS("tests/build/threeme4x4_model.rds")
##   data  <- readRDS("tests/data3me_4x4.rds")
##   res   <- thor_solver_sparse(model, 2016, 2050, data, index_time = "year")
##
## Several models can be loaded at once; each keeps its own compiled code.
##
## ---- On compilation and caching -------------------------------------------
##
## Compiled models are cached on disk between sessions, so a given version of
## a model is compiled once per machine rather than once per session:
##
##   first ever build of ThreeME 4x4       ~31 s  (19 s symbolic + 14 s compile)
##   rebuilding it unchanged               ~20 s  (compile served from cache)
##   reloading and solving in a new session  0.5 s
##
## The cache invalidates itself when the equations change, and also on a new
## platform, compiler, R version or Rcpp version, so a stale object cannot be
## picked up. Editing the model rebuilds normally.
##
##   tresthor_cache_dir()          # where it lives
##   clear_model_cache()           # empty it
##   options(tresthor.cache.dir = "/path")   # move it, e.g. onto a shared disk
##   options(tresthor.cache.dir = FALSE)     # switch caching off
##
## Passing cache = FALSE to create_model_sparse() or thor_solver_sparse() does
## the same for a single call.
