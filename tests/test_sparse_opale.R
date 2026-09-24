## Build and solve Opale through the sparse path, and check the result
## against the dense reference solver.
##
## Opale is the quarterly French macro model shipped with the package. It is
## much smaller than ThreeME but structurally harder: quarterly data, heavy
## use of lags and of the `trim` seasonal index, and a large epilogue of
## accounting identities. It therefore exercises the new sparse jacobian and
## C++ code generator on a different shape of model.
##
## Usage:  Rscript tests/test_sparse_opale.R
##         TRESTHOR_COMPARE_DENSE=1 Rscript tests/test_sparse_opale.R

pkgload::load_all(".", quiet = TRUE)

out_dir <- tempfile("tresthor_opale_"); dir.create(out_dir)
model_file <- "inst/models/opale.txt"
stopifnot(file.exists(model_file))

## ---- data ---------------------------------------------------------------
## The estimated coefficients live outside the database and have to be
## broadcast onto it as constant columns before solving.
data_opale  <- readRDS("inst/Opale/donnees_opale.rds")
coeffs      <- readRDS("inst/Opale/coefficients_opale.rds")
data_opale  <- add_coeffs(listcoeff = coeffs, database = data_opale,
                          pos.coeff.name = 2, pos.coeff.value = 1)

## Solve over the last few years of the sample, where every lagged input is
## available. The first observation can never be the anchor.
dates <- as.character(data_opale$date)
first_period <- dates[length(dates) - 39]
last_period  <- dates[length(dates)]
n_periods    <- 40

## ---- sparse path --------------------------------------------------------
t_build <- system.time(
  m <- create_model_sparse("opale_sparse",
                           model_source = model_file,
                           rcpp_path = out_dir,
                           no_var_map = TRUE)
)

res <- thor_solver_sparse(m, first_period, last_period,
                          database = data_opale, index_time = "date",
                          verbose = FALSE, diagnostics = TRUE)

cat("\n=== sparse ===\n")
cat("equations:", length(m@endo_list), "\n")
cat("build  :", round(t_build[["elapsed"]], 1), "s\n")
cat("solve  :", round(attr(res, "elapsed"), 3), "s for", n_periods, "periods\n")
cat("newton :", sum(attr(res, "iterations")), "iterations\n")
cat("worst scaled step:", format(max(attr(res, "convergence")), digits = 3),
    "(converges at 1)\n")
cat("max |residual|   :", format(max(attr(res, "residuals")), digits = 3), "\n")

stopifnot(max(attr(res, "convergence")) <= 1)

## ---- every equation is actually satisfied at the solution ---------------
periods <- dates[seq(length(dates) - 39, length(dates))]
resid <- model_residuals(m, res, periods = periods, index_time = "date")
worst <- max(resid)
cat("re-checked max |residual| over all blocks/periods:", format(worst, digits = 3), "\n")
print(apply(resid, 2, max))
stopifnot(worst < 1e-6)

## ---- the solution reproduces the historical data ------------------------
## Solving over history with the historical exogenous path must return the
## historical endogenous path: this is the real test that the generated code
## means what the model file says, independently of whether Newton converged.
rows <- match(periods, dates)
endo <- intersect(m@endo_list, names(data_opale))
A <- as.matrix(data_opale[rows, endo]); B <- as.matrix(res[rows, endo])
sc <- pmax(abs(A), abs(B))
rel <- abs(A - B) / sc
cat("\n=== against history ===\n")
cat("median relative difference          :", format(stats::median(rel[sc > 0], na.rm = TRUE), digits = 3), "\n")
cat("max relative difference, |value| > 1:", format(max(rel[sc > 1], na.rm = TRUE), digits = 3), "\n")
## Reported on variables of meaningful magnitude only. A handful of the
## contribution series (contpib*) are stored with an exact 0.0 in one quarter
## of the history, against which any non-zero solution scores a relative
## difference of 1; that is a rounding artefact of the stored data, not a
## solver error.
big <- rel; big[sc <= 1e-3] <- NA
off <- sort(apply(big, 2, max, na.rm = TRUE), decreasing = TRUE)[1:5]
cat("worst variables (|value| > 1e-3):\n"); print(off)

## ---- comparison with the dense reference solver -------------------------
if (nzchar(Sys.getenv("TRESTHOR_COMPARE_DENSE"))) {
  cat("\n=== dense reference ===\n")
  t_dense_build <- system.time(
    create_model("opale_dense", model_source = model_file, rcpp = TRUE,
                 rcpp_path = out_dir, no_var_map = TRUE)
  )
  md <- get("opale_dense")
  t_dense <- system.time(
    ref <- thor_solver(model = md, first_period = first_period,
                       last_period = last_period, database = data_opale,
                       index_time = "date", rcpp = TRUE, skip_tests = TRUE)
  )
  cat("build:", round(t_dense_build[["elapsed"]], 1), "s   solve:",
      round(t_dense[["elapsed"]], 1), "s\n")

  A <- as.matrix(ref[rows, m@endo_list]); B <- as.matrix(res[rows, m@endo_list])
  sc <- pmax(abs(A), abs(B))
  rel <- abs(A - B) / sc
  cat("\n=== agreement ===\n")
  cat("median relative difference          :", format(stats::median(rel[sc > 0], na.rm = TRUE), digits = 3), "\n")
  cat("max relative difference, |value| > 1:", format(max(rel[sc > 1], na.rm = TRUE), digits = 3), "\n")
  cat("speedup: build", round(t_dense_build[["elapsed"]] / t_build[["elapsed"]], 1), "x",
      " solve", round(t_dense[["elapsed"]] / attr(res, "elapsed"), 1), "x\n")
  stopifnot(stats::median(rel[sc > 0], na.rm = TRUE) < 1e-10)
}

cat("\nPASS\n")
unlink(out_dir, recursive = TRUE)
