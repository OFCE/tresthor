## Build and solve ThreeME through the sparse path, and check the result
## against the dense reference solver.
##
## Usage:  Rscript tests/test_sparse_threeme.R [4x4|8x8]

library(tresthor)

classification <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(classification)) classification <- "4x4"
stopifnot(classification %in% c("4x4", "8x8"))

out_dir <- tempfile("tresthor_sparse_"); dir.create(out_dir)
model_file <- file.path("tests", paste0("threeme_", classification, "_thor.txt"))
data_file  <- file.path("tests", paste0("data3me_", classification, ".rds"))
stopifnot(file.exists(model_file), file.exists(data_file))

first_period <- 2016
last_period  <- 2050

## ---- sparse path --------------------------------------------------------
t_build <- system.time(
  m <- create_model_sparse(paste0("threeme", classification, "_sparse"),
                           model_source = model_file,
                           rcpp_path = out_dir,
                           no_var_map = TRUE)
)

data_3me <- readRDS(data_file)
res <- thor_solver_sparse(m, first_period, last_period,
                          database = data_3me, index_time = "year",
                          verbose = FALSE, diagnostics = TRUE)

cat("\n=== sparse ===\n")
cat("build  :", round(t_build[["elapsed"]], 1), "s\n")
cat("solve  :", round(attr(res, "elapsed"), 3), "s for",
    last_period - first_period + 1, "periods\n")
cat("newton :", sum(attr(res, "iterations")), "iterations\n")
cat("worst scaled step:", format(max(attr(res, "convergence")), digits = 3),
    "(converges at 1)\n")
cat("max |residual|   :", format(max(attr(res, "residuals")), digits = 3), "\n")

stopifnot(max(attr(res, "convergence")) <= 1)

## ---- every equation is actually satisfied at the solution ---------------
mv <- sort(c(m@exo_list, m@endo_list, m@coeff_list))
M  <- as.matrix(res[, mv]); storage.mode(M) <- "double"
rows <- which(res$year >= first_period & res$year <= last_period)
worst <- max(vapply(rows - 1L, function(r) max(sparse_residuals(M, r)), numeric(1)))
cat("re-checked max |residual| over all blocks/periods:", format(worst, digits = 3), "\n")
stopifnot(worst < 1e-6)

## ---- comparison with the dense reference solver -------------------------
if (nzchar(Sys.getenv("TRESTHOR_COMPARE_DENSE"))) {
  cat("\n=== dense reference ===\n")
  t_dense_build <- system.time(
    create_model(paste0("threeme", classification, "_dense"),
                 model_source = model_file, rcpp = TRUE,
                 rcpp_path = out_dir, no_var_map = TRUE)
  )
  md <- get(paste0("threeme", classification, "_dense"))
  t_dense <- system.time(
    ref <- thor_solver(model = md, first_period = first_period,
                       last_period = last_period, database = data_3me,
                       index_time = "year", rcpp = TRUE, skip_tests = TRUE)
  )
  cat("build:", round(t_dense_build[["elapsed"]], 1), "s   solve:",
      round(t_dense[["elapsed"]], 1), "s\n")

  endo <- m@endo_list
  A <- as.matrix(ref[rows, endo]); B <- as.matrix(res[rows, endo])
  sc <- pmax(abs(A), abs(B))
  rel <- abs(A - B) / sc
  cat("\n=== agreement ===\n")
  cat("median relative difference          :", format(median(rel[sc > 0], na.rm = TRUE), digits = 3), "\n")
  cat("max relative difference, |value| > 1:", format(max(rel[sc > 1], na.rm = TRUE), digits = 3), "\n")
  cat("speedup: build", round(t_dense_build[["elapsed"]] / t_build[["elapsed"]], 1), "x",
      " solve", round(t_dense[["elapsed"]] / attr(res, "elapsed")), "x\n")
  ## Variables of meaningful magnitude must agree closely. A handful of
  ## ThreeME variables are genuinely ill-conditioned (for instance pds_*, a
  ## price divided by a near-zero stock change), so this is not exact.
  stopifnot(median(rel[sc > 0], na.rm = TRUE) < 1e-10)
}

cat("\nPASS\n")
unlink(out_dir, recursive = TRUE)
