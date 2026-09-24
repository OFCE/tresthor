## ---------------------------------------------------------------------------
## Regression check: the original (dense) API still behaves as it did.
##
## The sparse solver added in 1.9.0 is a parallel code path -- no file under
## R/ that existed before was modified. This script confirms that in practice
## rather than by inspection.
##
##   Rscript tests/regression_legacy.R
##
## Run from the root of the package. Takes a few minutes: it builds ThreeME
## 4x4 and Opale through the legacy path.
## ---------------------------------------------------------------------------

library(tresthor)

OUT <- file.path(tempdir(), "tresthor_regression"); dir.create(OUT, showWarnings = FALSE)
failures <- character(0)

ok <- function(label, expr) {
  r <- try(force(expr), silent = TRUE)
  bad <- inherits(r, "try-error")
  if (bad) failures <<- c(failures, label)
  cat(sprintf("  %-46s %s\n", label,
              if (bad) paste("FAIL:", sub("\n.*", "", attr(r, "condition")$message)) else "ok"))
  invisible(r)
}

cat("tresthor", as.character(packageVersion("tresthor")), "\n")

## ---- 1. model building and the utility functions --------------------------

cat("\n== model building ==\n")
src <- system.file("models", "model_type.txt", package = "tresthor")

ok("create_model(rcpp = FALSE)",
   create_model("base_m", model_source = src, rcpp = FALSE,
                no_var_map = FALSE, env = globalenv()))
ok("create_model(rcpp = TRUE)",
   create_model("base_c", model_source = src, rcpp = TRUE, rcpp_path = OUT,
                no_var_map = TRUE, env = globalenv()))

cat("\n== model modification ==\n")
ok("create_equation()",
   create_equation("newone", formula = "endovar5 = endovar1 + exovar1",
                   endogenous = "endovar5", coefflist = NULL, env = globalenv()))
ok("model_equations_add()",
   model_equations_add(base_model = base_m, new_model_name = "m_add",
                       thor_equations_add = list(newone),
                       rcpp = FALSE, algo = TRUE, env = globalenv()))
ok("model_equations_remove()",
   model_equations_remove(base_model = base_m, new_model_name = "m_rm",
                          equations_to_remove = "equation_var3",
                          endos_to_remove = "endovar3",
                          rcpp = FALSE, algo = TRUE, env = globalenv()))
ok("model_endo_exo_switch()",
   model_endo_exo_switch(base_model = base_m, new_model_name = "m_sw",
                         new_endo = "exovar1", new_exo = "endovar4",
                         rcpp = FALSE, algo = TRUE, env = globalenv()))

cat("\n== import / export ==\n")
ok("export_model()", export_model(base_m, file.path(OUT, "base_m.txt")))
ok("re-import the exported model",
   create_model("reimported", model_source = file.path(OUT, "base_m.txt"),
                rcpp = FALSE, env = globalenv()))
ok("round trip preserves the model",
   stopifnot(identical(sort(base_m@endo_list), sort(reimported@endo_list)),
             nrow(base_m@equation_list) == nrow(reimported@equation_list)))
ok("save_model() + load_model()",
   { save_model(base_m, OUT); rm(base_m, envir = globalenv())
     load_model(model = "base_m", folder_path = OUT)
     stopifnot(inherits(get("base_m", globalenv()), "thoR.model")); TRUE })

cat("\n== utilities ==\n")
ok("var_info_model()", var_info_model("endovar1", get("base_m", globalenv())))
ok("quick_solve()",    quick_solve(formula = "5*x=210", endogenous = "x",
                                   init = 1, quiet = TRUE))
ok("formula_latex()",  formula_latex("delta(1,log(y))=a*x"))
ok("delta() / newdiff()",
   stopifnot(all.equal(delta(1, c(1, 3, 6))[2:3], c(2, 3))))

## ---- 2. Opale, end to end, on both solver backends ------------------------

cat("\n== Opale, both solver backends ==\n")
opale_src <- system.file("models", "opale.txt", package = "tresthor")
dat <- readRDS(system.file("Opale", "donnees_opale.rds", package = "tresthor"))
cf  <- readRDS(system.file("Opale", "coefficients_opale.rds", package = "tresthor"))

ok("create_model(opale, rcpp = TRUE)",
   create_model("opale_c", model_source = opale_src, rcpp = TRUE, rcpp_path = OUT,
                no_var_map = TRUE, env = globalenv()))
ok("create_model(opale, rcpp = FALSE)",
   create_model("opale_r", model_source = opale_src, rcpp = FALSE,
                no_var_map = TRUE, env = globalenv()))
ok("add_coeffs()",
   assign("dat2", add_coeffs(listcoeff = cf, database = dat,
                             pos.coeff.name = "name", pos.coeff.value = "value"),
          envir = globalenv()))

t0    <- which(dat2$date == as.Date("2015-01-01"))
first <- dat2$date[t0]; last <- dat2$date[t0 + 7]

ok("thor_solver(rcpp = TRUE)",
   assign("r_c", thor_solver(opale_c, first_period = first, last_period = last,
                             database = dat2, index_time = "date",
                             rcpp = TRUE, skip_tests = TRUE), envir = globalenv()))
ok("thor_solver(rcpp = FALSE)",
   assign("r_r", thor_solver(opale_r, first_period = first, last_period = last,
                             database = dat2, index_time = "date",
                             rcpp = FALSE, skip_tests = TRUE), envir = globalenv()))

rows <- t0:(t0 + 7)
A <- as.matrix(r_c[rows, opale_c@endo_list])
B <- as.matrix(r_r[rows, opale_c@endo_list])
ok("both backends produce a complete solution",
   stopifnot(!anyNA(A), !anyNA(B)))
cat(sprintf("  %-46s %.2e\n", "max relative difference between backends",
            max(abs(A - B) / pmax(abs(A), 1e-8), na.rm = TRUE)))

## ---- 3. ThreeME 4x4 through the legacy path -------------------------------
##
## The heaviest exercise of the old code: symbolic jacobian, Rcpp generation
## and the dense solver, on a 1729 equation model. Skipped unless the test
## data is present (it lives in tests/, not in the installed package).

if (file.exists("tests/threeme_4x4_thor.txt")) {
  cat("\n== ThreeME 4x4, legacy dense path ==\n")
  ok("create_model(threeme 4x4, rcpp = TRUE)",
     create_model("t4_legacy", model_source = "tests/threeme_4x4_thor.txt",
                  rcpp = TRUE, rcpp_path = OUT, no_var_map = TRUE, env = globalenv()))
  d4 <- readRDS("tests/data3me_4x4.rds")
  ok("thor_solver(threeme 4x4, rcpp = TRUE)",
     assign("t4res", thor_solver(model = t4_legacy, first_period = 2016,
                                 last_period = 2050, database = d4,
                                 index_time = "year", rcpp = TRUE,
                                 skip_tests = TRUE), envir = globalenv()))
  ok("solution is complete",
     stopifnot(!anyNA(as.matrix(t4res[t4res$year >= 2016 & t4res$year <= 2050,
                                      t4_legacy@endo_list]))))
} else {
  cat("\n(skipping ThreeME: tests/threeme_4x4_thor.txt not found)\n")
}

## ---- verdict --------------------------------------------------------------

cat("\n", strrep("-", 60), "\n", sep = "")
if (length(failures)) {
  cat("FAILED:\n"); cat(paste0("  - ", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("All legacy checks passed.\n")
