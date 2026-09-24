## Compilation of generated model sources.
##
## Generated model code is unlike hand-written code: a few thousand very large
## arithmetic expressions with no control flow. Two of R's default compiler
## settings behave pathologically on it.
##
## Measured on ThreeME 4x4 (1729 equations, 370 KB of generated C++):
##
##   R's default flags (-g -O2)      299 s
##   same without -g                  12 s
##
## `-g` makes clang emit debug metadata for every subexpression, which is
## worthless here -- nobody steps through generated code in a debugger -- and
## costs 25x the compile time. The effect grows with model size, so on a 25k
## equation model it is the difference between minutes and an hour.

#' Compile a generated model source file
#'
#' Compiles with the platform's usual flags minus `-g`, by way of a temporary
#' `R_MAKEVARS_USER` file. The user's own `~/.R/Makevars` is left untouched.
#'
#' @param path path to the generated .cpp file
#' @param debug boolean. TRUE to keep `-g` (much slower; only useful when
#'   debugging the code generator itself). Default FALSE.
#' @param quiet boolean. TRUE to suppress compiler output. Default TRUE.
#' @return the elapsed time in seconds, invisibly
#' @keywords internal
compile_model_cpp <- function(path, debug = FALSE, quiet = TRUE) {

  stopifnot(file.exists(path))

  if (debug) {
    t <- system.time(Rcpp::sourceCpp(path, rebuild = TRUE, verbose = !quiet))
    return(invisible(t[["elapsed"]]))
  }

  ## Start from the platform's own flags so we keep -arch, -falign-functions
  ## and anything else the build was configured with, and only drop -g.
  cxxflags <- tryCatch(system2("R", c("CMD", "config", "CXXFLAGS"),
                               stdout = TRUE, stderr = FALSE),
                       error = function(e) character(0))
  cxxflags <- paste(cxxflags, collapse = " ")
  if (!nzchar(trimws(cxxflags))) cxxflags <- "-O2"
  ## drop -g / -ggdb / -g3 ... but keep -g0 if it is already there
  flags <- strsplit(trimws(cxxflags), "\\s+")[[1]]
  flags <- flags[!grepl("^-g([0-9]|gdb.*)?$", flags)]
  cxxflags <- paste(c(flags, "-g0"), collapse = " ")

  mk <- tempfile(pattern = "tresthor_makevars_")
  writeLines(c(paste0("CXXFLAGS = ", cxxflags),
               paste0("CXX17FLAGS = ", cxxflags),
               paste0("CXX20FLAGS = ", cxxflags)), mk)

  old <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
  Sys.setenv(R_MAKEVARS_USER = mk)
  on.exit({
    if (is.na(old)) Sys.unsetenv("R_MAKEVARS_USER") else Sys.setenv(R_MAKEVARS_USER = old)
    unlink(mk)
  }, add = TRUE)

  t <- system.time(Rcpp::sourceCpp(path, rebuild = TRUE, verbose = !quiet))
  invisible(t[["elapsed"]])
}
