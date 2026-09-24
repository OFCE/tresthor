## Compilation and loading of generated model sources.
##
## Generated model code is unlike hand-written code: a few thousand very large
## arithmetic expressions with no control flow. One of R's default compiler
## settings behaves pathologically on it.
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

## Compiled models, keyed by the absolute path of their source file. Every
## generated file exports functions with the same names (`sparse_solver`,
## `sparse_residuals`), so each model is loaded into its own environment
## rather than into the global one: otherwise a second model in the same
## session would silently be solved with the first one's code.
.tresthor_compiled <- new.env(parent = emptyenv())

#' Directory used to cache compiled models between sessions
#'
#' Without a persistent cache, `Rcpp::sourceCpp()` builds into `tempdir()` and
#' so recompiles every model once per R session -- around 14 s for ThreeME
#' 4x4 and 30 s for 8x8, and considerably more for larger models. With one,
#' the compilation happens once per machine and per version of the model.
#'
#' Rcpp keys its cache on the hash of the source, the platform and its own
#' version, so an edited model or a changed toolchain rebuilds by itself. It
#' does not key on the R version, and R binaries are not compatible across
#' minor releases, so that is added here.
#'
#' Set `options(tresthor.cache.dir = "...")` to move it, or
#' `options(tresthor.cache.dir = FALSE)` to disable caching entirely.
#'
#' @return the cache directory, created if needed, or NULL if caching is off
#' @export
tresthor_cache_dir <- function() {

  d <- getOption("tresthor.cache.dir", NULL)
  if (isFALSE(d)) return(NULL)

  if (is.null(d)) {
    d <- tryCatch(tools::R_user_dir("tresthor", "cache"),
                  error = function(e) file.path(tempdir(), "tresthor-cache"))
  }

  d <- file.path(path.expand(d), paste0("R-", getRversion()))
  if (!dir.exists(d)) {
    ok <- dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(d)) {
      warning("Could not create the tresthor cache directory '", d,
              "'. Models will be recompiled each session.", call. = FALSE)
      return(NULL)
    }
  }
  d
}

#' Remove cached compiled models
#'
#' @param confirm boolean. FALSE to skip the interactive confirmation.
#' @return the number of files removed, invisibly
#' @export
clear_model_cache <- function(confirm = interactive()) {
  d <- tresthor_cache_dir()
  if (is.null(d) || !dir.exists(d)) {
    cat("Nothing cached.\n")
    return(invisible(0L))
  }
  files <- list.files(d, recursive = TRUE, full.names = TRUE)
  size <- sum(file.info(files)$size, na.rm = TRUE)
  cat("Cache: ", d, "\n", length(files), " files, ",
      round(size / 1024^2, 1), " MB\n", sep = "")
  if (confirm && !isTRUE(utils::askYesNo("Delete?"))) return(invisible(0L))
  unlink(d, recursive = TRUE)
  ## also drop anything loaded in this session, so the next solve rebuilds
  rm(list = ls(.tresthor_compiled, all.names = TRUE), envir = .tresthor_compiled)
  invisible(length(files))
}

#' Compile a generated model source file and return its functions
#'
#' Compiles with the platform's usual flags minus `-g`, by way of a temporary
#' `R_MAKEVARS_USER` file. The user's own `~/.R/Makevars` is left untouched.
#'
#' The result is cached per source file, so solving the same model repeatedly
#' compiles once per session.
#'
#' @param path path to the generated .cpp file
#' @param rebuild boolean. TRUE to compile from scratch, ignoring both the
#'   in-session and the on-disk cache.
#' @param cache directory in which to cache the compiled object between
#'   sessions, FALSE to disable, or NULL (the default) for
#'   `tresthor_cache_dir()`.
#' @param debug boolean. TRUE to keep `-g` (much slower; only useful when
#'   debugging the code generator itself). Default FALSE.
#' @param quiet boolean. TRUE to suppress compiler output. Default TRUE.
#' @return an environment holding the model's compiled functions
#' @keywords internal
compile_model_cpp <- function(path, rebuild = FALSE, cache = NULL,
                              debug = FALSE, quiet = TRUE) {

  stopifnot(file.exists(path))
  key <- normalizePath(path)

  if (!rebuild && !is.null(.tresthor_compiled[[key]])) {
    return(.tresthor_compiled[[key]])
  }

  env <- new.env(parent = globalenv())

  if (!debug) {
    ## Start from the platform's own flags so we keep -arch, -falign-functions
    ## and anything else the build was configured with, and only drop -g.
    cxxflags <- tryCatch(system2("R", c("CMD", "config", "CXXFLAGS"),
                                 stdout = TRUE, stderr = FALSE),
                         error = function(e) character(0))
    cxxflags <- paste(cxxflags, collapse = " ")
    if (!nzchar(trimws(cxxflags))) cxxflags <- "-O2"
    flags <- strsplit(trimws(cxxflags), "\\s+")[[1]]
    flags <- flags[!grepl("^-g([0-9]|gdb.*)?$", flags)]   # keep an explicit -g0
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
  }

  cache_dir <- if (is.null(cache)) tresthor_cache_dir() else if (isFALSE(cache)) NULL else cache
  if (is.null(cache_dir)) cache_dir <- tempdir()

  Rcpp::sourceCpp(path, env = env, rebuild = rebuild,
                  cacheDir = cache_dir, verbose = !quiet)

  if (!exists("sparse_solver", envir = env, inherits = FALSE)) {
    stop("'", basename(path), "' does not define sparse_solver(). ",
         "Was it generated by create_model_sparse()?")
  }

  .tresthor_compiled[[key]] <- env
  env
}

#' Residuals of a model's equations at one observation
#'
#' Evaluates every equation of the model at a given row of the data and
#' returns the largest absolute residual per block. Useful to check that a
#' solution really does satisfy the model, and to compare solvers.
#'
#' @param model a `thoR.model` built by `create_model_sparse()`
#' @param database data.frame holding the data
#' @param periods the periods to check, as found in `index_time`. Default: all
#'   rows except the first.
#' @param index_time name of the time column in `database`
#' @return a matrix of maximum absolute residuals, periods by block
#' @export
model_residuals <- function(model, database, periods = NULL, index_time = "date") {

  env <- compile_model_cpp(model@rcpp_source)

  key <- as.character(database[[index_time]])
  if (is.null(periods)) {
    rows <- seq_along(key)[-1L]
  } else {
    rows <- match(as.character(periods), key)
    if (anyNA(rows)) {
      stop("Periods not found in '", index_time, "': ",
           paste(as.character(periods)[is.na(rows)], collapse = ", "))
    }
  }

  mv <- sort(c(model@exo_list, model@endo_list, model@coeff_list))
  M <- as.matrix(database[, mv, drop = FALSE])
  storage.mode(M) <- "double"

  f <- env$sparse_residuals
  if (!is.function(f)) {
    stop("'", basename(model@rcpp_source), "' does not define sparse_residuals(). ",
         "Rebuild the model with create_model_sparse().")
  }

  first <- f(M, as.integer(rows[1] - 1L))
  out <- t(vapply(rows, function(r) f(M, as.integer(r - 1L)), first))
  rownames(out) <- key[rows]
  out
}
