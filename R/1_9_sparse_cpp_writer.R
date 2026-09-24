## Generate the model-specific C++ solver, sparse from end to end.
##
## What this replaces (1_7_1 / 1_7_2 / 1_7):
##   * `mat Jacobian_n(n, n, fill::zeros)` allocated on every Newton iteration
##     of every period, then scanned element-by-element to build a sparse copy.
##   * the constant entries of the jacobian passed as a runtime-parsed
##     `arma::uvec("0 5 19 ...")` string literal.
##   * `spsolve()`, which redoes the fill-reducing ordering and the symbolic
##     analysis on every call.
##
## The generated code instead holds one Eigen sparse matrix per block whose
## pattern is fixed at build time. Each Newton iteration only overwrites the
## numerical values of the entries that actually vary, and the fill-reducing
## ordering is computed once per block for the whole simulation.

## Statements per generated function. Compilers are superlinear in the size of
## a single function body, so large blocks are split across several.
CPP_CHUNK <- 1500L

#' Format an integer vector as a C array initialiser
#' @keywords internal
cpp_int_array <- function(x, per_line = 20L) {
  if (length(x) == 0L) return("0")
  grp <- split(x, ceiling(seq_along(x) / per_line))
  paste(vapply(grp, function(g) paste(g, collapse = ","), character(1)),
        collapse = ",\n")
}

#' Format a double vector as a C array initialiser
#' @keywords internal
cpp_dbl_array <- function(x, per_line = 12L) {
  if (length(x) == 0L) return("0.0")
  s <- format(x, digits = 17L, scientific = FALSE, trim = TRUE)
  grp <- split(s, ceiling(seq_along(s) / per_line))
  paste(vapply(grp, function(g) paste(g, collapse = ","), character(1)),
        collapse = ",\n")
}

#' Split a vector of C++ statements into chunked functions
#'
#' @param stmts character vector of statement lines
#' @param fname base name of the generated function
#' @param sig argument list of the generated function
#' @param callargs arguments forwarded to each chunk
#' @return character vector of C++ lines defining `fname` and its chunks
#' @keywords internal
cpp_chunked_function <- function(stmts, fname, sig, callargs) {
  if (length(stmts) == 0L) {
    return(c(paste0("static void ", fname, "(", sig, "){ (void)M; (void)t; }"), ""))
  }
  idx <- split(seq_along(stmts), ceiling(seq_along(stmts) / CPP_CHUNK))
  out <- character(0)
  for (c_i in seq_along(idx)) {
    out <- c(out,
             paste0("static void ", fname, "_", c_i - 1L, "(", sig, "){"),
             stmts[idx[[c_i]]],
             "}", "")
  }
  out <- c(out,
           paste0("static void ", fname, "(", sig, "){"),
           paste0("  ", fname, "_", seq_along(idx) - 1L, "(", callargs, ");"),
           "}", "")
  out
}

#' Emit the C++ for one block of the model
#'
#' @param name block name ("prologue", "heart", "epilogue")
#' @param jac a `thoR.sparse_jacobian`
#' @param formulas character vector of the block's residual formulas, in the
#'   same order as `jac$equations`
#' @param vidx named integer vector: variable name -> 0-based column in M
#' @keywords internal
cpp_emit_block <- function(name, jac, formulas, vidx) {

  n   <- jac$n
  nnz <- length(jac$i)

  ## --- CSC pattern ---------------------------------------------------------
  ## jac$i / jac$j are already sorted by (column, row).
  Ai <- jac$i - 1L                                   # 0-based row indices
  Ap <- c(0L, cumsum(tabulate(jac$j, nbins = n)))    # column pointers

  ## --- constant vs varying entries ----------------------------------------
  num <- suppressWarnings(as.numeric(jac$expr))
  is_const <- !is.na(num)
  cst_k <- which(is_const) - 1L                      # 0-based slot in valuePtr
  cst_v <- num[is_const]
  vary  <- which(!is_const)

  ## --- statements ----------------------------------------------------------
  jac_stmts <- vapply(vary, function(k) {
    paste0("v[", k - 1L, "]=", cpp_from_formula(jac$expr[k], vidx), ";")
  }, character(1))

  res_stmts <- vapply(seq_along(formulas), function(k) {
    paste0("f[", k - 1L, "]=", cpp_from_formula(formulas[k], vidx), ";")
  }, character(1))

  ## --- block endogenous columns in the full data matrix --------------------
  endo_cols <- vidx[jac$endo]
  if (anyNA(endo_cols)) {
    stop("Block '", name, "': endogenous variables missing from the variable map.")
  }

  c(
    paste0("// ===================== block: ", name, " ====================="),
    sprintf("static const int %s_n   = %d;", name, n),
    sprintf("static const int %s_nnz = %d;", name, nnz),
    sprintf("static const int %s_ncst = %d;", name, length(cst_k)),
    sprintf("static const int %s_Ap[] = {\n%s};", name, cpp_int_array(Ap)),
    sprintf("static const int %s_Ai[] = {\n%s};", name, cpp_int_array(Ai)),
    sprintf("static const int %s_cst_k[] = {\n%s};", name, cpp_int_array(cst_k)),
    sprintf("static const double %s_cst_v[] = {\n%s};", name, cpp_dbl_array(cst_v)),
    sprintf("static const int %s_endo[] = {\n%s};", name, cpp_int_array(as.integer(endo_cols))),
    "",
    cpp_chunked_function(jac_stmts, paste0(name, "_jac"),
                         "const MapMat& M, int t, double* v", "M,t,v"),
    cpp_chunked_function(res_stmts, paste0(name, "_res"),
                         "const MapMat& M, int t, double* f", "M,t,f"),
    ""
  )
}

## The block-invariant part of the generated file: sparse Newton with a
## once-per-block symbolic factorisation, residual-based convergence,
## backtracking damping and a hard iteration cap.
CPP_RUNTIME <- '
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]
//
// C++17 is pinned deliberately. R >= 4.5 compiles Rcpp sources as gnu++20 by
// default, and on clang that makes this file take ~25x longer to compile
// (394 s versus 16 s for ThreeME 4x4) for no benefit here.
//
// Only the Eigen headers actually needed are included: pulling in the whole
// RcppEigen.h umbrella drags in the unsupported modules and costs seconds per
// translation unit.
#include <Rcpp.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include <cstdio>

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::SparseMatrix<double>  SpMat;

struct Block {
  const char* name;
  int n, nnz, ncst;
  const int* Ap; const int* Ai;
  const int* cst_k; const double* cst_v;
  const int* endo;
  void (*jac)(const MapMat&, int, double*);
  void (*res)(const MapMat&, int, double*);

  SpMat A;
  Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> > lu;
  bool ready;

  Block() : ready(false) {}

  // Build the sparsity pattern once, plant the constant entries, and run the
  // fill-reducing ordering + symbolic factorisation a single time.
  void setup(){
    std::vector<double> z(nnz, 0.0);
    A = Eigen::Map<const SpMat>(n, n, nnz, Ap, Ai, z.data());
    if (A.nonZeros() != nnz)
      Rcpp::stop(std::string("block ") + name + ": sparsity pattern was not preserved.");
    double* v = A.valuePtr();
    for (int c = 0; c < ncst; ++c) v[cst_k[c]] = cst_v[c];
    // NB: Eigen only initialises SparseLU::m_info in factorize(), so info()
    // must not be consulted after analyzePattern() alone -- it would read
    // uninitialised memory. Failures surface on the first factorize().
    lu.analyzePattern(A);
    ready = true;
  }
};

// Scale-free convergence measure. Macro models mix variables of very
// different magnitude -- in ThreeME 4x4 the heart block spans 0 to 1.2e7, and
// 338 of its variables are identically zero in a given scenario -- so neither
// a purely absolute nor a purely relative step test works.
//
// This is the standard mixed criterion: a step is acceptable when
//     |dx_i| <= rtol*|x_i| + atol
// which is returned in normalised form, so <= 1 means converged. Writing it
// instead as |dx_i|/(atol+|x_i|) <= rtol would be 1/rtol times stricter on
// variables sitting at zero, where it would measure nothing but rounding
// noise (dx ~ 5e-16 against a 1e-8 floor reads as 5e-8).
static double step_measure(const Eigen::VectorXd& dx, const Eigen::VectorXd& x,
                           double rtol, double atol)
{
  double c = 0.0;
  for (int i = 0; i < dx.size(); ++i) {
    double d = std::fabs(dx[i]) / (rtol * std::fabs(x[i]) + atol);
    if (d > c) c = d;
  }
  return c;
}

// Newton on one block at one period. Returns the number of iterations used,
// or -1 if it did not converge.
static int newton_block(Block& B, MapMat& M, int t,
                        double rtol, double atol, int max_iter, bool damping,
                        double& final_resid, double& final_conv)
{
  const int n = B.n;
  Eigen::VectorXd x(n), f(n), dx(n), xtry(n);

  for (int i = 0; i < n; ++i) x[i] = M(t, B.endo[i]);

  B.res(M, t, f.data());
  double rnorm = f.lpNorm<Eigen::Infinity>();
  double conv = std::numeric_limits<double>::infinity();

  for (int it = 0; it < max_iter; ++it) {

    if (!(rnorm == rnorm))                      // NaN in the residual
      { final_resid = rnorm; final_conv = conv; return -1; }

    B.jac(M, t, B.A.valuePtr());
    B.lu.factorize(B.A);
    if (B.lu.info() != Eigen::Success) { final_resid = rnorm; final_conv = conv; return -1; }

    dx = B.lu.solve(f);
    if (B.lu.info() != Eigen::Success) { final_resid = rnorm; final_conv = conv; return -1; }

    // Full Newton step, backtracked only if it makes the residual worse.
    double lambda = 1.0;
    double new_rnorm = rnorm;
    int tries = damping ? 12 : 1;
    bool ok = false;
    for (int b = 0; b < tries; ++b) {
      xtry = x - lambda * dx;
      for (int i = 0; i < n; ++i) M(t, B.endo[i]) = xtry[i];
      B.res(M, t, f.data());
      new_rnorm = f.lpNorm<Eigen::Infinity>();
      if (!damping || (new_rnorm == new_rnorm && new_rnorm <= rnorm)) { ok = true; break; }
      lambda *= 0.5;
    }

    conv = step_measure(lambda * dx, x, rtol, atol);

    if (!ok) {
      // No descent direction. Near the solution this just means the residual
      // has bottomed out in floating point, which is a success, not a failure.
      if (conv <= 1.0) { final_resid = rnorm; final_conv = conv; return it + 1; }
      for (int i = 0; i < n; ++i) M(t, B.endo[i]) = x[i];
      final_resid = rnorm; final_conv = conv;
      return -1;
    }

    x = xtry;
    rnorm = new_rnorm;

    if (conv <= 1.0) { final_resid = rnorm; final_conv = conv; return it + 1; }
  }

  final_resid = rnorm; final_conv = conv;
  return -1;
}
'

#' Generate the sparse C++ solver source for a model
#'
#' @param model_name name of the model
#' @param blocks named list of blocks; each element is
#'   `list(jac = <thoR.sparse_jacobian>, formulas = <character>)`, named
#'   "prologue" / "heart" / "epilogue". Blocks that are absent are skipped.
#' @param all_model_vars character vector of every model variable, sorted
#' @param rcpp_path directory in which to write the source file
#' @return the path of the generated file
#' @keywords internal
create_model_rcpp_sparse <- function(model_name, blocks, all_model_vars, rcpp_path) {

  rcpp_path <- normalizePath(rcpp_path, mustWork = TRUE)
  all_model_vars <- sort(all_model_vars)
  vidx <- seq_along(all_model_vars) - 1L          # 0-based columns of M
  names(vidx) <- all_model_vars

  f <- file.path(rcpp_path, paste0(model_name, "_sparse_newton.cpp"))

  active <- names(blocks)[vapply(blocks, function(b) !is.null(b) && b$jac$n > 0L, logical(1))]

  lines <- c(
    paste0("// Generated by tresthor for model '", model_name, "'. Do not edit."),
    paste0("// ", length(all_model_vars), " variables; blocks: ",
           paste(active, collapse = ", ")),
    CPP_RUNTIME,
    ""
  )

  for (nm in active) {
    b <- blocks[[nm]]
    cat("   - generating block '", nm, "' (", b$jac$n, " equations, ",
        length(b$jac$i), " jacobian entries)\n", sep = "")
    lines <- c(lines, cpp_emit_block(nm, b$jac, b$formulas, vidx))
  }

  ## --- registration and the time loop -------------------------------------
  reg <- unlist(lapply(active, function(nm) c(
    sprintf('  B[%d].name="%s"; B[%d].n=%s_n; B[%d].nnz=%s_nnz; B[%d].ncst=%s_ncst;',
            match(nm, active) - 1L, nm, match(nm, active) - 1L, nm,
            match(nm, active) - 1L, nm, match(nm, active) - 1L, nm),
    sprintf('  B[%d].Ap=%s_Ap; B[%d].Ai=%s_Ai; B[%d].cst_k=%s_cst_k; B[%d].cst_v=%s_cst_v;',
            match(nm, active) - 1L, nm, match(nm, active) - 1L, nm,
            match(nm, active) - 1L, nm, match(nm, active) - 1L, nm),
    sprintf('  B[%d].endo=%s_endo; B[%d].jac=&%s_jac; B[%d].res=&%s_res;',
            match(nm, active) - 1L, nm, match(nm, active) - 1L, nm,
            match(nm, active) - 1L, nm)
  )))

  lines <- c(lines, c(
    sprintf("static const int NB = %d;", length(active)),
    sprintf("static const char* BLOCK_NAMES[%d] = {%s};", length(active),
            paste0('"', active, '"', collapse = ",")),
    "",
    "static void tresthor_register(Block* B){",
    reg,
    "}",
    "",
    "// Maximum absolute residual of each block at one observation. Lets the",
    "// caller check that a solution really does satisfy the equations, and",
    "// compare solvers on the same footing.",
    "// [[Rcpp::export]]",
    "Rcpp::NumericVector sparse_residuals(Rcpp::NumericMatrix data, int row)",
    "{",
    "  MapMat M(data.begin(), data.nrow(), data.ncol());",
    "  Block B[NB];",
    "  tresthor_register(B);",
    "  Rcpp::NumericVector out(NB);",
    "  for (int b = 0; b < NB; ++b) {",
    "    std::vector<double> f(B[b].n);",
    "    B[b].res(M, row, f.data());",
    "    double mx = 0.0;",
    "    for (int i = 0; i < B[b].n; ++i) { double a = std::fabs(f[i]); if (a > mx) mx = a; }",
    "    out[b] = mx;",
    "  }",
    "  out.attr(\"names\") = Rcpp::CharacterVector(BLOCK_NAMES, BLOCK_NAMES + NB);",
    "  return out;",
    "}",
    "",
    "// [[Rcpp::export]]",
    "Rcpp::List sparse_solver(Rcpp::NumericMatrix data,",
    "                         int first_date, int last_date,",
    "                         double rtol, double atol, int max_iter,",
    "                         bool damping, bool verbose)",
    "{",
    "  MapMat M(data.begin(), data.nrow(), data.ncol());",
    "  // Block holds an Eigen::SparseLU, which is not copyable, so this has to",
    "  // be a fixed array rather than a std::vector.",
    "  Block B[NB];",
    "  tresthor_register(B);",
    "  for (int b = 0; b < NB; ++b) B[b].setup();",
    "",
    "  Rcpp::IntegerVector iters(last_date - first_date + 1);",
    "  Rcpp::NumericVector resid(last_date - first_date + 1);",
    "  Rcpp::NumericVector conv(last_date - first_date + 1);",
    "",
    "  for (int t = first_date; t <= last_date; ++t) {",
    "    // carry the previous period forward as the starting point",
    "    for (int b = 0; b < NB; ++b)",
    "      for (int i = 0; i < B[b].n; ++i)",
    "        M(t, B[b].endo[i]) = M(t - 1, B[b].endo[i]);",
    "",
    "    int tot = 0; double worst = 0.0; double worstc = 0.0;",
    "    for (int b = 0; b < NB; ++b) {",
    "      double r = 0.0, c = 0.0;",
    "      int it = newton_block(B[b], M, t, rtol, atol, max_iter, damping, r, c);",
    "      if (it < 0) {",
    "        char buf[512];",
    "        std::snprintf(buf, sizeof(buf),",
    "          \"Newton did not converge on block '%s' at row %d after %d iterations: \"",
    "          \"scaled step %.3e (converges at 1.0, rtol %.3e), max |residual| %.3e. \"",
    "          \"Try a looser convergence_criteria, a higher max_iter, or check the data at that period.\",",
    "          B[b].name, t + 1, max_iter, c, rtol, r);",
    "        Rcpp::stop(buf);",
    "      }",
    "      tot += it; if (r > worst) worst = r; if (c > worstc) worstc = c;",
    "    }",
    "    iters[t - first_date] = tot;",
    "    resid[t - first_date] = worst;",
    "    conv[t - first_date] = worstc;",
    "    if (verbose) Rcpp::Rcout << \"  \" << (t + 1) << \" (\" << tot << \" it) \";",
    "  }",
    "  if (verbose) Rcpp::Rcout << std::endl;",
    "",
    "  return Rcpp::List::create(Rcpp::_[\"data\"] = data,",
    "                            Rcpp::_[\"iterations\"] = iters,",
    "                            Rcpp::_[\"residuals\"] = resid,",
    "                            Rcpp::_[\"convergence\"] = conv);",
    "}"
  ))

  writeLines(lines, f)
  f
}
