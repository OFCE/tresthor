## Sparse symbolic jacobian.
##
## `symbolic_jacobian()` (1_5) materialises an n x n character matrix and the
## downstream writers then loop over all n^2 cells. On ThreeME 8x8 the heart
## block alone is 2517 x 2517 = 6.3M cells holding 9k real entries; at 25k
## equations that representation cannot be built at all.
##
## Here the jacobian is kept as triplets (i, j, expr) from the start. Nothing
## proportional to n^2 is ever allocated.

#' Compute a block's symbolic jacobian in triplet form
#'
#' Only the partial derivatives that can be non-zero are computed: for each
#' equation, we differentiate with respect to the contemporaneous endogenous
#' variables that actually occur in it.
#'
#' @param equations_list_df data.frame of equations, as built by `create_model`
#' @param eqns_vars_list named list: equation id -> variables occurring in it
#' @param endo_vec character vector of the block's endogenous variables
#' @param equations_subset character vector of the block's equation ids
#' @param id_col column of `equations_list_df` holding the equation id
#' @param formula_col column of `equations_list_df` holding the formula
#'
#' @return a `thoR.sparse_jacobian`: list with `n`, `equations`, `endo`,
#'   `i`, `j` (1-based, CSC order) and `expr`.
#' @import Deriv
#' @keywords internal
symbolic_jacobian_sparse <- function(equations_list_df, eqns_vars_list, endo_vec,
                                     equations_subset, id_col = "id",
                                     formula_col = "new_formula") {

  ## `delta` is opaque to Deriv unless we declare its derivative rule: it is
  ## linear in its second argument, so d/dx delta(n, x) = 1 (the lagged part
  ## is carried by `lag.*` symbols, which Deriv treats as unrelated atoms).
  delta <- function(n, x, y = TRUE) NULL
  drule[["delta"]] <- alist(x = 1, y = NULL)

  eqs  <- sort(equations_subset)
  endo <- sort(endo_vec)

  if (length(eqs) != length(endo)) {
    stop(sprintf("Block is not square: %d equations for %d endogenous variables.",
                 length(eqs), length(endo)))
  }

  ids <- as.character(equations_list_df[[id_col]])
  pos <- match(eqs, ids)
  if (anyNA(pos)) {
    stop("Some equations of the block were not found in the equation list: ",
         paste(eqs[is.na(pos)], collapse = ", "))
  }
  formulas <- as.character(equations_list_df[[formula_col]])[pos]

  ## column lookup: endogenous variable -> index within the block
  endo_col <- seq_along(endo)
  names(endo_col) <- endo

  ## Upper bound on the number of entries, so we fill pre-allocated vectors
  ## instead of growing them equation by equation.
  cand <- lapply(eqs, function(id) {
    v <- eqns_vars_list[[id]]
    intersect(unique(as.character(v)), endo)
  })
  cap <- sum(lengths(cand))

  out_i <- integer(cap)
  out_j <- integer(cap)
  out_e <- character(cap)
  k <- 0L

  for (r in seq_along(eqs)) {
    vars <- cand[[r]]
    if (length(vars) == 0L) next
    f <- parse(text = formulas[r], keep.source = FALSE)
    for (v in vars) {
      d <- paste(Deriv::Deriv(f, v, cache.exp = FALSE))
      ## Deriv returns "0" for variables that only appear inside lag.* atoms
      if (identical(d, "0") || identical(d, "0L")) next
      k <- k + 1L
      out_i[k] <- r
      out_j[k] <- endo_col[[v]]
      out_e[k] <- d
    }
  }

  length(out_i) <- k
  length(out_j) <- k
  length(out_e) <- k

  ## canonical CSC order: by column, then by row. This is the layout Eigen
  ## uses, so the generated code can write straight into the value array.
  o <- order(out_j, out_i)

  structure(list(
    n         = length(endo),
    equations = eqs,
    endo      = endo,
    i         = out_i[o],
    j         = out_j[o],
    expr      = out_e[o]
  ), class = "thoR.sparse_jacobian")
}

#' @export
print.thoR.sparse_jacobian <- function(x, ...) {
  nnz <- length(x$i)
  cat(sprintf("Sparse symbolic jacobian: %d x %d, %d non-zero (%.4f%% dense), %.2f per row\n",
              x$n, x$n, nnz, 100 * nnz / (x$n^2), nnz / x$n))
  invisible(x)
}

#' An empty sparse jacobian, for blocks that are not used
#' @keywords internal
empty_sparse_jacobian <- function() {
  structure(list(n = 0L, equations = character(0), endo = character(0),
                 i = integer(0), j = integer(0), expr = character(0)),
            class = "thoR.sparse_jacobian")
}
