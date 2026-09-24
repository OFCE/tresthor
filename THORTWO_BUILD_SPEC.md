# thortwo — build instructions

Instructions to myself for building `thortwo`, a new R package extracted from
`tresthor` (branch `improvements_2026`, v1.9.0) containing **only** the model
solver: parse a `.txt` model file, build a model object, solve it against a
database, store the results. Nothing else.

Read this first, then `git log --oneline` in the tresthor repo for the five
commits that introduced the sparse path — they are the reference implementation
for half of what follows.

---

## 1. Scope

**In.** The pipeline, end to end:

```
model .txt  ──parse──>  equations table + endo/exo/coeff vectors
            ──decompose──>  prologue / heart / epilogue blocks
            ──differentiate──>  symbolic jacobian (dense matrix OR sparse triplets)
            ──generate──>  C++ source  (OR R closures, for the dense-R path)
            ──compile──>  .so, cached or not
            ──solve──>  Newton per block per period
            ──store──>  data.frame + diagnostics, saveable
```

Both solver generations, side by side and interchangeable:

- **classic** — dense jacobian, `create_model()` / `thor_solver()`, with an R
  backend and an RcppArmadillo backend (`rcpp = TRUE`).
- **sparse** — `create_model_sparse()` / `thor_solver_sparse()`, RcppEigen,
  triplet jacobians, generated C++ only.

Plus, as first-class and independently switchable: the **compile cache**
(on/off, location) and **persistence** (save/load a built model, or rebuild
from source every time).

**Out.** Everything else in tresthor: model modification (`4_0`, 901 lines),
analysis (`6_0`), single-equation solving and estimation (`7_*`, `9_0`),
plots and contributions (`8_*`), the dictionary/var_map, the Opale HTML
guides, the shipped datasets. Do not port these. Do not port "just the small
useful bit" of them either — the point of the new package is a boundary.

The one exception worth arguing about is `var_map` / `build_dico` (`1_6`,
53 lines): it is already optional behind `no_var_map`. **Drop it.** It is a
model-inspection feature, not a solver feature, and keeping it drags
`splitstackshape` in.

---

## 2. What to port, file by file

Source paths are in the tresthor repo. "Port" means copy and clean; "rewrite"
means the existing code is not a good starting point.

| New file | From | Notes |
|---|---|---|
| `R/aaa-classes.R` | `0_thor_classes.R` | Rewrite. See §4 — the model class needs to change. |
| `R/parse-file.R` | `1_2_check_model_files.R` | Port. The eight `check_*` / `is_in_formulas` / `acceptable_var_name` / `parser_lag_delta_check` functions. |
| `R/parse-equations.R` | `1_1_create_model_subfunctions.R` | Port. `create_equations_list`, `table_contemporaneous_endos`, `formatting_formulas`, `get_variables_from_string`. |
| `R/decompose.R` | `1_3_decomposition_algo.R` | Port as is. 127 lines, self-contained, works. |
| `R/jacobian-dense.R` | `1_5_symbolic_jacobian.R` | Port as is. |
| `R/jacobian-sparse.R` | `1_5b_sparse_jacobian.R` | Port as is. |
| `R/codegen-expr.R` | `1_8_cpp_codegen.R` | Port as is. The R-expression → C++ translator; shared by any generated backend. |
| `R/codegen-sparse.R` | `1_9_sparse_cpp_writer.R` | Port as is. |
| `R/codegen-dense.R` | `1_7_2_create_rcpp_source_nonsuperlu.R` | Port. **Drop** `1_7_1_create_rcpp.R` and `1_7_rcpp_source_builder.R` (the SuperLU variant) — dead weight, and `use.superlu` is a support burden for a library almost nobody has installed. |
| `R/codegen-rfuns.R` | `1_4_function_command_writers.R` | Port, but see defect D1: it must stop writing `.R` files to the working directory. |
| `R/build.R` | `1_0_create_model.R` + `1_10_create_model_sparse.R` | **Rewrite as one function.** See §3. |
| `R/compile.R` | `1_11_compile_model.R` | Port. The best-designed file in the set; keep its structure and its comments. |
| `R/solve-dense.R` | `2_0_solver.R` | Rewrite. 308 lines of two solvers interleaved; split them. |
| `R/solve-sparse.R` | `2_1_sparse_solver.R` | Port as is. |
| `R/residuals.R` | `1_11_compile_model.R:179` | Port `model_residuals`. Extend it to work on classic models too — currently sparse-only. |
| `R/persist.R` | `3_0_files_and_data_management.R` | Rewrite. See §5. |
| `R/data-checks.R` | `3_1_data_checks.R`, `3_2_NA_check.R` | Port. |
| `R/operators.R` | `5_0_side_functions.R` | Port `delta`, `newdiff`, `add_coeffs`. These are needed **at solve time**, not just at build time — `delta()` appears in model formulas. Easy to forget and the failure is obscure. |

**DESCRIPTION.** Measured against the ported scope only:

```
Imports: methods, stats, utils, tools,
         Deriv, Rcpp, RcppEigen, RcppArmadillo, Matrix,
         assertthat, stringr, purrr
LinkingTo: Rcpp, RcppEigen, RcppArmadillo
```

Dropped relative to tresthor: `ggplot2`, `tidyr`, `scales`, `cointReg`,
`gsubfn`, `splitstackshape`, `dplyr`. Two of those need work:

- **`dplyr`** is only there for `%>%`, in seven ported files. Replace with `|>`
  (needs R >= 4.1, fine) and set `Depends: R (>= 4.1)`. Do not add `magrittr`.
- **`purrr`** is 16 calls, almost all `purrr::map`. Replace with `lapply`/
  `vapply` and drop it too. `purrr::is_empty(x)` → `length(x) == 0`.

If `RcppArmadillo` can also go — by regenerating the classic C++ backend on
Eigen, reusing `codegen-expr.R` — do that; it halves the compiled-dependency
surface. Decide after the classic path is ported and passing, not before.

---

## 3. The build function

tresthor has two near-identical 200–350 line builders that diverge only from
step 4 onwards. Do not copy that. One function, one `backend` argument:

```r
thor_model(name,
           source      = NULL,        # path to .txt
           endogenous = NULL, exogenous = NULL,
           coefficients = NULL, equations = NULL,   # manual input instead
           backend  = c("sparse", "dense-cpp", "dense-r"),
           decompose = TRUE,
           workdir  = NULL,           # where generated C++ lands; see below
           compile  = TRUE,
           cache    = NULL)           # NULL = default dir, FALSE = off, or a path
```

Steps 1–3 (read, check, index, decompose) are backend-independent and run
once. Step 4 branches: dense matrix jacobian for `dense-*`, triplets for
`sparse`. Step 5 branches: `codegen-sparse` / `codegen-dense` / `codegen-rfuns`.
Step 6 builds the object. Step 7 compiles if asked.

`backend` is stored on the object, and `thor_solve()` dispatches on it, so
there is exactly one solve entry point:

```r
thor_solve(model, from, to, data,
           index_time = "date",
           rtol = 1e-10, atol = 1e-8, max_iter = 100L,
           damping = TRUE, cache = NULL,
           verbose = TRUE, diagnostics = FALSE)
```

Keep `thor_solver()` and `thor_solver_sparse()` as thin deprecated wrappers if
migration matters; otherwise don't.

**`workdir` default matters more than it looks.** See defect D5.

---

## 4. The model object

`thoR.model` has 24 slots, six of which are `function` slots holding
`print("none")` stubs whenever the backend is not `dense-r`, and three dense
jacobian slots holding `matrix()` whenever the backend is sparse. The sparse
path then smuggles its real jacobians through `attr(model, "sparse_jacobians")`
because there is no slot for them. That is the seam to fix.

Design: a small common class, with backend-specific payload in one slot.

```r
setClass("thor_model", slots = c(
  name        = "character",
  backend     = "character",
  equations   = "data.frame",   # id, name, formula, new_formula, part
  vars        = "list",         # endo, exo, coeff (all lower-case, sorted)
  blocks      = "list",         # per block: name, present, endo, equation ids
  jacobian    = "list",         # backend-specific: dense matrices or triplets
  generated   = "list",         # see §5 on paths
  meta        = "list"          # version, built-at, source hash, digest
))
```

Verified today: `attr()` on an S4 object *does* survive `saveRDS`, so the
current smuggling works — but it is invisible to `validObject()`, to `show()`
and to anyone reading the class definition. A real slot costs nothing.

Add a `setValidity()` and a `show()` method. tresthor has neither, which is
why a malformed model surfaces as an error three steps later.

---

## 5. Caching and persistence — the two axes to get right

These are the parts the user explicitly wants configurable, and both are
subtly broken in tresthor today.

### Compile cache

Port `tresthor_cache_dir()` / `clear_model_cache()` / `compile_model_cpp()`
from `1_11` essentially unchanged, including the `-g` stripping (that comment
block documents a measured 25x compile-time effect; keep it verbatim). Three
states, per the `cache` argument: default dir (`tools::R_user_dir`), a given
path, or `FALSE` for off.

**But fix D5 first**, or the cache does nothing in practice.

### Persistence

`save_model()` writes the model with `saveRDS` and `load_model()` reads it back
and `sourceCpp`s `model@rcpp_source`, an **absolute path** recorded at build
time. Verified today: that path points into whatever temp directory the model
was built in. So a saved model is machine-local and usually session-local — it
cannot be mailed to a colleague, committed, or reused after a reboot clears
`/tmp`.

For thortwo, a saved model must be self-contained:

```r
thor_save(model, path)     # one .rds holding the object AND the generated
                           # C++ source as a string in `generated$code`
thor_load(path, compile = TRUE, cache = NULL)
                           # writes the code to a work dir, compiles (hitting
                           # the cache if the code is unchanged), returns ready
                           # to solve
```

`generated$path` becomes a hint, not the source of truth; `generated$code` and
`generated$hash` are authoritative. Then "saving or not saving" is genuinely
free: a built model and a loaded-from-disk model are the same thing, and
whether compilation actually happens is decided by the cache, not by which of
the two you have.

Keep `export_model()` (model object → `.txt`) — it round-trips with the parser
and is the cheapest possible regression test. See §7.

---

## 6. Defects in tresthor to fix during the port

Do not port these forward. Each is a real bug, not a style preference.

- **D1 — build writes to the current working directory.**
  `1_0_create_model.R:233` does `dir.create("temp_paprfn",)` (note the stray
  comma), writes six `.R` files into it, `source()`s them, and `unlink`s the
  directory at the end — so a failed build leaves litter, two concurrent builds
  in one directory race, and a read-only working directory breaks the build
  outright. Generate the closures into a `tempfile()` directory, or better,
  `eval(parse(text = ...))` them into a dedicated environment and skip the
  filesystem entirely (the commented-out code at `1_0_create_model.R:251-256`
  shows someone already tried this; find out why it was abandoned before
  copying it).

- **D2 — global state mutation at build time.**
  Both builders set `options(stringsAsFactors = FALSE)` (`1_0:91`, `1_10:49`)
  and write into `Deriv`'s `drule` table (`1_0:93-94`, `1_10:51-54`) without
  restoring either. `stringsAsFactors` has been a no-op default since R 4.0 —
  just delete it. The `drule` entries for `delta` and `abs` are necessary, but
  set them with `on.exit()` restore, or register them in the package's own
  `.onLoad`.

- **D3 — the classic C++ solver has no per-model isolation.**
  `2_0_solver.R:279` calls `Rcpp::sourceCpp(model@rcpp_source)` on *every
  solve*, loading `Rcpp_solver` into the global environment. Two classic models
  in one session: the second silently overwrites the first, and the first is
  then solved with the wrong code. The sparse path already solved this —
  `1_11`'s `.tresthor_compiled` registry keyed by source path, each model in its
  own environment. Route the classic backend through the same registry. This is
  the single highest-value fix in the list.

- **D4 — `database = t_data` as a default argument** (`2_0_solver.R:29`) makes
  the solver depend on a global named `t_data`. Make `data` required.

- **D5 — the cross-session cache is defeated by the usual call pattern.**
  Rcpp keys its cache on the source *path*; `create_model_sparse()` writes to
  whatever `rcpp_path` it is given, and every caller (including both test
  scripts in `tests/`) passes a fresh `tempfile()`. Result: a full recompile
  every run, plus a ~1.3 MB cache entry left behind each time. Measured: the
  same model built twice with a *stable* path costs 7.0 s then 0.3 s; with
  tempdir paths it is 7 s every time.
  Fix: default `workdir` to a **stable, model-derived** location —
  `file.path(tools::R_user_dir("thortwo","cache"), "src", name)` — so caching
  works by default and a tempdir is something you opt into.

- **D6 — silent truncation risk in the parser.** `1_0_create_model.R:98-101`
  reads the model file by *fixed line numbers* (`model_input[2]`, `[5]`, `[8]`,
  `[11:length]`), which is why `inst/models/opale.txt` carries three
  `####DO NOT ADD BLANK LINES OR REMOVE THIS LINE####` banners. One stray blank
  line silently shifts every section. Parse by section header instead, keep
  reading the legacy layout when the banners are present, and error loudly
  rather than mis-parse.

- **D7 — `thor_solver()` mutates its `database` argument's time column**
  (character conversion) and converts back at the end. Work on a copy.

---

## 7. Verification

Port the two existing scripts (`tests/test_sparse_threeme.R`,
`tests/test_sparse_opale.R`) as the acceptance tests, and extend them across
the new axes. Fixtures to copy from tresthor: `inst/models/opale.txt`,
`inst/Opale/donnees_opale.rds`, `inst/Opale/coefficients_opale.rds`,
`tests/threeme_{4x4,8x8}_thor.txt`, `tests/data3me_{4x4,8x8}.rds`.

The matrix that must pass:

| | classic-R | classic-C++ | sparse |
|---|---|---|---|
| Opale, 40 quarters | ✓ | ✓ | ✓ |
| ThreeME 4x4 | slow, mark skip | ✓ | ✓ |
| ThreeME 8x8 | skip | ✓ | ✓ |
| cache on / off | n/a | ✓ | ✓ |
| save → load → solve | ✓ | ✓ | ✓ |

Four invariants, in increasing order of strength:

1. **Converged** — worst scaled Newton step ≤ 1.
2. **Residuals** — `model_residuals()` over every block and period < 1e-6.
   This is the one that catches codegen bugs; convergence alone does not.
3. **Backends agree** — median relative difference between any two backends on
   the same model < 1e-10.
4. **History reproduces** — solving over history with historical exogenous
   inputs returns the historical endogenous path. On Opale: median relative
   difference 1.4e-15, max 1.8e-11 for variables with |value| > 1.

Baselines measured on Opale (496 equations, 201/98/197 blocks, jacobians
0.93–2.6% dense) on an M-series Mac, R 4.6.1:

```
sparse:      build 8.9 s (6.6 s compile)   solve 0.010 s / 40 periods
classic C++: build 14.7 s                  solve 0.3 s
             → build 1.8x, solve 34x
cache warm:  build 2.6 s (compile 0.3 s)
```

If the port comes in materially slower than these, something regressed.

Two known data quirks to encode in the tests rather than rediscover:

- Opale's `contpib7_td_p51p_d` is stored as exactly `0.0` in one historical
  quarter, so any non-zero solution scores relative difference 1 against it.
  Compare only where |value| > 1e-3.
- A handful of ThreeME variables (`pds_*`, a price over a near-zero stock
  change) are genuinely ill-conditioned. Compare medians, not maxima.

Also: `export_model(thor_model(f)) |> thor_model()` must produce an identical
equations table. Cheap, and it pins the parser.

---

## 8. Order of work

1. Skeleton, DESCRIPTION, class + validity + show. No logic.
2. Parse and decompose, with D6 fixed. Test: Opale and ThreeME 8x8 parse to the
   expected block sizes (201/98/197 for Opale).
3. Sparse path end to end — jacobian, codegen, compile, solve. This is the
   known-good code; getting it green first gives a reference for everything
   after. Test: Opale invariants 1, 2, 4.
4. Compile cache with D5 fixed. Test: second build of the same model compiles
   in < 1 s.
5. Persistence, self-contained. Test: save, restart R, load, solve.
6. Classic C++ backend, through the shared registry (D3). Test: invariant 3
   against sparse.
7. Classic R backend (D1, D2). Slowest and least used — last.
8. Drop `purrr`/`dplyr`, then reconsider `RcppArmadillo`.

Milestone 3 is the useful package. Everything after is compatibility.
