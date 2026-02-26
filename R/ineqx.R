# ============================================================================ #
# ineqx: Causal variance decomposition (main user-facing function)
# ============================================================================ #

#' Causal variance decomposition
#'
#' Decomposes the treatment effect on inequality into between-group and
#' within-group components. Supports both cross-sectional and longitudinal
#' designs, with automatic detection based on the input parameters.
#'
#' There are two usage modes:
#' \describe{
#'   \item{Mode A: Pre-extracted params}{Pass an \code{ineqx_params} object
#'     via \code{params}. Use \code{\link{ineqx_params}} to create one
#'     manually or from a fitted gamlss model.}
#'   \item{Mode B: Integrated fitting}{Pass \code{formula_mu},
#'     \code{formula_sigma}, \code{data}, \code{treat}, and \code{group}.
#'     The function fits a GAMLSS model, extracts parameters, and decomposes.}
#' }
#'
#' @param params An \code{ineqx_params} object for Mode A. If provided,
#'   formula/data arguments are ignored.
#' @param formula_mu Formula for the mean equation (Mode B). E.g.,
#'   \code{y ~ treat * group + controls}.
#' @param formula_sigma Formula for the log-SD equation (Mode B). One-sided,
#'   e.g., \code{~ treat * group + controls}.
#' @param data Data.frame for Mode B (individual-level data).
#' @param treat Character, treatment variable name (Mode B). Coded 0/1.
#' @param group Character, grouping variable name (Mode B).
#' @param time Character, time variable name (Mode B). NULL for cross-sectional.
#' @param weights Character, weight variable name (Mode B). NULL for equal weights.
#' @param post Character, pre/post indicator for DiD (Mode B). NULL for
#'   simple difference estimator.
#' @param ref Numeric, reference time period. Used in Mode B; in Mode A,
#'   taken from \code{params}.
#' @param ystat Character, \code{"Var"} or \code{"CV2"}. Default: \code{"Var"}.
#' @param order Decomposition ordering for longitudinal data. Either:
#'   \itemize{
#'     \item A character vector of length 3: a permutation of
#'       \code{c("behavioral", "compositional", "pretreatment")} for a
#'       single ordering.
#'     \item \code{"shapley"}: averages across all 6 possible orderings.
#'   }
#'   Ignored for cross-sectional data. Default:
#'   \code{c("behavioral", "compositional", "pretreatment")}.
#' @param se Standard error method. One of:
#'   \itemize{
#'     \item \code{"delta"} (default): delta method SEs (requires vcov in params)
#'     \item \code{"none"}: skip SE computation
#'     \item A \code{boot_config} object: bootstrap SEs
#'   }
#'
#' @return For cross-sectional data: an \code{ineqx_causal_cross} object.
#'   For longitudinal data with a specific ordering: an \code{ineqx_causal_longit}
#'   object. For \code{order = "shapley"}: an \code{ineqx_shapley} object.
#'
#' @examples
#' # Mode A: from manual params
#' params <- ineqx_params(
#'   data = data.frame(
#'     group = c("workers", "managers"),
#'     pi = c(0.5, 0.5),
#'     mu0 = c(500, 1000),
#'     sigma0 = c(200, 400),
#'     beta = c(60, 100),
#'     lambda = c(-0.1, -0.2)
#'   )
#' )
#' result <- ineqx(params = params)
#'
#' \dontrun{
#' # Mode B: integrated GAMLSS fitting
#' result <- ineqx(
#'   formula_mu = inc ~ x * factor(group),
#'   formula_sigma = ~ x * factor(group),
#'   data = incdat, treat = "x", group = "group"
#' )
#' }
#'
#' @export
ineqx <- function(params = NULL,
                  formula_mu = NULL, formula_sigma = NULL,
                  data = NULL, treat = NULL, group = NULL,
                  time = NULL, weights = NULL, post = NULL,
                  ref = NULL, ystat = "Var",
                  order = c("behavioral", "compositional", "pretreatment"),
                  se = "delta") {

  ystat <- match.arg(ystat, c("Var", "CV2"))

  # -------------------------------------------------------------------- #
  # Determine mode and obtain params
  # -------------------------------------------------------------------- #

  if (is.null(params)) {
    # Mode B: Integrated GAMLSS fitting
    if (is.null(formula_mu) || is.null(data) || is.null(treat) || is.null(group)) {
      stop("For Mode B (integrated fitting), 'formula_mu', 'data', 'treat', ",
           "and 'group' are required.\n",
           "Alternatively, pass pre-computed parameters via 'params'.")
    }
    if (is.null(formula_sigma)) {
      stop("'formula_sigma' is required for GAMLSS fitting.")
    }

    # Fit GAMLSS
    model <- fit_ineqx_model(
      formula_mu = formula_mu,
      formula_sigma = formula_sigma,
      data = data,
      weights = weights
    )

    # Determine whether to extract vcov (needed for delta method SEs)
    extract_vcov <- is.character(se) && se == "delta"

    # Extract params
    params <- ineqx_params(
      model = model, data = data,
      treat = treat, group = group,
      time = time, post = post,
      ystat = ystat, ref = ref,
      vcov = extract_vcov
    )

  } else {
    # Mode A: Pre-extracted params
    stopifnot(inherits(params, "ineqx_params"))
  }

  # -------------------------------------------------------------------- #
  # Parse se argument
  # -------------------------------------------------------------------- #

  boot <- NULL
  if (inherits(se, "ineqx_boot_config")) {
    boot <- se
    se_method <- "bootstrap"
  } else if (is.character(se)) {
    se <- match.arg(se, c("delta", "none"))
    se_method <- se
  } else {
    stop("'se' must be 'delta', 'none', or a boot_config() object.")
  }

  # For delta method, check vcov availability
  if (se_method == "delta" && is.null(params$vcov)) {
    message("vcov not available in params; falling back to se = 'none'.")
    se_method <- "none"
  }

  # Suppress vcov-based SEs in decomposition if not using delta method
  params_for_decomp <- params
  if (se_method != "delta") {
    params_for_decomp$vcov <- NULL
  }

  # -------------------------------------------------------------------- #
  # Parse order argument
  # -------------------------------------------------------------------- #

  use_shapley <- is.character(order) && length(order) == 1 && order == "shapley"

  if (!use_shapley) {
    valid_components <- c("behavioral", "compositional", "pretreatment")
    order <- match.arg(order, valid_components, several.ok = TRUE)
    if (length(order) != 3 || !setequal(order, valid_components)) {
      stop("'order' must be 'shapley' or a permutation of ",
           "c('behavioral', 'compositional', 'pretreatment')")
    }
  }

  # -------------------------------------------------------------------- #
  # Dispatch decomposition
  # -------------------------------------------------------------------- #

  if (params$type == "cross_sectional") {
    # Cross-sectional: order is ignored
    result <- causal_decompose_cross(params_for_decomp)

  } else if (use_shapley) {
    # Longitudinal with Shapley averaging
    result <- causal_shapley(params_for_decomp)

  } else {
    # Longitudinal with specific ordering
    result <- causal_decompose_longit(params_for_decomp, order = order)
  }

  # -------------------------------------------------------------------- #
  # Bootstrap SEs (if requested)
  # -------------------------------------------------------------------- #

  if (se_method == "bootstrap") {
    boot_group <- if (!is.null(boot$group)) boot$group else group
    boot_time  <- if (!is.null(boot$time))  boot$time  else time

    if (is.null(boot_group)) {
      stop("'group' must be specified either in boot_config() or ineqx().")
    }

    boot_result <- bootstrap_se(
      data = boot$data,
      formula_mu = boot$formula_mu,
      formula_sigma = boot$formula_sigma,
      treat = boot$treat,
      group = boot_group,
      time = boot_time,
      post = boot$post,
      ref = if (!is.null(params$ref)) params$ref else ref,
      ystat = params$ystat,
      order = if (use_shapley) c("behavioral", "compositional", "pretreatment") else order,
      B = boot$B,
      family = boot$family,
      parallel = boot$parallel,
      ncores = boot$ncores,
      seed = boot$seed,
      verbose = boot$verbose,
      gamlss_args = boot$gamlss_args
    )

    result$se <- boot_result$se
    result$se_method <- "bootstrap"
    result$boot <- boot_result
  } else if (se_method == "delta") {
    result$se_method <- "delta"
  } else {
    result$se_method <- "none"
  }

  result
}
