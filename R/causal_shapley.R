# ============================================================================ #
# Shapley value computation for longitudinal causal decomposition
# ============================================================================ #

#' Compute Shapley values for the longitudinal causal decomposition
#'
#' Computes the average decomposition across all 6 possible orderings of
#' the three component types (behavioral, compositional, pretreatment).
#' Since orderings are paired (between-group and within-group use the same
#' structural ordering), there are exactly 6 evaluations.
#'
#' Shapley values provide a robustness check against path dependence in the
#' sequential parameter-switching decomposition. When changes over time are
#' small relative to levels, the ordering has little practical impact. When
#' changes are large, Shapley values provide a useful robustness check.
#'
#' @param params An \code{ineqx_params} object with multiple time periods
#'
#' @return An object of class \code{"ineqx_causal_longit"} with
#'   \code{order = "shapley"} and additional fields:
#' \describe{
#'   \item{shapley}{data.frame with Shapley-averaged values for each component
#'     at each time period}
#'   \item{all_orderings}{List of 6 \code{ineqx_causal_longit} results,
#'     one per ordering}
#'   \item{ranges}{data.frame showing min/max for each component across orderings}
#' }
#'
#' @keywords internal
causal_shapley <- function(params, ref = NULL) {

  stopifnot(inherits(params, "ineqx_params"))
  stopifnot(params$type == "longitudinal")

  if (is.null(ref)) {
    stop("'ref' is required for longitudinal decomposition.")
  }

  # Generate all 6 permutations of the 3 component types
  perms <- .all_permutations(c("behavioral", "compositional", "pretreatment"))

  # Run decomposition for each ordering
  all_orderings <- lapply(perms, function(ord) {
    causal_decompose_longit(params, order = ord, ref = ref)
  })
  names(all_orderings) <- vapply(perms, paste, character(1), collapse = "->")

  # Extract component names
  # Split components (10): each parameter switch contributes to both tau_B and tau_W.
  # For V the off-diagonal parts are 0 by construction; for CV2 they capture
  # the shared-denominator coupling.
  component_names <- c(
    # Split (10)
    "Delta_beta_B", "Delta_beta_W",
    "Delta_lambda_B", "Delta_lambda_W",
    "Delta_pi_B", "Delta_pi_W",
    "Delta_mu_B", "Delta_mu_W",
    "Delta_sigma_B", "Delta_sigma_W",
    # Aggregate parameter components (4; pi already in split)
    "Delta_beta", "Delta_lambda", "Delta_mu", "Delta_sigma",
    # 3 combined components
    "Delta_behavioral", "Delta_compositional", "Delta_pretreatment",
    # 4-component grouping
    "Delta_pi", "Delta_pre",
    # Between/within and grand totals
    "Delta_B", "Delta_W", "Delta_total"
  )

  # Compute Shapley values (average) and ranges for each time period
  times <- setdiff(params$times, ref)
  shapley_list <- list()
  ranges_list <- list()
  results <- list()

  # Check if delta-method SEs are available (from any ordering)
  has_se <- !is.null(all_orderings[[1]]$se)
  se_names <- if (has_se) {
    paste0("se_", component_names)
  } else {
    character(0)
  }
  se_list <- list()

  for (t in times) {
    t_char <- as.character(t)

    # Collect values across all orderings
    values <- matrix(NA, nrow = length(all_orderings), ncol = length(component_names))
    colnames(values) <- component_names

    # Collect SEs across orderings (if available)
    se_values <- if (has_se) {
      matrix(NA, nrow = length(all_orderings), ncol = length(se_names))
    } else {
      NULL
    }

    for (i in seq_along(all_orderings)) {
      r <- all_orderings[[i]]$results[[t_char]]
      for (cn in component_names) {
        values[i, cn] <- r[[cn]]
      }
      # Collect SEs
      if (has_se && !is.null(all_orderings[[i]]$se[[t_char]])) {
        se_t <- all_orderings[[i]]$se[[t_char]]
        for (si in seq_along(se_names)) {
          se_val <- se_t[[se_names[si]]]
          if (!is.null(se_val)) se_values[i, si] <- se_val
        }
      }
    }

    # Shapley = mean across orderings
    shapley_vals <- colMeans(values)
    shapley_row <- as.data.frame(as.list(shapley_vals))
    shapley_row$time <- t
    shapley_list[[length(shapley_list) + 1]] <- shapley_row

    # Build results entry (same structure as causal_decompose_longit)
    # Cross-sectional components are ordering-invariant; take from first ordering
    ref_result <- all_orderings[[1]]$results[[t_char]]
    result_t <- c(
      list(time = t),
      as.list(shapley_vals),
      list(
        tau_B_t0 = ref_result$tau_B_t0,
        tau_W_t0 = ref_result$tau_W_t0,
        tau_B_t  = ref_result$tau_B_t,
        tau_W_t  = ref_result$tau_W_t,
        components_t0 = ref_result$components_t0,
        components_t  = ref_result$components_t
      )
    )
    results[[t_char]] <- result_t

    # Average SEs across orderings
    if (has_se && !is.null(se_values)) {
      avg_se <- colMeans(se_values, na.rm = TRUE)
      se_row <- as.list(avg_se)
      names(se_row) <- se_names
      se_row$time <- t
      se_list[[length(se_list) + 1]] <- as.data.frame(se_row)
    }

    # Ranges
    range_row <- data.frame(time = t, stringsAsFactors = FALSE)
    for (cn in component_names) {
      range_row[[paste0(cn, "_min")]] <- min(values[, cn])
      range_row[[paste0(cn, "_max")]] <- max(values[, cn])
      range_row[[paste0(cn, "_range")]] <- max(values[, cn]) - min(values[, cn])
    }
    ranges_list[[length(ranges_list) + 1]] <- range_row
  }

  shapley_df <- do.call(rbind, shapley_list)
  ranges_df <- do.call(rbind, ranges_list)

  # Reorder columns to put time first
  shapley_df <- shapley_df[, c("time", component_names)]
  rownames(shapley_df) <- NULL
  rownames(ranges_df) <- NULL

  # Build SE named list (same structure as longit: list per time period)
  shapley_se <- NULL
  if (has_se && length(se_list) > 0) {
    se_df <- do.call(rbind, se_list)
    shapley_se <- list()
    for (i in seq_len(nrow(se_df))) {
      t_char <- as.character(se_df$time[i])
      shapley_se[[t_char]] <- as.list(se_df[i, setdiff(names(se_df), "time")])
    }
  }

  # Cross-sectional SEs and Wald tests are ordering-invariant — take from first ordering
  cross_se <- all_orderings[[1]]$cross_se
  lambda_tests <- all_orderings[[1]]$lambda_tests
  beta_tests <- all_orderings[[1]]$beta_tests

  structure(
    list(
      results = results,
      se = shapley_se,
      cross_se = cross_se,
      lambda_tests = lambda_tests,
      beta_tests = beta_tests,
      order = "shapley",
      ystat = params$ystat,
      ref = ref,
      params = params,
      shapley = shapley_df,
      all_orderings = all_orderings,
      ranges = ranges_df
    ),
    class = "ineqx_causal_longit"
  )
}


#' Generate all permutations of a vector
#' @keywords internal
.all_permutations <- function(x) {
  if (length(x) == 1) return(list(x))
  result <- list()
  for (i in seq_along(x)) {
    rest <- .all_permutations(x[-i])
    for (perm in rest) {
      result <- c(result, list(c(x[i], perm)))
    }
  }
  result
}
