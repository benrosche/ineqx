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
#' @return An object of class \code{"ineqx_shapley"} containing:
#' \describe{
#'   \item{shapley}{data.frame with Shapley-averaged values for each component
#'     at each time period}
#'   \item{all_orderings}{List of 6 \code{ineqx_causal_longit} results,
#'     one per ordering}
#'   \item{ranges}{data.frame showing min/max for each component across orderings}
#'   \item{params}{The input ineqx_params object}
#' }
#'
#' @keywords internal
causal_shapley <- function(params) {

  stopifnot(inherits(params, "ineqx_params"))
  stopifnot(params$type == "longitudinal")

  # Generate all 6 permutations of the 3 component types
  perms <- .all_permutations(c("behavioral", "compositional", "pretreatment"))

  # Run decomposition for each ordering
  all_orderings <- lapply(perms, function(ord) {
    causal_decompose_longit(params, order = ord)
  })
  names(all_orderings) <- vapply(perms, paste, character(1), collapse = "->")

  # Extract component names
  component_names <- c("Delta_beta", "Delta_lambda",
                       "Delta_pi_B", "Delta_pi_W",
                       "Delta_mu", "Delta_sigma",
                       "Delta_behavioral", "Delta_compositional",
                       "Delta_pretreatment",
                       "Delta_B", "Delta_W", "Delta_total")

  # Compute Shapley values (average) and ranges for each time period
  times <- setdiff(params$times, params$ref)
  shapley_list <- list()
  ranges_list <- list()

  for (t in times) {
    t_char <- as.character(t)

    # Collect values across all orderings
    values <- matrix(NA, nrow = length(all_orderings), ncol = length(component_names))
    colnames(values) <- component_names

    for (i in seq_along(all_orderings)) {
      r <- all_orderings[[i]]$results[[t_char]]
      for (cn in component_names) {
        values[i, cn] <- r[[cn]]
      }
    }

    # Shapley = mean across orderings
    shapley_vals <- colMeans(values)
    shapley_row <- as.data.frame(as.list(shapley_vals))
    shapley_row$time <- t
    shapley_list[[length(shapley_list) + 1]] <- shapley_row

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

  structure(
    list(
      shapley = shapley_df,
      all_orderings = all_orderings,
      ranges = ranges_df,
      ystat = params$ystat,
      ref = params$ref,
      params = params
    ),
    class = "ineqx_shapley"
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
