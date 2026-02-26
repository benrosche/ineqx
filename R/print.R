# ============================================================================ #
# Print and summary methods for ineqx result objects
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# ineqx_desc
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_desc <- function(x, ...) {
  cat("Descriptive variance decomposition\n")
  cat("Inequality measure:", x$ystat, "\n\n")

  cat("Totals by time:\n")
  if (x$ystat == "Var") {
    cols <- c("VarW", "VarB", "VarT")
  } else {
    cols <- c("CV2W", "CV2B", "CV2T")
  }
  if ("time" %in% names(x$totals)) cols <- c("time", cols)
  print(x$totals[, cols, drop = FALSE], row.names = FALSE, digits = 4)

  if (!is.null(x$deltas)) {
    cat("\nChange relative to ref =", x$ref, ":\n")
    # Show summary deltas (per-parameter totals)
    summary_cols <- c("time", "delta_mu_T", "delta_sigma_T", "delta_pi_T", "delta_T")
    summary_cols <- intersect(summary_cols, names(x$deltas))
    print(x$deltas[, summary_cols, drop = FALSE], row.names = FALSE, digits = 4)
  }

  invisible(x)
}

#' @export
summary.ineqx_desc <- function(object, ...) {
  cat("Descriptive variance decomposition\n")
  cat("Inequality measure:", object$ystat, "\n")
  if (!is.null(object$ref)) cat("Reference period:", object$ref, "\n")
  cat("\n--- Group-level statistics ---\n")
  print(object$wibe, row.names = FALSE, digits = 4)
  cat("\n--- Totals ---\n")
  print(object$totals, row.names = FALSE, digits = 4)
  if (!is.null(object$deltas)) {
    cat("\n--- Deltas (change from reference) ---\n")
    print(object$deltas, row.names = FALSE, digits = 4)
  }
  invisible(object)
}

# ---------------------------------------------------------------------------- #
# ineqx_causal_cross
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_causal_cross <- function(x, ...) {
  cat("Cross-sectional causal variance decomposition\n")
  cat("Inequality measure:", x$ystat, "\n\n")
  cat("Treatment effect on inequality:\n")

  if (!is.null(x$se)) {
    se_label <- if (!is.null(x$se_method)) x$se_method else "delta"
    cat("  (SEs via ", se_label, " method)\n", sep = "")
    .print_estimate_with_se("Between-group (delta_B)", x$delta_B, x$se$se_delta_B)
    .print_estimate_with_se("Within-group  (delta_W)", x$delta_W, x$se$se_delta_W)
    .print_estimate_with_se("Total         (delta_T)", x$delta_total, x$se$se_delta_total)
    cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  } else {
    cat("  Between-group (delta_B):", round(x$delta_B, 4), "\n")
    cat("  Within-group  (delta_W):", round(x$delta_W, 4), "\n")
    cat("  Total         (delta_T):", round(x$delta_total, 4), "\n")
  }

  cat("\nBetween-group sub-components:\n")
  cat("  Var_pi(beta):            ", round(x$components$Var_pi_beta, 4), "\n")
  cat("  2*Cov_pi(mu0, beta):     ", round(2 * x$components$Cov_pi_mu_beta, 4), "\n")

  cat("\nWithin-group sub-components:\n")
  cat("  mean(sigma0^2) * mean(f):", round(x$components$mean_sigma2_0 * x$components$mean_f, 4), "\n")
  cat("  Cov_pi(sigma0^2, f):     ", round(x$components$Cov_pi_sigma2_f, 4), "\n")

  invisible(x)
}

#' @export
summary.ineqx_causal_cross <- function(object, ...) {
  print(object, ...)
  cat("\n--- Group-level contributions ---\n")
  print(object$by_group, row.names = FALSE, digits = 4)
  invisible(object)
}

# ---------------------------------------------------------------------------- #
# ineqx_causal_longit
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_causal_longit <- function(x, ...) {
  cat("Longitudinal causal variance decomposition\n")
  cat("Inequality measure:", x$ystat, "\n")
  cat("Reference period:", x$ref, "\n")
  cat("Ordering:", paste(x$order, collapse = " -> "), "\n\n")

  has_se <- !is.null(x$se)

  # Build summary table
  times <- names(x$results)
  rows <- lapply(times, function(t) {
    r <- x$results[[t]]
    row <- data.frame(
      time = r$time,
      Delta_behavioral = round(r$Delta_behavioral, 4),
      Delta_compositional = round(r$Delta_compositional, 4),
      Delta_pretreatment = round(r$Delta_pretreatment, 4),
      Delta_total = round(r$Delta_total, 4),
      stringsAsFactors = FALSE
    )
    if (has_se && t %in% names(x$se)) {
      s <- x$se[[t]]
      row$se_behavioral <- round(s$se_Delta_behavioral, 4)
      row$se_compositional <- round(s$se_Delta_compositional, 4)
      row$se_pretreatment <- round(s$se_Delta_pretreatment, 4)
      row$se_total <- round(s$se_Delta_total, 4)
    }
    row
  })
  df <- do.call(rbind, rows)

  cat("Four-component decomposition:\n")
  print(df, row.names = FALSE)

  if (has_se) {
    se_label <- if (!is.null(x$se_method)) x$se_method else "delta"
    cat("(SEs via ", se_label, " method)\n", sep = "")
    cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  invisible(x)
}

#' @export
summary.ineqx_causal_longit <- function(object, ...) {
  cat("Longitudinal causal variance decomposition\n")
  cat("Inequality measure:", object$ystat, "\n")
  cat("Reference period:", object$ref, "\n")
  cat("Ordering:", paste(object$order, collapse = " -> "), "\n\n")

  times <- names(object$results)
  rows <- lapply(times, function(t) {
    r <- object$results[[t]]
    data.frame(
      time = r$time,
      Delta_beta = round(r$Delta_beta, 4),
      Delta_lambda = round(r$Delta_lambda, 4),
      Delta_pi_B = round(r$Delta_pi_B, 4),
      Delta_pi_W = round(r$Delta_pi_W, 4),
      Delta_mu = round(r$Delta_mu, 4),
      Delta_sigma = round(r$Delta_sigma, 4),
      Delta_total = round(r$Delta_total, 4),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)

  cat("Six-component decomposition:\n")
  print(df, row.names = FALSE)

  invisible(object)
}

# ---------------------------------------------------------------------------- #
# ineqx_shapley
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_shapley <- function(x, ...) {
  cat("Shapley values for longitudinal causal decomposition\n")
  cat("Inequality measure:", x$ystat, "\n")
  cat("Reference period:", x$ref, "\n")
  cat("(Averaged across all 6 orderings)\n\n")

  # Show 4-component Shapley values
  cols <- c("time", "Delta_behavioral", "Delta_compositional",
            "Delta_pretreatment", "Delta_total")
  df <- x$shapley[, cols]
  df[, -1] <- round(df[, -1], 4)
  cat("Four-component Shapley values:\n")
  print(df, row.names = FALSE)

  invisible(x)
}

#' @export
summary.ineqx_shapley <- function(object, ...) {
  print(object, ...)

  cat("\n--- Six-component Shapley values ---\n")
  cols6 <- c("time", "Delta_beta", "Delta_lambda",
             "Delta_pi_B", "Delta_pi_W",
             "Delta_mu", "Delta_sigma", "Delta_total")
  df6 <- object$shapley[, cols6]
  df6[, -1] <- round(df6[, -1], 4)
  print(df6, row.names = FALSE)

  cat("\n--- Ranges across orderings ---\n")
  range_cols <- c("time",
                  "Delta_behavioral_range", "Delta_compositional_range",
                  "Delta_pretreatment_range")
  if (all(range_cols %in% names(object$ranges))) {
    df_r <- object$ranges[, range_cols]
    df_r[, -1] <- round(df_r[, -1], 6)
    print(df_r, row.names = FALSE)
  }

  invisible(object)
}



# ---------------------------------------------------------------------------- #
# Internal helpers for SE display
# ---------------------------------------------------------------------------- #

.print_estimate_with_se <- function(label, estimate, se) {
  pval <- .pval_from_se(estimate, se)
  stars <- .signif_stars(pval)
  cat(sprintf("  %s: %10.4f  (SE = %.4f) %s\n", label, estimate, se, stars))
}
