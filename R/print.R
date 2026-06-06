# ============================================================================ #
# Print and summary methods for ineqx result objects
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# ineqx_desc
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_desc <- function(x, ...) {
  cat("Descriptive variance decomposition\n")
  cat("Inequality measure:", x$ystat, "\n")
  if (!is.null(x$ref)) cat("Reference period:", x$ref, "\n")
  if (!is.null(x$order)) cat("Ordering:", paste(x$order, collapse = " -> "), "\n")
  cat("\n")

  # --- Table 1: Totals by time ---
  cat("Totals by time:\n")
  if (x$ystat %in% c("Var", "VL")) {
    cols <- c("VarW", "VarB", "VarT")
  } else {
    cols <- c("CV2W", "CV2B", "CV2T")
  }
  if ("time" %in% names(x$totals)) cols <- c("time", cols)
  print(x$totals[, cols, drop = FALSE], row.names = FALSE, digits = 4)

  # --- Table 2: Per-time blocks ---
  if (!is.null(x$deltas) && !is.null(x$ref)) {
    if (!is.null(x$ref_params)) {
      cat(sprintf("\nDecomposition of changes in %s relative to counterfactual reference:\n", x$ystat))
    } else {
      cat(sprintf("\nDecomposition of changes in %s relative to ref = %s:\n", x$ystat, x$ref))
    }

    all_times <- sort(unique(x$deltas$time))
    for (t in all_times) {
      t_char <- as.character(t)
      if (t == x$ref) {
        cat(sprintf("\n  time %s: (reference)\n", t_char))
        next
      }
      row <- x$deltas[x$deltas$time == t, ]
      cat(sprintf("\n  time %s:\n", t_char))
      .print_desc_block(row, ystat = x$ystat)
    }
  }

  invisible(x)
}

#' @export
summary.ineqx_desc <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

# ---------------------------------------------------------------------------- #
# ineqx_causal_cross
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_causal_cross <- function(x, ...) {
  pad <- 40
  cat("Cross-sectional causal variance decomposition\n")
  cat("Inequality measure:", x$ystat, "\n")
  if (!is.null(x$se_method) && x$se_method != "none") {
    cat(sprintf("SEs via %s method\n", x$se_method))
  }
  stat_label <- if (x$ystat == "Var") "Var" else "CV2"
  stat_name <- if (x$ystat == "Var") "variance" else "CV2"
  cat(sprintf("\nTreatment effect on outcome %s:\n", stat_name))

  # NB: x[["se"]] avoids partial matching to x$se_method
  se <- x[["se"]]
  .print_line(sprintf("%s[Y | T = 0]", stat_label), x$ineq_control,
              se = if (!is.null(se)) se$se_ineq_control, pad = pad)
  .print_line(sprintf("%s[Y | T = 1]", stat_label), x$ineq_treated,
              se = if (!is.null(se)) se$se_ineq_treated, pad = pad)

  # Total, then B/W indented
  if (!is.null(se)) {
    .print_line("Total effect (tau_T)", x$tau_total,
                se = se$se_tau_total, pad = pad)
    .print_line("Between-group (tau_B)", x$tau_B,
                se = se$se_tau_B, indent = 6, pad = pad)
    .print_line("Within-group (tau_W)", x$tau_W,
                se = se$se_tau_W, indent = 6, pad = pad)
  } else {
    .print_line("Total effect (tau_T)", x$tau_total, pad = pad)
    .print_line("Between-group (tau_B)", x$tau_B, indent = 6, pad = pad)
    .print_line("Within-group (tau_W)", x$tau_W, indent = 6, pad = pad)
  }

  # Wald test for lambda = 0
  if (!is.null(x$lambda_test)) {
    lt <- x$lambda_test
    p_fmt <- if (lt$p_value < 0.0001) "p < 0.0001" else sprintf("p = %.4f", lt$p_value)
    cat(sprintf("      -> Wald test H0: lambda = 0: chi2(%d) = %.4f, %s\n",
                lt$df, lt$statistic, p_fmt))
    if (lt$p_value < 0.05) {
      cat("         Significant: evidence of treatment effect heterogeneity\n")
    } else {
      cat("         Not significant: no evidence of treatment effect heterogeneity\n")
    }
  }

  # Wald test for beta homogeneity (effect on the mean equal across groups).
  # Null is scale-dependent: identity -> constant absolute effect; log ->
  # constant proportional (geometric-mean) effect.
  if (!is.null(x$beta_test)) {
    bt <- x$beta_test
    null_lab <- if (identical(bt$scale, "log")) "constant proportional effect"
                else "constant absolute effect"
    p_fmt <- if (bt$p_value < 0.0001) "p < 0.0001" else sprintf("p = %.4f", bt$p_value)
    cat(sprintf("      -> Wald test H0: beta homogeneous (%s): chi2(%d) = %.4f, %s\n",
                null_lab, bt$df, bt$statistic, p_fmt))
    if (bt$p_value < 0.05) {
      cat("         Significant: effect on the mean differs across groups\n")
    } else {
      cat("         Not significant: no evidence the mean effect differs across groups\n")
    }
  }

  comps <- x$components
  se_sub <- if (!is.null(se)) se$se_sub else NULL
  if (x$ystat == "Var") {
    cat("\nBetween-group sub-components:\n")
    .print_line("Var_pi(beta)", comps$Var_pi_beta,
                se = if (!is.null(se_sub)) se_sub$se_het_B, pad = pad)
    .print_line("2*Cov_pi(mu0, beta)", 2 * comps$Cov_pi_mu_beta,
                se = if (!is.null(se_sub)) se_sub$se_cov_B, pad = pad)
    cat("\nWithin-group sub-components:\n")
    .print_line("mean(sigma0^2) * mean(f)", comps$mean_sigma2_0 * comps$mean_f,
                se = if (!is.null(se_sub)) se_sub$se_het_W, pad = pad)
    .print_line("Cov_pi(sigma0^2, f)", comps$Cov_pi_sigma2_f,
                se = if (!is.null(se_sub)) se_sub$se_cov_W, pad = pad)
  } else {
    cat("\nBetween-group sub-components:\n")
    .print_line("Treatment heterogeneity", comps$het_B,
                se = if (!is.null(se_sub)) se_sub$se_het_B, pad = pad)
    .print_line("Scale effect", comps$cov_B,
                se = if (!is.null(se_sub)) se_sub$se_cov_B, pad = pad)
    cat("\nWithin-group sub-components:\n")
    .print_line("Treatment heterogeneity", comps$het_W,
                se = if (!is.null(se_sub)) se_sub$se_het_W, pad = pad)
    .print_line("Scale effect", comps$cov_W,
                se = if (!is.null(se_sub)) se_sub$se_cov_W, pad = pad)
  }

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
# ineqx_causal_longit (includes Shapley-averaged results)
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_causal_longit <- function(x, ...) {
  is_shapley <- identical(x$order, "shapley")

  cat("Longitudinal causal variance decomposition\n")
  cat("Inequality measure:", x$ystat, "\n")
  if (!is.null(x$se_method) && x$se_method != "none") {
    cat(sprintf("SEs via %s method\n", x$se_method))
  }
  cat("Reference period:", x$ref, "\n")
  if (is_shapley) {
    cat("Ordering: Shapley (averaged across all 6 orderings)\n\n")
  } else {
    cat("Ordering:", paste(x$order, collapse = " -> "), "\n\n")
  }

  # --- Table 1: Sub-component table ---
  # For Shapley, cross-sectional components are ordering-invariant
  comps_results <- if (is_shapley) x$all_orderings[[1]]$results else x$results
  cat("Treatment effect on variance by time:\n")
  cat("(tau_B = Var(beta) + 2Cov(mu,beta);  tau_W = E(sigma^2)*E(f) + Cov(sigma^2,f);  tau_T = tau_B + tau_W)\n\n")
  .build_longit_subcomponents(comps_results, x$ref,
                               cross_se = x$cross_se,
                               lambda_tests = x$lambda_tests,
                               beta_tests = x$beta_tests)

  # --- Table 2: Per-time blocks ---
  # NB: x[["se"]] avoids partial matching to x$se_method
  se <- x[["se"]]
  has_se <- !is.null(se)

  cat(sprintf("\nDecomposition of changes in group-specific treatment effects on variance relative to ref = %s:\n", x$ref))

  all_times <- .get_all_times(x$results, x$ref)
  for (t in all_times) {
    t_char <- as.character(t)
    if (t == x$ref) {
      cat(sprintf("\n  time %s: (reference)\n", t_char))
      next
    }
    r <- x$results[[t_char]]
    s <- if (has_se && t_char %in% names(se)) se[[t_char]] else NULL
    cat(sprintf("\n  time %s:\n", t_char))
    .print_causal_block(r, s)
  }

  invisible(x)
}

#' @export
summary.ineqx_causal_longit <- function(object, ...) {
  print(object, ...)
  invisible(object)
}


# ---------------------------------------------------------------------------- #
# Internal helpers
# ---------------------------------------------------------------------------- #

#' Print a single causal per-time block
#' @keywords internal
.print_causal_block <- function(r, se = NULL) {
  pad <- 52  # label width for alignment

  # Effects on means (delta_beta) — single component, SE available
  .print_line("Effects on means (delta_beta)", r$Delta_beta,
              se = if (!is.null(se)) se$se_Delta_beta, pad = pad)

  # Effects on SDs (delta_lambda) — single component, SE available
  .print_line("Effects on SDs (delta_lambda)", r$Delta_lambda,
              se = if (!is.null(se)) se$se_Delta_lambda, pad = pad)

  # Distribution of treatment (delta_pi) — total on header line, B/W below
  .print_line("Distribution of treatment (delta_pi)", r$Delta_pi,
              se = if (!is.null(se)) se$se_Delta_compositional, pad = pad)
  .print_line("Between-group", r$Delta_pi_B, indent = 6, pad = pad)
  .print_line("Within-group", r$Delta_pi_W, indent = 6, pad = pad)

  # Pre-treatment inequality (delta_pre) — total on header line, sub-components below
  .print_line("Pre-treatment inequality (delta_pre)", r$Delta_pre,
              se = if (!is.null(se)) se$se_Delta_pretreatment, pad = pad)
  .print_line("Means", r$Delta_mu, indent = 6, pad = pad)
  .print_line("Variances", r$Delta_sigma, indent = 6, pad = pad)

  # Grand total
  .print_line("Total", r$Delta_total,
              se = if (!is.null(se)) se$se_Delta_total, pad = pad)
}

#' Print a single descriptive per-time block
#' @keywords internal
.print_desc_block <- function(row, ystat = "Var") {
  pad <- 44  # label width for alignment

  # Between-group (delta_mu): changing group means
  # For Var: purely between-group (VarW doesn't depend on mu)
  # For CV2: has both B/W sub-components (CV2W depends on mu via grand mean denominator)
  if (ystat == "CV2") {
    .print_line("Between-group (delta_mu)", row$delta_mu, pad = pad)
    .print_line("Between-group", row$delta_mu_B, indent = 6, pad = pad)
    .print_line("Within-group", row$delta_mu_W, indent = 6, pad = pad)
  } else {
    .print_line("Between-group (delta_mu)", row$delta_mu, pad = pad)
  }

  # Within-group (delta_sigma): changing group dispersions
  # Always purely within-group (neither VarB nor CV2B depends on sigma)
  .print_line("Within-group (delta_sigma)", row$delta_sigma, pad = pad)

  # Compositional (delta_pi): changing group sizes
  .print_line("Compositional (delta_pi)", row$delta_pi, pad = pad)
  .print_line("Between-group", row$delta_pi_B, indent = 6, pad = pad)
  .print_line("Within-group", row$delta_pi_W, indent = 6, pad = pad)

  # Grand total
  .print_line("Total", row$delta_T, pad = pad)
}

#' Print a single line: label + value + optional SE
#' @keywords internal
.print_line <- function(label, value, se = NULL, indent = 4, pad = 52) {
  spaces <- paste(rep(" ", indent), collapse = "")
  label_with_colon <- paste0(label, ":")
  if (!is.null(se)) {
    cat(sprintf("%s%-*s %10.4f  (SE = %.4f)\n", spaces, pad - indent, label_with_colon, value, se))
  } else {
    cat(sprintf("%s%-*s %10.4f\n", spaces, pad - indent, label_with_colon, value))
  }
}

#' Print a data.frame with indentation
#' @keywords internal
.print_indented_df <- function(df, indent = 4, digits = 4) {
  txt <- capture.output(print(df, row.names = FALSE, digits = digits))
  spaces <- paste(rep(" ", indent), collapse = "")
  cat(paste0(spaces, txt, "\n"), sep = "")
}

#' Get sorted list of all times (ref + result times)
#' @keywords internal
.get_all_times <- function(results, ref) {
  result_times <- as.numeric(names(results))
  sort(unique(c(ref, result_times)))
}

#' Build and print cross-sectional sub-component table from longitudinal results
#' @keywords internal
.build_longit_subcomponents <- function(results, ref, cross_se = NULL,
                                         lambda_tests = NULL, beta_tests = NULL) {
  times <- names(results)
  first_result <- results[[1]]

  # Collect all time periods (ref + result times) in order
  all_times <- sort(unique(c(ref, sapply(results, function(r) r$time))))

  # Extract components for each time
  comps_list <- list()
  # Reference: use components_t0 from first result
  c0 <- first_result$components_t0
  comps_list[[as.character(ref)]] <- c0
  for (t in times) {
    comps_list[[t]] <- results[[t]]$components_t
  }

  # Column widths
  w <- 10  # width for numeric columns
  has_se <- !is.null(cross_se)
  has_wald <- !is.null(lambda_tests)
  has_beta <- !is.null(beta_tests)

  # Header
  hdr <- sprintf(" %4s %*s %*s %*s %*s %*s %*s %*s",
                  "time", w, "Var(beta)", w + 4, "2Cov(mu,beta)", w, "tau_B",
                  w + 4, "E(sigma2)*E(f)", w + 4, "Cov(sigma2,f)", w, "tau_W", w, "tau_T")
  if (has_wald) hdr <- paste0(hdr, sprintf(" %*s", w, "lambda=0"))
  if (has_beta) hdr <- paste0(hdr, sprintf(" %*s", w, "beta=hom"))
  cat(hdr, "\n")

  for (t in all_times) {
    tc <- as.character(t)
    ct <- comps_list[[tc]]

    vals <- c(
      ct$Var_pi_beta,
      2 * ct$Cov_pi_mu_beta,
      ct$Var_pi_beta + 2 * ct$Cov_pi_mu_beta,
      ct$mean_sigma2_0 * ct$mean_f,
      ct$Cov_pi_sigma2_f,
      ct$mean_sigma2_0 * ct$mean_f + ct$Cov_pi_sigma2_f,
      ct$Var_pi_beta + 2 * ct$Cov_pi_mu_beta +
        ct$mean_sigma2_0 * ct$mean_f + ct$Cov_pi_sigma2_f
    )

    # Value row
    val_line <- sprintf(" %4s %*s %*s %*s %*s %*s %*s %*s",
                        tc,
                        w, formatC(vals[1], format = "f", digits = 4, width = w),
                        w + 4, formatC(vals[2], format = "f", digits = 4, width = w + 4),
                        w, formatC(vals[3], format = "f", digits = 4, width = w),
                        w + 4, formatC(vals[4], format = "f", digits = 4, width = w + 4),
                        w + 4, formatC(vals[5], format = "f", digits = 4, width = w + 4),
                        w, formatC(vals[6], format = "f", digits = 4, width = w),
                        w, formatC(vals[7], format = "f", digits = 4, width = w))
    if (has_wald) {
      lt <- lambda_tests[[tc]]
      if (!is.null(lt)) {
        p_fmt <- if (lt$p_value < 0.001) "<0.001" else formatC(lt$p_value, format = "f", digits = 4)
        val_line <- paste0(val_line, sprintf(" %*s", w, p_fmt))
      } else {
        val_line <- paste0(val_line, sprintf(" %*s", w, ""))
      }
    }
    if (has_beta) {
      bt <- beta_tests[[tc]]
      if (!is.null(bt)) {
        p_fmt <- if (bt$p_value < 0.001) "<0.001" else formatC(bt$p_value, format = "f", digits = 4)
        val_line <- paste0(val_line, sprintf(" %*s", w, p_fmt))
      } else {
        val_line <- paste0(val_line, sprintf(" %*s", w, ""))
      }
    }
    cat(val_line, "\n")

    # SE row (if available)
    if (has_se && tc %in% names(cross_se)) {
      se <- cross_se[[tc]]
      se_vals <- c(
        se$se_sub$se_het_B,
        se$se_sub$se_cov_B,
        se$se_tau_B,
        se$se_sub$se_het_W,
        se$se_sub$se_cov_W,
        se$se_tau_W,
        se$se_tau_total
      )
      # Format as (X.XXXX)
      se_strs <- vapply(se_vals, function(v) {
        sprintf("(%s)", formatC(v, format = "f", digits = 4))
      }, character(1))

      se_line <- sprintf(" %4s %*s %*s %*s %*s %*s %*s %*s",
                         "SE",
                         w, se_strs[1], w + 4, se_strs[2], w, se_strs[3],
                         w + 4, se_strs[4], w + 4, se_strs[5],
                         w, se_strs[6], w, se_strs[7])
      if (has_wald) se_line <- paste0(se_line, sprintf(" %*s", w, ""))
      if (has_beta) se_line <- paste0(se_line, sprintf(" %*s", w, ""))
      cat(se_line, "\n")
    }
  }
}

#' Print estimate with SE (legacy, used by cross-sectional)
#' @keywords internal
.print_val_with_se <- function(label, estimate, se) {
  cat(sprintf("    %-28s %10.4f  (SE = %.4f)\n", paste0(label, ":"), estimate, se))
}
