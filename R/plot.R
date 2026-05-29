# ============================================================================ #
# Plot methods for ineqx result objects
# ============================================================================ #

#' @importFrom ggplot2 .data
NULL


# ---------------------------------------------------------------------------- #
# Internal: Resolve ci argument to a standard form
# ---------------------------------------------------------------------------- #

#' @keywords internal
.resolve_ci <- function(ci, x) {
  if (isFALSE(ci) || identical(ci, "none")) {
    return(list(method = "none", boot = NULL))
  }
  if (isTRUE(ci)) {
    return(list(method = "delta", boot = NULL))
  }
  if (identical(ci, "delta")) {
    return(list(method = "delta", boot = NULL))
  }
  if (identical(ci, "boot")) {
    return(list(method = "boot", boot = boot_config()))
  }
  if (inherits(ci, "ineqx_boot_config")) {
    return(list(method = "boot", boot = ci))
  }
  stop("'ci' must be FALSE, 'none', TRUE, 'delta', 'boot', or boot_config().")
}

# ---------------------------------------------------------------------------- #
# Internal: Resolve ci for params plots (delta only; boot falls back to delta)
# ---------------------------------------------------------------------------- #

#' @keywords internal
.resolve_ci_params <- function(ci) {
  if (is.list(ci)) {
    if (ci$method == "boot") {
      warning("Bootstrap CIs not available for params plots. ",
              "Run ineqx() with se = boot_config() first. ",
              "Falling back to delta-method CIs.", call. = FALSE)
    }
    ci$method != "none"
  } else {
    isTRUE(ci)
  }
}

# ---------------------------------------------------------------------------- #
# Internal: Normalize the `stats` argument for type = "wibe"
# ---------------------------------------------------------------------------- #
#
# Shared by .plot_causal_cross_bar and .plot_causal_longit_wibe so both methods
# accept the same vocabulary. Expands shorthands and validates.
#
# Canonical: c("tau", "tau_b", "tau_w",
#              "het_b", "cov_b", "rescale_b",
#              "het_w", "cov_w", "rescale_w")
# Shorthands:
#   "het_cov"   -> c("het_b", "cov_b", "het_w", "cov_w")           [no rescale]
#   "het_cov_b" -> c("het_b", "cov_b")
#   "het_cov_w" -> c("het_w", "cov_w")
#   "subs"      -> c("het_b","cov_b","rescale_b","het_w","cov_w","rescale_w")
#   "subs_b"    -> c("het_b","cov_b","rescale_b")
#   "subs_w"    -> c("het_w","cov_w","rescale_w")

#' @keywords internal
.normalize_wibe_stats <- function(stats,
                                   default = c("tau", "tau_b", "tau_w")) {
  if (is.null(stats)) return(default)
  expansions <- list(
    het_cov   = c("het_b", "cov_b", "het_w", "cov_w"),
    het_cov_b = c("het_b", "cov_b"),
    het_cov_w = c("het_w", "cov_w"),
    subs      = c("het_b", "cov_b", "rescale_b",
                  "het_w", "cov_w", "rescale_w"),
    subs_b    = c("het_b", "cov_b", "rescale_b"),
    subs_w    = c("het_w", "cov_w", "rescale_w")
  )
  out <- unlist(lapply(stats, function(s) {
    if (s %in% names(expansions)) expansions[[s]] else s
  }))
  canonical <- c("tau", "tau_b", "tau_w",
                 "het_b", "cov_b", "rescale_b",
                 "het_w", "cov_w", "rescale_w")
  bad <- setdiff(out, canonical)
  if (length(bad) > 0) {
    stop("Invalid 'stats' values: ", paste(bad, collapse = ", "),
         ". Valid options: ",
         paste(c(canonical, names(expansions)), collapse = ", "),
         call. = FALSE)
  }
  unique(out)
}

# ---------------------------------------------------------------------------- #
# Internal: Decomposition component registry
# ---------------------------------------------------------------------------- #

#' @keywords internal
.decomp_registry <- function() {
  list(
    behavioral    = list(label = "Behavioral",
                         field = "Delta_behavioral",
                         se_field = "se_Delta_behavioral",
                         color = "#008b94"),
    compositional = list(label = "Compositional",
                         field = "Delta_compositional",
                         se_field = "se_Delta_compositional",
                         color = "#7b2d8e"),
    pretreatment  = list(label = "Pre-treatment",
                         field = "Delta_pretreatment",
                         se_field = "se_Delta_pretreatment",
                         color = "#e85d00"),
    total         = list(label = "Total",
                         field = "Delta_total",
                         se_field = "se_Delta_total",
                         color = "black"),
    delta_beta    = list(label = "\u0394\u03b2",
                         field = "Delta_beta",
                         se_field = "se_Delta_beta",
                         color = "#008b94"),
    delta_beta_b  = list(label = "\u0394\u03b2_b",
                         field = "Delta_beta_B",
                         se_field = "se_Delta_beta_B",
                         color = "#008b94"),
    delta_beta_w  = list(label = "\u0394\u03b2_w",
                         field = "Delta_beta_W",
                         se_field = "se_Delta_beta_W",
                         color = "#5cc4cb"),
    delta_lambda  = list(label = "\u0394\u03bb",
                         field = "Delta_lambda",
                         se_field = "se_Delta_lambda",
                         color = "#cc6677"),
    delta_pi      = list(label = "\u0394\u03c0",
                         field = "Delta_compositional",
                         se_field = "se_Delta_compositional",
                         color = "#7b2d8e"),
    delta_pi_b    = list(label = "\u0394\u03c0_b",
                         field = "Delta_pi_B",
                         se_field = "se_Delta_pi_B",
                         color = "#7b2d8e"),
    delta_pi_w    = list(label = "\u0394\u03c0_w",
                         field = "Delta_pi_W",
                         se_field = "se_Delta_pi_W",
                         color = "#b47cc5"),
    delta_mu      = list(label = "\u0394\u03bc",
                         field = "Delta_mu",
                         se_field = "se_Delta_mu",
                         color = "#e85d00"),
    delta_mu_b    = list(label = "\u0394\u03bc_b",
                         field = "Delta_mu_B",
                         se_field = "se_Delta_mu_B",
                         color = "#e85d00"),
    delta_mu_w    = list(label = "\u0394\u03bc_w",
                         field = "Delta_mu_W",
                         se_field = "se_Delta_mu_W",
                         color = "#ffae66"),
    delta_sigma   = list(label = "\u0394\u03c3",
                         field = "Delta_sigma",
                         se_field = "se_Delta_sigma",
                         color = "#f5a623")
  )
}

# ---------------------------------------------------------------------------- #
# plot.ineqx_desc
# ---------------------------------------------------------------------------- #

#' Plot descriptive variance decomposition
#'
#' @param x An \code{ineqx_desc} object
#' @param type Character, plot type:
#'   \code{"wibe"} for within/between levels over time,
#'   \code{"decomp"} for change contributions over time from
#'   delta-mu, delta-sigma, and delta-pi
#'   relative to the reference period (requires ref),
#'   \code{"params"} for group-level means and SDs over time,
#'   \code{"ineq"} for overall inequality statistics over time,
#'   \code{"ineq.group"} for per-group inequality statistics over time.
#' @param stats For \code{type = "ineq"} or \code{"ineq.group"}: a list of
#'   inequality statistics to compute. Can be character strings
#'   (\code{"V"}, \code{"VL"}, \code{"CV2"}, \code{"Gini"}, \code{"Theil"})
#'   or custom functions with signature \code{f(y, w = NULL)}.
#'   Default: \code{c("V")}.
#' @param ci How to compute confidence intervals. Options:
#'   \code{FALSE} or \code{"none"} (default): no CIs.
#'   \code{TRUE} or \code{"delta"}: delta method CIs (reads stored SEs from
#'   \code{ineqx()}).
#'   \code{"boot"}: bootstrap CIs with default settings.
#'   \code{boot_config()}: bootstrap CIs with custom settings.
#' @param style Character: \code{"line"} (default) for connected lines with
#'   ribbon CIs, or \code{"point"} for points with error bar CIs.
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot2 object
#' @export
plot.ineqx_desc <- function(x, type = "decomp", stats = NULL,
                            ci = FALSE, style = "line", ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  style <- match.arg(style, c("line", "point"))

  if (type == "wibe") {
    .plot_desc_wibe(x, ci = ci, style = style)
  } else if (type == "decomp") {
    if (is.null(x$deltas)) {
      stop("Decomposition not available. Set 'ref' in ineqx() to compute deltas.")
    }
    .plot_desc_deltas(x, ci = ci, style = style)
  } else if (type == "params") {
    .plot_desc_params(x, ci = ci, style = style)
  } else if (type == "ineq") {
    if (is.null(stats)) stats <- list("V")
    .plot_desc_ineq(x, stats, ci = ci, style = style)
  } else if (type == "ineq.group") {
    if (is.null(stats)) stats <- list("V")
    .plot_desc_ineq_group(x, stats, ci = ci, style = style)
  } else {
    stop("Unknown plot type '", type,
         "'. Use 'wibe', 'decomp', 'params', 'ineq', or 'ineq.group'.")
  }
}

#' @keywords internal
.plot_desc_wibe <- function(x, ci = FALSE, style = "line") {
  totals <- x$totals
  if (!"time" %in% names(totals)) {
    stop("Time variable required for plotting")
  }

  if (x$ystat %in% c("Var", "VL")) {
    plot_df <- data.frame(
      time = rep(totals$time, 3),
      Component = rep(c("Within", "Between", "Total"), each = nrow(totals)),
      value = c(totals$VarW, totals$VarB, totals$VarT)
    )
  } else {
    plot_df <- data.frame(
      time = rep(totals$time, 3),
      Component = rep(c("Within", "Between", "Total"), each = nrow(totals)),
      value = c(totals$CV2W, totals$CV2B, totals$CV2T)
    )
  }

  plot_df$Component <- factor(plot_df$Component,
                              levels = c("Total", "Within", "Between"))

  # CIs
  ci_info <- .resolve_ci(ci, x)
  if (ci_info$method %in% c("stored", "delta") && !is.null(x[["se"]])) {
    # Read stored delta method SEs
    se_prefix <- if (x$ystat %in% c("Var", "VL")) "se_Var" else "se_CV2"
    se_fields <- paste0(se_prefix, c("W", "B", "T"))
    se_vec <- numeric(nrow(plot_df))
    for (i in seq_along(totals$time)) {
      t_key <- as.character(totals$time[i])
      se_t <- x$se$totals[[t_key]]
      if (!is.null(se_t)) {
        for (j in seq_along(se_fields)) {
          row_idx <- (j - 1) * nrow(totals) + i
          se_val <- se_t[[se_fields[j]]]
          if (!is.null(se_val)) se_vec[row_idx] <- se_val
        }
      }
    }
    plot_df$ymin <- plot_df$value - 1.96 * se_vec
    plot_df$ymax <- plot_df$value + 1.96 * se_vec
  } else if (ci_info$method == "boot" && !is.null(x$raw_data)) {
    B <- ci_info$boot$B
    ystat <- if (x$ystat == "VL") "Var" else x$ystat
    ci_df <- .bootstrap_desc_ci(x$raw_data, B = B, compute_fn = function(d) {
      wibe <- .compute_wibe_from_d(d)
      totals <- .compute_totals_n(wibe, sort(unique(d$time)))
      if (ystat %in% c("Var", "VL")) {
        data.frame(time = rep(totals$time, 3),
                   Component = rep(c("Within", "Between", "Total"),
                                   each = nrow(totals)),
                   value = c(totals$VarW, totals$VarB, totals$VarT),
                   stringsAsFactors = FALSE)
      } else {
        data.frame(time = rep(totals$time, 3),
                   Component = rep(c("Within", "Between", "Total"),
                                   each = nrow(totals)),
                   value = c(totals$CV2W, totals$CV2B, totals$CV2T),
                   stringsAsFactors = FALSE)
      }
    })
    plot_df <- merge(plot_df, ci_df, by = c("time", "Component"))
  }

  colors <- c("Total" = "black", "Within" = "#008b94", "Between" = "#e85d00")
  fill_colors <- c("Total" = "grey50", "Within" = "#008b94", "Between" = "#e85d00")

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              color = .data$Component,
                                              group = .data$Component))

  has_ci <- "ymin" %in% names(plot_df)
  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos)
    }
    p <- p + ggplot2::geom_point(size = 2, position = pos$pos)
  } else {
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$Component),
        alpha = 0.15, color = NA) +
        ggplot2::scale_fill_manual(values = fill_colors)
    }
    p <- p + ggplot2::geom_line(linewidth = 1)
  }

  p + ggplot2::labs(x = "Time", y = x$ystat,
                    title = "Within/Between Decomposition") +
    ggplot2::scale_color_manual(values = colors) +
    theme_ineqx
}

#' @keywords internal
.plot_desc_deltas <- function(x, ci = FALSE, style = "line") {
  deltas <- x$deltas

  # Plot per-parameter total deltas (mu, sigma, pi contributions)
  plot_df <- data.frame(
    time = rep(deltas$time, 4),
    Component = rep(c("\u0394\u03bc", "\u0394\u03c3",
                      "\u0394\u03c0", "Total"),
                    each = nrow(deltas)),
    value = c(deltas$delta_mu, deltas$delta_sigma,
              deltas$delta_pi, deltas$delta_T)
  )

  plot_df$Component <- factor(plot_df$Component,
                              levels = c("Total", "\u0394\u03bc",
                                         "\u0394\u03c3", "\u0394\u03c0"))

  # CIs
  ci_info <- .resolve_ci(ci, x)
  if (ci_info$method == "boot" && !is.null(x$raw_data) && !is.null(x$ref)) {
    B <- ci_info$boot$B
    ref <- x$ref
    ystat <- if (x$ystat == "VL") "Var" else x$ystat
    order <- x$order
    ci_df <- .bootstrap_desc_ci(x$raw_data, B = B, compute_fn = function(d) {
      wibe <- .compute_wibe_from_d(d)
      time_levels <- sort(unique(d$time))
      dd <- .compute_desc_deltas(wibe, time_levels, ref, ystat, order)
      data.frame(
        time = rep(dd$time, 4),
        Component = rep(c("\u0394\u03bc", "\u0394\u03c3",
                          "\u0394\u03c0", "Total"),
                        each = nrow(dd)),
        value = c(dd$delta_mu, dd$delta_sigma, dd$delta_pi, dd$delta_T),
        stringsAsFactors = FALSE)
    })
    plot_df <- merge(plot_df, ci_df, by = c("time", "Component"))
  } else if (ci_info$method %in% c("stored", "delta") && !is.null(x[["se"]]$deltas)) {
    # Read stored delta method SEs
    se_names <- c(
      "\u0394\u03bc" = "se_delta_mu", "\u0394\u03c3" = "se_delta_sigma",
      "\u0394\u03c0" = "se_delta_pi", "Total" = "se_delta_T"
    )
    se_vec <- numeric(nrow(plot_df))
    has_real <- FALSE
    for (i in seq_len(nrow(plot_df))) {
      t_key <- as.character(plot_df$time[i])
      comp  <- as.character(plot_df$Component[i])
      se_val <- x[["se"]]$deltas[[t_key]][[se_names[comp]]]
      if (!is.null(se_val) && !is.na(se_val)) {
        se_vec[i] <- se_val
        has_real <- TRUE
      }
    }
    if (has_real) {
      plot_df$ymin <- plot_df$value - 1.96 * se_vec
      plot_df$ymax <- plot_df$value + 1.96 * se_vec
    } else {
      warning("Delta method CIs for decomposition components are not available. ",
              "Re-run ineqx() with se = 'delta', or use ci = 'boot'.",
              call. = FALSE)
    }
  } else if (ci_info$method != "none") {
    warning("No stored standard errors found. ",
            "Run ineqx() with se = 'delta' or use ci = 'boot'.",
            call. = FALSE)
  }

  colors <- c("Total" = "black", "\u0394\u03bc" = "#e85d00",
              "\u0394\u03c3" = "#008b94", "\u0394\u03c0" = "#7b2d8e")
  fill_colors <- colors
  fill_colors["Total"] <- "grey50"

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              color = .data$Component,
                                              group = .data$Component))

  has_ci <- "ymin" %in% names(plot_df)
  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos)
    }
    p <- p + ggplot2::geom_point(size = 2, position = pos$pos)
  } else {
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$Component),
        alpha = 0.15, color = NA) +
        ggplot2::scale_fill_manual(values = fill_colors)
    }
    p <- p + ggplot2::geom_line(linewidth = 1)
  }

  p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Time",
                  y = paste0("Change in ", x$ystat, " (from ref = ", x$ref, ")"),
                  title = "Descriptive Decomposition of Change") +
    ggplot2::scale_color_manual(values = colors) +
    theme_ineqx
}

# ---------------------------------------------------------------------------- #
# plot.ineqx_causal_cross
# ---------------------------------------------------------------------------- #

#' Plot cross-sectional causal decomposition
#'
#' @param x An \code{ineqx_causal_cross} object
#' @param type Character, the plot type:
#'   \describe{
#'     \item{\code{"wibe"}}{Within/between treatment effects (bar chart).}
#'     \item{\code{"wibe.group"}}{Group-level contributions to within/between.}
#'     \item{\code{"treat"}}{Predicted treatment effect distributions by group.}
#'     \item{\code{"treat.params"}}{Treatment effect parameters (beta, lambda)
#'       by group (bar chart).}
#'     \item{\code{"outcome"}}{Predicted outcome distributions (control vs
#'       treated) by group.}
#'     \item{\code{"outcome.params"}}{Predicted means and SDs under control vs
#'       treatment. For simple-difference models, a bar chart per group. For
#'       DiD models, a netted-out line chart anchored at the pre-period
#'       treated level: the post-period gap between observed and
#'       counterfactual lines equals the DiD ATT.}
#'   }
#' @param ci Whether to show confidence intervals. Accepts \code{FALSE} or
#'   \code{"none"} (no CIs), \code{TRUE} or \code{"delta"} (use stored SEs),
#'   \code{"boot"} (bootstrap with defaults), or \code{boot_config()} (bootstrap
#'   with custom settings). Default \code{FALSE}.
#' @param stats Character vector. For \code{type = "wibe"}: components to
#'   display. The bar plot recognizes the same canonical vocabulary as the
#'   longitudinal \code{type = "wibe"} method:
#'   \code{"tau"} (Total τ bar), \code{"tau_b"} (Between τ bar),
#'   \code{"tau_w"} (Within τ bar), \code{"het_b"}, \code{"cov_b"},
#'   \code{"rescale_b"} (Between sub-components: heterogeneity, sorting,
#'   rescaling), and the within mirrors \code{"het_w"}, \code{"cov_w"},
#'   \code{"rescale_w"}. Single-segment bars are allowed. The legend labels
#'   \code{"cov_*"} as ``Sorting'' and \code{"rescale_*"} as ``Rescaling''.
#'   Shorthands: \code{"het_cov"} = the four heterogeneity / sorting sub-components
#'   (no rescaling),
#'   \code{"het_cov_b"} = \code{c("het_b","cov_b")},
#'   \code{"het_cov_w"} = \code{c("het_w","cov_w")},
#'   \code{"subs"} = all six sub-components (het + cov + rescale, between and within),
#'   \code{"subs_b"} = \code{c("het_b","cov_b","rescale_b")},
#'   \code{"subs_w"} = \code{c("het_w","cov_w","rescale_w")}.
#'   When \code{tau} and any sub-component appear at the same column, the bars
#'   are dodged. Default \code{"tau"} preserves the classic three-bar
#'   (Total/Between/Within) layout; \code{c("tau","het_cov")} preserves the
#'   dodged classic layout.
#' @param trim Numeric between 0 and 1. For distribution plots
#'   (\code{"treat"}, \code{"outcome"}): quantile at which to trim tails.
#'   Default \code{0.995} (trims 0.5\% from each tail).
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot2 object
#' @export
plot.ineqx_causal_cross <- function(x, type = "wibe", ci = FALSE,
                                     stats = NULL, trim = 0.995, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  # Resolve ci: accept "none", "delta", TRUE/FALSE, boot_config()
  ci_info <- .resolve_ci(ci, x)
  show_ci   <- ci_info$method != "none" && !is.null(x[["se"]])
  ci_params <- ci_info

  if (type == "wibe") {
    if (is.null(stats)) stats <- "tau"
    .plot_causal_cross_bar(x, ci = show_ci, stats = stats)
  } else if (type == "wibe.group") {
    bg <- x$by_group
    plot_df <- data.frame(
      group = factor(rep(bg$group, 2)),
      Component = rep(c("Within", "Between"), each = nrow(bg)),
      value = c(bg$contrib_W, bg$contrib_B)
    )

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$group, y = .data$value,
                                           fill = .data$Component)) +
      ggplot2::geom_col(position = "dodge", width = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "Group", y = paste0("Contribution to ", x$ystat),
                    title = "Group-Level Contributions") +
      ggplot2::scale_fill_manual(values = c("Within" = "#008b94",
                                             "Between" = "#e85d00")) +
      theme_ineqx
  } else if (type == "treat") {
    .plot_causal_cross_dist(x, show = "treat", trim = trim)
  } else if (type == "treat.params") {
    .plot_causal_cross_treat(x, ci = ci_params)
  } else if (type == "outcome") {
    .plot_causal_cross_dist(x, show = "outcome", trim = trim)
  } else if (type == "outcome.params") {
    .plot_causal_cross_outcome(x, ci = ci_params)
  } else if (type == "pretrends") {
    stop("Plot type 'pretrends' requires a longitudinal model with multiple ",
         "pre-period observations across time. For a single cross-section, ",
         "the parallel-trends assumption cannot be visualized.")
  } else {
    stop("Unknown plot type '", type,
         "'. Use 'wibe', 'wibe.group', 'treat', 'treat.params', ",
         "'outcome', or 'outcome.params'.")
  }
}

# ---------------------------------------------------------------------------- #
# plot.ineqx_causal_longit
# ---------------------------------------------------------------------------- #

#' Plot longitudinal causal decomposition
#'
#' @param x An \code{ineqx_causal_longit} object (including Shapley-averaged results)
#' @param type Character, the plot type:
#'   \describe{
#'     \item{\code{"decomp"}}{Decomposition of changes relative to reference.
#'       Use \code{stats} to select components (default: 3 aggregate + total).}
#'     \item{\code{"wibe"}}{Cross-sectional treatment effects on within/between
#'       inequality at each time (levels, not changes).}
#'     \item{\code{"treat"}}{Predicted treatment effect distributions.
#'       When \code{time} is specified, shows per-group distributions.
#'       When \code{time = NULL}, shows pi-weighted marginals across all times.}
#'     \item{\code{"treat.params"}}{Treatment effect parameters (beta/lambda)
#'       over time (line chart).}
#'     \item{\code{"outcome"}}{Predicted outcome distributions (control vs
#'       treated). When \code{time} is specified, shows per-group distributions.
#'       When \code{time = NULL}, shows pi-weighted marginals.}
#'     \item{\code{"outcome.params"}}{Predicted means and SDs under control vs
#'       treatment over time (line chart). For DiD models, the "Control" line
#'       represents the DiD-implied counterfactual untreated level for the
#'       treated post-period subpopulation, so the gap between the lines
#'       equals the DiD ATT.}
#'     \item{\code{"pretrends"}}{Pre-period predicted Treated and Control
#'       levels over time. Available only for DiD models. Under parallel
#'       trends the two lines should evolve with the same slope; their
#'       (constant) vertical gap reflects pre-existing selection. Diverging
#'       slopes indicate a violation of parallel trends.}
#'     \item{\code{"shapley"}}{Shapley averages with ordering ranges
#'       (only when \code{order = "shapley"}). Use \code{stats} to select
#'       components.}
#'   }
#' @param ci Whether to show confidence intervals. Accepts \code{FALSE} or
#'   \code{"none"} (no CIs), \code{TRUE} or \code{"delta"} (use stored SEs),
#'   \code{"boot"} (bootstrap with defaults), or \code{boot_config()} (bootstrap
#'   with custom settings). Default \code{FALSE}.
#' @param style Character: \code{"line"} (default) for connected lines with
#'   ribbon CIs, or \code{"point"} for points with error bar CIs.
#' @param stats Character vector of components to display. Meaning depends on
#'   \code{type}:
#'   \describe{
#'     \item{For \code{"decomp"} and \code{"shapley"}:}{Any combination of
#'       aggregate components (\code{"behavioral"}, \code{"compositional"},
#'       \code{"pretreatment"}, \code{"total"}) and/or individual deltas
#'       (\code{"delta_beta"}, \code{"delta_lambda"}, \code{"delta_pi"},
#'       \code{"delta_pi_b"}, \code{"delta_pi_w"}, \code{"delta_mu"},
#'       \code{"delta_sigma"}).
#'       Default for \code{"decomp"}: \code{c("behavioral", "compositional",
#'       "pretreatment", "total")}.
#'       Default for \code{"shapley"}: \code{c("behavioral", "compositional",
#'       "pretreatment")}.}
#'     \item{For \code{"wibe"}:}{Components to display,
#'       default \code{c("tau", "tau_b", "tau_w")}. Canonical options:
#'       \code{"tau"} (total), \code{"tau_b"}, \code{"tau_w"},
#'       \code{"het_b"}, \code{"cov_b"}, \code{"rescale_b"},
#'       \code{"het_w"}, \code{"cov_w"}, \code{"rescale_w"}.
#'       Legend labels \code{"cov_*"} as ``Sorting'' and \code{"rescale_*"} as
#'       ``Rescaling'' (rescaling values are zero under \code{ystat = "Var"}).
#'       Shorthands: \code{"het_cov"} = the four het/cov sub-components
#'       (no rescaling),
#'       \code{"het_cov_b"} = \code{c("het_b","cov_b")},
#'       \code{"het_cov_w"} = \code{c("het_w","cov_w")},
#'       \code{"subs"} = all six sub-components,
#'       \code{"subs_b"} = \code{c("het_b","cov_b","rescale_b")},
#'       \code{"subs_w"} = \code{c("het_w","cov_w","rescale_w")}.}
#'   }
#' @param time Numeric, time point for \code{type = "treat"} or
#'   \code{type = "outcome"}. If \code{NULL} (default), pi-weighted marginal
#'   distributions are shown across all time points.
#' @param trim Numeric between 0 and 1. For distribution plots
#'   (\code{"treat"}, \code{"outcome"}): quantile at which to trim tails.
#'   Default \code{0.995}.
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot2 object
#' @export
plot.ineqx_causal_longit <- function(x, type = "decomp", ci = FALSE,
                                      style = "line", stats = NULL,
                                      time = NULL, trim = 0.995, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  style <- match.arg(style, c("line", "point"))

  # Resolve ci: accept "none", "delta", TRUE/FALSE, boot_config()
  ci_info <- .resolve_ci(ci, x)
  # For decomp/wibe/shapley: CIs come from stored decomposition SEs (x$se)
  show_ci       <- ci_info$method != "none" && !is.null(x[["se"]])
  # For params plots: pass the full ci_info so the helper can choose delta vs boot
  ci_params <- ci_info

  is_shapley <- identical(x$order, "shapley")

  if (type == "decomp") {
    if (is.null(stats)) {
      stats <- c("behavioral", "compositional", "pretreatment", "total")
    }
    return(.plot_causal_longit_decomp(x, stats = stats, ci = show_ci,
                                       style = style))
  } else if (type == "wibe") {
    return(.plot_causal_longit_wibe(x, ci = show_ci, style = style, stats = stats))

  } else if (type == "shapley") {
    if (!is_shapley) {
      stop("Plot type 'shapley' requires order = 'shapley'. ",
           "Re-run ineqx() with order = 'shapley'.")
    }
    if (is.null(stats)) {
      stats <- c("behavioral", "compositional", "pretreatment")
    }
    return(.plot_shapley_ranges(x, stats = stats, style = style))
  } else if (type == "treat") {
    return(.plot_causal_longit_dist(x, time = time, show = "treat", trim = trim))
  } else if (type == "treat.params") {
    return(.plot_causal_longit_treat(x, ci = ci_params, style = style))
  } else if (type == "outcome") {
    return(.plot_causal_longit_dist(x, time = time, show = "outcome", trim = trim))
  } else if (type == "outcome.params") {
    return(.plot_causal_longit_outcome(x, ci = ci_params, style = style))
  } else if (type == "pretrends") {
    if (!isTRUE(x$params$is_did)) {
      stop("Plot type 'pretrends' requires a DiD model (post argument). ",
           "Refit with ineqx(..., post = '<your post column>').")
    }
    return(.plot_causal_longit_pretrends(x, ci = ci_params, style = style))
  } else {
    stop("Unknown plot type '", type,
         "'. Use 'decomp', 'wibe', 'shapley', 'treat', ",
         "'treat.params', 'outcome', 'outcome.params', or 'pretrends'.")
  }
}

# ---------------------------------------------------------------------------- #
# Unified decomp helper (replaces old decomp + decomp.six)
# ---------------------------------------------------------------------------- #

#' @keywords internal
.plot_causal_longit_decomp <- function(x, stats, ci = FALSE, style = "line") {
  registry <- .decomp_registry()
  invalid <- setdiff(stats, names(registry))
  if (length(invalid) > 0) {
    stop("Unknown stats values: ", paste(invalid, collapse = ", "),
         ". Valid: ", paste(names(registry), collapse = ", "))
  }

  is_shapley <- identical(x$order, "shapley")

  # Build data from results
  times <- as.numeric(names(x$results))
  ref_time <- x$ref
  if (!is.null(ref_time) && !ref_time %in% times) {
    times <- sort(c(ref_time, times))
  }

  # Extract values for each requested component
  comp_info <- registry[stats]
  labels <- vapply(comp_info, `[[`, character(1), "label")
  fields <- vapply(comp_info, `[[`, character(1), "field")
  se_fields <- vapply(comp_info, `[[`, character(1), "se_field")
  colors <- vapply(comp_info, `[[`, character(1), "color")
  names(colors) <- labels

  n_comp <- length(stats)
  vals <- matrix(0, nrow = length(times), ncol = n_comp)
  for (i in seq_along(times)) {
    t_char <- as.character(times[i])
    r <- x$results[[t_char]]
    if (!is.null(r)) {
      for (k in seq_len(n_comp)) {
        vals[i, k] <- r[[fields[k]]]
      }
    }
  }

  plot_df <- data.frame(
    time = rep(times, n_comp),
    Component = rep(labels, each = length(times)),
    value = as.vector(vals)
  )
  # Factor levels in user-specified order (reversed so first appears on top)
  plot_df$Component <- factor(plot_df$Component, levels = rev(labels))

  # Add CIs from existing $se
  if (ci) {
    se_vec <- numeric(nrow(plot_df))
    for (i in seq_len(n_comp)) {
      for (j in seq_along(times)) {
        t_key <- as.character(times[j])
        row_idx <- (i - 1) * length(times) + j
        se_val <- x$se[[t_key]][[se_fields[i]]]
        if (!is.null(se_val)) {
          se_vec[row_idx] <- se_val
        } else if (!is.null(ref_time) && times[j] == ref_time) {
          se_vec[row_idx] <- 0
        } else {
          se_vec[row_idx] <- NA_real_
        }
      }
    }
    plot_df$se <- se_vec
    plot_df$ymin <- plot_df$value - 1.96 * plot_df$se
    plot_df$ymax <- plot_df$value + 1.96 * plot_df$se
  }

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              color = .data$Component,
                                              group = .data$Component))

  has_ci <- "ymin" %in% names(plot_df)
  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos)
    }
    p <- p + ggplot2::geom_point(size = 2, position = pos$pos)
  } else {
    if (has_ci) {
      fill_colors <- colors
      if ("Total" %in% names(fill_colors)) fill_colors["Total"] <- "grey50"
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$Component),
        alpha = 0.15, color = NA) +
        ggplot2::scale_fill_manual(values = fill_colors)
    }
    p <- p + ggplot2::geom_line(linewidth = 1)
  }

  if (is_shapley) {
    y_label <- paste0("Shapley value (", x$ystat, ")")
    title <- "Longitudinal Causal Decomposition (Shapley)"
  } else {
    y_label <- paste0("Change in treatment effect on ", x$ystat)
    title <- "Longitudinal Causal Decomposition"
  }

  p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Time", y = y_label, title = title) +
    ggplot2::scale_color_manual(values = colors) +
    theme_ineqx
}


# ---------------------------------------------------------------------------- #
# type = "shapley": point + range across orderings
# ---------------------------------------------------------------------------- #

#' @keywords internal
.plot_shapley_ranges <- function(x, stats = c("behavioral", "compositional",
                                               "pretreatment"),
                                 style = "line") {
  registry <- .decomp_registry()
  invalid <- setdiff(stats, names(registry))
  if (length(invalid) > 0) {
    stop("Unknown stats values: ", paste(invalid, collapse = ", "),
         ". Valid: ", paste(names(registry), collapse = ", "))
  }

  s <- x$shapley
  r <- x$ranges
  times <- s$time

  # Add reference time point (all components = 0 at baseline)
  ref_time <- x$ref
  if (!is.null(ref_time) && !ref_time %in% times) {
    times <- sort(c(ref_time, times))
    ref_row_s <- as.data.frame(c(list(time = ref_time),
      setNames(as.list(rep(0, ncol(s) - 1)), setdiff(names(s), "time"))))
    s <- rbind(ref_row_s, s)
    ref_row_r <- as.data.frame(c(list(time = ref_time),
      setNames(as.list(rep(0, ncol(r) - 1)), setdiff(names(r), "time"))))
    r <- rbind(ref_row_r, r)
  }

  comp_info <- registry[stats]
  labels <- vapply(comp_info, `[[`, character(1), "label")
  fields <- vapply(comp_info, `[[`, character(1), "field")
  colors <- vapply(comp_info, `[[`, character(1), "color")
  names(colors) <- labels

  plot_df <- data.frame(
    time = rep(times, length(labels)),
    Component = rep(labels, each = length(times)),
    value = unlist(lapply(fields, function(f) s[[f]])),
    ymin = unlist(lapply(fields, function(f) r[[paste0(f, "_min")]])),
    ymax = unlist(lapply(fields, function(f) r[[paste0(f, "_max")]]))
  )
  plot_df$Component <- factor(plot_df$Component, levels = labels)

  fill_colors <- colors

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              color = .data$Component,
                                              group = .data$Component))

  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos) +
      ggplot2::geom_point(size = 2, position = pos$pos)
  } else {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$Component),
        alpha = 0.15, color = NA) +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2)
  }

  p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Time",
                  y = paste0("Shapley value (", x$ystat, ")"),
                  title = "Ordering Dependence (Shapley Range)") +
    ggplot2::scale_color_manual(values = colors) +
    theme_ineqx
}


# ---------------------------------------------------------------------------- #
# type = "wibe": cross-sectional treatment effects (levels) over time
# ---------------------------------------------------------------------------- #

#' @keywords internal
.plot_causal_longit_wibe <- function(x, ci = FALSE, style = "line",
                                      stats = NULL) {
  stats <- .normalize_wibe_stats(stats,
                                  default = c("tau", "tau_b", "tau_w"))

  is_shapley <- identical(x$order, "shapley")

  # For Shapley, cross-sectional levels are ordering-invariant
  results <- if (is_shapley) x$all_orderings[[1]]$results else x$results

  # Build full time series (ref + non-ref times)
  non_ref_times <- as.numeric(names(results))
  first_r <- results[[1]]
  ref_time <- x$ref
  all_times <- sort(c(ref_time, non_ref_times))
  n_times <- length(all_times)

  # Assemble level values at each time
  tau_B_all     <- numeric(n_times)
  tau_W_all     <- numeric(n_times)
  het_B_all     <- numeric(n_times)
  cov_B_all     <- numeric(n_times)
  rescale_B_all <- numeric(n_times)
  het_W_all     <- numeric(n_times)
  cov_W_all     <- numeric(n_times)
  rescale_W_all <- numeric(n_times)

  for (i in seq_along(all_times)) {
    t <- all_times[i]
    if (t == ref_time) {
      tau_B_all[i] <- first_r$tau_B_t0
      tau_W_all[i] <- first_r$tau_W_t0
      comps <- first_r$components_t0
    } else {
      r <- results[[as.character(t)]]
      tau_B_all[i] <- r$tau_B_t
      tau_W_all[i] <- r$tau_W_t
      comps <- r$components_t
    }
    het_B_all[i]     <- comps$het_B
    cov_B_all[i]     <- comps$cov_B
    rescale_B_all[i] <- comps$rescale_B %||% 0
    het_W_all[i]     <- comps$het_W
    cov_W_all[i]     <- comps$cov_W
    rescale_W_all[i] <- comps$rescale_W %||% 0
  }

  tau_total_all <- tau_B_all + tau_W_all
  add_total <- "tau" %in% stats
  has_subs <- any(c("het_b", "cov_b", "rescale_b",
                    "het_w", "cov_w", "rescale_w") %in% stats)

  # Map stats values to color group (Between/Within) and style group
  # (Effect / Heterogeneity / Sorting / Rescaling)
  comp_map <- list(
    tau_b     = list(color_grp = "Between", style_grp = "Effect",
                     values = tau_B_all),
    tau_w     = list(color_grp = "Within",  style_grp = "Effect",
                     values = tau_W_all),
    het_b     = list(color_grp = "Between", style_grp = "Heterogeneity",
                     values = het_B_all),
    cov_b     = list(color_grp = "Between", style_grp = "Sorting",
                     values = cov_B_all),
    rescale_b = list(color_grp = "Between", style_grp = "Rescaling",
                     values = rescale_B_all),
    het_w     = list(color_grp = "Within",  style_grp = "Heterogeneity",
                     values = het_W_all),
    cov_w     = list(color_grp = "Within",  style_grp = "Sorting",
                     values = cov_W_all),
    rescale_w = list(color_grp = "Within",  style_grp = "Rescaling",
                     values = rescale_W_all)
  )

  # Build plot data
  rows <- list()
  for (s in setdiff(stats, "tau")) {
    m <- comp_map[[s]]
    rows[[length(rows) + 1]] <- data.frame(
      time = all_times, color_grp = m$color_grp, style_grp = m$style_grp,
      value = m$values, stringsAsFactors = FALSE)
  }

  if (add_total) {
    rows[[length(rows) + 1]] <- data.frame(
      time = all_times, color_grp = "Total", style_grp = "Effect",
      value = tau_total_all, stringsAsFactors = FALSE)
  }

  plot_df <- do.call(rbind, rows)

  color_levels <- intersect(c("Total", "Between", "Within"),
                             unique(plot_df$color_grp))
  style_levels <- intersect(c("Effect", "Heterogeneity", "Sorting",
                              "Rescaling"),
                             unique(plot_df$style_grp))
  plot_df$color_grp <- factor(plot_df$color_grp, levels = color_levels)
  plot_df$style_grp <- factor(plot_df$style_grp, levels = style_levels)

  color_vals <- c("Total" = "black", "Between" = "#e85d00", "Within" = "#008b94")
  lty_vals <- c("Effect" = "solid", "Heterogeneity" = "dashed",
                "Sorting" = "dotted", "Rescaling" = "twodash")

  # CI: delta-method SEs (preferred), fall back to bootstrap percentiles
  has_vcov <- !is.null(x$params) && !is.null(x$params$vcov)
  has_boot <- !is.null(x$boot) && !is.null(x$boot$replicates)
  if (ci && !has_vcov && !has_boot) {
    warning("CIs requested but no vcov or bootstrap replicates available. ",
            "Use se = 'delta' or se = boot_config() in ineqx() to enable CIs.",
            call. = FALSE)
  }

  # Map (color_grp, style_grp) -> SE/replicate key
  #   Effect/Total = tau_total, Effect/Between = tau_B, Effect/Within = tau_W,
  #   Heterogeneity/{Between,Within} = het_{B,W},
  #   Sorting/{Between,Within}      = cov_{B,W},
  #   Rescaling/{Between,Within}    = rescale_{B,W}
  se_key_for <- function(color_grp, style_grp) {
    if (style_grp == "Effect") {
      switch(color_grp, "Total" = "tau_total",
             "Between" = "tau_B", "Within" = "tau_W", NA_character_)
    } else if (style_grp == "Heterogeneity") {
      switch(color_grp, "Between" = "het_B", "Within" = "het_W", NA_character_)
    } else if (style_grp == "Sorting") {
      switch(color_grp, "Between" = "cov_B", "Within" = "cov_W", NA_character_)
    } else if (style_grp == "Rescaling") {
      switch(color_grp,
             "Between" = "rescale_B",
             "Within"  = "rescale_W",
             NA_character_)
    } else NA_character_
  }

  if (ci && has_vcov) {
    # Pre-compute cross-sectional delta SEs per time
    se_per_time <- vector("list", n_times)
    for (i in seq_along(all_times)) {
      se_per_time[[i]] <- tryCatch(
        delta_method_se(x$params, type = "cross", ref = all_times[i]),
        error = function(e) NULL)
    }

    se_vec <- rep(NA_real_, nrow(plot_df))
    for (i in seq_len(nrow(plot_df))) {
      key <- se_key_for(as.character(plot_df$color_grp[i]),
                        as.character(plot_df$style_grp[i]))
      if (is.na(key)) next
      time_idx <- which(all_times == plot_df$time[i])
      se_t <- se_per_time[[time_idx]]
      if (is.null(se_t)) next
      se_nm <- paste0("se_", key)
      # tau_* keys live at top level; het_*/cov_* live in se_sub
      if (se_nm %in% names(se_t)) {
        se_vec[i] <- se_t[[se_nm]]
      } else if (!is.null(se_t$se_sub) && se_nm %in% names(se_t$se_sub)) {
        se_vec[i] <- se_t$se_sub[[se_nm]]
      }
    }
    plot_df$se <- se_vec
    plot_df$ymin <- plot_df$value - 1.96 * plot_df$se
    plot_df$ymax <- plot_df$value + 1.96 * plot_df$se
  } else if (ci && has_boot) {
    reps <- x$boot$replicates
    ymin_vec <- rep(NA_real_, nrow(plot_df))
    ymax_vec <- rep(NA_real_, nrow(plot_df))
    se_vec   <- rep(NA_real_, nrow(plot_df))
    for (i in seq_len(nrow(plot_df))) {
      key <- se_key_for(as.character(plot_df$color_grp[i]),
                        as.character(plot_df$style_grp[i]))
      if (is.na(key)) next
      t  <- plot_df$time[i]
      col <- paste0(key, "_", t)
      if (!(col %in% colnames(reps))) next
      rep_vec <- reps[, col]
      ymin_vec[i] <- stats::quantile(rep_vec, 0.025)
      ymax_vec[i] <- stats::quantile(rep_vec, 0.975)
      se_vec[i]   <- stats::sd(rep_vec)
    }
    plot_df$se   <- se_vec
    plot_df$ymin <- ymin_vec
    plot_df$ymax <- ymax_vec
  }

  # Build plot
  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              color = .data$color_grp,
                                              linetype = .data$style_grp,
                                              group = interaction(.data$color_grp,
                                                                  .data$style_grp)))

  has_ci <- "ymin" %in% names(plot_df)
  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    if (has_ci) {
      ci_df <- plot_df[!is.na(plot_df$se), ]
      p <- p + ggplot2::geom_errorbar(
        data = ci_df,
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos)
    }
    p <- p + ggplot2::geom_point(size = 2, position = pos$pos)
  } else {
    if (has_ci) {
      ci_df <- plot_df[!is.na(plot_df$se), ]
      fill_colors <- color_vals
      if ("Total" %in% names(fill_colors)) fill_colors["Total"] <- "grey50"
      p <- p + ggplot2::geom_ribbon(
        data = ci_df,
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$color_grp),
        alpha = 0.15, color = NA, linetype = "solid") +
        ggplot2::scale_fill_manual(values = fill_colors, guide = "none")
    }
    p <- p + ggplot2::geom_line(linewidth = 1)
  }

  # Suppress linetype legend if only one style level
  guide_lty <- if (length(style_levels) > 1) {
    ggplot2::guide_legend()
  } else {
    "none"
  }

  p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Time",
                  y = paste0("Treatment effect on ", x$ystat),
                  title = "Treatment Effect on Inequality Over Time",
                  color = "", linetype = "", fill = "") +
    ggplot2::scale_color_manual(values = color_vals) +
    ggplot2::scale_linetype_manual(values = lty_vals) +
    ggplot2::guides(linetype = guide_lty) +
    theme_ineqx
}


# ============================================================================ #
# type = "params" helpers
# ============================================================================ #

#' @keywords internal
.plot_desc_params <- function(x, ci = FALSE, style = "line") {
  wibe <- x$wibe
  if (!"time" %in% names(wibe)) {
    stop("Time variable required for 'params' plot")
  }

  ci_info <- .resolve_ci(ci, x)

  plot_df <- data.frame(
    time  = rep(wibe$time, 2),
    group = rep(wibe$group, 2),
    param = rep(c("Mean (\u03bc)", "SD (\u03c3)"), each = nrow(wibe)),
    value = c(wibe$mu, wibe$sigma),
    stringsAsFactors = FALSE
  )

  if (ci_info$method != "none") {
    se_mu <- wibe$sigma / sqrt(wibe$n)
    se_sigma <- wibe$sigma / sqrt(2 * (wibe$n - 1))
    plot_df$se <- c(se_mu, se_sigma)
    plot_df$ymin <- plot_df$value - 1.96 * plot_df$se
    plot_df$ymax <- plot_df$value + 1.96 * plot_df$se
  }

  plot_df$param <- factor(plot_df$param, levels = c("Mean (\u03bc)", "SD (\u03c3)"))
  plot_df$group <- factor(plot_df$group)

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              color = .data$group,
                                              group = .data$group))

  has_ci <- "ymin" %in% names(plot_df)
  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos)
    }
    p <- p + ggplot2::geom_point(position = pos$pos)
  } else {
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$group),
        alpha = 0.15, color = NA)
    }
    p <- p + ggplot2::geom_line() + ggplot2::geom_point()
  }

  p <- p + ggplot2::facet_wrap(~ .data$param, scales = "free_y") +
    ggplot2::labs(x = "Time", y = "", color = "Group") +
    theme_ineqx

  if (has_ci && style == "line") {
    p <- p + ggplot2::labs(fill = "Group")
  }
  p
}

#' @keywords internal
.plot_causal_cross_bar <- function(x, ci = FALSE,
                                    stats = "tau") {

  stats <- .normalize_wibe_stats(stats, default = "tau")

  # Legacy renderings preserved verbatim:
  #   "tau" alone         -> three τ bars (Total/Between/Within)
  #   c("tau","het_cov")  -> three τ bars + three stacked het+cov bars
  if (setequal(stats, "tau")) {
    stats <- c("tau", "tau_b", "tau_w")
  } else if (setequal(stats,
                       c("tau", "het_b", "cov_b", "het_w", "cov_w"))) {
    stats <- c("tau", "tau_b", "tau_w",
               "het_b", "cov_b", "het_w", "cov_w")
  }

  # When all four het/cov sub-stats are present, also render a Total-column
  # stacked sub bar with summed values (legacy "het_cov" three-position view).
  is_full_het_cov <- all(c("het_b", "cov_b", "het_w", "cov_w") %in% stats)

  comps     <- x$components
  het_B     <- comps$het_B
  cov_B     <- comps$cov_B
  rescale_B <- comps$rescale_B %||% 0
  het_W     <- comps$het_W
  cov_W     <- comps$cov_W
  rescale_W <- comps$rescale_W %||% 0

  fill_levels <- c("Between", "Within", "Total",
                   "Heterogeneity", "Sorting", "Rescaling")
  fill_vals <- c("Between" = "#e85d00", "Within" = "#008b94",
                 "Total" = "#1a4e66",
                 "Heterogeneity" = "#bdbdbd", "Sorting" = "#252525",
                 "Rescaling" = "#888888")

  # Each canonical stat -> bar specification
  stat_def <- list(
    tau       = list(xpos = 1, fill = "Total",
                     value = x$tau_total, side = "tau"),
    tau_b     = list(xpos = 2, fill = "Between",
                     value = x$tau_B,     side = "tau"),
    tau_w     = list(xpos = 3, fill = "Within",
                     value = x$tau_W,     side = "tau"),
    het_b     = list(xpos = 2, fill = "Heterogeneity",
                     value = het_B,       side = "sub"),
    cov_b     = list(xpos = 2, fill = "Sorting",
                     value = cov_B,       side = "sub"),
    rescale_b = list(xpos = 2, fill = "Rescaling",
                     value = rescale_B,   side = "sub"),
    het_w     = list(xpos = 3, fill = "Heterogeneity",
                     value = het_W,       side = "sub"),
    cov_w     = list(xpos = 3, fill = "Sorting",
                     value = cov_W,       side = "sub"),
    rescale_w = list(xpos = 3, fill = "Rescaling",
                     value = rescale_W,   side = "sub")
  )

  side_positions <- function(side) {
    unique(unlist(lapply(stats, function(s) {
      if (stat_def[[s]]$side == side) stat_def[[s]]$xpos else integer(0)
    })))
  }
  tau_positions <- side_positions("tau")
  sub_positions <- side_positions("sub")
  if (is_full_het_cov) sub_positions <- unique(c(sub_positions, 1))
  dodge_positions <- intersect(tau_positions, sub_positions)
  has_dodge <- length(dodge_positions) > 0

  if (has_dodge) {
    off <- 0.32
    bw_tau <- 0.35
    bw_sub <- 0.22
  } else {
    off <- 0
    bw_tau <- 0.6
    bw_sub <- 0.6
  }

  # τ-side bars
  tau_stats <- intersect(c("tau", "tau_b", "tau_w"), stats)
  tau_df <- NULL
  if (length(tau_stats) > 0) {
    tau_df <- data.frame(
      xpos = vapply(tau_stats, function(s) {
        xb <- stat_def[[s]]$xpos
        if (xb %in% dodge_positions) xb - off / 2 else xb
      }, numeric(1)),
      fill_var = factor(
        vapply(tau_stats, function(s) stat_def[[s]]$fill, character(1)),
        levels = fill_levels),
      value = vapply(tau_stats,
                     function(s) stat_def[[s]]$value, numeric(1)),
      stringsAsFactors = FALSE
    )

    if (ci && !is.null(x[["se"]])) {
      se_keys <- c(tau = "se_tau_total",
                   tau_b = "se_tau_B",
                   tau_w = "se_tau_W")
      ses <- vapply(tau_stats, function(s) {
        v <- x$se[[se_keys[[s]]]]
        if (is.null(v)) NA_real_ else v
      }, numeric(1))
      tau_df$ymin <- tau_df$value - 1.96 * ses
      tau_df$ymax <- tau_df$value + 1.96 * ses
    }
  }

  # Sub-component bars (heterogeneity / sorting / rescaling), stacked
  sub_stats <- intersect(c("het_b", "cov_b", "rescale_b",
                           "het_w", "cov_w", "rescale_w"), stats)
  sub_df <- NULL
  if (length(sub_stats) > 0 || is_full_het_cov) {
    rows <- list()
    for (s in sub_stats) {
      sd <- stat_def[[s]]
      xb <- sd$xpos
      xpos <- if (xb %in% dodge_positions) xb + off / 2 else xb
      rows[[length(rows) + 1L]] <- data.frame(
        xpos = xpos, fill_var = sd$fill, value = sd$value,
        stringsAsFactors = FALSE)
    }
    if (is_full_het_cov) {
      xb <- 1
      xpos <- if (xb %in% dodge_positions) xb + off / 2 else xb
      rows[[length(rows) + 1L]] <- data.frame(
        xpos = xpos, fill_var = "Heterogeneity",
        value = het_B + het_W, stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(
        xpos = xpos, fill_var = "Sorting",
        value = cov_B + cov_W, stringsAsFactors = FALSE)
    }
    sub_df <- do.call(rbind, rows)
    sub_df$fill_var <- factor(sub_df$fill_var, levels = fill_levels)
  }

  p <- ggplot2::ggplot()

  if (!is.null(tau_df)) {
    p <- p + ggplot2::geom_col(
      data = tau_df,
      ggplot2::aes(x = .data$xpos, y = .data$value,
                   fill = .data$fill_var),
      width = bw_tau)

    if (ci && "ymin" %in% names(tau_df)) {
      ci_df <- tau_df[!is.na(tau_df$ymin), , drop = FALSE]
      if (nrow(ci_df) > 0) {
        p <- p + ggplot2::geom_errorbar(
          data = ci_df,
          ggplot2::aes(x = .data$xpos,
                       ymin = .data$ymin, ymax = .data$ymax),
          width = 0.15)
      }
    }
  }

  if (!is.null(sub_df)) {
    p <- p + ggplot2::geom_col(
      data = sub_df,
      ggplot2::aes(x = .data$xpos, y = .data$value,
                   fill = .data$fill_var),
      width = bw_sub, position = "stack")
  }

  p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")

  fills_present <- character(0)
  if (!is.null(tau_df)) {
    fills_present <- c(fills_present, as.character(tau_df$fill_var))
  }
  if (!is.null(sub_df)) {
    fills_present <- c(fills_present, as.character(sub_df$fill_var))
  }
  brks <- intersect(c("Total", "Between", "Within",
                       "Heterogeneity", "Sorting", "Rescaling"),
                    unique(fills_present))

  p + ggplot2::scale_x_continuous(
        breaks = c(1, 2, 3),
        labels = c("Total", "Between", "Within")) +
    ggplot2::scale_fill_manual(
      values = fill_vals,
      breaks = brks,
      labels = brks) +
    ggplot2::labs(x = "", y = paste0("Treatment effect on ", x$ystat),
                  title = "Treatment Effect on Inequality",
                  fill = "") +
    theme_ineqx
}

#' @keywords internal
.plot_causal_cross_treat <- function(x, ci = FALSE) {
  use_ci <- .resolve_ci_params(ci)
  bg <- x$by_group

  plot_df <- data.frame(
    group = rep(bg$group, 2),
    param = rep(c("beta", "lambda"), each = nrow(bg)),
    value = c(bg$beta, bg$lambda),
    stringsAsFactors = FALSE
  )

  # Extract SEs: prefer vcov (delta method), fall back to bootstrap replicates
  has_boot <- !is.null(x$boot) && !is.null(x$boot$replicates)
  if (use_ci && is.null(x$params$vcov) && !has_boot) {
    warning("CIs requested but no vcov available in params. ",
            "Use se = 'delta' or se = boot_config() in ineqx() to enable parameter CIs.",
            call. = FALSE)
  }
  if (use_ci && !is.null(x$params$vcov)) {
    V <- x$params$vcov
    J <- x$params$n_groups
    se_beta <- sqrt(diag(V)[1:J])
    se_lambda <- sqrt(diag(V)[(2 * J + 1):(3 * J)])
    plot_df$se <- c(se_beta, se_lambda)
    plot_df$ymin <- plot_df$value - 1.96 * plot_df$se
    plot_df$ymax <- plot_df$value + 1.96 * plot_df$se
  } else if (use_ci && is.null(x$params$vcov) && has_boot) {
    reps <- x$boot$replicates
    groups <- bg$group
    ci_lo_beta   <- numeric(length(groups))
    ci_hi_beta   <- numeric(length(groups))
    ci_lo_lambda <- numeric(length(groups))
    ci_hi_lambda <- numeric(length(groups))
    for (i in seq_along(groups)) {
      g <- groups[i]
      col_b <- paste0("beta_",   g)
      col_l <- paste0("lambda_", g)
      if (col_b %in% colnames(reps)) {
        ci_lo_beta[i] <- stats::quantile(reps[, col_b], 0.025)
        ci_hi_beta[i] <- stats::quantile(reps[, col_b], 0.975)
      }
      if (col_l %in% colnames(reps)) {
        ci_lo_lambda[i] <- stats::quantile(reps[, col_l], 0.025)
        ci_hi_lambda[i] <- stats::quantile(reps[, col_l], 0.975)
      }
    }
    plot_df$ymin <- c(ci_lo_beta, ci_lo_lambda)
    plot_df$ymax <- c(ci_hi_beta, ci_hi_lambda)
  }

  plot_df$param <- factor(plot_df$param,
    levels = c("beta", "lambda"),
    labels = c("Effect on mean (\u03b2)", "Effect on SD (\u03bb)"))
  plot_df$group <- factor(plot_df$group)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$group, y = .data$value,
                                              fill = .data$group))

  if (use_ci && "ymin" %in% names(plot_df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
      width = 0.2)
  }

  p + ggplot2::facet_wrap(~ .data$param, scales = "free_y") +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Group", y = "", fill = "Group") +
    ggplot2::guides(fill = "none") +
    theme_ineqx
}


# ============================================================================ #
# type = "outcome": predicted means & SDs under control vs treatment
# ============================================================================ #

#' @keywords internal
.plot_causal_cross_outcome <- function(x, ci = FALSE) {
  ci <- .resolve_ci_params(ci)

  # Branch: DiD models get the netted-out classic DiD line chart anchored
  # at the pre-period treated level. Simple-difference models keep the
  # original bar-chart form.
  if (isTRUE(x$params$is_did)) {
    return(.plot_causal_cross_outcome_did(x, ci = ci))
  }

  bg <- x$by_group
  mu1    <- bg$mu0 + bg$beta
  sigma1 <- bg$sigma0 * exp(bg$lambda)
  J <- nrow(bg)

  plot_df <- data.frame(
    group     = rep(bg$group, 4),
    param     = rep(c("mu", "mu", "sigma", "sigma"), each = J),
    Condition = rep(c("Control", "Treated", "Control", "Treated"), each = J),
    value     = c(bg$mu0, mu1, bg$sigma0, sigma1),
    stringsAsFactors = FALSE
  )
  plot_df$group     <- factor(plot_df$group)
  plot_df$param     <- factor(plot_df$param,
    levels = c("mu", "sigma"),
    labels = c("Mean (\u03bc)", "SD (\u03c3)"))
  plot_df$Condition <- factor(plot_df$Condition, levels = c("Control", "Treated"))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$group, y = .data$value,
                                              fill = .data$Condition))

  has_boot <- !is.null(x$boot) && !is.null(x$boot$replicates)
  if (ci && is.null(x$params$vcov) && !has_boot) {
    warning("CIs requested but no vcov available in params. ",
            "Use se = 'delta' or se = boot_config() in ineqx() to enable parameter CIs.",
            call. = FALSE)
  }

  if (ci && !is.null(x$params$vcov)) {
    V <- x$params$vcov
    # vcov layout: (beta_1..J, mu0_1..J, lambda_1..J, log_sigma0_1..J)
    se_beta       <- sqrt(diag(V)[1:J])
    se_mu0        <- sqrt(diag(V)[(J + 1):(2 * J)])
    se_lambda     <- sqrt(diag(V)[(2 * J + 1):(3 * J)])
    se_log_sigma0 <- sqrt(diag(V)[(3 * J + 1):(4 * J)])
    se_sigma0     <- bg$sigma0 * se_log_sigma0  # delta method: sigma0 = exp(log_sigma0)

    # Treated mu = mu0 + beta: se = sqrt(se_mu0^2 + se_beta^2) (approx, ignoring cov)
    se_mu1    <- sqrt(se_mu0^2 + se_beta^2)
    # Treated sigma = sigma0 * exp(lambda): delta method
    se_sigma1 <- sqrt((exp(bg$lambda) * se_sigma0)^2 +
                       (bg$sigma0 * exp(bg$lambda) * se_lambda)^2)

    plot_df$se <- c(se_mu0, se_mu1, se_sigma0, se_sigma1)
    plot_df$ymin <- plot_df$value - 1.96 * plot_df$se
    plot_df$ymax <- plot_df$value + 1.96 * plot_df$se
  } else if (ci && is.null(x$params$vcov) && has_boot) {
    # Percentile CIs from bootstrap replicates
    reps <- x$boot$replicates
    ymin_vec <- rep(NA_real_, 4 * J)
    ymax_vec <- rep(NA_real_, 4 * J)
    for (i in seq_len(J)) {
      g <- bg$group[i]
      col_mu0    <- paste0("mu0_",    g)
      col_sigma0 <- paste0("sigma0_", g)
      col_beta   <- paste0("beta_",   g)
      col_lambda <- paste0("lambda_", g)
      if (!all(c(col_mu0, col_sigma0, col_beta, col_lambda) %in% colnames(reps))) next
      mu0_rep    <- reps[, col_mu0]
      sigma0_rep <- reps[, col_sigma0]
      mu1_rep    <- mu0_rep + reps[, col_beta]
      sigma1_rep <- sigma0_rep * exp(reps[, col_lambda])
      # plot_df row order: mu/Control (1..J), mu/Treated (J+1..2J),
      #                    sigma/Control (2J+1..3J), sigma/Treated (3J+1..4J)
      ymin_vec[i]           <- stats::quantile(mu0_rep,    0.025)
      ymax_vec[i]           <- stats::quantile(mu0_rep,    0.975)
      ymin_vec[J + i]       <- stats::quantile(mu1_rep,    0.025)
      ymax_vec[J + i]       <- stats::quantile(mu1_rep,    0.975)
      ymin_vec[2 * J + i]   <- stats::quantile(sigma0_rep, 0.025)
      ymax_vec[2 * J + i]   <- stats::quantile(sigma0_rep, 0.975)
      ymin_vec[3 * J + i]   <- stats::quantile(sigma1_rep, 0.025)
      ymax_vec[3 * J + i]   <- stats::quantile(sigma1_rep, 0.975)
    }
    plot_df$ymin <- ymin_vec
    plot_df$ymax <- ymax_vec
  }

  if (ci && "ymin" %in% names(plot_df) && any(!is.na(plot_df$ymin))) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$group, y = .data$value,
                                                fill = .data$Condition))
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                   group = .data$Condition),
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2)
  }

  p + ggplot2::facet_wrap(~ .data$param, scales = "free_y") +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.6),
                      width = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_fill_manual(values = c("Control" = "grey40",
                                           "Treated" = "#008b94")) +
    ggplot2::labs(x = "Group", y = "",
                  title = "Predicted Levels by Group",
                  fill = "") +
    theme_ineqx
}

#' Cross-sectional DiD outcome.params: netted-out classic DiD line chart
#'
#' Two lines per group, both anchored at the pre-period treated level
#' (\code{mu1_pre} / \code{sigma1_pre}). The "Treated" line ends at the
#' observed treated post-period level (\code{mu1}); the "Counterfactual"
#' line ends at the DiD-implied counterfactual untreated post-period level
#' (\code{mu0}). Lines coincide at "pre" by construction; the post-period
#' gap equals the ATT (beta) on the mu panel and lambda on the log-SD
#' panel.
#' @keywords internal
.plot_causal_cross_outcome_did <- function(x, ci = FALSE) {
  pdata <- x$params$data
  groups <- unique(pdata$group)
  J <- length(groups)

  # Build long-format data: 4 rows per (group, param) — pre and post for each
  # of Treated and Counterfactual.
  rows <- list()
  for (g in groups) {
    row <- pdata[pdata$group == g, ]
    # mu panel
    rows[[length(rows) + 1]] <- data.frame(
      group = g, param = "mu", Condition = "Treated",
      period = "pre",  value = row$mu1_pre, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      group = g, param = "mu", Condition = "Treated",
      period = "post", value = row$mu1, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      group = g, param = "mu", Condition = "Counterfactual",
      period = "pre",  value = row$mu1_pre, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      group = g, param = "mu", Condition = "Counterfactual",
      period = "post", value = row$mu0, stringsAsFactors = FALSE)
    # sigma panel
    rows[[length(rows) + 1]] <- data.frame(
      group = g, param = "sigma", Condition = "Treated",
      period = "pre",  value = row$sigma1_pre, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      group = g, param = "sigma", Condition = "Treated",
      period = "post", value = row$sigma1, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      group = g, param = "sigma", Condition = "Counterfactual",
      period = "pre",  value = row$sigma1_pre, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      group = g, param = "sigma", Condition = "Counterfactual",
      period = "post", value = row$sigma0, stringsAsFactors = FALSE)
  }
  plot_df <- do.call(rbind, rows)
  plot_df$group     <- factor(plot_df$group, levels = groups)
  plot_df$param     <- factor(plot_df$param,
                              levels = c("mu", "sigma"),
                              labels = c("Mean (μ)", "SD (σ)"))
  plot_df$Condition <- factor(plot_df$Condition,
                              levels = c("Treated", "Counterfactual"))
  plot_df$period    <- factor(plot_df$period, levels = c("pre", "post"))

  ggplot2::ggplot(plot_df, ggplot2::aes(
      x = .data$period, y = .data$value,
      color = .data$group,
      linetype = .data$Condition,
      group = interaction(.data$group, .data$Condition))) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(~ .data$param, scales = "free_y") +
    ggplot2::scale_linetype_manual(values = c("Treated"        = "solid",
                                                "Counterfactual" = "dashed")) +
    ggplot2::labs(x = NULL, y = "",
                  color = "Group", linetype = "",
                  title = "Predicted Levels: Observed vs DiD Counterfactual",
                  subtitle = "Lines anchored at pre-period treated level; post-period gap = ATT") +
    theme_ineqx
}

#' Generate dist plot title
#' @keywords internal
.dist_title <- function(show) {
  if (show == "outcome") {
    "Predicted outcome distribution"
  } else {
    "Predicted treatment effect distribution"
  }
}

#' Build treatment effect dist data (one density per group)
#' @keywords internal
.build_te_dist_data <- function(groups, mu0, sigma0, beta, lambda,
                                 trim = 0.995) {
  sigma1 <- sigma0 * exp(lambda)
  J <- length(groups)
  # Compute common x-range across all groups
  te_sds <- sqrt(sigma0^2 + sigma1^2)
  x_lo <- min(stats::qnorm(1 - trim, mean = beta, sd = te_sds))
  x_hi <- max(stats::qnorm(trim, mean = beta, sd = te_sds))
  x_seq <- seq(x_lo, x_hi, length.out = 200)
  rows <- list()
  for (i in seq_len(J)) {
    rows[[length(rows) + 1]] <- data.frame(
      group = groups[i], x = x_seq,
      density = stats::dnorm(x_seq, mean = beta[i], sd = te_sds[i]),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' Build outcome dist data (Control + Treated per group)
#' @keywords internal
.build_outcome_dist_data <- function(groups, mu0, sigma0, beta, lambda,
                                      trim = 0.995) {
  mu1    <- mu0 + beta
  sigma1 <- sigma0 * exp(lambda)
  # Compute common x-range across all groups and conditions
  all_mu <- c(mu0, mu1)
  all_sd <- c(sigma0, sigma1)
  x_lo <- min(stats::qnorm(1 - trim, mean = all_mu, sd = all_sd))
  x_hi <- max(stats::qnorm(trim, mean = all_mu, sd = all_sd))
  x_seq <- seq(x_lo, x_hi, length.out = 200)
  J <- length(groups)
  rows <- list()
  for (i in seq_len(J)) {
    rows[[length(rows) + 1]] <- data.frame(
      group = groups[i], Condition = "Control", x = x_seq,
      density = stats::dnorm(x_seq, mean = mu0[i], sd = sigma0[i]),
      stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      group = groups[i], Condition = "Treated", x = x_seq,
      density = stats::dnorm(x_seq, mean = mu1[i], sd = sigma1[i]),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' @keywords internal
.plot_causal_cross_dist <- function(x, show = "treat", trim = 0.995) {
  bg <- x$by_group
  title <- .dist_title(show)

  if (show == "treat") {
    plot_df <- .build_te_dist_data(bg$group, bg$mu0, bg$sigma0,
                                    bg$beta, bg$lambda, trim = trim)
    plot_df$group <- factor(plot_df$group, levels = bg$group)

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$density,
                                           color = .data$group,
                                           fill = .data$group)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_area(alpha = 0.15, position = "identity") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "grey50") +
      ggplot2::labs(x = x$y, y = "Density", title = title,
                    color = "Group", fill = "Group") +
      theme_ineqx
  } else {
    plot_df <- .build_outcome_dist_data(bg$group, bg$mu0, bg$sigma0,
                                         bg$beta, bg$lambda, trim = trim)
    plot_df$group <- factor(plot_df$group, levels = bg$group)
    plot_df$Condition <- factor(plot_df$Condition,
                                 levels = c("Control", "Treated"))

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$density,
                                           color = .data$Condition,
                                           fill = .data$Condition)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_area(alpha = 0.15, position = "identity") +
      ggplot2::facet_wrap(~ .data$group, ncol = 1, scales = "free_y") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "grey50") +
      ggplot2::scale_color_manual(values = c("Control" = "grey40",
                                               "Treated" = "#008b94")) +
      ggplot2::scale_fill_manual(values = c("Control" = "grey40",
                                              "Treated" = "#008b94")) +
      ggplot2::labs(x = x$y, y = "Density", title = title,
                    color = "", fill = "") +
      theme_ineqx
  }
}

#' @keywords internal
.plot_causal_longit_outcome <- function(x, ci = FALSE, style = "line") {
  ci <- .resolve_ci_params(ci)
  is_shapley <- identical(x$order, "shapley")
  results <- if (is_shapley) x$all_orderings[[1]]$results else x$results

  non_ref_times <- as.numeric(names(results))
  ref_time <- x$ref
  all_times <- sort(c(ref_time, non_ref_times))

  pdata <- x$params$data
  groups <- unique(pdata$group)
  J <- x$params$n_groups

  # Check CI availability: prefer vcov (delta), fall back to bootstrap replicates
  has_boot <- !is.null(x$boot) && !is.null(x$boot$replicates)
  if (ci && is.null(x$params$vcov) && !has_boot) {
    warning("CIs requested but no vcov available in params. ",
            "Use se = 'delta' or se = boot_config() in ineqx() to enable parameter CIs.",
            call. = FALSE)
  }
  has_vcov <- ci && !is.null(x$params$vcov)
  use_boot <- ci && !has_vcov && has_boot
  reps <- if (use_boot) x$boot$replicates else NULL

  rows <- list()
  for (t in all_times) {
    # Extract vcov for this time period
    t_key <- as.character(t)
    V <- NULL
    if (has_vcov) {
      V <- if (is.list(x$params$vcov)) x$params$vcov[[t_key]] else x$params$vcov
    }

    for (g_i in seq_along(groups)) {
      g <- groups[g_i]
      idx <- pdata$time == t & pdata$group == g
      if (!any(idx)) next
      row <- pdata[idx, ]
      # Use stored mu1/sigma1 if available (DiD models post-fix), else
      # reconstruct from mu0+beta. Equivalent for simple-difference; for
      # DiD this guarantees the plotted gap equals the ATT.
      mu1    <- if (!is.null(row$mu1)    && !is.na(row$mu1))    row$mu1
                else row$mu0 + row$beta
      sigma1 <- if (!is.null(row$sigma1) && !is.na(row$sigma1)) row$sigma1
                else row$sigma0 * exp(row$lambda)

      ymin_mu0 <- ymax_mu0 <- ymin_mu1 <- ymax_mu1 <- NA_real_
      ymin_sigma0 <- ymax_sigma0 <- ymin_sigma1 <- ymax_sigma1 <- NA_real_

      if (!is.null(V)) {
        # SEs via delta method
        # vcov layout: (beta_1..J, mu0_1..J, lambda_1..J, log_sigma0_1..J)
        se_beta       <- sqrt(diag(V)[g_i])
        se_mu0        <- sqrt(diag(V)[J + g_i])
        se_lambda     <- sqrt(diag(V)[2 * J + g_i])
        se_log_sigma0 <- sqrt(diag(V)[3 * J + g_i])
        se_sigma0     <- row$sigma0 * se_log_sigma0  # delta: sigma0 = exp(log_sigma0)
        se_mu1        <- sqrt(se_mu0^2 + se_beta^2)
        se_sigma1     <- sqrt((exp(row$lambda) * se_sigma0)^2 +
                              (row$sigma0 * exp(row$lambda) * se_lambda)^2)
        ymin_mu0    <- row$mu0 - 1.96 * se_mu0;   ymax_mu0    <- row$mu0 + 1.96 * se_mu0
        ymin_mu1    <- mu1     - 1.96 * se_mu1;   ymax_mu1    <- mu1     + 1.96 * se_mu1
        ymin_sigma0 <- row$sigma0 - 1.96 * se_sigma0
        ymax_sigma0 <- row$sigma0 + 1.96 * se_sigma0
        ymin_sigma1 <- sigma1  - 1.96 * se_sigma1
        ymax_sigma1 <- sigma1  + 1.96 * se_sigma1
      } else if (use_boot) {
        # Percentile CIs from bootstrap replicates
        col_mu0    <- paste0("mu0_",    g, "_", t)
        col_sigma0 <- paste0("sigma0_", g, "_", t)
        col_beta   <- paste0("beta_",   g, "_", t)
        col_lambda <- paste0("lambda_", g, "_", t)
        if (all(c(col_mu0, col_sigma0, col_beta, col_lambda) %in% colnames(reps))) {
          mu0_rep    <- reps[, col_mu0]
          sigma0_rep <- reps[, col_sigma0]
          mu1_rep    <- mu0_rep + reps[, col_beta]
          sigma1_rep <- sigma0_rep * exp(reps[, col_lambda])
          ymin_mu0    <- stats::quantile(mu0_rep,    0.025); ymax_mu0    <- stats::quantile(mu0_rep,    0.975)
          ymin_mu1    <- stats::quantile(mu1_rep,    0.025); ymax_mu1    <- stats::quantile(mu1_rep,    0.975)
          ymin_sigma0 <- stats::quantile(sigma0_rep, 0.025); ymax_sigma0 <- stats::quantile(sigma0_rep, 0.975)
          ymin_sigma1 <- stats::quantile(sigma1_rep, 0.025); ymax_sigma1 <- stats::quantile(sigma1_rep, 0.975)
        }
      }

      rows[[length(rows) + 1]] <- data.frame(
        time = t, group = g, param = "mu", Condition = "Control",
        value = row$mu0, ymin = ymin_mu0, ymax = ymax_mu0,
        stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        time = t, group = g, param = "mu", Condition = "Treated",
        value = mu1, ymin = ymin_mu1, ymax = ymax_mu1,
        stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        time = t, group = g, param = "sigma", Condition = "Control",
        value = row$sigma0, ymin = ymin_sigma0, ymax = ymax_sigma0,
        stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        time = t, group = g, param = "sigma", Condition = "Treated",
        value = sigma1, ymin = ymin_sigma1, ymax = ymax_sigma1,
        stringsAsFactors = FALSE)
    }
  }
  plot_df <- do.call(rbind, rows)
  plot_df$group     <- factor(plot_df$group)
  plot_df$param     <- factor(plot_df$param,
    levels = c("mu", "sigma"),
    labels = c("Mean (\u03bc)", "SD (\u03c3)"))
  plot_df$Condition <- factor(plot_df$Condition, levels = c("Control", "Treated"))

  has_ci <- (has_vcov || use_boot) && any(!is.na(plot_df$ymin))

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data$time, y = .data$value,
    color = .data$group,
    linetype = .data$Condition,
    group = interaction(.data$group, .data$Condition)))

  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos)
    }
    p <- p + ggplot2::geom_point(size = 2, position = pos$pos)
  } else {
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$group),
        alpha = 0.15, color = NA)
    }
    p <- p + ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2)
  }

  is_did <- isTRUE(x$params$is_did)
  plot_title <- if (is_did) {
    "Predicted Post-Period Levels Over Time (Gap = DiD ATT)"
  } else {
    "Predicted Levels Over Time"
  }

  p <- p + ggplot2::facet_wrap(~ .data$param, scales = "free_y") +
    ggplot2::scale_linetype_manual(values = c("Control" = "dashed",
                                                "Treated" = "solid")) +
    ggplot2::labs(x = "Time", y = "",
                  color = "Group", linetype = "",
                  title = plot_title) +
    theme_ineqx

  if (has_ci && style == "line") {
    p <- p + ggplot2::labs(fill = "Group")
  }
  p
}

#' Longitudinal pretrends diagnostic plot
#'
#' For DiD models, plots the pre-period predicted Treated and Control levels
#' across calendar time. Under parallel trends, the two lines should evolve
#' with the same slope; their (constant) vertical gap reflects pre-existing
#' selection (\eqn{\beta_D}). Diverging slopes indicate a violation of
#' parallel trends.
#'
#' Both panels (mu and sigma) are shown.
#'
#' @keywords internal
.plot_causal_longit_pretrends <- function(x, ci = FALSE, style = "line") {
  pdata  <- x$params$data
  groups <- unique(pdata$group)
  times  <- sort(unique(pdata$time))

  rows <- list()
  for (t in times) {
    for (g in groups) {
      idx <- pdata$time == t & pdata$group == g
      if (!any(idx)) next
      row <- pdata[idx, ]
      if (any(is.na(c(row$mu0_pre, row$mu1_pre,
                      row$sigma0_pre, row$sigma1_pre)))) next

      rows[[length(rows) + 1]] <- data.frame(
        time = t, group = g, param = "mu", Condition = "Control",
        value = row$mu0_pre, stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        time = t, group = g, param = "mu", Condition = "Treated",
        value = row$mu1_pre, stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        time = t, group = g, param = "sigma", Condition = "Control",
        value = row$sigma0_pre, stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        time = t, group = g, param = "sigma", Condition = "Treated",
        value = row$sigma1_pre, stringsAsFactors = FALSE)
    }
  }

  if (length(rows) == 0) {
    stop("No pre-period anchor values (mu0_pre/mu1_pre/sigma0_pre/sigma1_pre) ",
         "found in params$data. The DiD model may need to be re-fit.")
  }

  plot_df <- do.call(rbind, rows)
  plot_df$group     <- factor(plot_df$group)
  plot_df$param     <- factor(plot_df$param,
                              levels = c("mu", "sigma"),
                              labels = c("Mean (μ)", "SD (σ)"))
  plot_df$Condition <- factor(plot_df$Condition, levels = c("Control", "Treated"))
  plot_df$time      <- .time_to_factor(plot_df$time)

  ggplot2::ggplot(plot_df, ggplot2::aes(
      x = .data$time, y = .data$value,
      color = .data$group,
      linetype = .data$Condition,
      group = interaction(.data$group, .data$Condition))) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ .data$param, scales = "free_y") +
    ggplot2::scale_linetype_manual(values = c("Control" = "dashed",
                                                "Treated" = "solid")) +
    ggplot2::labs(x = "Time", y = "",
                  color = "Group", linetype = "",
                  title = "Pre-Period Predicted Levels (Parallel-Trends Diagnostic)",
                  subtitle = "Lines should evolve with the same slope under parallel trends") +
    theme_ineqx
}

#' @keywords internal
.plot_causal_longit_dist <- function(x, time = NULL, show = "treat",
                                     trim = 0.995) {
  is_shapley <- identical(x$order, "shapley")
  pdata <- x$params$data

  if (!is.null(time)) {
    # --- Single time point: same as cross-sectional ---
    sub <- pdata[pdata$time == time, ]
    if (nrow(sub) == 0) {
      stop("No data for time = ", time, ". Available: ",
           paste(sort(unique(pdata$time)), collapse = ", "))
    }
    title <- .dist_title(show)

    if (show == "treat") {
      plot_df <- .build_te_dist_data(sub$group, sub$mu0, sub$sigma0,
                                      sub$beta, sub$lambda, trim = trim)
      plot_df$group <- factor(plot_df$group, levels = sub$group)

      return(ggplot2::ggplot(plot_df, ggplot2::aes(
          x = .data$x, y = .data$density,
          color = .data$group, fill = .data$group)) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_area(alpha = 0.15, position = "identity") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                             color = "grey50") +
        ggplot2::labs(x = x$y, y = "Density", title = title,
                      color = "Group", fill = "Group") +
        theme_ineqx)
    } else {
      plot_df <- .build_outcome_dist_data(sub$group, sub$mu0, sub$sigma0,
                                           sub$beta, sub$lambda, trim = trim)
      plot_df$group <- factor(plot_df$group, levels = sub$group)
      plot_df$Condition <- factor(plot_df$Condition,
                                   levels = c("Control", "Treated"))

      return(ggplot2::ggplot(plot_df, ggplot2::aes(
          x = .data$x, y = .data$density,
          color = .data$Condition, fill = .data$Condition)) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_area(alpha = 0.15, position = "identity") +
        ggplot2::facet_wrap(~ .data$group, ncol = 1, scales = "free_y") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                             color = "grey50") +
        ggplot2::scale_color_manual(values = c("Control" = "grey40",
                                                 "Treated" = "#008b94")) +
        ggplot2::scale_fill_manual(values = c("Control" = "grey40",
                                                "Treated" = "#008b94")) +
        ggplot2::labs(x = x$y, y = "Density", title = title,
                      color = "", fill = "") +
        theme_ineqx)
    }
  }

  # --- All time points: pi-weighted marginal distributions, stacked ---
  results <- if (is_shapley) x$all_orderings[[1]]$results else x$results
  non_ref_times <- as.numeric(names(results))
  ref_time <- x$ref
  all_times <- sort(c(ref_time, non_ref_times))
  title <- .dist_title(show)

  # Pre-compute common x-range across ALL time points for axis alignment
  global_lo <- Inf
  global_hi <- -Inf
  for (t in all_times) {
    sub <- pdata[pdata$time == t, ]
    if (nrow(sub) == 0) next
    mu1    <- sub$mu0 + sub$beta
    sigma1 <- sub$sigma0 * exp(sub$lambda)
    if (show == "outcome") {
      all_mu <- c(sub$mu0, mu1)
      all_sd <- c(sub$sigma0, sigma1)
      lo <- min(stats::qnorm(1 - trim, mean = all_mu, sd = all_sd))
      hi <- max(stats::qnorm(trim, mean = all_mu, sd = all_sd))
    } else {
      te_sd_j <- sqrt(sub$sigma0^2 + sigma1^2)
      lo <- min(stats::qnorm(1 - trim, mean = sub$beta, sd = te_sd_j))
      hi <- max(stats::qnorm(trim, mean = sub$beta, sd = te_sd_j))
    }
    global_lo <- min(global_lo, lo)
    global_hi <- max(global_hi, hi)
  }
  x_seq <- seq(global_lo, global_hi, length.out = 300)

  rows <- list()
  for (t in all_times) {
    sub <- pdata[pdata$time == t, ]
    if (nrow(sub) == 0) next

    mu1    <- sub$mu0 + sub$beta
    sigma1 <- sub$sigma0 * exp(sub$lambda)
    pi_j   <- sub$pi

    if (show == "outcome") {
      dens_ctrl <- rowSums(vapply(seq_len(nrow(sub)), function(j) {
        pi_j[j] * stats::dnorm(x_seq, mean = sub$mu0[j], sd = sub$sigma0[j])
      }, numeric(length(x_seq))))
      dens_treat <- rowSums(vapply(seq_len(nrow(sub)), function(j) {
        pi_j[j] * stats::dnorm(x_seq, mean = mu1[j], sd = sigma1[j])
      }, numeric(length(x_seq))))

      rows[[length(rows) + 1]] <- data.frame(
        time_label = as.character(t), Condition = "Control",
        x = x_seq, density = dens_ctrl, stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        time_label = as.character(t), Condition = "Treated",
        x = x_seq, density = dens_treat, stringsAsFactors = FALSE)
    } else {
      te_sd_j <- sqrt(sub$sigma0^2 + sigma1^2)

      dens_te <- rowSums(vapply(seq_len(nrow(sub)), function(j) {
        pi_j[j] * stats::dnorm(x_seq, mean = sub$beta[j], sd = te_sd_j[j])
      }, numeric(length(x_seq))))

      rows[[length(rows) + 1]] <- data.frame(
        time_label = as.character(t), Condition = "Effect",
        x = x_seq, density = dens_te, stringsAsFactors = FALSE)
    }
  }

  plot_df <- do.call(rbind, rows)
  time_lvls <- as.character(sort(all_times))
  plot_df$time_label <- factor(plot_df$time_label, levels = time_lvls)

  if (show == "outcome") {
    plot_df$Condition <- factor(plot_df$Condition,
                                 levels = c("Control", "Treated"))
    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$density,
                                           color = .data$Condition,
                                           fill = .data$Condition)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_area(alpha = 0.15, position = "identity") +
      ggplot2::facet_wrap(~ .data$time_label, ncol = 1, scales = "free_y") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "grey50") +
      ggplot2::scale_color_manual(values = c("Control" = "grey40",
                                               "Treated" = "#008b94")) +
      ggplot2::scale_fill_manual(values = c("Control" = "grey40",
                                              "Treated" = "#008b94")) +
      ggplot2::labs(x = x$y, y = "Density", title = title,
                    color = "", fill = "") +
      theme_ineqx
  } else {
    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$density)) +
      ggplot2::geom_line(linewidth = 0.8, color = "#1a4e66") +
      ggplot2::geom_area(alpha = 0.15, fill = "#1a4e66") +
      ggplot2::facet_wrap(~ .data$time_label, ncol = 1, scales = "free_y") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "grey50") +
      ggplot2::labs(x = x$y, y = "Density", title = title) +
      theme_ineqx
  }
}

#' @keywords internal
.plot_causal_longit_treat <- function(x, ci = FALSE, style = "line") {
  use_ci <- .resolve_ci_params(ci)

  pdata <- x$params$data

  plot_df <- data.frame(
    time  = rep(pdata$time, 2),
    group = rep(pdata$group, 2),
    param = rep(c("beta", "lambda"), each = nrow(pdata)),
    value = c(pdata$beta, pdata$lambda),
    stringsAsFactors = FALSE
  )

  # Extract SEs: prefer vcov (delta method), fall back to bootstrap replicates
  has_boot <- !is.null(x$boot) && !is.null(x$boot$replicates)
  if (use_ci && is.null(x$params$vcov) && !has_boot) {
    warning("CIs requested but no vcov available in params. ",
            "Use se = 'delta' or se = boot_config() in ineqx() to enable parameter CIs.",
            call. = FALSE)
  }
  if (use_ci && !is.null(x$params$vcov)) {
    J <- x$params$n_groups
    time_levels <- sort(unique(pdata$time))
    se_beta_all <- numeric(nrow(pdata))
    se_lambda_all <- numeric(nrow(pdata))

    for (i in seq_along(time_levels)) {
      t <- time_levels[i]
      t_key <- as.character(t)
      V <- if (is.list(x$params$vcov)) x$params$vcov[[t_key]] else x$params$vcov
      if (!is.null(V)) {
        idx <- pdata$time == t
        se_beta_all[idx] <- sqrt(diag(V)[1:J])
        se_lambda_all[idx] <- sqrt(diag(V)[(2 * J + 1):(3 * J)])
      }
    }

    plot_df$se <- c(se_beta_all, se_lambda_all)
    plot_df$ymin <- plot_df$value - 1.96 * plot_df$se
    plot_df$ymax <- plot_df$value + 1.96 * plot_df$se
  } else if (use_ci && is.null(x$params$vcov) && has_boot) {
    # Bootstrap-based CIs from replicate columns named beta_<g>_<t> / lambda_<g>_<t>
    reps <- x$boot$replicates
    se_beta_all   <- numeric(nrow(pdata))
    se_lambda_all <- numeric(nrow(pdata))
    ci_lo_beta    <- numeric(nrow(pdata))
    ci_hi_beta    <- numeric(nrow(pdata))
    ci_lo_lambda  <- numeric(nrow(pdata))
    ci_hi_lambda  <- numeric(nrow(pdata))

    for (i in seq_len(nrow(pdata))) {
      g <- pdata$group[i]
      t <- pdata$time[i]
      col_b <- paste0("beta_",   g, "_", t)
      col_l <- paste0("lambda_", g, "_", t)
      if (col_b %in% colnames(reps)) {
        se_beta_all[i]  <- stats::sd(reps[, col_b])
        ci_lo_beta[i]   <- stats::quantile(reps[, col_b], 0.025)
        ci_hi_beta[i]   <- stats::quantile(reps[, col_b], 0.975)
      }
      if (col_l %in% colnames(reps)) {
        se_lambda_all[i] <- stats::sd(reps[, col_l])
        ci_lo_lambda[i]  <- stats::quantile(reps[, col_l], 0.025)
        ci_hi_lambda[i]  <- stats::quantile(reps[, col_l], 0.975)
      }
    }

    plot_df$ymin <- c(ci_lo_beta, ci_lo_lambda)
    plot_df$ymax <- c(ci_hi_beta, ci_hi_lambda)
  }

  plot_df$param <- factor(plot_df$param,
    levels = c("beta", "lambda"),
    labels = c("Effect on mean (\u03b2)", "Effect on SD (\u03bb)"))
  plot_df$group <- factor(plot_df$group)

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              color = .data$group,
                                              group = .data$group))

  has_ci <- use_ci && "ymin" %in% names(plot_df)
  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos)
    }
    p <- p + ggplot2::geom_point(position = pos$pos)
  } else {
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$group),
        alpha = 0.15, color = NA)
    }
    p <- p + ggplot2::geom_line() + ggplot2::geom_point()
  }

  p <- p + ggplot2::facet_wrap(~ .data$param, scales = "free_y") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Time", y = "", color = "Group") +
    theme_ineqx

  # Only set fill label when ribbon CI is active (avoids ggplot2 warning)
  if (has_ci && style == "line") {
    p <- p + ggplot2::labs(fill = "Group")
  }
  p
}


# ============================================================================ #
# type = "ineq" and "ineq.group" helpers
# ============================================================================ #

#' Resolve stats list: strings to registry functions, keep functions as-is
#' @keywords internal
.resolve_stats <- function(stats) {
  resolved <- list()
  for (i in seq_along(stats)) {
    s <- stats[[i]]
    if (is.character(s)) {
      if (!s %in% names(.ineq_stat_registry)) {
        stop("Unknown inequality stat '", s, "'. ",
             "Available: ", paste(names(.ineq_stat_registry), collapse = ", "))
      }
      resolved[[i]] <- list(name = s, fn = .ineq_stat_registry[[s]])
    } else if (is.function(s)) {
      nm <- names(stats)[i]
      if (is.null(nm) || nm == "") {
        nm <- paste0("custom_", i)
      }
      resolved[[i]] <- list(name = nm, fn = s)
    } else {
      stop("Each element of 'stats' must be a character string or function")
    }
  }
  resolved
}

#' @keywords internal
.plot_desc_ineq <- function(x, stats, ci = FALSE, style = "line") {
  if (is.null(x$raw_data)) {
    stop("Raw data not available. Re-run ineqx() to store individual-level data.")
  }

  resolved <- .resolve_stats(stats)
  d <- x$raw_data
  time_levels <- sort(unique(d$time))

  rows <- list()
  for (t in time_levels) {
    idx <- d$time == t
    y_t <- d$y[idx]
    w_t <- d$w[idx]
    w_arg <- if (all(w_t == 1)) NULL else w_t
    for (s in resolved) {
      rows[[length(rows) + 1]] <- data.frame(
        time = t,
        stat = s$name,
        value = s$fn(y_t, w = w_arg),
        stringsAsFactors = FALSE
      )
    }
  }

  plot_df <- do.call(rbind, rows)
  stat_names <- vapply(resolved, `[[`, character(1), "name")
  plot_df$stat <- factor(plot_df$stat, levels = stat_names)

  # CIs
  ci_info <- .resolve_ci(ci, x)
  ci_z <- stats::qnorm(0.975)
  if (ci_info$method == "delta") {
    ci_rows <- list()
    for (t in time_levels) {
      idx   <- d$time == t
      y_t   <- d$y[idx]
      w_t   <- d$w[idx]
      w_arg <- if (all(w_t == 1)) NULL else w_t
      w_se  <- if (is.null(w_arg)) rep(1, length(y_t)) else w_arg
      for (s in resolved) {
        se_fn <- .ineq_se_registry[[s$name]]
        val   <- plot_df$value[plot_df$time == t & plot_df$stat == s$name]
        if (is.null(se_fn)) {
          warning("Delta CI not available for custom stat '", s$name,
                  "'; omitting CI for this statistic.")
          next
        }
        se <- se_fn(y_t, w = w_se, val = val)
        ci_rows[[length(ci_rows) + 1]] <- data.frame(
          time = t, stat = s$name,
          ymin = val - ci_z * se, ymax = val + ci_z * se,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(ci_rows) > 0)
      plot_df <- merge(plot_df, do.call(rbind, ci_rows), by = c("time", "stat"))
  } else if (ci_info$method == "boot" && !is.null(x$raw_data)) {
    B <- ci_info$boot$B
    ci_df <- .bootstrap_desc_ci(d, B = B, compute_fn = function(d_boot) {
      time_lvls <- sort(unique(d_boot$time))
      r <- list()
      for (tt in time_lvls) {
        ii <- d_boot$time == tt
        y_tt <- d_boot$y[ii]
        w_tt <- d_boot$w[ii]
        w_a <- if (all(w_tt == 1)) NULL else w_tt
        for (ss in resolved) {
          r[[length(r) + 1]] <- data.frame(
            time = tt, stat = ss$name,
            value = ss$fn(y_tt, w = w_a),
            stringsAsFactors = FALSE)
        }
      }
      do.call(rbind, r)
    })
    plot_df <- merge(plot_df, ci_df, by = c("time", "stat"))
  }

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              group = 1L))

  has_ci <- "ymin" %in% names(plot_df)
  if (style == "point") {
    if (has_ci) {
      pos <- .dodge_pos(plot_df$time)
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4)
    }
    p <- p + ggplot2::geom_point(size = 2)
  } else {
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        alpha = 0.15)
    }
    p <- p + ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point()
  }

  p + ggplot2::facet_wrap(~ .data$stat, scales = "free_y") +
    ggplot2::labs(x = "Time", y = "") +
    theme_ineqx
}

#' @keywords internal
.plot_desc_ineq_group <- function(x, stats, ci = FALSE, style = "line") {
  if (is.null(x$raw_data)) {
    stop("Raw data not available. Re-run ineqx() to store individual-level data.")
  }

  resolved <- .resolve_stats(stats)
  d <- x$raw_data
  time_levels <- sort(unique(d$time))
  group_levels <- sort(unique(d$group))

  rows <- list()
  for (t in time_levels) {
    for (g in group_levels) {
      idx <- d$time == t & d$group == g
      y_tg <- d$y[idx]
      w_tg <- d$w[idx]
      w_arg <- if (all(w_tg == 1)) NULL else w_tg
      for (s in resolved) {
        rows[[length(rows) + 1]] <- data.frame(
          time = t,
          group = g,
          stat = s$name,
          value = s$fn(y_tg, w = w_arg),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  plot_df <- do.call(rbind, rows)
  stat_names <- vapply(resolved, `[[`, character(1), "name")
  plot_df$stat <- factor(plot_df$stat, levels = stat_names)
  plot_df$group <- factor(plot_df$group)

  # CIs
  ci_info <- .resolve_ci(ci, x)
  ci_z <- stats::qnorm(0.975)
  if (ci_info$method == "delta") {
    ci_rows <- list()
    for (t in time_levels) {
      for (g in group_levels) {
        idx   <- d$time == t & d$group == g
        y_tg  <- d$y[idx]
        w_tg  <- d$w[idx]
        w_arg <- if (all(w_tg == 1)) NULL else w_tg
        w_se  <- if (is.null(w_arg)) rep(1, length(y_tg)) else w_arg
        for (s in resolved) {
          se_fn <- .ineq_se_registry[[s$name]]
          val   <- plot_df$value[
            plot_df$time == t & plot_df$group == g & plot_df$stat == s$name
          ]
          if (is.null(se_fn)) {
            warning("Delta CI not available for custom stat '", s$name,
                    "'; omitting CI for this statistic.")
            next
          }
          se <- se_fn(y_tg, w = w_se, val = val)
          ci_rows[[length(ci_rows) + 1]] <- data.frame(
            time = t, group = g, stat = s$name,
            ymin = val - ci_z * se, ymax = val + ci_z * se,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(ci_rows) > 0)
      plot_df <- merge(plot_df, do.call(rbind, ci_rows),
                       by = c("time", "group", "stat"))
  } else if (ci_info$method == "boot" && !is.null(x$raw_data)) {
    B <- ci_info$boot$B
    ci_df <- .bootstrap_desc_ci(d, B = B, compute_fn = function(d_boot) {
      time_lvls <- sort(unique(d_boot$time))
      grp_lvls <- sort(unique(d_boot$group))
      r <- list()
      for (tt in time_lvls) {
        for (gg in grp_lvls) {
          ii <- d_boot$time == tt & d_boot$group == gg
          y_tg <- d_boot$y[ii]
          w_tg <- d_boot$w[ii]
          w_a <- if (all(w_tg == 1)) NULL else w_tg
          for (ss in resolved) {
            r[[length(r) + 1]] <- data.frame(
              time = tt, group = gg, stat = ss$name,
              value = ss$fn(y_tg, w = w_a),
              stringsAsFactors = FALSE)
          }
        }
      }
      do.call(rbind, r)
    })
    plot_df <- merge(plot_df, ci_df,
                     by = c("time", "group", "stat"))
  }

  plot_df$time <- .time_to_factor(plot_df$time)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                              color = .data$group,
                                              group = .data$group))

  has_ci <- "ymin" %in% names(plot_df)
  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        width = pos$width * 0.4, position = pos$pos)
    }
    p <- p + ggplot2::geom_point(position = pos$pos)
  } else {
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$group),
        alpha = 0.15, color = NA)
    }
    p <- p + ggplot2::geom_line() + ggplot2::geom_point()
  }

  p <- p + ggplot2::facet_wrap(~ .data$stat, scales = "free_y") +
    ggplot2::labs(x = "Time", y = "", color = "Group") +
    theme_ineqx

  if (has_ci && style == "line") {
    p <- p + ggplot2::labs(fill = "Group")
  }
  p
}


# ============================================================================ #
# Style helper: dodge position for "point" style
# ============================================================================ #

#' Convert time column to ordered factor for discrete x-axis
#'
#' Ensures time points are displayed as evenly-spaced categories, which avoids
#' large gaps when counterfactual time values (e.g. 0) are mixed with calendar
#' years (e.g. 1980, 1985).
#' @keywords internal
.time_to_factor <- function(time_vals) {
  lvls <- as.character(sort(unique(as.numeric(as.character(time_vals)))))
  factor(time_vals, levels = lvls)
}

#' Compute dodge position from time values
#' @keywords internal
.dodge_pos <- function(time_vals) {
  list(
    width = 0.6,
    pos = ggplot2::position_dodge(width = 0.6)
  )
}


# ============================================================================ #
# Bootstrap CI helper for descriptive plots
# ============================================================================ #

#' Bootstrap CIs for descriptive plot types
#'
#' Resamples rows of raw_data (stratified by time), applies compute_fn
#' to each replicate, and returns quantile-based CIs.
#'
#' @param raw_data Data frame with columns y, group, time, w
#' @param B Number of bootstrap replicates
#' @param compute_fn Function taking a resampled data frame and returning
#'   a data frame with at least a 'value' column plus grouping columns
#'   (time, Component, stat, group, etc.)
#' @param level Confidence level (default 0.95)
#' @return Data frame with grouping columns + ymin, ymax
#' @keywords internal
.bootstrap_desc_ci <- function(raw_data, B, compute_fn, level = 0.95) {
  alpha <- 1 - level
  time_levels <- sort(unique(raw_data$time))

  # Pre-split indices by time for stratified resampling
  time_idx <- lapply(time_levels, function(t) which(raw_data$time == t))
  names(time_idx) <- as.character(time_levels)

  # Progress bar
  show_progress <- interactive() && B > 1
  if (show_progress) {
    cat(sprintf("Bootstrap CI (B=%d): ", B))
    bar_width <- 40
    progress_at <- unique(round(seq(1, B, length.out = bar_width + 1)))
  }

  # Collect all replicates
  reps <- vector("list", B)
  for (b in seq_len(B)) {
    # Stratified resample: within each time period
    boot_rows <- unlist(lapply(time_idx, function(idx) {
      sample(idx, length(idx), replace = TRUE)
    }))
    d_boot <- raw_data[boot_rows, , drop = FALSE]
    rep_df <- tryCatch(compute_fn(d_boot), error = function(e) NULL)
    if (!is.null(rep_df)) {
      rep_df$.rep <- b
      reps[[b]] <- rep_df
    }

    # Update progress bar
    if (show_progress && b %in% progress_at) {
      pct <- round(100 * b / B)
      filled <- round(bar_width * b / B)
      bar <- paste0(strrep("=", filled), strrep(" ", bar_width - filled))
      cat(sprintf("\rBootstrap CI (B=%d): [%s] %3d%%", B, bar, pct))
      utils::flush.console()
    }
  }
  if (show_progress) cat("\n")

  all_reps <- do.call(rbind, reps)
  if (is.null(all_reps) || nrow(all_reps) == 0) {
    return(data.frame())
  }

  # Identify grouping columns (everything except value and .rep)
  group_cols <- setdiff(names(all_reps), c("value", ".rep"))

  # Compute quantiles per group
  ci_list <- list()
  if (length(group_cols) > 0) {
    groups <- unique(all_reps[, group_cols, drop = FALSE])
    for (i in seq_len(nrow(groups))) {
      mask <- rep(TRUE, nrow(all_reps))
      for (gc in group_cols) {
        mask <- mask & all_reps[[gc]] == groups[i, gc]
      }
      vals <- all_reps$value[mask]
      q <- stats::quantile(vals, probs = c(alpha / 2, 1 - alpha / 2),
                            na.rm = TRUE)
      row <- groups[i, , drop = FALSE]
      row$ymin <- q[1]
      row$ymax <- q[2]
      ci_list[[i]] <- row
    }
  }

  ci_df <- do.call(rbind, ci_list)
  rownames(ci_df) <- NULL
  ci_df
}


#' Recompute wibe summary from raw data (for bootstrap)
#' @keywords internal
.compute_wibe_from_d <- function(d) {
  time_levels <- sort(unique(d$time))
  group_levels <- sort(unique(d$group))

  wibe_list <- list()
  for (t in time_levels) {
    for (g in group_levels) {
      idx <- d$time == t & d$group == g
      yi <- d$y[idx]
      wi <- d$w[idx]
      ni <- sum(idx)

      if (ni == 0) {
        wibe_list[[length(wibe_list) + 1]] <- data.frame(
          time = t, group = g, n_raw = 0, sw = 0,
          mu = NA_real_, sigma2 = NA_real_,
          stringsAsFactors = FALSE)
        next
      }

      sw <- sum(wi)
      mu_g <- sum(wi * yi) / sw
      sigma2_g <- if (ni > 1) ni / (sw * (ni - 1)) * sum(wi * (yi - mu_g)^2) else 0

      wibe_list[[length(wibe_list) + 1]] <- data.frame(
        time = t, group = g, n_raw = ni, sw = sw,
        mu = mu_g, sigma2 = sigma2_g,
        stringsAsFactors = FALSE)
    }
  }

  wibe_df <- do.call(rbind, wibe_list)
  wibe_df$sigma <- sqrt(wibe_df$sigma2)

  # Compute n and pi
  wibe_out_list <- list()
  for (t in time_levels) {
    idx <- wibe_df$time == t
    sub <- wibe_df[idx, ]
    total_sw <- sum(sub$sw, na.rm = TRUE)
    total_n <- sum(sub$n_raw, na.rm = TRUE)
    sub$n <- (sub$sw / total_sw) * total_n
    sub$pi <- sub$n / sum(sub$n)
    wibe_out_list[[length(wibe_out_list) + 1]] <- sub
  }

  wibe_df <- do.call(rbind, wibe_out_list)
  wibe_df[, c("time", "group", "n", "pi", "mu", "sigma", "sigma2")]
}
