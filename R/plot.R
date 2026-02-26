# ============================================================================ #
# Plot methods for ineqx result objects
# ============================================================================ #

#' @importFrom ggplot2 .data
NULL

# ---------------------------------------------------------------------------- #
# plot.ineqx_desc
# ---------------------------------------------------------------------------- #

#' Plot descriptive variance decomposition
#'
#' @param x An \code{ineqx_desc} object
#' @param type Character, plot type:
#'   \code{"wibe"} for within/between levels over time,
#'   \code{"deltas"} for change components over time (requires ref).
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot2 object
#' @export
plot.ineqx_desc <- function(x, type = "deltas", ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  if (type == "wibe") {
    .plot_desc_wibe(x)
  } else if (type == "deltas") {
    if (is.null(x$deltas)) {
      stop("Deltas not available. Set 'ref' in ineq() to compute deltas.")
    }
    .plot_desc_deltas(x)
  } else {
    stop("Unknown plot type '", type, "'. Use 'wibe' or 'deltas'.")
  }
}

#' @keywords internal
.plot_desc_wibe <- function(x) {
  totals <- x$totals
  if (!"time" %in% names(totals)) {
    stop("Time variable required for plotting")
  }

  if (x$ystat == "Var") {
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

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                         color = .data$Component,
                                         linetype = .data$Component)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(x = "Time", y = x$ystat, title = "Within/Between Decomposition") +
    ggplot2::scale_color_manual(values = c("Total" = "black",
                                            "Within" = "#008b94",
                                            "Between" = "#e85d00")) +
    theme_ineqx
}

#' @keywords internal
.plot_desc_deltas <- function(x) {
  deltas <- x$deltas

  # Plot per-parameter total deltas (mu, sigma, pi contributions)
  plot_df <- data.frame(
    time = rep(deltas$time, 4),
    Component = rep(c("Means (mu)", "Dispersions (sigma)",
                      "Composition (pi)", "Total"),
                    each = nrow(deltas)),
    value = c(deltas$delta_mu_T, deltas$delta_sigma_T,
              deltas$delta_pi_T, deltas$delta_T)
  )

  plot_df$Component <- factor(plot_df$Component,
                              levels = c("Total", "Means (mu)",
                                         "Dispersions (sigma)",
                                         "Composition (pi)"))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                         color = .data$Component,
                                         linetype = .data$Component)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Time",
                  y = paste0("Change in ", x$ystat, " (from ref = ", x$ref, ")"),
                  title = "Descriptive Decomposition of Change") +
    ggplot2::scale_color_manual(values = c("Total" = "black",
                                            "Means (mu)" = "#e85d00",
                                            "Dispersions (sigma)" = "#008b94",
                                            "Composition (pi)" = "#7b2d8e")) +
    ggplot2::scale_linetype_manual(values = c("Total" = "solid",
                                               "Means (mu)" = "solid",
                                               "Dispersions (sigma)" = "solid",
                                               "Composition (pi)" = "solid")) +
    theme_ineqx
}

# ---------------------------------------------------------------------------- #
# plot.ineqx_causal_cross
# ---------------------------------------------------------------------------- #

#' Plot cross-sectional causal decomposition
#'
#' @param x An \code{ineqx_causal_cross} object
#' @param type Character: \code{"bar"} for bar chart, \code{"group"} for
#'   group-level contributions.
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot2 object
#' @export
plot.ineqx_causal_cross <- function(x, type = "bar", ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  if (type == "bar") {
    plot_df <- data.frame(
      Component = factor(c("Between", "Within", "Total"),
                         levels = c("Total", "Between", "Within")),
      value = c(x$delta_B, x$delta_W, x$delta_total)
    )

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Component, y = .data$value,
                                           fill = .data$Component)) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "", y = paste0("Treatment effect on ", x$ystat),
                    title = "Treatment Effect on Inequality") +
      ggplot2::scale_fill_manual(values = c("Total" = "black",
                                             "Between" = "#e85d00",
                                             "Within" = "#008b94")) +
      ggplot2::guides(fill = "none") +
      theme_ineqx
  } else if (type == "group") {
    bg <- x$by_group
    plot_df <- data.frame(
      group = rep(bg$group, 2),
      Component = rep(c("Within", "Between"), each = nrow(bg)),
      value = c(bg$contrib_W,
                bg$pi * bg$beta^2)  # simplified group-level B contribution
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
  } else {
    stop("Unknown plot type '", type, "'. Use 'bar' or 'group'.")
  }
}

# ---------------------------------------------------------------------------- #
# plot.ineqx_causal_longit
# ---------------------------------------------------------------------------- #

#' Plot longitudinal causal decomposition
#'
#' @param x An \code{ineqx_causal_longit} object
#' @param type Character:
#'   \code{"four"} for four-component decomposition,
#'   \code{"six"} for six-component decomposition,
#'   \code{"wb"} for within/between totals.
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot2 object
#' @export
plot.ineqx_causal_longit <- function(x, type = "four", ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  # Build data from results
  times <- as.numeric(names(x$results))

  if (type == "four") {
    vals <- t(vapply(x$results, function(r) {
      c(r$Delta_behavioral, r$Delta_compositional,
        r$Delta_pretreatment, r$Delta_total)
    }, numeric(4)))

    plot_df <- data.frame(
      time = rep(times, 4),
      Component = rep(c("Behavioral", "Compositional",
                        "Pre-treatment", "Total"), each = length(times)),
      value = as.vector(vals)
    )
    plot_df$Component <- factor(plot_df$Component,
                                levels = c("Total", "Behavioral",
                                           "Compositional", "Pre-treatment"))

    colors <- c("Total" = "black", "Behavioral" = "#008b94",
                "Compositional" = "#7b2d8e", "Pre-treatment" = "#e85d00")

  } else if (type == "six") {
    vals <- t(vapply(x$results, function(r) {
      c(r$Delta_beta, r$Delta_lambda, r$Delta_pi_B, r$Delta_pi_W,
        r$Delta_mu, r$Delta_sigma)
    }, numeric(6)))

    comp_names <- c("Delta_beta", "Delta_lambda", "Delta_pi_B",
                    "Delta_pi_W", "Delta_mu", "Delta_sigma")
    plot_df <- data.frame(
      time = rep(times, 6),
      Component = rep(comp_names, each = length(times)),
      value = as.vector(vals)
    )
    plot_df$Component <- factor(plot_df$Component, levels = comp_names)

    colors <- c("Delta_beta" = "#008b94", "Delta_lambda" = "#00bfc4",
                "Delta_pi_B" = "#7b2d8e", "Delta_pi_W" = "#b47cc5",
                "Delta_mu" = "#e85d00", "Delta_sigma" = "#f5a623")

  } else if (type == "wb") {
    vals <- t(vapply(x$results, function(r) {
      c(r$Delta_B, r$Delta_W, r$Delta_total)
    }, numeric(3)))

    plot_df <- data.frame(
      time = rep(times, 3),
      Component = rep(c("Between", "Within", "Total"),
                      each = length(times)),
      value = as.vector(vals)
    )
    plot_df$Component <- factor(plot_df$Component,
                                levels = c("Total", "Between", "Within"))

    colors <- c("Total" = "black", "Between" = "#e85d00",
                "Within" = "#008b94")

  } else {
    stop("Unknown plot type '", type, "'. Use 'four', 'six', or 'wb'.")
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                         color = .data$Component)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Time",
                  y = paste0("Change in treatment effect on ", x$ystat),
                  title = "Longitudinal Causal Decomposition") +
    ggplot2::scale_color_manual(values = colors) +
    theme_ineqx
}

# ---------------------------------------------------------------------------- #
# plot.ineqx_shapley
# ---------------------------------------------------------------------------- #

#' Plot Shapley values
#'
#' @param x An \code{ineqx_shapley} object
#' @param type Character: \code{"four"} or \code{"six"}
#' @param ... Additional arguments (currently unused)
#'
#' @return A ggplot2 object
#' @export
plot.ineqx_shapley <- function(x, type = "four", ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  s <- x$shapley

  if (type == "four") {
    comp_names <- c("Behavioral", "Compositional", "Pre-treatment")
    plot_df <- data.frame(
      time = rep(s$time, 3),
      Component = rep(comp_names, each = nrow(s)),
      value = c(s$Delta_behavioral, s$Delta_compositional,
                s$Delta_pretreatment)
    )
    plot_df$Component <- factor(plot_df$Component, levels = comp_names)

    colors <- c("Behavioral" = "#008b94", "Compositional" = "#7b2d8e",
                "Pre-treatment" = "#e85d00")
  } else if (type == "six") {
    comp_names <- c("Delta_beta", "Delta_lambda", "Delta_pi_B",
                    "Delta_pi_W", "Delta_mu", "Delta_sigma")
    plot_df <- data.frame(
      time = rep(s$time, 6),
      Component = rep(comp_names, each = nrow(s)),
      value = c(s$Delta_beta, s$Delta_lambda, s$Delta_pi_B,
                s$Delta_pi_W, s$Delta_mu, s$Delta_sigma)
    )
    plot_df$Component <- factor(plot_df$Component, levels = comp_names)

    colors <- c("Delta_beta" = "#008b94", "Delta_lambda" = "#00bfc4",
                "Delta_pi_B" = "#7b2d8e", "Delta_pi_W" = "#b47cc5",
                "Delta_mu" = "#e85d00", "Delta_sigma" = "#f5a623")
  } else {
    stop("Unknown plot type '", type, "'. Use 'four' or 'six'.")
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$value,
                                         color = .data$Component)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Time",
                  y = paste0("Shapley value (", x$ystat, ")"),
                  title = "Shapley Values (Averaged Across Orderings)") +
    ggplot2::scale_color_manual(values = colors) +
    theme_ineqx
}

