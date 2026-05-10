# ============================================================================ #
# Compare multiple ineqx result objects
# ============================================================================ #

#' Compare multiple ineqx results
#'
#' Stacks two or more ineqx result objects into a single comparison object
#' with a scenario column. Useful for comparing observed vs counterfactual
#' decompositions, or multiple counterfactual scenarios.
#'
#' @param ... Named ineqx result objects, or a single named list of them.
#'   All objects must be the same class (e.g., all \code{ineqx_causal_longit}).
#' @param names Optional character vector of scenario labels. If provided,
#'   overrides names from \code{...}.
#'
#' @return An object of class \code{ineqx_compare} containing:
#' \describe{
#'   \item{data}{Long-format data.frame with a \code{scenario} column}
#'   \item{class_type}{Character: the shared class of the input objects}
#'   \item{scenarios}{Character vector of scenario names}
#'   \item{ystat}{Character: shared inequality measure}
#'   \item{ref}{Numeric: shared reference period (for longitudinal)}
#'   \item{objects}{List of original input objects}
#' }
#'
#' @examples
#' \dontrun{
#' comp <- compare(Observed = res1, Counterfactual = res2)
#' print(comp)
#' plot(comp)
#' }
#'
#' @export
compare <- function(..., names = NULL) {
  dots <- list(...)


  # Allow single named list

  if (length(dots) == 1 && is.list(dots[[1]]) && !inherits(dots[[1]], "ineqx_desc") &&
      !inherits(dots[[1]], "ineqx_causal_cross") && !inherits(dots[[1]], "ineqx_causal_longit")) {
    dots <- dots[[1]]
  }

  if (length(dots) < 2) {
    stop("compare() requires at least 2 ineqx result objects.")
  }

  # Apply names
  if (!is.null(names)) {
    if (length(names) != length(dots)) {
      stop("Length of 'names' must match the number of objects.")
    }
    names(dots) <- names
  }

  if (is.null(names(dots)) || any(names(dots) == "")) {
    stop("All objects must be named. Use named arguments or the 'names' parameter.")
  }

  scenarios <- names(dots)

  # Validate: same class
  classes <- vapply(dots, function(x) class(x)[1], character(1))
  if (length(unique(classes)) > 1) {
    stop("All objects must be the same class. Found: ",
         paste(unique(classes), collapse = ", "))
  }
  class_type <- classes[1]

  valid_classes <- c("ineqx_desc", "ineqx_causal_cross", "ineqx_causal_longit")
  if (!class_type %in% valid_classes) {
    stop("compare() only works with ineqx result objects (ineqx_desc, ",
         "ineqx_causal_cross, ineqx_causal_longit).")
  }

  # Validate: same ystat
  ystats <- vapply(dots, function(x) x$ystat, character(1))
  if (length(unique(ystats)) > 1) {
    stop("All objects must have the same ystat. Found: ",
         paste(unique(ystats), collapse = ", "))
  }

  # Validate: same ref for longitudinal
  ref <- NULL
  if (class_type == "ineqx_causal_longit") {
    refs <- vapply(dots, function(x) x$ref, numeric(1))
    if (length(unique(refs)) > 1) {
      warning("Objects have different reference periods: ",
              paste(unique(refs), collapse = ", "),
              ". Using ref = ", refs[1])
    }
    ref <- refs[1]
  } else if (class_type == "ineqx_desc") {
    refs <- vapply(dots, function(x) if (!is.null(x$ref)) x$ref else NA_real_, numeric(1))
    if (any(!is.na(refs))) ref <- refs[!is.na(refs)][1]
  }

  # Extract data
  extract_fn <- switch(class_type,
    ineqx_desc = .extract_desc_data,
    ineqx_causal_cross = .extract_cross_data,
    ineqx_causal_longit = .extract_longit_data
  )

  data_list <- mapply(extract_fn, dots, scenarios, SIMPLIFY = FALSE)
  data <- do.call(rbind, data_list)
  data$scenario <- factor(data$scenario, levels = scenarios)
  rownames(data) <- NULL

  structure(list(
    data       = data,
    class_type = class_type,
    scenarios  = scenarios,
    ystat      = ystats[1],
    ref        = ref,
    objects    = dots
  ), class = "ineqx_compare")
}


# ---------------------------------------------------------------------------- #
# Data extraction helpers
# ---------------------------------------------------------------------------- #

#' @keywords internal
.extract_longit_data <- function(x, scenario) {
  is_shapley <- identical(x$order, "shapley")

  # All component names (show-value keys) and their data fields
  comp_names <- c("behavioral", "compositional", "pretreatment", "total",
                  "delta_beta", "delta_lambda", "delta_pi_b", "delta_pi_w",
                  "delta_mu", "delta_sigma")
  data_fields <- c("Delta_behavioral", "Delta_compositional",
                    "Delta_pretreatment", "Delta_total",
                    "Delta_beta", "Delta_lambda", "Delta_pi_B", "Delta_pi_W",
                    "Delta_mu", "Delta_sigma")

  if (is_shapley && !is.null(x$shapley)) {
    sh <- x$shapley
    n_comp <- length(comp_names)
    df <- data.frame(
      scenario  = scenario,
      time      = rep(sh$time, n_comp),
      component = rep(comp_names, each = nrow(sh)),
      value     = unlist(lapply(data_fields, function(f) sh[[f]])),
      detail    = "four",
      stringsAsFactors = FALSE
    )
    # Add reference period as zero anchor (excluded from x$shapley by construction)
    ref_time <- x$ref
    if (!is.null(ref_time) && !ref_time %in% sh$time) {
      ref_rows <- data.frame(
        scenario  = scenario,
        time      = ref_time,
        component = comp_names,
        value     = 0,
        detail    = "four",
        stringsAsFactors = FALSE
      )
      df <- rbind(ref_rows, df)
    }
    df
  } else {
    times <- as.numeric(names(x$results))
    ref_time <- x$ref
    if (!is.null(ref_time) && !ref_time %in% times) {
      times <- sort(c(ref_time, times))
    }

    n_comp <- length(comp_names)
    vals <- matrix(0, nrow = length(times), ncol = n_comp)
    for (i in seq_along(times)) {
      r <- x$results[[as.character(times[i])]]
      if (!is.null(r)) {
        for (k in seq_len(n_comp)) {
          vals[i, k] <- r[[data_fields[k]]]
        }
      }
    }

    data.frame(
      scenario  = scenario,
      time      = rep(times, n_comp),
      component = rep(comp_names, each = length(times)),
      value     = as.vector(vals),
      detail    = "four",
      stringsAsFactors = FALSE
    )
  }
}

#' @keywords internal
.extract_cross_data <- function(x, scenario) {
  # Tau-level data
  tau_df <- data.frame(
    scenario  = scenario,
    component = c("tau_B", "tau_W", "tau_total"),
    value     = c(x$tau_B, x$tau_W, x$tau_total),
    se        = if (!is.null(x$se)) c(x$se$se_tau_B, x$se$se_tau_W, x$se$se_tau_total)
                else rep(NA_real_, 3),
    detail    = "tau",
    stringsAsFactors = FALSE
  )

  # Sub-component data
  comps <- x$components
  sub_df <- data.frame(
    scenario  = scenario,
    component = c("het_B", "cov_B", "het_W", "cov_W"),
    value     = c(comps$het_B, comps$cov_B, comps$het_W, comps$cov_W),
    se        = NA_real_,
    detail    = "subcomp",
    stringsAsFactors = FALSE
  )

  rbind(tau_df, sub_df)
}

#' @keywords internal
.extract_desc_data <- function(x, scenario) {
  # Totals (wibe levels)
  totals <- x$totals
  if (x$ystat %in% c("Var", "VL")) {
    wibe_df <- data.frame(
      scenario  = scenario,
      time      = rep(totals$time, 3),
      component = rep(c("Within", "Between", "Total"), each = nrow(totals)),
      value     = c(totals$VarW, totals$VarB, totals$VarT),
      detail    = "wibe",
      stringsAsFactors = FALSE
    )
  } else {
    wibe_df <- data.frame(
      scenario  = scenario,
      time      = rep(totals$time, 3),
      component = rep(c("Within", "Between", "Total"), each = nrow(totals)),
      value     = c(totals$CV2W, totals$CV2B, totals$CV2T),
      detail    = "wibe",
      stringsAsFactors = FALSE
    )
  }

  # Deltas (if ref exists)
  if (!is.null(x$deltas)) {
    deltas <- x$deltas
    delta_df <- data.frame(
      scenario  = scenario,
      time      = rep(deltas$time, 4),
      component = rep(c("Between-group (mu)", "Within-group (sigma)",
                        "Compositional (pi)", "Total"),
                      each = nrow(deltas)),
      value     = c(deltas$delta_mu, deltas$delta_sigma,
                    deltas$delta_pi, deltas$delta_T),
      detail    = "decomp",
      stringsAsFactors = FALSE
    )
    return(rbind(wibe_df, delta_df))
  }

  wibe_df
}


# ---------------------------------------------------------------------------- #
# Print method
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_compare <- function(x, ...) {
  cat(sprintf("Comparison of %d scenarios (%s, %s)\n",
              length(x$scenarios), x$class_type, x$ystat))
  if (!is.null(x$ref)) cat("Reference:", x$ref, "\n")
  cat("Scenarios:", paste(x$scenarios, collapse = ", "), "\n\n")

  if (x$class_type == "ineqx_causal_longit") {
    .print_compare_longit(x)
  } else if (x$class_type == "ineqx_causal_cross") {
    .print_compare_cross(x)
  } else if (x$class_type == "ineqx_desc") {
    .print_compare_desc(x)
  }

  invisible(x)
}

#' @keywords internal
.print_compare_longit <- function(x) {
  d <- x$data[x$data$detail == "four", ]
  d <- d[d$component != "total", ]

  wide <- stats::reshape(d[, c("scenario", "time", "component", "value")],
                         idvar = c("time", "component"),
                         timevar = "scenario",
                         direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  wide <- wide[order(wide$time, wide$component), ]

  # Round numeric columns
  num_cols <- setdiff(names(wide), c("time", "component"))
  for (col in num_cols) wide[[col]] <- round(wide[[col]], 4)

  cat("Decomposition by time and scenario:\n")
  print(wide, row.names = FALSE)
}

#' @keywords internal
.print_compare_cross <- function(x) {
  d <- x$data[x$data$detail == "tau", ]
  wide <- stats::reshape(d[, c("scenario", "component", "value")],
                         idvar = "component",
                         timevar = "scenario",
                         direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))

  num_cols <- setdiff(names(wide), "component")
  for (col in num_cols) wide[[col]] <- round(wide[[col]], 4)

  cat("Treatment effects by scenario:\n")
  print(wide, row.names = FALSE)
}

#' @keywords internal
.print_compare_desc <- function(x) {
  d <- x$data[x$data$detail == "wibe", ]
  wide <- stats::reshape(d[, c("scenario", "time", "component", "value")],
                         idvar = c("time", "component"),
                         timevar = "scenario",
                         direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  wide <- wide[order(wide$time, wide$component), ]

  num_cols <- setdiff(names(wide), c("time", "component"))
  for (col in num_cols) wide[[col]] <- round(wide[[col]], 4)

  cat(sprintf("%s by time and scenario:\n", x$ystat))
  print(wide, row.names = FALSE)
}


# ---------------------------------------------------------------------------- #
# Plot method
# ---------------------------------------------------------------------------- #

#' Plot comparison of ineqx results
#'
#' Visualizes differences between scenarios. Color encodes the decomposition
#' component, linetype encodes the scenario.
#'
#' @param x An \code{ineqx_compare} object from \code{\link{compare}}.
#' @param type Character: plot type. Defaults depend on class:
#'   \itemize{
#'     \item \code{ineqx_causal_longit}: \code{"decomp"}, \code{"wibe"}
#'     \item \code{ineqx_causal_cross}: \code{"wibe"}
#'     \item \code{ineqx_desc}: \code{"wibe"}, \code{"decomp"}
#'   }
#' @param style Character: \code{"line"} or \code{"point"} (longitudinal only).
#' @param stats Character vector of decomposition components to display (for
#'   \code{type = "decomp"} with longitudinal data). See
#'   \code{\link{plot.ineqx_causal_longit}} for valid values. Default:
#'   \code{c("behavioral", "compositional", "pretreatment", "total")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object
#' @export
plot.ineqx_compare <- function(x, type = NULL, style = "line",
                                stats = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  # Default type per class
  if (is.null(type)) {
    type <- switch(x$class_type,
      ineqx_causal_longit = "decomp",
      ineqx_causal_cross  = "wibe",
      ineqx_desc          = "wibe"
    )
  }

  switch(x$class_type,
    ineqx_causal_longit = .plot_compare_longit(x, type = type, style = style,
                                                stats = stats),
    ineqx_causal_cross  = .plot_compare_cross(x, type = type),
    ineqx_desc          = .plot_compare_desc(x, type = type, style = style)
  )
}

#' @keywords internal
.plot_compare_longit <- function(x, type = "decomp", style = "line",
                                  stats = NULL) {
  style <- match.arg(style, c("line", "point"))

  if (type == "decomp") {
    if (is.null(stats)) {
      stats <- c("behavioral", "compositional", "pretreatment", "total")
    }
    registry <- .decomp_registry()
    invalid <- setdiff(stats, names(registry))
    if (length(invalid) > 0) {
      stop("Unknown stats values: ", paste(invalid, collapse = ", "))
    }

    comp_info <- registry[stats]
    labels <- vapply(comp_info, `[[`, character(1), "label")
    colors <- vapply(comp_info, `[[`, character(1), "color")
    names(colors) <- labels

    # Filter to requested components and map to display labels
    plot_df <- x$data[x$data$component %in% stats, ]
    label_map <- stats::setNames(labels, stats)
    plot_df$component <- label_map[plot_df$component]
    plot_df$component <- factor(plot_df$component, levels = rev(labels))

    y_label <- paste0("Change in ", x$ystat)
    title <- "Longitudinal Causal Decomposition"

  } else if (type == "wibe") {
    plot_df <- .extract_compare_wibe_longit(x)
    colors <- c("Total" = "black", "Between" = "#e85d00", "Within" = "#008b94")
    y_label <- paste0("Treatment effect on ", x$ystat)
    title <- "Within/Between Treatment Effects"

  } else {
    stop("Unknown plot type '", type,
         "'. Use 'decomp' or 'wibe'.")
  }

  plot_df$time <- .time_to_factor(plot_df$time)
  n_scenarios <- length(x$scenarios)
  lty_vals <- stats::setNames(seq_len(n_scenarios), x$scenarios)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data$time, y = .data$value,
    color = .data$component,
    linetype = .data$scenario,
    group = interaction(.data$component, .data$scenario)))

  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    p <- p + ggplot2::geom_point(size = 2, position = pos$pos)
  } else {
    p <- p + ggplot2::geom_line(linewidth = 1)
  }

  p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = colors, name = "Component") +
    ggplot2::scale_linetype_manual(values = lty_vals, name = "Scenario") +
    ggplot2::labs(x = "Time", y = y_label, title = title) +
    theme_ineqx
}

#' @keywords internal
.extract_compare_wibe_longit <- function(x) {
  rows <- list()
  for (scen in x$scenarios) {
    obj <- x$objects[[scen]]
    is_shapley <- identical(obj$order, "shapley")
    results <- if (is_shapley) obj$all_orderings[[1]]$results else obj$results

    non_ref_times <- as.numeric(names(results))
    ref_time <- obj$ref
    all_times <- sort(c(ref_time, non_ref_times))
    first_r <- results[[1]]

    for (i in seq_along(all_times)) {
      t <- all_times[i]
      if (t == ref_time) {
        tau_B <- first_r$tau_B_t0
        tau_W <- first_r$tau_W_t0
      } else {
        r <- results[[as.character(t)]]
        tau_B <- r$tau_B_t
        tau_W <- r$tau_W_t
      }
      rows[[length(rows) + 1]] <- data.frame(
        scenario  = scen,
        time      = t,
        component = c("Between", "Within", "Total"),
        value     = c(tau_B, tau_W, tau_B + tau_W),
        stringsAsFactors = FALSE
      )
    }
  }
  plot_df <- do.call(rbind, rows)
  plot_df$component <- factor(plot_df$component,
                              levels = c("Total", "Between", "Within"))
  plot_df$scenario <- factor(plot_df$scenario, levels = x$scenarios)
  plot_df
}

#' @keywords internal
.plot_compare_cross <- function(x, type = "wibe") {
  if (type != "wibe") {
    stop("Unknown plot type '", type, "'. Use 'wibe'.")
  }

  d <- x$data[x$data$detail == "tau", ]
  d$component <- factor(d$component,
                        levels = c("tau_total", "tau_B", "tau_W"),
                        labels = c("Total", "Between", "Within"))

  fill_vals <- c("Between" = "#e85d00", "Within" = "#008b94", "Total" = "#1a4e66")

  ggplot2::ggplot(d, ggplot2::aes(x = .data$component, y = .data$value,
                                   fill = .data$component,
                                   alpha = .data$scenario)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = fill_vals, name = "Component") +
    ggplot2::scale_alpha_manual(
      values = stats::setNames(seq(1, 0.4, length.out = length(x$scenarios)),
                               x$scenarios),
      name = "Scenario") +
    ggplot2::labs(x = "", y = paste0("Treatment effect on ", x$ystat),
                  title = "Treatment Effect on Inequality") +
    theme_ineqx
}

#' @keywords internal
.plot_compare_desc <- function(x, type = "wibe", style = "line") {
  style <- match.arg(style, c("line", "point"))

  if (type == "wibe") {
    plot_df <- x$data[x$data$detail == "wibe", ]
    plot_df$component <- factor(plot_df$component,
                                levels = c("Total", "Within", "Between"))
    colors <- c("Total" = "black", "Within" = "#008b94", "Between" = "#e85d00")
    y_label <- x$ystat
    title <- "Variance Decomposition"

  } else if (type == "decomp") {
    plot_df <- x$data[x$data$detail == "decomp", ]
    if (nrow(plot_df) == 0) {
      stop("No delta data available. Did you set a reference period?")
    }
    plot_df$component <- factor(plot_df$component,
                                levels = c("Total", "Between-group (mu)",
                                           "Within-group (sigma)",
                                           "Compositional (pi)"))
    colors <- c("Total" = "black", "Between-group (mu)" = "#e85d00",
                "Within-group (sigma)" = "#008b94", "Compositional (pi)" = "#7b2d8e")
    y_label <- paste0("Change in ", x$ystat)
    title <- "Descriptive Decomposition of Change"

  } else {
    stop("Unknown plot type '", type, "'. Use 'wibe' or 'decomp'.")
  }

  plot_df$time <- .time_to_factor(plot_df$time)
  n_scenarios <- length(x$scenarios)
  lty_vals <- stats::setNames(seq_len(n_scenarios), x$scenarios)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data$time, y = .data$value,
    color = .data$component,
    linetype = .data$scenario,
    group = interaction(.data$component, .data$scenario)))

  if (style == "point") {
    pos <- .dodge_pos(plot_df$time)
    p <- p + ggplot2::geom_point(size = 2, position = pos$pos)
  } else {
    p <- p + ggplot2::geom_line(linewidth = 1)
  }

  p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = colors, name = "Component") +
    ggplot2::scale_linetype_manual(values = lty_vals, name = "Scenario") +
    ggplot2::labs(x = "Time", y = y_label, title = title) +
    theme_ineqx
}
