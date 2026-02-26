# ============================================================================ #
# ineq: Descriptive within/between variance decomposition
# ============================================================================ #

#' Descriptive variance decomposition
#'
#' Decomposes total variance (or CV^2) into within-group and between-group
#' components. For longitudinal data with a reference period, further decomposes
#' the change over time into components attributable to changing means (mu),
#' dispersions (sigma), and group composition (pi), using sequential parameter
#' switching with optional Shapley averaging.
#'
#' @param y Character, name of the outcome variable in \code{data}
#' @param group Character, name of the grouping variable in \code{data}
#' @param time Character, name of the time variable in \code{data}. If NULL,
#'   a single cross-section is assumed.
#' @param weights Character, name of the weight variable in \code{data}. If NULL,
#'   equal weights are used.
#' @param data A data.frame containing the variables
#' @param ref Numeric, reference time period for computing deltas. If NULL,
#'   deltas are not computed.
#' @param ystat Character, either \code{"Var"} (default) or \code{"CV2"}.
#' @param order Decomposition ordering for the counterfactual deltas. Either:
#'   \itemize{
#'     \item A character vector of length 3: a permutation of
#'       \code{c("mu", "sigma", "pi")} for a single ordering.
#'     \item \code{"shapley"} (default): averages across all 6 possible orderings.
#'   }
#'   Only used when \code{time} and \code{ref} are provided.
#'
#' @return An object of class \code{"ineqx_desc"} containing:
#' \describe{
#'   \item{wibe}{data.frame by group and time: n, pi, mu, sigma, sigma2}
#'   \item{totals}{data.frame by time: W, B, Total (and CV2 variants)}
#'   \item{deltas}{data.frame by time: counterfactual delta components
#'     relative to ref (NULL if ref is not specified)}
#'   \item{ystat}{The inequality measure used}
#'   \item{order}{The ordering used for counterfactual decomposition}
#' }
#'
#' @examples
#' data(incdat)
#' # Single cross-section
#' ineq(y = "inc", group = "group", data = incdat[incdat$t == 0, ])
#'
#' # Over time with reference period
#' ineq(y = "inc", group = "group", time = "year", data = incdat,
#'      ref = 1, ystat = "Var")
#'
#' @export
ineq <- function(y, group, time = NULL, weights = NULL,
                 data, ref = NULL, ystat = "Var", order = "shapley") {

  ystat <- match.arg(ystat, c("Var", "CV2"))

  # Validate inputs
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!y %in% names(data)) stop("'", y, "' not found in data")
  if (!group %in% names(data)) stop("'", group, "' not found in data")
  if (!is.null(time) && !time %in% names(data)) stop("'", time, "' not found in data")
  if (!is.null(weights) && !weights %in% names(data)) stop("'", weights, "' not found in data")

  # Extract and rename columns for internal use
  d <- data.frame(
    y = data[[y]],
    group = data[[group]],
    stringsAsFactors = FALSE
  )
  if (!is.null(time)) {
    d$time <- data[[time]]
  } else {
    d$time <- 1L
  }
  if (!is.null(weights)) {
    d$w <- data[[weights]]
  } else {
    d$w <- 1
  }

  # Drop NAs
  d <- d[complete.cases(d), ]
  if (nrow(d) == 0) stop("No complete cases in data")

  # Get levels
  group_levels <- sort(unique(d$group))
  time_levels <- sort(unique(d$time))

  # Validate ref
  if (!is.null(ref) && !(ref %in% time_levels)) {
    stop("'ref' = ", ref, " not found in time variable. ",
         "Available: ", paste(time_levels, collapse = ", "))
  }

  # ------------------------------------------------------------ #
  # Compute group-level statistics (wibe)
  # ------------------------------------------------------------ #

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
          stringsAsFactors = FALSE
        )
        next
      }

      sw <- sum(wi)
      mu_g <- sum(wi * yi) / sw

      if (ni > 1) {
        sigma2_g <- ni / (sw * (ni - 1)) * sum(wi * (yi - mu_g)^2)
      } else {
        sigma2_g <- 0
      }

      wibe_list[[length(wibe_list) + 1]] <- data.frame(
        time = t, group = g, n_raw = ni, sw = sw,
        mu = mu_g, sigma2 = sigma2_g,
        stringsAsFactors = FALSE
      )
    }
  }

  wibe_df <- do.call(rbind, wibe_list)
  wibe_df$sigma <- sqrt(wibe_df$sigma2)

  # Compute weighted n and pi per time period
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

  wibe_out <- wibe_df[, c("time", "group", "n", "pi", "mu", "sigma", "sigma2")]
  rownames(wibe_out) <- NULL

  # ------------------------------------------------------------ #
  # Compute totals by time
  # ------------------------------------------------------------ #

  totals_list <- list()
  for (t in time_levels) {
    idx <- wibe_out$time == t
    sub <- wibe_out[idx, ]
    n_vec <- sub$n
    mu_vec <- sub$mu
    sigma_vec <- sub$sigma

    W_var <- VarW_n(n_vec, sigma_vec)
    B_var <- VarB_n(n_vec, mu_vec)
    T_var <- W_var + B_var

    W_cv2 <- CV2W_n(n_vec, mu_vec, sigma_vec)
    B_cv2 <- CV2B_n(n_vec, mu_vec)
    T_cv2 <- W_cv2 + B_cv2

    N_total <- sum(n_vec)
    gmu <- sum((n_vec / N_total) * mu_vec)

    totals_list[[length(totals_list) + 1]] <- data.frame(
      time = t, N = N_total, grand_mean = gmu,
      VarW = W_var, VarB = B_var, VarT = T_var,
      CV2W = W_cv2, CV2B = B_cv2, CV2T = T_cv2,
      stringsAsFactors = FALSE
    )
  }

  totals <- do.call(rbind, totals_list)
  rownames(totals) <- NULL

  # ------------------------------------------------------------ #
  # Compute counterfactual deltas (change relative to ref)
  # ------------------------------------------------------------ #

  deltas <- NULL
  if (!is.null(ref) && length(time_levels) > 1) {
    deltas <- .compute_desc_deltas(wibe_out, time_levels, ref, ystat, order)
  }

  # ------------------------------------------------------------ #
  # Remove time column if single cross-section
  # ------------------------------------------------------------ #

  if (is.null(time)) {
    wibe_out$time <- NULL
    totals$time <- NULL
  }

  structure(
    list(
      wibe = wibe_out,
      totals = totals,
      deltas = deltas,
      ystat = ystat,
      order = order,
      ref = ref,
      call = match.call()
    ),
    class = "ineqx_desc"
  )
}


# ---------------------------------------------------------------------------- #
# Internal: Counterfactual delta decomposition
# ---------------------------------------------------------------------------- #

.compute_desc_deltas <- function(wibe_out, time_levels, ref, ystat, order) {

  use_shapley <- is.character(order) && length(order) == 1 && order == "shapley"

  if (!use_shapley) {
    valid_params <- c("mu", "sigma", "pi")
    order <- match.arg(order, valid_params, several.ok = TRUE)
    if (length(order) != 3 || !setequal(order, valid_params)) {
      stop("'order' must be 'shapley' or a permutation of c('mu', 'sigma', 'pi')")
    }
  }

  ref_wibe <- wibe_out[wibe_out$time == ref, ]
  pi_0 <- ref_wibe$pi
  mu_0 <- ref_wibe$mu
  sigma_0 <- ref_wibe$sigma

  # Inequality evaluation function
  .eval_ineq <- function(pi, mu, sigma) {
    if (ystat == "Var") {
      W <- VarW_pi(pi, sigma)
      B <- VarB_pi(pi, mu)
    } else {
      W <- CV2W_pi(pi, mu, sigma)
      B <- CV2B_pi(pi, mu)
    }
    list(W = W, B = B, Total = W + B)
  }

  # Single ordering decomposition
  .decompose_one_order <- function(pi_0, mu_0, sigma_0, pi_t, mu_t, sigma_t, ord) {
    param_map <- list(
      mu    = list(param = "mu",    from = mu_0,    to = mu_t),
      sigma = list(param = "sigma", from = sigma_0, to = sigma_t),
      pi    = list(param = "pi",    from = pi_0,    to = pi_t)
    )

    current <- list(pi = pi_0, mu = mu_0, sigma = sigma_0)
    prev <- .eval_ineq(current$pi, current$mu, current$sigma)
    baseline <- prev

    components <- list()
    for (step in ord) {
      p <- param_map[[step]]$param
      current[[p]] <- param_map[[step]]$to
      new_val <- .eval_ineq(current$pi, current$mu, current$sigma)
      components[[step]] <- list(
        delta_W = new_val$W - prev$W,
        delta_B = new_val$B - prev$B,
        delta_T = new_val$Total - prev$Total
      )
      prev <- new_val
    }
    components
  }

  deltas_list <- list()
  for (t in time_levels) {
    sub <- wibe_out[wibe_out$time == t, ]
    pi_t <- sub$pi
    mu_t <- sub$mu
    sigma_t <- sub$sigma

    if (use_shapley) {
      # Average over all 6 orderings
      perms <- .all_permutations(c("mu", "sigma", "pi"))

      # Accumulate
      accum <- list(
        mu    = list(delta_W = 0, delta_B = 0, delta_T = 0),
        sigma = list(delta_W = 0, delta_B = 0, delta_T = 0),
        pi    = list(delta_W = 0, delta_B = 0, delta_T = 0)
      )

      for (perm in perms) {
        comp <- .decompose_one_order(pi_0, mu_0, sigma_0, pi_t, mu_t, sigma_t, perm)
        for (p in c("mu", "sigma", "pi")) {
          accum[[p]]$delta_W <- accum[[p]]$delta_W + comp[[p]]$delta_W
          accum[[p]]$delta_B <- accum[[p]]$delta_B + comp[[p]]$delta_B
          accum[[p]]$delta_T <- accum[[p]]$delta_T + comp[[p]]$delta_T
        }
      }

      n_perms <- length(perms)
      for (p in c("mu", "sigma", "pi")) {
        accum[[p]]$delta_W <- accum[[p]]$delta_W / n_perms
        accum[[p]]$delta_B <- accum[[p]]$delta_B / n_perms
        accum[[p]]$delta_T <- accum[[p]]$delta_T / n_perms
      }
      components <- accum

    } else {
      components <- .decompose_one_order(pi_0, mu_0, sigma_0, pi_t, mu_t, sigma_t, order)
    }

    deltas_list[[length(deltas_list) + 1]] <- data.frame(
      time = t,
      # Per-parameter deltas
      delta_mu_W = components$mu$delta_W,
      delta_mu_B = components$mu$delta_B,
      delta_mu_T = components$mu$delta_T,
      delta_sigma_W = components$sigma$delta_W,
      delta_sigma_B = components$sigma$delta_B,
      delta_sigma_T = components$sigma$delta_T,
      delta_pi_W = components$pi$delta_W,
      delta_pi_B = components$pi$delta_B,
      delta_pi_T = components$pi$delta_T,
      # Aggregated
      delta_W = components$mu$delta_W + components$sigma$delta_W + components$pi$delta_W,
      delta_B = components$mu$delta_B + components$sigma$delta_B + components$pi$delta_B,
      delta_T = components$mu$delta_T + components$sigma$delta_T + components$pi$delta_T,
      stringsAsFactors = FALSE
    )
  }

  deltas <- do.call(rbind, deltas_list)
  rownames(deltas) <- NULL
  deltas
}


# ---------------------------------------------------------------------------- #
# Backward compatibility wrapper
# ---------------------------------------------------------------------------- #

#' @keywords internal
desc_decompose <- function(y, group, time = NULL, weights = NULL,
                            data, ref = NULL, ystat = "Var") {
  ineq(y = y, group = group, time = time, weights = weights,
       data = data, ref = ref, ystat = ystat, order = "shapley")
}
