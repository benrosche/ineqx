# ============================================================================ #
# Internal: Compute group-level within/between statistics from raw data
# Shared helper used by .ineq_descriptive() and .ineq_descriptive_with_ref()
# ============================================================================ #

.compute_wibe <- function(y, group, time = NULL, weights = NULL, data) {

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

  # Compute group-level statistics
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

  list(wibe_out = wibe_out, time_levels = time_levels, group_levels = group_levels,
       raw_data = d)
}


# ============================================================================ #
# Internal: Compute totals from wibe (n-based, from raw data)
# ============================================================================ #

.compute_totals_n <- function(wibe_out, time_levels) {
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
  totals
}


# ============================================================================ #
# Internal: Compute totals from wibe (pi-based, from counterfactual params)
# ============================================================================ #

.compute_totals_pi <- function(wibe_out, time_levels) {
  totals_list <- list()
  for (t in time_levels) {
    idx <- wibe_out$time == t
    sub <- wibe_out[idx, ]
    pi_vec <- sub$pi
    mu_vec <- sub$mu
    sigma_vec <- sub$sigma

    W_var <- VarW_pi(pi_vec, sigma_vec)
    B_var <- VarB_pi(pi_vec, mu_vec)
    T_var <- W_var + B_var

    W_cv2 <- CV2W_pi(pi_vec, mu_vec, sigma_vec)
    B_cv2 <- CV2B_pi(pi_vec, mu_vec)
    T_cv2 <- W_cv2 + B_cv2

    N_total <- sum(sub$n)
    gmu <- sum(pi_vec * mu_vec)

    totals_list[[length(totals_list) + 1]] <- data.frame(
      time = t, N = N_total, grand_mean = gmu,
      VarW = W_var, VarB = B_var, VarT = T_var,
      CV2W = W_cv2, CV2B = B_cv2, CV2T = T_cv2,
      stringsAsFactors = FALSE
    )
  }

  totals <- do.call(rbind, totals_list)
  rownames(totals) <- NULL
  totals
}


# ============================================================================ #
# Internal: Descriptive within/between variance decomposition (from raw data)
# Called by ineqx() when treat is NULL and params is NULL
# ============================================================================ #

.ineq_descriptive <- function(y, group, time = NULL, weights = NULL,
                               data, ref = NULL, ystat = "Var",
                               order = "shapley") {

  # Compute group-level statistics
  wibe_result <- .compute_wibe(y, group, time, weights, data)
  wibe_out <- wibe_result$wibe_out
  time_levels <- wibe_result$time_levels

  # Validate ref
  if (!is.null(ref) && !(ref %in% time_levels)) {
    stop("'ref' = ", ref, " not found in time variable. ",
         "Available: ", paste(time_levels, collapse = ", "))
  }

  # Compute totals by time
  totals <- .compute_totals_n(wibe_out, time_levels)

  # Compute counterfactual deltas (change relative to ref)
  deltas <- NULL
  if (!is.null(ref) && length(time_levels) > 1) {
    deltas <- .compute_desc_deltas(wibe_out, time_levels, ref, ystat, order)
  }

  # Remove time column if single cross-section
  if (is.null(time)) {
    wibe_out$time <- NULL
    totals$time <- NULL
  }

  structure(
    list(
      wibe = wibe_out,
      totals = totals,
      deltas = deltas,
      raw_data = wibe_result$raw_data,
      ystat = ystat,
      order = order,
      ref = ref,
      call = match.call()
    ),
    class = "ineqx_desc"
  )
}


# ============================================================================ #
# Internal: Descriptive decomposition with counterfactual reference params
# Called by ineqx() when params is an ineqx_desc_params object
# ============================================================================ #

.ineq_descriptive_with_ref <- function(y, group, time = NULL, weights = NULL,
                                        data, ref_params, ref,
                                        ystat = "Var", order = "shapley") {

  # Compute observed group-level statistics from raw data
  wibe_result <- .compute_wibe(y, group, time, weights, data)
  obs_wibe <- wibe_result$wibe_out
  obs_time_levels <- wibe_result$time_levels
  obs_groups <- wibe_result$group_levels

  # Validate group alignment
  ref_groups <- ref_params$groups
  if (!setequal(obs_groups, ref_groups)) {
    stop("Groups in 'params' (", paste(ref_groups, collapse = ", "),
         ") do not match groups in data (", paste(obs_groups, collapse = ", "), ")")
  }

  # Build wibe rows from params
  p_data <- ref_params$data
  p_data <- p_data[order(p_data$group), ]

  if (is.null(ref_params$times)) {
    # Single cross-section params: auto-assign time=0
    params_times <- 0L
    p_data$time <- 0L
  } else {
    params_times <- ref_params$times
  }

  params_wibe <- data.frame(
    time   = p_data$time,
    group  = p_data$group,
    n      = p_data$pi,
    pi     = p_data$pi,
    mu     = p_data$mu,
    sigma  = p_data$sigma,
    sigma2 = p_data$sigma^2,
    stringsAsFactors = FALSE
  )

  # Merge: params periods override observed periods
  # Keep observed periods that don't overlap with params
  obs_only_times <- setdiff(obs_time_levels, params_times)
  obs_wibe_keep <- obs_wibe[obs_wibe$time %in% obs_only_times, ]

  # Combine: params periods + non-overlapping observed periods
  wibe_out <- rbind(params_wibe, obs_wibe_keep)
  wibe_out <- wibe_out[order(wibe_out$time, wibe_out$group), ]
  rownames(wibe_out) <- NULL
  all_time_levels <- sort(unique(wibe_out$time))

  # Validate ref
  if (!(ref %in% all_time_levels)) {
    stop("'ref' = ", ref, " not found in merged time periods. ",
         "Available: ", paste(all_time_levels, collapse = ", "))
  }

  # Compute totals: pi-based for params periods, n-based for observed periods
  totals_list <- list()
  for (t in all_time_levels) {
    if (t %in% params_times) {
      totals_list[[length(totals_list) + 1]] <-
        .compute_totals_pi(wibe_out[wibe_out$time == t, ], t)
    } else {
      totals_list[[length(totals_list) + 1]] <-
        .compute_totals_n(wibe_out[wibe_out$time == t, ], t)
    }
  }
  totals <- do.call(rbind, totals_list)
  rownames(totals) <- NULL

  # Compute counterfactual deltas
  deltas <- .compute_desc_deltas(wibe_out, all_time_levels, ref, ystat, order)

  structure(
    list(
      wibe = wibe_out,
      totals = totals,
      deltas = deltas,
      raw_data = wibe_result$raw_data,
      ystat = ystat,
      order = order,
      ref = ref,
      ref_params = ref_params,
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

  # CV2 requires nonzero grand mean at reference (CV2 = Var/mean^2)
  if (ystat == "CV2" && sum(pi_0 * mu_0) == 0) {
    stop("CV2 decomposition requires a nonzero grand mean at the reference period. ",
         "A zero-mean reference is not supported because CV2 = Var/mean^2 is undefined when mean = 0.")
  }

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
      delta_mu = components$mu$delta_T,
      delta_sigma_W = components$sigma$delta_W,
      delta_sigma_B = components$sigma$delta_B,
      delta_sigma = components$sigma$delta_T,
      delta_pi_W = components$pi$delta_W,
      delta_pi_B = components$pi$delta_B,
      delta_pi = components$pi$delta_T,
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
