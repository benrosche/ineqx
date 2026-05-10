# ============================================================================ #
# Variance decomposition formulas
# Pure base-R functions, no external dependencies
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Pi-based variants (for causal decomposition)
# pi = normalized group shares (sum to 1)
# ---------------------------------------------------------------------------- #

#' Within-group variance (pi-based)
#' @param pi Numeric vector of group shares (must sum to 1)
#' @param sigma Numeric vector of group standard deviations
#' @return Scalar within-group variance: sum(pi * sigma^2)
#' @keywords internal
VarW_pi <- function(pi, sigma) {
  sum(pi * sigma^2)
}

#' Between-group variance (pi-based)
#' @param pi Numeric vector of group shares (must sum to 1)
#' @param mu Numeric vector of group means
#' @return Scalar between-group variance: sum(pi * (mu - grand_mean)^2)
#' @keywords internal
VarB_pi <- function(pi, mu) {
  gmu <- sum(pi * mu)
  sum(pi * (mu - gmu)^2)
}

#' Total variance (pi-based)
#' @param pi Numeric vector of group shares
#' @param mu Numeric vector of group means
#' @param sigma Numeric vector of group standard deviations
#' @return Scalar total variance: W + B
#' @keywords internal
VarT_pi <- function(pi, mu, sigma) {
  VarW_pi(pi, sigma) + VarB_pi(pi, mu)
}

#' Within-group CV^2 (pi-based)
#' @param pi Numeric vector of group shares
#' @param mu Numeric vector of group means
#' @param sigma Numeric vector of group standard deviations
#' @return Scalar within-group CV^2: sum(pi * sigma^2) / (sum(pi * mu))^2
#' @keywords internal
CV2W_pi <- function(pi, mu, sigma) {
  gmu <- sum(pi * mu)
  num <- sum(pi * sigma^2)
  if (gmu == 0) return(if (num == 0) 0 else Inf)
  num / gmu^2
}

#' Between-group CV^2 (pi-based)
#' @param pi Numeric vector of group shares
#' @param mu Numeric vector of group means
#' @return Scalar between-group CV^2: sum(pi * (mu - grand_mean)^2) / grand_mean^2
#' @keywords internal
CV2B_pi <- function(pi, mu) {
  gmu <- sum(pi * mu)
  num <- sum(pi * (mu - gmu)^2)
  if (gmu == 0) return(if (num == 0) 0 else Inf)
  num / gmu^2
}

#' Total CV^2 (pi-based)
#' @param pi Numeric vector of group shares
#' @param mu Numeric vector of group means
#' @param sigma Numeric vector of group standard deviations
#' @return Scalar total CV^2: CV2W + CV2B
#' @keywords internal
CV2T_pi <- function(pi, mu, sigma) {
  CV2W_pi(pi, mu, sigma) + CV2B_pi(pi, mu)
}

# ---------------------------------------------------------------------------- #
# Count-based variants (for descriptive decomposition)
# n = group counts (not necessarily normalized)
# ---------------------------------------------------------------------------- #

#' Within-group variance (count-based)
#' @param n Numeric vector of group counts
#' @param sigma Numeric vector of group standard deviations
#' @return Scalar within-group variance
#' @keywords internal
VarW_n <- function(n, sigma) {
  pi <- n / sum(n)
  sum(pi * sigma^2)
}

#' Between-group variance (count-based)
#' @param n Numeric vector of group counts
#' @param mu Numeric vector of group means
#' @return Scalar between-group variance
#' @keywords internal
VarB_n <- function(n, mu) {
  pi <- n / sum(n)
  gmu <- sum(pi * mu)
  sum(pi * (mu - gmu)^2)
}

#' Total variance (count-based)
#' @param n Numeric vector of group counts
#' @param mu Numeric vector of group means
#' @param sigma Numeric vector of group standard deviations
#' @return Scalar total variance
#' @keywords internal
VarT_n <- function(n, mu, sigma) {
  VarW_n(n, sigma) + VarB_n(n, mu)
}

#' Within-group CV^2 (count-based)
#' @param n Numeric vector of group counts
#' @param mu Numeric vector of group means
#' @param sigma Numeric vector of group standard deviations
#' @return Scalar within-group CV^2
#' @keywords internal
CV2W_n <- function(n, mu, sigma) {
  pi <- n / sum(n)
  CV2W_pi(pi, mu, sigma)
}

#' Between-group CV^2 (count-based)
#' @param n Numeric vector of group counts
#' @param mu Numeric vector of group means
#' @return Scalar between-group CV^2
#' @keywords internal
CV2B_n <- function(n, mu) {
  pi <- n / sum(n)
  CV2B_pi(pi, mu)
}

#' Total CV^2 (count-based)
#' @param n Numeric vector of group counts
#' @param mu Numeric vector of group means
#' @param sigma Numeric vector of group standard deviations
#' @return Scalar total CV^2
#' @keywords internal
CV2T_n <- function(n, mu, sigma) {
  CV2W_n(n, mu, sigma) + CV2B_n(n, mu)
}

# ---------------------------------------------------------------------------- #
# Helper: Compute delta_B^D (treatment effect on between-group component)
# Used by causal decomposition
# ---------------------------------------------------------------------------- #

#' Compute treatment effect on between-group inequality
#' @param pi Group shares among treated
#' @param mu0 Baseline group means mu_j(0)
#' @param beta Treatment effects on mean
#' @param ystat "Var" or "CV2"
#' @return Scalar delta_B^D
#' @keywords internal
compute_delta_B <- function(pi, mu0, beta, ystat) {
  if (ystat == "Var") {
    VarB_pi(pi, mu0 + beta) - VarB_pi(pi, mu0)
  } else {
    CV2B_pi(pi, mu0 + beta) - CV2B_pi(pi, mu0)
  }
}

#' Compute treatment effect on within-group inequality
#' @param pi Group shares among treated
#' @param mu0 Baseline group means (needed for CV2)
#' @param sigma0 Baseline group SDs sigma_j(0)
#' @param lambda Treatment effects on log-SD
#' @param ystat "Var" or "CV2"
#' @return Scalar delta_W^D
#' @keywords internal
compute_delta_W <- function(pi, mu0, sigma0, lambda, ystat) {
  sigma1 <- sigma0 * exp(lambda)
  if (ystat == "Var") {
    # delta_W^D = W(1) - W(0) = sum(pi * sigma1^2) - sum(pi * sigma0^2)
    # = sum(pi * sigma0^2 * (exp(2*lambda) - 1))
    VarW_pi(pi, sigma1) - VarW_pi(pi, sigma0)
  } else {
    # For CV2: total(1) - B(1) - [total(0) - B(0)]
    # = [W_CV2(1) - W_CV2(0)]
    # Need full CV2 with the post-treatment means for treated
    mu1 <- mu0 + 0  # Note: within-group CV2 uses the SAME means
    # Actually for delta_W^D with CV2, we need:
    # CV2W under treatment - CV2W without treatment
    # But B changes too due to beta... Let's compute total - B for each
    # Under d=1: means are mu0+beta, sigmas are sigma1
    # Under d=0: means are mu0, sigmas are sigma0
    # But wait - delta_W is specifically holding beta constant and only
    # looking at sigma changes. Actually no - delta_W^D = W(d=1) - W(d=0)
    # where W(d) = Total(d) - B(d), and both mu and sigma change with d.
    #
    # For the sequential switching in the longitudinal case, we need
    # the full function. Here for cross-sectional:
    # delta_W uses the actual treated means mu0+beta and untreated means mu0
    # Let me just compute it directly:
    stop("CV2 within-group treatment effect requires beta. Use compute_delta_WB instead.")
  }
}

#' Compute treatment effect on both W and B components
#' @param pi Group shares among treated
#' @param mu0 Baseline group means
#' @param sigma0 Baseline group SDs
#' @param beta Treatment effects on mean
#' @param lambda Treatment effects on log-SD
#' @param ystat "Var" or "CV2"
#' @return List with delta_B, delta_W, delta_total
#' @keywords internal
compute_delta_WB <- function(pi, mu0, sigma0, beta, lambda, ystat) {
  sigma1 <- sigma0 * exp(lambda)
  mu1 <- mu0 + beta

  if (ystat == "Var") {
    B1 <- VarB_pi(pi, mu1)
    B0 <- VarB_pi(pi, mu0)
    W1 <- VarW_pi(pi, sigma1)
    W0 <- VarW_pi(pi, sigma0)
  } else {
    B1 <- CV2B_pi(pi, mu1)
    B0 <- CV2B_pi(pi, mu0)
    W1 <- CV2W_pi(pi, mu1, sigma1)
    W0 <- CV2W_pi(pi, mu0, sigma0)
  }

  list(
    tau_B = B1 - B0,
    tau_W = W1 - W0,
    tau_total = (B1 + W1) - (B0 + W0)
  )
}


# ============================================================================ #
# Individual-level inequality statistics (for plot type = "ineq")
# All take (y, w) where w = weights (NULL or numeric vector)
# ============================================================================ #

#' Weighted variance
#' @keywords internal
.ineq_V <- function(y, w = NULL) {
  if (is.null(w)) return(stats::var(y))
  n <- length(y)
  sw <- sum(w)
  mu <- sum(w * y) / sw
  n / (sw * (n - 1)) * sum(w * (y - mu)^2)
}

#' Variance of log
#' @keywords internal
.ineq_VL <- function(y, w = NULL) {
  y_pos <- y[y > 0]
  if (is.null(w)) {
    w_pos <- NULL
  } else {
    w_pos <- w[y > 0]
  }
  .ineq_V(log(y_pos), w_pos)
}

#' Coefficient of variation squared
#' @keywords internal
.ineq_CV2 <- function(y, w = NULL) {
  v <- .ineq_V(y, w)
  if (is.null(w)) {
    mu <- mean(y)
  } else {
    mu <- sum(w * y) / sum(w)
  }
  if (mu == 0) return(if (v == 0) 0 else Inf)
  v / mu^2
}

#' Weighted Gini coefficient
#' @keywords internal
.ineq_Gini <- function(y, w = NULL) {
  n <- length(y)
  if (n <= 1) return(0)
  if (is.null(w)) w <- rep(1, n)

  ord <- order(y)
  y <- y[ord]
  w <- w[ord]

  sw <- sum(w)
  mu <- sum(w * y) / sw
  if (mu == 0) return(0)

  cum_w <- cumsum(w)
  # Weighted Gini using covariance formula
  # G = (2 / (mu * sw^2)) * sum(w_i * y_i * (cum_w_i - w_i/2)) - 1
  rank_w <- cum_w - w / 2
  2 * sum(w * y * rank_w) / (mu * sw^2) - 1
}

#' Weighted Theil index (GE(1))
#' @keywords internal
.ineq_Theil <- function(y, w = NULL) {
  y_pos <- y[y > 0]
  if (is.null(w)) {
    w_pos <- rep(1, length(y_pos))
  } else {
    w_pos <- w[y > 0]
  }
  if (length(y_pos) <= 1) return(0)

  sw <- sum(w_pos)
  mu <- sum(w_pos * y_pos) / sw
  if (mu == 0) return(0)

  s <- y_pos / mu
  sum(w_pos * s * log(s)) / sw
}

#' Registry mapping stat names to functions
#' @keywords internal
.ineq_stat_registry <- list(
  V     = .ineq_V,
  VL    = .ineq_VL,
  CV2   = .ineq_CV2,
  Gini  = .ineq_Gini,
  Theil = .ineq_Theil
)


# ============================================================================ #
# Influence-function SEs for inequality statistics (for delta CIs)
# All take (y, w, val) where val is the already-computed statistic value.
# SE = sqrt(sum(w^2 * psi^2)) / sum(w)  (linearisation estimator)
# ============================================================================ #

#' @keywords internal
.ineq_se_V <- function(y, w = NULL, val) {
  if (is.null(w)) w <- rep(1, length(y))
  sw <- sum(w)
  mu <- sum(w * y) / sw
  psi <- (y - mu)^2 - val
  sqrt(sum(w^2 * psi^2)) / sw
}

#' @keywords internal
.ineq_se_VL <- function(y, w = NULL, val) {
  keep <- y > 0
  y_k  <- log(y[keep])
  w_k  <- if (is.null(w)) rep(1, sum(keep)) else w[keep]
  .ineq_se_V(y_k, w_k, val)
}

#' @keywords internal
.ineq_se_CV2 <- function(y, w = NULL, val) {
  if (is.null(w)) w <- rep(1, length(y))
  sw <- sum(w)
  mu <- sum(w * y) / sw
  V  <- val * mu^2
  psi <- ((y - mu)^2 - V) / mu^2 - 2 * V * (y - mu) / mu^3
  sqrt(sum(w^2 * psi^2)) / sw
}

#' Linearised SE for weighted Gini (Monti 1991 / Deville 1999)
#' @keywords internal
.ineq_se_Gini <- function(y, w = NULL, val) {
  if (is.null(w)) w <- rep(1, length(y))
  ord  <- order(y)
  y    <- y[ord]
  w    <- w[ord]
  sw   <- sum(w)
  mu   <- sum(w * y) / sw
  if (mu == 0) return(0)
  # weighted CDF at each point (mid-point convention)
  cum_w <- cumsum(w)
  F_i   <- (cum_w - w / 2) / sw          # centred weighted CDF
  L_i   <- cumsum(w * y) / (mu * sw)     # weighted Lorenz ordinates
  # influence function: psi_i = (2/mu)*(F_i*y_i - L_i*mu) - val
  psi   <- (2 / mu) * (F_i * y - L_i * mu) - val
  sqrt(sum(w^2 * psi^2)) / sw
}

#' @keywords internal
.ineq_se_Theil <- function(y, w = NULL, val) {
  keep <- y > 0
  y_k  <- y[keep]
  w_k  <- if (is.null(w)) rep(1, sum(keep)) else w[keep]
  sw   <- sum(w_k)
  mu   <- sum(w_k * y_k) / sw
  if (mu == 0) return(0)
  psi  <- (y_k / mu) * log(y_k / mu) - val
  sqrt(sum(w_k^2 * psi^2)) / sw
}

#' Registry mapping stat names to SE functions (for delta CIs)
#' @keywords internal
.ineq_se_registry <- list(
  V     = .ineq_se_V,
  VL    = .ineq_se_VL,
  CV2   = .ineq_se_CV2,
  Gini  = .ineq_se_Gini,
  Theil = .ineq_se_Theil
)
