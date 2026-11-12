# ---------------------------------------------------------------------------- #
# Create `earnings_panel`: a synthetic two-wave panel of women's weekly earnings
# ---------------------------------------------------------------------------- #
#
# The data are fully simulated. The data-generating process is a stylized
# version of a CPS (IPUMS) panel of women observed in two consecutive years,
# some of whom have their first child between the two waves:
#   - many zero earners; motherhood lowers employment
#   - log earnings rise with education, SES, and over time; the college and
#     SES premiums grow over time
#   - within-group dispersion of log earnings rises over time
#   - motherhood lowers mean log earnings and raises their dispersion, more so
#     for less-educated women, and these penalties shrink over time
#   - the share of mothers declines and the share with a BA rises over time
#
# To keep results smooth with a modest sample size, draws are stratified:
#   - group shares (mother, educ, ses, race) and ages are allocated exactly
#     within each year
#   - employment counts are fixed at their expected value within each
#     year x wave x mother x educ cell
#   - log-earnings residuals are evenly spaced normal quantiles within each
#     decade x wave x mother x educ x ses cell, spread evenly across years and
#     winsorized at the 1st/99th percentiles (akin to CPS top-coding)
# Within-person persistence across waves is induced through rank correlation.
#
# Run from the package root:  source("create/create_earnings_panel.R")
# ---------------------------------------------------------------------------- #

make_earnings_panel <- function(seed = 4, n_per_year = 200,
                                years = 1982:2019) {

  set.seed(seed)

  educ_lab <- c("HS or less", "Some college", "BA+")
  ses_lab  <- c("Low", "Middle", "High")
  race_lab <- c("White", "Black", "Hispanic", "Other")

  # Exact allocation of n units to categories (largest remainder), shuffled
  alloc <- function(n, p) {
    k <- floor(n * p)
    r <- n * p - k
    short <- n - sum(k)
    if (short > 0) k[order(-r)[seq_len(short)]] <- k[order(-r)[seq_len(short)]] + 1
    sample(rep(seq_along(p), k))
  }

  # Evenly spaced standard-normal quantiles, winsorized and standardized
  strat_z <- function(m, cap = Inf) {
    z <- qnorm((seq_len(m) - 0.5) / m)
    z <- pmin(pmax(z, -cap), cap)
    if (m > 2) z <- (z - mean(z)) / sd(z)
    z
  }

  # -------------------------------------------------------------------------- #
  # Person-level characteristics (wave 1)
  # -------------------------------------------------------------------------- #

  persons <- do.call(rbind, lapply(years, function(y) {
    t <- (y - min(years)) / (max(years) - min(years))   # 0 ... 1
    n <- n_per_year

    p_mother <- 0.50 - 0.12 * t
    mother <- alloc(n, c(1 - p_mother, p_mother)) - 1L

    p_educ <- c(0.55 - 0.30 * t, 0.30 + 0.05 * t, 0.15 + 0.25 * t)
    educ <- integer(n)
    for (m in 0:1) educ[mother == m] <- alloc(sum(mother == m), p_educ)

    p_ses <- rbind(c(0.36, 0.52, 0.12),   # HS or less
                   c(0.16, 0.55, 0.29),   # Some college
                   c(0.07, 0.43, 0.50))   # BA+
    ses <- integer(n)
    for (e in 1:3) ses[educ == e] <- alloc(sum(educ == e), p_ses[e, ])

    p_race <- c(0.80 - 0.15 * t, 0.10 + 0.02 * t, 0.06 + 0.11 * t, 0.04 + 0.02 * t)
    race <- alloc(n, p_race)

    # Mothers: first birth in the mid/late 20s, later over time and with
    # education; non-mothers spread evenly over ages 20-48
    age <- integer(n)
    for (e in 1:3) {
      ix <- which(mother == 1 & educ == e)
      if (length(ix)) age[ix] <- sample(pmin(pmax(round(
        26.5 + 3 * t + 1.2 * (e - 1) + 4 * strat_z(length(ix))), 18), 45))
    }
    ix <- which(mother == 0)
    age[ix] <- (20:48)[alloc(length(ix), rep(1 / 29, 29))]

    p_married <- ifelse(mother == 1, 0.85 - 0.20 * t, 0.55 - 0.10 * t)
    married <- as.integer(runif(n) < p_married)

    data.frame(year1 = y, t = t, mother, educ, ses, race, age, married)
  }))

  persons$id <- seq_len(nrow(persons))
  N <- nrow(persons)

  # Person-level latent traits that persist across waves
  persons$u_emp  <- rnorm(N)   # employment propensity
  persons$u_earn <- rnorm(N)   # earnings rank

  # -------------------------------------------------------------------------- #
  # Expand to two waves and simulate earnings
  # -------------------------------------------------------------------------- #

  d <- rbind(transform(persons, post = 0L), transform(persons, post = 1L))
  d$year   <- d$year1 + d$post
  d$age    <- d$age + d$post
  d$t      <- d$t + d$post / (max(years) - min(years))
  d$decade <- 10L * (d$year1 %/% 10L)   # decade of wave 1 (same for both waves)

  treated <- d$mother == 1 & d$post == 1
  a <- d$age - 35

  # Motherhood effects, larger for less-educated women, shrinking over time
  shrink   <- 1 - 0.8 * d$t
  beta_emp <- c(-0.90, -0.70, -0.45)[d$educ] * shrink   # logit employment
  beta_mu  <- c(-0.44, -0.26, -0.11)[d$educ] * shrink   # mean log earnings
  beta_sig <- c( 0.55,  0.40,  0.18)[d$educ] * shrink   # log SD of log earnings

  # Employment: exact expected count per year x wave x mother x educ cell;
  # who is employed follows a persistent propensity
  lp_emp <- 0.60 + c(0, 0.35, 0.90)[d$educ] - 0.002 * a^2 +
    0.25 * d$mother + 0.15 * d$t + treated * beta_emp
  d$employed <- 0L
  emp_cells <- split(seq_len(nrow(d)), list(d$year1, d$post, d$mother, d$educ),
                     drop = TRUE)
  for (ix in emp_cells) {
    k <- round(sum(plogis(lp_emp[ix])))
    score <- lp_emp[ix] + 1.5 * d$u_emp[ix] + rlogis(length(ix), scale = 0.5)
    d$employed[ix[order(-score)[seq_len(k)]]] <- 1L
  }

  # Log earnings among the employed
  mu <- 6.00 + c(0, 0.15, 0.55)[d$educ] +
    c(0, 0.08, 0.18)[d$ses] + c(0, 0.05, 0.15)[d$ses] * d$t +
    (0.15 + c(0, 0.05, 0.30)[d$educ]) * d$t +
    0.01 * a + 0.05 * d$mother + 0.01 * d$post +
    treated * beta_mu
  sig <- (0.35 + c(0, 0.03, 0.06)[d$educ] + 0.15 * d$t) * exp(treated * beta_sig)

  # Residuals: evenly spaced quantiles within each decade x wave x mother x
  # educ x ses cell. Within-year ranks are spread evenly over the cell's range
  # so each year gets a balanced set of quantiles; ranks follow a persistent
  # person-level trait (rank correlation 0.85 across waves).
  rho <- 0.85
  rank_score <- ifelse(d$post == 0, d$u_earn,
                       rho * d$u_earn + sqrt(1 - rho^2) * rnorm(nrow(d)))
  d$z <- 0
  emp <- which(d$employed == 1)
  z_cells <- split(emp, d[emp, c("decade", "post", "mother", "educ", "ses")],
                   drop = TRUE)
  for (ix in z_cells) {
    u <- ave(rank_score[ix], d$year1[ix],
             FUN = function(s) (rank(s) - 0.5 + runif(length(s), -0.05, 0.05)) / length(s))
    d$z[ix[order(u)]] <- strat_z(length(ix), cap = qnorm(0.99))
  }

  d$earnweek <- ifelse(d$employed == 1, round(exp(mu + d$z * sig), 2), 0)

  # Survey weights (CPS-like magnitude, mildly variable, constant within person)
  d$earnwt <- round(exp(rnorm(N, log(2500), 0.15)))[d$id]

  # -------------------------------------------------------------------------- #
  # Finalize
  # -------------------------------------------------------------------------- #

  d$educ    <- factor(educ_lab[d$educ], levels = educ_lab)
  d$ses     <- factor(ses_lab[d$ses],   levels = ses_lab)
  d$race    <- factor(race_lab[d$race], levels = race_lab)
  d$age     <- as.integer(d$age)
  d$year    <- as.integer(d$year)
  d$mother  <- as.integer(d$mother)
  d$married <- as.integer(d$married)

  out <- d[order(d$id, d$post),
           c("id", "year", "decade", "post", "mother", "earnweek",
             "earnwt", "educ", "ses", "age", "race", "married")]
  rownames(out) <- NULL
  out
}

if (sys.nframe() == 0L || !exists("no_save")) {
  earnings_panel <- make_earnings_panel()
  save(earnings_panel, file = "data/earnings_panel.rda", compress = "xz")
}
