#' Synthetic panel of women's weekly earnings
#'
#' A simulated two-wave panel of women's weekly earnings, used to demonstrate
#' the descriptive and causal variance decompositions. Its structure and
#' variable names follow the Current Population Survey (CPS): each woman is
#' observed in two consecutive years, and some women (\code{mother = 1}) have
#' their first child between the two waves. All values are simulated; no real
#' survey records are included.
#'
#' The data-generating process is stylized: motherhood lowers employment and
#' mean earnings and raises earnings dispersion, more so for less-educated
#' women, and these penalties shrink over time; the share of mothers falls,
#' the share with a BA rises, and within-group dispersion grows. Draws are
#' stratified so that results are smooth despite the modest sample size. See
#' \code{create/create_earnings_panel.R} in the package source for the full
#' data-generating process.
#'
#' @format A data frame with 15,200 rows (7,600 women x 2 waves) and 12
#'   variables:
#' \describe{
#'   \item{id}{Person identifier.}
#'   \item{year}{Survey year (1982--2020).}
#'   \item{decade}{Decade of the first wave (1980, 1990, 2000, 2010); the same
#'     for both waves of a person.}
#'   \item{post}{Wave indicator: 0 = first wave, 1 = second wave (after the
#'     birth, for mothers).}
#'   \item{mother}{1 if the woman has her first child between the two waves,
#'     0 otherwise.}
#'   \item{earnweek}{Weekly earnings in dollars (0 if not employed).}
#'   \item{earnwt}{Earnings weight.}
#'   \item{educ}{Education: factor with levels \code{"HS or less"},
#'     \code{"Some college"}, \code{"BA+"}.}
#'   \item{ses}{Socioeconomic status: factor with levels \code{"Low"},
#'     \code{"Middle"}, \code{"High"}.}
#'   \item{age}{Age in years.}
#'   \item{race}{Race/ethnicity: factor with levels \code{"White"},
#'     \code{"Black"}, \code{"Hispanic"}, \code{"Other"}.}
#'   \item{married}{1 if married, 0 otherwise.}
#' }
#'
#' @docType data
#' @usage data(earnings_panel)
#' @keywords datasets
#'
#' @examples
#' data(earnings_panel)
#' head(earnings_panel)
"earnings_panel"

#' Income data
#'
#' A simulated income dataset for demonstrating the variance decomposition
#' functions.
#'
#' @docType data
#' @usage data(incdat)
#' @keywords datasets
#'
#' @examples
#' data(incdat)
"incdat"
