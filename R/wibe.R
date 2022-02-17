#' @title Descriptive within/between decomposition
#'
#' @description [...]
#'
#' @param y Dependent variable
#' @param group Grouping variable to decompose variance into within- and between-group components
#' @param time Time variable to analyze change over time
#' @param weights Probability weights
#' @param dat Dataframe
#' @param smoothDat Logical. Should data be smoothed?
#' @param ref Number or FALSE. Should values be reported in reference to a specific time?
#' @param long Logical. Should output be in long format?
#'
#' @return List of length 2. Element 1 returns the decomposition by group and time. Elements 2 returns the decomposition by time.
#'
#' @examples data(incdat)
#' wibe1 <- wibe(y="inc", group="SES", time="year", dat=dat1)
#'
#' @export wibe
#'
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#'
#' @details
#' ...

wibe <- function(y=NULL, group=NULL, time=NULL, weights=NULL, dat, smoothDat=F, ref=F, long=F) {

  # dat <- dat.dWB1 %>% dplyr::filter(x==1); time <- "year"; group <- "i.group"; y <- "c.inc"; weights=NULL; ref <- 1; smoothDat <- F

  # ============================================================================================== #
  # Dissect input ----
  # ============================================================================================== #

  if(isTRUE(ref)) stop("ref must be numeric")

  # Y
  y <- dissectVar(y, cicheck="c")[[1]]

  # Group
  if(is.null(rlang::enexpr(group))) {
    dat <- dat %>% dplyr::mutate(group=1)
    group <- substitute(group)
  }
  group <- dissectVar(group, cicheck="i")[[1]]

  # Time
  if(is.null(rlang::enexpr(time))) {
    dat <- dat %>% dplyr::mutate(time=1)
    time <- substitute(time)
  }
  time <- dissectVar(time, cicheck="i")[[1]]

  # Check whether variables are in the dataset
  if(!as.character(y) %in% names(dat)) stop(paste0(y, " not in dataset"))
  if(!as.character(group) %in% names(dat)) stop(paste0(group, " not in dataset"))
  if(!as.character(time) %in% names(dat)) stop(paste0(time, " not in dataset"))
  if(!is.null(weights)) if(!as.character(weights) %in% names(dat)) stop(paste0(weights, " not in dataset"))

  # ============================================================================================== #
  # Rename ----
  # ============================================================================================== #

  # Rename variables and filter NAs to avoid creating another grouping level
  dat <-
    dat %>%
    dplyr::select({{ group }}, {{ time }}, {{ y }}, {{ weights }}) %>%
    dplyr::rename(group:={{ group }}, time:={{ time }}, y:={{ y }}, w:={{ weights }}) %>%
    drop_na()

  # Weights
  if(is.null(weights)) dat <- dat %>% dplyr::mutate(w=1)

  # Number of levels of group and time var
  group_levels <- dat %>% .$group %>% unique() %>% sort()
  time_levels <- dat %>% .$time %>% unique() %>% sort()

  # ============================================================================================== #
  # Create dat.between ----
  # ============================================================================================== #

  dat.between <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    left_join(
      dat %>%
        group_by(time, group) %>%
        dplyr::summarise(
          n=n(),
          sw=sum(w),
          mu=1/sw * sum(w * y, na.rm = T),
          sigma2=n/(sw*(n-1)) * sum(w*(y-mu)^2),
          sigma=sqrt(sigma2)
        ) %>%
        ungroup() %>%
        dplyr::mutate(
          n=sw/sum(sw)*sum(n), # weighted n
          CV2=sigma2/(mu^2) # squared CV
        ) %>%
        dplyr::select(time, group, n, mu, sigma, sigma2, CV2),
      by=c("time", "group"))

  # Weights follow Stata's sum command:
  # https://www.stata.com/support/faqs/statistics/weights-and-summary-statistics/

  # Replace n, mu, sigma2, and CV2 with smoothed versions?
  if(smoothDat==T) {

    # Creates a smoothed version of "var"
    loess.pred <- function(.dat, var) {
      # Function arguments
      # var: string
      # Returns predictions of the same size as .dat
      return(predict(loess(as.formula(paste(var, "~ time + group + time*group")), control=loess.control(surface="direct"), .dat), .dat))
    }

    dat.between <-
      dat.between %>%
      mutate(n=loess.pred(., "n")) %>% # I replace the actual values with a smoothed version from loess
      mutate(mu=loess.pred(., "mu")) %>%
      mutate(sigma=loess.pred(., "sigma")) %>%
      mutate(sigma2=loess.pred(., "sigma2")) %>%
      mutate(CV2=loess.pred(., "CV2"))
  }

  # ============================================================================================== #
  # Create dat.total ----
  # ============================================================================================== #

  dat.total <-
    dat.between %>%
    group_by(time) %>%
    dplyr::summarise(N=sum(n, na.rm = T),
                     gmu=sum(n*mu, na.rm = T)/N,
                     VarW=sum(n*sigma2, na.rm = T)/sum(n, na.rm = T),
                     VarB=sum(n*(mu-gmu)^2, na.rm = T)/sum(n, na.rm = T),
                     CV2W=sum((n/N)*(mu/gmu)^2*CV2, na.rm = T),
                     CV2B=(sum((n/N)*(((mu/gmu)^2)-1), na.rm = T))) %>%
    ungroup() %>%
    mutate(VarWBRatio=VarW/VarB, CV2WBRatio=CV2W/CV2B) %>%
    # mutate(VarT=VarW+VarB, CV2T=CV2W+CV2B) %>%
    inner_join(
      dat %>%
        group_by(time) %>%
        dplyr::summarise(
          VarT = var(y, na.rm = T),
          CV2T = var(y, na.rm = T)/(mean(y, na.rm = T)^2)
        ) %>%
        ungroup(), by="time")

  # ============================================================================================== #
  # Prepare output ----
  # ============================================================================================== #

  if(!isFALSE(ref)) {

    # Relative to year == ref

    dat.between <-
      dat.between %>%
      group_by(group) %>%
      mutate(
        mu1= mu[time==ref],
        mu = (mu - mu1)/mu1*100+100,
        sigma1 = sigma[time==ref],
        sigma  = (sigma - sigma1)/sigma1*100+100,
        sigma21 = sigma2[time==ref],
        sigma2  = (sigma2 - sigma21)/sigma21*100+100,
        CV21 = CV2[time==ref],
        CV2  = (CV2 - CV21)/CV21*100+100
      ) %>%
      ungroup() %>%
      dplyr::select(-mu1, -sigma1, -sigma21, -CV21)

    dat.total <-
      dat.total %>%
      mutate(
        VarB=(VarB-VarB[time==ref])/VarB[time==ref]*100+100,
        VarW=(VarW-VarW[time==ref])/VarW[time==ref]*100+100,
        VarT=(VarT-VarT[time==ref])/VarT[time==ref]*100+100,
        VarWBRatio=(VarWBRatio-VarWBRatio[time==ref])/VarWBRatio[time==ref]*100+100,
        CV2B=(CV2B-CV2B[time==ref])/CV2B[time==ref]*100+100,
        CV2W=(CV2W-CV2W[time==ref])/CV2W[time==ref]*100+100,
        CV2T=(CV2T-CV2T[time==ref])/CV2T[time==ref]*100+100,
        CV2WBRatio=(CV2WBRatio-CV2WBRatio[time==ref])/CV2WBRatio[time==ref]*100+100
      )

  }

  if(long==T) dat.between <- dat.between %>% pivot_longer(cols=c(-time, -group), names_to = "variable", values_to = "value")
  if(long==T) dat.total <- dat.total %>% pivot_longer(cols=c(-time), names_to = "variable", values_to = "value")

  # Rename variables back to their original names
  dat.between <-
    dat.between %>%
    dplyr::rename(
      !!enquo(time) := time,
      !!enquo(group) := group
    )

  dat.total <-
    dat.total %>%
    dplyr::rename(
      !!enquo(time) := time
    )

  return(list("By group and time"=dat.between, "By time"=dat.total))

}
