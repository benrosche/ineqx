#' @title Descriptive within/between decomposition
#'
#' @description [...]
#'
#' @param dat Dataframe
#' @param y Dependent variable
#' @param groupvar Grouping variable to decompose variance into within- and between-group components
#' @param timevar Time variable to analyze change over time
#' @param smoothDat Logical. Should data be smoothed?
#' @param rel Number or FALSE. Should values be reported relative to specified time?
#' @param long Logical. Should output be in long format?
#'
#' @return List of length 2. Element 1 returns the decomposition by group and time. Elements 2 returns the decomposition by time.
#'
#' @examples data(incdat)
#' wibe1 <- wibe(dat, y=inc, groupvar=SES, timevar=year)
#'
#' @export wibe
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>

wibe <- function(y=NULL, groupvar=NULL, timevar=NULL, dat, smoothDat=F, rel=F, long=F) {

  # dat <- incdat; timevar <- "i.year"; groupvar <- "i.group"; y <- "inc"; rel <- F; smoothDat <- F

  # ============================================================================================== #
  # Dissect input
  # ============================================================================================== #

  if(isTRUE(rel)) stop("rel must be numeric")

  # Y
  y <- dissectVar(y, cicheck="c")[[1]]

  # Group
  if(is.null(rlang::enexpr(groupvar))) {
    dat <- dat %>% dplyr::mutate(group=1)
    groupvar <- substitute(group)
  }
  groupvar <- dissectVar(groupvar, cicheck="i")[[1]]

  # Time
  if(is.null(rlang::enexpr(timevar))) {
    dat <- dat %>% dplyr::mutate(time=1)
    timevar <- substitute(time)
  }
  timevar <- dissectVar(timevar, cicheck="i")[[1]]

  # Check whether variables are in the dataset
  if(!as.character(y) %in% names(dat)) stop(paste0(y, " not in dataset"))
  if(!as.character(groupvar) %in% names(dat)) stop(paste0(groupvar, " not in dataset"))
  if(!as.character(timevar) %in% names(dat)) stop(paste0(timevar, " not in dataset"))

  # ============================================================================================== #
  # Rename
  # ============================================================================================== #

  # Rename variables and filter NAs to avoid creating another grouping level
  dat <- dat %>% dplyr::rename(group:={{ groupvar }}, time:={{ timevar }}, y:={{ y }}) %>% filter(!is.na(group))

  # Number of levels of group and time var
  group_levels <- dat %>% .$group %>% unique() %>% sort()
  time_levels <- dat %>% .$time %>% unique() %>% sort()

  # ============================================================================================== #
  # Create dat.between
  # ============================================================================================== #

  dat.between <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    left_join(
      dat %>%
        group_by(time, group) %>%
        dplyr::summarise(
          mu=mean(y, na.rm = T),
          sigma=sd(y, na.rm = T),
          sigma2=var(y, na.rm = T),
          n=n(),
        ) %>%
        ungroup() %>%
        mutate(CV2=sigma2/(mu^2)) %>% # Squared CV
        dplyr::select(time, group, n, mu, sigma, sigma2, CV2),
      by=c("time", "group"))

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
  # Create dat.total
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
  # Prepare output
  # ============================================================================================== #

  if(!isFALSE(rel)) {

    # Relative to year == rel

    dat.between <-
      dat.between %>%
      group_by(group) %>%
      mutate(
        mu1= mu[time==rel],
        mu = (mu - mu1)/mu1*100+100,
        sigma1 = sigma[time==rel],
        sigma  = (sigma - sigma1)/sigma1*100+100,
        sigma21 = sigma2[time==rel],
        sigma2  = (sigma2 - sigma21)/sigma21*100+100,
        CV21 = CV2[time==rel],
        CV2  = (CV2 - CV21)/CV21*100+100
      ) %>%
      ungroup() %>%
      dplyr::select(-mu1, -sigma1, -sigma21, -CV21)

    dat.total <-
      dat.total %>%
      mutate(
        VarB=(VarB-VarB[time==rel])/VarB[time==rel]*100+100,
        VarW=(VarW-VarW[time==rel])/VarW[time==rel]*100+100,
        VarT=(VarT-VarT[time==rel])/VarT[time==rel]*100+100,
        VarWBRatio=(VarWBRatio-VarWBRatio[time==rel])/VarWBRatio[time==rel]*100+100,
        CV2B=(CV2B-CV2B[time==rel])/CV2B[time==rel]*100+100,
        CV2W=(CV2W-CV2W[time==rel])/CV2W[time==rel]*100+100,
        CV2T=(CV2T-CV2T[time==rel])/CV2T[time==rel]*100+100,
        CV2WBRatio=(CV2WBRatio-CV2WBRatio[time==rel])/CV2WBRatio[time==rel]*100+100
      )

  }

  if(long==T) dat.between <- dat.between %>% pivot_longer(cols=c(-time, -group), names_to = "variable", values_to = "value")
  if(long==T) dat.total <- dat.total %>% pivot_longer(cols=c(-time), names_to = "variable", values_to = "value")

  # Rename variables back to their original names
  dat.between <-
    dat.between %>%
    dplyr::rename(
      !!enquo(timevar) := time,
      !!enquo(groupvar) := group
    )

  dat.total <-
    dat.total %>%
    dplyr::rename(
      !!enquo(timevar) := time
    )

  return(list("By group and time"=dat.between, "By time"=dat.total))

}
