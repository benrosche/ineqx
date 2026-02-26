# ============================================================================ #
# extract_params: Thin wrapper for backward compatibility
# Delegates to ineqx_params(model = ...)
# ============================================================================ #

#' @keywords internal
extract_params <- function(model, treat, group, time = NULL, post = NULL,
                           data, ref = NULL, ystat = "Var", vcov = TRUE) {
  ineqx_params(
    data = data, model = model,
    treat = treat, group = group,
    time = time, post = post,
    ystat = ystat, ref = ref, vcov = vcov
  )
}
