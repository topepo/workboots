#' Calculate bootstrap prediction intervals for new samples.
#'
#' @param object An object of class `"bootstrapped_models"` created by
#' [add_bootstrap_models()].
#' @param new_data A data frame of predictor values.
#' @param interval_width The coverage rate of the intervals (e.g. 0.95 for a
#' 95% interval).
#' @param ... Not used.
#' @export
#'
predict.bootstrapped_models <- function(object, new_data, ..., interval_width) {

  req_pkgs <- c("workboots", "parsnip", "workflows",
                generics::required_pkgs(object$.models[[1]]$fit))
  req_pkgs <- unique(req_pkgs)
  rlang::check_installed(req_pkgs)

  pred_res <-
    future.apply::future_lapply(
      object$.models,
      predict_new,
      new_data = new_data,
      future.packages = req_pkgs,
      future.seed = TRUE,
      future.stdout = TRUE
    )

  pred_res <- do.call("cbind", pred_res)
  pred_res <- as.matrix(pred_res)
  low_quant <- (1 - interval_width) / 2
  high_quant <- 1 - low_quant
  val_lower <- apply(pred_res, 1, function(x) unname(quantile(x, probs = low_quant, na.rm = TRUE)))
  val_median <- apply(pred_res, 1, function(x) stats::median(x, na.rm = TRUE))
  val_upper <- apply(pred_res, 1, function(x) unname(quantile(x, probs = high_quant, na.rm = TRUE)))
  tibble::tibble(.lower = val_lower, .pred = val_median, .upper = val_upper)
}


predict_new <- function(x, new_data) {
  # TODO unbundle?
  n <- nrow(new_data)
  res <- predict(x$fit, new_data)
  if (x$interval == "prediction") {
    res$.pred <- res$.pred + stats::rnorm(n, 0, x$sd_resid)
  }
  res
}
