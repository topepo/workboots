#' Calculate bootstrap prediction intervals for new samples
#' @export
#'
predict.bootstrapped_models <- function(x, new_data, ..., interval_width) {

  pred_res <-
    future.apply::future_lapply(
      x$.models,
      predict_new,
      new_data = new_data,
      future.packages = c("tidymodels", "earth"),
      future.seed = TRUE,
      future.stdout = TRUE
    )
  pred_res <- do.call("cbind", pred_res)
  low_quant <- (1 - interval_width) / 2
  high_quant <- 1 - low_quant
  val_lower <- purrr::map_dbl(pred_res, ~ unname(quantile(.x, probs = low_quant, na.rm = TRUE)))
  val_median <- purrr::map_dbl(pred_res, ~ median(.x, na.rm = TRUE))
  val_upper <- purrr::map_dbl(pred_res, ~ unname(quantile(.x, probs = high_quant, na.rm = TRUE)))
  tibble::tibble(.lower = val_lower, .pred = val_median, .upper = val_upper)
}


predict_new <- function(x, new_data) {
  n <- nrow(new_data)
  res <- predict(x$fit, new_data)
  if (x$interval == "prediction") {
    res$.pred <- res$.pred + stats::rnorm(n, 0, x$sd_resid)
  }
  res
}
