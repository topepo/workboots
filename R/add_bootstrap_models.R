#' Fit from a workflow using many bootstrap resamples.
#'
#' Generate a prediction interval from arbitrary model types using bootstrap
#' resampling. `add_bootstrap_models()` fits a model to each resample.
#'
#'
#' @param interval One of `prediction`, `confidence`. Specifies the interval
#' type to be generated.
#' @param remove_splits A logical to remove the `splits` column.
#' @inheritParams vi_boots
#'
#' @export
add_bootstrap_models <- function(resamples,
                                 workflow,
                                 interval = c("prediction", "confidence"),
                                 verbose = FALSE,
                                 remove_splits = TRUE,
                                 ...) {

  if (!inherits(resamples, "bootstraps")) {
    cli::cli_abort("{.arg resamples} should be generated from {.fn rsample::bootstraps}.")
  }
  # TODO check mode

  apparently <- resamples$id == "Apparent"
  if (any(apparently)) {
    resamples <- resamples[!apparently, ]
  }

  # convert interval type
  interval <- rlang::arg_match(interval)

  # # check arguments
  # workboots:::assert_workflow(workflow)
  # workboots:::assert_n(n)
  # workboots:::assert_pred_data(workflow, training_data, "training")

  # check apparent

  # warn if low n
  if (nrow(resamples) < 2000) {

    rlang::warn(
      paste0("At least 2000 resamples recommended for stable results.")
    )

  }

  req_pkgs <- c("workboots", "parsnip", "workflows", "rsample",
                generics::required_pkgs(workflow$.models[[1]]$fit))
  req_pkgs <- unique(req_pkgs)
  rlang::check_installed(req_pkgs)

  model_res <-
    future.apply::future_lapply(
      resamples$splits,
      fit_single_model,
      workflow = workflow,
      interval = interval,
      verbose = verbose,
      future.packages = req_pkgs,
      future.seed = TRUE,
      future.stdout = TRUE
    )

  # Check for failed models
  fit_failed <- purrr::map_lgl(model_res, is.null)
  if (any(fit_failed)) {
    resamples <- resamples[!fit_failed, ]
    model_res <- model_res[!fit_failed]
  }
  resamples$.models <- model_res

  if (remove_splits) {
    resamples$splits <- NULL
  }

  # TODO make control function?

  class(resamples) <- c("bootstrapped_models", class(resamples))
  resamples

}

fit_single_model <- function(split,
                             workflow,
                             interval,
                             verbose) {

  # get training data from bootstrap resample split
  boot_train <- rsample::analysis(split)

  # get oob sample
  boot_oob <- rsample::assessment(split)

  # fit workflow to training data
  model <- try(generics::fit(workflow, boot_train), silent = TRUE)
  if (inherits(model, "try-error")) {
    return(NULL)
  }

  # get predicted var name
  # TODO update from main
  pred_name <- names(model$pre$mold$blueprint$ptypes$outcomes)

  # TODO make this a parameter (default = "prediction") for a control object?

  # apply prediction interval using bootstrap 632+ estimate
  # if not, just returns absolute prediction (when summarised, this generates a confidence interval)
  if (interval == "prediction") {

    # get training residuals
    preds_train <- dplyr::pull(stats::predict(model, boot_train), .pred)
    actuals_train <- dplyr::pull(boot_train, rlang::sym(pred_name))
    resids_train <- actuals_train - preds_train
    resids_train <- resids_train - mean(resids_train)

    # get oob residuals
    preds_oob <- dplyr::pull(stats::predict(model, boot_oob), .pred)
    actuals_oob <- dplyr::pull(boot_oob, rlang::sym(pred_name))
    resids_oob <- actuals_oob - preds_oob
    resids_oob <- resids_oob - mean(resids_oob, na.rm = TRUE)

    # TODO lose a dependency by writing a simple RMSE function (that takes care of NA values too)
    # calculate no-information error rate (rmse_ni) with RMSE as loss function
    combos <- tidyr::crossing(actuals_train, preds_train)
    rmse_ni <- Metrics::rmse(combos$actuals_train, combos$preds_train)

    # calculate overfit rate
    rmse_oob <- Metrics::rmse(actuals_oob, preds_oob)
    rmse_train <- Metrics::rmse(actuals_train, preds_train)
    overfit <- (rmse_oob - rmse_train)/(rmse_ni - rmse_train)

    # calculate weight (if overfit = 0, weight = .632 & residual used will just be .632)
    # uses the actual proportion of distinct training/oob samples, rather than average of 0.632/0.368
    prop_368 <- nrow(boot_oob) / nrow(boot_train)
    prop_632 <- 1 - prop_368
    weight <- prop_632 / (1 - (prop_368 * overfit))

    # determine residual std.dev based on weight
    sd_oob <- stats::sd(resids_oob)
    sd_train <- stats::sd(resids_train)
    sd_resid <- weight * sd_oob + (1 - weight) * sd_train

  } else {
    sd_resid <- NA_real_
  }

  # use progressr?
  # print progress when verbose is set to TRUE
  # verbose_print(verbose, index, times)

  # bundle?
  res <- list(fit = model, sd_resid = sd_resid, interval = interval)
  class(res) <- "boot_model"
  res
}
