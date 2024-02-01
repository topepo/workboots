# Util function for checking workflow passed to predict/vi boots
assert_workflow <- function(workflow) {

  # assert that object is actually a workflow
  assertthat::assert_that(
    assertthat::are_equal(class(workflow), "workflow"),
    msg = "argument `workflow` must be of class \"workflow\"."
  )

  # Check for a regression model
  # wflow_mode <-
  #   workflow %>%
  #   hardhat:::extract_fit_parsnip() %>%
  #   purrr::pluck("spec") %>%
  #   purrr::pluck("mode")
  #
  # if (wflow_mode != "regression") {
  #   cli::cli_abort("workboots only works with regression models. The mode \\
  #                   of this model is {.val {wflow_mode}}.")
  # }

  # check that it is fitted

}

# Util function for checking that workflow argument is not set to tune()
assert_tune <- function(workflow, index) {

  # get model arg (need to know if it's NULL or not)
  model_arg <- as.character(rlang::eval_tidy(workflow$fit$actions$model$spec$args[[index]]))

  # don't check if NULL
  if (length(model_arg) != 0) {

    assertthat::assert_that(
      model_arg[1] != "tune",
      msg = paste0("all tuning parameters must be final before passing workflow to `predict_boots()`.")
    )

  }

}

# Util function for checking n param
assert_n <- function(n) {

  # >= 1
  assertthat::assert_that(
    n >= 1,
    msg = "argument `n` must be >= 1."
  )

  # don't pass double
  assertthat::assert_that(
    n == as.integer(n),
    msg = "argmuent `n` must be an integer."
  )

}

# Util function for checking training_data and new_data params
assert_pred_data <- function(workflow, data, type) {

  # get colnames from workflow
  var_info <- workflow$pre$actions$recipe$recipe$var_info

  if (type == "training") {

    # check that colnames include all predictors and outcomes
    var_info <- dplyr::filter(var_info, role %in% c("predictor", "outcome"))

  } else { # type == "new"

    # check that colnames include all predictors
    var_info <- dplyr::filter(var_info, role == "predictor")

  }

  # get colnames for comparison
  cols_wf <- var_info$variable
  cols_dat <- colnames(data)

  # check that all cols in wf appear in data
  # message displays any names that appear in cols_wf that don't appear in cols_dat
  assertthat::assert_that(
    sum(cols_wf %in% cols_dat) == length(cols_wf),
    msg = paste0("missing cols in ", type, "_data:\n",
                 paste(cols_wf[which(!cols_wf %in% cols_dat)],
                       collapse = ", "))
  )

}

# Util function for checking .data passed to summary functions
assert_pred_summary <- function(data) {

  # check for col .preds
  assertthat::assert_that(
    ".preds" %in% names(data),
    msg = "col `.preds` missing."

  )

  # check that .preds is list
  assertthat::assert_that(
    class(data$.preds) == "list",
    msg = "col `.preds` must be a list-col."
  )

  # check that model.pred is numeric
  assertthat::assert_that(
    class(data$.preds[[1]]$model.pred) == "numeric",
    msg = "col `model.pred` must be numeric."
  )

}

# Util function for checking .data passed to importance summary function
assert_importance_summary <- function(data) {

  # check for col .importances
  assertthat::assert_that(
    ".importances" %in% names(data),
    msg = "col `.importances` missing."
  )

  # check that .importances is a list
  assertthat::assert_that(
    class(data$.importances) == "list",
    msg = "col `.importances` must be a list-col."
  )

  # check that model.importance is numeric
  assertthat::assert_that(
    class(data$.importances[[1]]$model.importance) == "numeric",
    msg = "col `model.importance` must be numeric."
  )

}

# Util function for checkinginterval
assert_interval <- function(interval_width) {

  # numeric
  assertthat::assert_that(
    is.numeric(interval_width),
    msg = "argument `interval_width` must be numeric."
  )

  # must be between [0, 1]
  assertthat::assert_that(
    interval_width >= 0 && interval_width <= 1,
    msg = "argument `interval_width` must be between [0, 1]."
  )

}

check_installs_fitting <- function(x) {
  req_pkgs <- c("workboots", "parsnip", "workflows", "rsample",
                generics::required_pkgs(x))
  req_pkgs <- unique(req_pkgs)
  rlang::check_installed(req_pkgs)
  req_pkgs
}


check_installs_predicting <- function(x) {
  req_pkgs <- c("workboots", "parsnip", "workflows",
                generics::required_pkgs(x))
  req_pkgs <- unique(req_pkgs)
  rlang::check_installed(req_pkgs)
  req_pkgs
}

check_resamples <- function(x) {
  if (!inherits(x, "bootstraps")) {
    cli::cli_abort("{.arg resamples} should be generated from {.fn rsample::bootstraps}.")
  }

  # Exclude apparent samples, if any
  apparently <- x$id == "Apparent"
  if (any(x)) {
    x <- x[!apparently, ]
  }

  # warn if low n
  if (nrow(x) < 2000) {
    cli::cli_warn("At least 2000 resamples recommended for stable results.")
  }

  x
}
