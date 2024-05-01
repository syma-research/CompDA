test_oneX <- function(
    x, y, z, epsilon,
    family_y = "binomial",
    method_y = "glmnet", cv_folds = 5,
    m = 1e4, debug_file = NULL) {

  n <- length(x)

  # fit y vs. z model
  y_res <- get_y_res(x = z, y = y,
                     family = family_y,
                     method = method_y, nfolds = cv_folds)
  y_res <- y_res / sd(y_res)

  # fit x vs z model
  delta <- x > 0
  if(sum(delta) <= 2) {
    result <- c(0, 1)
    if(!is.null(debug_file)) {
      result_debug <- list(y_res = y_res,
                           x_res = rep(0, n),
                           xtilde_res = matrix(n, nrow = n, ncol = m))
      save(result_debug, file = debug_file)
    }
    return(result)
  }

  fitted_delta <- glmnet_wrapper(
    x = z, y = delta, family = "binomial",
    nfolds = cv_folds, alpha = 1,standardize = TRUE,thresh = 1e-04)
  fitted_mean <- glmnet_wrapper(
    x = z[delta, , drop = FALSE], y = log(x[delta] + epsilon),
    family = "gaussian", nfolds = cv_folds, alpha = 1,
    standardize = TRUE, thresh = 1e-04)
  predicted_mean <-
    predict_wrapper(fitted = fitted_mean, newx = z,
                    s = "lambda.min", y = log(x[delta] + epsilon))
  x_res_nonzero <-
    log(x[delta] + epsilon) - predicted_mean[delta]
  fitted_var <- glmnet_wrapper(
    x = z[delta, , drop = FALSE], y = x_res_nonzero^2,
    family = "gaussian",
    nfolds = cv_folds, alpha = 1, standardize = TRUE, thresh = 1e-04)
  x_res_sd <- sqrt(predict_wrapper(fitted = fitted_var, newx = z,
                                   s = "lambda.1se", y = x_res_nonzero^2,
                                   nn = TRUE))
  x_res_nonzero_stand <- x_res_nonzero / x_res_sd[delta]

  predicted_delta <- predict_wrapper(fitted = fitted_delta, newx = z,
                                     s = "lambda.min", y = delta)
  mean_x <- (1 - predicted_delta) * log(epsilon) +
    predicted_delta * predicted_mean
  var_x <- (1 - predicted_delta) * (log(epsilon) - mean_x)^2 +
    predicted_delta * ((predicted_mean - mean_x)^2 + x_res_sd^2)

  # generate randomized samples
  deltatilde <- matrix(rbinom(n = n * m, size = 1, prob = predicted_delta),
                       nrow = n)
  i_iter <- 0
  while(TRUE) {
    if(i_iter > 100)
      stop("Cannot generate deltatilde that is not all zero!")
    ind_allzero <- apply(deltatilde == 0, 2, all)
    if(all(!ind_allzero))
      break
    deltatilde[, ind_allzero] <- rbinom(n = n * sum(ind_allzero),
                                        size = 1,
                                        prob = predicted_delta)
    i_iter <- i_iter + 1
  }
  xtilde_nonzero <- seq(1, m) %>%
    vapply(
      function(i_m)
        sample(x_res_nonzero_stand, replace = TRUE, size = n),
      rep(0.0, n)) * x_res_sd +
    predicted_mean
  xtilde <- deltatilde * xtilde_nonzero  +
    (1 - deltatilde) * log(epsilon)

  x_res <- log(x + epsilon) - mean(log(x + epsilon)) - mean_x + mean(mean_x)
  xtilde_res <- t(t(xtilde) - apply(xtilde, 2, mean)) - mean_x + mean(mean_x)
  x_res <- x_res / sd(x_res)
  xtilde_res <- t(t(xtilde_res) / apply(xtilde_res, 2, sd))

  beta_original <- mean(y_res * x_res)
  beta_perm <- apply(y_res * xtilde_res, 2, mean)
  p <- (sum(abs(beta_original) <= abs(beta_perm)) + 1) / (m + 1)

  result <- c(beta_original, p)
  if(!is.null(debug_file)) {
    result_debug <- list(y_res = y_res,
                         x_res = x_res,
                         xtilde_res = xtilde_res)
    save(result_debug, file = debug_file)
  }

  return(result)
}

glmnet_wrapper <- function(
    x, y, family,
    alpha = 1, nfolds = 10, standardize = TRUE, thresh = 1e-4) {
  glmnet_options <- glmnet::glmnet.control()
  glmnet::glmnet.control(epsnr = 1.01e-4)
  err_fit <- tryCatch(
    fitted <- glmnet::cv.glmnet(x = x,
                                y = y,
                                family = family,
                                nfolds = nfolds,
                                alpha = alpha,
                                standardize = standardize,
                                thresh = thresh),
    error = function(e) e
  )
  do.call(glmnet::glmnet.control, glmnet_options)
  return(err_fit)
}

predict_wrapper <- function(
    fitted, newx, s, y, nn = FALSE, thresh = 1e-4
) {
  if("error" %in% class(fitted)) {
    return(rep(mean(y), length = nrow(newx)))
  } else {
    predicted <- glmnet:::predict.cv.glmnet(
      fitted,
      newx = newx,
      s = s)[, 1]
    if(is.character(family(fitted)))
      if(family(fitted) == "binomial")
        predicted <- exp(predicted) / (1 + exp(predicted))
    if(nn)
      predicted[predicted < thresh] <- thresh
    return(predicted)
  }
}

get_y_res <- function(x, y,
                      family = "binomial",
                      method = "glmnet", nfolds = cv_folds) {
  if(method == "glmnet") {
    fitted_y <-  glmnet_wrapper(
      x = x, y = y, family = family,
      nfolds = nfolds, alpha = 1,
      standardize = TRUE, thresh = 1e-04)
    yhat <- predict_wrapper(fitted = fitted_y, newx = x,
                            s = "lambda.min", y = y)
  } else {
    if(method == "gbt") {
      tag <- "xgbTree"
    } else if(method == "rf") {
      tag <- "rf"
    } else if (method == "svm") {
      tag <- "svmRadial"
    }
    x_qc <- x[, apply(x, 2, function(i_x) length(unique(i_x))) > 1,
              drop = FALSE]
    fit_control <-
      caret::trainControl(
        method = "cv",
        number = nfolds,
        classProbs =  TRUE)

    if(family == "binomial") {
      y_fct <- factor(paste0("X", y))
    }

    vb <- capture.output(
      fitted_y <-
        caret::train(x = x_qc,
                     y = {
                       if(family == "binomial")
                         y_fct
                       else
                         y},
                     method = tag,
                     metric = {
                       if(family == "binomial")
                         "Accuracy"
                       else
                         "RMSE"},
                     trControl = fit_control))
    if(family == "binomial")
      yhat <- caret::predict.train(fitted_y, newdata = x_qc, type = "prob")[, 2]
    else
      yhat <- caret::predict.train(fitted_y, newdata = x_qc, type = "raw")
  }
  y_res <- y - yhat
  return(y_res)
}
