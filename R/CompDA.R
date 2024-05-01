#' Compositional Differential Abundance Analysis for Microbiome Data
#'
#' @param x A matrix of microbiome relative abundances. Rows are samples and columns are microbes.
#' @param y The outcome of interests. Currently only continuous and binary y are supported.
#' @param covariates Matrix of potential confounders.
#' @param family_y Specify the glm family for y. Only "binomial" (for binary y) and "gaussian" (for continuous y) are supported.
#' @param method_y Prediction model for y. Currently options "glmnet", "xgbTree", "rf", and "svm" are supported.
#' @param m Number of random CRT samples to generate for computing p-values.
#' @param cv_folds Number of CV folds in the modeling for y.
#' @param epsilon Pseudo count to add to relative abundances before transformation.
#' @param verbose Whether or not to print verbose information.
#' @param debug_dir Directory to save intermediate output (for debugging purpose).
#'
#' @return A matrix of two columns. Each row corresponds to one microbe's testing results. First column indicates conditional associations with the outcome. Second column indicates p-values.
#' @export
CompDA <-
  function(x, y,
           covariates = NULL,
           family_y = "gaussian",
           method_y = "glmnet",
           m = 1e4,
           cv_folds = 5,
           epsilon = NULL,
           verbose = TRUE,
           debug_dir = NULL) {

    if(!(all(x <= 1) & all( x>= 0)))
      stop("Currently only supports relative abundance values for x!")
    if(!(nrow(x) == length(y)))
      stop("Numbers of observations in x and y do not agree!")
    if(!is.null(covariates))
      if(!(nrow(x) == nrow(covariates)))
        stop("Covariates should be provided as a matrix, with rows corresponding to samples!")

    if(!(family_y %in% c("gaussian", "binomial")))
      stop("Currently only gaussian and binomial families for y are accepted!")
    if(!(method_y %in% c("glmnet", "xgbTree", "rf", "svm")))
       stop("Currently only supports glmnet, xgbTree, rf, and svm prediction models!")

    if(is.null(epsilon))
      epsilon <- min(setdiff(x, 0)) / 2

    if(!exists(".Random.seed")) set.seed(NULL)
    old <- .Random.seed
    on.exit({ .Random.seed <<- old })
    set.seed(0)

    if(!is.null(debug_dir))
      dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)
    features <- seq(1, ncol(x))

    progressr::with_progress({
      p <- progressr::progressor(along = features)
      results <- future.apply::future_vapply(
        seq(features),
        function(i_feature) {
          p(sprintf("i_feature=%g", i_feature))
          i_x <- x[, i_feature]
          i_z <- x[, -i_feature, drop = FALSE]
          i_z <- t(apply(i_z, 1, tss_withzero))
          i_z <- log(i_z + epsilon)
          if(!is.null(covariates))
            i_z <- cbind(i_z, covariates)

          i_fit <-
            test_oneX(
              y = y, x = i_x, z = i_z, epsilon = epsilon,
              family_y = family_y, method_y = method_y, cv_folds = cv_folds, m = m,
              debug_file = {
                if(!is.null(debug_dir))
                  paste0(debug_dir, "/", i_feature, ".RData")
                else
                  NULL
              })
          return(i_fit)
        },
        FUN.VALUE = c(0.0, 0.0),
        future.seed = TRUE) %>% t()
    },
    handlers = progressr::handlers("txtprogressbar"),
    enable = verbose)

    colnames(results) <- c("Conditional Assoc.", "p")
    if(!is.null(colnames(x)))
      rownames(results) <- colnames(x)
    else
      rownames(results) <- seq(1, ncol(x))
    return(results)
  }
