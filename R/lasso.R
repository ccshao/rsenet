#' Perform elastic net regression on scRNA-seq and spatial reference data
#'
#' This function takes two input matrix, x and y,  performs elastic net regression (via \code{glmnet}) on each column
#' of y with the x used the predictors.
#'
#' @param x a feature by sample matrix used as predictors.
#' @param y a feature by sample response matrix while each column will be used as response variables in separate lasso regressions.
#'        The row names of y MUST be identical to those in x. If family = "gaussian", it might be better to transform counts via logarithm (e.g., log2)
#' @param adaptive whether using adaptive lasso (Default TRUE)
#' @param hybrid whether to elastic net (hybrid of lasso and ridge regression) or only lasso.
#' @param tau weights in transforming the penalty. 1 by default, other choices are 0.5, 2.
#' @param nfolds the number of folds used in cross validation to estimate lambda (and alpha if sepecified).
#' @param n_run number of regression to run per cell.
#' @param ... additional parameters to \code{glmnet::cv.glmnet} or \code{glmnetUtils::glmnetUtils}.
#' @details
#' For each sample in the y, this function tries to find the weights to the bins in x, with non-negative lasso.  If \code{adaptive} is TRUE,
#' ridge regression will be first perfomred to estimate the \code{penalty.factor}. Optionaly, the alpha value is estiamted via \code{glmnetUtils::cva.glmnet}.
#' lambda value will always be choosen from \code{nfolds} cross validation.
#' The estmatied coefficients are normalied to the sum, and the normalized value from \code{n_run} are averaged.
#' @return A probility marix with bins in rownames and samples in columns, suggesting the probility of a sample assigned to a bin.
#' @export
enet <- function(x, y, adaptive = TRUE, hybrid = TRUE, tau = 1, nfolds = 10, n_run = 10, ...) {
  i <- j <- NULL

  if (adaptive) {
    message("(II) Using adaptive procedure.")
  } else {
    message("(II) **NOT** using adaptive procedure.")
  }

  if (hybrid) {
    message("(II) Using elastic net.")
  } else {
    message("(II) Using lasso.")
  }

  p <- progressr::progressor(steps = ncol(y))

  fn_penalty <- fn_adaptive_switch(adaptive)
  fn_fit     <- fn_hybrid_switch(hybrid)

  res <- foreach(i = seq_len(ncol(y))) %dorng% {
    #- Prediction of one cell with multiple runs to get consensus results.
    p("Running ...")
    res_one <- foreach(j = seq_len(n_run)) %do% {
      penalty_f <- fn_penalty(x, y, i, tau, nfolds, ...)
      bin_coef  <- fn_fit(x, y, i, nfolds, penalty_f, ...)

      return(bin_coef / sum(bin_coef))
    }

    #- A bins by runs matrix.
    res_one <- fn_list2mtx(res_one, rnames = colnames(x))
    avg_one <- rowMeans(res_one, na.rm = TRUE)
    #- The normalized weights of bins to one cell.
    return(avg_one)
  }

  #- A bins by cells matrix.
  res <- fn_list2mtx(res, rnames = colnames(x), cnames = colnames(y))
  return(res)
}

#- Convert list to matrix.
fn_list2mtx <- function(x, rnames = NULL, cnames = NULL) {
  res <- matrix(unlist(x), ncol = length(x)) %>%
    set_rownames(rnames) %>%
    set_colnames(cnames)
  return(res)
}

#- Swith adaptive.
# fn_adaptive_switch <- function(adaptive = TRUE, ...) {
fn_adaptive_switch <- function(adaptive = TRUE) {
  if (adaptive) {
    fn_adaptive_p
  } else {
    fn_fix_p
  }
}

fn_adaptive_p <- function(x, y, i, tau, nfolds, ...) {
  penalty_f <- rep_len(1, ncol(x))

  #- Get the coefficients of ridge regression.
  ridge_cv   <- glmnet::cv.glmnet(x, y[, i], nfolds = nfolds, alpha = 0, lower.limits = 0, ...)
  ridge_coef <- coef(ridge_cv, s = "lambda.1se")[-1, 1]

  #- Replace 0 with minimum non-zero value, if all zero then taking 1.
  if (any(ridge_coef %!=% 0)) {
    ridge_coef[ridge_coef %==% 0] <- min(ridge_coef[ridge_coef %!=% 0])
    penalty_f <- abs(ridge_coef)^(-tau)
  }

  return(penalty_f)
}

fn_fix_p <- function(x, ...) return(rep_len(1, ncol(x)))

#- Has to add ellipsis otherwise raise error.
# fn_hybrid_switch <- function(hybrid = TRUE, ...) {
fn_hybrid_switch <- function(hybrid = TRUE) {
  if (hybrid) {
    fn_enet
  } else {
    fn_lasso
  }
}

fn_enet <- function(x, y, i, nfolds, penalty_f, ...) {
  enet_cv   <- glmnetUtils::cva.glmnet(x, y[, i], nfolds = nfolds, lower.limits = 0, penalty.factor = penalty_f, ...)
  #- Find the alpha with lowest cvm of lambda.min.
  opt_alpha <- enet_cv$alpha[which.min(vapply(enet_cv$modlist, \(x) min(x$cvm), numeric(1)))]
  bin_coef  <- coef(enet_cv, alpha = opt_alpha, s = "lambda.1se")[-1, 1]
}

fn_lasso <- function(x, y, i, nfolds, penalty_f, ...) {
  lasso_cv <- glmnet::cv.glmnet(x, y[, i], nfolds = nfolds, lower.limits = 0, alpha = 1, penalty.factor = penalty_f, ...)
  bin_coef <- coef(lasso_cv, s = "lambda.1se")[-1, 1]
}
