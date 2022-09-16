#' Perform elastic net regression with scRNA-seq and spatial reference data
#'
#' This function takes two input matrix, x and y, to  performs elastic net regression (via \code{glmnet}) on each column
#' of y (cells) with the x used as the predictors.
#'
#' @param x a feature by bins matrix used as predictors.
#' @param y a feature by cells matrix while each column will be used as response variables separately in separate lasso regressions.
#'        The row names of y MUST be identical to those in x. If family = "gaussian" (default), it might be better to transform counts via logarithm (e.g., log2)
#' @param adaptive whether using the adaptive process (Default TRUE).
#' @param hybrid whether to elastic net or lasso (Default TRUE).
#' @param tau cofficient in transforming the weights to penalty as abs(weights)^(-tau). 1 by default, other choices could be 0.5, 2.
#' @param nfolds the number of folds used in cross validation to estimate lambda (lambda and alpha in elastic net).
#' @param n_run times of regresss to repeat per cell.
#' @param alpha list of alpha value used in evaluation. Note high alpha value tends to lasso regression.
#' @param ... additional parameters to \code{glmnet::cv.glmnet} or \code{glmnetUtils::glmnetUtils}.
#' @details
#'     For each sample in the y, \code{enet} tries to find the weights of bins in x using non-negative elastic net. If \code{adaptive} is TRUE,
#'     ridge regression will be first performed to estimate the \code{penalty.factor}. If \code{hybrid} is TRUE,
#'     elastic net regression will be used with alpha and lambda are both estimated via cross validation.
#'     The estmatied coefficients are normalized by the sum, and averaged value from repeated runs are returned.
#' @return A probility marix with bins in rownames and samples in columns, suggesting the probility of a sample assigned to a bin.
#' @export
enet <- function(x, y, adaptive = TRUE, hybrid = TRUE, tau = 1, nfolds = 10, n_run = 10, alpha = seq(0.2, 1.0, by = .1), ...) {
  stopifnot(identical(rownames(x), rownames(y)))
  stopifnot(is(x, "sparseMatrix") || is.matrix(x))
  stopifnot(is(y, "sparseMatrix") || is.matrix(y))

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
      bin_coef  <- fn_fit(x, y, i, nfolds, penalty_f, alpha, ...)

      return(bin_coef / sum(bin_coef))
    }

    #- A bins by runs matrix, return normalized weights of bins to one cell.
    avg_one <- do.call(cbind, res_one) %>% rowMeans(na.rm = TRUE)
    return(avg_one)
  }

  #- A bins by cells matrix.
  res <- do.call(cbind, res) %>% set_colnames(colnames(y))
  return(res)
}

#- Swith adaptive.
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

fn_hybrid_switch <- function(hybrid = TRUE) {
  if (hybrid) {
    fn_enet
  } else {
    fn_lasso
  }
}

fn_enet <- function(x, y, i, nfolds, penalty_f, alpha, ...) {
  enet_cv   <- glmnetUtils::cva.glmnet(x, y[, i], nfolds = nfolds, lower.limits = 0, penalty.factor = penalty_f, alpha = alpha, ...)
  #- Find the alpha with lowest cvm of lambda.min.
  opt_alpha <- enet_cv$alpha[which.min(vapply(enet_cv$modlist, \(x) min(x$cvm), numeric(1)))]
  write(opt_alpha, "alpha.txt", append = TRUE)
  bin_coef  <- coef(enet_cv, alpha = opt_alpha, s = "lambda.1se")[-1, 1]
}

fn_lasso <- function(x, y, i, nfolds, penalty_f, alpha, ...) {
  lasso_cv <- glmnet::cv.glmnet(x, y[, i], nfolds = nfolds, lower.limits = 0, alpha = 1, penalty.factor = penalty_f, ...)
  bin_coef <- coef(lasso_cv, s = "lambda.1se")[-1, 1]
}
