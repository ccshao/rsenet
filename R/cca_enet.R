#' Perform elastic net regression on scRNA-seq and spatial reference data based k nearest neighbors found in CCA sub-space.
#'
#' This function takes two input matrix, x and y, performs elastic net regression (via \code{glmnet}) on each column
#' of y with the k nearest in CCA sub-sapce x used as the predictors.
#'
#' @param x a feature by sample matrix used as predictors.
#' @param y a feature by sample response matrix while each column will be used as response variables in separate lasso regressions.
#'        The row names of y MUST be identical to those in x. If family = "gaussian", it might be better to transform counts via logarithm (e.g., log2)
#' @param adaptive whether using adaptive lasso (Default TRUE)
#' @param hybrid whether to elastic net (hybrid of lasso and ridge regression) or only lasso.
#' @param tau weights in transforming the penalty. 1 by default, other choices are 0.5, 2.
#' @param nfolds the number of folds used in cross validation to estimate lambda (and alpha if sepecified).
#' @param n_run number of regression to run per cell.
#' @param num_cc number of canonical correlation vector to use.
#' @param num_bin number of nearest bins to select.
#' @param ... additional parameters to \code{glmnet::cv.glmnet} or \code{glmnetUtils::glmnetUtils}.
#' @details
#'     Similar to \code{enet}, but with bins in x are pre-selected by CCA. Additially, the euclidean distance between cells and nearest bins
#'     are used as penalty.
#' @return A probility marix with bins in rownames and samples in columns, suggesting the probility of a sample assigned to a bin.
#' @export
enet_cca <- function(x, y, adaptive = TRUE, hybrid = TRUE, tau = 1, nfolds = 10, n_run = 5, num_cc = 20, num_bin = 50, ...) {
  i <- NULL

  cca_mtx <- knn_p(x, y, num_cc, num_bin)

  #- Equal penalty.
  if (!adaptive) cca_mtx[cca_mtx != -1] <- 1

  #- Skip the steps on calculating penalty.
  p <- progressr::progressor(steps = ncol(y))

  fn_fit <- fn_hybrid_switch(hybrid)

  res <- foreach(i = seq_len(ncol(y))) %dorng% {
    #- Prediction of one cell with multiple runs to get consensus results.
    p("Running ...")
    #- Vector with -1 and eucl dist.
    cca_nn <- cca_mtx[, i]

    x_new     <- x[, cca_nn != -1]
    penalty_f <- cca_nn[cca_nn != -1]

    res_one <- foreach(j = seq_len(n_run)) %do% {
      bin_coef <- fn_fit(x_new, y, i, nfolds, penalty_f, ...)
      return(bin_coef / sum(bin_coef))
    }

    #- A bins by runs matrix, return normalized weights of bins to one cell.
    avg_one <- do.call(cbind, res_one) %>% rowMeans(na.rm = TRUE)

    bin_vec <- rep(0, ncol(x)) %>% set_names(colnames(x))
    bin_vec[names(avg_one)] <- avg_one
    return(bin_vec)
  }

  #- A bins by cells matrix.
  res <- do.call(cbind, res) %>% set_colnames(colnames(y))
  return(res)
}

#- Diagonal CCA, codes are modified from Seurat.
dcca <- function(x, y, k) {
  #- featurs in rows, cells/bins in columns.
  x <- scale(x)
  y <- scale(y)

  res <- irlba::irlba(crossprod(x, y), k)

  #- location of cells/bins in latent clusters, cells by features.
  latent_x <- res$u %>% set_colnames(paste0("k_", seq_len(k))) %>% set_rownames(colnames(x))
  latent_y <- res$v %>% set_colnames(paste0("k_", seq_len(k))) %>% set_rownames(colnames(y))

  #- L2 norm per cells.
  norm_x <- sweep(latent_x, 1, apply(latent_x, 1, \(x) sqrt(sum(x^2))), "/")
  norm_y <- sweep(latent_y, 1, apply(latent_y, 1, \(x) sqrt(sum(x^2))), "/")
  return(list(x = norm_x, y = norm_y))
}

knn_p <- function(x, y, num_cc = 20, num_bin = 50) {
  cca_res <- dcca(x, y, num_cc)
  knn_res <- RANN::nn2(cca_res$x, cca_res$y, num_bin)

  #- cells by bins matrix.
  mtx <- matrix(-1, nrow = ncol(y), ncol = ncol(x))

  for (i in seq_len(nrow(mtx))) {
    mtx[i, knn_res$nn.idx[i, ]] <- knn_res$nn.dists[i, ]
  }

  #- A bins by cells matrix.
  return(t(mtx))
}
