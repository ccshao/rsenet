#' Perform elastic net regression wth scRNA-seq and spatial reference data based on k nearest neighbors found in CCA
#'
#' This function takes two input matrix, x and y, to  performs elastic net regression (via \code{glmnet}) on each column
#' of y (cells) with the x used as the predictors. Column of x are pre-selected by canonical correlation analysis.
#'
#' @param x a feature by bins matrix used as predictors.
#' @param y a feature by cells matrix while each column will be used as response variables separately in separate lasso regressions.
#'        The row names of y MUST be identical to those in x. If family = "gaussian" (default), it might be better to transform counts via logarithm (e.g., log2)
#' @param adaptive whether using adaptive process (Default TRUE).
#' @param hybrid whether to elastic net or lasso (Default TRUE).
#' @param nfolds the number of folds used in cross validation to estimate lambda (lambda and alpha in elastic net).
#' @param n_run times of regresss to repeat per cell.
#' @param num_cc number of canonical correlation vector to use.
#' @param num_bin number of nearest bins to select.
#' @param alpha list of alpha value used in evaluation. Note high alpha value tends to lasso regression.
#' @param ... additional parameters to \code{glmnet::cv.glmnet} or \code{glmnetUtils::glmnetUtils}.
#' @details
#'     Similar to \code{enet}, but with bins in x are pre-selected by diagonal CCA. The euclidean distance between cells and nearest bins
#'     are used as penalty (\code{adaptive} is TRUE). See Tim Stuart, et al, 2019; Andrew Butler, et al, 2018 for detail on CCA.
#'     The calculation is much faster than \code{enet} due to no individual estimation of penalty via ridge regression per cell.
#' @return A probility marix with bins in rownames and cells in columns, suggesting the probility of a cell assigned to a bin.
#' @export
enet_cca <- function(x,
  y,
  adaptive = TRUE,
  hybrid   = TRUE,
  nfolds   = 10,
  n_run    = 5,
  num_cc   = 20,
  num_bin  = 50,
  alpha    = seq(0.2, 0.8, by = .1),
  ...) {
  stopifnot(identical(rownames(x), rownames(y)))
  stopifnot(is(x, "sparseMatrix") || is.matrix(x))
  stopifnot(is(y, "sparseMatrix") || is.matrix(y))

  i <- NULL

  cca_mtx <- knn_p(x, y, num_cc, num_bin)

  #- Equal penalty.
  if (!adaptive) cca_mtx[cca_mtx != -1] <- 1

  p <- progressr::progressor(steps = ncol(y))

  fn_fit <- fn_hybrid_switch(hybrid)

  res <- foreach(i = seq_len(ncol(y))) %dorng% {
    p("Running ...")
    #- Vector with -1 and eucl dist.
    cca_nn <- cca_mtx[, i]

    #- Reduced x.
    x_new     <- x[, cca_nn != -1]
    penalty_f <- cca_nn[cca_nn != -1]

    res_one <- foreach(j = seq_len(n_run)) %do% {
      bin_res <- fn_fit(x_new, y, i, nfolds, penalty_f, alpha, ...)
      return(list(bin_res[[1]] / sum(bin_res[[1]]), bin_res[[2]]))
    }

    #- A bins by runs matrix, return normalized weights of bins to one cell.
    avg_one   <- do.call(cbind, lapply(res_one, "[[", 1)) %>% rowMeans(na.rm = TRUE)
    alpha_one <- sapply(res_one, "[[", 2) %>% mean

    bin_vec <- rep(0, ncol(x)) %>% set_names(colnames(x))
    bin_vec[names(avg_one)] <- avg_one
    return(list(bin_vec, alpha_one))
  }

  #- A bins by cells matrix.
  res_pred  <- do.call(cbind, lapply(res, "[[", 1)) %>% set_colnames(colnames(y))
  alpha_all <- sapply(res, "[[", 2)

  write.table(alpha_all, "alpha.txt", row.names = FALSE, col.names = FALSE)
  return(res_pred)
}

knn_p <- function(x, y, num_cc = 20, num_bin = 50) {
  cca_res <- dcca(x, y, num_cc)
  knn_res <- RANN::nn2(cca_res$x, cca_res$y, num_bin)

  #- cells by bins matrix.
  mtx <- matrix(-1, nrow = ncol(y), ncol = ncol(x)) %>%
    set_rownames(colnames(y)) %>%
    set_colnames(colnames(x))

  for (i in seq_len(nrow(mtx))) {
    mtx[i, knn_res$nn.idx[i, ]] <- knn_res$nn.dists[i, ]
  }

  #- A bins by cells matrix.
  return(t(mtx))
}

#- Diagonal CCA, codes are modified from Seurat.
dcca <- function(x, y, k) {
  #- featurs in rows, cells/bins in columns.
  x <- scale(x)
  y <- scale(y)

  res <- irlba::irlba(crossprod(x, y), k)

  #- Location of cells/bins in latent clusters, cells by features.
  latent_x <- res$u %>% set_colnames(paste0("k_", seq_len(k))) %>% set_rownames(colnames(x))
  latent_y <- res$v %>% set_colnames(paste0("k_", seq_len(k))) %>% set_rownames(colnames(y))

  #- Make the u and v unique to sign.
  sign_uv  <- apply(latent_x, 2, \(x) sign(x[1]))
  latent_x <- sweep(latent_x, 2, sign_uv, "*")
  latent_y <- sweep(latent_y, 2, sign_uv, "*")

  #- L2 norm per cells across latent dimensions.
  norm_x <- sweep(latent_x, 1, apply(latent_x, 1, \(x) sqrt(sum(x^2))), "/")
  norm_y <- sweep(latent_y, 1, apply(latent_y, 1, \(x) sqrt(sum(x^2))), "/")
  return(list(x = norm_x, y = norm_y))
}
