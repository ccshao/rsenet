dcca <- function(x, y, k) {
  #- featurs in rows, cells/bins in columns.
  x <- scale(x)
  y <- scale(y)

  res <- irlba(crossprod(x, y), k)

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
  knn_res <- nn2(cca_res$x, cca$y, num_bin)

  #- cells by bins matrix.
  mtx <- matrix(0, nrow = ncol(y), ncol = ncol(x))

  for (i in seq_len(nrow(mtx))) {
    mtx[i, knn_res$nn.idx[i, ]] <- knn_res$nn.dists[i, ]
  }

  return(mtx)
}

