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

