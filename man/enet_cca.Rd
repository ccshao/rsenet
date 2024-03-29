% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enet_cca.R
\name{enet_cca}
\alias{enet_cca}
\title{Perform elastic net regression wth scRNA-seq and spatial reference data based on k nearest neighbors found in CCA}
\usage{
enet_cca(
  x,
  y,
  adaptive = TRUE,
  hybrid = TRUE,
  nfolds = 10,
  n_run = 5,
  num_cc = 20,
  num_bin = 50,
  alpha = seq(0.2, 0.8, by = 0.1),
  ...
)
}
\arguments{
\item{x}{a feature by bins matrix used as predictors.}

\item{y}{a feature by cells matrix while each column will be used as response variables separately in separate lasso regressions.
The row names of y MUST be identical to those in x. If family = "gaussian" (default), it might be better to transform counts via logarithm (e.g., log2)}

\item{adaptive}{whether using adaptive process (Default TRUE).}

\item{hybrid}{whether to elastic net or lasso (Default TRUE).}

\item{nfolds}{the number of folds used in cross validation to estimate lambda (lambda and alpha in elastic net).}

\item{n_run}{times of regresss to repeat per cell.}

\item{num_cc}{number of canonical correlation vector to use.}

\item{num_bin}{number of nearest bins to select.}

\item{alpha}{list of alpha value used in evaluation. Note high alpha value tends to lasso regression.}

\item{...}{additional parameters to \code{glmnet::cv.glmnet} or \code{glmnetUtils::glmnetUtils}.}
}
\value{
A probility marix with bins in rownames and cells in columns, suggesting the probility of a cell assigned to a bin.
}
\description{
This function takes two input matrix, x and y, to  performs elastic net regression (via \code{glmnet}) on each column
of y (cells) with the x used as the predictors. Column of x are pre-selected by canonical correlation analysis.
}
\details{
Similar to \code{enet}, but with bins in x are pre-selected by diagonal CCA. The euclidean distance between cells and nearest bins
    are used as penalty (\code{adaptive} is TRUE). See Tim Stuart, et al, 2019; Andrew Butler, et al, 2018 for detail on CCA.
    The calculation is much faster than \code{enet} due to no individual estimation of penalty via ridge regression per cell.
}
