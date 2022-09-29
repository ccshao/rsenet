#' @import Matrix
#' @importFrom magrittr %>% set_rownames set_colnames set_names
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom fpCompare %==% %!=%
#' @importFrom doRNG %dorng%
#' @importFrom stats coef
#' @importFrom methods is
#' @importFrom utils write.table
NULL

#- deal with . in magrittr
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
