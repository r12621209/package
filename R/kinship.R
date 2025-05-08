
#' The function calculate additive matrix 
#'
#' major allele = 2, heterozygote = 1, minor allele = 0.
#'
#' @param X  A a set of SNP data.
#' @export
#' @examples
#' data(DST2_maize)
#' K.t2 <- kinship(snp_DST2_maize)
#'

kinship <- function(X){
  
  for (i in 1:ncol(X)) {
    sep <- table(X[, i])
    if(length(sep)==3){
      p <- (sep['2'] + (sep['1'] * 0.5)) / nrow(X)
    }else{
      p <- sep['2'] / nrow(X)
    }
    X[, i] <- X[, i] - (2 * p)
  }
  X=as.matrix(X)
  K=X%*%t(X)
  K=K/mean(diag(K))
  return(K)
}
