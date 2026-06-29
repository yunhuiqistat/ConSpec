#' Contrastive PCA
#'
#' This function performs contrastive PCA (cPCA) or sparse contrastive PCA (scPCA), in comparison, also perform target/ancillary group only PCA/sPCA
#' @import PMA
#' @importFrom stats cov
#' @importFrom stats dist 
#' @importFrom stats kmeans 
#' @param Xt data frame for target group, samples in rows, variables in columns.
#' @param Xa data frame for ancillary group, samples in rows, variables in columns.
#' @param eta tuning parameter controlling how much ancillary variation should be contrasted off from the target variation. It can be theo \code{eta_opt} from function eta_tuning_general()
#' @param sparse_ctst parameter controlling the sparsity of scPCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_trt parameter controlling the sparsity of target sparse PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_ctrl parameter controlling the sparsity of ancillary sparse PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @export
#' @return A list containing contrastive scores \code{cPC} and loadings \code{cU}, target scores \code{PCt} and loadings \code{Ut}, and ancillary scores \code{PCa} and loadings \code{Ua}.
#' @examples
#' library(ConSpec)
#' data(intro)
#' cpca_res <- cPCA(Xt = intro$cPCA$Xt, Xa = intro$cPCA$Xa, eta = 1.2, sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL)


cPCA <- function(Xt, Xa, eta, sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL){

  ## pre-processing
  Xt <- scale(Xt, center = TRUE, scale = FALSE)
  Xa <- scale(Xa, center = TRUE, scale = FALSE)

  ## contrastive
  if(is.null(sparse_ctst)){
    # cPCA
    cU <- eigen(cov(Xt)-eta*cov(Xa))$vectors
  }
  else{
    # scPCA
    if(length(sparse_ctst) > 1){
      # use CV to find penalty parameters
      sumabs_ctst <- PMA::PMD.cv(cov(Xt)-eta*cov(Xa), type="standard",
                       sumabss = sparse_ctst, center = FALSE, trace = FALSE)$bestsumabs
    }
    else{sumabs_ctst <- sparse_ctst}
    cV <- PMA::PMD(cov(Xt)-eta*cov(Xa), K=5, type = "standard", sumabs = sumabs_ctst, center = FALSE, trace = FALSE)$v
    cU <- PMA::PMD(cov(Xt)-eta*cov(Xa), K=5, type = "standard", sumabs = sumabs_ctst, center = FALSE, trace = FALSE)$u
    cU <- cU[,order(diag(t(cU)%*%cV), decreasing = TRUE)[1:2]]
  }
  cPC <- Xt%*%cU

  ## target only
  if(!is.null(sparse_trt)){
    # use sPCA for target group
    if(length(sparse_trt) == 1){sumabs_trt <- sparse_trt} # with provided sparsity
    else{# cv to choose sparsity
      sumabs_trt <- PMA::PMD.cv(cov(Xt), type="standard", sumabss = sparse_trt, center = FALSE, trace = FALSE)$bestsumabs
    }
    Ut <- PMA::PMD(cov(Xt), K=5, type = "standard", sumabs = sumabs_trt, center = FALSE, trace = FALSE)$v}
  else{Ut <- eigen(cov(Xt))$vectors}  # use PCA for target

  PCt <- Xt%*%Ut

  ## ancillary only
  if(!is.null(sparse_ctrl)){
    # use sPCA for ancillary group
    if(length(sparse_ctrl) == 1){sumabs_ctrl <- sparse_ctrl} # with provided sparsity
    else{# cv to choose sparsity
      sumabs_ctrl <- PMA::PMD.cv(cov(Xa), type="standard", sumabss = sparse_ctrl, center = FALSE, trace = FALSE)$bestsumabs
    }
    Ua <- PMA::PMD(cov(Xa), K=5, type = "standard", sumabs = sumabs_ctrl, center = FALSE, trace = FALSE)$v}
  else{Ua <- eigen(cov(Xa))$vectors}  # use PCA for control

  PCa <- Xa%*%Ua

  return(list(PCt = PCt, Ut = Ut, PCa = PCa, Ua = Ua, cPC = cPC, cU = cU))
}
