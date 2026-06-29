#' Simulate Contrastive Cross Covariance Analysis Data
#'
#' Simulates data for a contrastive cross covariance analysis scenario. Computes clustering performance using standard CCA and
#' contrastive CCA approaches.
#'
#' @param nt Integer. Number of target samples (default 200).
#' @param na Integer. Number of ancillary samples (default 200).
#' @param p Integer. Number of X features (default 100).
#' @param q Integer. Number of Y features (default 100).
#' @param mu Numeric. Signal strength parameter (default 1.0).
#' @param theta Numeric. Noise variance parameter (default 400).
#' @param eta_true Numeric. True signal-to-noise ratio (default 2.5).
#' @param repetition Integer. Number of simulation repetitions (default 50).
#'
#' @return A data frame with columns: Iteration, ARI_Target, ARI_Contrastive, Eta_Est
#'
#' @importFrom mclust adjustedRandIndex
#' @importFrom stats cov
#' @importFrom stats dist 
#' @importFrom stats kmeans 
#' @export
#' @examples
#' # See vignettes/ContrastiveSpectralClustering.Rmd for example usage.
sim_3CA <- function(nt = 200, na = 200, p = 100, q = 100, mu = 1.0, 
                            theta = 400, eta_true = 2.5, repetition = 50){
  res_df <- data.frame()
  
  # REMOVED SCALING: Just use the raw mu for every feature
  M <- matrix(rep(c(mu, -mu), each = p), nrow = 2, byrow = TRUE)
  N <- matrix(rep(c(mu, -mu), each = q), nrow = 2, byrow = TRUE)
  
  # Dense nuisance directions
  wx <- runif(p, -1, 1); wx <- wx / sqrt(sum(wx^2))
  wy <- runif(q, -1, 1); wy <- wy / sqrt(sum(wy^2))
  
  for (repe in 1:repetition) {
    labels <- sample(1:2, nt, replace = TRUE)
    L <- unname(model.matrix(~ factor(labels) - 1))
    
    Za <- rnorm(na, 0, sqrt(theta))
    Zt <- rnorm(nt, 0, sqrt(eta_true * theta))
    
    Xa <- Za %*% t(wx) + matrix(rnorm(na * p), na, p)
    Ya <- Za %*% t(wy) + matrix(rnorm(na * q), na, q)
    
    Xt <- L %*% M + Zt %*% t(wx) + matrix(rnorm(nt * p), nt, p)
    Yt <- L %*% N + Zt %*% t(wy) + matrix(rnorm(nt * q), nt, q)
    
    Xt <- scale(Xt, scale=FALSE); Yt <- scale(Yt, scale=FALSE)
    Xa <- scale(Xa, scale=FALSE); Ya <- scale(Ya, scale=FALSE)
    
    eta_est <- svd(cov(Xt, Yt))$d[1] / svd(cov(Xa, Ya))$d[1]
    
    CCA_trt <- svd(cov(Xt, Yt))
    score_trt <- cbind(Xt %*% CCA_trt$u[,1], Yt %*% CCA_trt$v[,1])
    ari_trt <- adjustedRandIndex(kmeans(score_trt, centers = 2, nstart = 10)$cluster, labels)
    
    C <- cov(Xt, Yt) - eta_est * cov(Xa, Ya)
    CCA_ctst <- svd(C)
    score_ctst <- cbind(Xt %*% CCA_ctst$u[,1], Yt %*% CCA_ctst$v[,1])
    ari_ctst <- adjustedRandIndex(kmeans(score_ctst, centers = 2, nstart = 10)$cluster, labels)
    
    res_df <- rbind(res_df, data.frame(
      Iteration = repe, ARI_Target = ari_trt, ARI_Contrastive = ari_ctst, Eta_Est = eta_est
    ))
  }
  return(res_df)
}


