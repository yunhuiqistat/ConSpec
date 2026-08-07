#' Simulate Sparse Contrastive Cross Covariance Analysis Data
#'
#' Simulates high-dimensional data for a sparse contrastive cross covariance analysis (s3CA) 
#' scenario with target and ancillary groups. Evaluates 
#' both clustering performance and feature selection accuracy using standard sparse CCA 
#' and contrastive sparse CCA approaches. The continuous nuisance variation is also 
#' simulated as a sparse signal to adequately challenge the feature selection algorithm.
#'
#' @param nt Integer. Number of target samples (default 200).
#' @param na Integer. Number of ancillary samples (default 200).
#' @param p Integer. Total number of features in modality X (default 300).
#' @param q Integer. Total number of features in modality Y (default 300).
#' @param s Integer. Number of true active signal features in each modality (sparsity level, default 20).
#' @param mu Numeric. Signal strength parameter for the cluster means (default 1.0).
#' @param theta Numeric. Variance of the dominant continuous nuisance factor (default 80).
#' @param eta_true Numeric. True scaling parameter for the nuisance variation between target and ancillary groups (default 2.5).
#' @param repetition Integer. Number of simulation repetitions (default 50).
#' @param eta Contrast parameter rule. A numeric vector of length > 1 is the candidate grid
#'   handed to \code{\link{eta_tuning_bimodal}}, which is the default. A single number is used
#'   directly. \code{"ratio"} instead uses the spectral-scale estimator
#'   \eqn{\lambda_1(S_t)/\lambda_1(S_a)}, which is recorded in the \code{Eta_Ratio} column of
#'   the output either way.
#' @param km_cluster Integer. Number of reference (nuisance) clusters used by the tuning algorithm.
#' @param ckm_cluster Integer. Number of contrastive clusters used by the tuning algorithm.
#'
#' @return A data frame with columns: Iteration, ARI_Target, ARI_Contrastive, TPR_Target, FPR_Target, TPR_Contrastive, FPR_Contrastive, Eta_Est, Eta_Ratio.
#'
#' @importFrom mclust adjustedRandIndex
#' @importFrom stats cov rnorm model.matrix kmeans
#' @importFrom PMA PMD
#' @export
#' @examples
#' # See vignettes/ContrastiveSpectralClustering.Rmd for example usage.
sim_s3CA <- function(nt = 200, na = 200, p = 300, q = 300, s = 20, mu = 1.0,
                             theta = 80, eta_true = 2.5, repetition = 50,
                             eta = seq(0.1, 6, 0.1), km_cluster = 2, ckm_cluster = 2){

  res_df <- data.frame()

  # true sparsity handed to PMD. PMD's scalar `sumabs` maps to sumabsu = sumabs*sqrt(nrow)
  # and sumabsv = sumabs*sqrt(ncol); with p == q this reproduces sumabsu = sumabsv = sqrt(s),
  # i.e. exactly the true sparsity used for the analyses below.
  sumabs_true <- sqrt(s / p)
  
  M <- matrix(0, nrow = 2, ncol = p)
  N <- matrix(0, nrow = 2, ncol = q)
  M[1, 1:s] <- mu; M[2, 1:s] <- -mu
  N[1, 1:s] <- mu; N[2, 1:s] <- -mu
  
  wx <- rep(0, p); wy <- rep(0, q)
  wx[(s+1):(2*s)] <- runif(s, 0.5, 1) 
  wy[(s+1):(2*s)] <- runif(s, 0.5, 1)
  wx <- wx / sqrt(sum(wx^2)); wy <- wy / sqrt(sum(wy^2))
  
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
    
    # spectral-scale estimator, recorded in every run for reference
    eta_ratio <- svd(cov(Xt, Yt))$d[1] / svd(cov(Xa, Ya))$d[1]

    # contrast parameter actually used
    if(identical(eta, "ratio")){eta_est <- eta_ratio}
    else if(length(eta) > 1){
      eta_est <- eta_tuning_bimodal(Xt = Xt, Yt = Yt, Xa = Xa, Ya = Ya,
                                    eta_range = eta, plot = FALSE,
                                    sparse_ctst = sumabs_true, sparse_trt = sumabs_true,
                                    reference = "target", num_comps = 1,
                                    km_cluster = km_cluster,
                                    ckm_cluster = ckm_cluster)$eta_opt
    }
    else{eta_est <- eta}

    sCCA_trt <- PMD(cov(Xt, Yt), type = "standard", center = FALSE,
                    sumabsu=sqrt(s), sumabsv=sqrt(s), K=1, trace=FALSE)
    u_trt <- as.numeric(sCCA_trt$u); v_trt <- as.numeric(sCCA_trt$v)
    
    score_trt <- cbind(Xt %*% u_trt, Yt %*% v_trt)
    ari_trt <- adjustedRandIndex(kmeans(score_trt, centers = 2, nstart = 10)$cluster, labels)
    
    tpr_trt <- (sum(u_trt[1:s] != 0) + sum(v_trt[1:s] != 0)) / (2*s)
    fpr_trt <- (sum(u_trt[(s+1):p] != 0) + sum(v_trt[(s+1):q] != 0)) / (p + q - 2*s)
    
    C <- cov(Xt, Yt) - eta_est * cov(Xa, Ya)
    sCCA_ctst <- PMD(C, type = "standard", center = FALSE, 
                     sumabsu=sqrt(s), sumabsv=sqrt(s), K=1, trace=FALSE)
    u_ctst <- as.numeric(sCCA_ctst$u); v_ctst <- as.numeric(sCCA_ctst$v)
    
    score_ctst <- cbind(Xt %*% u_ctst, Yt %*% v_ctst)
    ari_ctst <- adjustedRandIndex(kmeans(score_ctst, centers = 2, nstart = 10)$cluster, labels)
    
    tpr_ctst <- (sum(u_ctst[1:s] != 0) + sum(v_ctst[1:s] != 0)) / (2*s)
    fpr_ctst <- (sum(u_ctst[(s+1):p] != 0) + sum(v_ctst[(s+1):q] != 0)) / (p + q - 2*s)
    
    res_df <- rbind(res_df, data.frame(
      Iteration = repe,
      ARI_Target = ari_trt, ARI_Contrastive = ari_ctst,
      TPR_Target = tpr_trt, FPR_Target = fpr_trt,
      TPR_Contrastive = tpr_ctst, FPR_Contrastive = fpr_ctst,
      Eta_Est = eta_est, Eta_Ratio = eta_ratio
    ))
  }
  return(res_df)
}