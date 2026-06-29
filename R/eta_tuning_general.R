#' Select \eqn{\eta} for contrastive spectral clustering by cPCA or scPCA
#' @import PMA
#' @import ggplot2
#' @import ggsci
#' @importFrom aricode ARI
#' @importFrom cluster silhouette
#' @importFrom stats cov
#' @importFrom stats dist 
#' @importFrom stats kmeans 
#' @param Xt data frame for target group, samples in rows, variables in columns. If PCA is applied to correlation matrix, Xt should be scaled and centered.
#' @param Xa data frame for ancillary group, samples in rows, variables in columns. If cPCA is applied to correlation matrix, Xa should also be scaled and centered.
#' @param eta_range a vector containing candidate eta range to be tuning from.
#' @param plot a logical quantity indicating whether plot the ARI(nuisance) and silhouette score.
#' @param sparse_trt parameter controlling the sparsity of target sparse PCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_ctst parameter controlling the sparsity of scPCA, if NULL, use cPCA, else use scPCA. For use of scPCA, if vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param num_comps if scPCA, specify the number of contrastive components that will be estimated by PMD, default is 20.
#' @param km_cluster number showing number of nuisance clusters defined from kmeans clustering on target PCA
#' @param ckm_cluster number of interesting clusters defined from kmeans clustering on cPCA/scPCA
#' @export
#' @return A list containing selected value of \eqn{\eta}  \code{eta_opt}, the vector of ARI(nuisance) \code{nuisance_ARI} and the vector of Silhouette \code{silhouette}, and the ggplot object \code{p}.
#' @examples
#' library(ConSpec)
#' # A simple and balanced data example
#' data(intro)
#' eta_range <- seq(0.1, 5, 0.05)
#' eta_opt <- eta_tuning_general(Xt = intro$cPCA$Xt, Xa = intro$cPCA$Xa, eta_range = eta_range, plot = TRUE, sparse_trt = NULL, sparse_ctst = NULL)$eta_opt
#' # An unbalanced data example
#' data(sim_data_CSC)
#' eta_tuning_general(Xt = sim_data_CSC$Xt, Xa = sim_data_CSC$Xa, eta_range = seq(0, 100, 0.2))

eta_tuning_general <- function(Xt, Xa, eta_range = seq(0.5, 10, 0.5), plot = TRUE,
                       sparse_trt = NULL, sparse_ctst = NULL, num_comps = 20,
                       km_cluster = 2, ckm_cluster = 2){

  ## pre-processing
  Xt <- scale(Xt, center = TRUE, scale = FALSE)
  Xa <- scale(Xa, center = TRUE, scale = FALSE)

  ## trtPCA/sPCA to define nuisance labels
  if(!is.null(sparse_trt)){
    # use sPCA
    if(length(sparse_trt) == 1){sumabs_trt <- sparse_trt} # with provided sparsity
    else{# cv to choose sparsity
      sumabs_trt <- PMD.cv(cov(Xt), type="standard", sumabss = sparse_trt, center = FALSE, trace = FALSE)$bestsumabs
    }
    Ut <- PMD(cov(Xt), K=2, type = "standard", sumabs = sumabs_trt, center = FALSE, trace = FALSE)$v}
  else{Ut <- eigen(cov(Xt))$vectors[,1:2]} # use PCA

  PCt <- Xt%*%Ut
  km <- kmeans(PCt[,1], centers = km_cluster)
  nuisance_labels <- km$cluster

  ## cPCA/scPCA for a range of eta
  nuisance_ARI <- c()
  silhouette <- c()
  for (eta in eta_range){
    if(!is.null(sparse_ctst)){ # use scPCA
      if(length(sparse_ctst) == 1){sumabs_ctst <- sparse_ctst}
      else{# cv to choose best sparsity
        sumabs_ctst <- PMD.cv(cov(Xt)-eta*cov(Xa), type="standard",
                              sumabss = sparse_ctst, center = FALSE, trace = FALSE)$bestsumabs
      }
      # since the contrast matrix is symmetric, its fine to use either cU or cV, they are the same in theory, close in practice.
      cU <- PMD(cov(Xt)-eta*cov(Xa), type="standard", sumabs=sumabs_ctst, K=num_comps, center = FALSE, trace = FALSE)$u
      cV <- PMD(cov(Xt)-eta*cov(Xa), type="standard", sumabs=sumabs_ctst, K=num_comps, center = FALSE, trace = FALSE)$v
      cU <- matrix(cU[,which.max(diag(t(cU)%*%cV))], ncol=1)
    }
    else{cU <- eigen(cov(Xt)-eta*cov(Xa))$vectors[,1:2]} # use cPCA
    cPC <- Xt%*%cU
    ckm <- kmeans(cPC[,1], centers = ckm_cluster)
    interesting_labels <- ckm$cluster
    silhouette_score <- silhouette(interesting_labels, dist(cPC[,1]))
    silhouette <- c(silhouette, mean(silhouette_score[, "sil_width"]))
    nuisance_ARI <- c(nuisance_ARI, ARI(nuisance_labels, interesting_labels))
  }
  eta_opt <- eta_range[which.max(silhouette - nuisance_ARI)]
  p <- NULL
  if(plot){
    plot_df <- data.frame(eta = rep(eta_range,2), metrics = c(nuisance_ARI, silhouette), type = rep(c("ARI","Silhouette"), each = length(eta_range)))
    p <- ggplot(plot_df)+
      geom_point(aes(x = eta, y = metrics, color = type, shape = type))+
      theme_bw()+
      geom_vline(xintercept = eta_opt, linetype = "dashed")+
      labs(color = "", shape = "", title = "", x = expression(eta), y = "")+
      scale_shape_manual(values=c(15,16))+
      scale_color_cosmic()+
      theme(
        text = element_text(size = 14),  # Adjust text size
        axis.title = element_text(size = 16),  # Adjust axis title size
        axis.text = element_text(size = 12),   # Adjust axis text size
        legend.position = "bottom")
    print(p)
  }
  return(list(eta_opt = eta_opt, eta_range = eta_range, nuisance_ARI = nuisance_ARI, silhouette = silhouette, plot=p))
}
