#' Select \eqn{\eta} for cross-modality contrastive spectral clustering by 3CA or s3CA
#'
#' Bimodal counterpart of \code{\link{eta_tuning_general}}. The rationale is the
#' same as in the single-modality algorithm: a feasible \eqn{\eta} should produce
#' a contrastive arrangement of the target samples that (i) disagrees with the
#' reference (nuisance) arrangement obtained without contrast, measured by a small
#' ARI, and (ii) is itself well separated, measured by a large average silhouette
#' score. The selected value maximizes \eqn{\mbox{Sil}(\eta) - \mbox{ARI}(\eta)}.
#'
#'
#' @import PMA
#' @import ggplot2
#' @import ggsci
#' @importFrom aricode ARI
#' @importFrom cluster silhouette
#' @importFrom stats cov
#' @importFrom stats dist
#' @importFrom stats kmeans
#' @param Xt data frame for one data type in target group, samples in rows, variables in columns. To tune based on cross correlation matrices, each feature should be centered and standardized.
#' @param Yt data frame for another data type in target group, samples in rows, variables in columns.
#' @param Xa data frame for one data type in ancillary group, samples in rows, variables in columns.
#' @param Ya data frame for another data type in ancillary group, samples in rows, variables in columns.
#' @param eta_range a vector containing candidate eta range to be tuned from.
#' @param plot a logical quantity indicating whether plot the ARI(nuisance) and silhouette score.
#' @param sparse_ctst parameter controlling the sparsity of s3CA, if NULL, use 3CA, else use s3CA. If vector, refer sumabss in PMD.cv from package PMA, if single number, refer sumabs in PMD in package PMA.
#' @param sparse_trt parameter controlling the sparsity of the target-only (sparse) cross covariance analysis used to define the reference labels when \code{reference = "target"}. NULL means non-sparse.
#' @param sparse_ctrl parameter controlling the sparsity of the ancillary-only (sparse) cross covariance analysis used to define the reference labels when \code{reference = "ancillary"}. NULL means non-sparse.
#' @param reference which non-contrastive analysis defines the reference (nuisance) arrangement of the target samples. \code{"target"} (default) mimics \code{eta_tuning_general} and clusters the target samples on the target-only leading component pair. \code{"ancillary"} clusters the target samples after projecting them on the ancillary-only leading component pair, which targets the nuisance direction directly.
#' @param num_comps number of leading contrastive singular pairs used to build the score matrix for clustering. Default is 1, i.e. the top pair, giving two score columns (\eqn{X_t\hat u_1}, \eqn{Y_t\hat v_1}).
#' @param km_cluster number of nuisance clusters defined from kmeans clustering on the reference score pair.
#' @param ckm_cluster number of interesting clusters defined from kmeans clustering on the contrastive score pair.
#' @param nstart number of random starts passed to \code{\link[stats]{kmeans}}.
#' @export
#' @return A list containing the selected value of \eqn{\eta} \code{eta_opt}, the candidate grid \code{eta_range}, the vector of ARI(nuisance) \code{nuisance_ARI}, the vector of silhouette scores \code{silhouette}, the reference labels \code{nuisance_labels}, and the ggplot object \code{plot}.
#' @examples
#' library(ConSpec)
#' data(intro)
#' eta_tuning_bimodal(Xt = intro$cCCA$Xt, Yt = intro$cCCA$Yt,
#'                    Xa = intro$cCCA$Xa, Ya = intro$cCCA$Ya,
#'                    eta_range = seq(0.1, 5, 0.1), plot = TRUE)

eta_tuning_bimodal <- function(Xt, Yt, Xa, Ya, eta_range = seq(0.5, 10, 0.5), plot = TRUE,
                               sparse_ctst = NULL, sparse_trt = NULL, sparse_ctrl = NULL,
                               reference = c("target", "ancillary"),
                               num_comps = 1,
                               km_cluster = 2, ckm_cluster = 2, nstart = 25){

  reference <- match.arg(reference)

  ## pre-processing: center only, never scale (scaling, if wanted, is done by the caller)
  Xt <- scale(Xt, center = TRUE, scale = FALSE)
  Yt <- scale(Yt, center = TRUE, scale = FALSE)
  Xa <- scale(Xa, center = TRUE, scale = FALSE)
  Ya <- scale(Ya, center = TRUE, scale = FALSE)

  ## helper: leading num_comps singular pairs of a cross-covariance matrix,
  ## sparse (PMD) if a sparsity parameter is supplied, dense (SVD) otherwise.
  lead_pair <- function(Omega, sparse, K){
    if(is.null(sparse)){ # dense SVD
      sv <- svd(Omega)
      return(list(u = sv$u[, 1:K, drop = FALSE], v = sv$v[, 1:K, drop = FALSE]))
    }
    if(length(sparse) == 1){sumabs_use <- sparse} # with provided sparsity
    else{# cv to choose sparsity
      sumabs_use <- PMD.cv(Omega, type = "standard", sumabss = sparse,
                           center = FALSE, trace = FALSE)$bestsumabs
    }
    pmd <- PMD(Omega, type = "standard", sumabs = sumabs_use, K = K,
               center = FALSE, trace = FALSE)
    list(u = pmd$u[, 1:K, drop = FALSE], v = pmd$v[, 1:K, drop = FALSE])
  }

  ## reference (nuisance) labels for the target samples, no contrast involved
  if(reference == "target"){
    ref <- lead_pair(cov(Xt, Yt), sparse_trt, num_comps)
  }
  else{ # project the target samples on the ancillary-only directions
    ref <- lead_pair(cov(Xa, Ya), sparse_ctrl, num_comps)
  }
  ref_score <- cbind(Xt %*% ref$u, Yt %*% ref$v)
  km <- kmeans(ref_score, centers = km_cluster, nstart = nstart)
  nuisance_labels <- km$cluster

  ## 3CA/s3CA for a range of eta
  nuisance_ARI <- c()
  silhouette <- c()
  for (eta in eta_range){
    ctst <- lead_pair(cov(Xt, Yt) - eta * cov(Xa, Ya), sparse_ctst, num_comps)
    # the contrastive score pair: metabolome-side score against microbiome-side score
    cScore <- cbind(Xt %*% ctst$u, Yt %*% ctst$v)
    ckm <- kmeans(cScore, centers = ckm_cluster, nstart = nstart)
    interesting_labels <- ckm$cluster
    silhouette_score <- silhouette(interesting_labels, dist(cScore))
    silhouette <- c(silhouette, mean(silhouette_score[, "sil_width"]))
    nuisance_ARI <- c(nuisance_ARI, ARI(nuisance_labels, interesting_labels))
  }
  eta_opt <- eta_range[which.max(silhouette - nuisance_ARI)]
  p <- NULL
  if(plot){
    plot_df <- data.frame(eta = rep(eta_range, 2), metrics = c(nuisance_ARI, silhouette),
                          type = rep(c("ARI", "Silhouette"), each = length(eta_range)))
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
  return(list(eta_opt = eta_opt, eta_range = eta_range, nuisance_ARI = nuisance_ARI,
              silhouette = silhouette, nuisance_labels = nuisance_labels, plot = p))
}
