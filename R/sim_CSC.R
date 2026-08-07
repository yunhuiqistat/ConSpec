#' Simulation function for contrastive spectral clustering (CSC)
#'
#' This function performs simulation to prove the interesting of hidden clusters by using cPCA and eta tuning algorithm.
#' Same as in the example of Section 2, we consider two target clusters and two nuisance clusters
#'
#' @import PMA
#' @import MASS
#' @importFrom stats cov
#' @importFrom stats dist 
#' @importFrom stats kmeans 
#' @importFrom stats runif
#' @importFrom stats model.matrix 
#' @importFrom stats rnorm 
#' @importFrom ggpubr ggarrange
#' @param n11 number of target samples in interesting cluster 1 and nuisance cluster 1.
#' @param n12 number of target samples in interesting cluster 1 and nuisance cluster 2.
#' @param n21 number of target samples in interesting cluster 2 and nuisance cluster 1.
#' @param n22 number of target samples in interesting cluster 2 and nuisance cluster 2.
#' @param m1tilde number of ancillary samples in nuisance cluster 1.
#' @param m2tilde number of ancillary samples in nuisance cluster 2.
#' @param mu1 non-zero values in mu_t1.
#' @param mu2 non-zero values in mu_t2.
#' @param theta1 non-zero values in theta_1.
#' @param theta2 non-zero values in theta_2
#' @param p1 number of non-zero elements in mu_t1 and mu_t2.
#' @param p2 number of non-zero elements in theta_1 and theta_2.
#' @param sigma_t standard deviation of the uncorrelated errors in target group.
#' @param sigma_a standard deviation of the uncorrelated errors in ancillary group.
#' @param Sigma_t the error covariance matrix in target group, if is NULL, use standard deviation specified and assume uncorrelation.
#' @param Sigma_a the error covariance matrix in ancillary group, if is NULL, use standard deviation specified and assume uncorrelation.
#' @param eta if eta is a vector, eta_tuning algorithm will be called to select optimal eta; if a number, eta will be used.
#' @param repetition number of replications of the simulation.
#' @param save_plot logical, whether save the score plot and loading plot for the first replication.
#' @param clust_comps number of leading component scores used for kmeans clustering, both inside
#'   \code{\link{eta_tuning_general}} and for the reported ARIs. \code{"auto"} applies the
#'   \code{ckm_cluster - 1} rule, which is 1 for the two interesting clusters in this design.
#' @param nstart number of random starts passed to \code{\link[stats]{kmeans}}.
#' @export
#' @return A data frame containing the selected optimal eta, the four metrics values.
#' @examples
#' # See vignettes/ContrastiveSpectralClustering.Rmd for example usage.
sim_CSC <- function(n11, n12, n21, n22, m1tilde, m2tilde,
                    mu1, mu2, theta1, theta2, p1, p2,
                    sigma_t = 1, sigma_a = 1,
                    Sigma_t = NULL, Sigma_a = NULL,
                    eta, repetition = 10, save_plot = FALSE, clust_comps = "auto", nstart = 25){

  ## ckm_cluster is 2 in this design, so the "auto" rule gives 1. A numeric value is used as is.
  if(identical(clust_comps, "auto")){clust_comps <- 1}

  if(length(eta) > 1){
    eta_range <- eta
  }
  else{eta_use <- eta}

  ARI_trt_nuisance <- c()
  ARI_trt_interesting <- c()
  ARI_ctst_nuisance <- c()
  ARI_ctst_interesting <- c()
  eta_opt_set <- c()

  p <- p1+p2

  mu_t1 <- matrix(rep(c(mu1,0), c(p1,p2)), ncol=1)
  mu_t2 <- matrix(rep(c(mu2,0), c(p1,p2)), ncol=1)
  theta_1 <- matrix(rep(c(0,theta1), c(p1,p2)), ncol=1)
  theta_2 <- matrix(rep(c(0,theta2), c(p1,p2)), ncol=1)

  if(is.null(Sigma_t)){
    Sigma_t <- diag(sigma_t, p)
  }
  if(is.null(Sigma_a)){
    Sigma_a <- diag(sigma_a, p)
  }

  covariates_t <- data.frame(dataset = rep("target", n11+n12+n21+n22),
                             interesting = rep(c("int=1","int=2"), c(n11+n12, n21+n22)),
                             nuisance = rep(c("nui=1","nui=2","nui=1","nui=2"), c(n11, n12, n21, n22)))
  covariates_a <- data.frame(dataset = rep("ancillary", m1tilde+m2tilde),
                             nuisance = rep(c("nui=1","nui=2"), c(m1tilde, m2tilde)))
  for(repe in 1:repetition){
    # generate data
    Xt <- rbind(matrix(mu_t1, nrow = n11, ncol = p, byrow = TRUE) + matrix(theta_1, nrow = n11, ncol = p, byrow = TRUE) + mvrnorm(n = n11, mu = rep(0, p), Sigma = Sigma_t),
                matrix(mu_t1, nrow = n12, ncol = p, byrow = TRUE) + matrix(theta_2, nrow = n12, ncol = p, byrow = TRUE) + mvrnorm(n = n12, mu = rep(0, p), Sigma = Sigma_t),
                matrix(mu_t2, nrow = n21, ncol = p, byrow = TRUE) + matrix(theta_1, nrow = n21, ncol = p, byrow = TRUE) + mvrnorm(n = n21, mu = rep(0, p), Sigma = Sigma_t),
                matrix(mu_t2, nrow = n22, ncol = p, byrow = TRUE) + matrix(theta_2, nrow = n22, ncol = p, byrow = TRUE) + mvrnorm(n = n22, mu = rep(0, p), Sigma = Sigma_t))

    Xa <- rbind(matrix(theta_1, nrow = m1tilde, ncol = p, byrow = TRUE) + mvrnorm(n = m1tilde, mu = rep(0, p), Sigma = Sigma_a),
                matrix(theta_2, nrow = m2tilde, ncol = p, byrow = TRUE) + mvrnorm(n = m2tilde, mu = rep(0, p), Sigma = Sigma_a))


    ## pre-processing
    Xt <- scale(Xt, center = TRUE, scale = FALSE)
    Xa <- scale(Xa, center = TRUE, scale = FALSE)

    if(length(eta) > 1){
      eta_use <- eta_tuning_general(Xt = Xt, Xa = Xa, eta_range = eta_range,
                                    plot = FALSE, sparse_trt = NULL, sparse_ctst = NULL,
                                    num_comps = 20,
                                    km_cluster = 2, ckm_cluster = 2,
                                    clust_comps = clust_comps, nstart = nstart)$eta_opt
      eta_opt_set <- c(eta_opt_set, eta_use)
      cat(paste("Use optimal eta:", eta_use, "\n"))
    }



    # cPCA
    Ut <- eigen(cov(Xt))$vectors
    cU <- eigen(cov(Xt)-eta_use*cov(Xa))$vectors

    PCt <- Xt%*%Ut[,1:2]
    cPC <- Xt%*%cU[,1:2]

    # kmeans clustering using PCs
    PCt_use <- PCt[, 1:clust_comps, drop = FALSE]
    cPC_use <- cPC[, 1:clust_comps, drop = FALSE]
    ARI_trt_nuisance <- c(ARI_trt_nuisance, spec_clust(PCt_use, group = covariates_t$nuisance)$ari)
    ARI_trt_interesting <- c(ARI_trt_interesting, spec_clust(PCt_use, group = covariates_t$interesting)$ari)
    ARI_ctst_nuisance <- c(ARI_ctst_nuisance, spec_clust(cPC_use, group = covariates_t$nuisance)$ari)
    ARI_ctst_interesting <- c(ARI_ctst_interesting, spec_clust(cPC_use, group = covariates_t$interesting)$ari)

    if(repe == 1 & save_plot){

      plot_data <- data.frame(cPC1 = cPC[,1], cPC2 = cPC[,2],
                              PC1 = PCt[,1], PC2 = PCt[,2],
                              covariates_t)
      cP <- ggplot(plot_data)+
        geom_point(aes(x = cPC1, y = cPC2, color = interesting, shape = nuisance))+
        labs(title = "Contrastive")+
        theme_bw()+
        theme(panel.grid = element_blank(),
              text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      Pt <- ggplot(plot_data)+
        geom_point(aes(x=PC1, y=PC2, color= interesting, shape = nuisance), show.legend = FALSE)+
        labs(title = "Treatment")+
        theme_bw()+
        theme(panel.grid = element_blank(),
              text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      PUt1 <- ggplot(data.frame(index = 1:nrow(Ut), Ut1 = Ut[,1]))+
        geom_point(aes(x=index, y=Ut1))+
        ylim(-0.8,0.8)+
        labs(title="Treatment PCA eigenvector 1", x = "Variable index")+
        theme_bw()+
        theme(text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      PUt2 <- ggplot(data.frame(index = 1:nrow(Ut), Ut2 = Ut[,2]))+
        geom_point(aes(x=index, y=Ut2))+
        ylim(-0.8,0.8)+
        labs(title = "Treatment PCA eigenvector 2", x = "Variable index")+
        theme_bw()+
        theme(text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      PcU1 <- ggplot(data.frame(index = 1:nrow(cU), cU1 = cU[,1]))+
        geom_point(aes(x=index, y=cU1))+
        ylim(-0.8,0.8)+
        labs(title = "Contrastive PCA eigenvector 1", x = "Variable index")+
        theme_bw()+
        theme(text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      PcU2 <- ggplot(data.frame(index = 1:nrow(cU), cU2 = cU[,2]))+
        geom_point(aes(x=index, y=cU2))+
        ylim(-0.8,0.8)+
        labs(title = "Contrastive PCA eigenvector 2", x = "Variable index")+
        theme_bw()+
        theme(text = element_text(size = 14),  # Adjust text size
              axis.title = element_text(size = 16),  # Adjust axis title size
              axis.text = element_text(size = 12))
      p_all <- ggarrange(Pt,cP, PUt1, PcU1, PUt2, PcU2, nrow=3, ncol=2,
                         common.legend = TRUE, legend = "bottom", align = "v")
      ggsave(paste("../inst/SavedFigures/", "mu2", mu2,"_p1", p1,"_etaopt.pdf", sep = ""), plot = p_all, width = 8, height = 10, units = "in")
    }
  }
  return(data.frame(ARI_ctst_interesting = ARI_ctst_interesting,
                    ARI_ctst_nuisance = ARI_ctst_nuisance,
                    ARI_trt_interesting = ARI_trt_interesting,
                    ARI_trt_nuisance = ARI_trt_nuisance,
                    eta = if(length(eta_opt_set)) eta_opt_set else rep(eta, repetition)))
}
