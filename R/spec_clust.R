#' Cluster samples using ConSpec component scores
#'
#' Performs k-means clustering on component scores obtained from a ConSpec
#' analysis, such as contrastive PCA or contrastive CCA. If reference labels
#' are supplied, the adjusted Rand index (ARI) is also calculated.
#'
#' @param score_mat A numeric matrix or data frame containing component scores,
#'   with samples in rows and components in columns.
#' @param k Number of clusters. If `NULL` and `group` is supplied, `k` is set
#'   to the number of unique non-missing reference labels.
#' @param group Optional vector of reference group labels. Must have one value
#'   per row of `score_mat`.
#' @param components Integer vector specifying the score dimensions used for
#'   clustering. By default, all columns of `score_mat` are used.
#' @param nstart Number of random initializations passed to
#'   [stats::kmeans()]. Default is 25.
#' @param iter.max Maximum number of k-means iterations.
#' @param seed Optional random seed for reproducibility.
#' @param ... Additional arguments passed to [stats::kmeans()].
#'
#' @return A list containing:
#' \describe{
#'   \item{cluster}{Integer vector of cluster assignments.}
#'   \item{ari}{Adjusted Rand index, or `NULL` if `group` is not supplied.}
#'   \item{kmeans}{The fitted object returned by [stats::kmeans()].}
#'   \item{scores}{The score matrix used for clustering.}
#'   \item{components}{The component indices used for clustering.}
#' }
#'
#' @importFrom stats kmeans
#' @importFrom aricode ARI
#' @export
#'
#' @examples
#' set.seed(1)
#' scores <- matrix(rnorm(100), nrow = 20)
#'
#' # Unsupervised clustering
#' fit <- spec_clust(scores, k = 2)
#' fit$cluster
#'
#' # Clustering evaluation using known labels
#' group <- rep(c("A", "B"), each = 10)
#' fit <- spec_clust(scores, k = 2, group = group)
#' fit$ari
spec_clust <- function(
    score_mat,
    k = NULL,
    group = NULL,
    components = NULL,
    nstart = 25,
    iter.max = 100,
    seed = NULL,
    ...
) {
    score_mat <- as.matrix(score_mat)

    if (!is.numeric(score_mat)) {
        stop("`score_mat` must be a numeric matrix or data frame.",
             call. = FALSE)
    }

    if (nrow(score_mat) < 2L) {
        stop("`score_mat` must contain at least two samples.",
             call. = FALSE)
    }

    if (ncol(score_mat) < 1L) {
        stop("`score_mat` must contain at least one score dimension.",
             call. = FALSE)
    }

    if (anyNA(score_mat) || any(!is.finite(score_mat))) {
        stop("`score_mat` cannot contain missing or non-finite values.",
             call. = FALSE)
    }

    if (is.null(components)) {
        components <- seq_len(ncol(score_mat))
    }

    if (
        !is.numeric(components) ||
        anyNA(components) ||
        any(components < 1L) ||
        any(components > ncol(score_mat))
    ) {
        stop(
            "`components` must contain valid column indices of `score_mat`.",
            call. = FALSE
        )
    }

    components <- unique(as.integer(components))
    scores_used <- score_mat[, components, drop = FALSE]

    if (!is.null(group)) {
        if (length(group) != nrow(score_mat)) {
            stop(
                "`group` must have one label for each row of `score_mat`.",
                call. = FALSE
            )
        }

        if (anyNA(group)) {
            stop("`group` cannot contain missing values.",
                 call. = FALSE)
        }
    }

    if (is.null(k)) {
        if (is.null(group)) {
            stop(
                "Supply `k`, or supply `group` so that `k` can be inferred.",
                call. = FALSE
            )
        }

        k <- length(unique(group))
    }

    if (
        length(k) != 1L ||
        is.na(k) ||
        k != as.integer(k) ||
        k < 2L ||
        k >= nrow(score_mat)
    ) {
        stop(
            "`k` must be a single integer between 2 and nrow(score_mat) - 1.",
            call. = FALSE
        )
    }

    k <- as.integer(k)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    km_fit <- stats::kmeans(
        x = scores_used,
        centers = k,
        nstart = nstart,
        iter.max = iter.max,
        ...
    )

    ari <- NULL

    if (!is.null(group)) {
        ari <- unname(aricode::ARI(km_fit$cluster, group))
    }

    structure(
        list(
            cluster = km_fit$cluster,
            ari = ari,
            kmeans = km_fit,
            scores = scores_used,
            components = components
        ),
        class = "conspec_clust"
    )
}