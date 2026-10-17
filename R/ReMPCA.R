#' Regularized (Smooth and Sparse) Principal Component Analysis for Hybrid Data
#'
#' Fits regularized principal components to a hybrid data object (functional,
#' image, and non-functional blocks observed on a common set of rows) by a
#' regularized power algorithm. Smoothness and sparsity can be imposed on the
#' row direction (scores, \eqn{u}) and, block by block, on the column
#' direction (loadings, \eqn{v}). Smoothing parameters are selected by
#' (conditional) generalized cross-validation and sparsity parameters by
#' \eqn{K}-fold cross-validation whenever more than one candidate value is
#' supplied.
#'
#' @param hd An hdClass object.
#' @param centerhds A logical; if True, it demeans the data before calculating the principal components.
#' @param num_pcs An integer. The number of principal components.
#' @param sparse_tuning_type A character string specifying the sparse calculation method. Must be one of "soft" (default), "hard", or "SCAD".
#' @param smoothness_type A character string specifying the method used in smoothing u and/or v,
#' must be one of "Second_order" (default), "First_order" or "Indicator".
#' @param nfolds_u  An integer. It's used in cross validation approach for tuning the level of sparsity of rows.
#' @param nfolds_v An integer vector of length \code{p} (the number of variables), where
#' the \code{i}-th element specifies the number of cross-validation folds to use for the
#' \code{i}-th variable. If not provided (\code{NULL}), a value of 5 will be assigned to
#' all variables by default.
#' @param tuning_order A character string representing the tuning order.
#' If set to 'Sparsity', sparsity parameters are tuned first, followed by smoothness.
#'  If set to 'Smoothness', the order is reversed.
#' @param cv.pick A character string specifying the rule used to select the optimal
#' tuning parameter during cross-validation for sparsity. If set to `'min'`, the tuning
#' parameter corresponding to the minimum cross-validation error is chosen. If set to `'1se'`
#' (the default), the 1-standard-error rule is applied, selecting the most regularized
#' model whose error is within one standard error of the minimum.
#' @param parallel Logical; if \code{TRUE}, the cross-validation fits used to tune the
#' sparsity parameters are distributed over a local PSOCK cluster created internally
#' (at most 6 workers, and never more than \code{parallel::detectCores()}). The cluster
#' is stopped before the function returns. Default is \code{FALSE}.
#' @param tuning_iter Integer specifying the number of iterations to perform during the tuning process for conditional smoothing and sparsity parameters.
#' @param sparse_tuning_u Optional. Specifies the sparsity level(s) for rows:
#' \itemize{
#'   \item A single non-negative integer for fixed sparsity.
#'   \item A numeric vector of candidate values, to be selected via cross-validation (CV).
#'   \item Set to \code{0} for no sparsity.
#'   \item If \code{NULL}, the function defaults to the \code{Sparsity_parameter} attribute from the \code{hdClass} object.
#' }
#' @param sparse_tuning_v Optional. A list of length \code{p} (number of variables), where the \code{i}-th element is either a numeric value or a vector specifying candidate sparsity levels for the \code{i}-th variable.
#' If set to \code{NULL}, the function will default to using the \code{Sparsity_parameter_col} attribute from the input object of class \code{hdClass}.
#' @param smooth_tuning_u Optional. Specifies smoothing parameter for rows:
#' \itemize{
#'   \item A single number for fixed smoothness.
#'   \item A numeric vector of candidate values, to be selected via generalized cross-validation (GCV).
#'   \item Set to \code{0} for no smoothness
#'   \item If \code{NULL}, the function defaults to the \code{Smoothing_parameter} attribute from the \code{hdClass} object.
#' }
#' @param smooth_tuning_v Optional. A list of length \code{p} (number of variables), where the \code{i}-th element is either a numeric value or a vector specifying candidate smoothing parameters for the \code{i}-th variable.
#' If set to \code{NULL}, the function will default to using the \code{Smoothing_parameter_col} attribute from the input object of class \code{hdClass}.
#' @param weights Optional numeric vector of scaling weights.
#' If `NULL`, the function automatically computes weights based on the inverse square root
#' of the average variance of each variable. If set to `0`, no scaling is applied.
#' If a numeric vector is provided, its length must match the number of variables, and
#' each element is used to scale the corresponding variable.
#' This is used to adjust for scale differences across variables in the hybrid data object.
#' @param thresh The convergence threshold in power algorithm.
#' @param maxit Maximum number of iterations in power algorithm.
#' @param verbose Logical; if \code{TRUE} (default), progress bars and progress
#' messages are printed. Set to \code{FALSE} for silent fitting.
#'
#' @details
#' Components are extracted sequentially. After the \eqn{j}-th regularized
#' component \eqn{\hat u_j \hat v_j^\top} is estimated (with \eqn{\|\hat v_j\| = 1}),
#' it is subtracted from the (centered and weighted) data matrix and the next
#' component is estimated from the residual matrix.
#'
#' Candidate sparsity levels are interpreted as the number of entries set to zero
#' and should be supplied in increasing order (least sparse first), because the
#' one-standard-error rule (\code{cv.pick = "1se"}) selects the largest candidate
#' whose CV error is within one standard error of the minimum. The cross-validation
#' folds are formed with a fixed internal seed so that results are reproducible;
#' the global random number generator state of the user is not modified.
#'
#' @references
#' Huang, J. Z., Shen, H. and Buja, A. (2008). Functional principal components
#' analysis via penalized rank one approximation. \emph{Electronic Journal of
#' Statistics}, 2, 678--695. \doi{10.1214/08-EJS218}
#'
#' Shen, H. and Huang, J. Z. (2008). Sparse principal component analysis via
#' regularized low rank matrix approximation. \emph{Journal of Multivariate
#' Analysis}, 99(6), 1015--1034. \doi{10.1016/j.jmva.2007.06.007}
#'
#' Huang, J. Z., Shen, H. and Buja, A. (2009). The analysis of two-way functional
#' data using two-way regularized singular value decompositions. \emph{Journal of
#' the American Statistical Association}, 104(488), 1609--1620.
#' \doi{10.1198/jasa.2009.tm08024}
#'
#' @importFrom utils  txtProgressBar setTxtProgressBar
#' @importFrom Matrix bdiag
#' @importFrom stats var
#'
#' @return A list with components
#' \describe{
#'   \item{ReconstructedData}{List (one element per component) of the fitted
#'     rank-one matrices \eqn{\hat u_j \hat v_j^\top}, on the centered and weighted scale.}
#'   \item{PCFunctions}{List (one element per component) of lists (one element per
#'     variable) holding the estimated loadings \eqn{\hat v_j}.}
#'   \item{PCScores}{Matrix of scores \eqn{\hat u_j}, one column per component.}
#'   \item{OptimalAlphaV, OptimalAlphaU}{Selected smoothing parameters, per component.}
#'   \item{OptimalGammaV, OptimalGammaU}{Selected sparsity parameters, per component.}
#'   \item{GCVResultsV, GCVResultsU}{GCV scores over the candidate smoothing parameters.}
#'   \item{CVResultsV, CVResultsU}{CV errors and standard errors over the candidate
#'     sparsity parameters.}
#'   \item{variable_types}{Type of each variable: \code{"hd"} (functional) or \code{"rd"} (regular).}
#'   \item{VarianceExplained}{Proportion of the total variation of the (centered and
#'     weighted) data explained by each estimated component, computed by projecting
#'     the data onto the span of the estimated loadings.}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 40
#' u <- sin(seq(0, pi, length.out = n))
#' X1 <- outer(u, cos(seq(0, pi, length.out = 15))) + matrix(rnorm(n * 15, sd = 0.2), n)
#' X2 <- outer(u, c(1, -1, 0, 0, 0.5)) + matrix(rnorm(n * 5, sd = 0.2), n)
#' hd <- hdClass(list(fdClass(X1, Smoothing_parameter = c(0.01, 0.1, 1)),
#'                    rdClass(X2, Sparsity_parameter = 0:2)))
#' fit <- ReMPCA(hd, num_pcs = 1, nfolds_u = 2, nfolds_v = c(2, 2),
#'               cv.pick = "min", verbose = FALSE)
#' fit$OptimalAlphaV
#' fit$OptimalGammaV
#' plot_pc_functions(fit)
#' @export
################### Smooth and Sparse Multivariate PCA ###################
ReMPCA <- function(hd,
                   centerhds = TRUE,
                   num_pcs = 1,
                   nfolds_u = 5,
                   nfolds_v = NULL,
                   thresh = 1e-10,
                   maxit = 100,
                   tuning_iter = 1,
                   parallel = FALSE,
                   weights = NULL,
                   smoothness_type = "Second_order",
                   sparse_tuning_type = "soft",
                   tuning_order = "Sparsity",
                   cv.pick = "1se",
                   sparse_tuning_u = NULL,
                   sparse_tuning_v = NULL,
                   smooth_tuning_u = NULL,
                   smooth_tuning_v = NULL,
                   verbose = TRUE) {

  # Check if hd is a hd object
  if (!inherits(hd, "hdClass")) {
    stop("hd must be of class 'hdClass'!")
  }
  # Validate character options
  smoothness_type <- match.arg(smoothness_type,
                               c("Second_order", "First_order", "Indicator"))
  sparse_tuning_type <- match.arg(sparse_tuning_type, c("soft", "hard", "SCAD"))
  tuning_order <- match.arg(tuning_order, c("Sparsity", "Smoothness"))
  cv.pick <- match.arg(cv.pick, c("1se", "min"))
  if (!is.numeric(num_pcs) || length(num_pcs) != 1 || num_pcs < 1) {
    stop("'num_pcs' must be a positive integer.")
  }
  n_var <- attr(hd, "n_var")

  # weights
  if (is.null(weights)) {
    scale_hd_out <- scale_hd(hd)
    weights <- scale_hd_out$weights
    hddata <- scale_hd_out$scaled_hd

  } else if (identical(weights, 0)) {
    hddata <- hd$matrix

  } else if (is.vector(weights) &&
             is.numeric(weights) &&
             length(weights) == n_var) {
    scaled_matrix <- hd$matrix
    start_idx <- 1
    ncol_vec <- as.numeric(data.frame(attr(hd, "ncol")))
    for (i in 1:n_var) {
      end_idx <- as.numeric(start_idx + ncol_vec[i] - 1)
      mat <- scaled_matrix[, start_idx:end_idx, drop = FALSE]
      scaled_matrix[, start_idx:end_idx] <- weights[i] * mat
      start_idx <- end_idx + 1
    }
    hddata <- scaled_matrix

  } else {
    stop("The 'weights' must be NULL, 0, or a numeric vector of length equal to the number of variables.")
  }


  # Combine matrices side by side
  n <- nrow(hd)
  n_var <- attr(hd, "n_var") # Number of variables (# of matrices in hd)
  ncol <- data.frame(attr(hd, "ncol")) # Number of columns of each matrix
  # nfolds_v with no default
  if (is.null(nfolds_v)) {
    nfolds_v <- rep(5, attr(hd, "n_var"))
  }

  ####### Smoothing Parameter ##########
  # Smoothing parameters for column
  # Generate all combinations alphas (one row per combination): A matrix
  if(!is.null(smooth_tuning_v)){
    smooth_tuning_col <- expand.grid(smooth_tuning_v)
    }else{
      smooth_tuning_col <- expand.grid(attr(hd, "Smoothing_parameter_col"))
    }

  # Smoothing parameters for row: A vector
  if(!is.null(smooth_tuning_u)){
    smooth_tuning_row <- smooth_tuning_u
  }else{
    smooth_tuning_row <- attr(hd, "Smoothing_parameter")
  }

  ####### level of sparsity #######
  if(!is.null(sparse_tuning_u)){
    sparsity_row_list <- sparse_tuning_u
  }else{
    sparsity_row_list <- attr(hd, "Sparsity_parameter")
  }
  if(!is.null(sparse_tuning_v)){
    sparsity_col_list <- sparse_tuning_v
  }else{
    sparsity_col_list <- attr(hd, "Sparsity_parameter_col")
  }

  ####### Pre-processing: Centralizing the data #######
  X <- data.frame()
  if (centerhds) {
    X <- apply(hddata, 2, function(x) x - mean(x))
  }else{
    X <- hddata
  }

  ####### Grid Points #######
  GridPoints_u <- attr(hd, "GridPoints_u")
  GridPoints_v <- attr(hd, "GridPoints_v")
  variable_types <- attr(hd, "variable_types")

  ####### S_alpha for all alphas #######
  S_alpha_list_u <- S_alpha_list_v <- Omegas_v <- list()
  index <- 0
  if (verbose) message("Preprocessing ...")
  n_iter1 <- nrow(smooth_tuning_col)     # The number of alphas
  if (verbose) pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter1,# Maximum value of the progress bar
                       style = 3,    # Progress bar style
                       width = 50,   # Progress bar width.
                       char = "=")   # Character used to create the bar

  # S_alpha for v
  for (alpha_index in 1:nrow(smooth_tuning_col)) {
    index <- index + 1
    S_u <- S_v <- Omega_v <- Omegas_u <- list()
    for (i in 1:n_var) {
      alpha <- as.numeric(smooth_tuning_col[alpha_index,i])
      if(is.null(GridPoints_v[[i]])){
        S_v[[i]] <- diag(ncol[,i])
        Omega_v[[i]] <- diag(ncol[,i])
      }else{
        get.pen.result <- get.pen(td = GridPoints_v[[i]],
                                  alpha = alpha,
                                  type = smoothness_type)
        S_v[[i]] <- get.pen.result$S.alpha
        Omega_v[[i]] <- get.pen.result$ Omega
      }
    }
    S_alpha_list_v[[index]] <- S_v #as.matrix(bdiag(S_v))
    Omegas_v[[index]] <- Omega_v
    if (verbose) setTxtProgressBar(pb, index)
  }

  # S_alpha for u
  if(!is.null(smooth_tuning_row)){
    for (alpha_index in 1:length(smooth_tuning_row)) {
      alpha <- as.numeric(smooth_tuning_row[alpha_index])
      get.pen.result <- get.pen(td = GridPoints_u,
                                alpha = alpha,
                                type = smoothness_type)

      S_alpha_list_u[[alpha_index]] <- get.pen.result$S.alpha
      Omegas_u[[alpha_index]] <- get.pen.result$Omega
    }
  }
  if (verbose) close(pb)

  # Initializing the lists
  opt_S_v <- opt_S_u <- list()
  smooth_result_u <- smooth_result_v <- list()
  sparse_result_u <- sparse_result_v <- list()
  GCV_v <- GCV_u <- CV_v <- CV_u <- list()
  lsv <- lsu <- c()
  funcs <- PCs <- list()
  V_hat <- matrix(0, nrow = ncol(X), ncol = num_pcs)  # estimated loadings
  X_orig <- X; X_temp <- X
  pve <- numeric(num_pcs)

  ####### ReMPCA Implementation #######
  for (j in 1:num_pcs) {

    if (verbose) message(sprintf("Computing the %s PC ...", ordinal(j)))
    # For j > 1, X_temp is the residual after removing the previously
    # estimated regularized components (deflation step of the algorithm).

    # Tuning Parameters
    param_result <- list()
    param_result <- parameter_selection(X_temp =  X_temp,
                                        n_var = n_var,
                                        ncol = ncol,
                                        n = n,
                                        GridPoints_u = GridPoints_u,
                                        GridPoints_v = GridPoints_v,
                                        smooth_tuning_v = smooth_tuning_col,
                                        smooth_tuning_u = smooth_tuning_row,
                                        sparse_tuning_u = sparsity_row_list,
                                        sparse_tuning_v = sparsity_col_list,
                                        sparse_tuning_type = sparse_tuning_type,
                                        nfolds_u = nfolds_u,
                                        nfolds_v = nfolds_v,
                                        parallel = parallel,
                                        S_alpha_list_v = S_alpha_list_v,
                                        S_alpha_list_u = S_alpha_list_u,
                                        Omegas_u = Omegas_u,
                                        Omegas_v = Omegas_v,
                                        tuning_iter = tuning_iter,
                                        tuning_order = tuning_order,
                                        thresh = thresh,
                                        maxit = maxit,
                                        cv.pick = cv.pick,
                                        smoothness_type = smoothness_type,
                                        verbose = verbose)

    # Optimal parameters
    sparse_result_u[[j]] <- param_result$sparse_tuning_selection_u
    sparse_result_v[[j]] <- param_result$sparse_tuning_selection_v
    smooth_result_u[[j]] <- param_result$smooth_tuning_selection_u
    smooth_result_v[[j]] <- param_result$smooth_tuning_selection_v

    # Optimal Smoothing Matrices
    opt_S_u[[j]] <- param_result$last_opt_S_u
    opt_S_v[[j]] <- param_result$last_opt_S_v

    # Generalized cross-validation scores
    GCV_v[[j]] <- param_result$last_gcv_result_v
    GCV_u[[j]] <- param_result$last_gcv_result_u

    # Cross-validation scores
    CV_v[[j]] <- param_result$last_cv_result_v
    CV_u[[j]] <- param_result$last_cv_result_u

    # Extracting v and u having the optimal parameters
    test_result <- power_algo(data = X_temp,
                              n_var = n_var,
                              ncol = ncol,
                              conditional = FALSE,
                              thresh = thresh,
                              maxit = maxit,
                              sparse_tuning_result_u = sparse_result_u[[j]],
                              sparse_tuning_result_v = sparse_result_v[[j]],
                              S_alpha_v = opt_S_v[[j]],
                              S_alpha_u = opt_S_u[[j]],
                              alpha_Omega_v=param_result$last_opt_alpha_omega_v,
                              alpha_Omega_u=param_result$last_opt_alpha_omega_u,
                              sparse_tuning_type = sparse_tuning_type)


    v <- test_result[[1]]
    u <- test_result[[2]]
    lsv <- cbind(lsv, v)
    lsu <- cbind(lsu, u)
    funcs[[j]] <- u%*%t(v)
    V_hat[, j] <- v

    # Deflation: remove the estimated regularized component
    X_temp <- X_temp - funcs[[j]]

    # Splitting v for variables
    new_PC <- list()
    for (i in 1:n_var) {
      lsv <- data.frame(lsv)
      rows_to_extract <- 1:as.integer(ncol[i])
      new_PC[[i]] <- lsv[rows_to_extract,]
      lsv <- lsv[-rows_to_extract,]
    }
    PCs[[j]] <- new_PC
  }

  # Variance explained by each PC
  V_list <- lapply(seq_len(num_pcs), function(k) V_hat[, k])
  Variance <- compute_variance_explained(X_orig, V_list = V_list)

  return(list(
    ReconstructedData = funcs,              # Reconstructed hybrid data
    PCFunctions = PCs,                      # Estimated PC (v, per var/PC)
    PCScores = lsu,                         # PC scores u for each component
    OptimalAlphaV = smooth_result_v,        # Alph_v (per var/PC)
    OptimalAlphaU = smooth_result_u,        # Alph_u (per var/PC)
    OptimalGammaV = sparse_result_v,        # gamma_v (per var/PC)
    OptimalGammaU = sparse_result_u,        # gamma_u (per var/PC)
    GCVResultsV = GCV_v,                    # GCV for v (per variable/component)
    GCVResultsU = GCV_u,                    # GCV for u (per component)
    CVResultsV = CV_v,                      # CV for v (per var/PC)
    CVResultsU = CV_u,                      # CV for u (per component)
    variable_types = variable_types,        # Variable types
    VarianceExplained = Variance$AdjPercVar # Variance explained by each PC
  ))
}


################ Percentage of variance explained by each PCA ################
compute_variance_explained <- function(X, V_list) {
  K <- length(V_list)
  n <- nrow(X)
  m <- ncol(X)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  # PC matrix V_k
  Vmat_list <- lapply(1:K, function(k) do.call(cbind, V_list[1:k]))

  # Compute projection matrices H_k
  H_list <- lapply(Vmat_list, function(Vk) {
    solve_term <- solve(t(Vk) %*% Vk)
    Hk <- Vk %*% solve_term %*% t(Vk)
    return(Hk)
  })

  # Project the data matrix onto V
  Xk_list <- lapply(H_list, function(Hk) X %*% Hk)

  # Total variance of the original data
  TotalVar <- sum(X^2)

  # Compute variance explained and adjusted variance
  VarExplained <- sapply(Xk_list, function(Xk) sum(Xk^2))
  AdjVar <- c(VarExplained[1], diff(VarExplained))

  # Convert to percentages
  CPEV <- VarExplained / TotalVar
  AdjPercVar <- AdjVar / TotalVar

  return(list(
    CPEV = CPEV,
    AdjPercVar = AdjPercVar
  ))
}
