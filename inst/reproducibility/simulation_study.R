###############################################################################
# Simulation study, benchmark comparison and convergence diagnostics
#
# Pourmoshir, M. and Maadooliat, M. "Regularized Principal Component Analysis
# for Multivariate Two-way Hybrid Data", Journal of Data Science (revision 2).
#
# Reproduces
#   * Table 2 / Figure 2 : relative reconstruction error of the six ReMPCA
#                          configurations over the 3 x 3 noise-by-sparsity
#                          scenarios (100 runs each),
#   * Figure 3           : one representative run (medium noise, medium sparsity),
#   * Table 3            : benchmark of ReHPCA against competing hybrid and
#                          regularized PCA methods on the same data sets,
#   * Table S1 / Fig. S1 : convergence and initialization diagnostics of
#                          Algorithm 1 (SVD start versus random starts).
#
# Usage
#   Sys.setenv(REMPCA_OUTPUT_DIR = "ReMPCA_output",   # output folder
#              REMPCA_CORES      = "8",               # forked workers (1 on Windows)
#              REMPCA_NREP       = "100")             # runs per scenario
#   source(system.file("reproducibility", "simulation_study.R", package = "ReMPCA"))
#
# Every data set is generated from its own seed (scenario id and run number),
# so the results do not depend on the number of workers.
#
# Required packages: ReMPCA, MFPCA (>= 1.3), funData, PMA, fda, ggplot2,
# patchwork, Matrix, parallel (base R).
###############################################################################

suppressPackageStartupMessages({
  library(ReMPCA)
  library(funData)
  library(MFPCA)
  library(PMA)
  library(fda)
  library(ggplot2)
  library(patchwork)
  library(parallel)
})

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

out_dir <- Sys.getenv("REMPCA_OUTPUT_DIR", "ReMPCA_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
n_rep   <- as.integer(Sys.getenv("REMPCA_NREP", "100"))
n_cores <- as.integer(Sys.getenv("REMPCA_CORES",
                                 max(1L, parallel::detectCores() - 1L)))
if (.Platform$OS.type == "windows") n_cores <- 1L
n_conv  <- as.integer(Sys.getenv("REMPCA_NCONV", "20"))   # runs per scenario
n_start <- as.integer(Sys.getenv("REMPCA_NSTART", "20"))  # random starts per run
tuning_order <- Sys.getenv("REMPCA_TUNING_ORDER", "Smoothness")
cv_pick      <- Sys.getenv("REMPCA_CV_PICK", "min")
scen_subset  <- Sys.getenv("REMPCA_SCENARIOS", "")  # e.g. "5" or "1,5,9"; "" = all
stage        <- Sys.getenv("REMPCA_STAGE", "all")   # "all", "fit", "conv", "traces", "figs"
# runs whose iteration traces are drawn in Figure S1 ("scenario:run"): a
# medium-noise run and two high-noise runs (slow convergence / oscillation)
trace_runs   <- strsplit(Sys.getenv("REMPCA_TRACE_RUNS", "5:1,7:9,7:18"), ",")[[1]]

###############################################################################
# 1. Simulation design (Section 3)
###############################################################################
# Rank-one two-way hybrid signal X = d * u v^T with three column blocks:
#   X1: functional (m1 = 40 grid points), X2: non-functional vector (m2 = 20),
#   X3: functional (m3 = 30). Rows (n = 150) are ordered on [0, 1].
# For sparsity level c (fraction of inactive entries):
#   u(s)  = 0 for s <= c and -sin{pi (s - c) / (1 - c)} for s > c,
#   v1(t) = -0.5 cos{2 pi t / (1 - c)} for t < 1 - c and 0 otherwise,
#   v2    = N(0, 1) entries with round(c * m2) of them set to zero,
#   v3(t) = cos(pi t + pi / 2).
# With c = 1/2 these reduce to the formulas displayed in Section 3. u and
# v = (v1, v2, v3) are scaled to unit length and d = 4.4, which places the
# three noise levels sigma = 0.12, 0.40, 0.85 below, near and above the
# detection threshold sigma * (n p)^(1/4) of the unpenalized rank-one SVD.
n <- 150
m <- c(X1 = 40, X2 = 20, X3 = 30)
p <- sum(m)
idx <- split(seq_len(p), rep(names(m), m))[names(m)]
d_signal <- 4.4
noise_levels    <- c(Low = 0.12, Medium = 0.40, High = 0.85)
sparsity_levels <- c(Low = 0.10, Medium = 0.45, High = 0.72)
scenarios <- expand.grid(sparsity = names(sparsity_levels),
                         noise = names(noise_levels),
                         stringsAsFactors = FALSE)[, c("noise", "sparsity")]
scenarios$id <- seq_len(nrow(scenarios))
scenarios$sigma <- unname(noise_levels[scenarios$noise])
scenarios$c_sp  <- unname(sparsity_levels[scenarios$sparsity])
run_ids <- if (nzchar(scen_subset)) as.integer(strsplit(scen_subset, ",")[[1]]) else scenarios$id

grid_s  <- seq(1 / n, 1, length.out = n)          # same grids as hdClass()/fdClass()
grid_t1 <- seq(1 / m[["X1"]], 1, length.out = m[["X1"]])
grid_t3 <- seq(1 / m[["X3"]], 1, length.out = m[["X3"]])

true_components <- function(c_sp) {
  a  <- 1 - c_sp
  u  <- ifelse(grid_s <= c_sp, 0, -sin(pi * (grid_s - c_sp) / a))
  v1 <- ifelse(grid_t1 < a, -0.5 * cos(2 * pi * grid_t1 / a), 0)
  v2 <- rnorm(m[["X2"]])
  v2[sample.int(m[["X2"]], round(c_sp * m[["X2"]]))] <- 0
  v3 <- cos(pi * grid_t3 + pi / 2)
  v  <- c(v1, v2, v3)
  list(u = u / sqrt(sum(u^2)), v = v / sqrt(sum(v^2)))
}

simulate_run <- function(scn, rep) {
  seed <- 100000L * scn$id + rep
  set.seed(seed)
  tc <- true_components(scn$c_sp)
  X  <- d_signal * tcrossprod(tc$u, tc$v)
  Y  <- X + matrix(rnorm(n * p, sd = scn$sigma), n, p)
  list(X = X, Y = Y, u = tc$u, v = tc$v, seed = seed)
}

###############################################################################
# 2. Evaluation
###############################################################################
sin_angle <- function(a, b) {
  cs <- abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2))
  sqrt(max(0, 1 - cs^2))
}
# Relative reconstruction error ||Xhat - X||_F / ||X||_F, the sine of the
# angle between the true u (v) and the leading left (right) singular vector of
# the fitted matrix (for rank-one fits these are the fitted u and v), and the
# F1 score for recovering the zero entries of u and v.
evaluate <- function(Xhat, sim) {
  sv <- svd(Xhat, nu = 1, nv = 1)
  vh <- sv$v[, 1]; uh <- sv$u[, 1]
  f1 <- function(truth, est) {
    tp <- sum(truth & est); if (tp == 0) return(0)
    2 * tp / (sum(truth) + sum(est))
  }
  c(rel_err = norm(Xhat - sim$X, "F") / norm(sim$X, "F"),
    sin_u = sin_angle(uh, sim$u),
    sin_v = sin_angle(vh, sim$v),
    f1_zero_u = f1(sim$u == 0, abs(uh) < 1e-10),
    f1_zero_v = f1(sim$v == 0, abs(vh) < 1e-10))
}
na_row <- function(time) c(rel_err = NA, sin_u = NA, sin_v = NA, f1_zero_u = NA,
                           f1_zero_v = NA, time = time)

###############################################################################
# 3. ReMPCA configurations (Table 2 / Figure 2)
###############################################################################
alpha_grid   <- 2^seq(-20, 0, by = 2)               # smoothing candidates (GCV)
gamma_u_grid <- seq(0, 135, by = 15)                # number of zeroed rows (CV)
gamma_v_grid <- list(seq(0, 36, by = 4),            # zeroed entries, X1
                     seq(0, 18, by = 2),            # zeroed entries, X2
                     seq(0, 27, by = 3))            # zeroed entries, X3
configs <- list(
  "No penalty"            = c(su = FALSE, pu = FALSE, sv = FALSE, pv = FALSE),
  "Smooth+Sparse u"       = c(su = TRUE,  pu = TRUE,  sv = FALSE, pv = FALSE),
  "Smooth+Sparse v"       = c(su = FALSE, pu = FALSE, sv = TRUE,  pv = TRUE),
  "Two-way Smooth"        = c(su = TRUE,  pu = FALSE, sv = TRUE,  pv = FALSE),
  "Two-way Sparse"        = c(su = FALSE, pu = TRUE,  sv = FALSE, pv = TRUE),
  "Two-way Smooth+Sparse" = c(su = TRUE,  pu = TRUE,  sv = TRUE,  pv = TRUE))

make_hd <- function(Y) {
  hdClass(list(fdClass(Y[, idx$X1, drop = FALSE]),
               rdClass(Y[, idx$X2, drop = FALSE]),
               fdClass(Y[, idx$X3, drop = FALSE])))
}

fit_rempca <- function(Y, cfg) {
  fit <- ReMPCA(make_hd(Y), centerhds = FALSE, num_pcs = 1, weights = 0,
                nfolds_u = 5, nfolds_v = c(5, 5, 5), thresh = 1e-10, maxit = 100,
                smoothness_type = "Second_order", sparse_tuning_type = "soft",
                tuning_order = tuning_order, cv.pick = cv_pick,
                sparse_tuning_u = if (cfg[["pu"]]) gamma_u_grid else 0,
                sparse_tuning_v = if (cfg[["pv"]]) gamma_v_grid else list(0, 0, 0),
                smooth_tuning_u = if (cfg[["su"]]) alpha_grid else 0,
                smooth_tuning_v = if (cfg[["sv"]]) list(alpha_grid, 0, alpha_grid)
                                  else list(0, 0, 0),
                verbose = FALSE)
  list(Xhat = as.matrix(fit$ReconstructedData[[1]]),
       u = as.numeric(fit$PCScores[, 1]),
       v = unlist(lapply(fit$PCFunctions[[1]], as.numeric)),
       alpha_u = as.numeric(unlist(fit$OptimalAlphaU[[1]])),
       alpha_v = as.numeric(unlist(fit$OptimalAlphaV[[1]])),
       gamma_u = as.numeric(unlist(fit$OptimalGammaU[[1]])),
       gamma_v = as.numeric(unlist(fit$OptimalGammaV[[1]])))
}

###############################################################################
# 4. Competing methods (Table 3)
###############################################################################
# (a) Pre-smoothing hybrid PCA (Ramsay and Silverman, 2005): each functional
#     curve is smoothed by a cubic smoothing spline with a common, GCV-selected
#     penalty per block; the smoothed curves are concatenated with the vector
#     block and an ordinary (centered) PCA is applied.
presmooth_block <- function(Xb, tgrid, loglam = seq(-12, 0, by = 0.5)) {
  basis <- fda::create.bspline.basis(range(tgrid), norder = 4, breaks = tgrid)
  gcv <- vapply(loglam, function(ll) {
    sum(fda::smooth.basis(tgrid, t(Xb), fda::fdPar(basis, 2, 10^ll))$gcv)
  }, numeric(1))
  sm <- fda::smooth.basis(tgrid, t(Xb), fda::fdPar(basis, 2, 10^loglam[which.min(gcv)]))
  t(fda::eval.fd(tgrid, sm$fd))
}
fit_presmooth_hpca <- function(Y) {
  Z <- cbind(presmooth_block(Y[, idx$X1], grid_t1), Y[, idx$X2],
             presmooth_block(Y[, idx$X3], grid_t3))
  mu <- colMeans(Z)
  sv <- svd(sweep(Z, 2, mu), nu = 1, nv = 1)
  list(Xhat = matrix(mu, n, p, byrow = TRUE) + sv$d[1] * tcrossprod(sv$u, sv$v))
}

# (b) HFV-PCA (Jang, 2021): (1) MFPCA of the functional blocks (Happ and
#     Greven, 2018; univariate PACE expansions explaining 99% of the variance),
#     (2) PCA of the vector block, (3) concatenation of both score sets with the
#     weight omega = sum_i ||Y_i - Ybar||^2_L2 / sum_i ||x_i - xbar||^2 that
#     balances functional and vector variation, and (4) eigen-analysis of the
#     covariance of the concatenated scores. The rank-one fit is
#     mean + first hybrid component.
fit_hfvpca <- function(Y, pve = 0.99) {
  Y1 <- Y[, idx$X1]; Y2 <- Y[, idx$X2]; Y3 <- Y[, idx$X3]
  f1 <- funData(grid_t1, Y1); f3 <- funData(grid_t3, Y3)
  npc <- c(PACE(f1, pve = pve)$npc, PACE(f3, pve = pve)$npc)
  mf <- MFPCA(multiFunData(f1, f3), M = sum(npc),
              uniExpansions = list(list(type = "uFPCA", npc = npc[1]),
                                   list(type = "uFPCA", npc = npc[2])),
              fit = FALSE)
  L <- which(cumsum(mf$values) / sum(mf$values) > pve)[1]
  if (is.na(L)) L <- length(mf$values)
  rho  <- mf$scores[, seq_len(L), drop = FALSE]
  psi1 <- mf$functions[[1]]@X[seq_len(L), , drop = FALSE]   # L x m1
  psi3 <- mf$functions[[2]]@X[seq_len(L), , drop = FALSE]   # L x m3
  mu1 <- as.numeric(mf$meanFunction[[1]]@X); mu3 <- as.numeric(mf$meanFunction[[2]]@X)
  Y1c <- sweep(Y1, 2, colMeans(Y1)); Y3c <- sweep(Y3, 2, colMeans(Y3))
  mu2 <- colMeans(Y2); Y2c <- sweep(Y2, 2, mu2)
  omega <- (sum(Y1c^2) / m[["X1"]] + sum(Y3c^2) / m[["X3"]]) / sum(Y2c^2)
  sv2 <- svd(sqrt(omega) * Y2c)
  J <- which(cumsum(sv2$d^2) / sum(sv2$d^2) > pve)[1]
  zeta <- sv2$u[, seq_len(J), drop = FALSE] %*% diag(sv2$d[seq_len(J)], J)
  B  <- cbind(rho, zeta)
  e1 <- eigen(stats::cov(B), symmetric = TRUE)$vectors[, 1]
  s1 <- as.numeric(B %*% e1)                       # first hybrid PC scores
  c1 <- e1[seq_len(L)]; d1 <- e1[L + seq_len(J)]
  Xhat <- cbind(matrix(mu1, n, m[["X1"]], byrow = TRUE) + outer(s1, as.numeric(c1 %*% psi1)),
                matrix(mu2, n, m[["X2"]], byrow = TRUE) +
                  outer(s1, as.numeric(sv2$v[, seq_len(J), drop = FALSE] %*% d1) / sqrt(omega)),
                matrix(mu3, n, m[["X3"]], byrow = TRUE) + outer(s1, as.numeric(c1 %*% psi3)))
  list(Xhat = Xhat, L = L, J = J, omega = omega)
}

# (c) Penalized matrix decomposition (Witten, Tibshirani and Hastie, 2009):
#     rank-one L1-penalized SVD with sparsity in both directions; the bound is
#     selected by the package's 5-fold cross-validation.
fit_pmd <- function(Y) {
  # sumabs * sqrt(p) must be at least 1, hence the lower end of the grid
  cv <- PMD.cv(Y, type = "standard", sumabss = seq(0.12, 0.96, by = 0.06),
               center = FALSE, trace = FALSE)
  f <- PMD(Y, type = "standard", sumabs = cv$bestsumabs, K = 1,
           center = FALSE, trace = FALSE)
  list(Xhat = f$d * tcrossprod(f$u, f$v), sumabs = cv$bestsumabs)
}

# (d) Two-way regularized SVD (Huang, Shen and Buja, 2009) applied to the
#     concatenated matrix as if it were a single two-way functional data set:
#     one roughness penalty over all p columns (ignoring block boundaries and
#     the non-functional block) and one over the rows, both selected by GCV.
fit_twrsvd <- function(Y) {
  fit <- ReMPCA(hdClass(list(fdClass(Y))), centerhds = FALSE, num_pcs = 1,
                weights = 0, thresh = 1e-10, maxit = 100,
                smoothness_type = "Second_order", tuning_order = "Smoothness",
                sparse_tuning_u = 0, sparse_tuning_v = list(0),
                smooth_tuning_u = alpha_grid, smooth_tuning_v = list(alpha_grid),
                verbose = FALSE)
  list(Xhat = as.matrix(fit$ReconstructedData[[1]]))
}

competitors <- list(
  "Pre-smoothing HPCA" = fit_presmooth_hpca,
  "HFV-PCA"            = fit_hfvpca,
  "PMD"                = fit_pmd,
  "Two-way rSVD"       = fit_twrsvd)

###############################################################################
# 5. One simulation run: all methods on the same data set
###############################################################################
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- tryCatch(expr, error = function(e) e)
  list(val = val, time = proc.time()[["elapsed"]] - t0)
}

run_one <- function(job) {
  scn <- scenarios[job$scn, ]
  sim <- simulate_run(scn, job$rep)
  rows <- list(); sel <- list(); errs <- character(0)
  for (nm in names(configs)) {
    set.seed(sim$seed + 1L)
    tt <- timed(fit_rempca(sim$Y, configs[[nm]]))
    if (inherits(tt$val, "error")) {
      rows[[nm]] <- na_row(tt$time); errs <- c(errs, paste(nm, conditionMessage(tt$val))); next
    }
    rows[[nm]] <- c(evaluate(tt$val$Xhat, sim), time = tt$time)
    if (nm %in% c("Two-way Smooth", "Two-way Smooth+Sparse"))
      sel[[nm]] <- tt$val[c("alpha_u", "alpha_v", "gamma_u", "gamma_v")]
  }
  for (nm in names(competitors)) {
    set.seed(sim$seed + 2L)
    tt <- timed(competitors[[nm]](sim$Y))
    if (inherits(tt$val, "error")) {
      rows[[nm]] <- na_row(tt$time); errs <- c(errs, paste(nm, conditionMessage(tt$val))); next
    }
    rows[[nm]] <- c(evaluate(tt$val$Xhat, sim), time = tt$time)
  }
  res <- data.frame(scenario = scn$id, noise = scn$noise, sparsity = scn$sparsity,
                    rep = job$rep, method = names(rows),
                    do.call(rbind, rows), row.names = NULL, check.names = FALSE)
  list(res = res, sel = sel, errs = errs)
}

###############################################################################
# 6. Convergence and initialization diagnostics (Section 2.3, Table S1)
###############################################################################
# power_iterate() mirrors the package's internal power algorithm (Algorithm 1,
# conditional = FALSE) but accepts a user-supplied starting vector v0 and
# records ||v^(k) - v^(k-1)||^2 at every iteration.
power_iterate <- function(Y, S_u, S_v, gamma_u, gamma_v, v0,
                          type = "soft", thresh = 1e-10, maxit = 100) {
  thr <- getFromNamespace("sparse_pen_fun", "ReMPCA")
  v <- v0 / sqrt(sum(v0^2)); delta <- numeric(0); fit_change <- numeric(0)
  F_old <- NULL
  for (k in seq_len(maxit)) {
    Xv <- Y %*% v
    ur <- thr(y = Xv, tuning_parameter = gamma_u, type = type)
    if (all(ur == 0)) ur <- Xv
    u <- S_u %*% ur
    vn <- unlist(lapply(seq_along(idx), function(i) {
      Xu <- crossprod(Y[, idx[[i]], drop = FALSE], u)
      vr <- thr(y = Xu, tuning_parameter = as.integer(gamma_v[i]), type = type)
      if (all(vr == 0)) vr <- Xu
      as.numeric(S_v[[i]] %*% as.vector(vr))
    }))
    vn <- vn / sqrt(sum(vn^2))
    delta <- c(delta, sum((vn - v)^2)); v <- vn
    F_new <- tcrossprod(as.numeric(u), v)
    if (!is.null(F_old)) fit_change <- c(fit_change, norm(F_new - F_old, "F") / norm(F_new, "F"))
    F_old <- F_new
    if (delta[k] <= thresh) break
  }
  list(u = as.numeric(u), v = v, fit = F_new, iter = length(delta), delta = delta,
       fit_change = fit_change, converged = tail(delta, 1) <= thresh)
}

smoother <- function(td, alpha) {
  getFromNamespace("get.pen", "ReMPCA")(td = td, alpha = alpha, type = "Second_order")$S.alpha
}
mat_sqrt <- function(S) {
  e <- eigen(S, symmetric = TRUE); e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
}

conv_one <- function(job) {
  scn <- scenarios[job$scn, ]
  sim <- simulate_run(scn, job$rep)
  out <- list()
  for (nm in c("Two-way Smooth", "Two-way Smooth+Sparse")) {
    set.seed(sim$seed + 1L)
    fit <- fit_rempca(sim$Y, configs[[nm]])
    S_u <- smoother(grid_s, fit$alpha_u)
    S_v <- list(smoother(grid_t1, fit$alpha_v[1]), diag(m[["X2"]]),
                smoother(grid_t3, fit$alpha_v[3]))
    g_u <- fit$gamma_u; g_v <- fit$gamma_v
    v_svd <- svd(sim$Y, nu = 0, nv = 1)$v[, 1]
    ref <- power_iterate(sim$Y, S_u, S_v, g_u, g_v, v_svd)
    set.seed(sim$seed + 3L)
    rnd <- lapply(seq_len(n_start), function(r)
      power_iterate(sim$Y, S_u, S_v, g_u, g_v, rnorm(p)))
    ang <- vapply(rnd, function(z) sin_angle(z$v, ref$v), numeric(1))
    resid <- function(z) norm(sim$Y - z$fit, "F")
    relerr <- function(z) norm(z$fit - sim$X, "F") / norm(sim$X, "F")
    res_ref <- resid(ref); res_rnd <- vapply(rnd, resid, numeric(1))
    err_ref <- relerr(ref); err_rnd <- vapply(rnd, relerr, numeric(1))
    last10 <- function(z) if (length(z$fit_change)) max(tail(z$fit_change, 10)) else 0
    closed <- NA_real_; gap <- NA_real_
    if (nm == "Two-way Smooth") {
      # Smoothing only: Algorithm 1 is the power method for S_v X^T S_u X, whose
      # leading eigenvector is S_v^{1/2} b_1, with b_1 the leading right singular
      # vector of S_u^{1/2} X S_v^{1/2}; the convergence rate is (s_2 / s_1)^2.
      Sv_half <- as.matrix(Matrix::bdiag(lapply(S_v, mat_sqrt)))
      sb <- svd(mat_sqrt(S_u) %*% sim$Y %*% Sv_half, nu = 0, nv = 1)
      closed <- sin_angle(Sv_half %*% sb$v[, 1], ref$v)
      gap <- (sb$d[2] / sb$d[1])^2
    }
    z <- data.frame(
      scenario = scn$id, noise = scn$noise, sparsity = scn$sparsity, rep = job$rep,
      method = nm, iter_svd = ref$iter, conv_svd = ref$converged,
      agree_fit = sin_angle(ref$v, fit$v),
      iter_rand_med = median(vapply(rnd, `[[`, numeric(1), "iter")),
      iter_rand_max = max(vapply(rnd, `[[`, numeric(1), "iter")),
      conv_rand = mean(vapply(rnd, `[[`, logical(1), "converged")),
      same_sol = mean(ang < 1e-4), max_sin_rand = max(ang),
      rand_better_fit = mean(res_rnd < res_ref - 1e-8 * res_ref),
      resid_svd = res_ref, resid_rand_min = min(res_rnd),
      relerr_svd = err_ref, relerr_rand_med = median(err_rnd), relerr_rand_max = max(err_rnd),
      fit_change_last10 = last10(ref),
      closed_form = closed, rate = gap)
    if (paste(scn$id, job$rep, sep = ":") %in% trace_runs)
      attr(z, "trace") <- list(scenario = scn$id, noise = scn$noise, rep = job$rep,
                               converged = ref$converged, iter = ref$iter, svd = ref$delta,
                               rand = lapply(rnd, `[[`, "delta"))
    out[[nm]] <- z
  }
  out
}

###############################################################################
# 7. Run
###############################################################################
jobs <- unlist(lapply(run_ids, function(s) lapply(seq_len(n_rep), function(r)
  list(scn = s, rep = r))), recursive = FALSE)

if (stage %in% c("all", "fit")) {
  message(sprintf("Fitting %d data sets x %d methods on %d workers ...",
                  length(jobs), length(configs) + length(competitors), n_cores))
  t0 <- Sys.time()
  fits <- mclapply(jobs, run_one, mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- vapply(fits, inherits, logical(1), "try-error")
  if (any(bad)) warning(sum(bad), " runs failed: ", paste(unique(unlist(fits[bad])), collapse = "; "))
  errs <- unique(unlist(lapply(fits[!bad], `[[`, "errs")))
  if (length(errs)) writeLines(errs, file.path(out_dir, "simulation_errors.txt"))
  results <- do.call(rbind, lapply(fits[!bad], `[[`, "res"))
  sel <- do.call(rbind, lapply(fits[!bad], function(f) {
    s <- f$sel[["Two-way Smooth+Sparse"]]
    if (is.null(s)) return(NULL)
    data.frame(alpha_u = s$alpha_u, alpha_v1 = s$alpha_v[1], alpha_v3 = s$alpha_v[3],
               gamma_u = s$gamma_u, gamma_v1 = s$gamma_v[1], gamma_v2 = s$gamma_v[2],
               gamma_v3 = s$gamma_v[3])
  }))
  write.csv(results, file.path(out_dir, "simulation_results.csv"), row.names = FALSE)
  if (!is.null(sel)) write.csv(sel, file.path(out_dir, "simulation_selected_tuning_TwoWaySmoothSparse.csv"),
                               row.names = FALSE)
  message("  done in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
}

if (stage %in% c("all", "conv")) {
  cjobs <- unlist(lapply(run_ids, function(s) lapply(seq_len(min(n_conv, n_rep)),
    function(r) list(scn = s, rep = r))), recursive = FALSE)
  message(sprintf("Convergence diagnostics: %d data sets x %d random starts ...",
                  length(cjobs), n_start))
  conv <- mclapply(cjobs, conv_one, mc.cores = n_cores, mc.preschedule = FALSE)
  conv <- conv[!vapply(conv, inherits, logical(1), "try-error")]
  traces <- list()
  for (cc in conv) for (z in cc) if (!is.null(attr(z, "trace")))
    traces[[paste(z$method, z$scenario, z$rep, sep = "|")]] <- attr(z, "trace")
  conv_df <- do.call(rbind, lapply(conv, function(cc) do.call(rbind, cc)))
  write.csv(conv_df, file.path(out_dir, "convergence_diagnostics.csv"), row.names = FALSE)
  saveRDS(traces, file.path(out_dir, "convergence_traces.rds"))
}

if (stage == "traces") {
  tj <- lapply(strsplit(trace_runs, ":"), function(a) list(scn = as.integer(a[1]), rep = as.integer(a[2])))
  conv <- mclapply(tj, conv_one, mc.cores = min(length(tj), n_cores))
  traces <- list()
  for (cc in conv) for (z in cc) traces[[paste(z$method, z$scenario, z$rep, sep = "|")]] <- attr(z, "trace")
  saveRDS(traces, file.path(out_dir, "convergence_traces.rds"))
}

# Table S1: convergence and initialization summary by configuration and noise level
if (file.exists(file.path(out_dir, "convergence_diagnostics.csv"))) {
  cd <- read.csv(file.path(out_dir, "convergence_diagnostics.csv"))
  cd$noise <- factor(cd$noise, levels = names(noise_levels))
  tabS1 <- do.call(rbind, lapply(split(cd, list(cd$method, cd$noise), drop = TRUE), function(d) {
    nc <- d[!d$conv_svd, ]
    data.frame(method = d$method[1], noise = as.character(d$noise[1]), n_datasets = nrow(d),
               pct_converged_svd = 100 * mean(d$conv_svd),
               median_iter_svd = median(d$iter_svd), max_iter_svd = max(d$iter_svd),
               pct_random_same_solution = 100 * mean(d$same_sol),
               pct_runs_random_better_fit = 100 * mean(d$resid_rand_min < d$resid_svd * (1 - 1e-3)),
               nonconv_fit_change_min = if (nrow(nc)) min(nc$fit_change_last10) else NA,
               nonconv_fit_change_max = if (nrow(nc)) max(nc$fit_change_last10) else NA,
               median_rate = if (all(is.na(d$rate))) NA else median(d$rate, na.rm = TRUE),
               max_sin_closed_form_converged = if (all(is.na(d$closed_form))) NA
                                               else max(d$closed_form[d$conv_svd], na.rm = TRUE))
  }))
  write.csv(tabS1, file.path(out_dir, "tableS1_convergence_summary.csv"), row.names = FALSE)
}

###############################################################################
# 8. Tables
###############################################################################
method_levels <- c(names(configs), names(competitors))
tex_label <- c("No penalty" = "No penalty (HPCA)",
               "Smooth+Sparse u" = "Smooth+Sparse $\\bm{u}$",
               "Smooth+Sparse v" = "Smooth+Sparse $\\bm{v}$",
               "Two-way Smooth" = "Two-way Smooth", "Two-way Sparse" = "Two-way Sparse",
               "Two-way Smooth+Sparse" = "Two-way Smooth+Sparse",
               "Pre-smoothing HPCA" = "Pre-smoothing HPCA",
               "HFV-PCA" = "HFV-PCA", "PMD" = "PMD", "Two-way rSVD" = "Two-way rSVD")
fmt <- function(x) sprintf("%.3f", x)

if (stage %in% c("all", "fit", "figs") &&
    file.exists(file.path(out_dir, "simulation_results.csv"))) {
  results <- read.csv(file.path(out_dir, "simulation_results.csv"), check.names = FALSE)
  results$noise    <- factor(results$noise, levels = names(noise_levels))
  results$sparsity <- factor(results$sparsity, levels = names(sparsity_levels))
  results$method   <- factor(results$method, levels = method_levels)

  summ <- aggregate(cbind(rel_err, sin_u, sin_v, f1_zero_u, f1_zero_v, time) ~
                      noise + sparsity + method, data = results, na.action = na.omit,
                    FUN = function(x) c(mean = mean(x), sd = sd(x), median = median(x)))
  summ <- do.call(data.frame, summ)
  write.csv(summ, file.path(out_dir, "simulation_summary.csv"), row.names = FALSE)

  latex_rows <- function(methods, metric = "rel_err", labels = tex_label) {
    lines <- character(0)
    mu <- paste0(metric, ".mean"); s <- paste0(metric, ".sd")
    for (nz in names(noise_levels)) {
      block <- summ[summ$noise == nz & summ$method %in% methods, ]
      best <- tapply(block[[mu]], as.character(block$sparsity), min)
      for (k in seq_along(methods)) {
        mt <- methods[k]
        cells <- vapply(names(sparsity_levels), function(sp) {
          r <- block[block$method == mt & block$sparsity == sp, ]
          if (!nrow(r)) return("--")
          txt <- sprintf("%s (%s)", fmt(r[[mu]]), fmt(r[[s]]))
          if (abs(r[[mu]] - best[[sp]]) < 5e-4) txt <- paste0("\\textbf{", txt, "}")
          txt
        }, character(1))
        lab <- if (k == 1) sprintf("\\multirow{%d}{*}{%s noise}", length(methods), nz) else ""
        lines <- c(lines, sprintf("%s & %s & %s \\\\", lab, labels[[mt]],
                                  paste(cells, collapse = " & ")))
      }
      lines <- c(lines, "\\midrule")
    }
    lines[length(lines)] <- "\\bottomrule"
    lines
  }
  bench <- c("No penalty", names(competitors), "Two-way Smooth", "Two-way Smooth+Sparse")
  writeLines(latex_rows(names(configs)), file.path(out_dir, "table2_rows.tex"))
  bench_label <- tex_label
  bench_label[c("No penalty", "Two-way Smooth", "Two-way Smooth+Sparse")] <-
    c("HPCA (no penalty)", "ReHPCA, two-way smooth", "ReHPCA, two-way smooth+sparse")
  writeLines(latex_rows(bench, labels = bench_label), file.path(out_dir, "table3_benchmark_rows.tex"))
  writeLines(latex_rows(bench, "sin_u", labels = bench_label),
             file.path(out_dir, "tableS_benchmark_sinu_rows.tex"))
  writeLines(latex_rows(bench, "sin_v", labels = bench_label),
             file.path(out_dir, "tableS_benchmark_sinv_rows.tex"))
  # median computing time per fit (seconds), including tuning
  tm <- aggregate(time ~ method, data = results[results$method %in% bench, ], FUN = median)
  write.csv(tm, file.path(out_dir, "tableS_benchmark_time.csv"), row.names = FALSE)
}

###############################################################################
# 9. Figures (drawn at final print size: 6.3 in = JDS text width)
###############################################################################
pal_cfg <- c("No penalty" = "#999999", "Smooth+Sparse u" = "#E69F00",
             "Smooth+Sparse v" = "#56B4E9", "Two-way Smooth" = "#009E73",
             "Two-way Sparse" = "#CC79A7", "Two-way Smooth+Sparse" = "#0072B2")
theme_paper <- theme_bw(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey88"),
        strip.background = element_rect(fill = "grey94", colour = "grey60", linewidth = 0.3),
        strip.text = element_text(size = 8.5, face = "bold"),
        axis.text = element_text(size = 7.5, colour = "grey20"),
        axis.title = element_text(size = 8.5),
        legend.text = element_text(size = 8), legend.title = element_blank(),
        legend.key.size = grid::unit(0.35, "cm"), legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(3, 4, 3, 3))

if (stage %in% c("all", "fit", "figs") && exists("results")) {
  r2 <- results[results$method %in% names(configs), ]
  r2$method <- droplevels(r2$method)
  r2$noise  <- factor(paste(r2$noise, "noise"), levels = paste(names(noise_levels), "noise"))
  g2 <- ggplot(r2, aes(x = sparsity, y = rel_err, fill = method)) +
    geom_boxplot(width = 0.78, position = position_dodge(0.86), linewidth = 0.25,
                 outlier.size = 0.35, outlier.stroke = 0.2) +
    facet_wrap(~ noise, nrow = 1, scales = "free_y") +
    scale_y_log10() +
    scale_fill_manual(values = pal_cfg,
                      labels = c("No penalty",
                                 expression("Smooth+Sparse " * bolditalic(u)),
                                 expression("Smooth+Sparse " * bolditalic(v)),
                                 "Two-way Smooth", "Two-way Sparse", "Two-way Smooth+Sparse")) +
    labs(x = "Sparsity level", y = "Relative reconstruction error (log scale)") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_paper + theme(legend.text = element_text(size = 7.5),
                        legend.key.width = grid::unit(0.3, "cm"))
  ggsave(file.path(out_dir, "fig2_reconstruction_error_boxplot.pdf"), g2,
         width = 6.3, height = 3.2, units = "in")
}

if (stage %in% c("all", "figs")) {
  # Figure 3: representative run (medium noise, medium sparsity, run 1)
  # The displayed run is the one with the median ReHPCA (Two-way Smooth+Sparse)
  # error among the runs of the medium-noise / medium-sparsity scenario.
  scn <- scenarios[scenarios$noise == "Medium" & scenarios$sparsity == "Medium", ]
  rep_show <- 1L
  if (file.exists(file.path(out_dir, "simulation_results.csv"))) {
    rr <- read.csv(file.path(out_dir, "simulation_results.csv"))
    rr <- rr[rr$scenario == scn$id & rr$method == "Two-way Smooth+Sparse" & !is.na(rr$rel_err), ]
    rr <- rr[order(rr$rel_err), ]
    rep_show <- rr$rep[ceiling(nrow(rr) / 2)]
  }
  message("Figure 3 shows run ", rep_show)
  writeLines(as.character(rep_show), file.path(out_dir, "fig3_run_shown.txt"))
  sim <- simulate_run(scn, rep_show)
  set.seed(sim$seed + 1L)
  f_re <- fit_rempca(sim$Y, configs[["Two-way Smooth+Sparse"]])
  sv0 <- svd(sim$Y, nu = 1, nv = 1)
  align <- function(x, ref) { x <- x / sqrt(sum(x^2)); if (sum(x * ref) < 0) -x else x }
  long_mat <- function(M, lab) data.frame(row = rep(seq_len(n), p), col = rep(seq_len(p), each = n),
                                          value = as.vector(M), panel = lab)
  lim <- max(abs(sim$X))
  hm <- rbind(long_mat(sim$X, "True signal"), long_mat(sim$Y, "Noisy data"),
              long_mat(f_re$Xhat, "ReHPCA fit"))
  hm$panel <- factor(hm$panel, levels = c("True signal", "Noisy data", "ReHPCA fit"))
  hm$value <- pmax(pmin(hm$value, lim), -lim)
  block_mid <- c(mean(idx$X1), mean(idx$X2), mean(idx$X3))
  g3a <- ggplot(hm, aes(col, row, fill = value)) +
    geom_raster() +
    geom_vline(xintercept = c(40.5, 60.5), colour = "black", linewidth = 0.3) +
    facet_wrap(~ panel, nrow = 1) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         limits = c(-lim, lim), name = NULL) +
    scale_x_continuous(breaks = block_mid, labels = c("X1", "X2", "X3"), expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0), breaks = c(1, 50, 100, 150)) +
    labs(x = NULL, y = "Row index") +
    theme_paper + theme(legend.position = "right", legend.key.height = grid::unit(0.8, "cm"),
                        legend.key.width = grid::unit(0.25, "cm"), panel.grid.major = element_blank(),
                        axis.ticks.x = element_blank())
  vec_df <- function(est_svd, est_re, truth) {
    data.frame(x = rep(seq_along(truth), 3),
               value = c(truth, align(est_svd, truth), align(est_re, truth)),
               which = factor(rep(c("Truth", "HPCA (no penalty)", "ReHPCA"), each = length(truth)),
                              levels = c("Truth", "HPCA (no penalty)", "ReHPCA")))
  }
  pal_vec <- c("Truth" = "black", "HPCA (no penalty)" = "#999999", "ReHPCA" = "#0072B2")
  lw_vec  <- c("Truth" = 0.45, "HPCA (no penalty)" = 0.35, "ReHPCA" = 0.6)
  du <- vec_df(sv0$u[, 1], f_re$u, sim$u)
  g3b <- ggplot(du, aes(x, value, colour = which, linewidth = which)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
    geom_line() +
    scale_colour_manual(values = pal_vec) +
    scale_linewidth_manual(values = lw_vec, guide = "none") +
    labs(x = "Row index", y = expression("Row component " * bolditalic(u))) +
    theme_paper
  dv <- vec_df(sv0$v[, 1], f_re$v, sim$v)
  dv$block <- factor(rep(rep(names(m), m), 3), levels = names(m))
  g3c <- ggplot(dv, aes(x, value, colour = which)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
    geom_vline(xintercept = c(40.5, 60.5), colour = "grey50", linewidth = 0.25, linetype = 2) +
    geom_line(data = dv[dv$block != "X2", ],
              aes(group = interaction(which, block), linewidth = which)) +
    geom_point(data = dv[dv$block == "X2", ], size = 0.7, show.legend = FALSE) +
    annotate("text", x = block_mid, y = Inf, vjust = 1.25, size = 2.6, lineheight = 0.85,
             label = c("X1\n(functional)", "X2\n(vector)", "X3\n(functional)")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    scale_colour_manual(values = pal_vec) +
    scale_linewidth_manual(values = lw_vec, guide = "none") +
    labs(x = "Column index", y = expression("Loading " * bolditalic(v))) +
    theme_paper
  g3 <- g3a / ((g3b | g3c) + plot_layout(widths = c(1, 1.6), guides = "collect") &
                 theme(legend.position = "bottom")) +
    plot_layout(heights = c(1, 0.9))
  ggsave(file.path(out_dir, "fig3_representative_run.pdf"), g3, width = 6.3, height = 4.6, units = "in")
}

if (stage %in% c("all", "conv", "figs") &&
    file.exists(file.path(out_dir, "convergence_traces.rds"))) {
  traces <- readRDS(file.path(out_dir, "convergence_traces.rds"))
  if (length(traces)) {
    tr <- do.call(rbind, lapply(names(traces), function(key) {
      z <- traces[[key]]; nm <- sub("\\|.*$", "", key)
      lab <- sprintf("%s noise, run %d", z$noise, z$rep)
      rbind(data.frame(method = nm, runlab = lab, start = "SVD", run = 0,
                       iter = seq_along(z$svd), delta = z$svd),
            do.call(rbind, lapply(seq_along(z$rand), function(r)
              data.frame(method = nm, runlab = lab, start = "Random", run = r,
                         iter = seq_along(z$rand[[r]]), delta = z$rand[[r]]))))
    }))
    status <- do.call(rbind, lapply(names(traces), function(key) {
      z <- traces[[key]]
      data.frame(method = sub("\\|.*$", "", key), runlab = sprintf("%s noise, run %d", z$noise, z$rep),
                 txt = if (z$converged) sprintf("SVD start: converged\nin %d iterations", z$iter)
                       else "SVD start: not converged\nin 100 iterations")
    }))
    run_levels <- unique(tr$runlab[order(tr$runlab != "Medium noise, run 1")])
    tr$runlab <- factor(tr$runlab, levels = run_levels)
    status$runlab <- factor(status$runlab, levels = run_levels)
    tr$delta <- pmax(tr$delta, 1e-17)
    tr$start <- factor(tr$start, levels = c("Random", "SVD"))
    gS <- ggplot(tr, aes(iter, delta, group = interaction(start, run), colour = start)) +
      geom_hline(yintercept = 1e-10, linetype = 2, linewidth = 0.25) +
      geom_line(data = tr[tr$start == "Random", ], linewidth = 0.3) +
      geom_line(data = tr[tr$start == "SVD", ], linewidth = 0.55) +
      geom_text(data = status, aes(x = 100, y = 1e-15, label = txt), inherit.aes = FALSE,
                hjust = 1, vjust = 0, size = 2.4, lineheight = 0.9) +
      facet_grid(method ~ runlab) +
      scale_y_log10(limits = c(1e-17, 10), breaks = 10^c(-15, -10, -5, 0)) +
      scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
      scale_colour_manual(values = c(Random = "grey65", SVD = "#D55E00"),
                          labels = c("Random start", "SVD start")) +
      labs(x = "Iteration",
           y = expression("||" * bolditalic(v)^(k) - bolditalic(v)^(k-1) * "||"^2)) +
      theme_paper + theme(strip.text = element_text(size = 7.5))
    ggsave(file.path(out_dir, "figS_convergence_traces.pdf"), gS, width = 6.3, height = 3.8, units = "in")
  }
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_simulation.txt"))
