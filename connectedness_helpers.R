# =============================================================================
# connectedness_helpers.R
# Shared helper functions for Diebold-Yilmaz total connectedness computation.
#
# Sourced by:
#   Compute_Ct.R         (housing-only,    51-variable VAR)
#   Compute_Ct_macro.R   (housing + macro, 57-variable VAR)
#
# Functions exported:
#   LVAR()              fit elastic-net VAR(p) via cv.glmnet (alpha = 0.5)
#   Bcoef()             stack VAR equation coefficients into matrix B
#   Acoef()             split B into list of K x K companion matrices A_1,...,A_p
#   Phi()               compute MA(inf) coefficient matrices Phi_0,...,Phi_{nstep-1}
#   res_fn()            compute VAR residuals
#   fevd_generalised()  Pesaran-Shin (1998) generalised forecast error variance
#                       decomposition; returns row-normalised K x K matrix
#   sptable2012()       Diebold-Yilmaz (2012) spillover table
#                       (used by Compute_Ct.R only)
#
# IMPORTANT:
#   sptable2012() reads a global object 'data' (a data.frame) to label its
#   rows/columns. Callers MUST set
#       data <- <data.frame whose colnames are the K series names>
#   in the global environment before calling sptable2012(). This pattern
#   matches the Lee & Ma (2025) replication code and is preserved here for
#   bit-for-bit numerical equivalence.
# =============================================================================

library(glmnet)


# -----------------------------------------------------------------------------
# LVAR: equation-by-equation Elastic-Net VAR (alpha = 0.5)
# -----------------------------------------------------------------------------
LVAR <- function(y, p = 1) {

  y    <- scale(y)
  colnames(y) <- make.names(colnames(y))
  y.orig <- y

  obs    <- nrow(y)
  K      <- ncol(y)
  sample <- obs - p

  # Build lagged regressors matrix
  ylags <- embed(y, dimension = p + 1)[, -(1:K)]
  temp1 <- NULL
  for (i in 1:p) {
    temp1 <- c(temp1, paste(colnames(y), ".l", i, sep = ""))
  }
  colnames(ylags) <- temp1

  yend     <- y[-c(1:p), ]
  datamat  <- ylags
  equation <- list()
  n_coef   <- ncol(datamat)

  # Equation-by-equation Elastic-Net (alpha = 0.5), following Lee (2025)
  for (i in 1:K) {
    yi <- yend[, i]
    set.seed(2)
    cv.out  <- cv.glmnet(datamat, yi, family = "gaussian", alpha = 0.5,
                         intercept = FALSE)
    bestlam <- cv.out$lambda.min
    equation[[colnames(yend)[i]]] <-
      predict(cv.out, type = "coefficients", s = bestlam,
              intercept = FALSE)[2:(n_coef + 1), ]
  }

  list(varresult = equation,
       datamat   = data.frame(cbind(yend, ylags)),
       datamat1  = ylags,
       y         = y.orig,
       p         = p,
       K         = K,
       obs       = sample,
       totobs    = sample + p,
       type      = "ELASSO")
}


# -----------------------------------------------------------------------------
# Bcoef: stack equation-level coefficients into a single coefficient matrix B
# -----------------------------------------------------------------------------
Bcoef <- function(x) {
  y.names <- colnames(x$datamat[, 1:x$K])
  Z  <- x$datamat[, -(1:x$K)]
  B  <- matrix(0, nrow = x$K, ncol = ncol(Z))
  for (i in 1:x$K) {
    B[i, ] <- as.vector(x$varresult[[i]])
  }
  colnames(B) <- colnames(Z)
  rownames(B) <- y.names
  B
}


# -----------------------------------------------------------------------------
# Acoef: split B into list of companion matrices A_1, ..., A_p
# -----------------------------------------------------------------------------
Acoef <- function(x) {
  K <- x$K
  p <- x$p
  A <- Bcoef(x)[, 1:(K * p)]
  As <- list()
  start <- seq(1, p * K, K)
  end   <- seq(K, p * K, K)
  for (i in 1:p) {
    As[[i]] <- matrix(A[, start[i]:end[i]], nrow = K, ncol = K)
    rownames(As[[i]]) <- rownames(A)
    colnames(As[[i]]) <- colnames(A[, start[i]:end[i]])
  }
  As
}


# -----------------------------------------------------------------------------
# Phi: VAR(p) -> MA(inf) recursion. Returns array of Phi_0, ..., Phi_{nstep-1}.
# -----------------------------------------------------------------------------
Phi <- function(x, nstep = 10) {
  nstep <- abs(as.integer(nstep))
  K <- x$K
  p <- x$p
  A <- as.array(Acoef(x))

  As <- array(0, dim = c(K, K, max(nstep, p)))
  for (i in 1:p) As[, , i] <- A[[i]]

  Phi_arr <- array(0, dim = c(K, K, nstep))
  if (nstep >= 1) Phi_arr[, , 1] <- diag(K)
  if (nstep >= 2) Phi_arr[, , 2] <- Phi_arr[, , 1] %*% As[, , 1]
  if (nstep > 2) {
    for (i in 3:nstep) {
      tmp1 <- Phi_arr[, , 1] %*% As[, , i - 1]
      tmp2 <- matrix(0, nrow = K, ncol = K)
      idx  <- (i - 2):1
      for (j in 1:(i - 2)) {
        tmp2 <- tmp2 + Phi_arr[, , j + 1] %*% As[, , idx[j]]
      }
      Phi_arr[, , i] <- tmp1 + tmp2
    }
  }
  Phi_arr
}


# -----------------------------------------------------------------------------
# res_fn: VAR residuals (named res_fn to avoid shadowing R's res in plot calls)
# -----------------------------------------------------------------------------
res_fn <- function(x) {
  y <- x$datamat[, 1:x$K]
  as.matrix(y) - as.matrix(x$datamat1) %*% t(as.matrix(Bcoef(x)))
}


# -----------------------------------------------------------------------------
# fevd_generalised: Pesaran-Shin (1998) generalised FEVD
#   Returns a K x K row-normalised matrix when normalize = TRUE (default).
#   Set normalize = FALSE to inspect the un-normalised matrix.
# -----------------------------------------------------------------------------
fevd_generalised <- function(model, n.ahead = 10, normalize = TRUE) {
  A       <- Phi(model, n.ahead)
  epsilon <- as.matrix(res_fn(model))
  Sigma   <- t(epsilon) %*% epsilon / model$obs

  gi     <- array(0, dim(A))
  sigmas <- sqrt(diag(Sigma))
  for (j in 1:dim(A)[3]) {
    gi[, , j] <- t(t(A[, , j] %*% Sigma) / sigmas)
  }

  d <- array(0, dim(A)[c(2, 3)])
  for (j in 1:ncol(d)) {
    d[, j] <- diag(A[, , j] %*% Sigma %*% t(A[, , j]))
  }

  num  <- apply(gi^2, 1:2, sum)
  den  <- c(apply(d, 1, sum))
  fevd <- num / den

  if (normalize) {
    fevd / apply(fevd, 1, sum)
  } else {
    fevd
  }
}


# -----------------------------------------------------------------------------
# sptable2012: Diebold-Yilmaz (2012) spillover table
#   sp : K x K row-normalised GFEVD matrix scaled by 100
#
# Returns: (K+2) x (K+1) spillover table:
#   - Top K x K block:    GFEVD entries
#   - Last column "From": row-sum off-diagonal -> spillover received by row i
#   - Row K+1 "To":       column-sum off-diagonal -> spillover sent by col j
#   - Row K+2 "Net":      To - From (with totSP in the last cell)
#   - Bottom-right cell:  total connectedness Ct (%)  [matspill index]
#
# NOTE: Uses global 'data' for column names. See file header.
# -----------------------------------------------------------------------------
sptable2012 <- function(sp) {
  own    <- diag(sp)
  from   <- apply(sp, 1, sum) - own
  to     <- apply(sp, 2, sum) - own
  net    <- to - from
  offsum <- sum(to)
  to     <- c(to, offsum)
  totSP  <- (sum(sp) - sum(diag(sp))) / sum(sp) * 100
  net    <- c(net, totSP)

  ab    <- cbind(sp, from)
  to    <- rbind(ab, to)
  sptab <- rbind(to, net)

  dimnames(sptab) <- list(
    c(colnames(data), "To", "Net"),
    c(colnames(data), "From")
  )
  sptab <- round(sptab, 3)
  return(sptab)
}
