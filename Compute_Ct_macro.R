# =============================================================================
# Compute_Ct_macro.R
# Connectedness index WITH macroeconomic controls
#
# Augments the rolling elastic-net VAR with macro variables already present
# in EI_DATA.csv (columns 52:end of sym block).  The total connectedness
# index is extracted from the K_house × K_house (51 × 51) housing sub-block
# of the full-system GFEVD, following the approach in Lee (2025).
#
# OUTPUT: Ct_macro_monthly.csv
# =============================================================================
rm(list = ls())
library(glmnet)

# ---------------------------------------------------------------------------
# LVAR: elastic-net VAR  (same as Compute_Ct.R)
# ---------------------------------------------------------------------------
LVAR <- function(y, p = 1) {
  y    <- scale(y)
  colnames(y) <- make.names(colnames(y))
  y.orig <- y

  obs    <- nrow(y);  K <- ncol(y)
  sample <- obs - p

  ylags <- embed(y, dimension = p + 1)[, -(1:K)]
  temp1 <- NULL
  for (i in 1:p) temp1 <- c(temp1, paste(colnames(y), ".l", i, sep = ""))
  colnames(ylags) <- temp1

  yend    <- y[-c(1:p), ]
  datamat <- ylags
  equation <- list()
  n_coef   <- ncol(datamat)

  for (i in 1:K) {
    yi      <- yend[, i]
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
       p         = p,  K = K,
       obs       = sample,  totobs = sample + p,
       type      = "ELASSO")
}

Bcoef <- function(x) {
  y.names <- colnames(x$datamat[, 1:x$K])
  Z  <- x$datamat[, -(1:x$K)]
  B  <- matrix(0, nrow = x$K, ncol = ncol(Z))
  for (i in 1:x$K) B[i, ] <- as.vector(x$varresult[[i]])
  colnames(B) <- colnames(Z);  rownames(B) <- y.names
  B
}

Acoef <- function(x) {
  K <- x$K;  p <- x$p
  A <- Bcoef(x)[, 1:(K * p)]
  As <- list()
  start <- seq(1, p * K, K);  end <- seq(K, p * K, K)
  for (i in 1:p) {
    As[[i]] <- matrix(A[, start[i]:end[i]], nrow = K, ncol = K)
    rownames(As[[i]]) <- rownames(A)
    colnames(As[[i]]) <- colnames(A[, start[i]:end[i]])
  }
  As
}

Phi <- function(x, nstep = 10) {
  nstep <- abs(as.integer(nstep))
  K <- x$K;  p <- x$p
  A <- as.array(Acoef(x))

  As <- array(0, dim = c(K, K, max(nstep, p)))
  for (i in 1:p) As[, , i] <- A[[i]]

  Phi_arr <- array(0, dim = c(K, K, nstep))
  if (nstep >= 1) Phi_arr[, , 1] <- diag(K)
  if (nstep >= 2) Phi_arr[, , 2] <- Phi_arr[, , 1] %*% As[, , 1]
  if (nstep > 2) {
    for (i in 3:nstep) {
      tmp1 <- Phi_arr[, , 1] %*% As[, , i - 1]
      tmp2 <- matrix(0, K, K)
      idx  <- (i - 2):1
      for (j in 1:(i - 2))
        tmp2 <- tmp2 + Phi_arr[, , j + 1] %*% As[, , idx[j]]
      Phi_arr[, , i] <- tmp1 + tmp2
    }
  }
  Phi_arr
}

res_fn <- function(x) {
  y <- x$datamat[, 1:x$K]
  as.matrix(y) - as.matrix(x$datamat1) %*% t(Bcoef(x))
}

fevd_generalised <- function(model, n.ahead = 10) {
  A       <- Phi(model, n.ahead)
  epsilon <- res_fn(model)
  Sigma   <- t(epsilon) %*% epsilon / model$obs

  gi     <- array(0, dim(A))
  sigmas <- sqrt(diag(Sigma))
  for (j in 1:dim(A)[3])
    gi[, , j] <- t(t(A[, , j] %*% Sigma) / sigmas)

  d <- array(0, dim(A)[c(2, 3)])
  for (j in 1:ncol(d))
    d[, j] <- diag(A[, , j] %*% Sigma %*% t(A[, , j]))

  num  <- apply(gi^2, 1:2, sum)
  den  <- c(apply(d, 1, sum))
  fevd <- num / den
  fevd / apply(fevd, 1, sum)   # row-normalise
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
raw  <- read.csv("EI_DATA.csv", header = TRUE)
sym  <- raw[, c(4:ncol(raw))]      # all variable columns

K_house <- 51                       # number of housing series (always first 51)
dat_all <- as.matrix(sym)           # full system (housing + macro)
K_total <- ncol(dat_all)

cat("Data loaded:", nrow(dat_all), "monthly obs |",
    K_house, "housing +", K_total - K_house, "macro variables\n")

# Rolling window parameters  (identical to benchmark)
window  <- 120
p       <- 1
n.ahead <- 10

last <- nrow(dat_all) - window - p + 1
cat("Rolling windows:", last, "\n\n")

matspill <- matrix(0, last, 1)

cat("Starting computation with macro controls...\n\n")
start_time <- proc.time()

for (i in seq_len(last)) {

  data1       <- dat_all[i:(window - 1 + p + i), ]
  Rolling_VAR <- LVAR(data1, p = p)

  # Full-system GFEVD (K_total × K_total, row-normalised)
  fevd_full <- fevd_generalised(Rolling_VAR, n.ahead = n.ahead) * 100

  # Extract housing sub-block (first K_house rows/cols) and re-normalise rows
  sp_sub <- fevd_full[1:K_house, 1:K_house]
  sp_sub <- sp_sub / rowSums(sp_sub) * 100   # re-normalise to sum to 100

  # Total connectedness from housing sub-block (DY 2012 definition)
  matspill[i] <- (sum(sp_sub) - sum(diag(sp_sub))) / sum(sp_sub) * 100

  if (i %% 10 == 0 || i == 1 || i == last) {
    elapsed <- (proc.time() - start_time)[3]
    eta     <- if (i > 1) elapsed / i * (last - i) else NA
    cat(sprintf("Window %d / %d  |  Ct_macro = %.3f  |  Elapsed: %.1f min  |  ETA: %.1f min\n",
                i, last, matspill[i], elapsed / 60,
                ifelse(is.na(eta), NA, eta / 60)))
  }
}

cat("\nDone! Total:", round((proc.time() - start_time)[3] / 60, 1), "min\n")

write.csv(matspill, file = "Ct_macro_monthly.csv", row.names = TRUE)
cat("Saved: Ct_macro_monthly.csv\n")
