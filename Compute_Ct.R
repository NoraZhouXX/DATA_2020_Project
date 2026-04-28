rm(list = ls())

library(glmnet)

LVAR <- function(y, p = 1) {

  y    <- scale(y)
  colnames(y) <- make.names(colnames(y))
  y.orig <- y

  obs <- dim(y)[1]
  K   <- dim(y)[2]
  sample <- obs - p

  # Build lagged regressors matrix
  ylags <- embed(y, dimension = p + 1)[, -(1:K)]
  temp1 <- NULL
  for (i in 1:p) {
    temp  <- paste(colnames(y), ".l", i, sep = "")
    temp1 <- c(temp1, temp)
  }
  colnames(ylags) <- temp1

  yend     <- y[-c(1:p), ]
  rhs      <- ylags
  datamat  <- ylags
  equation <- list()
  n_coef   <- dim(datamat)[2]

  # Elastic-Net (alpha = 0.5), following Lee (2025)
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

  result <- list(
    varresult = equation,
    datamat   = data.frame(cbind(yend, rhs)),
    datamat1  = rhs,
    y         = y.orig,
    p         = p,
    K         = K,
    obs       = sample,
    totobs    = sample + p,
    type      = "ELASSO"
  )
  return(result)
}


# Bcoef: Extract coefficient matrix B from a fitted LVAR object
Bcoef <- function(x) {
  y.names <- colnames(x$datamat[, c(1:x$K)])
  Z  <- x$datamat[, -c(1:x$K)]
  Z1 <- as.numeric(dim(Z)[2])
  B  <- matrix(0, nrow = x$K, ncol = Z1)

  for (i in 1:x$K) {
    B[i, ] <- t(as.vector(x$varresult[[i]]))
  }
  colnames(B) <- colnames(Z)
  rownames(B) <- y.names
  return(B)
}


# Acoef: Extract list of companion matrices A_1, ..., A_p
Acoef <- function(x) {
  K   <- x$K
  p   <- x$p
  A   <- Bcoef(x)[, 1:(K * p)]
  As  <- list()
  start <- seq(1,   p * K, K)
  end   <- seq(K,   p * K, K)

  for (i in 1:p) {
    As[[i]] <- matrix(A[, start[i]:end[i]], nrow = K, ncol = K)
    rownames(As[[i]]) <- rownames(A)
    colnames(As[[i]]) <- colnames(A[, start[i]:end[i]])
  }
  return(As)
}


# Phi: Compute MA coefficient matrices Phi_0, ..., Phi_{nstep}
Phi <- function(x, nstep = 10, ...) {
  nstep <- abs(as.integer(nstep))
  K <- x$K
  p <- x$p
  A <- as.array(Acoef(x))

  if (nstep > p) {
    As <- array(0, dim = c(K, K, nstep))
    for (i in (p + 1):(nstep)) {
      As[, , i] <- matrix(0, nrow = K, ncol = K)
    }
  } else {
    As <- array(0, dim = c(K, K, p))
  }
  for (i in 1:p) {
    As[, , i] <- A[[i]]
  }

  Phi_arr <- array(0, dim = c(K, K, nstep))

  if (nstep == 1) {
    Phi_arr[, , 1] <- diag(K)
  }
  if (nstep == 2) {
    Phi_arr[, , 1] <- diag(K)
    Phi_arr[, , 2] <- Phi_arr[, , 1] %*% As[, , 1]
  }
  if (nstep > 2) {
    Phi_arr[, , 1] <- diag(K)
    Phi_arr[, , 2] <- Phi_arr[, , 1] %*% As[, , 1]
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
  return(Phi_arr)
}


# res: Compute residuals from a fitted LVAR object
res <- function(x) {
  k   <- x$K
  y   <- x$datamat[, 1:k]
  X   <- x$datamat1
  res <- y - as.matrix(X) %*% t(as.matrix(Bcoef(x)))
  return(res)
}


# fevd_generalised: Generalised Forecast Error Variance Decomposition
# Koop, Pesaran & Potter (1996); Pesaran & Shin (1998)
# Returns a K x K matrix (row-normalised by default)
fevd_generalised <- function(model, n.ahead = 10, normalize = TRUE) {
  A       <- Phi(model, n.ahead)
  epsilon <- as.matrix(res(model))
  Sigma   <- t(epsilon) %*% epsilon / (model$obs)

  gi     <- array(0, dim(A))
  sigmas <- sqrt(diag(Sigma))
  for (j in 1:dim(A)[3]) {
    gi[, , j] <- t(t(A[, , j] %*% Sigma) / sigmas)
  }

  d <- array(0, dim(A)[c(2, 3)])
  for (j in 1:dim(d)[2]) {
    d[, j] <- diag(A[, , j] %*% Sigma %*% t(A[, , j]))
  }

  num  <- apply(gi^2, 1:2, sum)
  den  <- c(apply(d, 1, sum))
  fevd <- num / den

  if (normalize) {
    return(fevd / apply(fevd, 1, sum))
  } else {
    return(fevd)
  }
}


# sptable2012: Build Diebold-Yilmaz (2012) spillover table
# NOTE: uses global variable 'data' for column names
# Args:
#   sp : K x K GFEVD matrix (row-normalised, scaled by 100)
# Returns: (K+2) x (K+1) spillover table; bottom-right cell = Ct
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


# MAIN: COMPUTE ROLLING CONNECTEDNESS INDEX

# Load data
raw        <- read.csv("EI_DATA.csv", header = TRUE)

# Columns 4-54: 51 state housing return series (AK through WY)
sym        <- raw[, c(4:60)]
data_House <- as.matrix(sym[, c(1:51)])

# 'data' referenced globally by sptable2012() for column names
dat  <- data_House
data <- sym[, c(1:51)]

cat("Data loaded:", nrow(dat), "monthly observations,", ncol(dat), "states\n")
cat("Sample: Year", raw$Year[1], "Month", raw$Month[1],
    "to Year", raw$Year[nrow(raw)], "Month", raw$Month[nrow(raw)], "\n\n")

# Rolling window parameters
window  <- 120   # 10-year rolling window (months)
p       <- 1     # VAR lag order
n.ahead <- 10    # GFEVD forecast horizon

last <- nrow(dat) - window - p + 1
cat("Rolling windows:", last, "\n")
cat("Ct dated at center of each window; series spans 1981Q1 to 2016Q4\n\n")

# Storage
matspill <- matrix(0, last, 1)

# Rolling window loop
cat("Starting computation (est. 1-3 hours)...\n\n")
start_time <- proc.time()

for (i in seq_len(last)) {

  data1       <- dat[i:(window - 1 + p + i), ]
  Rolling_VAR <- LVAR(data1, p = p)
  rolling_fe  <- fevd_generalised(Rolling_VAR, n.ahead = n.ahead) * 100
  th0         <- sptable2012(rolling_fe)

  # Total connectedness = bottom-right cell (Net row, From column)
  matspill[i] <- th0[(Rolling_VAR$K + 2), (Rolling_VAR$K + 1)]

  if (i %% 10 == 0 || i == 1 || i == last) {
    elapsed <- (proc.time() - start_time)[3]
    eta     <- if (i > 1) elapsed / i * (last - i) else NA
    cat(sprintf("Window %d / %d  |  Ct = %.3f  |  Elapsed: %.1f min  |  ETA: %.1f min\n",
                i, last, matspill[i], elapsed / 60,
                ifelse(is.na(eta), NA, eta / 60)))
  }
}

cat("\nDone! Total:", round((proc.time() - start_time)[3] / 60, 1), "min\n")

# Save results
write.csv(matspill, file = "Ct_monthly.csv", row.names = TRUE)
cat("Saved: Ct_monthly.csv\n")

