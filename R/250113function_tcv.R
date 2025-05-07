library(foreach)
library(doParallel)
library(irlba)
library(doSNOW)
library(openxlsx)
library(countsplit)
library(bcv)
library(Matrix)
library(tidyr)
library(invgamma)
library(mvtnorm)
library(pracma)
library(pcaPP)
library(irlba)
library(Rfast)
library(dplyr)
library(foreach)
library(doParallel)
library(irlba)
library(doSNOW)
library(openxlsx)
library(countsplit)
library(bcv)
library(tidyr)
library(dplyr)
library(invgamma)
library(mvtnorm)
library(pracma)
library(pcaPP)
library(Rfast)
#Rcpp::sourceCpp("/Users/chengjunwang/Downloads/confirm_jmle_omp_poisson_with_intercept_missing.cpp")

##################multiDT
#' Thinning Cross-Validation for Determining the Number of Factors in Poisson Factor Models.
#'
#' @param x n by p matrix containing count responses, where n is the number of samples, p is the number of features.
#' @param K The number of folds for conducting thinning cross-validation. Default as 5.
#' @param rmax The maximum number of factors for conducting thinning cross-validation. Default as 8.
#'
#' @return The TCV error and TICV error. The number of factors corresponding to the smallest result is TCV estimator.
#' @author Zhijing Wang
multiDT <- function(x, K = 5, rmax = 8) {
  # Define required packages
  required_packages <- c(
    "foreach", "doParallel", "irlba", "doSNOW", "openxlsx",
    "countsplit", "bcv", "Matrix", "tidyr", "invgamma",
    "mvtnorm", "pracma", "pcaPP", "Rfast", "dplyr", "GFM"
  )

  # Check if packages are installed, if not prompt to install
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Please install the '", pkg, "' package: install.packages('", pkg, "')"))
    }
  }

  # Load packages while suppressing messages and warnings
  suppressPackageStartupMessages({
    for (pkg in required_packages) {
      library(pkg, character.only = TRUE)
    }
  })

  n <- nrow(x)
  p <- ncol(x)

  # Suppress the message from countsplit
  suppressMessages({
    xsplit <- countsplit::countsplit(x, folds = K)
  })
  obj_mat <- matrix(0, nrow = rmax, ncol = K)
  co_mat <- matrix(0, nrow = rmax, ncol = K)
  theta0 <- cbind(rep(1, n), matrix(runif(n*rmax, 0, 1), n, rmax))
  A0 <- matrix(runif(p*(rmax+1), -1, 1), p, (rmax+1))

  for (r in (1:rmax)) {
    QB <- matrix()
    objB <- c()
    co_es <- c()

    for (k in (1:K)) {
      XA <- x - as.matrix(xsplit[[k]])
      XB <- x - XA

      lamB <- tryCatch({
        tra <- suppressMessages(GFM::gfm(list(XA), types = c("poisson"),
                                         maxIter = 30, q = r, algorithm = "AM"))
        lamA <- exp(cbind(tra$hH, 1) %*% t(cbind(tra$hB, tra$hmu)))
        lamA/(K-1)
      }, error = function(e) {
        true_matrix <- matrix(TRUE, nrow = n, ncol = p)
        jic.res <- confirm_CJMLE_poisson_cpp(XA, true_matrix,
                                             theta0[, 1:(r+1)],
                                             A0[, 1:(r+1)],
                                             matrix(TRUE, p, (r+1)),
                                             C = 100,
                                             tol = .01/(n*p))
        lamA <- exp((jic.res$theta) %*% t(jic.res$A))
        lamA/(K-1)
      })

      eps1 <- 1e-20
      objB[k] <- sum(-log(dpois(XB, lamB) + eps1))/(n*p)
      cat("Completed fold", k, "estimation.\n")
    }

    obj_mat[r, ] <- objB
    cat("Factor number r =", r, "\n")
  }

  cv <- rowSums(obj_mat)
  icv <- log(cv)

  return(list(TCV = cv, TICV = icv))
}

add_identifiability <- function(H, B, mu){
  # Load the irlba library
  #library(irlba)

  # Perform SVD decomposition with rank k = 10

  mu <- mu + B %*% colMeans(H)
  q <- ncol(H); n <- nrow(H)
  svdHB <- irlba((H- matrix(colMeans(H), n, q, byrow = TRUE)) %*% t(B), nv= q)
  #signB1 <- sign(svdHB$v[1,])
  signB1 = sign(apply(svdHB$v, 2, function(x) x[which.max(abs(x))]))
  H <- sqrt(n) * svdHB$u %*% Diag(signB1)

  B <- svdHB$v %*% Diag(svdHB$d[1:q]*signB1) / sqrt(n)

  return(list(H=H, B=B, mu=mu))
}




#' Estimating the Number of Factor by Eigenvalue Ratio of Natural Parameter Matrix in Generalized Factor Model.
#'
#' @param XList A list that containing an n by p matrix, where n is the number of samples, p is the number of features.
#' @param types The type of data. In Poisson factor models, the type is "poisson".
#' @param q_set The maximum number of factors for conducting ratio methods. Default as 10.
#' @param select_method The methods to conduct GFM. Default as AM.
#' @param offset Default as FALSE.
#' @param dc_eps The tolerance for convergence. Default as 1e-4.
#' @param maxIter The maximum iteration times. Defualt as 30.
#' @param verbose Default as TRUE.
#' @param parallelList Whether to use parallel. Default as NULL.
#'
#' @return The number of factors estimated by ratio methods.
#' @author Zhijing Wang

chooseFacNumber_ratio <- function(XList, types, q_set = 2: 10, select_method = c("SVR", "IC"),offset=FALSE, dc_eps=1e-4, maxIter=30,
                                  verbose = TRUE, parallelList=NULL){

  required_packages <- c(
    "foreach", "doParallel", "irlba", "doSNOW", "openxlsx",
    "countsplit", "bcv", "Matrix", "tidyr", "invgamma",
    "mvtnorm", "pracma", "pcaPP", "Rfast", "dplyr", "GFM"
  )

  # Check if packages are installed, if not prompt to install
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Please install the '", pkg, "' package: install.packages('", pkg, "')"))
    }
  }

  # Load packages while suppressing messages and warnings
  suppressPackageStartupMessages({
    for (pkg in required_packages) {
      library(pkg, character.only = TRUE)
    }
  })

  q_max <- max(q_set)
  gfm2 <- GFM::gfm(XList, types, q=q_max, algorithm = 'AM', offset=offset, dc_eps=dc_eps, maxIter=maxIter,verbose = verbose)


  svalues <- svd(gfm2$hB)$d
  svalues_use <- svalues[svalues>1e-2]
  q_max <- length(svalues_use)
  q <- which.max(svalues[-q_max] / svalues[-1])

  message('The singular value ratio  method estimates the factor number q  as ', q, '. \n')

  return(q)

}




SR <- function(eigvals, Kmax, dim, size) {
  m_prime_value <- rep(0, Kmax + 1)
  for (j in 1:(Kmax + 1)) {
    cj <- (dim - j) / size
    m_prime_value[j] <- mean(na.omit(1 / ((eigvals[-(1:j)] - eigvals[j])^2+dim^{-4/3}) ) )
  }
  m_prime_ratio <- m_prime_value[1:Kmax + 1] / m_prime_value[1:Kmax]
  Khat_SR <- which.max(m_prime_ratio)
  return(Khat_SR)
}


# ------------------------------------------------------------------------------------
# Function: Estimate population spiked eigenvalues from sample eigenvalues
# Reference: Bai and Ding (2012, RMTA), Fan et al. (2020, JASA)
# Input:
#         1. eigvals: sample eigenvalues of Pearson's correlation matrix,
#         2. Kmax: predetermined upper bound of the true number of spiked eigenvalues
# Output:
#         1. eigvals_C: Kmax estimated spiked eigenvalues
# ------------------------------------------------------------------------------------
BaiDing12_RMTA <- function(eigvals, Kmax) {
  mz <- rep(0, Kmax)
  for (kk in 1:Kmax) {
    qu <- 3 / 4
    eigvals_1 <- eigvals[-seq(1, kk, 1)]
    z0 <- qu * eigvals[kk] + (1 - qu) * eigvals[kk + 1]
    ttt <- c(eigvals_1 - eigvals[kk], z0 - eigvals[kk])
    y0 <- (length(eigvals_1)) / (n - 1)
    mz[kk] <- -(1 - y0) / eigvals[kk] + y0 * (sum(na.omit(1 / ttt)) / (p - kk))
  }
  eigvals_C <- -1 / mz
  return(eigvals_C)
}
# ------------------------------------------------------------------------
# Function: ACT method
# Reference: Fan et al. (2020, JASA)
# Input:
#         1. eigvals: sample eigenvalues of Pearson's correlation matrix
#         2. threshold: threshold for adjusted eigenvalues
#         3. Kmax: predetermined upper bound of the true number of factors
# Output:
#         1. Khat_ACT: estimated number of common factors
# ------------------------------------------------------------------------
ACT <- function(eigvals, threshold, Kmax, n ,p) {
  mz <- rep(0, Kmax)
  for (kk in 1:Kmax) {
    qu <- 3 / 4
    eigvals_1 <- eigvals[-seq(1, kk, 1)]
    z0 <- qu * eigvals[kk] + (1 - qu) * eigvals[kk + 1]
    ttt <- c(eigvals_1 - eigvals[kk], z0 - eigvals[kk])
    y0 <- (length(eigvals_1)) / (n - 1)
    mz[kk] <- -(1 - y0) / eigvals[kk] + y0 * (sum(na.omit(1 / ttt)) / (p - kk))
  }
  eigvals_C <- -1 / mz
  #eigvals_C <- BaiDing12_RMTA(eigvals, Kmax)
  temp <- eigvals_C - threshold
  temp1 <- seq(1, Kmax, 1)
  temp2 <- cbind(temp1, temp)
  Khat_ACT <- max((temp2[, 1][temp2[, 2] > 0]), 0)
  return(Khat_ACT)
}


getCurrentPath <- function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  # If script is not run from command line, get path from RStudio
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}


library(invgamma)
library(mvtnorm)
library(pracma)
library(pcaPP)
library(Rfast)


# -------------------------------------------------------------------
# Function: Calculate Multivariate Kendall's tau matrix (MKen)
# Input:
#         1. X: p x n, data matrix (p = dimension, n = sample size)
# Output:
#         1. TK: p x p, MKen of X
# -------------------------------------------------------------------
SK <- function(X) {
  p <- nrow(X)
  n <- ncol(X)
  TK <- matrix(0, p, p)
  for (i in 1:(n - 2)) {
    TT <- X[, i] - X[, (i + 1):n]
    V <- sqrt(colsums(TT^2))
    TT_spatial_sign <- sweep(TT, 2, V, FUN = "/")
    TT <- tcrossprod(TT_spatial_sign)
    TK <- TK + TT
  }
  TT <- X[, n - 1] - X[, n]
  TK <- TK + TT %*% t(TT) / sum(TT^2)
  TK <- 2 / (n * (n - 1)) * TK
  return(TK)
}



# ---------------------------------------------------------------------------
# Function: ER method
# Reference: Lam and Yao (2012, AoS), Ahn and Horenstein (2013, Econometrica)
# Input:
#         1. eigvals: sample eigenvalues of some kinds of matrices, e.g.,
#                     sample covariance matrix
#         2. Kmax: predetermined upper bound of the true number of factors
# Output:
#         1. Khat_ER: estimated number of common factors
# ---------------------------------------------------------------------------
ER <- function(eigvals, Kmax,n ,p) {
  eigvals_lag1 <- c(eigvals[-1], 0)
  ER <- (eigvals_lag1 / eigvals)[-n]
  Khat_ER <- which.min(ER[1:Kmax])
  return(Khat_ER)
}


# ---------------------------------------------------------------------------
# Function: GR method
# Reference: Ahn and Horenstein (2013, Econometrica)
# Input:
#         1. eigvals: sample eigenvalues of some kinds of matrices, e.g.,
#                     sample covariance matrix
#         2. Kmax: predetermined upper bound of the true number of factors
# Output:
#         1. Khat_GR: estimated number of common factors
# ---------------------------------------------------------------------------
GR <- function(eigvals, Kmax, dim, size) {
  hm1 <- hm2 <- hm3 <- matrix(rep(0, dim * Kmax), dim, Kmax)
  for (i in 1:Kmax) {
    for (j in i:dim) {
      hm1[j, i] <- 1
    }
  }
  for (i in 1:Kmax) {
    for (j in (i + 1):dim) {
      hm2[j, i] <- 1
    }
  }
  for (i in 1:Kmax) {
    for (j in (i + 2):dim) {
      hm3[j, i] <- 1
    }
  }
  gr1 <- eigvals %*% hm1
  gr2 <- eigvals %*% hm2
  gr3 <- eigvals %*% hm3
  GR <- log(gr1 / gr2) / log(gr2 / gr3)
  Khat_GR <- which.max(GR)
  return(Khat_GR)
}



# ------------------------------------------------------------------------
# Function: ED method
# Reference: Onatski (2010, RES)
# Input:
#         1. eigvals: sample eigenvalues of sample covariance matrix
#         2. Kmax: predetermined upper bound of the true number of factors
#         3. iter: number of iteration
# Output:
#         1. Khat_ED: estimated number of common factors
# ------------------------------------------------------------------------
ED <- function(eigvals, Kmax, iter) {
  iter <- 4
  ED_seq <- eigvals[1:Kmax] - eigvals[2:(Kmax + 1)]
  noj <- Kmax + 1
  for (j in 1:iter) {
    y <- eigvals[noj:(noj + iter)]
    x <- (noj + seq(-1, (iter - 1), 1))^(2 / 3)
    beta <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
    delta <- 2 * abs(beta)
    # min factor number is 1
    # Khat_ED <- ifelse(max(ED_seq - delta) < 0, 0, max(which(ED_seq >= delta)))
    Khat_ED <- ifelse(max(ED_seq - delta) < 0, 1, max(which(ED_seq >= delta)))
    noj <- Khat_ED + 1
  }
  return(Khat_ED)
}


# ------------------------------------------------------------------------
# Function: NE method
# Reference: Nadakuditi and Edelman (2008, IEEE Transactions on Signal Processing)
# Input:
#         1. eigvals: sample eigenvalues of Pearson's correlation matrix
#         2. Kmax: predetermined upper bound of the true number of factors
#         3. dim: dimension
#         4. size: sample size
# Output:
#         1. Khat_SR: estimated number of common factors
# ------------------------------------------------------------------------
NE <- function(eigvals, dim, size, maxfactor=min(p, n)) {
  p <- dim
  n <- size
  t <- rep(0, maxfactor)
  ## min number of factor is one
  for (i in 1:(maxfactor)) {
    temp <- p * ((p - i) * (sum(eigvals[(i + 1):p]^2) / sum(eigvals[(i + 1):p])^2) - 1 - p / n) - p / n
    t[i] <- temp^2 * (n / p)^2 / 4 + 2 * (i + 1)
  }
  Khat_NE <- which.min(t)
  # return(t)
  return(Khat_NE)
}










