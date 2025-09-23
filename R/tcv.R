

#' Enforce Identifiability Constraints on Factor Model Components
#'
#' Post-processes the factor scores (H), loadings (B), and intercept (mu)
#' to ensure a unique solution by applying SVD-based rotation.
#' This typically enforces orthogonality constraints.
#'
#' @param H A numeric matrix of factor scores (n x q).
#' @param B A numeric matrix of factor loadings (p x q).
#' @param mu A numeric vector for the intercept/mean term.
#' @return A list containing the transformed H, B, and mu that satisfy identifiability constraints.
add_identifiability <- function(H, B, mu){
  # Load the irlba library
  #library(irlba)

  # Perform SVD decomposition with rank k = 10

  mu <- mu + B %*% colMeans(H)
  q <- ncol(H); n <- nrow(H)
  svdHB <- irlba::irlba((H- matrix(colMeans(H), n, q, byrow = TRUE)) %*% t(B), nv= q)
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
#' @param q_set The maximum number of factors for conducting ratio methods. Default as 5.
#' @param select_method The methods to conduct GFM. Default as AM.
#' @param offset Default as FALSE.
#' @param dc_eps The tolerance for convergence. Default as 1e-4.
#' @param maxIter The maximum iteration times. Defualt as 30.
#' @param verbose Default as FALSE
#' @param parallelList Whether to use parallel. Default as NULL.
#'
#' @return The number of factors estimated by ratio methods.
chooseFacNumber_ratio <- function(XList, types, q_set = 1: 5, select_method = c("SVR", "IC"),offset=FALSE, dc_eps=1e-4, maxIter=30,
                                  verbose = FALSE, parallelList=NULL){


  q_max <- max(q_set)
  gfm2 <- GFM::gfm(XList, types, q=q_max, algorithm = 'AM', offset=offset, dc_eps=dc_eps, maxIter=maxIter,verbose = verbose)


  svalues <- svd(gfm2$hB)$d
  svalues_use <- svalues[svalues>1e-2]
  q_max <- length(svalues_use)
  q <- which.max(svalues[-q_max] / svalues[-1])

  message('The singular value ratio (SVR) method estimates the factor number q  as ', q, '. \n')

  return(q)

}



#' Create a Diagonal Matrix Robustly
#'
#' @param vec A numeric vector to be placed on the diagonal of the matrix.
#'
#' @return A square matrix where the main diagonal contains the elements of `vec`.
#' @noRd
Diag <- function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  } else {
    y <- matrix(vec, 1, 1)
  }
  return(y)
}



################## multiDT Function ##################

#' Perform Thinning Cross-Validation to Select Factor Number
#'
#' This function implements a K-fold cross-validation scheme based on data thinning
#' (count splitting) to determine the optimal number of factors for a Poisson matrix factorization model.
#'
#' @param x A numeric matrix of count data (n x p).
#' @param K An integer, the number of folds for cross-validation. Default is 5.
#' @param rmax An integer, the maximum number of factors to test. Default is 8.
#' @return A list containing two elements:
#'         - TCV: A numeric vector of total cross-validation error for each number of factors.
#'         - TICV: A numeric vector of the natural logarithm of TCV.
#' @importFrom stats runif dpois
#' @examples
#' # 1. Set parameters for data generation
#' # Use smaller dimensions for a quick example
#' n <- 50 # Number of samples
#' p <- 30 # Number of features
#' true_q <- 2  # True number of factors
#'
#' # 2. Generate data from a Poisson factor model
#' set.seed(123) # For reproducibility
#'
#' # Factor matrix (scores)
#' FF <- matrix(rnorm(n * true_q), nrow = n, ncol = true_q)
#'
#' # Loading matrix
#' BB <- matrix(runif(p * true_q, min = -1, max = 1), nrow = p, ncol = true_q)
#'
#' # Intercept term
#' a <- runif(p, min = 0, max = 1)
#'
#' # Enforce identifiability for a unique generating model
#' FF0 <- add_identifiability(FF, BB, a)$H
#' BB0 <- add_identifiability(FF, BB, a)$B
#' alpha <- add_identifiability(FF, BB, a)$mu
#'
#' # Calculate the mean matrix (lambda) with some noise
#' lambda <- exp(FF0 %*% t(BB0) + rep(1, n) %*% t(alpha) + matrix(rnorm(n*p, 0, 0.5), n, p))
#'
#' # Generate the final count data matrix 'x'
#' x <- matrix(rpois(n * p, lambda = as.vector(lambda)), nrow = n, ncol = p)
#'
#' # 3. Run multiDT to find the best number of factors
#' # Use small K and rmax for a quick example run
#' cv_results <- multiDT(x, K = 2, rmax = 4)
#'
#' # 4. Print results and select the best 'r' based on the minimum TCV
#' print(cv_results$TCV)
#' best_r <- which.min(cv_results$TCV)
#'
multiDT=function(x, K = 5, rmax=8){
  n=nrow(x)
  p=ncol(x)
  xsplit=countsplit::countsplit(x, folds = K)
  obj_mat=matrix(0, nrow = rmax, ncol = K)
  co_mat=matrix(0, nrow = rmax, ncol = K)
  theta0 = cbind(rep(1, n), matrix(runif(n*rmax, 0,1), n,rmax))
  A0 = matrix(runif(p*(rmax+1), -1,1), p, (rmax+1))
  for (r in (1:rmax)) {
    QB=matrix()
    objB=c()
    co_es=c()
    for (k in (1:K)) {
      XA=x-as.matrix(xsplit[[k]])
      XB=x-XA

      lamB = tryCatch({
        tra=suppressMessages(GFM::gfm(list(XA), types = c("poisson"), maxIter = 30, q = r, algorithm = "AM"))
        lamA=exp(cbind(tra$hH,1)%*%t(cbind(tra$hB, tra$hmu)))
        lamA/(K-1)
      }, error = function(e){
        true_matrix <- matrix(TRUE, nrow = n, ncol = p)
        jic.res= confirm_CJMLE_poisson_cpp(XA,  true_matrix, theta0[,1:(r+1)], A0[,1:(r+1)], matrix(TRUE, p, (r+1)), C = Inf, tol = .01/(n*p))
        lamA=exp((jic.res$theta)%*%t(jic.res$A))
        lamB=lamA/(K-1)
      })

      eps1 <- 1e-20
      objB[k] <- sum(-log(dpois(XB, lamB)+eps1))/(n*p)
    }
    obj_mat[r,]=objB
    #cat("factor number r = ", r, "\n")
  }
  cv=rowSums(obj_mat)
  icv=log(cv)
  return(list(TCV=cv, TICV=icv))
}






