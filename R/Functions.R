#' Make a p-Matrix
#'
#' This function creates a 2I x 2I matrix where the diagonals equal zero, and the
#' off-diagonal elements (i, j) contain the probability the ith observation has
#' Z = max{Z_i, Z_j} and the jth observation has Z = min{Z_i, Z_j}. 
#'
#' @param Z a I-length vector of treatment values
#' @param X a I x k matrix of covariate values
#' @param method currently takes either "linCDE" or "RFCDE" as values. Determines
#' which conditional density estimation method to use to estimate the generalized
#' propensity score. Defaults to "linCDE"
#' @return a 2I x 2I numeric matrix
#' @examples
#' set.seed(12345)
#' X <- rnorm(1000, 0, 5)
#' Z <- X + rnorm(1000, 0, (1+sqrt(abs(X))))
#' make_pmatrix(Z, X)
#' @export
make_pmatrix <- function(Z, X, method = "linCDE"){
  if(method %in% c("linCDE", "RFCDE") == FALSE){
    stop("Failed to select valid CDE method")
  }
  
  X <- data.frame(X)
  p_matrix <- matrix(0, nrow = length(Z), ncol = length(Z))

  if(method == "linCDE"){
    
    propensity_matrix <- matrix(0, nrow = length(Z), ncol = length(Z))
    model <- LinCDE::LinCDE.boost(y = Z, X = X, terminalSize = 20, verbose = FALSE, centering = TRUE, centeringMethod = "linearRegression")
  
    for(i in 1:nrow(propensity_matrix)){
      propensity_matrix[i, ] <- predict(model, X = X[i, ], y = Z)
    }
    
    for(i in 1:(length(Z)-1)){
      for(j in (i+1):(length(Z))){
        
        observed_treatment <- propensity_matrix[i, i] * propensity_matrix[j, j]
        opposing_treatment <- propensity_matrix[j, i] * propensity_matrix[i, j]
        p_observed <- observed_treatment / (observed_treatment + opposing_treatment)
        p_matrix[i, j] <- ifelse(Z[i] > Z[j], p_observed, 1-p_observed)
        p_matrix[j, i] <- p_matrix[i, j]
        
      }
    }  
  }
  
  if(method == "RFCDE"){
    
    model <- RFCDE(data[, 3:ncol(data)], data[, 1])
    predicted <- predict(model, as.matrix(data[, 3:ncol(data)]), "CDE", data[, 1])
    same_propensity_vector <- 1:nrow(data)
    
    for(i in 1:nrow(data)){
      same_propensity_vector[i] <- predicted[i, i]
    }
    
    observed_treatment_matrix <- same_propensity_vector %*% base::t(same_propensity_vector)
    p_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
    
    for(i in 1:(nrow(data)-1)){
      for(j in (i+1):nrow(data)){
        
        opposing_treatment <- predicted[j, i] * predicted[i, j]
        p_observed <- observed_treatment_matrix[i, j] / (observed_treatment_matrix[i, j] + opposing_treatment)
        p_matrix[i, j] <- ifelse(Z[i] > Z[j], p_observed, 1-p_observed)
        p_matrix[j, i] <- p_matrix[i, j]
      }
    }
  }

  return(p_matrix)
}

#' non-bipartite matching with treatment assignment caliper
#'
#' This function creates a I x 2 dataframe containing the indices of observations
#' that form our set of matched pairs. It uses the nbpMatch package along with 
#' a p-atrix in order to create I matched pairs using a treatment assignment
#' caliper. 
#'
#' @param Z a I-length vector of treatment values
#' @param X a I x k matrix of covariate values
#' @param pmat a 2I x 2I matrix where the diagonals equal zero, and the
#' off-diagonal elements (i, j) contain the probability the ith observation has
#' Z = max{Z_i, Z_j} and the jth observation has Z = min{Z_i, Z_j}. We can create
#' a p-matrix using the make_pmatrix function.
#' @param xi a number in the range [0, 0.5], the cutoff related to the treatment
#' assignment probability caliper.
#' @param M an integer determining the penalty of the treatment assignment
#' probability caliper. If a potential matched pair between observations i and j 
#' has treatment assignment probability less than xi or greater than 1-xi, we add 
#' M to the distance matrix in the (i, j) and (j, i) entry.
#' @return I x 2 dataframe
#' @examples
#' set.seed(12345)
#' X <- rnorm(1000, 0, 5)
#' Z <- X + rnorm(1000, 0, (1+sqrt(abs(X))))
#' pmat <- make_pmatrix(Z, X)
#' nbp_caliper(Z, X, pmat, xi = 0.1, M = 10000)
#' @export
nbp_caliper <- function(Z, X, pmat, xi, M){
  
  X <- data.frame(X)
  
  # matrix encoding our caliper penalty
  caliper_matrix <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  
  for(i in 1:(nrow(X)-1)){
    for(j in (i+1):nrow(X)){
      
      # if probability of observed treatment assignment is extreme, enforce
      # infinite penalty on distance matrix
      caliper_matrix[i, j] <- ifelse(pmat[i, j] > (1-xi) | pmat[i, j] < xi, Inf, 0)
      caliper_matrix[j, i] <- caliper_matrix[i, j]
      
    }
  }
  
  # making mahalanobis distance matrix
  
  Matching_data <- cbind(Z, X)
  test.dist <- gendistance(Matching_data, idcol = 1)
  test.mdm <- distancematrix(test.dist)
  
  # applying caliper to mahalnobis distance matrix
  
  caliper_dist_matrix <- distancematrix(test.mdm + caliper_matrix)
  test.match <- nonbimatch(caliper_dist_matrix)
  matches <- test.match$matches
  matches$group <- ifelse(as.numeric(matches$Group1.ID) < as.numeric(matches$Group2.ID), "L", "H")
  high_groups <- matches[which(matches$group == "H"), 2]
  low_groups <- matches[which(matches$group == "L"), 2]
  
  balanced_pairs <- as.data.frame(matrix(nrow = 0, ncol = 2))
  high_data <- matches[which(matches$group == "H"), ]
  balanced_pairs <- high_data[, c(2, 4)]
  names(balanced_pairs) <- c("High", "Low")
  
  return(balanced_pairs)
  
}

#' Classic Neyman Sample Average Treatment Effect Estimator
#'
#' This function estimates the sample average treatment effect for a set of
#' matched pairs using the classic Neyman estimator.
#'
#' @param Y a I-length vector of outcome values
#' @param Z a I-length vector of treatment values
#' @param pairs a I x 2 dataframe containing the indices of observations
#' that form our set of matched pairs.
#' @return I x 2 dataframe
#' @examples
#' set.seed(12345)
#' X <- rnorm(1000, 0, 5)
#' Z <- X + rnorm(1000, 0, (1+sqrt(abs(X))))
#' Y <- X + Z + rnorm(1000, 0, 0.5)
#' pmat <- make_pmatrix(Z, X)
#' pairs <- nbp_caliper(Z, X, pmat, xi = 0.1, M = 10000)
#' classic_Neyman(Y, Z, pairs)
#' @export
classic_Neyman <- function(Y, Z, pairs){
  
  numerator <- ifelse(Z[pairs[, 1]] > Z[pairs[, 2]], Y[pairs[, 1]] - Y[pairs[, 2]], Y[pairs[, 2]] - Y[pairs[, 1]])
  denominator <- ifelse(Z[pairs[, 1]] > Z[pairs[, 2]], Z[pairs[, 1]] - Z[pairs[, 2]], Z[pairs[, 2]] - Z[pairs[, 1]])
  return(sum(numerator) / sum(denominator))
  
}

#' Bias-corrected Neyman Sample Average Treatment Effect Estimator
#'
#' This function estimates the sample average treatment effect for a set of
#' matched pairs using the bias-corrected Neyman estimator.
#'
#' @param Y a I-length vector of outcome values
#' @param Z a I-length vector of treatment values
#' @param pairs a I x 2 dataframe containing the indices of observations
#' that form our set of matched pairs.
#' @param pmat a 2I x 2I matrix where the diagonals equal zero, and the
#' off-diagonal elements (i, j) contain the probability the ith observation has
#' Z = max{Z_i, Z_j} and the jth observation has Z = min{Z_i, Z_j}. We can create
#' a p-matrix using the make_pmatrix function.
#' @param xi a number in the range [0, 0.5], the cutoff related to the treatment
#' assignment probability caliper.
#' @return I x 2 dataframe
#' @examples
#' set.seed(12345)
#' X <- rnorm(1000, 0, 5)
#' Z <- X + rnorm(1000, 0, (1+sqrt(abs(X))))
#' Y <- X + Z + rnorm(1000, 0, 0.5)
#' pmat <- make_pmatrix(Z, X)
#' pairs <- nbp_caliper(Z, X, pmat, xi = 0.1, M = 10000)
#' biasCorrected_Neyman(Y, Z, pairs, pmat, xi = 0.1)
#' @export
biasCorrected_Neyman <- function(Y, Z, pairs, pmat, xi){

  p <- diag(pmat[pairs[, 1], pairs[, 2]])
  p[p < xi] <- xi
  p[p > (1-xi)] <- 1-xi
  numerator <- ifelse(Z[pairs[, 1]] > Z[pairs[, 2]], (p^(-1))*(Y[pairs[, 1]] - Y[pairs[, 2]]), ((1-p)^(-1))*(Y[pairs[, 2]] - Y[pairs[, 1]]))
  denominator <- ifelse(Z[pairs[, 1]] > Z[pairs[, 2]], 2*(Z[pairs[, 1]] - Z[pairs[, 2]]), 2*(Z[pairs[, 2]] - Z[pairs[, 1]]))
  return(sum(numerator) / sum(denominator))
  
}

#' Generate some example data
#'
#' This function creates some example data using the data generation process described
#' in simulation 1 of Frazier, Heng, Zhou (2024). The dataframe contains a treatment
#' variable Z, outcome variable Y, and five covariate X1,...,X5.
#'
#' @param N number of observations to simulate
#' @return an N x 7 matrix containing treatment, outcome, and covariates.
#' @examples
#' generate_data_dose(N = 100)
#' @export
generate_data_dose <- function(N){
  # integer N: determines number of observations generated 
  
  # Generate covariates
  x1 = rnorm(N, 0, 1); x2 = rnorm(N, 0, 1); x3 = rnorm(N, 0, 1); x4 = rnorm(N, 0, 1); x5 = rnorm(N, 0, 1);
  
  # Generate treatment variable
  error1 = rnorm(N, 0, 1)
  z = x1 + x2^2 + abs(x3 * x4) + ifelse(x4 > 0, 1, 0) + log(1 + abs(x5)) + error1
  
  # Generate outcome variable
  error2 = rnorm(N, 0, 3)
  y = z + 0.3*x1*z + 0.2*x3^3*z + exp(abs(x2 - x4)) - sin(x5) + error2
  
  # Place variables in a dataframe
  data = data.frame("Z" = z, "Y" = y, "X1" = x1, "X2" = x2, 
                    "X3" = x3, "X4" = x4, "X5" = x5)
  
  return(data)
}

#' Covariate-Adjusted Variance Estimation
#'
#' This function calculates the covariate-adjusted conservative variance estimator
#' For the (classic or bias-corrected) Neyman estimator
#'
#' @param Y a I-length vector of outcome values
#' @param Z a I-length vector of treatment values
#' @param X a I x k matrix of covariate values
#' @param pairs a I x 2 dataframe containing the indices of observations
#' that form our set of matched pairs.
#' @param pmat a 2I x 2I matrix where the diagonals equal zero, and the
#' off-diagonal elements (i, j) contain the probability the ith observation has
#' Z = max{Z_i, Z_j} and the jth observation has Z = min{Z_i, Z_j}. We can create
#' a p-matrix using the make_pmatrix function.
#' @param xi a number in the range [0, 0.5], the cutoff related to the treatment
#' assignment probability caliper.
#' @param Q an arbitrary I x L matrix, where L < I
#' @return a 2I x 2I numeric matrix
#' @examples
#' set.seed(12345)
#' X <- rnorm(1000, 0, 5)
#' Z <- X + rnorm(1000, 0, (1+sqrt(abs(X))))
#' Y <- X + Z + rnorm(1000, 0, 0.5)
#' pmat <- make_pmatrix(Z, X)
#' pairs <- nbp_caliper(Z, X, pmat, xi = 0.1, M = 10000)
#' covAdj_variance(Y, Z, X, pairs, pmat, xi = 0.1)
#' @export
covAdj_Variance <- function(Y, Z, X, pairs, pmat, xi, Q){
  
  X <- data.frame(X)
  if(length(unique(c(length(Z), length(Y), nrow(X)))) != 1){
    stop("One of X, Y, or Z is not the right length.")
  }
  
  if(missing(Q) == TRUE){
    if(missing(X) == TRUE){
      Q <- rep(1, length(Y))
    }
    if(missing(X) == FALSE){
      Q <- (X[pairs[, 1], ] + X[pairs[, 2], ])/2
    }
  }
  Q <- matrix(Q)
  singular_columns <- which(colSums(abs(Q)) == 0)
  
  if(length(singular_columns) != 0){
    Q <- Q[, -singular_columns] 
  }
  
  if(dim(Q)[1] < dim(Q)[2]){
    stop("Q matrix has more columns than rows.")
  }
  
  hatQ <- Q %*% solve(base::t(Q)%*%Q) %*% base::t(Q)
  
  if(missing(pmat)){
    pmat <- matrix(0.5, nrow = length(Y), ncol = length(Z))
  }
  
  pvec <- diag(pmat[pairs[, 1], pairs[, 2]])
  pvec[pvec < xi] <- xi
  pvec[pvec > (1-xi)] <- (1-xi)
  yvec <- ifelse(Z[pairs[, 1]] > Z[pairs[, 2]], (pvec^(-1))*(Y[pairs[, 1]] - Y[pairs[, 2]])/sqrt(1-diag(hatQ)), ((1-pvec)^(-1))*(Y[pairs[, 2]] - Y[pairs[, 1]])/sqrt(1-diag(hatQ)))
  denominator <- sum(ifelse(Z[pairs[, 1]] > Z[pairs[, 2]], 4*(Z[pairs[, 1]] - Z[pairs[, 2]])^2, 4*(Z[pairs[, 2]] - Z[pairs[, 1]])^2))
  cov_adj_var <- (yvec %*% (diag(1, nrow(hatQ)) - hatQ) %*% yvec)/denominator
  return(cov_adj_var)
}
