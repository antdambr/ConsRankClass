#'Tau_x rank correlation coefficient for vectors
#'
#'Tau_x rank correlation coefficient for large vectors  
#'
#' @param X A vector of length n
#' @param Y A vector of length n
#' 
#' @return The tau_x rank correlation coefficient between the 2 vectors
#' 
#' @export


vecTaux <- function(X,Y){
  
  #Tau_X rank correlation coefficient for vectorized vectors
  
  n <- length(X)
  if (n != length(Y)) stop("X and Y must have the same length")
  
  # Matrices (no loop)
  X_mat <- outer(X, X, "<=")
  Y_mat <- outer(Y, Y, "<=")
  
  # Sum (avoiding diagonal)
  sum_pairs <- sum((2*X_mat - 1) * (2*Y_mat - 1) * !diag(n))
  
  Tx <- sum_pairs / (n * (n - 1))
  return(Tx)
  
}